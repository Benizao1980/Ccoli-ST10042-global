#!/usr/bin/env python3
"""Apply the validated INNUca-like post-assembly filter to one Azevedo genome.

This is the scale-up worker used by the Slurm array. It deliberately performs
only local assembly filtering/QC; PubMLST REST queries are kept out of the array
so we do not send a large burst of simultaneous requests.

Validated workflow:
  SPAdes contigs
  -> retain contigs >=200 bp with SPAdes k-mer coverage >=2
  -> map the same fastp-trimmed paired reads back with Bowtie2
  -> retain contigs with mean mapped depth >= max(10x, assembly mean / 3)
  -> write filtered FASTA + per-contig coverage + per-sample QC JSON
"""

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path

import pandas as pd

SPADES_COV_RE=re.compile(r"(?:^|_)cov_([0-9.eE+-]+)(?:_|$)")


def run(cmd, stdout=None):
    print("+"," ".join(map(str,cmd)),flush=True)
    subprocess.run(cmd,check=True,stdout=stdout)


def read_fasta(path):
    records=[]
    header=None
    seq=[]
    with open(path,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header,"".join(seq)))
                header=line[1:].split()[0]
                seq=[]
            elif header is not None:
                seq.append(line.strip())
        if header is not None:
            records.append((header,"".join(seq)))
    return records


def spades_cov(header):
    m=SPADES_COV_RE.search(header)
    return float(m.group(1)) if m else None


def write_fasta(records,path):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,"w",encoding="utf-8") as out:
        for h,s in records:
            out.write(f">{h}\n")
            for i in range(0,len(s),80):
                out.write(s[i:i+80]+"\n")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",default="results/azevedo_spades_manifest.tsv")
    p.add_argument("--index",type=int,required=True)
    p.add_argument("--assembly-dir",default="data/europe_assemblies")
    p.add_argument("--trim-dir",default="data/europe_trimmed")
    p.add_argument("--out-dir",default="data/europe_assemblies_innuca_like_trimmed")
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--min-length",type=int,default=200)
    p.add_argument("--min-spades-cov",type=float,default=2.0)
    p.add_argument("--keep-mapping-files",action="store_true")
    a=p.parse_args()

    for exe in ("bowtie2","bowtie2-build","samtools"):
        if not shutil.which(exe):
            raise SystemExit(f"ERROR {exe} not found; activate st10042-assembly")

    m=pd.read_csv(a.manifest,sep="\t",dtype=str).fillna("")
    if a.index < 0 or a.index >= len(m):
        raise SystemExit(f"ERROR index {a.index} outside manifest range 0..{len(m)-1}")

    r=m.iloc[a.index]
    sample=r["Strain_ID"]
    source=Path(a.assembly_dir)/sample/"contigs.fasta"
    r1=Path(a.trim_dir)/sample/f"{sample}_R1.fastq.gz"
    r2=Path(a.trim_dir)/sample/f"{sample}_R2.fastq.gz"

    for path,label in ((source,"SPAdes assembly"),(r1,"trimmed R1"),(r2,"trimmed R2")):
        if not path.exists() or path.stat().st_size==0:
            raise SystemExit(f"ERROR missing/empty {label}: {path}")

    od=Path(a.out_dir)/sample
    od.mkdir(parents=True,exist_ok=True)
    filtered=od/"contigs.innuca_like.fasta"
    coverage_tsv=od/"contig_coverage.tsv"
    qc_json=od/"qc.json"

    # A completed JSON + non-empty FASTA is the resumability marker.
    if qc_json.exists() and filtered.exists() and filtered.stat().st_size>0:
        qc=json.loads(qc_json.read_text(encoding="utf-8"))
        print(
            f"SKIP complete {sample}: filtered_contigs={qc.get('filtered_contigs')} "
            f"filtered_bp={qc.get('filtered_bp')}"
        )
        return

    print(f"Sample {sample} (array index {a.index})",flush=True)

    records=read_fasta(source)
    pre=[]
    no_cov_header=0
    for h,s in records:
        cov=spades_cov(h)
        if cov is None:
            no_cov_header += 1
        if len(s) < a.min_length:
            continue
        if cov is not None and cov < a.min_spades_cov:
            continue
        pre.append((h,s))

    if not pre:
        raise SystemExit(f"ERROR no contigs remain after stage-1 filtering for {sample}")

    stage1=od/"contigs.kmercov_filtered.fasta"
    write_fasta(pre,stage1)
    print(
        f"stage1: {len(pre)}/{len(records)} contigs "
        f"(len>={a.min_length}, SPAdes_cov>={a.min_spades_cov}); "
        f"headers_without_cov={no_cov_header}",
        flush=True
    )

    prefix=str(od/"bt2")
    run(["bowtie2-build","--quiet",str(stage1),prefix])

    sam=od/"reads.sam"
    bam=od/"reads.sorted.bam"
    with open(sam,"w") as sfh:
        run([
            "bowtie2","--quiet","--very-sensitive-local",
            "--threads",str(a.threads),"--fr","-I","0","-X","2000",
            "-x",prefix,"-1",str(r1),"-2",str(r2)
        ],stdout=sfh)
    run(["samtools","sort","-@",str(a.threads),"-o",str(bam),str(sam)])
    run(["samtools","index",str(bam)])
    sam.unlink(missing_ok=True)

    lengths={h:len(s) for h,s in pre}
    covsum={h:0 for h,_ in pre}

    # Stream depth rather than writing ~1.6 Mb of per-base depth per sample.
    proc=subprocess.Popen(
        ["samtools","depth","-a",str(bam)],
        stdout=subprocess.PIPE,
        text=True
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        chrom,_pos,dep=line.rstrip("\n").split("\t")[:3]
        if chrom in covsum:
            covsum[chrom] += int(dep)
    rc=proc.wait()
    if rc != 0:
        raise SystemExit(f"ERROR samtools depth failed for {sample} with exit {rc}")

    total_positions=sum(lengths.values())
    total_depth=sum(covsum.values())
    assembly_mean=(total_depth/total_positions) if total_positions else 0.0
    threshold=max(10.0,assembly_mean/3.0)

    keep=[]
    per_contig=[]
    for h,s in pre:
        mean=(covsum[h]/lengths[h]) if lengths[h] else 0.0
        kept=mean >= threshold
        per_contig.append({
            "contig":h,
            "length":lengths[h],
            "spades_kmer_cov":spades_cov(h),
            "mapped_mean_depth":mean,
            "kept":kept,
        })
        if kept:
            keep.append((h,s))

    if not keep:
        raise SystemExit(f"ERROR no contigs remain after mapped-depth filtering for {sample}")

    write_fasta(keep,filtered)

    pd.DataFrame(per_contig).to_csv(coverage_tsv,sep="\t",index=False)

    flag=subprocess.run(
        ["samtools","flagstat",str(bam)],
        check=True,text=True,capture_output=True
    ).stdout.splitlines()
    mapline=next((x for x in flag if " mapped (" in x),"")

    qc={
        "array_index":int(a.index),
        "Strain_ID":sample,
        "Country":r.get("Country",""),
        "run_accession":r.get("run_accession",""),
        "Figure4_cluster":r.get("Figure4_cluster",""),
        "source_contigs":len(records),
        "source_bp":sum(len(s) for _,s in records),
        "stage1_contigs":len(pre),
        "stage1_bp":sum(len(s) for _,s in pre),
        "filtered_contigs":len(keep),
        "filtered_bp":sum(len(s) for _,s in keep),
        "assembly_mean_depth":assembly_mean,
        "coverage_threshold":threshold,
        "mapped_flagstat_line":mapline.strip(),
        "min_length":a.min_length,
        "min_spades_cov":a.min_spades_cov,
        "filtered_fasta":str(filtered),
        "coverage_tsv":str(coverage_tsv),
    }
    qc_json.write_text(json.dumps(qc,indent=2)+"\n",encoding="utf-8")

    print(
        f"filtered: {len(keep)}/{len(pre)} contigs; "
        f"{qc['filtered_bp']} bp; mean_depth={assembly_mean:.2f}x; "
        f"threshold={threshold:.2f}x",
        flush=True
    )
    print(mapline.strip(),flush=True)
    print(f"Wrote {filtered}",flush=True)
    print(f"Wrote {qc_json}",flush=True)

    if not a.keep_mapping_files:
        bam.unlink(missing_ok=True)
        Path(str(bam)+".bai").unlink(missing_ok=True)
        for suffix in (
            ".1.bt2",".2.bt2",".3.bt2",".4.bt2",
            ".rev.1.bt2",".rev.2.bt2",
            ".1.bt2l",".2.bt2l",".3.bt2l",".4.bt2l",
            ".rev.1.bt2l",".rev.2.bt2l",
        ):
            Path(prefix+suffix).unlink(missing_ok=True)


if __name__=="__main__":
    main()
