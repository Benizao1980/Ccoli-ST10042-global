#!/usr/bin/env python3
"""Pilot an INNUca-like post-assembly coverage filter on known overlap isolates.

INNUca's assembly workflow filters low-kmer-coverage contigs, maps reads back to
the assembly, and keeps contigs whose mean mapped depth is at least max(10x,
one-third of the assembly-wide mean depth). This script approximates that
post-assembly filtering on our fastp+SPAdes assemblies, then re-queries PubMLST
cgMLST-v2 scheme 8.

This is a validation pilot on the known Azevedo/PubMLST overlaps before scaling
the filter to all Azevedo assemblies.
"""

import argparse
import base64
import csv
import json
import re
import shutil
import subprocess
import time
from pathlib import Path

import pandas as pd
import requests

BASE="https://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef"
SCHEME_ID=8
EXPECTED_LOCI=1142
SPADES_COV_RE=re.compile(r"(?:^|_)cov_([0-9.eE+-]+)(?:_|$)")


def run(cmd, stdout=None):
    print("+"," ".join(map(str,cmd)),flush=True)
    subprocess.run(cmd,check=True,stdout=stdout)


def get_json(url,retries=4):
    last=None
    for i in range(1,retries+1):
        try:
            r=requests.get(url,timeout=120)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            last=e
            if i<retries:
                time.sleep(2**(i-1))
    raise RuntimeError(f"GET failed {url}: {last}")


def post_json(url,payload,retries=3):
    last=None
    for i in range(1,retries+1):
        try:
            r=requests.post(url,json=payload,timeout=900)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            last=e
            if i<retries:
                time.sleep(2**i)
    raise RuntimeError(f"POST failed {url}: {last}")


def scheme_loci():
    obj=get_json(f"{BASE}/schemes/{SCHEME_ID}/loci")
    return {str(u).rstrip("/").split("/")[-1] for u in (obj.get("loci",[]) or [])}


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


def scheme_query(path,cache):
    cache=Path(cache)
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    data=base64.b64encode(Path(path).read_bytes()).decode("ascii")
    res=post_json(
        f"{BASE}/schemes/{SCHEME_ID}/sequence",
        {"base64":True,"details":False,"sequence":data},
    )
    cache.parent.mkdir(parents=True,exist_ok=True)
    cache.write_text(json.dumps(res,indent=2),encoding="utf-8")
    return res


def summarise_scheme(res,loci):
    hits=res.get("exact_matches",{}) or {}
    any_loci=set(hits)
    unique=0
    ambiguous=0
    extra=0
    for locus,vals in hits.items():
        ids={str(v.get("allele_id")) for v in (vals or []) if v.get("allele_id") is not None}
        if len(ids)==1:
            unique += 1
        elif len(ids)>1:
            ambiguous += 1
            extra += len(ids)-1
    return {
        "exact_loci_any":len(any_loci),
        "unique_exact_loci":unique,
        "ambiguous_loci":ambiguous,
        "extra_allele_hits":extra,
        "loci_without_exact_hit":len(loci-any_loci),
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--combined",default="results/combined_st10042_provisional.tsv")
    p.add_argument("--spades-manifest",default="results/azevedo_spades_manifest.tsv")
    p.add_argument("--assembly-dir",default="data/europe_assemblies")
    p.add_argument("--trim-dir",default="data/europe_trimmed")
    p.add_argument("--work-dir",default="data/europe_assemblies_innuca_like_trimmed")
    p.add_argument("--cache-dir",default="results/pubmlst_v2_innuca_like_trimmed_raw")
    p.add_argument("--out",default="results/pubmlst_v2_innuca_like_trimmed_overlap.tsv")
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--min-length",type=int,default=200)
    p.add_argument("--min-spades-cov",type=float,default=2.0)
    a=p.parse_args()

    for exe in ("bowtie2","bowtie2-build","samtools"):
        if not shutil.which(exe):
            raise SystemExit(f"ERROR {exe} not found; update/activate st10042-assembly")

    loci=scheme_loci()
    print(f"PubMLST scheme {SCHEME_ID} loci: {len(loci)}")
    if len(loci)!=EXPECTED_LOCI:
        print(f"WARNING expected {EXPECTED_LOCI}")

    c=pd.read_csv(a.combined,sep="\t",dtype=str).fillna("")
    overlaps=c[c["provenance"].eq("Azevedo+PubMLST")].copy()

    m=pd.read_csv(a.spades_manifest,sep="\t",dtype=str).fillna("")
    reads=m.set_index("Strain_ID").to_dict("index")

    rows=[]
    for _,r in overlaps.iterrows():
        aid=r["azevedo_strain_id"]
        pid=str(r["pubmlst_id"]).split(";")[0]
        source=Path(a.assembly_dir)/aid/"contigs.fasta"
        if not source.exists() or aid not in reads:
            continue

        rr=reads[aid]
        # Match the reads used for our SPAdes assembly. INNUca likewise updates
        # fastq_files to the trimmed paired reads before assembly mapping.
        r1=Path(a.trim_dir)/aid/f"{aid}_R1.fastq.gz"
        r2=Path(a.trim_dir)/aid/f"{aid}_R2.fastq.gz"
        if not r1.exists() or not r2.exists():
            raise SystemExit(
                f"ERROR missing fastp-trimmed reads for {aid}: {r1} {r2}"
            )

        print(f"\n{aid} / PubMLST {pid}")
        od=Path(a.work_dir)/aid
        od.mkdir(parents=True,exist_ok=True)

        # Stage 1: approximate INNUca's SPAdes contig filters.
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

        stage1=od/"contigs.kmercov_filtered.fasta"
        write_fasta(pre,stage1)
        print(
            f"  stage1 contigs: {len(pre)}/{len(records)} "
            f"(min_len={a.min_length}, SPAdes_cov>={a.min_spades_cov}); "
            f"headers_without_cov={no_cov_header}"
        )

        # Map processed paired reads back to the stage-1 assembly.
        prefix=str(od/"bt2")
        if not Path(prefix+".1.bt2").exists() and not Path(prefix+".1.bt2l").exists():
            run(["bowtie2-build","--quiet",str(stage1),prefix])

        bam=od/"reads.sorted.bam"
        if not bam.exists():
            sam=od/"reads.sam"
            with open(sam,"w") as sfh:
                run([
                    "bowtie2","--quiet","--very-sensitive-local","--threads",str(a.threads),
                    "--fr","-I","0","-X","2000",
                    "-x",prefix,"-1",str(r1),"-2",str(r2)
                ],stdout=sfh)
            run(["samtools","sort","-@",str(a.threads),"-o",str(bam),str(sam)])
            run(["samtools","index",str(bam)])
            sam.unlink(missing_ok=True)

        # Compute depth at every position and per-contig mean depth.
        lengths={h:len(s) for h,s in pre}
        covsum={h:0 for h,_ in pre}
        mapped_positions={h:0 for h,_ in pre}
        depth=od/"depth.tsv"
        if not depth.exists():
            with open(depth,"w") as dfh:
                run(["samtools","depth","-a",str(bam)],stdout=dfh)

        with open(depth,encoding="utf-8") as fh:
            for line in fh:
                chrom,pos,dep=line.rstrip("\n").split("\t")[:3]
                if chrom in covsum:
                    covsum[chrom] += int(dep)
                    mapped_positions[chrom] += 1

        total_positions=sum(lengths.values())
        total_depth=sum(covsum.values())
        assembly_mean=(total_depth/total_positions) if total_positions else 0.0
        threshold=max(10.0,assembly_mean/3.0)

        keep=[]
        per_contig=[]
        for h,s in pre:
            denom=lengths[h]
            mean=(covsum[h]/denom) if denom else 0.0
            per_contig.append((h,denom,spades_cov(h),mean))
            if mean >= threshold:
                keep.append((h,s))

        filtered=od/"contigs.innuca_like.fasta"
        write_fasta(keep,filtered)

        # Mapping rate as a diagnostic.
        flag=subprocess.run(
            ["samtools","flagstat",str(bam)],
            check=True,text=True,capture_output=True
        ).stdout.splitlines()
        mapline=next((x for x in flag if " mapped (" in x),"")

        print(
            f"  mapped assembly mean depth={assembly_mean:.2f}x; "
            f"filter threshold={threshold:.2f}x"
        )
        print(
            f"  kept contigs={len(keep)}/{len(pre)}; "
            f"bp={sum(len(s) for _,s in keep)}; {mapline.strip()}"
        )

        # Save per-contig coverage report.
        with open(od/"contig_coverage.tsv","w",newline="",encoding="utf-8") as fh:
            w=csv.writer(fh,delimiter="\t")
            w.writerow(["contig","length","spades_kmer_cov","mapped_mean_depth","kept"])
            keep_names={h for h,_ in keep}
            for h,L,kcov,mcov in per_contig:
                w.writerow([h,L,"" if kcov is None else kcov,f"{mcov:.4f}","Yes" if h in keep_names else "No"])

        res=scheme_query(
            filtered,
            Path(a.cache_dir)/f"{aid}.innuca_like.sequence.json"
        )
        ss=summarise_scheme(res,loci)
        print(
            f"  scheme8 exact={ss['exact_loci_any']}/{len(loci)} "
            f"unique={ss['unique_exact_loci']} ambiguous={ss['ambiguous_loci']} "
            f"no_exact={ss['loci_without_exact_hit']}"
        )

        rows.append({
            "azevedo_strain_id":aid,
            "pubmlst_id":pid,
            "published_cluster":r["published_figure4_cluster"],
            "published_cluster21":r["published_cluster21"],
            "source_contigs":len(records),
            "stage1_contigs":len(pre),
            "filtered_contigs":len(keep),
            "filtered_bp":sum(len(s) for _,s in keep),
            "assembly_mean_depth":assembly_mean,
            "coverage_threshold":threshold,
            "flagstat_mapped_line":mapline.strip(),
            **ss,
            "filtered_fasta":str(filtered),
        })

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    print("\nSUMMARY")
    if len(out):
        cols=[
            "azevedo_strain_id","pubmlst_id","published_cluster",
            "filtered_contigs","filtered_bp","assembly_mean_depth","coverage_threshold",
            "exact_loci_any","unique_exact_loci","ambiguous_loci","loci_without_exact_hit"
        ]
        print(out[cols].to_string(index=False,float_format=lambda x:f"{x:.2f}"))
    else:
        print("No overlap samples processed.")

    print(f"\nWrote: {a.out}")
    print(f"Filtered assemblies: {a.work_dir}")
    print("This is an INNUca-like filter, not a byte-for-byte reproduction of the paper workflow.")


if __name__=="__main__":
    main()
