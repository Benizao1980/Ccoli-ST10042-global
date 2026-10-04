#!/usr/bin/env python3
"""Conservatively trim and assemble one Azevedo Illumina pair with SPAdes --isolate."""

import argparse
import csv
import json
import shutil
import subprocess
from pathlib import Path


def get_row(path,index):
    with open(path,encoding="utf-8") as fh:
        rows=list(csv.DictReader(fh,delimiter="\t"))
    if index<0 or index>=len(rows):
        raise SystemExit(f"ERROR index {index} outside manifest range 0..{len(rows)-1}")
    return rows[index]


def run(cmd):
    print("+"," ".join(str(x) for x in cmd),flush=True)
    subprocess.run(cmd,check=True)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",default="results/azevedo_spades_manifest.tsv")
    p.add_argument("--index",type=int,required=True)
    p.add_argument("--threads",type=int,default=4)
    p.add_argument("--memory-gb",type=int,default=15)
    p.add_argument("--assembly-dir",default="data/europe_assemblies")
    p.add_argument("--trim-dir",default="data/europe_trimmed")
    p.add_argument("--qc-dir",default="results/fastp")
    a=p.parse_args()

    r=get_row(a.manifest,a.index)
    sample=r["Strain_ID"]
    r1=Path(r["read1_path"])
    r2=Path(r["read2_path"])
    if not r1.exists() or not r2.exists():
        raise SystemExit(f"ERROR missing reads for {sample}: {r1} {r2}")

    fastp=shutil.which("fastp")
    spades=shutil.which("spades.py") or shutil.which("spades")
    if not fastp:
        raise SystemExit("ERROR fastp not found in PATH")
    if not spades:
        raise SystemExit("ERROR spades.py/spades not found in PATH")

    tdir=Path(a.trim_dir)/sample
    qdir=Path(a.qc_dir)
    odir=Path(a.assembly_dir)/sample
    tdir.mkdir(parents=True,exist_ok=True)
    qdir.mkdir(parents=True,exist_ok=True)
    odir.parent.mkdir(parents=True,exist_ok=True)

    tr1=tdir/f"{sample}_R1.fastq.gz"
    tr2=tdir/f"{sample}_R2.fastq.gz"
    fastp_json=qdir/f"{sample}.json"
    fastp_html=qdir/f"{sample}.html"

    if not tr1.exists() or not tr2.exists():
        run([
            fastp,
            "--in1",str(r1),"--in2",str(r2),
            "--out1",str(tr1),"--out2",str(tr2),
            "--detect_adapter_for_pe",
            "--thread",str(a.threads),
            "--json",str(fastp_json),
            "--html",str(fastp_html),
        ])
    else:
        print(f"Using existing trimmed reads for {sample}")

    contigs=odir/"contigs.fasta"
    if contigs.exists() and contigs.stat().st_size>0:
        print(f"Assembly already complete: {contigs}")
        return

    if odir.exists():
        # SPAdes refuses a partially populated output directory without --continue,
        # but --continue can preserve a failed command state. Remove only the
        # generated SPAdes directory; raw/trimmed reads live elsewhere.
        shutil.rmtree(odir)

    run([
        spades,
        "--isolate",
        "-1",str(tr1),"-2",str(tr2),
        "-o",str(odir),
        "-t",str(a.threads),
        "-m",str(a.memory_gb),
    ])

    if not contigs.exists() or contigs.stat().st_size==0:
        raise SystemExit(f"ERROR SPAdes finished without contigs: {sample}")
    print(f"OK assembly: {contigs}")


if __name__=="__main__":
    main()
