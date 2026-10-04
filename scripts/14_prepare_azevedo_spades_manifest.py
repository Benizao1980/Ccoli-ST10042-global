#!/usr/bin/env python3
"""Create a 0-based Slurm manifest for the 110 paired Azevedo FASTQ runs."""

import argparse
import re
from pathlib import Path
import pandas as pd


def split(v):
    return [x.strip() for x in str(v).split(";") if x.strip()]


def mate_number(url):
    name=Path(url).name
    pats=[
        (r"(?:^|[_\.])R?1(?:[_\.]|\.f(?:ast)?q)",1),
        (r"(?:^|[_\.])R?2(?:[_\.]|\.f(?:ast)?q)",2),
    ]
    for pat,mate in pats:
        if re.search(pat,name,re.I):
            return mate
    return None


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--plan",default="results/azevedo_accessible_fastq_plan.tsv")
    p.add_argument("--out",default="results/azevedo_spades_manifest.tsv")
    p.add_argument("--read-dir",default="data/europe_reads")
    a=p.parse_args()

    df=pd.read_csv(a.plan,sep="\t",dtype=str).fillna("")
    df=df[df["n_fastq_files"].eq("2")].copy()

    rows=[]
    ambiguous=0
    for _,r in df.iterrows():
        urls=split(r["fastq_ftp"])
        md5s=split(r["fastq_md5"])
        sizes=split(r["fastq_bytes"])
        if not (len(urls)==len(md5s)==len(sizes)==2):
            raise SystemExit(f"ERROR inconsistent paired metadata for {r['Strain_ID']}")

        mates=[mate_number(x) for x in urls]
        if mates==[1,2]:
            order=[0,1]; method="filename"
        elif mates==[2,1]:
            order=[1,0]; method="filename"
        else:
            order=[0,1]; method="ENA_order"
            ambiguous += 1

        u1,u2=urls[order[0]],urls[order[1]]
        m1,m2=md5s[order[0]],md5s[order[1]]
        b1,b2=sizes[order[0]],sizes[order[1]]
        sample=r["Strain_ID"]
        sd=Path(a.read_dir)/sample

        rows.append({
            "array_index":len(rows),
            "Strain_ID":sample,
            "Country":r["Country"],
            "run_accession":r["run_accession"],
            "Figure4_cluster":r["Figure4_cluster"],
            "read1_url":u1,
            "read2_url":u2,
            "read1_md5":m1,
            "read2_md5":m2,
            "read1_bytes":b1,
            "read2_bytes":b2,
            "read1_path":str(sd/Path(u1).name),
            "read2_path":str(sd/Path(u2).name),
            "pairing_method":method,
        })

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    print(f"Paired Azevedo runs in SPAdes manifest: {len(out)}")
    print(f"Filename-explicit R1/R2 pairing: {(out['pairing_method']=='filename').sum()}")
    print(f"Pairs using ENA file order fallback: {ambiguous}")
    if ambiguous:
        print(out[out["pairing_method"].eq("ENA_order")][
            ["Strain_ID","run_accession","read1_url","read2_url"]
        ].to_string(index=False))
    print(f"Wrote: {a.out}")


if __name__=="__main__":
    main()
