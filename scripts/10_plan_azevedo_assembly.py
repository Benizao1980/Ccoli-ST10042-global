#!/usr/bin/env python3
"""Plan the 111 currently accessible Azevedo FASTQ downloads before assembly."""

import argparse
from pathlib import Path
import pandas as pd


def split_field(value):
    if pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--ena", default="data/europe_ena_fastq_manifest.tsv")
    p.add_argument("--out", default="results/azevedo_accessible_fastq_plan.tsv")
    a=p.parse_args()

    df=pd.read_csv(a.ena, sep="\t", dtype=str).fillna("")
    df=df[df["status"].eq("resolved")].copy()

    rows=[]
    for _,r in df.iterrows():
        urls=split_field(r["fastq_ftp"])
        md5s=split_field(r["fastq_md5"])
        sizes=split_field(r["fastq_bytes"])
        sizes_int=[]
        for x in sizes:
            try:
                sizes_int.append(int(x))
            except ValueError:
                sizes_int.append(0)

        rows.append({
            "Strain_ID":r["Strain_ID"],
            "Country":r["Country"],
            "requested_run_accession":r["requested_run_accession"],
            "run_accession":r["run_accession"],
            "Figure4_cluster":r["Figure4_cluster"],
            "n_fastq_files":len(urls),
            "total_fastq_bytes":sum(sizes_int),
            "total_fastq_GiB":sum(sizes_int)/(1024**3),
            "fastq_ftp":r["fastq_ftp"],
            "fastq_md5":r["fastq_md5"],
            "fastq_bytes":r["fastq_bytes"],
            "metadata_length_consistent":len(urls)==len(md5s)==len(sizes_int),
        })

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    total_bytes=int(out["total_fastq_bytes"].sum())
    print(f"Accessible Azevedo runs: {len(out)}")
    print(f"FASTQ files: {int(out['n_fastq_files'].sum())}")
    print(f"Runs with 2 FASTQ files: {(out['n_fastq_files']==2).sum()}")
    print(f"Runs with 1 FASTQ file: {(out['n_fastq_files']==1).sum()}")
    print(f"Runs with other FASTQ file counts: {(~out['n_fastq_files'].isin([1,2])).sum()}")
    print(f"ENA metadata length-consistent: {out['metadata_length_consistent'].sum()}/{len(out)}")
    print(f"Compressed download size: {total_bytes/1e9:.2f} GB ({total_bytes/(1024**3):.2f} GiB)")

    print("\nBY COUNTRY")
    country=out.groupby("Country").agg(
        runs=("run_accession","size"),
        GiB=("total_fastq_GiB","sum"),
    ).sort_values("runs",ascending=False)
    print(country.to_string(float_format=lambda x:f"{x:.2f}"))

    print("\n10 LARGEST RUNS")
    cols=["Strain_ID","Country","run_accession","n_fastq_files","total_fastq_GiB"]
    print(
        out.sort_values("total_fastq_bytes",ascending=False)[cols]
        .head(10)
        .to_string(index=False,float_format=lambda x:f"{x:.2f}")
    )

    print(f"\nWrote: {a.out}")


if __name__=="__main__":
    main()
