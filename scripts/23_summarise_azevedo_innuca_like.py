#!/usr/bin/env python3
"""Summarise completion/QC of the full Azevedo INNUca-like filtered set."""

import argparse
import json
from pathlib import Path

import pandas as pd


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",default="results/azevedo_spades_manifest.tsv")
    p.add_argument("--filtered-dir",default="data/europe_assemblies_innuca_like_trimmed")
    p.add_argument("--out",default="results/azevedo_innuca_like_trimmed_qc.tsv")
    a=p.parse_args()

    m=pd.read_csv(a.manifest,sep="\t",dtype=str).fillna("")
    rows=[]
    missing=[]

    for _,r in m.iterrows():
        sample=r["Strain_ID"]
        q=Path(a.filtered_dir)/sample/"qc.json"
        fa=Path(a.filtered_dir)/sample/"contigs.innuca_like.fasta"
        if not q.exists() or not fa.exists() or fa.stat().st_size==0:
            missing.append({
                "array_index":r["array_index"],
                "Strain_ID":sample,
                "run_accession":r.get("run_accession",""),
            })
            continue
        x=json.loads(q.read_text(encoding="utf-8"))
        x["filtered_fasta_present"]=True
        rows.append(x)

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    print(f"Expected Azevedo genomes: {len(m)}")
    print(f"Complete filtered genomes: {len(out)}")
    print(f"Missing/incomplete: {len(missing)}")

    if len(out):
        for col in [
            "source_contigs","source_bp","stage1_contigs","stage1_bp",
            "filtered_contigs","filtered_bp","assembly_mean_depth","coverage_threshold"
        ]:
            x=pd.to_numeric(out[col],errors="coerce")
            print(
                f"{col}: median={x.median():.2f} "
                f"min={x.min():.2f} max={x.max():.2f}"
            )

        print("\nLargest filtered assemblies")
        cols=["array_index","Strain_ID","Country","Figure4_cluster",
              "filtered_contigs","filtered_bp","assembly_mean_depth","coverage_threshold"]
        print(
            out.sort_values("filtered_bp",ascending=False)[cols]
               .head(12).to_string(index=False)
        )

        print("\nMost fragmented filtered assemblies")
        print(
            out.sort_values("filtered_contigs",ascending=False)[cols]
               .head(12).to_string(index=False)
        )

    if missing:
        print("\nMISSING/INCOMPLETE")
        print(pd.DataFrame(missing).to_string(index=False))

    print(f"\nWrote: {a.out}")


if __name__=="__main__":
    main()
