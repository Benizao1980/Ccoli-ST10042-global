#!/usr/bin/env python3
"""Prepare the known Azevedo/PubMLST overlaps for a chewBBACA allele-call pilot."""

import argparse
from pathlib import Path
import pandas as pd


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--combined",default="results/combined_st10042_provisional.tsv")
    p.add_argument("--cgmlst-manifest",default="results/cgmlst_v2_input_manifest.tsv")
    p.add_argument("--out-list",default="results/chewbbaca_overlap_inputs.txt")
    p.add_argument("--out-meta",default="results/chewbbaca_overlap_inputs.tsv")
    a=p.parse_args()

    combined=pd.read_csv(a.combined,sep="\t",dtype=str).fillna("")
    manifest=pd.read_csv(a.cgmlst_manifest,sep="\t",dtype=str).fillna("")

    ov=combined[combined["provenance"].eq("Azevedo+PubMLST")].copy()
    keep=set("AZE_"+x for x in ov["azevedo_strain_id"] if x)
    q=manifest[manifest["analysis_id"].isin(keep)].copy()

    if len(q)!=9:
        raise SystemExit(f"ERROR expected 9 known overlap controls, found {len(q)}")

    missing=[]
    paths=[]
    for _,r in q.iterrows():
        f=Path(r["chewbbaca_fasta"])
        if not f.exists() or f.stat().st_size==0:
            missing.append((r["analysis_id"],str(f)))
        paths.append(str(f.absolute()))
    if missing:
        raise SystemExit("ERROR missing pilot FASTAs:\n" + "\n".join(f"{x}\t{p}" for x,p in missing))

    q=q.sort_values("analysis_id")
    Path(a.out_meta).parent.mkdir(parents=True,exist_ok=True)
    q.to_csv(a.out_meta,sep="\t",index=False)
    Path(a.out_list).write_text(
        "\n".join(str(Path(x).absolute()) for x in q["chewbbaca_fasta"])+"\n",
        encoding="utf-8",
    )

    print(f"Overlap controls: {len(q)}")
    print(q[[
        "analysis_id","azevedo_strain_id","pubmlst_id",
        "published_figure4_cluster","published_cluster21"
    ]].to_string(index=False))
    print(f"Wrote: {a.out_list}")
    print(f"Wrote: {a.out_meta}")


if __name__=="__main__":
    main()
