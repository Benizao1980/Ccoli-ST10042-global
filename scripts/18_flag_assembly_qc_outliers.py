#!/usr/bin/env python3
"""Flag Azevedo assemblies that merit manual review before cgMLST-v2 typing.

These are REVIEW flags, not automatic exclusions. The primary downstream gate is
PubMLST cgMLST-v2 quality (including <=25 missing loci), but gross size/GC/
fragmentation outliers should be inspected first because mixed or contaminated
assemblies can still produce misleading allele calls.
"""

import argparse
from pathlib import Path
import pandas as pd


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--qc",default="results/azevedo_spades_qc.tsv")
    p.add_argument("--out",default="results/azevedo_spades_qc_review.tsv")
    p.add_argument("--min-bp",type=int,default=1400000)
    p.add_argument("--max-bp",type=int,default=2200000)
    p.add_argument("--min-gc",type=float,default=29.0)
    p.add_argument("--max-gc",type=float,default=34.0)
    p.add_argument("--max-contigs",type=int,default=200)
    p.add_argument("--min-n50",type=int,default=30000)
    a=p.parse_args()

    x=pd.read_csv(a.qc,sep="\t")
    numcols=["total_bp_ge500","n_contigs_ge500","N50_ge500","largest_contig","GC_percent_all"]
    for c in numcols:
        x[c]=pd.to_numeric(x[c],errors="coerce")

    def flags(r):
        f=[]
        if not bool(r["assembly_present"]):
            f.append("missing_assembly")
            return ";".join(f)
        if r["total_bp_ge500"] < a.min_bp:
            f.append("small_genome")
        if r["total_bp_ge500"] > a.max_bp:
            f.append("large_genome")
        if r["GC_percent_all"] < a.min_gc or r["GC_percent_all"] > a.max_gc:
            f.append("GC_outlier")
        if r["n_contigs_ge500"] > a.max_contigs:
            f.append("fragmented")
        if r["N50_ge500"] < a.min_n50:
            f.append("low_N50")
        return ";".join(f)

    x["review_flags"]=x.apply(flags,axis=1)
    x["review_required"]=x["review_flags"].ne("")

    # Robust cohort context without making robust-z a pass/fail rule.
    for c in ["total_bp_ge500","n_contigs_ge500","N50_ge500","GC_percent_all"]:
        med=x[c].median()
        mad=(x[c]-med).abs().median()
        x[f"{c}_median"]=med
        x[f"{c}_MAD"]=mad
        x[f"{c}_robust_z"]=(0.6745*(x[c]-med)/mad) if mad and mad>0 else 0.0

    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    x.to_csv(a.out,sep="\t",index=False)

    review=x[x["review_required"]].copy()
    print(f"Assemblies assessed: {len(x)}")
    print(f"Assemblies flagged for review: {len(review)}")
    print(f"Assemblies without review flags: {len(x)-len(review)}")
    print()
    if len(review):
        cols=[
            "array_index","Strain_ID","Country","run_accession","Figure4_cluster",
            "total_bp_ge500","n_contigs_ge500","N50_ge500","largest_contig",
            "GC_percent_all","review_flags"
        ]
        print("FLAGGED ASSEMBLIES")
        print(review[cols].sort_values(
            ["total_bp_ge500","n_contigs_ge500"],ascending=[False,False]
        ).to_string(index=False))
    else:
        print("No assemblies crossed the review thresholds.")

    print("\nMOST EXTREME BY ASSEMBLY SIZE")
    print(x.sort_values("total_bp_ge500",ascending=False)[
        ["array_index","Strain_ID","run_accession","total_bp_ge500",
         "n_contigs_ge500","N50_ge500","GC_percent_all","review_flags"]
    ].head(10).to_string(index=False))

    print("\nMOST FRAGMENTED")
    print(x.sort_values("n_contigs_ge500",ascending=False)[
        ["array_index","Strain_ID","run_accession","total_bp_ge500",
         "n_contigs_ge500","N50_ge500","GC_percent_all","review_flags"]
    ].head(10).to_string(index=False))

    print(f"\nWrote: {a.out}")
    print("NOTE: review flags are diagnostic, not exclusion criteria.")


if __name__=="__main__":
    main()
