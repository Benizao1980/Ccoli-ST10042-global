#!/usr/bin/env python3
"""Join LINwalker placement results to ST10042/Azevedo metadata and summarise biology.

This is a project-specific post-processing layer. LINwalker provides conservative
reference-anchored placement; this script adds published Azevedo cluster labels and
current PubMLST/Peru anchor context.

It does not convert a supported prefix into an official LINcode.
"""

import argparse
import re
from pathlib import Path

import pandas as pd


def split_ids(value):
    return [x for x in str(value).split(";") if x and x.lower() not in {"nan", "none"}]


def first_nonempty(values):
    vals=[str(x) for x in values if str(x).strip() and str(x).lower()!="nan"]
    return ";".join(sorted(set(vals)))


def extract_analysis_id(value):
    s=str(value)
    m=re.search(r"(AZE_[A-Za-z0-9_.-]+)",s)
    if not m:
        return s
    x=m.group(1)
    for suffix in (".fasta.gz",".fa.gz",".fna.gz",".fasta",".fna",".fa"):
        if x.endswith(suffix):
            x=x[:-len(suffix)]
            break
    return x


def main():
    p=argparse.ArgumentParser()
    p.add_argument(
        "--placement",
        default="results/linwalker_st10042/tables/placement_summary.tsv",
    )
    p.add_argument(
        "--azevedo-metadata",
        default="results/pubmlst_genome_comparator/azevedo_110_metadata.tsv",
    )
    p.add_argument(
        "--anchors",
        default="results/pubmlst_genome_comparator/current_pubmlst_st10042_anchors.tsv",
    )
    p.add_argument("--out",default="results/st10042_v2_lin_placement.tsv")
    a=p.parse_args()

    place=pd.read_csv(a.placement,sep="\t",dtype=str).fillna("")
    az=pd.read_csv(a.azevedo_metadata,sep="\t",dtype=str).fillna("")
    anchors=pd.read_csv(a.anchors,sep="\t",dtype=str).fillna("")

    if "query_analysis_id" not in place.columns:
        place["query_analysis_id"]=place["query_id"].map(extract_analysis_id)
    else:
        place["query_analysis_id"]=place["query_analysis_id"].map(extract_analysis_id)

    if len(place)!=110:
        print(f"WARNING: expected 110 Azevedo placement rows, found {len(place)}")

    out=az.merge(
        place,
        left_on="analysis_id",
        right_on="query_analysis_id",
        how="left",
        validate="one_to_one",
    )

    missing=out[out["query_id"].fillna("").eq("")]
    if len(missing):
        print(f"WARNING: {len(missing)} Azevedo genomes have no LINwalker placement row")
        print("; ".join(missing["analysis_id"].head(20)))

    amap=anchors.set_index("pubmlst_id").to_dict("index")

    nearest_countries=[]
    nearest_sources=[]
    nearest_peru_n=[]
    nearest_c21_n=[]
    nearest_azevedo_ids=[]

    for _,r in out.iterrows():
        ids=split_ids(r.get("nearest_reference_id",""))
        meta=[amap[x] for x in ids if x in amap]
        nearest_countries.append(first_nonempty(x.get("country","") for x in meta))
        nearest_sources.append(first_nonempty(x.get("source","") for x in meta))
        nearest_peru_n.append(sum(x.get("peru_focal","")=="Yes" for x in meta))
        nearest_c21_n.append(sum(x.get("published_cluster21","")=="Yes" for x in meta))
        nearest_azevedo_ids.append(
            first_nonempty(x.get("azevedo_strain_id","") for x in meta)
        )

    out["nearest_reference_country"]=nearest_countries
    out["nearest_reference_source"]=nearest_sources
    out["nearest_reference_peru_n"]=nearest_peru_n
    out["nearest_reference_published_cluster21_n"]=nearest_c21_n
    out["nearest_reference_azevedo_strain_id"]=nearest_azevedo_ids

    # Ask whether supported fixed-threshold placements are represented in Peru and/or
    # among the current published cluster-21 PubMLST anchors.
    for threshold in (200,100,50,25,10,5):
        pcol=f"Cjc_cgc2_{threshold}_placement"
        if pcol not in out.columns or f"Cjc_cgc2_{threshold}" not in anchors.columns:
            continue

        peru_groups={}
        c21_groups={}
        for _,x in anchors.iterrows():
            group=str(x.get(f"Cjc_cgc2_{threshold}","")).strip()
            if not group:
                continue
            for g in split_ids(group):
                if x.get("peru_focal","")=="Yes":
                    peru_groups[g]=peru_groups.get(g,0)+1
                if x.get("published_cluster21","")=="Yes":
                    c21_groups[g]=c21_groups.get(g,0)+1

        out[f"{pcol}_peru_anchor_n"]=[
            peru_groups.get(str(g),0) if str(g) else 0 for g in out[pcol]
        ]
        out[f"{pcol}_published_cluster21_anchor_n"]=[
            c21_groups.get(str(g),0) if str(g) else 0 for g in out[pcol]
        ]

    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    print(f"Azevedo genomes in joined table: {len(out)}")
    print(f"Placement rows present: {out['query_id'].fillna('').ne('').sum()}")

    print("\nPLACEMENT STATUS")
    print(out["placement_status"].replace("",pd.NA).value_counts(dropna=False).to_string())

    print("\nPUBLISHED CLUSTER 21: DEEPEST SUPPORTED LIN THRESHOLD")
    c21=out[out["published_cluster21"].eq("Yes")].copy()
    if len(c21):
        print(
            c21["deepest_supported_difference_threshold"]
            .replace("",pd.NA)
            .value_counts(dropna=False)
            .sort_index()
            .to_string()
        )

        if "Cjc_cgc2_5_placement" in c21.columns:
            print("\nPUBLISHED CLUSTER 21: SUPPORTED cgc2_5 GROUPS")
            print(
                c21["Cjc_cgc2_5_placement"]
                .replace("",pd.NA)
                .value_counts(dropna=False)
                .to_string()
            )

        show=[
            "azevedo_strain_id","country","year","source",
            "nearest_reference_id","nearest_reference_country",
            "nearest_reference_peru_n","nearest_reference_published_cluster21_n",
            "nearest_normalised_AD","nearest_cgST",
            "deepest_supported_difference_threshold",
            "deepest_supported_LIN_prefix",
            "Cjc_cgc2_5_placement","Cjc_cgc2_5_status",
        ]
        show=[x for x in show if x in c21.columns]
        print("\nPUBLISHED CLUSTER 21 DETAIL")
        print(c21[show].sort_values("azevedo_strain_id").to_string(index=False))

    print(f"\nWrote: {a.out}")
    print("Supported LIN prefixes/cgc2 placements are reference-anchored interpretations,")
    print("not newly minted official PubMLST cgSTs or LINcodes.")


if __name__=="__main__":
    main()
