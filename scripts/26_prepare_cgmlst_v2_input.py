#!/usr/bin/env python3
"""Build the 142-genome cgMLST-v2 input manifest/list.

Uses:
- 110 coverage-filtered Azevedo assemblies from the assembly-ready manifest;
- 32 exact-ID-unmatched PubMLST assemblies from the provisional combined manifest.

No genome is excluded here on the basis of the exact-hit QC. The QC columns are
carried forward so PT-41/PT-58 and the known poor UK genomes can be adjudicated
after full allele calling.
"""

import argparse
from pathlib import Path
import pandas as pd

EXPECTED_AZE=110
EXPECTED_PUB=32
EXPECTED_TOTAL=142


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--combined",default="results/combined_st10042_provisional.tsv")
    p.add_argument("--spades-manifest",default="results/azevedo_spades_manifest.tsv")
    p.add_argument("--aze-qc",default="results/pubmlst_v2_all_azevedo_exact_qc.tsv")
    p.add_argument("--aze-dir",default="data/europe_assemblies_innuca_like_trimmed")
    p.add_argument("--pub-dir",default="data/pubmlst_st10042_assemblies")
    p.add_argument("--out",default="results/cgmlst_v2_input_manifest.tsv")
    p.add_argument("--list-out",default="results/cgmlst_v2_input_fastas.txt")
    a=p.parse_args()

    combined=pd.read_csv(a.combined,sep="\t",dtype=str).fillna("")
    spades=pd.read_csv(a.spades_manifest,sep="\t",dtype=str).fillna("")
    qc=pd.read_csv(a.aze_qc,sep="\t",dtype=str).fillna("")
    qc=qc.set_index("Strain_ID")

    combined_aze=combined.set_index("azevedo_strain_id")
    rows=[]

    for _,r in spades.iterrows():
        sid=r["Strain_ID"]
        fasta=Path(a.aze_dir)/sid/"contigs.innuca_like.fasta"
        meta=combined_aze.loc[sid] if sid in combined_aze.index else pd.Series(dtype=str)
        q=qc.loc[sid] if sid in qc.index else pd.Series(dtype=str)
        rows.append({
            "analysis_id":f"AZE_{sid}",
            "provenance":"Azevedo",
            "azevedo_strain_id":sid,
            "pubmlst_id":meta.get("pubmlst_id",""),
            "country":meta.get("country",r.get("Country","")),
            "town_or_city":meta.get("town_or_city",""),
            "year":meta.get("year",""),
            "source":meta.get("source",""),
            "peru_focal":"No",
            "published_figure4_cluster":meta.get("published_figure4_cluster",""),
            "published_cluster21":meta.get("published_cluster21",""),
            "prelim_exact_loci":q.get("exact_loci_any",""),
            "prelim_ambiguous_loci":q.get("ambiguous_loci",""),
            "prelim_loci_without_exact_hit":q.get("loci_without_exact_hit",""),
            "fasta":str(fasta.resolve()),
        })

    pub=combined[combined["provenance"].eq("PubMLST")].copy()
    for _,r in pub.iterrows():
        pid=str(r["pubmlst_id"]).strip()
        fasta=Path(a.pub_dir)/f"{pid}.fasta"
        rows.append({
            "analysis_id":f"PUB_{pid}",
            "provenance":"PubMLST",
            "azevedo_strain_id":"",
            "pubmlst_id":pid,
            "country":r.get("country",""),
            "town_or_city":r.get("town_or_city",""),
            "year":r.get("year",""),
            "source":r.get("source",""),
            "peru_focal":r.get("peru_focal",""),
            "published_figure4_cluster":"",
            "published_cluster21":"",
            "prelim_exact_loci":"",
            "prelim_ambiguous_loci":"",
            "prelim_loci_without_exact_hit":"",
            "fasta":str(fasta.resolve()),
        })

    out=pd.DataFrame(rows)
    missing=[]
    for _,r in out.iterrows():
        f=Path(r["fasta"])
        if not f.exists() or f.stat().st_size==0:
            missing.append((r["analysis_id"],str(f)))

    n_aze=(out["provenance"]=="Azevedo").sum()
    n_pub=(out["provenance"]=="PubMLST").sum()
    if n_aze!=EXPECTED_AZE or n_pub!=EXPECTED_PUB or len(out)!=EXPECTED_TOTAL:
        raise SystemExit(
            f"ERROR unexpected counts: Azevedo={n_aze} PubMLST={n_pub} total={len(out)}"
        )
    if missing:
        raise SystemExit("ERROR missing FASTAs:\n" + "\n".join(f"{x}\t{p}" for x,p in missing))

    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)
    Path(a.list_out).write_text("\n".join(out["fasta"])+"\n",encoding="utf-8")

    print(f"Azevedo filtered assemblies: {n_aze}")
    print(f"PubMLST-only assemblies: {n_pub}")
    print(f"Total cgMLST input genomes: {len(out)}")
    print(f"Peru focal genomes: {(out['peru_focal']=='Yes').sum()}")
    print(f"Published accessible cluster-21 anchors: {(out['published_cluster21']=='Yes').sum()}")
    print(f"Wrote: {a.out}")
    print(f"Wrote: {a.list_out}")


if __name__=="__main__":
    main()
