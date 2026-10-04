#!/usr/bin/env python3
"""Package the 110 filtered Azevedo genomes for PubMLST Genome Comparator.

Creates a real tar.gz archive (not symlinks) with unique FASTA names plus a
metadata table. Also writes the current ST10042 PubMLST isolate IDs that have
local contig FASTAs so they can be selected/checked as database anchors.
"""

import argparse
import tarfile
from pathlib import Path

import pandas as pd


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",default="results/cgmlst_v2_input_manifest.tsv")
    p.add_argument("--pubmlst",default="results/pubmlst_st10042_v2_lincodes.tsv")
    p.add_argument("--pubmlst-dir",default="data/pubmlst_st10042_assemblies")
    p.add_argument("--out-dir",default="results/pubmlst_genome_comparator")
    a=p.parse_args()

    outdir=Path(a.out_dir)
    outdir.mkdir(parents=True,exist_ok=True)

    m=pd.read_csv(a.manifest,sep="\t",dtype=str).fillna("")
    az=m[m["provenance"].eq("Azevedo")].copy()
    if len(az)!=110:
        raise SystemExit(f"ERROR expected 110 Azevedo genomes, found {len(az)}")

    archive=outdir/"azevedo_110_filtered_fastas.tar.gz"
    meta_rows=[]
    with tarfile.open(archive,"w:gz") as tf:
        for _,r in az.sort_values("analysis_id").iterrows():
            src=Path(r["fasta"])
            if not src.exists() or src.stat().st_size==0:
                raise SystemExit(f"ERROR missing FASTA: {src}")
            name=f"{r['analysis_id']}.fasta"
            tf.add(src,arcname=name,recursive=False)
            meta_rows.append({
                "filename":name,
                "analysis_id":r["analysis_id"],
                "azevedo_strain_id":r["azevedo_strain_id"],
                "country":r["country"],
                "year":r["year"],
                "source":r["source"],
                "published_figure4_cluster":r["published_figure4_cluster"],
                "published_cluster21":r["published_cluster21"],
            })

    meta=pd.DataFrame(meta_rows)
    meta_path=outdir/"azevedo_110_metadata.tsv"
    meta.to_csv(meta_path,sep="\t",index=False)

    pm=pd.read_csv(a.pubmlst,sep="\t",dtype=str).fillna("")
    anchor_rows=[]
    for _,r in pm.iterrows():
        pid=str(r["pubmlst_id"]).strip()
        fa=Path(a.pubmlst_dir)/f"{pid}.fasta"
        if fa.exists() and fa.stat().st_size>0:
            anchor_rows.append({
                "pubmlst_id":pid,
                "isolate":r.get("isolate",""),
                "country":r.get("country",""),
                "town_or_city":r.get("town_or_city",""),
                "year":r.get("year",""),
                "source":r.get("source",""),
                "cgST_v2":r.get("cgST_v2",""),
                "LINcode_v2":r.get("LINcode_v2",""),
                "Cjc_cgc2_5":r.get("Cjc_cgc2_5",""),
            })

    anchors=pd.DataFrame(anchor_rows)
    anchors_path=outdir/"current_pubmlst_st10042_anchors.tsv"
    anchors.to_csv(anchors_path,sep="\t",index=False)

    ids_path=outdir/"current_pubmlst_st10042_anchor_ids.txt"
    ids_path.write_text("\n".join(anchors["pubmlst_id"])+"\n",encoding="utf-8")

    print(f"Azevedo genomes packaged: {len(meta)}")
    print(f"Archive: {archive}")
    print(f"Archive size: {archive.stat().st_size/1024/1024:.1f} MB")
    print(f"Metadata: {meta_path}")
    print(f"Current PubMLST ST10042 anchors with contigs: {len(anchors)}")
    print(f"Anchor metadata: {anchors_path}")
    print(f"Anchor IDs: {ids_path}")
    print()
    print("Use the archive as uploaded genomes in PubMLST Genome Comparator with cgMLST v2 (scheme 8).")
    print("Use current database ST10042 isolates as reference anchors where the interface permits.")
    print("Uploaded genomes do not automatically acquire official cgST/LINcodes if they contain novel profiles.")


if __name__=="__main__":
    main()
