#!/usr/bin/env python3
"""Summarise completion and basic assembly statistics for Azevedo SPAdes outputs."""

import argparse
import csv
from pathlib import Path

import pandas as pd


def fasta_lengths(path):
    lengths=[]
    gc=0
    bases=0
    seq=[]
    with open(path,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            line=line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    s="".join(seq).upper()
                    lengths.append(len(s))
                    gc += s.count("G")+s.count("C")
                    bases += len(s)
                    seq=[]
            else:
                seq.append(line)
        if seq:
            s="".join(seq).upper()
            lengths.append(len(s))
            gc += s.count("G")+s.count("C")
            bases += len(s)
    return lengths,gc,bases


def n50(lengths):
    if not lengths:
        return 0
    xs=sorted(lengths,reverse=True)
    half=sum(xs)/2
    acc=0
    for x in xs:
        acc += x
        if acc >= half:
            return x
    return 0


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",default="results/azevedo_spades_manifest.tsv")
    p.add_argument("--assembly-dir",default="data/europe_assemblies")
    p.add_argument("--fastp-dir",default="results/fastp")
    p.add_argument("--out",default="results/azevedo_spades_qc.tsv")
    a=p.parse_args()

    m=pd.read_csv(a.manifest,sep="\t",dtype=str).fillna("")
    rows=[]
    for _,r in m.iterrows():
        sample=r["Strain_ID"]
        contigs=Path(a.assembly_dir)/sample/"contigs.fasta"
        fastp=Path(a.fastp_dir)/f"{sample}.json"

        row={
            "array_index":r["array_index"],
            "Strain_ID":sample,
            "Country":r["Country"],
            "run_accession":r["run_accession"],
            "Figure4_cluster":r["Figure4_cluster"],
            "assembly_present":contigs.exists() and contigs.stat().st_size>0,
            "contigs_path":str(contigs),
            "fastp_json_present":fastp.exists() and fastp.stat().st_size>0,
        }

        if row["assembly_present"]:
            lengths,gc,bases=fasta_lengths(contigs)
            ge500=[x for x in lengths if x>=500]
            # GC is reported over all assembled bases; for Campylobacter this is
            # mainly a sanity descriptor rather than a pass/fail criterion.
            row.update({
                "n_contigs_all":len(lengths),
                "total_bp_all":sum(lengths),
                "n_contigs_ge500":len(ge500),
                "total_bp_ge500":sum(ge500),
                "N50_ge500":n50(ge500),
                "largest_contig":max(lengths) if lengths else 0,
                "GC_percent_all":(100*gc/bases) if bases else "",
            })
        else:
            row.update({
                "n_contigs_all":"",
                "total_bp_all":"",
                "n_contigs_ge500":"",
                "total_bp_ge500":"",
                "N50_ge500":"",
                "largest_contig":"",
                "GC_percent_all":"",
            })
        rows.append(row)

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    complete=out[out["assembly_present"]]
    missing=out[~out["assembly_present"]]

    print(f"Expected paired Azevedo assemblies: {len(out)}")
    print(f"Assemblies with non-empty contigs.fasta: {len(complete)}")
    print(f"Missing/failed assemblies: {len(missing)}")
    print(f"fastp JSON reports present: {out['fastp_json_present'].sum()}/{len(out)}")

    if len(complete):
        for col in ["total_bp_ge500","n_contigs_ge500","N50_ge500","largest_contig","GC_percent_all"]:
            x=pd.to_numeric(complete[col],errors="coerce")
            print(f"{col}: median={x.median():.1f} min={x.min():.1f} max={x.max():.1f}")

    if len(missing):
        print("\nMISSING/FAILED")
        print(missing[["array_index","Strain_ID","run_accession"]].to_string(index=False))

    print(f"\nWrote: {a.out}")


if __name__=="__main__":
    main()
