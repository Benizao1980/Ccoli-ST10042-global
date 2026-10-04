#!/usr/bin/env python3
"""Build the currently sequence-accessible ST10042 analysis manifest.

Inputs:
- provisional 249-isolate combined manifest;
- ENA resolution inventory for Azevedo runs;
- local PubMLST assembly directory.

Rules:
- Azevedo records are sequence-accessible only when ENA currently exposes FASTQ;
- PubMLST-only records are sequence-accessible when their downloaded FASTA exists;
- exact Azevedo+PubMLST overlaps remain represented once by the Azevedo record.

This produces a current-access snapshot, not a biological/QC inclusion decision.
"""

import argparse
from pathlib import Path
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--combined", default="results/combined_st10042_provisional.tsv")
    p.add_argument("--ena", default="data/europe_ena_fastq_manifest.tsv")
    p.add_argument("--pubmlst-dir", default="data/pubmlst_st10042_assemblies")
    p.add_argument("--out", default="results/combined_st10042_currently_accessible.tsv")
    a = p.parse_args()

    x = pd.read_csv(a.combined, sep="\t", dtype=str).fillna("")
    en = pd.read_csv(a.ena, sep="\t", dtype=str).fillna("")

    ena_status = en.set_index("Strain_ID")["status"].to_dict()
    ena_fastq = en.set_index("Strain_ID")["fastq_ftp"].to_dict()

    availability = []
    sources = []
    locations = []

    for _, r in x.iterrows():
        provenance = r["provenance"]
        aid = r["azevedo_strain_id"]
        pid = r["pubmlst_id"]

        if provenance in ("Azevedo", "Azevedo+PubMLST"):
            ok = ena_status.get(aid, "") == "resolved"
            availability.append(ok)
            sources.append("ENA_FASTQ" if ok else "currently_unavailable")
            locations.append(ena_fastq.get(aid, "") if ok else "")
        elif provenance == "PubMLST":
            fasta = Path(a.pubmlst_dir) / f"{pid}.fasta"
            ok = fasta.exists()
            availability.append(ok)
            sources.append("PubMLST_contigs" if ok else "missing_local_PubMLST_contigs")
            locations.append(str(fasta) if ok else "")
        else:
            availability.append(False)
            sources.append("unknown")
            locations.append("")

    x["sequence_accessible_now"] = availability
    x["sequence_source_now"] = sources
    x["sequence_location_now"] = locations

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    x.to_csv(a.out, sep="\t", index=False)

    accessible = x[x["sequence_accessible_now"]]
    unavailable = x[~x["sequence_accessible_now"]]

    print(f"Provisional unique candidate isolates: {len(x)}")
    print(f"Sequence-accessible now: {len(accessible)}")
    print(f"Currently unavailable: {len(unavailable)}")
    print(f"Accessible Azevedo backbone: {(accessible['provenance'].isin(['Azevedo','Azevedo+PubMLST'])).sum()}")
    print(f"Accessible PubMLST-only additions: {(accessible['provenance'].eq('PubMLST')).sum()}")
    print(f"Accessible Peru focal isolates: {((accessible['peru_focal'].eq('Yes'))).sum()}")
    print()
    print("ACCESSIBLE PUBLISHED CLUSTER 21")
    c21 = accessible[accessible["published_cluster21"].eq("Yes")]
    print(f"n = {len(c21)}")
    if len(c21):
        print(c21["country"].value_counts().to_string())
    print()
    print("CURRENTLY UNAVAILABLE BY COUNTRY")
    print(unavailable["country"].value_counts().to_string())
    print()
    print(f"Wrote: {a.out}")


if __name__ == "__main__":
    main()
