#!/usr/bin/env python3
"""Build the provisional combined ST10042 manifest after exact-ID de-duplication.

Canonical rule:
- retain all 217 Azevedo records as the reference backbone;
- annotate any exact PubMLST overlap onto the corresponding Azevedo record;
- append PubMLST records with no exact Azevedo identifier match.

This is an exact-identifier de-duplication, not proof that all unmatched records are
different physical isolates. Sequence-/metadata-level duplicate checks remain a later QC step.
"""

import argparse
from pathlib import Path
import pandas as pd


def s(x):
    if pd.isna(x):
        return ""
    return str(x).strip()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--azevedo", default="data/europe_217_manifest.tsv")
    p.add_argument("--pubmlst", default="results/pubmlst_st10042_curated.tsv")
    p.add_argument("--overlap", default="results/azevedo_pubmlst_exact_overlap.tsv")
    p.add_argument("--out", default="results/combined_st10042_provisional.tsv")
    a = p.parse_args()

    az = pd.read_csv(a.azevedo, sep="\t", dtype=str).fillna("")
    pm = pd.read_csv(a.pubmlst, sep="\t", dtype=str).fillna("")
    ov = pd.read_csv(a.overlap, sep="\t", dtype=str).fillna("")

    # Collapse overlap evidence to one PubMLST -> Azevedo mapping.
    pairs = (
        ov[["pubmlst_id", "azevedo_strain_id"]]
        .drop_duplicates()
        .groupby("pubmlst_id")["azevedo_strain_id"]
        .agg(lambda x: ";".join(sorted(set(x))))
        .to_dict()
    )
    reverse = {}
    for pid, aid in pairs.items():
        for x in aid.split(";"):
            reverse.setdefault(x, []).append(pid)

    rows = []

    # Azevedo backbone.
    for _, r in az.iterrows():
        aid = s(r["Strain_ID"])
        pids = sorted(reverse.get(aid, []))
        rows.append({
            "canonical_id": f"AZE_{aid}",
            "provenance": "Azevedo+PubMLST" if pids else "Azevedo",
            "azevedo_strain_id": aid,
            "pubmlst_id": ";".join(pids),
            "country": s(r["Country"]),
            "town_or_city": "",
            "year": s(r["Isolation year"]),
            "source": s(r["Source"]),
            "run_accession": s(r["Run Accession Number"]),
            "bioproject": s(r["Bioproject"]),
            "published_figure4_included": s(r["Included in Figure 4"]),
            "published_figure4_cluster": s(r["Figure 4 cluster clean"]),
            "published_cluster21": s(r["Cluster 21 member"]),
            "peru_focal": "No",
            "sequence_source_preference": "Azevedo_run",
            "exact_id_dedup_status": "exact_overlap_retained_as_Azevedo" if pids else "Azevedo_only",
        })

    # PubMLST records not exactly matched to Azevedo.
    overlap_pids = set(pairs)
    added = 0
    for _, r in pm.iterrows():
        pid = s(r["pubmlst_id"])
        if pid in overlap_pids:
            continue
        city = s(r.get("town_or_city_curated", "")) or s(r.get("town_or_city", ""))
        rows.append({
            "canonical_id": f"PUB_{pid}",
            "provenance": "PubMLST",
            "azevedo_strain_id": "",
            "pubmlst_id": pid,
            "country": s(r.get("country", "")),
            "town_or_city": city,
            "year": s(r.get("year", "")),
            "source": s(r.get("source", "")),
            "run_accession": "",
            "bioproject": "",
            "published_figure4_included": "",
            "published_figure4_cluster": "",
            "published_cluster21": "",
            "peru_focal": "Yes" if s(r.get("country", "")) == "Peru" else "No",
            "sequence_source_preference": f"data/pubmlst_st10042_assemblies/{pid}.fasta",
            "exact_id_dedup_status": "PubMLST_no_exact_Azevedo_match",
        })
        added += 1

    out = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)

    print(f"Azevedo backbone: {len(az)}")
    print(f"Exact-overlap PubMLST records absorbed into Azevedo backbone: {len(overlap_pids)}")
    print(f"PubMLST records appended after exact-ID de-duplication: {added}")
    print(f"Provisional combined candidate isolates: {len(out)}")
    print(f"Peru focal isolates: {(out['peru_focal'] == 'Yes').sum()}")
    print()
    print("PROVENANCE")
    print(out["provenance"].value_counts().to_string())
    print()
    print("PUBLISHED QC FAILURES AMONG EXACT OVERLAPS")
    q = out[
        out["pubmlst_id"].ne("") &
        out["published_figure4_included"].eq("No")
    ][["canonical_id","pubmlst_id","published_figure4_cluster"]]
    print(q.to_string(index=False) if len(q) else "None")
    print()
    print(f"Wrote: {a.out}")
    print("NOTE: 249 is provisional until sequence/metadata duplicate QC of unmatched records is complete.")


if __name__ == "__main__":
    main()
