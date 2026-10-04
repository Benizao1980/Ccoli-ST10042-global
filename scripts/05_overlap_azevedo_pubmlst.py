#!/usr/bin/env python3
"""Find exact identifier overlap between Azevedo 217 and current PubMLST records.

The PubMLST JSON schema can vary between records. Instead of assuming one external-
accession field, this script recursively indexes every scalar value in each full
PubMLST JSON record, then tests the stable identifiers reported in the Azevedo
appendix (strain, ENA/SRA, BioProject and run accession).

Matches by run accession / ENA-SRA accession are considered high-confidence.
Strain-ID-only matches are reported separately and should be inspected manually.
"""

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path

def norm(x):
    if x is None:
        return ""
    return str(x).strip()

def scalar_index(obj, path="$", out=None):
    if out is None:
        out = defaultdict(list)
    if isinstance(obj, dict):
        for k, v in obj.items():
            scalar_index(v, f"{path}.{k}", out)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            scalar_index(v, f"{path}[{i}]", out)
    elif obj is not None:
        s = norm(obj)
        if s:
            out[s].append(path)
    return out

def pubmlst_id(record):
    prov = record.get("provenance", {}) or {}
    if prov.get("id") not in (None, ""):
        return str(prov["id"])
    for value in scalar_index(record):
        m = re.search(r"/isolates/(\d+)", value)
        if m:
            return m.group(1)
    return ""

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--azevedo", default="data/europe_217_manifest.tsv")
    p.add_argument("--pubmlst-json", default="data/pubmlst_st10042_all.json")
    p.add_argument("--out", default="results/azevedo_pubmlst_exact_overlap.tsv")
    p.add_argument("--unmatched", default="results/pubmlst_not_exactly_matched_to_azevedo.tsv")
    a = p.parse_args()

    with open(a.azevedo, encoding="utf-8") as fh:
        az = list(csv.DictReader(fh, delimiter="\t"))
    with open(a.pubmlst_json, encoding="utf-8") as fh:
        pub = json.load(fh)

    # Stable/semistable identifiers available in the published appendix.
    fields = [
        ("Run Accession Number", "run_accession", "high"),
        ("ENA/SRA_ID", "ena_sra", "high"),
        ("Bioproject", "bioproject", "supporting"),
        ("Strain_ID", "strain_id", "inspect"),
    ]

    pub_indexes = []
    for rec in pub:
        pub_indexes.append((pubmlst_id(rec), scalar_index(rec)))

    matches = []
    matched_pub = set()

    for row in az:
        for field, kind, confidence in fields:
            value = norm(row.get(field))
            if not value or value == "-":
                continue
            for pid, idx in pub_indexes:
                if value in idx:
                    matches.append({
                        "pubmlst_id": pid,
                        "azevedo_strain_id": norm(row.get("Strain_ID")),
                        "azevedo_country": norm(row.get("Country")),
                        "match_kind": kind,
                        "matched_value": value,
                        "confidence": confidence,
                        "pubmlst_json_paths": ";".join(idx[value]),
                    })
                    matched_pub.add(pid)

    # Deduplicate identical evidence rows.
    uniq = {}
    for m in matches:
        key = tuple(m[k] for k in m)
        uniq[key] = m
    matches = list(uniq.values())

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    fields_out = [
        "pubmlst_id","azevedo_strain_id","azevedo_country","match_kind",
        "matched_value","confidence","pubmlst_json_paths"
    ]
    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields_out, delimiter="\t")
        w.writeheader()
        w.writerows(sorted(matches, key=lambda x: (x["pubmlst_id"], x["match_kind"])))

    unmatched = []
    for rec in pub:
        pid = pubmlst_id(rec)
        if pid in matched_pub:
            continue
        prov = rec.get("provenance", {}) or {}
        unmatched.append({
            "pubmlst_id": pid,
            "isolate": norm(prov.get("isolate") or prov.get("isolate_name") or prov.get("strain") or prov.get("strain_id")),
            "country": norm(prov.get("country")),
            "year": norm(prov.get("year")),
            "source": norm(prov.get("source")),
        })

    with open(a.unmatched, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["pubmlst_id","isolate","country","year","source"],
            delimiter="\t"
        )
        w.writeheader()
        w.writerows(unmatched)

    high_pub = {m["pubmlst_id"] for m in matches if m["confidence"] == "high"}
    inspect_pub = {m["pubmlst_id"] for m in matches if m["confidence"] == "inspect"}
    print(f"Azevedo rows: {len(az)}")
    print(f"PubMLST records: {len(pub)}")
    print(f"PubMLST records with high-confidence accession overlap: {len(high_pub)}")
    print(f"Additional PubMLST records with strain-ID-only overlap: {len(inspect_pub - high_pub)}")
    print(f"PubMLST records without any exact identifier match: {len(pub) - len(matched_pub)}")
    print(f"Wrote: {a.out}")
    print(f"Wrote: {a.unmatched}")

if __name__ == "__main__":
    main()
