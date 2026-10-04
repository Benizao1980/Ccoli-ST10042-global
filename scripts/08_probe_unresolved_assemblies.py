#!/usr/bin/env python3
"""Probe ENA for public assemblies corresponding to Azevedo runs lacking public FASTQ.

This does not assume that raw-read availability and assembly availability are the same.
It queries the ENA assembly domain by BioProject/study and then performs conservative
exact string matching against Azevedo strain/ENA identifiers and run references.

Outputs:
- results/unresolved_study_assembly_inventory.tsv : every public assembly returned
  for studies represented among unresolved Azevedo runs.
- results/unresolved_azevedo_assembly_matches.tsv : candidate exact identifier matches
  back to unresolved Azevedo isolates.
"""

import argparse
import csv
import time
from pathlib import Path

import pandas as pd
import requests

ENA_SEARCH = "https://www.ebi.ac.uk/ena/portal/api/search"
ASSEMBLY_FIELDS = [
    "accession",
    "assembly_name",
    "assembly_title",
    "run_ref",
    "sample_accession",
    "secondary_sample_accession",
    "study_accession",
    "strain",
    "scientific_name",
]


def query_assemblies(study, retries=4):
    params = {
        "result": "assembly",
        "query": f'study_accession="{study}"',
        "fields": ",".join(ASSEMBLY_FIELDS),
        "format": "tsv",
        "limit": "0",
    }
    last = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(ENA_SEARCH, params=params, timeout=90)
            r.raise_for_status()
            text = r.text.strip()
            if not text:
                return []
            lines = [x for x in text.splitlines() if x]
            if len(lines) < 2:
                return []
            reader = csv.DictReader(lines, delimiter="\t")
            return list(reader)
        except requests.RequestException as e:
            last = e
            if attempt < retries:
                time.sleep(2 ** (attempt - 1))
    raise RuntimeError(f"ENA assembly query failed for {study}: {last}")


def exact_values(row):
    vals = set()
    for k, v in row.items():
        if v is None:
            continue
        s = str(v).strip()
        if not s:
            continue
        vals.add(s)
        # run_ref can contain multiple accessions separated by commas/semicolons/spaces
        if k == "run_ref":
            for sep in [",", ";", " "]:
                if sep in s:
                    vals.update(x.strip() for x in s.split(sep) if x.strip())
    return vals


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--azevedo", default="data/europe_217_manifest.tsv")
    p.add_argument("--ena-inventory", default="data/europe_ena_fastq_manifest.tsv")
    p.add_argument("--assembly-inventory", default="results/unresolved_study_assembly_inventory.tsv")
    p.add_argument("--matches", default="results/unresolved_azevedo_assembly_matches.tsv")
    a = p.parse_args()

    az = pd.read_csv(a.azevedo, sep="\t", dtype=str).fillna("")
    en = pd.read_csv(a.ena_inventory, sep="\t", dtype=str).fillna("")
    unresolved_ids = set(en.loc[en["status"].ne("resolved"), "Strain_ID"])
    u = az[az["Strain_ID"].isin(unresolved_ids)].copy()

    studies = sorted(x for x in u["Bioproject"].unique() if x)
    print(f"Unresolved Azevedo isolates: {len(u)}")
    print(f"Distinct studies/BioProjects to probe: {len(studies)}")

    all_assemblies = []
    for study in studies:
        recs = query_assemblies(study)
        print(f"{study}: {len(recs)} public assemblies")
        for rec in recs:
            rec["queried_study"] = study
            all_assemblies.append(rec)

    Path(a.assembly_inventory).parent.mkdir(parents=True, exist_ok=True)
    inv_fields = ["queried_study"] + ASSEMBLY_FIELDS
    with open(a.assembly_inventory, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=inv_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_assemblies)

    matches = []
    for _, ar in u.iterrows():
        wanted = {
            ar["Strain_ID"].strip(),
            ar["ENA/SRA_ID"].strip(),
            ar["Run Accession Number"].strip(),
        }
        wanted.discard("")
        for rec in all_assemblies:
            if rec.get("queried_study", "") != ar["Bioproject"].strip():
                continue
            vals = exact_values(rec)
            hit = sorted(wanted & vals)
            if hit:
                matches.append({
                    "Strain_ID": ar["Strain_ID"],
                    "Country": ar["Country"],
                    "Bioproject": ar["Bioproject"],
                    "Run Accession Number": ar["Run Accession Number"],
                    "ENA/SRA_ID": ar["ENA/SRA_ID"],
                    "Figure 4 cluster clean": ar["Figure 4 cluster clean"],
                    "Cluster 21 member": ar["Cluster 21 member"],
                    "assembly_accession": rec.get("accession", ""),
                    "assembly_name": rec.get("assembly_name", ""),
                    "assembly_title": rec.get("assembly_title", ""),
                    "assembly_strain": rec.get("strain", ""),
                    "assembly_run_ref": rec.get("run_ref", ""),
                    "matched_identifier": ";".join(hit),
                })

    match_fields = [
        "Strain_ID","Country","Bioproject","Run Accession Number","ENA/SRA_ID",
        "Figure 4 cluster clean","Cluster 21 member","assembly_accession",
        "assembly_name","assembly_title","assembly_strain","assembly_run_ref",
        "matched_identifier",
    ]
    with open(a.matches, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=match_fields, delimiter="\t")
        w.writeheader()
        w.writerows(matches)

    matched_ids = {m["Strain_ID"] for m in matches}
    c21 = u[u["Cluster 21 member"].eq("Yes")]
    c21_matched = set(c21["Strain_ID"]) & matched_ids

    print()
    print(f"Public assemblies returned across unresolved studies: {len(all_assemblies)}")
    print(f"Unresolved Azevedo isolates with exact assembly identifier match: {len(matched_ids)}")
    print(f"Unresolved cluster 21 isolates: {len(c21)}")
    print(f"Unresolved cluster 21 isolates rescued by exact assembly match: {len(c21_matched)}")
    print(f"Wrote: {a.assembly_inventory}")
    print(f"Wrote: {a.matches}")


if __name__ == "__main__":
    main()
