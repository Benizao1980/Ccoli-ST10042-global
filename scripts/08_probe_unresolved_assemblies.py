#!/usr/bin/env python3
"""Probe ENA for public assemblies corresponding to Azevedo runs lacking public FASTQ.

The Azevedo supplement contains one malformed BioProject value (FR-1 has an ERR run
accession in the BioProject column). Project-wide assembly queries are therefore made
only for syntactically valid PRJ*/ERP/DRP/SRP study accessions. Invalid project values
are reported and, where possible, probed once by exact run_ref instead.

Outputs:
- results/unresolved_study_assembly_inventory.tsv
- results/unresolved_azevedo_assembly_matches.tsv
"""

import argparse
import csv
import re
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

PROJECT_RE = re.compile(r"^(?:PRJ[A-Z]{2}\d+|[ESD]RP\d+)$")
RUN_RE = re.compile(r"^[ESD]RR\d+$")


def portal_query(query, retries=4):
    params = {
        "result": "assembly",
        "query": query,
        "fields": ",".join(ASSEMBLY_FIELDS),
        "format": "tsv",
        "limit": "0",
    }
    last = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(ENA_SEARCH, params=params, timeout=90)
            # A malformed/unsupported query is not transient: report it to caller.
            if r.status_code == 400:
                return [], f"HTTP400: {r.text.strip()[:300]}"
            r.raise_for_status()
            text = r.text.strip()
            if not text:
                return [], None
            lines = [x for x in text.splitlines() if x]
            if len(lines) < 2:
                return [], None
            return list(csv.DictReader(lines, delimiter="\t")), None
        except requests.RequestException as e:
            last = e
            if attempt < retries:
                time.sleep(2 ** (attempt - 1))
    return [], f"{type(last).__name__}: {last}"


def exact_values(row):
    vals = set()
    for k, v in row.items():
        if v is None:
            continue
        s = str(v).strip()
        if not s:
            continue
        vals.add(s)
        if k == "run_ref":
            # ENA may represent multiple run references in one field.
            vals.update(x for x in re.split(r"[,;\s]+", s) if x)
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

    project_values = sorted(x.strip() for x in u["Bioproject"].unique() if x.strip())
    valid_projects = [x for x in project_values if PROJECT_RE.match(x)]
    invalid_projects = [x for x in project_values if not PROJECT_RE.match(x)]

    print(f"Unresolved Azevedo isolates: {len(u)}")
    print(f"Valid studies/BioProjects to probe: {len(valid_projects)}")
    print(f"Invalid/non-project BioProject values: {len(invalid_projects)}")
    for x in invalid_projects:
        affected = u.loc[u["Bioproject"].eq(x), ["Strain_ID","Run Accession Number"]]
        print(f"WARNING invalid BioProject value {x!r}:")
        print(affected.to_string(index=False))

    all_assemblies = []
    query_errors = []

    # Bulk project queries: only a handful of API calls.
    for study in valid_projects:
        recs, err = portal_query(f'study_accession="{study}"')
        if err:
            print(f"WARNING assembly query failed for {study}: {err}")
            query_errors.append((study, err))
            continue
        print(f"{study}: {len(recs)} public assemblies")
        for rec in recs:
            rec["query_kind"] = "study_accession"
            rec["queried_value"] = study
            all_assemblies.append(rec)

    # Rescue malformed project metadata with a single exact run-ref query per affected row.
    # This avoids trying to use an ERR accession as a study accession.
    for _, ar in u[u["Bioproject"].isin(invalid_projects)].iterrows():
        run = ar["Run Accession Number"].strip()
        if not RUN_RE.match(run):
            print(f"WARNING cannot run-ref probe {ar['Strain_ID']}: invalid run {run!r}")
            continue
        recs, err = portal_query(f'run_ref="{run}"')
        if err:
            print(f"WARNING run_ref assembly query failed for {run}: {err}")
            query_errors.append((run, err))
            continue
        print(f"{run}: {len(recs)} public assemblies by exact run_ref")
        for rec in recs:
            rec["query_kind"] = "run_ref"
            rec["queried_value"] = run
            all_assemblies.append(rec)

    # De-duplicate assembly records returned by multiple query routes.
    uniq = {}
    for rec in all_assemblies:
        key = (
            rec.get("accession", ""),
            rec.get("assembly_name", ""),
            rec.get("run_ref", ""),
            rec.get("study_accession", ""),
        )
        uniq[key] = rec
    all_assemblies = list(uniq.values())

    Path(a.assembly_inventory).parent.mkdir(parents=True, exist_ok=True)
    inv_fields = ["query_kind", "queried_value"] + ASSEMBLY_FIELDS
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
            vals = exact_values(rec)
            hit = sorted(wanted & vals)
            if not hit:
                continue

            # If this came from a study query, require the Azevedo project to agree.
            qkind = rec.get("query_kind", "")
            if qkind == "study_accession":
                if ar["Bioproject"].strip() != rec.get("queried_value", ""):
                    continue

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
                "assembly_study_accession": rec.get("study_accession", ""),
                "matched_identifier": ";".join(hit),
            })

    match_fields = [
        "Strain_ID","Country","Bioproject","Run Accession Number","ENA/SRA_ID",
        "Figure 4 cluster clean","Cluster 21 member","assembly_accession",
        "assembly_name","assembly_title","assembly_strain","assembly_run_ref",
        "assembly_study_accession","matched_identifier",
    ]
    with open(a.matches, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=match_fields, delimiter="\t")
        w.writeheader()
        w.writerows(matches)

    matched_ids = {m["Strain_ID"] for m in matches}
    c21 = u[u["Cluster 21 member"].eq("Yes")]
    c21_matched = set(c21["Strain_ID"]) & matched_ids

    print()
    print(f"Public assemblies returned: {len(all_assemblies)}")
    print(f"Unresolved Azevedo isolates with exact assembly identifier match: {len(matched_ids)}")
    print(f"Unresolved cluster 21 isolates: {len(c21)}")
    print(f"Unresolved cluster 21 isolates rescued by exact assembly match: {len(c21_matched)}")
    print(f"Assembly-query errors: {len(query_errors)}")
    print(f"Wrote: {a.assembly_inventory}")
    print(f"Wrote: {a.matches}")


if __name__ == "__main__":
    main()
