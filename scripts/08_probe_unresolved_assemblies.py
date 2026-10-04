#!/usr/bin/env python3
"""Probe ENA for public assembly products corresponding to unresolved Azevedo runs.

ENA's Portal API fields vary by result type and can change over time. This script
therefore discovers current return/search fields at runtime rather than hard-coding
legacy fields such as assembly.run_ref.

Routes used:
1) result=assembly, queried by study_accession where possible;
2) result=analysis, restricted to SEQUENCE_ASSEMBLY where supported;
3) conservative exact matching to Azevedo isolate identifiers using returned
   strain/name/title/sample metadata.

The source Azevedo manifest is never rewritten.
"""

import argparse
import csv
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests

BASE = "https://www.ebi.ac.uk/ena/portal/api"
PROJECT_RE = re.compile(r"^(?:PRJ[A-Z]{2}\d+|[ESD]RP\d+)$")
RUN_RE = re.compile(r"^[ESD]RR\d+$")

PREFERRED_RETURN_FIELDS = {
    "assembly": [
        "accession",
        "study_accession",
        "sample_accession",
        "secondary_sample_accession",
        "assembly_name",
        "assembly_title",
        "strain",
        "scientific_name",
        "assembly_level",
        "genome_representation",
    ],
    "analysis": [
        "accession",
        "analysis_accession",
        "study_accession",
        "sample_accession",
        "secondary_sample_accession",
        "analysis_type",
        "analysis_title",
        "submitted_ftp",
        "submitted_md5",
        "scientific_name",
        "strain",
    ],
}


def get_metadata_response(endpoint, params, retries=4):
    """Return ENA metadata endpoint response as text plus content type.

    ENA Portal metadata endpoints such as returnFields/searchFields are not
    guaranteed to return JSON by default; deployments may return newline/TSV text.
    """
    last = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(
                f"{BASE}/{endpoint}",
                params={**params, "format": "json"},
                timeout=90,
            )
            r.raise_for_status()
            return r.text, r.headers.get("content-type", "")
        except requests.RequestException as e:
            last = e
            if attempt < retries:
                time.sleep(2 ** (attempt - 1))
    raise RuntimeError(f"ENA {endpoint} failed: {last}")


def _fields_from_json(data):
    fields = set()
    if isinstance(data, list):
        for item in data:
            if isinstance(item, str):
                fields.add(item)
            elif isinstance(item, dict):
                for key in ("columnId", "fieldName", "name", "id"):
                    if item.get(key):
                        fields.add(str(item[key]))
                        break
    elif isinstance(data, dict):
        for value in data.values():
            if isinstance(value, list):
                fields |= _fields_from_json(value)
    return fields


def _fields_from_text(text):
    """Parse ENA field metadata returned as newline or TSV text."""
    fields = set()
    header_tokens = {
        "field", "fieldname", "field_name", "name", "id",
        "description", "type", "searchfield", "returnfield", "columnid",
    }
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        # HTML here means the metadata endpoint did not return its advertised data.
        if line.startswith("<"):
            raise RuntimeError(
                "ENA field-discovery endpoint returned HTML rather than field metadata: "
                + line[:120]
            )
        first = line.split("\t", 1)[0].strip()
        # Some deployments emit a comma-separated header/list.
        if "\t" not in line and "," in first:
            candidates = [x.strip() for x in first.split(",")]
        else:
            candidates = [first]
        for value in candidates:
            if not value:
                continue
            if value.lower() in header_tokens:
                continue
            if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", value):
                fields.add(value)
    return fields


def discover_fields(result, endpoint):
    text, content_type = get_metadata_response(endpoint, {"result": result})

    # Prefer JSON when ENA honours format=json.
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        data = None

    fields = _fields_from_json(data) if data is not None else _fields_from_text(text)
    if not fields:
        snippet = text[:300].replace("\n", "\\n")
        raise RuntimeError(
            f"Could not parse any fields from ENA {endpoint} for result={result}; "
            f"content-type={content_type!r}; response={snippet!r}"
        )
    return fields


def portal_search(result, query, fields, retries=4):
    params = {
        "result": result,
        "query": query,
        "fields": ",".join(fields),
        "format": "tsv",
        "limit": "0",
    }
    last = None
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(f"{BASE}/search", params=params, timeout=90)
            if r.status_code == 400:
                return [], f"HTTP400: {r.text.strip()[:500]}"
            r.raise_for_status()
            lines = [x for x in r.text.strip().splitlines() if x]
            if len(lines) < 2:
                return [], None
            return list(csv.DictReader(lines, delimiter="\t")), None
        except requests.RequestException as e:
            last = e
            if attempt < retries:
                time.sleep(2 ** (attempt - 1))
    return [], f"{type(last).__name__}: {last}"


def tokens(value):
    if value is None:
        return set()
    s = str(value).strip()
    if not s:
        return set()
    out = {s}
    out.update(x for x in re.split(r"[,;|\s]+", s) if x)
    return out


def record_values(rec):
    vals = set()
    for value in rec.values():
        vals |= tokens(value)
    return vals


def exact_text_identifier_hits(rec, wanted):
    """Conservative identifier matching.

    Exact scalar/token equality is accepted. For free-text title fields we only
    accept identifiers as delimiter-bounded tokens, avoiding DE-1 matching DE-10.
    """
    hits = record_values(rec) & wanted
    for key in ("assembly_title", "analysis_title"):
        text = str(rec.get(key, "") or "")
        for ident in wanted:
            if not ident:
                continue
            pat = rf"(?<![A-Za-z0-9]){re.escape(ident)}(?![A-Za-z0-9])"
            if re.search(pat, text):
                hits.add(ident)
    return hits


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--azevedo", default="data/europe_217_manifest.tsv")
    p.add_argument("--ena-inventory", default="data/europe_ena_fastq_manifest.tsv")
    p.add_argument("--assembly-inventory", default="results/unresolved_study_assembly_inventory.tsv")
    p.add_argument("--analysis-inventory", default="results/unresolved_study_analysis_inventory.tsv")
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

    # Discover the API schema currently exposed by ENA.
    schemas = {}
    for result in ("assembly", "analysis"):
        return_fields = discover_fields(result, "returnFields")
        search_fields = discover_fields(result, "searchFields")
        selected = [f for f in PREFERRED_RETURN_FIELDS[result] if f in return_fields]
        # accession is useful and normally returned even when not requested, but ask if allowed.
        if not selected:
            raise RuntimeError(f"No usable ENA return fields discovered for result={result}")
        schemas[result] = {
            "return": return_fields,
            "search": search_fields,
            "selected": selected,
        }
        print(
            f"ENA {result}: discovered {len(return_fields)} return fields, "
            f"{len(search_fields)} search fields; requesting {','.join(selected)}"
        )

    query_errors = []
    assembly_records = []
    analysis_records = []

    for study in valid_projects:
        # Genome assembly catalogue.
        if "study_accession" in schemas["assembly"]["search"]:
            recs, err = portal_search(
                "assembly",
                f'study_accession="{study}"',
                schemas["assembly"]["selected"],
            )
            if err:
                print(f"WARNING assembly query failed for {study}: {err}")
                query_errors.append(("assembly", study, err))
            else:
                print(f"{study}: {len(recs)} public genome-assembly records")
                for rec in recs:
                    rec["queried_study"] = study
                    assembly_records.append(rec)
        else:
            print("WARNING ENA assembly schema does not expose study_accession as searchable")
            query_errors.append(("assembly", study, "study_accession not searchable"))

        # Submitted analysis objects can also carry sequence assemblies/files.
        if "study_accession" in schemas["analysis"]["search"]:
            q = f'study_accession="{study}"'
            if "analysis_type" in schemas["analysis"]["search"]:
                q += ' AND analysis_type="SEQUENCE_ASSEMBLY"'
            recs, err = portal_search("analysis", q, schemas["analysis"]["selected"])
            if err:
                print(f"WARNING analysis query failed for {study}: {err}")
                query_errors.append(("analysis", study, err))
            else:
                print(f"{study}: {len(recs)} public sequence-assembly analysis records")
                for rec in recs:
                    rec["queried_study"] = study
                    analysis_records.append(rec)
        else:
            print("WARNING ENA analysis schema does not expose study_accession as searchable")
            query_errors.append(("analysis", study, "study_accession not searchable"))

    # FR-1 has an ERR value in the published BioProject column. We cannot safely
    # translate that to a project. If current assembly/analysis schemas happen to
    # expose a run field, try it; otherwise report that the route is unavailable.
    for _, ar in u[u["Bioproject"].isin(invalid_projects)].iterrows():
        run = ar["Run Accession Number"].strip()
        if not RUN_RE.match(run):
            continue
        attempted = False
        for result, bucket in (("assembly", assembly_records), ("analysis", analysis_records)):
            candidate_run_fields = [
                f for f in ("run_accession", "run_ref", "run") if f in schemas[result]["search"]
            ]
            if not candidate_run_fields:
                continue
            attempted = True
            field = candidate_run_fields[0]
            recs, err = portal_search(
                result,
                f'{field}="{run}"',
                schemas[result]["selected"],
            )
            if err:
                print(f"WARNING {result} run query failed for {run}: {err}")
                query_errors.append((result, run, err))
            else:
                print(f"{run}: {len(recs)} public {result} records by {field}")
                for rec in recs:
                    rec["queried_study"] = ""
                    bucket.append(rec)
        if not attempted:
            print(
                f"NOTE {run}: current ENA assembly/analysis search schemas expose no "
                "run-linked field; FR-1 cannot be rescued by exact run query here."
            )

    # De-duplicate within each domain.
    def dedup(records):
        uniq = {}
        for rec in records:
            key = tuple(sorted((k, str(v)) for k, v in rec.items() if k != "queried_study"))
            uniq[key] = rec
        return list(uniq.values())

    assembly_records = dedup(assembly_records)
    analysis_records = dedup(analysis_records)

    Path(a.assembly_inventory).parent.mkdir(parents=True, exist_ok=True)

    def write_inventory(path, records, result):
        fields = ["queried_study"] + schemas[result]["selected"]
        # Preserve only columns that can actually be returned plus our provenance column.
        fields = list(dict.fromkeys(fields))
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(records)

    write_inventory(a.assembly_inventory, assembly_records, "assembly")
    write_inventory(a.analysis_inventory, analysis_records, "analysis")

    matches = []
    for _, ar in u.iterrows():
        wanted = {
            ar["Strain_ID"].strip(),
            ar["ENA/SRA_ID"].strip(),
            ar["Run Accession Number"].strip(),
        }
        wanted.discard("")
        study = ar["Bioproject"].strip()

        for result, records in (("assembly", assembly_records), ("analysis", analysis_records)):
            for rec in records:
                # For valid projects, keep matches inside the same queried study.
                if PROJECT_RE.match(study) and rec.get("queried_study", "") != study:
                    continue
                hit = sorted(exact_text_identifier_hits(rec, wanted))
                if not hit:
                    continue
                matches.append({
                    "Strain_ID": ar["Strain_ID"],
                    "Country": ar["Country"],
                    "Bioproject": ar["Bioproject"],
                    "Run Accession Number": ar["Run Accession Number"],
                    "ENA/SRA_ID": ar["ENA/SRA_ID"],
                    "Figure 4 cluster clean": ar["Figure 4 cluster clean"],
                    "Cluster 21 member": ar["Cluster 21 member"],
                    "ena_result_type": result,
                    "record_accession": rec.get("accession", rec.get("analysis_accession", "")),
                    "sample_accession": rec.get("sample_accession", ""),
                    "record_name": rec.get("assembly_name", ""),
                    "record_title": rec.get("assembly_title", rec.get("analysis_title", "")),
                    "strain": rec.get("strain", ""),
                    "submitted_ftp": rec.get("submitted_ftp", ""),
                    "matched_identifier": ";".join(hit),
                })

    match_fields = [
        "Strain_ID","Country","Bioproject","Run Accession Number","ENA/SRA_ID",
        "Figure 4 cluster clean","Cluster 21 member","ena_result_type",
        "record_accession","sample_accession","record_name","record_title",
        "strain","submitted_ftp","matched_identifier",
    ]
    with open(a.matches, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=match_fields, delimiter="\t")
        w.writeheader()
        w.writerows(matches)

    matched_ids = {m["Strain_ID"] for m in matches}
    c21 = u[u["Cluster 21 member"].eq("Yes")]
    c21_matched = set(c21["Strain_ID"]) & matched_ids

    print()
    print(f"Public genome-assembly records returned: {len(assembly_records)}")
    print(f"Public SEQUENCE_ASSEMBLY analysis records returned: {len(analysis_records)}")
    print(f"Unresolved Azevedo isolates with conservative exact identifier match: {len(matched_ids)}")
    print(f"Unresolved cluster 21 isolates: {len(c21)}")
    print(f"Unresolved cluster 21 isolates rescued by exact match: {len(c21_matched)}")
    print(f"Query errors: {len(query_errors)}")
    print(f"Wrote: {a.assembly_inventory}")
    print(f"Wrote: {a.analysis_inventory}")
    print(f"Wrote: {a.matches}")


if __name__ == "__main__":
    main()
