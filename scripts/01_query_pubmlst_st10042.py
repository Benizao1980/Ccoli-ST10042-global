#!/usr/bin/env python3
"""Query PubMLST/BIGSdb for current Campylobacter coli ST10042 records.

Authentication
--------------
PubMLST requires authentication to access records added after 31 Dec 2024.
Create a personal API key in your PubMLST/BIGSdb profile and export it as:

    export PUBMLST_API_KEY='...'

The key is sent as the X-API-Key header on search, record, and FASTA requests.
"""

import argparse
import csv
import json
import os
import re
from pathlib import Path

import requests

DB = "https://rest.pubmlst.org/db/pubmlst_campylobacter_isolates"
SEARCH = f"{DB}/isolates/search?return_all=1"


def auth_headers(json_content=False):
    headers = {}
    if json_content:
        headers["Content-Type"] = "application/json"
    key = os.environ.get("PUBMLST_API_KEY", "").strip()
    if key:
        headers["X-API-Key"] = key
    return headers


def pick(d, *names):
    for n in names:
        if n in d and d[n] not in (None, ""):
            return d[n]
    return ""


def isolate_id(url):
    m = re.search(r"/isolates/(\d+)", url)
    return int(m.group(1)) if m else ""


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--country")
    p.add_argument("--out", default="pubmlst_st10042.tsv")
    p.add_argument("--json-out", default="pubmlst_st10042_records.json")
    p.add_argument("--download-contigs")
    a = p.parse_args()

    if not os.environ.get("PUBMLST_API_KEY", "").strip():
        print(
            "WARNING: PUBMLST_API_KEY is not set. "
            "Post-2024 records may be hidden by PubMLST access policy."
        )

    session = requests.Session()

    r = session.post(
        SEARCH,
        headers=auth_headers(json_content=True),
        json={"scheme.1.ST": 10042},
        timeout=120,
    )
    if r.status_code == 401:
        raise SystemExit(
            "ERROR: PubMLST returned HTTP 401 Unauthorized. "
            "PUBMLST_API_KEY is present but was not accepted. "
            "Check/regenerate the personal API key in your PubMLST/BIGSdb profile, "
            "then export the replacement key before rerunning."
        )
    r.raise_for_status()
    result = r.json()
    urls = result.get("isolates", [])
    print(f"PubMLST returned {result.get('records', len(urls))} ST10042 records")

    out = []
    selected_raw = []

    for n, url in enumerate(urls, 1):
        rr = session.get(url, headers=auth_headers(), timeout=60)
        rr.raise_for_status()
        obj = rr.json()

        prov = obj.get("provenance", {}) or {}
        iid = prov.get("id") or isolate_id(url)
        country = str(pick(prov, "country")).strip()

        if a.country and country.lower() != a.country.lower():
            continue

        selected_raw.append(obj)

        row = {
            "pubmlst_id": iid,
            "isolate": pick(prov, "isolate", "isolate_name", "strain", "strain_id"),
            "country": country,
            "continent": pick(prov, "continent"),
            "region": pick(prov, "region", "state"),
            "town_or_city": pick(prov, "town_or_city", "city"),
            "year": pick(prov, "year"),
            "month": pick(prov, "month"),
            "source": pick(prov, "source"),
            "species": pick(prov, "species"),
            "ST": 10042,
            "record_url": url,
        }
        out.append(row)

        if a.download_contigs:
            d = Path(a.download_contigs)
            d.mkdir(parents=True, exist_ok=True)
            fa = session.get(
                f"{DB}/isolates/{iid}/contigs_fasta?header=original_designation",
                headers=auth_headers(),
                timeout=120,
            )
            if fa.ok and fa.text.startswith(">"):
                (d / f"{iid}.fasta").write_text(fa.text, encoding="utf-8")

        if n % 25 == 0:
            print(f"Visited {n}/{len(urls)} returned records")

    fields = [
        "pubmlst_id",
        "isolate",
        "country",
        "continent",
        "region",
        "town_or_city",
        "year",
        "month",
        "source",
        "species",
        "ST",
        "record_url",
    ]

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    Path(a.json_out).parent.mkdir(parents=True, exist_ok=True)

    with open(a.out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out)

    with open(a.json_out, "w", encoding="utf-8") as fh:
        json.dump(selected_raw, fh, indent=2)

    label = f" for country={a.country}" if a.country else ""
    print(f"Wrote {len(out)} ST10042 records{label} to {a.out}")


if __name__ == "__main__":
    main()
