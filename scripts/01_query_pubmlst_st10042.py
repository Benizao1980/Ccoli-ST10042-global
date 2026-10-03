#!/usr/bin/env python3
"""Query PubMLST/BIGSdb for current Campylobacter coli ST10042 records.

Full contemporary PubMLST searches require OAuth because isolate searches use
HTTP POST and PubMLST does not permit personal API keys for POST requests.

One-time OAuth setup (interactive, on the login node):

    bigsdb-downloader \
      --key_name PubMLST \
      --site PubMLST \
      --db pubmlst_campylobacter_isolates \
      --token_dir ~/.bigsdb_tokens \
      --setup

The setup asks for the PubMLST OAuth client key and client secret, then opens an
authorization step and stores the resulting access token outside this repository.
Subsequent runs obtain/renew short-lived session tokens automatically.
"""

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import requests

DB = "https://rest.pubmlst.org/db/pubmlst_campylobacter_isolates"
SEARCH = f"{DB}/isolates/search?return_all=1"


def pick(d, *names):
    for n in names:
        if n in d and d[n] not in (None, ""):
            return d[n]
    return ""


def isolate_id(url):
    m = re.search(r"/isolates/(\d+)", url)
    return int(m.group(1)) if m else ""


def oauth_fetch(url, method="GET", json_body=None, key_name="PubMLST", token_dir=None):
    """Fetch a BIGSdb resource using the official OAuth-aware downloader."""
    exe = shutil.which("bigsdb-downloader")
    if not exe:
        raise SystemExit(
            "ERROR: bigsdb-downloader is not installed. "
            "Run: conda env update -f environment.yml"
        )

    token_dir = Path(token_dir or "~/.bigsdb_tokens").expanduser()
    if not token_dir.exists():
        raise SystemExit(
            f"ERROR: OAuth token directory does not exist: {token_dir}\n"
            "Run the one-time bigsdb-downloader --setup command first."
        )

    fd, tmp_name = tempfile.mkstemp(prefix="bigsdb_", suffix=".out")
    os.close(fd)
    tmp = Path(tmp_name)

    cmd = [
        exe,
        "--key_name", key_name,
        "--site", "PubMLST",
        "--token_dir", str(token_dir),
        "--cron",
        "--url", url,
        "--output_file", str(tmp),
    ]
    if method.upper() == "POST":
        cmd += ["--method", "POST", "--json_body", json.dumps(json_body or {})]

    try:
        p = subprocess.run(cmd, text=True, capture_output=True)
        if p.returncode != 0:
            detail = (p.stderr or p.stdout or "").strip()
            raise SystemExit(
                "ERROR: authenticated PubMLST request failed via "
                f"bigsdb-downloader (exit {p.returncode}).\n{detail}"
            )
        return tmp.read_text(encoding="utf-8")
    finally:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass


def anonymous_fetch(url, method="GET", json_body=None):
    """Anonymous mode for public/pre-policy data only."""
    if method.upper() == "POST":
        r = requests.post(
            url,
            headers={"Content-Type": "application/json"},
            json=json_body or {},
            timeout=120,
        )
    else:
        r = requests.get(url, timeout=120)
    r.raise_for_status()
    return r.text


def fetch(url, method, json_body, auth, key_name, token_dir):
    if auth == "oauth":
        return oauth_fetch(
            url,
            method=method,
            json_body=json_body,
            key_name=key_name,
            token_dir=token_dir,
        )
    return anonymous_fetch(url, method=method, json_body=json_body)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--country")
    p.add_argument("--out", default="pubmlst_st10042.tsv")
    p.add_argument("--json-out", default="pubmlst_st10042_records.json")
    p.add_argument("--download-contigs")
    p.add_argument(
        "--auth",
        choices=["oauth", "anonymous"],
        default="oauth",
        help="OAuth is required for the complete contemporary dataset.",
    )
    p.add_argument(
        "--oauth-key-name",
        default=os.environ.get("PUBMLST_KEY_NAME", "PubMLST"),
    )
    p.add_argument(
        "--oauth-token-dir",
        default=os.environ.get("PUBMLST_TOKEN_DIR", "~/.bigsdb_tokens"),
    )
    a = p.parse_args()

    if a.auth == "anonymous":
        print(
            "WARNING: anonymous mode may omit records hidden by PubMLST's "
            "current access policy. Do not use it for the final global dataset."
        )

    search_text = fetch(
        SEARCH,
        method="POST",
        json_body={"scheme.1.ST": 10042},
        auth=a.auth,
        key_name=a.oauth_key_name,
        token_dir=a.oauth_token_dir,
    )
    result = json.loads(search_text)
    urls = result.get("isolates", [])
    print(f"PubMLST returned {result.get('records', len(urls))} ST10042 records")

    out = []
    selected_raw = []

    for n, url in enumerate(urls, 1):
        obj = json.loads(
            fetch(
                url,
                method="GET",
                json_body=None,
                auth=a.auth,
                key_name=a.oauth_key_name,
                token_dir=a.oauth_token_dir,
            )
        )

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
            fasta = fetch(
                f"{DB}/isolates/{iid}/contigs_fasta?header=original_designation",
                method="GET",
                json_body=None,
                auth=a.auth,
                key_name=a.oauth_key_name,
                token_dir=a.oauth_token_dir,
            )
            if fasta.startswith(">"):
                (d / f"{iid}.fasta").write_text(fasta, encoding="utf-8")

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
