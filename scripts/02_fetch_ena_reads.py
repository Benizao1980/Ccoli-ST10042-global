#!/usr/bin/env python3
import argparse
import csv
import hashlib
import time
from pathlib import Path

import requests

ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def ena_lookup(run, retries=4):
    params = {
        "accession": run,
        "result": "read_run",
        "fields": "run_accession,fastq_ftp,fastq_md5,fastq_bytes",
        "format": "tsv",
    }
    last_error = ""
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(ENA, params=params, timeout=60)
            r.raise_for_status()
            lines = [x for x in r.text.strip().splitlines() if x]
            if len(lines) < 2:
                return None, "no_ENA_FASTQ"
            rec = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
            if not rec.get("run_accession"):
                return None, "no_run_accession"
            return rec, "resolved"
        except requests.RequestException as e:
            last_error = f"{type(e).__name__}: {e}"
            if attempt < retries:
                time.sleep(2 ** (attempt - 1))
    return None, f"request_failed: {last_error}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", default="data/europe_217_manifest.tsv")
    p.add_argument("--resolved", default="data/europe_ena_fastq_manifest.tsv")
    p.add_argument("--unresolved", default="results/europe_ena_unresolved.tsv")
    p.add_argument("--download", action="store_true")
    p.add_argument("--outdir", default="data/europe_reads")
    a = p.parse_args()

    with open(a.manifest, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    all_rows = []
    unresolved = []

    for n, row in enumerate(rows, 1):
        run = row["Run Accession Number"].strip()
        rec, status = ena_lookup(run)

        base = {
            "Strain_ID": row["Strain_ID"],
            "Country": row["Country"],
            "Figure4_cluster": row["Figure 4 cluster clean"],
            "published_figure4_included": row["Included in Figure 4"],
            "requested_run_accession": run,
            "status": status,
        }

        if rec is None:
            print(f"WARNING {status}: {run}")
            unresolved.append(base.copy())
            all_rows.append({
                **base,
                "run_accession": "",
                "fastq_ftp": "",
                "fastq_md5": "",
                "fastq_bytes": "",
            })
            continue

        rec.update(base)
        all_rows.append(rec)

        if a.download:
            sd = Path(a.outdir) / row["Strain_ID"]
            sd.mkdir(parents=True, exist_ok=True)
            ftps = [x for x in rec.get("fastq_ftp", "").split(";") if x]
            md5s = [x for x in rec.get("fastq_md5", "").split(";") if x]
            for j, ftp in enumerate(ftps):
                url = "https://" + ftp if not ftp.startswith("http") else ftp
                dest = sd / Path(ftp).name
                if not dest.exists():
                    with requests.get(url, stream=True, timeout=600) as rr:
                        rr.raise_for_status()
                        with open(dest, "wb") as out:
                            for chunk in rr.iter_content(1024 * 1024):
                                if chunk:
                                    out.write(chunk)
                if j < len(md5s) and md5sum(dest) != md5s[j]:
                    raise RuntimeError(f"MD5 mismatch: {dest}")

        if n % 25 == 0:
            print(f"Resolved/checkpointed {n}/{len(rows)}")

    fields = [
        "Strain_ID",
        "Country",
        "Figure4_cluster",
        "published_figure4_included",
        "requested_run_accession",
        "status",
        "run_accession",
        "fastq_ftp",
        "fastq_md5",
        "fastq_bytes",
    ]

    Path(a.resolved).parent.mkdir(parents=True, exist_ok=True)
    with open(a.resolved, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)

    Path(a.unresolved).parent.mkdir(parents=True, exist_ok=True)
    with open(a.unresolved, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=[
                "Strain_ID",
                "Country",
                "Figure4_cluster",
                "published_figure4_included",
                "requested_run_accession",
                "status",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(unresolved)

    resolved_n = sum(x["status"] == "resolved" for x in all_rows)
    print()
    print(f"Input Azevedo isolates: {len(rows)}")
    print(f"ENA runs resolved: {resolved_n}")
    print(f"ENA runs unresolved: {len(unresolved)}")
    print(f"Wrote complete inventory: {a.resolved}")
    print(f"Wrote unresolved report: {a.unresolved}")


if __name__ == "__main__":
    main()
