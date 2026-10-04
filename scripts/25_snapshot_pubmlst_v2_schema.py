#!/usr/bin/env python3
"""Freeze a reproducible snapshot of the current PubMLST Campylobacter cgMLST-v2 schema.

Downloads:
- scheme 8 metadata and locus list;
- all 1,142 locus allele FASTAs;
- current scheme profiles TSV;
- native LINcode definitions/nicknames;
- Cjc_cgc2 classification-scheme metadata.

The snapshot is resumable and records SHA256 checksums for every locus FASTA.
It does not modify PubMLST and sends requests sequentially.
"""

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path
from urllib.parse import unquote, urlparse

import requests

BASE="https://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef"
SCHEME_ID=8
EXPECTED_LOCI=1142
CLASSIFICATION_IDS=[11,12,13,14,15,16]


def get(url, *, accept=None, retries=5):
    headers={}
    if accept:
        headers["Accept"]=accept
    last=None
    for i in range(1,retries+1):
        try:
            r=requests.get(url,headers=headers,timeout=180)
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            last=e
            if i<retries:
                time.sleep(min(30,2**(i-1)))
    raise RuntimeError(f"GET failed {url}: {last}")


def locus_name(url):
    return unquote(urlparse(url).path.rstrip("/").split("/")[-1])


def sha256(data):
    return hashlib.sha256(data).hexdigest()


def count_fasta_records(data):
    return sum(1 for line in data.splitlines() if line.startswith(b">"))


def main():
    p=argparse.ArgumentParser()
    p.add_argument(
        "--out-dir",
        default="schema/pubmlst_campy_cgmlst_v2_scheme8_2026-10-04",
    )
    p.add_argument("--sleep",type=float,default=0.05)
    a=p.parse_args()

    out=Path(a.out_dir)
    alleles=out/"alleles"
    out.mkdir(parents=True,exist_ok=True)
    alleles.mkdir(parents=True,exist_ok=True)

    scheme=get(f"{BASE}/schemes/{SCHEME_ID}").json()
    (out/"scheme.json").write_text(json.dumps(scheme,indent=2)+"\n",encoding="utf-8")

    loci_obj=get(f"{BASE}/schemes/{SCHEME_ID}/loci").json()
    loci_urls=loci_obj.get("loci",[]) or []
    if len(loci_urls)!=EXPECTED_LOCI:
        raise SystemExit(f"ERROR expected {EXPECTED_LOCI} loci, got {len(loci_urls)}")
    (out/"scheme_loci.json").write_text(json.dumps(loci_obj,indent=2)+"\n",encoding="utf-8")

    # Profiles: preserve server-returned tab-delimited content verbatim.
    profiles_url=scheme.get("profiles_csv") or f"{BASE}/schemes/{SCHEME_ID}/profiles_csv"
    profiles=get(profiles_url,accept="text/tab-separated-values").content
    (out/"profiles.tsv").write_bytes(profiles)

    # Native LIN-related resources are live BIGSdb routes used by the current scheme.
    for name in ("lincodes","lincode_nicknames"):
        url=f"{BASE}/schemes/{SCHEME_ID}/{name}"
        try:
            payload=get(url).json()
            (out/f"{name}.json").write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")
        except Exception as e:
            (out/f"{name}.ERROR.txt").write_text(str(e)+"\n",encoding="utf-8")

    cls_dir=out/"classification_schemes"
    cls_dir.mkdir(exist_ok=True)
    for cid in CLASSIFICATION_IDS:
        url=f"{BASE}/classification_schemes/{cid}"
        payload=get(url).json()
        (cls_dir/f"{cid}.json").write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")

    rows=[]
    for i,url in enumerate(loci_urls,1):
        locus=locus_name(url)
        dest=alleles/f"{locus}.fasta"

        if dest.exists() and dest.stat().st_size>0:
            data=dest.read_bytes()
            status="cached"
        else:
            data=get(url.rstrip("/")+"/alleles_fasta",accept="text/plain").content
            if not data.strip():
                raise RuntimeError(f"Empty allele FASTA returned for {locus}")
            dest.write_bytes(data)
            status="downloaded"

        rows.append({
            "locus":locus,
            "locus_url":url,
            "alleles_fasta":str(dest),
            "allele_records":count_fasta_records(data),
            "bytes":len(data),
            "sha256":sha256(data),
        })
        if i==1 or i%50==0 or i==len(loci_urls):
            print(f"[{i}/{len(loci_urls)}] {locus} ({status})",flush=True)
        time.sleep(a.sleep)

    with open(out/"locus_manifest.tsv","w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    snapshot={
        "database":BASE,
        "scheme_id":SCHEME_ID,
        "expected_loci":EXPECTED_LOCI,
        "downloaded_loci":len(rows),
        "scheme_description":scheme.get("description"),
        "scheme_last_added":scheme.get("last_added"),
        "scheme_last_updated":scheme.get("last_updated"),
        "profiles_sha256":sha256(profiles),
        "profiles_bytes":len(profiles),
        "classification_scheme_ids":CLASSIFICATION_IDS,
    }
    (out/"snapshot.json").write_text(json.dumps(snapshot,indent=2)+"\n",encoding="utf-8")

    print("\nSNAPSHOT COMPLETE")
    print(f"Loci: {len(rows)}")
    print(f"Profiles bytes: {len(profiles)}")
    print(f"Output: {out}")
    print("Use locus_manifest.tsv checksums to document the exact schema snapshot.")


if __name__=="__main__":
    main()
