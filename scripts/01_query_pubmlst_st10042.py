#!/usr/bin/env python3
import argparse, csv, json, os, re
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

def pick(d,*names):
    for n in names:
        if n in d and d[n] not in (None,""):
            return d[n]
    return ""

def isolate_id(url):
    m=re.search(r"/isolates/(\d+)",url)
    return int(m.group(1)) if m else ""

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--country")
    p.add_argument("--out",default="pubmlst_st10042.tsv")
    p.add_argument("--json-out",default="pubmlst_st10042_records.json")
    p.add_argument("--download-contigs")
    a=p.parse_args()

    r=requests.post(SEARCH,headers={"Content-Type":"application/json"},
                    json={"scheme.1.ST":10042},timeout=120)
    r.raise_for_status()
    result=r.json()
    urls=result.get("isolates",[])
    print(f"PubMLST returned {result.get('records',len(urls))} ST10042 records")

    out=[]
    raw=[]
    for n,url in enumerate(urls,1):
        rr=requests.get(url,timeout=60); rr.raise_for_status()
        obj=rr.json(); raw.append(obj)
        prov=obj.get("provenance",{}) or {}
        iid=prov.get("id") or isolate_id(url)
        country=str(pick(prov,"country")).strip()
        if a.country and country.lower()!=a.country.lower():
            continue
        row={
            "pubmlst_id":iid,
            "isolate":pick(prov,"isolate","isolate_name","strain","strain_id"),
            "country":country,
            "continent":pick(prov,"continent"),
            "region":pick(prov,"region","state"),
            "town_or_city":pick(prov,"town_or_city","city"),
            "year":pick(prov,"year"),
            "month":pick(prov,"month"),
            "source":pick(prov,"source"),
            "species":pick(prov,"species"),
            "ST":10042,
            "record_url":url,
        }
        out.append(row)
        if a.download_contigs:
            d=Path(a.download_contigs); d.mkdir(parents=True,exist_ok=True)
            fa=requests.get(f"{DB}/isolates/{iid}/contigs_fasta?header=original_designation",timeout=120)
            if fa.ok and fa.text.startswith(">"):
                (d/f"{iid}.fasta").write_text(fa.text)
        if n%25==0: print(f"Fetched {n}/{len(urls)}")

    fields=["pubmlst_id","isolate","country","continent","region","town_or_city","year","month","source","species","ST","record_url"]
    with open(a.out,"w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t"); w.writeheader(); w.writerows(out)
    with open(a.json_out,"w",encoding="utf-8") as fh:
        json.dump(raw,fh,indent=2)
    print(f"Wrote {len(out)} records to {a.out}")

if __name__=="__main__":
    main()
