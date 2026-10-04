#!/usr/bin/env python3
"""Download and MD5-verify one Azevedo paired-end FASTQ entry."""

import argparse
import csv
import hashlib
import time
from pathlib import Path
import requests


def md5sum(path):
    h=hashlib.md5()
    with open(path,"rb") as fh:
        for block in iter(lambda:fh.read(8*1024*1024),b""):
            h.update(block)
    return h.hexdigest()


def get_row(path,index):
    with open(path,encoding="utf-8") as fh:
        rows=list(csv.DictReader(fh,delimiter="\t"))
    if index<0 or index>=len(rows):
        raise SystemExit(f"ERROR index {index} outside manifest range 0..{len(rows)-1}")
    return rows[index]


def download(url,dest,expected_md5,retries=4):
    dest=Path(dest)
    dest.parent.mkdir(parents=True,exist_ok=True)

    if dest.exists():
        got=md5sum(dest)
        if got==expected_md5:
            print(f"OK existing {dest}")
            return
        print(f"WARNING existing MD5 mismatch; removing {dest}")
        dest.unlink()

    tmp=dest.with_suffix(dest.suffix+".part")
    for attempt in range(1,retries+1):
        if tmp.exists():
            tmp.unlink()
        try:
            http=url if url.startswith("http") else "https://"+url
            print(f"Downloading {http}")
            with requests.get(http,stream=True,timeout=(60,600)) as r:
                r.raise_for_status()
                with open(tmp,"wb") as out:
                    for chunk in r.iter_content(8*1024*1024):
                        if chunk:
                            out.write(chunk)
            got=md5sum(tmp)
            if got!=expected_md5:
                raise RuntimeError(f"MD5 mismatch: expected {expected_md5}, got {got}")
            tmp.replace(dest)
            print(f"OK {dest}")
            return
        except Exception as e:
            print(f"WARNING attempt {attempt}/{retries} failed: {e}")
            if tmp.exists():
                tmp.unlink()
            if attempt==retries:
                raise
            time.sleep(2**attempt)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",default="results/azevedo_spades_manifest.tsv")
    p.add_argument("--index",type=int,required=True)
    a=p.parse_args()

    r=get_row(a.manifest,a.index)
    print(f"Sample {r['Strain_ID']} run {r['run_accession']}")
    download(r["read1_url"],r["read1_path"],r["read1_md5"])
    download(r["read2_url"],r["read2_path"],r["read2_md5"])


if __name__=="__main__":
    main()
