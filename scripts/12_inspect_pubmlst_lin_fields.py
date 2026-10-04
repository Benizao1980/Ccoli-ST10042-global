#!/usr/bin/env python3
"""Inspect current ST10042 PubMLST JSON for cgMLST-v2/LIN-related fields.

This is deliberately schema-agnostic: BIGSdb isolate JSON can expose derived
scheme/classification fields in nested objects whose keys may change. We print
all paths containing cgST/cgMLST/LIN/cgc2/classification terms for the 41-record
snapshot, without altering the source JSON.
"""

import argparse
import json
import re
from pathlib import Path

PAT = re.compile(r"(lin|cgst|cgmlst|cgc2|classification)", re.I)


def walk(obj, path="$"):
    if isinstance(obj, dict):
        for k,v in obj.items():
            p=f"{path}.{k}"
            if PAT.search(str(k)):
                yield p,v
            yield from walk(v,p)
    elif isinstance(obj,list):
        for i,v in enumerate(obj):
            yield from walk(v,f"{path}[{i}]")


def isolate_id(rec):
    p=rec.get("provenance",{}) or {}
    return p.get("id","")


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--json",default="data/pubmlst_st10042_all.json")
    a=p.parse_args()
    records=json.loads(Path(a.json).read_text(encoding="utf-8"))

    hits=0
    records_with_hits=0
    for rec in records:
        found=list(walk(rec))
        if not found:
            continue
        records_with_hits += 1
        print(f"\nPUBMLST {isolate_id(rec)}")
        for path,value in found:
            hits += 1
            value_text=json.dumps(value,ensure_ascii=False)
            if len(value_text)>500:
                value_text=value_text[:500]+"..."
            print(f"{path} = {value_text}")

    print(f"\nRecords inspected: {len(records)}")
    print(f"Records with cgMLST/LIN/classification-related JSON fields: {records_with_hits}")
    print(f"Matching JSON paths printed: {hits}")


if __name__=="__main__":
    main()
