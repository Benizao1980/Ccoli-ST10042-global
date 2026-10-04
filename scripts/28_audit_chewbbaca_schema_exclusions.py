#!/usr/bin/env python3
"""Audit chewBBACA PrepExternalSchema allele exclusions against current ST10042 profiles.

PrepExternalSchema may remove external-schema sequences that are not valid complete
CDSs under chewBBACA rules. That does NOT mean the PubMLST allele definition is
biologically invalid. This audit asks the practical question we care about here:
do any current ST10042 official cgST profiles use one of those excluded alleles?
"""

import argparse
import csv
import re
from pathlib import Path

import pandas as pd

ALLELE_RE=re.compile(r"^(?P<locus>.+)_(?P<allele>[^_:]+)(?::|$)")


def parse_invalid(path):
    invalid={}
    for line in Path(path).read_text(encoding="utf-8",errors="replace").splitlines():
        line=line.strip()
        if not line:
            continue
        m=ALLELE_RE.match(line)
        if not m:
            continue
        invalid.setdefault(m.group("locus"),set()).add(m.group("allele"))
    return invalid


def main():
    p=argparse.ArgumentParser()
    p.add_argument(
        "--invalid",
        default="schema/chewbbaca_pubmlst_cgmlst_v2_scheme8_2026-10-04_invalid_alleles.txt",
    )
    p.add_argument(
        "--profiles",
        default="schema/pubmlst_campy_cgmlst_v2_scheme8_2026-10-04/profiles.tsv",
    )
    p.add_argument(
        "--official",
        default="results/pubmlst_st10042_v2_lincodes.tsv",
    )
    p.add_argument(
        "--out",
        default="results/chewbbaca_schema_invalid_alleles_st10042_audit.tsv",
    )
    a=p.parse_args()

    invalid=parse_invalid(a.invalid)
    invalid_total=sum(len(v) for v in invalid.values())
    print(f"Excluded allele IDs parsed: {invalid_total}")
    print(f"Loci containing excluded alleles: {len(invalid)}")

    off=pd.read_csv(a.official,sep="\t",dtype=str).fillna("")
    off=off[off["cgST_v2"].ne("")].copy()
    target_cgsts=set(off["cgST_v2"].astype(str))
    print(f"Current ST10042 records with official cgST: {len(off)}")
    print(f"Distinct ST10042 cgSTs to audit: {len(target_cgsts)}")

    rows=[]
    found_profiles=set()
    with open(a.profiles,encoding="utf-8",errors="replace",newline="") as fh:
        rd=csv.DictReader(fh,delimiter="\t")
        if not rd.fieldnames:
            raise SystemExit("ERROR profiles TSV has no header")
        pk=rd.fieldnames[0]
        lin_col="LINcode" if "LINcode" in rd.fieldnames else None
        locus_cols=[x for x in rd.fieldnames if x.startswith("CAMP")]
        print(f"Profile primary key column: {pk}")
        print(f"Profile locus columns: {len(locus_cols)}")
        print(f"LINcode column present: {'Yes' if lin_col else 'No'}")

        for rec in rd:
            cgst=str(rec.get(pk,""))
            if cgst not in target_cgsts:
                continue
            found_profiles.add(cgst)
            hits=[]
            for locus,alleles in invalid.items():
                val=str(rec.get(locus,""))
                if val in alleles:
                    hits.append(f"{locus}:{val}")
            rows.append({
                "cgST_v2":cgst,
                "LINcode_snapshot":rec.get(lin_col,"") if lin_col else "",
                "excluded_alleles_used":len(hits),
                "excluded_allele_ids":";".join(hits),
            })

    prof=pd.DataFrame(rows)
    ann=off[["pubmlst_id","isolate","country","town_or_city","year","source","cgST_v2","LINcode_v2"]]
    out=ann.merge(prof,on="cgST_v2",how="left")
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    affected=out[pd.to_numeric(out["excluded_alleles_used"],errors="coerce").fillna(0)>0]
    print(f"Official ST10042 cgST profiles found in snapshot: {len(found_profiles)}/{len(target_cgsts)}")
    print(f"Current ST10042 records using >=1 chewBBACA-excluded schema allele: {len(affected)}")
    if len(affected):
        print("\nAFFECTED CURRENT ST10042 RECORDS")
        print(
            affected[
                ["pubmlst_id","isolate","cgST_v2","LINcode_v2",
                 "excluded_alleles_used","excluded_allele_ids"]
            ].to_string(index=False)
        )
    else:
        print("No current ST10042 official profile uses any of the 163 excluded alleles.")

    print(f"\nWrote: {a.out}")


if __name__=="__main__":
    main()
