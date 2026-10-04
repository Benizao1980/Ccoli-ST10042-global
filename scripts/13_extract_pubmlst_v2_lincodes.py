#!/usr/bin/env python3
"""Extract official current PubMLST cgMLST-v2/LIN assignments for ST10042.

The isolate JSON stores scheme results as a list, so list position must NOT be
treated as a scheme ID. The cgMLST-v2 object is identified by its native LINcode
and/or Cjc_cgc2_* classification keys.

Outputs a compact table that can be joined to the global ST10042 manifest and
prints the focal Peru hierarchy, including any <=5-AD cgc2 group shared with a
published Azevedo cluster-21 anchor already present in PubMLST.
"""

import argparse
import json
from pathlib import Path

import pandas as pd

THRESHOLDS=[1119,1085,982,914,857,680,445,343,183,86,43,10,7,5,3,2,1,0]
CGC_THRESHOLDS=[200,100,50,25,10,5]


def sval(x):
    if x is None:
        return ""
    if isinstance(x,list):
        return ";".join(str(v) for v in x)
    return str(x)


def pick(d,*names):
    for n in names:
        if n in d and d[n] not in (None,""):
            return d[n]
    return ""


def v2_scheme(rec):
    for s in rec.get("schemes",[]) or []:
        if not isinstance(s,dict):
            continue
        cls=s.get("classification_schemes",{}) or {}
        if "LINcode" in s or any(str(k).startswith("Cjc_cgc2_") for k in cls):
            return s
    return {}


def group_ids(scheme, threshold):
    cls=(scheme.get("classification_schemes",{}) or {}).get(f"Cjc_cgc2_{threshold}",{}) or {}
    groups=cls.get("groups",[]) or []
    return ";".join(str(g.get("group")) for g in groups if g.get("group") is not None)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--json",default="data/pubmlst_st10042_all.json")
    p.add_argument("--combined",default="results/combined_st10042_provisional.tsv")
    p.add_argument("--out",default="results/pubmlst_st10042_v2_lincodes.tsv")
    a=p.parse_args()

    records=json.loads(Path(a.json).read_text(encoding="utf-8"))

    annotations={}
    cp=Path(a.combined)
    if cp.exists():
        c=pd.read_csv(cp,sep="\t",dtype=str).fillna("")
        for _,r in c.iterrows():
            for pid in str(r.get("pubmlst_id","")).split(";"):
                pid=pid.strip()
                if pid:
                    annotations[pid]={
                        "canonical_id":r.get("canonical_id",""),
                        "azevedo_strain_id":r.get("azevedo_strain_id",""),
                        "published_cluster21":r.get("published_cluster21",""),
                        "provenance":r.get("provenance",""),
                        "peru_focal":r.get("peru_focal",""),
                    }

    rows=[]
    for rec in records:
        prov=rec.get("provenance",{}) or {}
        pid=sval(prov.get("id"))
        s=v2_scheme(rec)
        lin=sval(s.get("LINcode"))
        parts=lin.split("_") if lin else []

        row={
            "pubmlst_id":pid,
            "isolate":sval(pick(prov,"isolate","isolate_name","strain","strain_id")),
            "country":sval(pick(prov,"country")),
            "town_or_city":sval(pick(prov,"town_or_city","city")),
            "year":sval(pick(prov,"year")),
            "source":sval(pick(prov,"source")),
            "cgST_v2":sval((s.get("fields",{}) or {}).get("cgST")),
            "LINcode_v2":lin,
            "has_v2_assignment":"Yes" if s else "No",
        }
        for i,t in enumerate(THRESHOLDS,1):
            row[f"LIN_L{i}_{t}AD"]=parts[i-1] if len(parts)>=i else ""
        for t in CGC_THRESHOLDS:
            row[f"Cjc_cgc2_{t}"]=group_ids(s,t) if s else ""

        row.update(annotations.get(pid,{
            "canonical_id":"","azevedo_strain_id":"","published_cluster21":"",
            "provenance":"","peru_focal":""
        }))
        rows.append(row)

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    peru=out[out["country"].eq("Peru")].copy()
    print(f"Current PubMLST ST10042 records: {len(out)}")
    print(f"Records with cgMLST-v2/LIN-related assignment: {(out['has_v2_assignment']=='Yes').sum()}")
    print(f"Peru records: {len(peru)}")
    print(f"Peru records with native LINcode: {peru['LINcode_v2'].ne('').sum()}")

    print("\nPERU cgc2 GROUPS")
    for t in CGC_THRESHOLDS:
        col=f"Cjc_cgc2_{t}"
        print(f"\n{t} AD")
        vc=peru[col].replace("",pd.NA).value_counts(dropna=False)
        print(vc.to_string())

    anchors=out[out["published_cluster21"].eq("Yes")]
    print("\nPUBLISHED CLUSTER-21 PUBMLST ANCHORS")
    if len(anchors):
        print(
            anchors[["pubmlst_id","azevedo_strain_id","LINcode_v2","Cjc_cgc2_5"]]
            .to_string(index=False)
        )
    else:
        print("None annotated in combined manifest")

    shared=set(x for x in anchors["Cjc_cgc2_5"] if x)
    if shared:
        hit=peru[peru["Cjc_cgc2_5"].isin(shared)]
        print("\nPERU SHARING A <=5-AD cgc2 GROUP WITH A PUBLISHED CLUSTER-21 ANCHOR")
        if len(hit):
            print(
                hit[["pubmlst_id","isolate","town_or_city","year","source",
                     "cgST_v2","LINcode_v2","Cjc_cgc2_5"]]
                .sort_values("pubmlst_id")
                .to_string(index=False)
            )
        else:
            print("None")
    else:
        print("\nNo cluster-21 anchor currently has a cgc2_5 assignment.")

    print(f"\nWrote: {a.out}")


if __name__=="__main__":
    main()
