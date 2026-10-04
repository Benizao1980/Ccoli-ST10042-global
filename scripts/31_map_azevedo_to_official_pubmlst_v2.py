#!/usr/bin/env python3
"""Map filtered Azevedo genomes onto *official* current PubMLST cgMLST-v2 profiles.

Primary principle:
- PubMLST/BIGSdb remains the nomenclature authority.
- Reuse the cached scheme-8 whole-genome query responses generated in script 24.
- Where BIGSdb identifies one existing cgST, retrieve that profile directly from
  PubMLST and export its official LINcode and Cjc_cgc2_* classifications.
- Keep genomes with no exact existing profile, or multiple profile matches, separate
  for nearest-profile / Genome Comparator follow-up. Do not invent cgST/LIN values.

This script does not perform local clustering.
"""

import argparse
import json
import re
import time
from pathlib import Path

import pandas as pd
import requests

BASE="https://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef"
SCHEME_ID=8
CGC=[200,100,50,25,10,5]


def get_json(url,cache,retries=4):
    cache=Path(cache)
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    last=None
    for i in range(1,retries+1):
        try:
            r=requests.get(url,timeout=180)
            r.raise_for_status()
            obj=r.json()
            cache.parent.mkdir(parents=True,exist_ok=True)
            cache.write_text(json.dumps(obj,indent=2)+"\n",encoding="utf-8")
            return obj
        except requests.RequestException as e:
            last=e
            if i<retries:
                time.sleep(2**(i-1))
    raise RuntimeError(f"GET failed {url}: {last}")


def split_profile_ids(value):
    if value is None:
        return []
    if isinstance(value,list):
        vals=value
    else:
        vals=re.split(r"[,;\s]+",str(value).strip())
    out=[]
    for x in vals:
        x=str(x).strip()
        if x:
            out.append(x)
    return sorted(set(out),key=lambda x:(not x.isdigit(),int(x) if x.isdigit() else x))


def cgc_group(profile,threshold):
    cs=profile.get("classification_schemes",{}) or {}
    key=f"Cjc_cgc2_{threshold}"
    obj=cs.get(key,{}) or {}
    group=(obj.get("group",{}) or {}).get("group")
    if group is not None:
        return str(group)

    # Defensive fallback if display names ever change slightly.
    for name,obj in cs.items():
        if str(name).lower().replace("-","_")==key.lower():
            group=((obj or {}).get("group",{}) or {}).get("group")
            return "" if group is None else str(group)
    return ""


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--qc",default="results/pubmlst_v2_all_azevedo_exact_qc.tsv")
    p.add_argument("--sequence-cache",default="results/pubmlst_v2_all_azevedo_raw")
    p.add_argument("--profile-cache",default="results/pubmlst_v2_official_profiles")
    p.add_argument("--official-current",default="results/pubmlst_st10042_v2_lincodes.tsv")
    p.add_argument("--out",default="results/azevedo_pubmlst_official_v2.tsv")
    p.add_argument("--sleep",type=float,default=0.15)
    a=p.parse_args()

    qc=pd.read_csv(a.qc,sep="\t",dtype=str).fillna("")
    current=pd.read_csv(a.official_current,sep="\t",dtype=str).fillna("")
    current_map=current.set_index("pubmlst_id").to_dict("index")

    rows=[]
    requested={}
    for i,r in qc.iterrows():
        sid=r["Strain_ID"]
        cache=Path(a.sequence_cache)/f"{sid}.scheme8.sequence.json"
        if not cache.exists():
            raise SystemExit(f"ERROR missing cached scheme-8 response for {sid}: {cache}")
        seq=json.loads(cache.read_text(encoding="utf-8"))

        fields=seq.get("fields",{}) or {}
        cgst_value=fields.get("cgST","")
        candidates=split_profile_ids(cgst_value)

        row={
            "Strain_ID":sid,
            "Country":r.get("Country",""),
            "Figure4_cluster":r.get("Figure4_cluster",""),
            "Cluster21_member":r.get("Cluster21_member",""),
            "pubmlst_overlap_id":r.get("pubmlst_id",""),
            "exact_loci_any":r.get("exact_loci_any",""),
            "ambiguous_loci":r.get("ambiguous_loci",""),
            "loci_without_exact_hit":r.get("loci_without_exact_hit",""),
            "pubmlst_cgST_candidates":";".join(candidates),
            "pubmlst_profile_match_status":(
                "single_official_profile" if len(candidates)==1 else
                "multiple_official_profiles" if len(candidates)>1 else
                "no_exact_official_profile"
            ),
            "official_cgST":"",
            "official_LINcode":"",
        }
        for t in CGC:
            row[f"Cjc_cgc2_{t}"]=""

        if len(candidates)==1:
            cgst=candidates[0]
            if cgst not in requested:
                url=f"{BASE}/schemes/{SCHEME_ID}/profiles/{cgst}?allele_ids_only=1"
                requested[cgst]=get_json(
                    url,
                    Path(a.profile_cache)/f"cgST_{cgst}.json"
                )
                time.sleep(a.sleep)
            prof=requested[cgst]
            row["official_cgST"]=str(prof.get("cgST",cgst))
            row["official_LINcode"]=str(prof.get("LINcode",""))
            for t in CGC:
                row[f"Cjc_cgc2_{t}"]=cgc_group(prof,t)

        # Internal check for Azevedo records that already have a known PubMLST ID.
        pid=str(r.get("pubmlst_id","")).split(";")[0].strip()
        off=current_map.get(pid,{}) if pid else {}
        row["known_overlap_official_cgST"]=off.get("cgST_v2","")
        row["known_overlap_official_LINcode"]=off.get("LINcode_v2","")
        if pid and row["official_cgST"] and off.get("cgST_v2",""):
            row["overlap_cgST_agrees"]="Yes" if row["official_cgST"]==off.get("cgST_v2","") else "No"
        else:
            row["overlap_cgST_agrees"]=""

        rows.append(row)
        print(
            f"[{i+1}/{len(qc)}] {sid}: "
            f"{row['pubmlst_profile_match_status']} "
            f"{row['official_cgST'] or row['pubmlst_cgST_candidates'] or '-'}"
        )

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    print("\nSUMMARY")
    print(f"Azevedo genomes: {len(out)}")
    print(out["pubmlst_profile_match_status"].value_counts().to_string())
    print(f"Official LINcodes recovered: {out['official_LINcode'].ne('').sum()}")

    ov=out[out["known_overlap_official_cgST"].ne("")]
    if len(ov):
        print("\nKNOWN PUBMLST OVERLAP VALIDATION")
        print(
            ov[[
                "Strain_ID","pubmlst_overlap_id",
                "official_cgST","known_overlap_official_cgST",
                "overlap_cgST_agrees",
                "official_LINcode","known_overlap_official_LINcode"
            ]].to_string(index=False)
        )
        disagree=ov[ov["overlap_cgST_agrees"].eq("No")]
        if len(disagree):
            print("\nWARNING: at least one exact profile assignment disagrees with the known PubMLST overlap.")
            print("Do not interpret those assignments until investigated.")

    print("\nPUBLISHED CLUSTER 21")
    c21=out[out["Cluster21_member"].str.lower().eq("yes")]
    cols=[
        "Strain_ID","official_cgST","official_LINcode",
        "Cjc_cgc2_25","Cjc_cgc2_10","Cjc_cgc2_5",
        "pubmlst_profile_match_status"
    ]
    print(c21[cols].to_string(index=False) if len(c21) else "None")

    print(f"\nWrote: {a.out}")
    print("Only values returned for an existing PubMLST profile are labelled official.")
    print("Unmatched genomes remain for PubMLST nearest-profile/Genome Comparator follow-up.")


if __name__=="__main__":
    main()
