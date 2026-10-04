#!/usr/bin/env python3
"""Validate chewBBACA allele numbering on known Azevedo/PubMLST overlap controls.

Primary validation:
- compare each chewBBACA numeric call to the unique exact allele IDs returned by
  the current PubMLST scheme-8 sequence endpoint for the same filtered assembly.

Secondary validation:
- for overlap records with an official current cgST, compare calls with the frozen
  PubMLST profile row and report its official frozen LINcode.

The goal is not to force every local assembly to equal its stored PubMLST profile:
our reassembly may recover loci that are missing in the official stored assembly.
The critical check is that a known exact PubMLST allele retains the same allele ID
after PrepExternalSchema/AlleleCall.
"""

import argparse
import csv
import json
import re
from pathlib import Path

import pandas as pd

NUMERIC=re.compile(r"^\d+$")


def find_results(path):
    p=Path(path)
    hits=list(p.rglob("results_alleles.tsv"))
    if len(hits)!=1:
        raise SystemExit(f"ERROR expected one results_alleles.tsv under {p}, found {len(hits)}")
    return hits[0]


def unique_rest_calls(path):
    obj=json.loads(Path(path).read_text(encoding="utf-8"))
    out={}
    for locus,hits in (obj.get("exact_matches",{}) or {}).items():
        ids={str(x.get("allele_id")) for x in (hits or []) if x.get("allele_id") is not None}
        if len(ids)==1:
            out[locus]=next(iter(ids))
    return out


def load_target_profiles(path,target):
    target=set(str(x) for x in target if str(x))
    rows={}
    if not target:
        return rows,None,None
    with open(path,encoding="utf-8",errors="replace",newline="") as fh:
        rd=csv.DictReader(fh,delimiter="\t")
        if not rd.fieldnames:
            raise SystemExit("ERROR frozen profiles TSV has no header")
        pk=rd.fieldnames[0]
        lin="LINcode" if "LINcode" in rd.fieldnames else None
        for rec in rd:
            key=str(rec.get(pk,""))
            if key in target:
                rows[key]=rec
    return rows,pk,lin


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--pilot-dir",default="results/chewbbaca_overlap_pilot")
    p.add_argument("--meta",default="results/chewbbaca_overlap_inputs.tsv")
    p.add_argument("--official",default="results/pubmlst_st10042_v2_lincodes.tsv")
    p.add_argument(
        "--profiles",
        default="schema/pubmlst_campy_cgmlst_v2_scheme8_2026-10-04/profiles.tsv",
    )
    p.add_argument(
        "--rest-cache",
        default="results/pubmlst_v2_innuca_like_trimmed_raw",
    )
    p.add_argument(
        "--invalid",
        default="schema/chewbbaca_pubmlst_cgmlst_v2_scheme8_2026-10-04_invalid_alleles.txt",
    )
    p.add_argument("--out",default="results/chewbbaca_overlap_validation.tsv")
    a=p.parse_args()

    calls_path=find_results(a.pilot_dir)
    calls=pd.read_csv(calls_path,sep="\t",dtype=str).fillna("")
    first=calls.columns[0]
    calls=calls.set_index(first)

    meta=pd.read_csv(a.meta,sep="\t",dtype=str).fillna("")
    official=pd.read_csv(a.official,sep="\t",dtype=str).fillna("")
    omap=official.set_index("pubmlst_id").to_dict("index")

    target_cgsts=[]
    for _,r in meta.iterrows():
        pid=str(r["pubmlst_id"]).split(";")[0]
        cgst=omap.get(pid,{}).get("cgST_v2","")
        if cgst:
            target_cgsts.append(cgst)
    profiles,pk,lin_col=load_target_profiles(a.profiles,target_cgsts)

    invalid=set()
    if Path(a.invalid).exists():
        for line in Path(a.invalid).read_text(encoding="utf-8",errors="replace").splitlines():
            token=line.split(":",1)[0].strip()
            if token:
                invalid.add(token)

    rows=[]
    total_exact=0
    total_match=0
    total_mismatch=0

    for _,m in meta.iterrows():
        analysis_id=m["analysis_id"]
        aid=m["azevedo_strain_id"]
        pid=str(m["pubmlst_id"]).split(";")[0]
        if analysis_id not in calls.index:
            raise SystemExit(f"ERROR pilot output missing {analysis_id}")
        rowcalls=calls.loc[analysis_id]

        rest=Path(a.rest_cache)/f"{aid}.innuca_like.sequence.json"
        if not rest.exists():
            raise SystemExit(f"ERROR missing REST validation cache for {aid}: {rest}")
        exact=unique_rest_calls(rest)

        matched=0
        mismatched=[]
        nonnumeric=[]
        excluded_exact=[]
        for locus,allele in exact.items():
            val=str(rowcalls.get(locus,""))
            if val==allele:
                matched+=1
            elif NUMERIC.match(val):
                mismatched.append(f"{locus}:{allele}>{val}")
            else:
                nonnumeric.append(f"{locus}:{allele}>{val}")
            if f"{locus}_{allele}" in invalid:
                excluded_exact.append(f"{locus}:{allele}")

        total_exact += len(exact)
        total_match += matched
        total_mismatch += len(mismatched)

        off=omap.get(pid,{})
        cgst=off.get("cgST_v2","")
        prof=profiles.get(cgst,{}) if cgst else {}
        profile_compared=0
        profile_match=0
        profile_mismatch=0
        profile_special=0
        if prof:
            for locus in (x for x in calls.columns if x.startswith("CAMP")):
                expected=str(prof.get(locus,""))
                if expected in ("","N","0"):
                    continue
                profile_compared+=1
                val=str(rowcalls.get(locus,""))
                if val==expected:
                    profile_match+=1
                elif NUMERIC.match(val):
                    profile_mismatch+=1
                else:
                    profile_special+=1

        rows.append({
            "analysis_id":analysis_id,
            "azevedo_strain_id":aid,
            "pubmlst_id":pid,
            "published_cluster":m.get("published_figure4_cluster",""),
            "published_cluster21":m.get("published_cluster21",""),
            "rest_unique_exact_loci":len(exact),
            "chew_matches_rest_exact_id":matched,
            "chew_numeric_mismatches_to_rest":len(mismatched),
            "chew_special_calls_at_rest_exact_loci":len(nonnumeric),
            "rest_exact_alleles_removed_by_prep_schema":len(excluded_exact),
            "rest_mismatch_details":";".join(mismatched[:50]),
            "rest_special_details":";".join(nonnumeric[:50]),
            "official_cgST_v2":cgst,
            "official_LINcode_v2":off.get("LINcode_v2",""),
            "snapshot_LINcode":prof.get(lin_col,"") if prof and lin_col else "",
            "official_profile_loci_compared":profile_compared,
            "chew_matches_official_profile":profile_match,
            "chew_numeric_mismatches_official_profile":profile_mismatch,
            "chew_special_calls_official_profile":profile_special,
        })

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    print(f"chewBBACA results: {calls_path}")
    print(f"Controls validated: {len(out)}")
    print(f"Unique PubMLST REST exact calls compared: {total_exact}")
    print(f"Exact allele IDs reproduced by chewBBACA: {total_match}")
    print(f"Numeric allele-ID mismatches: {total_mismatch}")
    print()

    cols=[
        "analysis_id","pubmlst_id","published_cluster",
        "rest_unique_exact_loci","chew_matches_rest_exact_id",
        "chew_numeric_mismatches_to_rest","chew_special_calls_at_rest_exact_loci",
        "rest_exact_alleles_removed_by_prep_schema",
        "official_cgST_v2","official_LINcode_v2","snapshot_LINcode",
    ]
    print(out[cols].to_string(index=False))

    bad=out[
        (pd.to_numeric(out["chew_numeric_mismatches_to_rest"],errors="coerce")>0)
        | (pd.to_numeric(out["rest_exact_alleles_removed_by_prep_schema"],errors="coerce")>0)
    ]
    print()
    if len(bad):
        print("VALIDATION WARNING")
        print("At least one known PubMLST exact allele was not preserved cleanly.")
        print("Do NOT scale to 142 genomes until investigated.")
    else:
        print("ALLELE-ID VALIDATION PASSED")
        print("All known PubMLST exact allele IDs were preserved by chewBBACA.")
        print("Special calls at other loci may represent genuine novel/missing/problematic loci.")

    print(f"\nWrote: {a.out}")


if __name__=="__main__":
    main()
