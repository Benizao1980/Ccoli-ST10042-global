#!/usr/bin/env python3
"""Pilot PubMLST cgMLST-v2 allele calling on Azevedo assemblies with known PubMLST overlaps.

Purpose:
1. validate our fastp+SPAdes assemblies against existing PubMLST records;
2. quantify exact scheme-8 allele recovery, missing loci, and ambiguous/multiple calls;
3. inspect whether PubMLST returns cgST/other scheme fields for the queried profile.

The BIGSdb scheme sequence endpoint is exact-match based. Therefore this is a
validation/pilot step, not yet the final solution for genuinely novel alleles.
A locus with multiple exact allele hits is recorded as ambiguous but is NOT
equated with a missing locus or automatic cgMLST failure.
Requests are sequential to avoid unnecessary load on PubMLST.
"""

import argparse
import base64
import json
import time
from pathlib import Path

import pandas as pd
import requests

BASE="https://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef"
SCHEME_ID=8
EXPECTED_LOCI=1142


def get_json(url, retries=4):
    last=None
    for attempt in range(1,retries+1):
        try:
            r=requests.get(url,timeout=120)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            last=e
            if attempt<retries:
                time.sleep(2**(attempt-1))
    raise RuntimeError(f"GET failed {url}: {last}")


def post_json(url,payload,retries=3):
    last=None
    for attempt in range(1,retries+1):
        try:
            r=requests.post(url,json=payload,timeout=900)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            last=e
            if attempt<retries:
                time.sleep(2**attempt)
    raise RuntimeError(f"POST failed {url}: {last}")


def scheme_loci():
    obj=get_json(f"{BASE}/schemes/{SCHEME_ID}/loci")
    urls=obj.get("loci",[]) or []
    loci={str(u).rstrip("/").split("/")[-1] for u in urls}
    return loci


def encode_fasta(path):
    return base64.b64encode(Path(path).read_bytes()).decode("ascii")


def unique_alleles(matches):
    out={}
    for locus,hits in (matches or {}).items():
        ids=sorted({str(h.get("allele_id")) for h in (hits or []) if h.get("allele_id") is not None})
        out[locus]=ids
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--combined",default="results/combined_st10042_provisional.tsv")
    p.add_argument("--official",default="results/pubmlst_st10042_v2_lincodes.tsv")
    p.add_argument("--assembly-dir",default="data/europe_assemblies")
    p.add_argument("--out",default="results/pubmlst_v2_overlap_validation.tsv")
    p.add_argument("--raw-dir",default="results/pubmlst_v2_overlap_raw")
    a=p.parse_args()

    loci=scheme_loci()
    print(f"PubMLST scheme {SCHEME_ID} loci returned: {len(loci)}")
    if len(loci)!=EXPECTED_LOCI:
        print(f"WARNING expected {EXPECTED_LOCI} loci; live API returned {len(loci)}")

    c=pd.read_csv(a.combined,sep="\t",dtype=str).fillna("")
    q=c[c["provenance"].eq("Azevedo+PubMLST")].copy()

    official=pd.read_csv(a.official,sep="\t",dtype=str).fillna("")
    omap=official.set_index("pubmlst_id").to_dict("index")

    rawdir=Path(a.raw_dir)
    rawdir.mkdir(parents=True,exist_ok=True)

    rows=[]
    candidates=[]
    for _,r in q.iterrows():
        aid=r["azevedo_strain_id"]
        fa=Path(a.assembly_dir)/aid/"contigs.fasta"
        if fa.exists() and fa.stat().st_size>0:
            candidates.append((r,fa))

    print(f"Known Azevedo/PubMLST overlaps with local SPAdes assembly: {len(candidates)}")
    print("Requests will run sequentially.")

    for i,(r,fa) in enumerate(candidates,1):
        aid=r["azevedo_strain_id"]
        pid=str(r["pubmlst_id"]).split(";")[0]
        print(f"\n[{i}/{len(candidates)}] {aid} / PubMLST {pid} / {fa}")

        payload={
            "base64":True,
            "details":False,
            "sequence":encode_fasta(fa),
        }
        seqres=post_json(f"{BASE}/schemes/{SCHEME_ID}/sequence",payload)
        (rawdir/f"{aid}.sequence.json").write_text(json.dumps(seqres,indent=2),encoding="utf-8")

        ua=unique_alleles(seqres.get("exact_matches",{}))
        unique={locus:ids[0] for locus,ids in ua.items() if len(ids)==1}
        ambiguous={locus:ids for locus,ids in ua.items() if len(ids)>1}
        exact_any=set(ua)
        known_loci=loci if loci else set(ua)
        missing=sorted(known_loci-exact_any)

        designations={locus:[{"allele":allele}] for locus,allele in unique.items()}
        dres={}
        if designations:
            dres=post_json(
                f"{BASE}/schemes/{SCHEME_ID}/designations",
                {"designations":designations},
            )
            (rawdir/f"{aid}.designations.json").write_text(
                json.dumps(dres,indent=2),encoding="utf-8"
            )

        off=omap.get(pid,{})
        fields=dres.get("fields",{}) if isinstance(dres,dict) else {}

        row={
            "azevedo_strain_id":aid,
            "pubmlst_id":pid,
            "published_cluster":r["published_figure4_cluster"],
            "published_cluster21":r["published_cluster21"],
            "assembly_bp":fa.stat().st_size,
            "scheme_loci":len(known_loci),
            "exact_loci_any":len(exact_any),
            "unique_exact_loci":len(unique),
            "ambiguous_loci":len(ambiguous),
            "missing_exact_loci":len(missing),
            "exact_call_fraction":len(unique)/len(known_loci) if known_loci else "",
            "within_25_loci_without_exact_hit":"Yes" if len(missing)<=25 else "No",
            "all_exact_hit_loci_unambiguous":"Yes" if len(ambiguous)==0 else "No",
            "official_cgST_v2":off.get("cgST_v2",""),
            "official_LINcode_v2":off.get("LINcode_v2",""),
            "official_Cjc_cgc2_5":off.get("Cjc_cgc2_5",""),
            "api_fields_json":json.dumps(fields,sort_keys=True),
            "sequence_response_keys":";".join(sorted(seqres.keys())),
            "designations_response_keys":";".join(sorted(dres.keys())) if isinstance(dres,dict) else "",
            "ambiguous_loci_json":json.dumps(ambiguous,sort_keys=True),
        }
        rows.append(row)

        print(
            f"  exact={len(exact_any)}/{len(known_loci)} "
            f"unique={len(unique)} ambiguous={len(ambiguous)} missing={len(missing)}"
        )
        print(f"  official LIN={row['official_LINcode_v2'] or '<none>'}")
        print(f"  API fields={row['api_fields_json']}")

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    print("\nSUMMARY")
    if len(out):
        cols=[
            "azevedo_strain_id","pubmlst_id","published_cluster",
            "exact_loci_any","unique_exact_loci","ambiguous_loci","missing_exact_loci",
            "within_25_loci_without_exact_hit","all_exact_hit_loci_unambiguous","official_cgST_v2",
            "official_Cjc_cgc2_5","official_LINcode_v2"
        ]
        print(out[cols].to_string(index=False))
    else:
        print("No overlap assemblies available.")

    print(f"\nWrote: {a.out}")
    print(f"Raw API responses: {rawdir}")


if __name__=="__main__":
    main()
