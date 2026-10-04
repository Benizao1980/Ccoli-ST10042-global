#!/usr/bin/env python3
"""Discover the current Campylobacter PubMLST cgMLST and classification schemes.

Do not hard-code the cgMLST v2 scheme ID or LIN classification-scheme ID. BIGSdb
exposes both through its REST API; this script records what PubMLST currently
advertises so the analysis remains auditable if IDs/descriptions change.
"""

import argparse
import json
from pathlib import Path
import requests

SEQDEF = "https://rest.pubmlst.org/db/pubmlst_campylobacter_seqdef"
ISOLATES = "https://rest.pubmlst.org/db/pubmlst_campylobacter_isolates"


def get_json(url):
    r=requests.get(url,timeout=120)
    r.raise_for_status()
    return r.json()


def normalise_links(obj, key):
    vals=obj.get(key, [])
    out=[]
    for x in vals:
        if isinstance(x,str):
            out.append({"href":x})
        elif isinstance(x,dict):
            out.append(x)
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/pubmlst_v2_lin_preflight.json")
    a=p.parse_args()

    schemes=get_json(f"{SEQDEF}/schemes")
    scheme_rows=[]
    for item in normalise_links(schemes,"schemes"):
        href=item.get("scheme") or item.get("href")
        if not href:
            continue
        detail=get_json(href)
        scheme_rows.append({
            "href":href,
            "id":detail.get("id"),
            "name":detail.get("name") or detail.get("description") or item.get("description"),
            "description":detail.get("description") or item.get("description"),
            "loci_count":detail.get("loci_count") or (
                len(detail.get("loci",[])) if isinstance(detail.get("loci"),list) else None
            ),
            "lincodes":detail.get("lincodes"),
            "raw":detail,
        })

    classifications=get_json(f"{ISOLATES}/classification_schemes")
    class_rows=[]
    for item in normalise_links(classifications,"classification_schemes"):
        href=item.get("classification_scheme") or item.get("scheme") or item.get("href")
        if not href:
            continue
        detail=get_json(href)
        class_rows.append({
            "href":href,
            "id":detail.get("id"),
            "name":detail.get("name") or detail.get("description") or item.get("description"),
            "description":detail.get("description") or item.get("description"),
            "scheme_id":detail.get("scheme_id"),
            "inclusion_threshold":detail.get("inclusion_threshold"),
            "relative_threshold":detail.get("relative_threshold"),
            "raw":detail,
        })

    payload={"schemes":scheme_rows,"classification_schemes":class_rows}
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    Path(a.out).write_text(json.dumps(payload,indent=2),encoding="utf-8")

    print("PUBMLST CAMPYLOBACTER SCHEMES")
    for x in scheme_rows:
        text=" ".join(str(x.get(k) or "") for k in ("name","description")).lower()
        flag="  <== cgMLST candidate" if "cgmlst" in text else ""
        print(f"id={x['id']} loci={x['loci_count']}  {x['name']}{flag}")

    v2 = next((x for x in scheme_rows if x["id"] == 8), None)
    print("\nNATIVE LINCODE DEFINITION FOR cgMLST v2")
    if v2 and v2.get("lincodes"):
        print(json.dumps(v2["lincodes"], indent=2))
    else:
        print("No lincodes object returned for scheme 8.")

    print("\nCLASSIFICATION SCHEMES (USEFUL cgc2 VIEWS; NOT THE FULL LINCODE)")
    for x in class_rows:
        print(
            f"id={x['id']} scheme_id={x.get('scheme_id')} "
            f"threshold={x.get('inclusion_threshold')}  {x['name']}"
        )

    print(f"\nWrote: {a.out}")
    print("Primary nomenclature: cgMLST v2 scheme 8 + its native lincodes object. "
          "Cjc_cgc2_200..5 are thresholded classification views, not the complete 18-level LIN code.")


if __name__=="__main__":
    main()
