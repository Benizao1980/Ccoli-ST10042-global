#!/usr/bin/env python3
"""Compare local SPAdes assemblies with official PubMLST contigs on overlap isolates.

This diagnoses whether multiple exact allele hits seen in the local assemblies are
introduced by our reassembly or are also present in the assembly currently stored
for the same isolate in PubMLST.

For each overlap with both local and PubMLST FASTA:
  - query scheme 8 on the full local SPAdes contigs;
  - query scheme 8 after filtering local contigs to >=500 bp;
  - query scheme 8 on the official PubMLST contigs;
  - compare exact-hit, ambiguous and no-exact-hit locus counts.

The API is queried sequentially and responses are cached under results/.
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


def get_json(url,retries=4):
    last=None
    for i in range(1,retries+1):
        try:
            r=requests.get(url,timeout=120)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            last=e
            if i<retries:
                time.sleep(2**(i-1))
    raise RuntimeError(f"GET failed {url}: {last}")


def post_json(url,payload,retries=3):
    last=None
    for i in range(1,retries+1):
        try:
            r=requests.post(url,json=payload,timeout=900)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            last=e
            if i<retries:
                time.sleep(2**i)
    raise RuntimeError(f"POST failed {url}: {last}")


def scheme_loci():
    obj=get_json(f"{BASE}/schemes/{SCHEME_ID}/loci")
    return {str(u).rstrip("/").split("/")[-1] for u in (obj.get("loci",[]) or [])}


def b64_bytes(data):
    return base64.b64encode(data).decode("ascii")


def filter_fasta_bytes(path,min_len=500):
    records=[]
    header=None
    seq=[]
    with open(path,encoding="utf-8",errors="replace") as fh:
        for line in fh:
            line=line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    s="".join(seq)
                    if len(s)>=min_len:
                        records.append((header,s))
                header=line
                seq=[]
            else:
                seq.append(line.strip())
        if header is not None:
            s="".join(seq)
            if len(s)>=min_len:
                records.append((header,s))
    text="".join(h+"\n"+s+"\n" for h,s in records)
    return text.encode("utf-8"),len(records)


def query_fasta_bytes(data,cache):
    cache=Path(cache)
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    res=post_json(
        f"{BASE}/schemes/{SCHEME_ID}/sequence",
        {"base64":True,"details":False,"sequence":b64_bytes(data)},
    )
    cache.parent.mkdir(parents=True,exist_ok=True)
    cache.write_text(json.dumps(res,indent=2),encoding="utf-8")
    return res


def summarise(res,loci):
    hits=res.get("exact_matches",{}) or {}
    any_loci=set(hits)
    ambiguous=0
    unique=0
    extra_hits=0
    for locus,vals in hits.items():
        ids={str(v.get("allele_id")) for v in (vals or []) if v.get("allele_id") is not None}
        if len(ids)==1:
            unique += 1
        elif len(ids)>1:
            ambiguous += 1
            extra_hits += len(ids)-1
    return {
        "exact_loci_any":len(any_loci),
        "unique_exact_loci":unique,
        "ambiguous_loci":ambiguous,
        "extra_allele_hits":extra_hits,
        "loci_without_exact_hit":len(loci-any_loci),
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--combined",default="results/combined_st10042_provisional.tsv")
    p.add_argument("--assembly-dir",default="data/europe_assemblies")
    p.add_argument("--pubmlst-dir",default="data/pubmlst_st10042_assemblies")
    p.add_argument("--cache-dir",default="results/pubmlst_v2_assembly_compare_raw")
    p.add_argument("--out",default="results/pubmlst_v2_assembly_compare.tsv")
    a=p.parse_args()

    loci=scheme_loci()
    print(f"PubMLST scheme {SCHEME_ID} loci: {len(loci)}")
    if len(loci)!=EXPECTED_LOCI:
        print(f"WARNING expected {EXPECTED_LOCI}")

    c=pd.read_csv(a.combined,sep="\t",dtype=str).fillna("")
    q=c[c["provenance"].eq("Azevedo+PubMLST")].copy()

    rows=[]
    for _,r in q.iterrows():
        aid=r["azevedo_strain_id"]
        pid=str(r["pubmlst_id"]).split(";")[0]
        local=Path(a.assembly_dir)/aid/"contigs.fasta"
        official=Path(a.pubmlst_dir)/f"{pid}.fasta"
        if not local.exists() or not official.exists():
            continue

        print(f"\n{aid} / PubMLST {pid}")
        local_bytes=local.read_bytes()
        local500,n500=filter_fasta_bytes(local,500)
        official_bytes=official.read_bytes()

        variants=[
            ("local_all",local_bytes),
            ("local_ge500",local500),
            ("pubmlst_official",official_bytes),
        ]
        for label,data in variants:
            cache=Path(a.cache_dir)/f"{aid}.{label}.json"
            res=query_fasta_bytes(data,cache)
            s=summarise(res,loci)
            rows.append({
                "azevedo_strain_id":aid,
                "pubmlst_id":pid,
                "published_cluster":r["published_figure4_cluster"],
                "published_cluster21":r["published_cluster21"],
                "assembly_variant":label,
                "fasta_bytes":len(data),
                "local_contigs_ge500":n500 if label=="local_ge500" else "",
                **s,
            })
            print(
                f"  {label:16s} exact={s['exact_loci_any']:4d} "
                f"unique={s['unique_exact_loci']:4d} "
                f"ambiguous={s['ambiguous_loci']:4d} "
                f"no_exact={s['loci_without_exact_hit']:3d}"
            )

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    print("\nCOMPARISON")
    if len(out):
        pivot=out.pivot(
            index=["azevedo_strain_id","pubmlst_id","published_cluster"],
            columns="assembly_variant",
            values=["exact_loci_any","ambiguous_loci","loci_without_exact_hit"],
        )
        print(pivot.to_string())
    else:
        print("No overlaps had both local and official PubMLST FASTA.")

    print(f"\nWrote: {a.out}")
    print("Interpretation: if local_ge500 approaches the official contig profile, small "
          "SPAdes contigs are driving ambiguity. If the official FASTA is similarly "
          "ambiguous, multiple allele hits are intrinsic to the stored isolate/assembly "
          "and should not be treated as automatic cgMLST failure.")


if __name__=="__main__":
    main()
