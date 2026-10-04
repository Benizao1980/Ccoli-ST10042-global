#!/usr/bin/env python3
"""Sequential PubMLST cgMLST-v2 exact-hit QC for all filtered Azevedo genomes.

This is deliberately NOT an array job. BIGSdb requests are sent sequentially and
cached locally. The scheme-8 sequence endpoint reports exact matches to currently
known PubMLST alleles. Therefore loci_without_exact_hit is a useful QC diagnostic,
but it is not yet the final missing-locus count because a genuine novel allele can
also lack an exact match.
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


def post_json(url,payload,retries=4):
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


def query_fasta(path,cache):
    cache=Path(cache)
    if cache.exists():
        return json.loads(cache.read_text(encoding="utf-8"))
    payload={
        "base64":True,
        "details":False,
        "sequence":base64.b64encode(Path(path).read_bytes()).decode("ascii"),
    }
    res=post_json(f"{BASE}/schemes/{SCHEME_ID}/sequence",payload)
    cache.parent.mkdir(parents=True,exist_ok=True)
    cache.write_text(json.dumps(res,indent=2),encoding="utf-8")
    return res


def summarise(res,loci):
    hits=res.get("exact_matches",{}) or {}
    any_loci=set(hits)
    unique=0
    ambiguous=0
    extra=0
    for locus,vals in hits.items():
        ids={str(v.get("allele_id")) for v in (vals or []) if v.get("allele_id") is not None}
        if len(ids)==1:
            unique += 1
        elif len(ids)>1:
            ambiguous += 1
            extra += len(ids)-1
    no_exact=len(loci-any_loci)
    return {
        "exact_loci_any":len(any_loci),
        "unique_exact_loci":unique,
        "ambiguous_loci":ambiguous,
        "extra_allele_hits":extra,
        "loci_without_exact_hit":no_exact,
        # Diagnostic only: novel alleles can also be counted here.
        "exact_unmatched_le25":"Yes" if no_exact<=25 else "No",
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",default="results/azevedo_spades_manifest.tsv")
    p.add_argument("--published",default="data/europe_217_manifest.tsv")
    p.add_argument("--qc",default="results/azevedo_innuca_like_trimmed_qc.tsv")
    p.add_argument("--filtered-dir",default="data/europe_assemblies_innuca_like_trimmed")
    p.add_argument("--cache-dir",default="results/pubmlst_v2_all_azevedo_raw")
    p.add_argument("--out",default="results/pubmlst_v2_all_azevedo_exact_qc.tsv")
    p.add_argument("--sleep",type=float,default=0.25)
    a=p.parse_args()

    loci=scheme_loci()
    print(f"PubMLST scheme {SCHEME_ID} loci: {len(loci)}")
    if len(loci)!=EXPECTED_LOCI:
        raise SystemExit(f"ERROR expected {EXPECTED_LOCI} loci")

    m=pd.read_csv(a.manifest,sep="\t",dtype=str).fillna("")
    pub=pd.read_csv(a.published,sep="\t",dtype=str).fillna("")
    pub=pub.set_index("Strain_ID")
    qc=pd.read_csv(a.qc,sep="\t",dtype=str).fillna("").set_index("Strain_ID")

    rows=[]
    for i,r in m.iterrows():
        sample=r["Strain_ID"]
        fasta=Path(a.filtered_dir)/sample/"contigs.innuca_like.fasta"
        if not fasta.exists() or fasta.stat().st_size==0:
            raise SystemExit(f"ERROR missing filtered FASTA for {sample}: {fasta}")

        print(f"[{i+1}/{len(m)}] {sample}",flush=True)
        res=query_fasta(
            fasta,
            Path(a.cache_dir)/f"{sample}.scheme8.sequence.json"
        )
        ss=summarise(res,loci)

        pmeta=pub.loc[sample] if sample in pub.index else pd.Series(dtype=str)
        qmeta=qc.loc[sample] if sample in qc.index else pd.Series(dtype=str)

        row={
            "array_index":r.get("array_index",""),
            "Strain_ID":sample,
            "run_accession":r.get("run_accession",""),
            "Country":pmeta.get("Country",r.get("Country","")),
            "Isolation_year":pmeta.get("Isolation year",""),
            "Source_Group":pmeta.get("Source Group clean",pmeta.get("Source Group","")),
            "Figure4_cluster":pmeta.get("Figure 4 cluster clean",""),
            "Cluster21_member":pmeta.get("Cluster 21 member",""),
            "filtered_contigs":qmeta.get("filtered_contigs",""),
            "filtered_bp":qmeta.get("filtered_bp",""),
            **ss,
            "filtered_fasta":str(fasta),
        }
        rows.append(row)

        print(
            f"  exact={ss['exact_loci_any']}/{len(loci)} "
            f"unique={ss['unique_exact_loci']} "
            f"ambiguous={ss['ambiguous_loci']} "
            f"no_exact={ss['loci_without_exact_hit']}",
            flush=True
        )
        time.sleep(a.sleep)

    out=pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(a.out,sep="\t",index=False)

    no=pd.to_numeric(out["loci_without_exact_hit"],errors="coerce")
    amb=pd.to_numeric(out["ambiguous_loci"],errors="coerce")

    print("\nSUMMARY")
    print(f"Genomes queried: {len(out)}")
    print(f"Exact-unmatched <=25 (diagnostic, not final allele-call QC): {(no<=25).sum()}")
    print(f"Exact-unmatched >25: {(no>25).sum()}")
    print(f"Zero ambiguous exact-hit loci: {(amb==0).sum()}")
    print(f"One or more ambiguous exact-hit loci: {(amb>0).sum()}")

    print("\nHighest loci_without_exact_hit")
    cols=[
        "Strain_ID","Country","Figure4_cluster","Cluster21_member",
        "filtered_contigs","filtered_bp","exact_loci_any",
        "unique_exact_loci","ambiguous_loci","loci_without_exact_hit"
    ]
    print(
        out.assign(_no=no)
           .sort_values(["_no","ambiguous_loci"],ascending=[False,False])
           .head(20)[cols].to_string(index=False)
    )

    print("\nPublished cluster 21")
    c21=out[out["Cluster21_member"].str.lower().eq("yes")]
    if len(c21):
        print(c21[cols].sort_values("loci_without_exact_hit",key=lambda x:pd.to_numeric(x,errors="coerce")).to_string(index=False))
    else:
        print("No accessible cluster-21 rows found.")

    print(f"\nWrote: {a.out}")
    print("NOTE: a locus without an exact known-allele hit may contain a genuine novel allele;")
    print("this table is the full-cohort exact-hit QC checkpoint, not the final allele matrix.")


if __name__=="__main__":
    main()
