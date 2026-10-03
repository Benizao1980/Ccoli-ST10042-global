#!/usr/bin/env python3
import argparse, csv, hashlib
from pathlib import Path
import requests

ENA="https://www.ebi.ac.uk/ena/portal/api/filereport"

def md5sum(path):
    h=hashlib.md5()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",default="data/europe_217_manifest.tsv")
    p.add_argument("--resolved",default="data/europe_ena_fastq_manifest.tsv")
    p.add_argument("--download",action="store_true")
    p.add_argument("--outdir",default="data/europe_reads")
    a=p.parse_args()
    with open(a.manifest,encoding="utf-8") as fh:
        rows=list(csv.DictReader(fh,delimiter="\t"))
    resolved=[]
    for n,row in enumerate(rows,1):
        run=row["Run Accession Number"].strip()
        params={"accession":run,"result":"read_run",
                "fields":"run_accession,fastq_ftp,fastq_md5,fastq_bytes","format":"tsv"}
        r=requests.get(ENA,params=params,timeout=60); r.raise_for_status()
        lines=[x for x in r.text.strip().splitlines() if x]
        if len(lines)<2:
            print("WARNING no ENA FASTQ:",run); continue
        rec=dict(zip(lines[0].split("\t"),lines[1].split("\t")))
        rec.update({"Strain_ID":row["Strain_ID"],"Country":row["Country"],
                    "Figure4_cluster":row["Figure 4 cluster clean"]})
        resolved.append(rec)
        if a.download:
            sd=Path(a.outdir)/row["Strain_ID"]; sd.mkdir(parents=True,exist_ok=True)
            ftps=[x for x in rec.get("fastq_ftp","").split(";") if x]
            md5s=[x for x in rec.get("fastq_md5","").split(";") if x]
            for j,ftp in enumerate(ftps):
                url="https://"+ftp if not ftp.startswith("http") else ftp
                dest=sd/Path(ftp).name
                if not dest.exists():
                    with requests.get(url,stream=True,timeout=600) as rr:
                        rr.raise_for_status()
                        with open(dest,"wb") as out:
                            for chunk in rr.iter_content(1024*1024):
                                if chunk: out.write(chunk)
                if j<len(md5s) and md5sum(dest)!=md5s[j]:
                    raise RuntimeError(f"MD5 mismatch: {dest}")
        if n%25==0: print(f"Resolved {n}/{len(rows)}")
    fields=["Strain_ID","Country","Figure4_cluster","run_accession","fastq_ftp","fastq_md5","fastq_bytes"]
    Path(a.resolved).parent.mkdir(parents=True,exist_ok=True)
    with open(a.resolved,"w",newline="",encoding="utf-8") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",extrasaction="ignore")
        w.writeheader(); w.writerows(resolved)
    print(f"Wrote {len(resolved)} resolved runs to {a.resolved}")

if __name__=="__main__":
    main()
