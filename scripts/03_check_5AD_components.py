#!/usr/bin/env python3
import argparse,csv
from collections import defaultdict,deque
BAD={"","-","INF","LNF","PLOT3","PLOT5","NIPH","NIPHEM","ALM","ASM"}

def dist(a,b):
    d=n=0
    for x,y in zip(a,b):
        if x in BAD or y in BAD: continue
        n+=1
        if x!=y: d+=1
    return d,n

def main():
    p=argparse.ArgumentParser()
    p.add_argument("alleles_tsv")
    p.add_argument("--threshold",type=int,default=5)
    p.add_argument("--out",default="clusters_le5.tsv")
    a=p.parse_args()
    with open(a.alleles_tsv,encoding="utf-8") as fh:
        r=csv.reader(fh,delimiter="\t"); next(r); rows=list(r)
    ids=[r[0] for r in rows]; prof=[r[1:] for r in rows]
    adj=defaultdict(set)
    for i in range(len(rows)):
        for j in range(i+1,len(rows)):
            d,n=dist(prof[i],prof[j])
            if n and d<=a.threshold: adj[i].add(j); adj[j].add(i)
    seen=set(); comps=[]
    for i in range(len(rows)):
        if i in seen: continue
        q=deque([i]); seen.add(i); c=[]
        while q:
            u=q.popleft(); c.append(u)
            for v in adj[u]:
                if v not in seen: seen.add(v); q.append(v)
        comps.append(c)
    comps.sort(key=len,reverse=True)
    mem={}
    for k,c in enumerate(comps,1):
        for i in c: mem[ids[i]]=(f"component_{k}",len(c))
    with open(a.out,"w",newline="",encoding="utf-8") as fh:
        w=csv.writer(fh,delimiter="\t"); w.writerow(["sample","component_le5","component_size"])
        for x in ids: w.writerow([x,*mem[x]])
    print(f"{len(rows)} genomes -> {len(comps)} components at <= {a.threshold} AD")
if __name__=="__main__": main()
