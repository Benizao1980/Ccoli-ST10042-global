#!/usr/bin/env python3
"""Summarise current public ENA availability of the Azevedo ST10042 dataset."""

import argparse
from pathlib import Path
import pandas as pd


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--azevedo", default="data/europe_217_manifest.tsv")
    p.add_argument("--ena", default="data/europe_ena_fastq_manifest.tsv")
    p.add_argument("--out", default="results/europe_ena_availability_by_isolate.tsv")
    a = p.parse_args()

    az = pd.read_csv(a.azevedo, sep="\t", dtype=str).fillna("")
    en = pd.read_csv(a.ena, sep="\t", dtype=str).fillna("")

    keep = [
        "Strain_ID", "Bioproject", "Run Accession Number", "Country",
        "Isolation year", "Source", "Included in Figure 4",
        "Figure 4 cluster clean", "Cluster 21 member"
    ]
    x = az[keep].merge(
        en[["Strain_ID","status","fastq_ftp"]],
        on="Strain_ID",
        how="left",
        validate="one_to_one",
    )
    x["ena_public_fastq"] = x["status"].eq("resolved")
    x["has_fastq_url"] = x["fastq_ftp"].ne("")

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    x.to_csv(a.out, sep="\t", index=False)

    print("OVERALL")
    print(x["ena_public_fastq"].value_counts().rename(index={True:"resolved",False:"unresolved"}).to_string())

    print("\nBY BIOPROJECT")
    tab = pd.crosstab(x["Bioproject"], x["ena_public_fastq"])
    tab = tab.rename(columns={False:"unresolved", True:"resolved"})
    for c in ["resolved","unresolved"]:
        if c not in tab:
            tab[c] = 0
    tab["total"] = tab["resolved"] + tab["unresolved"]
    tab = tab[["resolved","unresolved","total"]].sort_values(["unresolved","total"], ascending=False)
    print(tab.to_string())

    print("\nBY COUNTRY")
    tab = pd.crosstab(x["Country"], x["ena_public_fastq"]).rename(columns={False:"unresolved",True:"resolved"})
    for c in ["resolved","unresolved"]:
        if c not in tab:
            tab[c] = 0
    tab["total"] = tab["resolved"] + tab["unresolved"]
    print(tab[["resolved","unresolved","total"]].sort_values("total",ascending=False).to_string())

    print("\nPUBLISHED FIGURE 4")
    f4 = x[x["Included in Figure 4"].eq("Yes")]
    print(f"Figure 4 isolates: {len(f4)}")
    print(f"Figure 4 with public FASTQ: {f4['ena_public_fastq'].sum()}")
    print(f"Figure 4 without public FASTQ: {(~f4['ena_public_fastq']).sum()}")

    print("\nPUBLISHED CLUSTER 21")
    c21 = x[x["Cluster 21 member"].eq("Yes")]
    print(f"Cluster 21 isolates: {len(c21)}")
    print(f"Cluster 21 with public FASTQ: {c21['ena_public_fastq'].sum()}")
    print(f"Cluster 21 without public FASTQ: {(~c21['ena_public_fastq']).sum()}")
    print("\nCluster 21 availability by country")
    print(pd.crosstab(c21["Country"], c21["ena_public_fastq"]).rename(columns={False:"unresolved",True:"resolved"}).to_string())

    print("\nUNRESOLVED BIOPROJECT/RUN EXAMPLES")
    cols = ["Strain_ID","Country","Bioproject","Run Accession Number","Figure 4 cluster clean"]
    print(x.loc[~x["ena_public_fastq"], cols].head(30).to_string(index=False))

    print(f"\nWrote: {a.out}")


if __name__ == "__main__":
    main()
