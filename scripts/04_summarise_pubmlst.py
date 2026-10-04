#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", default="data/pubmlst_st10042_all.tsv")
    p.add_argument("--city-curation", default="data/peru_city_curation.tsv")
    p.add_argument("--out", default="results/pubmlst_st10042_curated.tsv")
    a = p.parse_args()

    df = pd.read_csv(a.manifest, sep="\t")
    cur = pd.read_csv(a.city_curation, sep="\t")
    cur["pubmlst_id"] = cur["pubmlst_id"].astype(str)
    df["pubmlst_id"] = df["pubmlst_id"].astype(str)

    df = df.merge(cur, on="pubmlst_id", how="left")
    df["town_or_city_curated"] = df["town_or_city"]
    m = df["curated_town_or_city"].notna()
    df.loc[m, "town_or_city_curated"] = df.loc[m, "curated_town_or_city"]

    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, sep="\t", index=False)

    print("GLOBAL COUNTRY")
    print(df["country"].value_counts(dropna=False).to_string())
    print("\nGLOBAL YEAR")
    print(df["year"].value_counts(dropna=False).sort_index().to_string())
    print("\nGLOBAL SOURCE")
    print(df["source"].value_counts(dropna=False).to_string())

    peru = df[df["country"].eq("Peru")].copy()
    print("\nPERU CURATED CITY")
    print(peru["town_or_city_curated"].value_counts(dropna=False).to_string())
    print("\nPERU YEAR")
    print(peru["year"].value_counts(dropna=False).sort_index().to_string())
    print("\nPERU SOURCE")
    print(peru["source"].value_counts(dropna=False).to_string())

    non_peru = df[~df["country"].eq("Peru")].copy()
    cols = ["pubmlst_id","isolate","country","town_or_city","year","source"]
    print("\nNON-PERU RECORDS")
    print(non_peru[cols].sort_values(["country","year","pubmlst_id"]).to_string(index=False))

    missing_year = df[df["year"].isna()][cols]
    if len(missing_year):
        print("\nMISSING YEAR")
        print(missing_year.to_string(index=False))

    print(f"\nWrote curated manifest: {a.out}")

if __name__ == "__main__":
    main()
