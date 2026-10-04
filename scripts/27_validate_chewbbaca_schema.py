#!/usr/bin/env python3
"""Validate the adapted chewBBACA copy of the frozen PubMLST cgMLST-v2 schema."""

import argparse
from pathlib import Path
import pandas as pd

EXPECTED=1142


def nonempty_lines(path):
    if not path.exists():
        return []
    return [x.strip() for x in path.read_text(encoding="utf-8",errors="replace").splitlines() if x.strip()]


def main():
    p=argparse.ArgumentParser()
    p.add_argument(
        "--schema",
        default="schema/chewbbaca_pubmlst_cgmlst_v2_scheme8_2026-10-04",
    )
    a=p.parse_args()

    schema=Path(a.schema)
    if not schema.exists():
        raise SystemExit(f"ERROR schema directory not found: {schema}")

    loci=sorted(x for x in schema.glob("*.fasta") if x.is_file())
    short=sorted((schema/"short").glob("*.fasta")) if (schema/"short").exists() else []

    parent=schema.parent
    stem=schema.name
    invalid_loci=nonempty_lines(parent/f"{stem}_invalid_loci.txt")
    invalid_alleles=nonempty_lines(parent/f"{stem}_invalid_alleles.txt")
    stats_path=parent/f"{stem}_summary_stats.tsv"

    print(f"Adapted schema loci: {len(loci)}")
    print(f"Representative-locus FASTAs: {len(short)}")
    print(f"Invalid loci: {len(invalid_loci)}")
    print(f"Invalid alleles: {len(invalid_alleles)}")

    if stats_path.exists():
        stats=pd.read_csv(stats_path,sep="\t",dtype=str).fillna("")
        print(f"Summary-stat rows: {len(stats)}")
        print(f"Summary stats: {stats_path}")
    else:
        raise SystemExit(f"ERROR missing PrepExternalSchema summary: {stats_path}")

    if invalid_loci:
        print("\nINVALID LOCI")
        print("\n".join(invalid_loci[:30]))

    if len(loci)!=EXPECTED:
        raise SystemExit(
            f"ERROR adapted schema has {len(loci)} loci; expected {EXPECTED}. "
            "Do not proceed to allele calling until investigated."
        )
    if len(short)!=EXPECTED:
        raise SystemExit(
            f"ERROR short/ representative directory has {len(short)} loci; expected {EXPECTED}."
        )

    print("\nSCHEMA VALIDATION PASSED")
    print("All 1,142 PubMLST cgMLST-v2 loci survived chewBBACA adaptation.")


if __name__=="__main__":
    main()
