# Authenticated PubMLST ST10042 snapshot — 2026-10-03/04

This document records the first successful authenticated PubMLST/BIGSdb retrieval
for the global *Campylobacter coli* ST10042 project.

## Retrieval checkpoint

Puma job `24088323` completed successfully using OAuth authentication.

- PubMLST ST10042 metadata records: **41**
- records with retrievable contigs: **40**
- records without contigs: **1**
- no-contig record: **PubMLST 151737**
- Peru ST10042 records: **23**
- non-Peru ST10042 records: **18**

This independently explains the original exploratory dataset: the PubMLST query
contains 41 ST10042 records, but only 40 can enter sequence-based analysis.

## Country composition

| Country | n |
|---|---:|
| Peru | 23 |
| UK | 8 |
| UK [Northern Ireland] | 1 |
| Portugal | 4 |
| Luxembourg | 2 |
| Spain | 2 |
| Vietnam | 1 |

The current PubMLST collection is therefore **not** a substitute for the broader
Azevedo et al. European surveillance collection. Several countries represented in
the 217-isolate Azevedo dataset are absent from the current PubMLST ST10042 query.

## Year composition

| Year | n |
|---|---:|
| 2019 | 1 |
| 2021 | 6 |
| 2022 | 4 |
| 2023 | 9 |
| 2024 | 7 |
| 2025 | 13 |
| missing | 1 |

## Source composition

| Source | n |
|---|---:|
| human stool | 22 |
| chicken offal or meat | 10 |
| goat | 3 |
| chicken | 2 |
| cattle faeces | 2 |
| environmental waters | 1 |
| beef offal or meat | 1 |

## Peru subset

The 23 Peru records comprise:

- **16 Lima** after applying project metadata curation to 11 PubMLST records with a
  blank town/city field;
- **7 Iquitos**;
- 2023: 4 records;
- 2024: 6 records;
- 2025: 13 records.

The Iquitos set is six human-stool isolates from 2024 plus one chicken
offal/meat isolate from 2025. The Peru collection is therefore not purely clinical:
it spans humans, poultry, goats, cattle and beef.

Raw PubMLST `town_or_city` values are never overwritten. The 11 additional Lima
assignments are stored separately in `data/peru_city_curation.tsv`.

## Interpretation

The authenticated PubMLST result is best treated as a **current public/database
snapshot**, not as an unbiased global surveillance sample. Peru contributes most of
the non-human diversity in the 41-record PubMLST set, whereas the Azevedo dataset
provides much broader European surveillance coverage.

The next formal step is to identify exact overlap between the 41 PubMLST records and
the Azevedo 217 using stable identifiers/accessions before concatenating sequence
sets.


## Exact overlap with Azevedo et al.

Exact stable-identifier matching found **9** current PubMLST records already represented
in the 217-isolate Azevedo collection:

| PubMLST | Azevedo | published Figure 4 status |
|---:|---|---|
| 112734 | LU-11 | singleton_68 |
| 118912 | UK-1 | singleton_64 |
| 118906 | UK-2 | excluded |
| 119062 | UK-3 | singleton_60 |
| 148932 | UK-4 | excluded |
| 148958 | UK-5 | singleton_82 |
| 148964 | UK-6 | excluded |
| 119231 | UK-7 | cluster_21 |
| 151737 | ES-1 | cluster_21 |

PubMLST 151737 has no contigs in BIGSdb, but the same isolate is represented by
Azevedo ES-1 with run accession SRR31589120. Thus its absence from the 40-genome
PubMLST set does not remove ES-1 from the combined sequence analysis.

After exact-identifier de-duplication, **32 PubMLST records are not exactly matched**
to the Azevedo set: 23 Peru and 9 non-Peru. This gives **249 provisional candidate
isolates** (217 + 32) before additional sequence-/metadata-level duplicate QC.

The 9 unmatched non-Peru PubMLST records are from Portugal (4), Luxembourg (1),
Spain (1), UK (1), UK Northern Ireland (1), and Vietnam (1).

A useful caution from the exploratory tree is that LU-11 (the Luxembourg 2021
environmental-water isolate) is a published **singleton_68**, not cluster 21.
Therefore proximity to LU-11 in the initial FastTree cannot be interpreted as
evidence that the Peru genomes belong to published cluster 21. UK-7 and ES-1 are
the two exact-overlap records assigned to published cluster 21; ES-1 is absent from
the PubMLST-contig tree because BIGSdb has no contigs for that record.
