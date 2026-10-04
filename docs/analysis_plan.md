# Analysis plan: global *Campylobacter coli* ST10042

## Primary question
Where do the Peru ST10042 genomes sit within the current global ST10042 population,
and specifically relative to the European multidrug-resistant lineage described by
Azevedo et al. (2026)?

## Stage 1 — exact reproduction of the European reference analysis
1. Start with the 217 Azevedo isolates in `data/europe_217_manifest.tsv`.
2. Resolve and obtain sequence data using the published run accessions.
3. Reproduce assembly/QC as closely as possible to the study.
4. Call PubMLST Campylobacter cgMLST v1 (1,343 loci) with chewBBACA 3.1.2.
5. Retain exact/inferred calls and exclude genomes with <95% loci called.
6. Reproduce the 211-genome dataset and published <=5 AD clustering.
7. Required checkpoint: cluster 21 = 44 isolates, cluster 20 = 15, cluster 15 = 10.

## Stage 2 — current PubMLST ST10042 population
1. Query all current PubMLST ST10042 isolate records.
2. Use authenticated access so post-2024 records are included.
3. Retain all metadata records for ascertainment/distribution analysis.
4. Retrieve sequence only for records with genome assemblies. Records without contigs remain in the metadata/distribution dataset but are excluded from cgMLST.
5. De-duplicate against the Azevedo dataset using run accession, BioSample/assembly,
   PubMLST ID, then strain/isolate metadata as fallback.
6. Preserve provenance flags (`Azevedo`, `PubMLST`, `Peru`) rather than discarding overlap.
7. Treat 249 (217 + 32 exact-unmatched PubMLST records) as a provisional pre-QC candidate count only; inspect unmatched European records for hidden duplicate submissions before finalising N.

## Stage 3 — focal Peru analysis
For each Peru genome calculate:
- cgMLST locus-call rate
- nearest global ST10042 neighbour
- minimum allelic distance to any Azevedo genome
- minimum allelic distance to published cluster 21
- <=5 AD cluster/component membership
- source, city and year context

## Stage 4 — AMR
Compare the Peru/global population for the lineage-defining markers reported by Azevedo:
- GyrA Thr86Ile
- `tet(O/32/O)` / `tet(O)`
- `blaOXA-61` and promoter -57G>T
- promoter -69delA
- `porA` allele
- 23S rRNA macrolide resistance mutations
- `cmeRABC`

## Interpretation guardrails
The seven-locus ST alone does not prove membership of the European MDR epidemic lineage.
FastTree topology alone should not be used to infer direction of transmission or a single
introduction. The first formal test is the matched cgMLST analysis; directionality requires
additional temporal/epidemiological evidence and ideally recombination-aware SNP analysis.
