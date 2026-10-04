# Analysis plan: global *Campylobacter coli* ST10042

## Primary question
Where do the Peru ST10042 genomes sit within the current global ST10042 population,
and specifically relative to the European multidrug-resistant lineage described by
Azevedo et al. (2026)?

## Stage 1 — published reference and current sequence inventory
1. Start with the 217 Azevedo isolates in `data/europe_217_manifest.tsv`.
2. Retain the published cgMLST v1 cluster assignments as reference annotations.
3. Resolve all currently public sequence data without pretending the unavailable
   106 Azevedo records can be reproduced.
4. Use the 23 currently accessible published cluster-21 genomes as labelled anchors.
5. Keep the published v1/<=5 AD definition separate from the primary modern analysis:
   v2/LIN results must not be described as a literal reproduction of Azevedo cluster 21.

## Stage 2 — current PubMLST ST10042 population
1. Query all current PubMLST ST10042 isolate records.
2. Use authenticated access so post-2024 records are included.
3. Retain all metadata records for ascertainment/distribution analysis.
4. Retrieve sequence only for records with genome assemblies. Records without contigs remain in the metadata/distribution dataset but are excluded from cgMLST.
5. De-duplicate against the Azevedo dataset using run accession, BioSample/assembly,
   PubMLST ID, then strain/isolate metadata as fallback.
6. Preserve provenance flags (`Azevedo`, `PubMLST`, `Peru`) rather than discarding overlap.
7. Treat 249 (217 + 32 exact-unmatched PubMLST records) as a provisional pre-QC candidate count only; inspect unmatched European records for hidden duplicate submissions before finalising N.

## Stage 3 — PubMLST cgMLST v2 + LIN analysis
Use the current PubMLST/BIGSdb Campylobacter cgMLST v2 scheme and PubMLST LIN
classification as the primary nomenclature. Discover the active scheme IDs from the
live API with `scripts/11_pubmlst_v2_lin_preflight.py`; do not hard-code IDs from
the older cgMLST v1 analysis.

For each Peru genome calculate/retain:
- cgMLST v2 locus-call rate
- PubMLST v2 allelic profile
- PubMLST LIN code / hierarchical classification
- nearest global ST10042 neighbour under v2
- minimum v2 allelic distance to any accessible Azevedo genome
- minimum v2 allelic distance to an accessible published cluster-21 anchor
- source, city and year context

Published Azevedo v1 cluster membership remains an annotation. A v2/LIN grouping is
not assumed to be identical to a published v1 <=5-AD cluster.

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

## Current archive-availability checkpoint (2026-10-04)

The ENA file-report API currently resolves generated FASTQ files for **111/217**
Azevedo run accessions and returns no public FASTQ record for **106/217**. This is
an archive-availability issue rather than a Slurm/resource failure: the resolution
job completed successfully and wrote an explicit per-run inventory.

The paper was released in September 2026 and states that raw reads were deposited in
ENA, so the unresolved accessions should be treated as **not currently retrievable
through the public ENA file-report endpoint**, not as absent or invalid data. Exact
reproduction of the full 217-isolate analysis is therefore blocked until those
records become public/retrievable or an alternative public sequence source is found.

Use `scripts/07_diagnose_ena_availability.py` to quantify the effect by BioProject,
country and published cluster (especially cluster 21).


### Assembly rescue step

Because ENA raw-read and assembly records are separate archive domains, query the
assembly domain for studies represented among the 106 unresolved runs before
concluding that those genomes are inaccessible. Use
`scripts/08_probe_unresolved_assemblies.py`.

Any assembly recovered this way should be matched back to an Azevedo isolate using
stable exact identifiers (run reference, strain/ENA identifier) before use. Do not
infer a match solely from country/year or approximate names.


### FR-1 BioProject anomaly

The Azevedo supplementary manifest records `ERR16782028` in the BioProject field
for FR-1. This is a run accession, not a BioProject/study accession. Preserve the
published value in the source manifest rather than silently correcting it. Assembly
probing therefore validates project-accession syntax and handles FR-1 separately by
its exact run accession.


### ENA assembly schema note

The current ENA Portal API does not expose the historical `run_ref` field in the
genome-assembly return schema used by older examples. The rescue script therefore
discovers `returnFields` and `searchFields` at runtime and queries both
`assembly` and `analysis` (sequence-assembly) domains using only fields that ENA
currently advertises. This avoids treating API-schema drift as evidence that an
assembly is absent.


## Public sequence availability conclusion (2026-10-04)

The archive rescue check is now complete for the 106 Azevedo records without public
FASTQ. Using ENA's current advertised assembly and analysis schemas, the three valid
unresolved BioProjects return **0 public genome-assembly records** and **0 public
SEQUENCE_ASSEMBLY analysis records**. FR-1 (published with an ERR accession in the
BioProject field) likewise returns no public assembly/analysis object by exact run
accession. Query errors are zero.

Therefore, as of this checkpoint:

- Azevedo isolates with publicly retrievable FASTQ: **111/217**;
- Azevedo isolates without currently retrievable FASTQ or assembly/analysis object:
  **106/217**;
- published cluster 21 with publicly retrievable sequence: **23/44**;
- published cluster 21 currently unavailable: **21/44**.

The 41-record PubMLST snapshot contains 40 assemblies; the sole no-contig record
(PubMLST 151737) is an exact duplicate of Azevedo ES-1 and is represented by the
public Azevedo run. All 32 PubMLST-only records remaining after exact-ID
de-duplication therefore have sequence.

This yields a **currently sequence-accessible candidate set of 143 isolates before
cgMLST QC**: 111 Azevedo + 32 PubMLST-only additions. This includes all 23 Peru
focal isolates and 23 published cluster-21 anchors.

Do not interpret distance from the 23 accessible cluster-21 anchors as a complete
negative test of cluster-21 membership, because 21/44 published members are
currently unavailable. A <=5-allele match to an accessible published cluster-21
genome is positive evidence; absence of such a match is not definitive exclusion.

Use `scripts/09_build_accessible_manifest.py` to generate the current-access
manifest reproducibly.


## Assembly strategy

The primary analysis is not an exact reconstruction of the 2026 INNUca/cgMLST-v1
workflow. It uses current PubMLST cgMLST v2 and LIN nomenclature.

For **Illumina-only reads**, assemble consistently with **SPAdes in isolate mode**.
For samples with both Illumina and Oxford Nanopore reads, use **Unicycler hybrid
assembly**. Do not use Unicycler merely as a wrapper around SPAdes for the
Illumina-only Azevedo set: direct SPAdes is simpler and avoids introducing another
assembly/graph-processing layer.

Existing PubMLST contigs are retained as the database assemblies for the 32
PubMLST-only additions. Assembly method is therefore not perfectly homogeneous;
control this with assembly QC and cgMLST-v2 completeness thresholds, and optionally
perform a sensitivity analysis on records for which raw reads are also available.

Before downloading the 111 FASTQ datasets, run
`scripts/10_plan_azevedo_assembly.py` to quantify compressed download volume and
read-file layout. The separate `environment-assembly.yml` contains SPAdes,
Unicycler (for hybrid data), fastp and QUAST.


## Current PubMLST v2/LIN checkpoint (2026-10-04)

Live PubMLST confirms **cgMLST v2 = scheme 8 (1,142 loci)**. Its native LINcode
uses 18 thresholds:

`1119,1085,982,914,857,680,445,343,183,86,43,10,7,5,3,2,1,0`

with `max_missing=25`. The fixed `Cjc_cgc2_200/100/50/25/10/5`
classification schemes are convenient views of the same v2 population structure,
but are not the full LINcode.

The ENA inventory contains 111 nominally resolved Azevedo runs, but **DK-2
(ERR10702984) has no FASTQ files**. The actual assembly-ready Azevedo set is
therefore **110 paired-end runs**, giving a current analysis-ready candidate set of
**142 genomes = 110 Azevedo + 32 PubMLST-only assemblies**, before assembly/cgMLST
quality filtering.

A preliminary read of the current PubMLST isolate JSON shows an important result:
published Azevedo cluster-21 anchor **UK-7 / PubMLST 119231** is in
`Cjc_cgc2_5 group 12905`. Four Peru ST10042 records — **149884, 149885,
149886 and 149892** — are also in group 12905. Their native LINcodes share the
same prefix through the 5-allele LIN level and split only at finer thresholds.
This is strong current-cgMLST-v2 evidence that at least four Peru isolates occupy
the same <=5-AD PubMLST group as a published cluster-21 anchor.

This does **not** make the v2 group synonymous with Azevedo's published v1
cluster 21. The full accessible Azevedo set must still be assembled and typed under
v2 to determine how the published v1 cluster maps onto current LIN/cgc2 space.

Use `scripts/13_extract_pubmlst_v2_lincodes.py` for the clean current-PubMLST
table and focal Peru summary. The earlier generic inspector was intentionally broad
and previously matched the substring "lin" inside fields such as tetracycline; this
has been corrected.


## Azevedo SPAdes completion checkpoint (2026-10-04)

All **110/110** paired Azevedo read sets completed fastp + SPAdes `--isolate`
successfully (Slurm exit 0 for every array task), and every sample has a non-empty
`contigs.fasta` plus a fastp JSON report.

Median assembly statistics are consistent with a typical *Campylobacter* genome
(1.654 Mb >=500 bp, 27 contigs >=500 bp, N50 163.7 kb, GC 31.6%), but the cohort
contains at least one obvious gross outlier: maxima include 4.332 Mb assembly size,
1,656 contigs >=500 bp and 47.3% GC, with minimum N50 8.45 kb. Do not pass the
entire set directly into biological interpretation without identifying these
sample(s).

Use `scripts/18_flag_assembly_qc_outliers.py` to list review candidates. Its
size/GC/fragmentation thresholds are diagnostic only, not automatic exclusions.
Final inclusion for the main analysis should combine assembly sanity checks with
PubMLST cgMLST-v2 typing quality, especially the native `max_missing=25` rule.


### QC interpretation after SPAdes

Eleven of 110 assemblies cross conservative review thresholds. The strongest
outliers are PT-32 and PT-33 (large, high-GC, highly fragmented) and UK-2/UK-4/UK-6
(large and highly fragmented). Notably, UK-2, UK-4 and UK-6 are the same three UK
isolates excluded from Azevedo Figure 4, which supports the utility of the QC flags.

Do not automatically discard all 11 flagged assemblies. Several flagged genomes
were retained in the published analysis, including cluster-21 anchors PT-26, PT-32,
PT-33 and UK-7. Instead, validate allele recovery first. UK-7 is especially useful
because it has an existing official PubMLST cgMLST-v2/LIN assignment.

Use `scripts/19_validate_pubmlst_v2_on_overlaps.py` to query the current PubMLST
scheme-8 sequence endpoint sequentially for local SPAdes assemblies that also have
known PubMLST records. Compare exact allele recovery and returned scheme fields to
the official current record before scaling allele calling to all 110 Azevedo
assemblies. BIGSdb's scheme sequence endpoint reports exact allele matches; genuine
novel alleles will require a different/local calling route for the final analysis.
