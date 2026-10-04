# Analysis plan: global *Campylobacter coli* ST10042

## Biological questions

1. **Are the Peru ST10042 isolates part of the same genomic lineage as the European
   isolates described by Azevedo et al. (2026)?**
2. **When did that lineage emerge?**
3. **Where did it most likely emerge, geographically and by source reservoir?**
4. **Are the Peru isolates also multidrug resistant, and are they associated with
   diarrhoeal disease rather than asymptomatic carriage in children?**

The immediate technical objective is to define the lineage robustly using current
PubMLST cgMLST v2 + LINcodes. Temporal, geographic/source and clinical inference
comes only after that genomic definition is established.

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

## Stage 5 — pangenome, AMR and virulence

After final cgMLST-v2 QC/de-duplication, use one fixed genome set for comparative
genomics so gene-content differences are not driven by assembly failures.

### Pangenome

1. Re-annotate all final assemblies consistently with the same annotation workflow.
2. Run **Panaroo** on the common annotation set.
3. Report:
   - core, soft-core, shell and cloud gene counts;
   - gene-presence/absence matrix;
   - core-gene alignment;
   - lineage-, country- and source-associated accessory genes;
   - sensitivity of accessory-gene calls to the remaining fragmented assemblies.
4. Overlay cgMLST/LIN groups, Azevedo cluster labels, country, source and Peru
   epidemiological metadata onto the pangenome structure.

Panaroo is used to characterise gene content, not to define the primary lineage:
the main nomenclature remains PubMLST cgMLST v2 + LINcodes.

### AMR

Run **AMRFinderPlus** across the final assembly set, retaining both gene and point-
mutation calls where supported. Then explicitly verify the resistance determinants
highlighted by Azevedo et al.:

- GyrA Thr86Ile;
- `tet(O/32/O)` / `tet(O)`;
- `blaOXA-61`;
- 23S rRNA macrolide-resistance mutations;
- `cmeRABC`.

AMRFinderPlus alone should not be assumed to recover every regulatory variant of
interest. Add targeted sequence checks for the `blaOXA-61` promoter -57G>T,
promoter -69delA and any other lineage-defining non-coding variants that are not
represented in the AMRFinderPlus catalogue.

Summarise AMR by LIN group, published cluster, country, source and Peru child
phenotype. Define MDR using the same antibiotic-class rule across all genomes.

### Virulence

Screen the same final genomes against **VFDB**, using a reproducible homology
threshold and recording identity/coverage rather than only binary presence/absence.
Confirm biologically important or borderline hits at protein level where needed.

For a closely related ST10042 population, many canonical *Campylobacter* virulence
genes may be core. Therefore also consider allelic variation, truncation/pseudogene
status and accessory-gene context rather than interpreting simple presence/absence
as a proxy for virulence.

## Stage 6 — recombination-aware phylogenomics

Recombination must be assessed before molecular-clock dating or directionality/origin
inference.

1. Define the final ST10042 analysis set and the focal cluster-21-related lineage
   using cgMLST v2/LINcodes.
2. Build a high-resolution core-genome/SNP alignment for the whole ST10042 set and,
   separately, for the focal close lineage where reference-based alignment is most
   defensible.
3. Use **Gubbins** (and, if useful, ClonalFrameML as a sensitivity analysis) to:
   - identify recombinant regions;
   - estimate recombination burden (including r/m where supported);
   - generate a recombination-masked clonal alignment/tree;
   - identify recurrent recombination hotspots.
4. Compare masked and unmasked topologies, especially for the four Lima 2023 animal
   isolates, the Iquitos human group and other Peru source-associated subclusters.

The pangenome tree and cgMLST/LIN structure provide complementary views; the
recombination-masked phylogeny is the preferred basis for dating.

## Stage 7 — temporal analysis and lineage dating

Do not fit a molecular clock automatically.

1. Use the recombination-masked tree/alignment to test temporal signal across the
   whole focal lineage.
2. Examine root-to-tip regression and perform date-randomisation/permutation tests.
3. If temporal signal is adequate, estimate:
   - substitution rate;
   - TMRCA of the focal ST10042/cluster-21-related lineage;
   - TMRCA of Peru-associated sublineages where the sampling dates contain enough
     temporal spread and substitutions to support a separate estimate.
4. Use **BactDating** as an efficient primary bacterial-dating framework, with a
   second method (for example BEAST2 or another clock framework) for key sensitivity
   analyses if warranted.
5. Report uncertainty intervals and avoid dating subclusters sampled within only one
   or two narrow time points.

The European collection spans 2018–2025, but several Peru subgroups are concentrated
in 2023–2025. Those short windows may be insufficient for independent subcluster
dating even if the broader lineage contains usable temporal signal.

## Stage 8 — geographic and source-reservoir history

Only after the recombination-aware dated lineage is established, reconstruct
geographic and source history.

- Treat country and source reservoir as separate discrete traits.
- Compare ancestral-state support rather than assigning an origin from the oldest
  sampled genome.
- Explicitly account for the strong surveillance imbalance and the 107 Azevedo
  genomes that are not currently sequence-accessible.
- Where appropriate, combine phylogenetic ancestral-state reconstruction with
  source-attribution models rather than relying on one method.

The desired outputs are estimates of the most plausible ancestral geography/source
states and the timing of major transitions, with uncertainty and sampling bias made
explicit.

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

The ENA inventory initially yielded 111 nominally resolved Azevedo runs. A later
file-level check showed that DK-2 (ERR10702984) has no FASTQ files, leaving **110
assembly-ready Azevedo isolates**. Together with the 32 PubMLST-only additions,
this yields a **currently sequence-accessible candidate set of 142 genomes before
cgMLST QC**. This includes all 23 Peru focal isolates and the currently accessible
published cluster-21 anchors.

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


### Exact-hit ambiguity checkpoint

The first scheme-8 API pilot recovered nearly all loci for several known overlaps
but found multiple exact allele hits in some assemblies. This should not be treated
as equivalent to PubMLST's `max_missing=25`: UK-7, for example, has an official
LINcode despite multiple cgST values in its current PubMLST record. The initial
pilot's requirement of zero ambiguous loci was therefore too conservative and has
been removed.

Before scaling, compare each local SPAdes overlap assembly with the official
PubMLST contigs for the same isolate using
`scripts/20_compare_overlap_assemblies.py`. The comparison also tests a >=500-bp
filtered version of the local SPAdes FASTA, because the first pilot queried every
SPAdes contig including small low-coverage graph fragments. This will establish
whether ambiguity is driven by our assembly representation or is intrinsic to the
stored isolate.


### INNUca-like coverage-filter pilot

The direct overlap comparison shows that >=500-bp length filtering does not solve
the multiple-exact-hit problem. UK-7 remains 75 ambiguous loci locally versus 0 in
the official PubMLST assembly; UK-3 and UK-5 show the same pattern. Conversely,
the official assemblies for UK-2/UK-4/UK-6 remain poor by scheme-8 coverage, which
supports their exclusion as sample-level failures rather than a local-assembly-only
artifact.

The published Azevedo workflow used INNUca. INNUca applies more than a simple
length filter: its documented defaults include SPAdes k-mer coverage filtering and
a read-mapping post-assembly filter that retains contigs with mean mapped coverage
>= max(10x, one-third of the assembly-wide mean coverage). We did not apply that
mapping filter in the first SPAdes pass.

Use `scripts/21_validate_innuca_like_filter.py` on the known overlaps before
scaling. It approximates the relevant INNUca post-assembly steps using the current
SPAdes contigs and fastp-processed reads: contig length >=200 bp, SPAdes k-mer
coverage >=2, Bowtie2 very-sensitive-local mapping, then the same dynamic mapped-
coverage threshold. It then re-queries PubMLST scheme 8. This is deliberately
described as INNUca-like rather than an exact reconstruction because our trimming
and software versions differ from the published pipeline.


#### Trimmed-read correction

The first successful INNUca-like overlap pilot strongly removed duplicate allele
hits, but inspection of the pilot code found that the mapping step used the raw
downloaded ENA FASTQs from the manifest rather than the fastp-trimmed pairs that
were actually supplied to SPAdes. INNUca itself updates `fastq_files` to its
trimmed paired reads before the assembly-mapping step. The pilot has therefore
been corrected to map `data/europe_trimmed/<sample>/<sample>_R1/R2.fastq.gz`
and writes to a new versioned output directory/cache so the previous raw-read
results cannot be silently reused.

The raw-read pilot remains useful evidence: after coverage filtering, LU-11,
ES-1, UK-1, UK-3 and UK-7 all had zero ambiguous scheme-8 loci and UK-5 had one;
UK-2/UK-4/UK-6 remained poor.

The **corrected trimmed-read pilot has now reproduced the same pattern**. LU-11,
ES-1, UK-1, UK-3 and UK-7 have zero ambiguous scheme-8 loci; UK-5 has one.
UK-7, the key published cluster-21 anchor, has 1,139/1,142 loci with exact hits,
all 1,139 uniquely called, and only 3 loci without an exact hit. UK-2, UK-4 and
UK-6 remain above the PubMLST v2 `max_missing=25` completeness threshold
(53, 87 and 66 loci without an exact hit, respectively), matching their poor
assembly profiles and published Figure-4 exclusion.

This validates the INNUca-like coverage-filtering strategy for scale-up to all 110
accessible Azevedo read sets.


## Full 110-genome coverage-filter completion (2026-10-04)

The validated INNUca-like coverage filter has now completed successfully for
**110/110 accessible Azevedo genomes** (all Slurm tasks COMPLETED; 0 incomplete
outputs).

Post-filter cohort summary:

- median filtered assembly size: **1,648,401 bp**;
- median filtered contig count: **27**;
- filtered size range: **1,632,875–2,517,580 bp**;
- filtered contig range: **18–404**;
- median mapped mean depth: **128.45x**.

UK-4, UK-6 and UK-2 remain the three most fragmented filtered UK assemblies
(404, 233 and 216 contigs respectively), consistent with the overlap pilot and
the published Figure-4 exclusions. Other larger/fragmented Portuguese genomes
remain in the dataset pending cgMLST-v2 QC; fragmentation alone is not used as
an exclusion criterion.

The next checkpoint is a **sequential scheme-8 exact-hit QC across all 110 filtered
assemblies** using `scripts/24_pubmlst_v2_qc_all_azevedo.py`. This step is
diagnostic only: a locus without an exact known-allele hit can represent a genuine
novel allele, so `loci_without_exact_hit` must not yet be treated as the final
missing-locus count.


## Full 110-genome scheme-8 exact-hit QC (2026-10-04)

Sequential PubMLST scheme-8 exact matching has completed for all 110 coverage-filtered
Azevedo assemblies.

- **105/110** have <=25 loci without an exact known-allele hit;
- **106/110** have zero ambiguous exact-hit loci;
- the five >25 genomes are UK-4 (87), UK-6 (66), UK-2 (53), PT-41 (43) and PT-58 (28);
- UK-4/UK-6/UK-2 remain obvious sample-level failures and match the published
  Figure-4 exclusions;
- PT-41/PT-58 have zero ambiguity and should not be excluded until full allele
  calling distinguishes true missing loci from novel alleles;
- **all 23 accessible published cluster-21 genomes pass the exact-hit diagnostic
  comfortably (1–10 loci without exact hits; 0 ambiguous loci).**

The exact-hit REST endpoint cannot classify novel alleles and therefore is not the
final cgMLST profile. The next reproducibility checkpoint is to freeze the live
PubMLST cgMLST-v2 scheme (all 1,142 locus allele FASTAs, profiles, LIN definitions,
and classification-scheme metadata) with checksums using
`scripts/25_snapshot_pubmlst_v2_schema.py`.

Then build the exact 142-genome input list with
`scripts/26_prepare_cgmlst_v2_input.py` and adapt the frozen external schema for
chewBBACA. chewBBACA can infer novel alleles while retaining exact matches to an
external BIGSdb schema, making it suitable for the final pairwise allelic-distance
matrix. Any LIN-like assignments calculated locally must be distinguished from
official PubMLST LINcodes unless PubMLST itself returns an assignment.
