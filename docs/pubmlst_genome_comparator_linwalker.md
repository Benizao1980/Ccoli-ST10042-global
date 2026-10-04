# PubMLST Genome Comparator → LINwalker workflow

This is the primary current-cgMLST workflow for the accessible ST10042 dataset.

## Inputs already prepared

Run:

```bash
python scripts/32_package_pubmlst_genome_comparator.py
```

This produces:

- `results/pubmlst_genome_comparator/azevedo_110_filtered_fastas.tar.gz`
- `results/pubmlst_genome_comparator/azevedo_110_metadata.tsv`
- `results/pubmlst_genome_comparator/current_pubmlst_st10042_anchors.tsv`
- `results/pubmlst_genome_comparator/current_pubmlst_st10042_anchor_ids.txt`

The archive contains the 110 filtered Azevedo assemblies. The anchor table/list
contains the 40 current PubMLST ST10042 records with contigs and their official
current cgMLST-v2/LIN metadata where available.

## PubMLST Genome Comparator

In the Campylobacter PubMLST isolate database:

1. Open **Genome Comparator**.
2. Paste/select the isolate IDs from
   `current_pubmlst_st10042_anchor_ids.txt`.
3. Upload `azevedo_110_filtered_fastas.tar.gz` as user genomes.
4. Analyse **against defined loci** and select
   **C. jejuni / C. coli cgMLST v2** (scheme 8; 1,142 loci).
5. Include the isolate/name field in displayed identifiers if offered. BIGSdb may
   internally label uploaded genomes as `u1`, `u2`, etc.; LINwalker can recover
   the `AZE_...` filename when it is included in the label.
6. Keep/download the Excel workbook. The `all` worksheet is the primary input for
   LINwalker because it retains per-locus allele calls (`X` = missing,
   `I` = incomplete, and `New#...` for alleles not already numbered).
7. The Genome Comparator distance matrix is useful as a cross-check, but LINwalker
   recalculates query-to-reference distances from the per-locus calls using the
   BIGSdb LIN missing-data normalisation.

Do not interpret uploaded genomes as having acquired official PubMLST cgSTs or
LINcodes simply because they fall close to database isolates.

## Run LINwalker

Create the pinned environment once:

```bash
conda env create -f environment-linwalker.yml
```

Then:

```bash
conda activate st10042-linwalker

bash scripts/34_run_linwalker_placement.sh \
  /path/to/GenomeComparator_output.xlsx
```

LINwalker uses the current 18-component Campylobacter v2 LIN thresholds:

```text
1119,1085,982,914,857,680,445,343,183,86,43,10,7,5,3,2,1,0
```

and reports:

- nearest official PubMLST reference genome(s);
- raw allele differences;
- BIGSdb-style missing-data-normalised allele distance;
- ST and clonal-complex context of the nearest references;
- an exact official cgST only when the uploaded profile is a complete exact match
  to a single official reference cgST;
- the deepest conservatively supported existing LIN prefix;
- supported/ambiguous `Cjc_cgc2_200/100/50/25/10/5` placements.

The local status labels distinguish `EXACT_REFERENCE_GENOME`,
`SUPPORTED_PREFIX` and `UNRESOLVED`. A supported prefix is not an official new
LINcode.

## Validation before interpretation

Inspect:

```text
results/linwalker_st10042/tables/reference_label_map.tsv
```

before using biological results. Every database anchor included in the Genome
Comparator run should map to exactly one PubMLST ID. Ambiguous or missing mappings
must be resolved before interpreting distances.

The project-specific joined output is:

```text
results/st10042_v2_lin_placement.tsv
```

and `scripts/33_summarise_linwalker_placement.py` prints the mapping of the 23
accessible published Azevedo cluster-21 genomes onto the current PubMLST LIN/cgc2
hierarchy, including whether their nearest current reference is a Peru isolate or
a published cluster-21 PubMLST anchor.
