# Data layout

Only small metadata/manifests are tracked in Git.

Tracked:
- `europe_217_manifest.tsv`: cleaned Azevedo et al. appendix (217 isolates)
- `cluster21_44.tsv`: published cluster 21 subset
- `excluded_from_figure4_6.tsv`: six isolates excluded from the paper's 211-genome Figure 4 analysis
- `europe_run_accessions.txt`: European run accessions
- `known_peru_2023_st10042.tsv`: four Peru ST10042 isolates identifiable in the older Peru manuscript metadata

Not tracked:
- raw FASTQ
- genome assemblies
- PubMLST contig downloads
- cgMLST schema files
- large intermediate/output files

The current global PubMLST ST10042 population should be generated at run time using
`scripts/01_query_pubmlst_st10042.py` rather than committed as a frozen source of truth.
