# Global population genomics of *Campylobacter coli* ST10042

This repository reproduces the 2026 Azevedo et al. European ST10042 analysis and
extends it to **all current PubMLST ST10042 genomes**, with a focal analysis of the
Peru population (Lima and Iquitos).

## Core questions

1. Can we reproduce the published 211-genome European cgMLST analysis and its major clusters?
2. What is the current global distribution of ST10042 in PubMLST?
3. Do the Peru ST10042 genomes fall within the published European <=5-allele clusters,
   particularly cluster 21, or do they form distinct global clusters?
4. Does the conserved European MDR genotype extend to the Peru/global population?

## Published reference analysis

Azevedo et al. analysed 217 ST10042 isolates from nine European countries. Their
Figure 4 retained 211 genomes after cgMLST completeness filtering. The analysis used:

- PubMLST Campylobacter cgMLST v1: 1,343 loci
- chewBBACA 3.1.2
- BSR threshold 0.6
- size threshold 0.2
- exact + inferred allele assignments for the profile matrix
- >=95% cgMLST loci called
- ReporTree 2.5.4
- GrapeTree MSTreeV2
- clustering at a maximum of 5 allelic differences

The cleaned appendix and published cluster assignments are included under `data/`.

## Repository layout

```text
config/     fixed analysis settings and reproduction targets
data/       small tracked manifests only; no raw sequence data
docs/       analysis plan
hpc/        HPC quick start + SLURM wrappers
scripts/    PubMLST/ENA retrieval and analysis helpers
results/    generated outputs (ignored by git)
schema/     local cgMLST schema (ignored by git)
```

## Quick start on the HPC

```bash
git clone git@github.com:Benizao1980/Ccoli-ST10042-global.git
cd Ccoli-ST10042-global
conda env create -f environment.yml
conda activate st10042
```

PubMLST requires authentication for records added after 31 December 2024.
Create a personal API key in your PubMLST/BIGSdb profile and keep it outside Git:

```bash
export PUBMLST_API_KEY='...'
```

Retrieve **all current ST10042 metadata and available PubMLST assemblies**, plus a
separate Peru manifest:

```bash
sbatch hpc/01_pubmlst.sbatch
```

Equivalent direct command for the global pull:

```bash
python scripts/01_query_pubmlst_st10042.py \
  --out data/pubmlst_st10042_all.tsv \
  --json-out data/pubmlst_st10042_all.json \
  --download-contigs data/pubmlst_st10042_assemblies
```

The exploratory query used to start this project contained 23 confirmed Peru ST10042
records (16 Lima, 7 Iquitos). This count is a validation checkpoint only: the analysis
always uses the live PubMLST query rather than hard-coding those isolates.

## Reproducibility checkpoint

Before interpreting the expanded global dataset, reproduce the published European result:

- 217 input isolates
- 211 genomes retained for Figure 4
- cluster 21: 44 genomes
- cluster 20: 15 genomes
- cluster 15: 10 genomes

Only after that checkpoint passes should the global/Peru population be added and
cluster labels recomputed.

## Data policy

Raw reads, assemblies, credentials, schemas and large generated outputs are not committed.
Small public metadata/manifests required to reproduce the analysis are tracked.

## Status

Initial analysis scaffold. Unpublished Peru/global results should be treated as preliminary.
