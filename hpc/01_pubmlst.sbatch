#!/bin/bash
#SBATCH --job-name=st10042_pubmlst
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=st10042_pubmlst-%j.out

set -euo pipefail
source ~/.bashrc
conda activate st10042

: "${PUBMLST_API_KEY:?Set PUBMLST_API_KEY before submitting the job}"

python scripts/01_query_pubmlst_st10042.py \
  --out data/pubmlst_st10042_all.tsv \
  --json-out data/pubmlst_st10042_all.json

python scripts/01_query_pubmlst_st10042.py \
  --country Peru \
  --out data/peru_st10042_current.tsv \
  --json-out data/peru_st10042_current.json \
  --download-contigs data/peru_assemblies
