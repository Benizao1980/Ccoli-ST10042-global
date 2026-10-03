#!/bin/bash
#SBATCH --job-name=st10042_ena
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=st10042_ena-%j.out

# Puma's /etc/bashrc is not safe under nounset during initialisation.
set -eo pipefail
source ~/.bashrc
set -u

conda activate st10042

python scripts/02_fetch_ena_reads.py \
  --manifest data/europe_217_manifest.tsv \
  --resolved data/europe_ena_fastq_manifest.tsv
