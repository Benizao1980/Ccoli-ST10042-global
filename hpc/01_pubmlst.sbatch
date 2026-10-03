#!/bin/bash
#SBATCH --job-name=st10042_pubmlst
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=st10042_pubmlst-%j.out
#SBATCH --error=st10042_pubmlst-%j.err

set -eo pipefail

echo "ST10042 PubMLST job starting"
echo "date: $(date)"
echo "host: $(hostname)"
echo "workdir: $(pwd)"
echo "conda env: ${CONDA_DEFAULT_ENV:-<not set>}"
echo "python: $(command -v python || true)"

# sbatch inherits the submitting shell environment by default on Puma.
# Submit this job from an activated 'st10042' environment; this avoids sourcing
# Puma's /etc/bashrc inside the batch shell.
if [[ "${CONDA_DEFAULT_ENV:-}" != "st10042" ]]; then
  echo "ERROR: activate the st10042 conda environment before sbatch." >&2
  exit 2
fi

if [[ -z "${PUBMLST_API_KEY:-}" ]]; then
  echo "ERROR: PUBMLST_API_KEY is not exported into the sbatch environment." >&2
  echo "Run: export PUBMLST_API_KEY='...' before sbatch." >&2
  exit 2
fi
echo "PubMLST API key: present"

set -u

# Current global ST10042 population + all available PubMLST assemblies.
python scripts/01_query_pubmlst_st10042.py \
  --out data/pubmlst_st10042_all.tsv \
  --json-out data/pubmlst_st10042_all.json \
  --download-contigs data/pubmlst_st10042_assemblies

# Focal Peru metadata subset (assemblies are already present above).
python scripts/01_query_pubmlst_st10042.py \
  --country Peru \
  --out data/peru_st10042_current.tsv \
  --json-out data/peru_st10042_current.json

echo
echo "Global ST10042 records:"
awk 'END{print NR-1}' data/pubmlst_st10042_all.tsv
echo "Downloaded PubMLST assemblies:"
find data/pubmlst_st10042_assemblies -name '*.fasta' | wc -l
echo "Peru ST10042 records:"
awk 'END{print NR-1}' data/peru_st10042_current.tsv
