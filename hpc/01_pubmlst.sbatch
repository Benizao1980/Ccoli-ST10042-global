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
echo "inherited conda env: ${CONDA_DEFAULT_ENV:-<not set>}"
echo "inherited CONDA_EXE: ${CONDA_EXE:-<not set>}"

# Re-initialise conda explicitly inside the batch shell.  Slurm can inherit
# CONDA_DEFAULT_ENV without inheriting the matching PATH, which previously
# left us with /usr/bin/python even though the environment name was st10042.
if [[ -z "${CONDA_EXE:-}" || ! -x "${CONDA_EXE}" ]]; then
  echo "ERROR: CONDA_EXE was not inherited. Submit from an activated conda shell." >&2
  exit 2
fi
eval "$("${CONDA_EXE}" shell.bash hook)"
conda activate st10042

echo "active conda env: ${CONDA_DEFAULT_ENV:-<not set>}"
echo "python: $(command -v python)"
python --version

if [[ -z "${PUBMLST_API_KEY:-}" ]]; then
  echo "ERROR: PUBMLST_API_KEY is not exported into the sbatch environment." >&2
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
