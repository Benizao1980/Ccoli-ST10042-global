#!/bin/bash
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --job-name=st10042_pubmlst
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --output=st10042_pubmlst-%j.out
#SBATCH --error=st10042_pubmlst-%j.err

set -eo pipefail

echo "ST10042 PubMLST job starting"
echo "date: $(date)"
echo "host: $(hostname)"
echo "workdir: $(pwd)"
echo "Slurm memory/cpu: ${SLURM_MEM_PER_CPU:-<not set>} MB"
echo "inherited conda env: ${CONDA_DEFAULT_ENV:-<not set>}"
echo "inherited CONDA_EXE: ${CONDA_EXE:-<not set>}"

if [[ -z "${CONDA_EXE:-}" || ! -x "${CONDA_EXE}" ]]; then
  echo "ERROR: CONDA_EXE was not inherited. Submit from an activated conda shell." >&2
  exit 2
fi

eval "$("${CONDA_EXE}" shell.bash hook)"
conda activate st10042

echo "active conda env: ${CONDA_DEFAULT_ENV:-<not set>}"
echo "python: $(command -v python)"
python --version
echo "bigsdb-downloader: $(command -v bigsdb-downloader || true)"

PUBMLST_TOKEN_DIR="${PUBMLST_TOKEN_DIR:-$HOME/.bigsdb_tokens}"
PUBMLST_KEY_NAME="${PUBMLST_KEY_NAME:-PubMLST}"
export PUBMLST_TOKEN_DIR PUBMLST_KEY_NAME

if ! command -v bigsdb-downloader >/dev/null 2>&1; then
  echo "ERROR: bigsdb-downloader is missing. Run: conda env update -f environment.yml" >&2
  exit 2
fi

if [[ ! -d "$PUBMLST_TOKEN_DIR" ]]; then
  echo "ERROR: OAuth tokens not found at $PUBMLST_TOKEN_DIR" >&2
  echo "Run the one-time bigsdb-downloader --setup step on the login node." >&2
  exit 2
fi

echo "PubMLST OAuth token directory: present"
set -u

python scripts/01_query_pubmlst_st10042.py \
  --auth oauth \
  --out data/pubmlst_st10042_all.tsv \
  --json-out data/pubmlst_st10042_all.json \
  --download-contigs data/pubmlst_st10042_assemblies

python scripts/01_query_pubmlst_st10042.py \
  --auth oauth \
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
