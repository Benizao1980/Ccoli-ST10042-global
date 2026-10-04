#!/bin/bash
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --job-name=st10042_ena
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=st10042_ena-%j.out
#SBATCH --error=st10042_ena-%j.err

set -eo pipefail

echo "ST10042 Azevedo ENA resolution starting"
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
set -u

echo "active conda env: ${CONDA_DEFAULT_ENV}"
echo "python: $(command -v python)"
python --version

python scripts/02_fetch_ena_reads.py \
  --manifest data/europe_217_manifest.tsv \
  --resolved data/europe_ena_fastq_manifest.tsv \
  --unresolved results/europe_ena_unresolved.tsv
