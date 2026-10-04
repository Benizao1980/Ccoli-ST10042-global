#!/bin/bash
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --job-name=st10042_spades
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --output=st10042_spades-%A_%a.out
#SBATCH --error=st10042_spades-%A_%a.err

set -euo pipefail

if [[ -z "${CONDA_EXE:-}" || ! -x "${CONDA_EXE}" ]]; then
  echo "ERROR: CONDA_EXE was not inherited. Submit from an activated conda shell." >&2
  exit 2
fi

eval "$("${CONDA_EXE}" shell.bash hook)"
conda activate st10042-assembly

python scripts/16_assemble_azevedo_spades.py \
  --manifest results/azevedo_spades_manifest.tsv \
  --index "${SLURM_ARRAY_TASK_ID}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --memory-gb 15
