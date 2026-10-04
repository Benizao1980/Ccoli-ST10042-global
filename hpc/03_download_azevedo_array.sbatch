#!/bin/bash
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --job-name=st10042_dl
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00
#SBATCH --output=st10042_dl-%A_%a.out
#SBATCH --error=st10042_dl-%A_%a.err

set -euo pipefail

if [[ -z "${CONDA_EXE:-}" || ! -x "${CONDA_EXE}" ]]; then
  echo "ERROR: CONDA_EXE was not inherited. Submit from activated st10042." >&2
  exit 2
fi

eval "$("${CONDA_EXE}" shell.bash hook)"
conda activate st10042

python scripts/15_download_azevedo_pair.py \
  --manifest results/azevedo_spades_manifest.tsv \
  --index "${SLURM_ARRAY_TASK_ID}"
