#!/bin/bash
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --job-name=st10042_covpilot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=03:00:00
#SBATCH --output=st10042_covpilot-%j.out
#SBATCH --error=st10042_covpilot-%j.err

set -euo pipefail

if [[ -z "${CONDA_EXE:-}" || ! -x "${CONDA_EXE}" ]]; then
  echo "ERROR: CONDA_EXE was not inherited. Submit from activated st10042-assembly." >&2
  exit 2
fi

eval "$("${CONDA_EXE}" shell.bash hook)"
conda activate st10042-assembly

python scripts/21_validate_innuca_like_filter.py \
  --threads "${SLURM_CPUS_PER_TASK}"
