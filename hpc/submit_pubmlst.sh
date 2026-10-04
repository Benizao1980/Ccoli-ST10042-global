#!/bin/bash
set -euo pipefail

# Slurm SBATCH_* environment variables override directives inside a batch script.
# Clear memory overrides so a stale shell setting cannot reduce the job to 1 MB.
unset SBATCH_MEM_PER_NODE SBATCH_MEM_PER_CPU SBATCH_MEM_PER_GPU

exec sbatch --account=cooperma --mem=5gb hpc/01_pubmlst.sbatch
