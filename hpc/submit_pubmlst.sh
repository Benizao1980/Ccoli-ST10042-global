#!/bin/bash
set -euo pipefail

# Known-good Puma configuration for the cooperma allocation:
# account=cooperma, partition=standard, memory via --mem-per-cpu.
unset SBATCH_MEM_PER_NODE SBATCH_MEM_PER_CPU SBATCH_MEM_PER_GPU
exec sbatch \
  --account=cooperma \
  --partition=standard \
  --cpus-per-task=1 \
  --mem-per-cpu=4G \
  hpc/01_pubmlst.sbatch
