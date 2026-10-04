#!/bin/bash
set -euo pipefail

script="hpc/11_map_azevedo_to_official_pubmlst_v2.sbatch"

unset SBATCH_MEM_PER_NODE SBATCH_MEM_PER_CPU SBATCH_MEM_PER_GPU

exec sbatch \
  --account=cooperma \
  --partition=standard \
  --cpus-per-task=1 \
  --mem-per-cpu=4G \
  "$script"
