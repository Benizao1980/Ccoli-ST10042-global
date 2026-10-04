#!/bin/bash
set -euo pipefail

script="hpc/10_chewbbaca_overlap_pilot.sbatch"
[[ "$(head -n 1 "$script")" == "#!/bin/bash" ]] || {
  echo "ERROR: malformed $script" >&2; exit 2;
}

unset SBATCH_MEM_PER_NODE SBATCH_MEM_PER_CPU SBATCH_MEM_PER_GPU

exec sbatch \
  --account=cooperma \
  --partition=standard \
  --cpus-per-task=8 \
  --mem-per-cpu=4G \
  "$script"
