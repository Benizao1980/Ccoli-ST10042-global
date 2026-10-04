#!/bin/bash
set -euo pipefail

script="hpc/05_innuca_like_overlap.sbatch"
[[ "$(head -n 1 "$script")" == "#!/bin/bash" ]] || {
  echo "ERROR: malformed $script" >&2; exit 2;
}

unset SBATCH_MEM_PER_NODE SBATCH_MEM_PER_CPU SBATCH_MEM_PER_GPU

exec sbatch \
  --account=cooperma \
  --partition=standard \
  --cpus-per-task=4 \
  --mem-per-cpu=4G \
  "$script"
