#!/bin/bash
set -euo pipefail

script="hpc/02_resolve_ena.sbatch"

if [[ "$(head -n 1 "$script")" != "#!/bin/bash" ]]; then
  echo "ERROR: $script is malformed: first line must be #!/bin/bash" >&2
  exit 2
fi

unset SBATCH_MEM_PER_NODE SBATCH_MEM_PER_CPU SBATCH_MEM_PER_GPU

exec sbatch \
  --account=cooperma \
  --partition=standard \
  --cpus-per-task=1 \
  --mem-per-cpu=4G \
  "$script"
