#!/bin/bash
set -euo pipefail

manifest="results/azevedo_spades_manifest.tsv"
script="hpc/03_download_azevedo_array.sbatch"

[[ "$(head -n 1 "$script")" == "#!/bin/bash" ]] || {
  echo "ERROR: malformed $script" >&2; exit 2;
}
[[ -s "$manifest" ]] || {
  echo "ERROR: missing $manifest; run scripts/14_prepare_azevedo_spades_manifest.py first" >&2; exit 2;
}

n=$(( $(wc -l < "$manifest") - 1 ))
(( n > 0 )) || { echo "ERROR: empty manifest" >&2; exit 2; }

unset SBATCH_MEM_PER_NODE SBATCH_MEM_PER_CPU SBATCH_MEM_PER_GPU

echo "Submitting $n download tasks (max 12 concurrent)"
exec sbatch \
  --account=cooperma \
  --partition=standard \
  --cpus-per-task=1 \
  --mem-per-cpu=4G \
  --array="0-$((n-1))%12" \
  "$script"
