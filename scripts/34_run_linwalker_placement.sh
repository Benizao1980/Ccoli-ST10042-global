#!/bin/bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <PubMLST_Genome_Comparator.xlsx> [outdir]" >&2
  exit 2
fi

xlsx="$1"
outdir="${2:-results/linwalker_st10042}"
anchors="results/pubmlst_genome_comparator/current_pubmlst_st10042_anchors.tsv"

[[ -s "$xlsx" ]] || { echo "ERROR missing Genome Comparator workbook: $xlsx" >&2; exit 2; }
[[ -s "$anchors" ]] || {
  echo "ERROR missing anchor metadata: $anchors" >&2
  echo "Run: python scripts/32_package_pubmlst_genome_comparator.py" >&2
  exit 2
}

python - <<'PY'
import linwalker
print(f"LINwalker {linwalker.__version__}")
PY

python -m linwalker place \
  --profiles "$xlsx" \
  --sheet all \
  --reference-metadata "$anchors" \
  --reference-id-col pubmlst_id \
  --lin-col LINcode_v2 \
  --cgst-col cgST_v2 \
  --st-col ST \
  --cc-col clonal_complex \
  --query-prefix AZE_ \
  --outdir "$outdir"

python scripts/33_summarise_linwalker_placement.py \
  --placement "$outdir/tables/placement_summary.tsv"

echo
echo "Primary project table:"
echo "  results/st10042_v2_lin_placement.tsv"
echo
echo "Before biological interpretation, inspect:"
echo "  $outdir/tables/reference_label_map.tsv"
echo "All 40 PubMLST anchor IDs should map unambiguously where present in the Genome Comparator run."
