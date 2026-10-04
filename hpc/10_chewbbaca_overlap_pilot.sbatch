#!/bin/bash
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --job-name=st10042_chewpilot
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --output=st10042_chewpilot-%j.out
#SBATCH --error=st10042_chewpilot-%j.err

set -euo pipefail

if [[ -z "${CONDA_EXE:-}" || ! -x "${CONDA_EXE}" ]]; then
  echo "ERROR: CONDA_EXE was not inherited. Submit from activated st10042-cgmlst." >&2
  exit 2
fi

eval "$("${CONDA_EXE}" shell.bash hook)"
conda activate st10042-cgmlst

SCHEMA="schema/chewbbaca_pubmlst_cgmlst_v2_scheme8_2026-10-04"
OUT="results/chewbbaca_overlap_pilot"

[[ -d "${SCHEMA}" ]] || { echo "ERROR missing adapted schema: ${SCHEMA}" >&2; exit 2; }

echo "=== Audit 163 PrepExternalSchema exclusions against current ST10042 cgSTs ==="
python scripts/28_audit_chewbbaca_schema_exclusions.py

echo
echo "=== Prepare nine known Azevedo/PubMLST controls ==="
python scripts/29_prepare_chewbbaca_overlap_pilot.py

if [[ -e "${OUT}" ]]; then
  echo "ERROR pilot output already exists: ${OUT}" >&2
  echo "Refusing to mix/reuse allele-call output." >&2
  exit 3
fi

echo
echo "=== chewBBACA AlleleCall pilot ==="
chewBBACA.py AlleleCall \
  -i results/chewbbaca_overlap_inputs.txt \
  -g "${SCHEMA}" \
  -o "${OUT}" \
  --cpu "${SLURM_CPUS_PER_TASK}" \
  --no-inferred \
  --output-masked \
  --output-novel \
  --output-missing

echo
echo "=== Validate PubMLST allele numbering ==="
python scripts/30_validate_chewbbaca_overlap_calls.py --pilot-dir "${OUT}"
