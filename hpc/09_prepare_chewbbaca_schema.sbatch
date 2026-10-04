#!/bin/bash
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --job-name=st10042_prepschema
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH --time=08:00:00
#SBATCH --output=st10042_prepschema-%j.out
#SBATCH --error=st10042_prepschema-%j.err

set -euo pipefail

if [[ -z "${CONDA_EXE:-}" || ! -x "${CONDA_EXE}" ]]; then
  echo "ERROR: CONDA_EXE was not inherited. Submit from activated st10042-cgmlst." >&2
  exit 2
fi

eval "$("${CONDA_EXE}" shell.bash hook)"
conda activate st10042-cgmlst

SNAPSHOT="schema/pubmlst_campy_cgmlst_v2_scheme8_2026-10-04"
SOURCE="${SNAPSHOT}/alleles"
OUT="schema/chewbbaca_pubmlst_cgmlst_v2_scheme8_2026-10-04"
TRAIN_DIR="schema/training"
TRAIN="${TRAIN_DIR}/Ccoli_ST10042_UK7.prodigal.trn"
TRAIN_GENOME="data/europe_assemblies_innuca_like_trimmed/UK-7/contigs.innuca_like.fasta"

[[ -d "${SOURCE}" ]] || { echo "ERROR missing frozen allele directory: ${SOURCE}" >&2; exit 2; }
[[ -s "${TRAIN_GENOME}" ]] || { echo "ERROR missing UK-7 training genome: ${TRAIN_GENOME}" >&2; exit 2; }

mkdir -p "${TRAIN_DIR}"

if [[ ! -s "${TRAIN}" ]]; then
  echo "Training Prodigal on clean UK-7 assembly..."
  prodigal -i "${TRAIN_GENOME}" -t "${TRAIN}" -g 11 -p single -o /dev/null
else
  echo "Using existing Prodigal training file: ${TRAIN}"
fi

if [[ -d "${OUT}" ]]; then
  echo "ERROR adapted schema output already exists: ${OUT}" >&2
  echo "Refusing to overwrite. Remove it only after checking any previous run." >&2
  exit 3
fi

echo "Adapting frozen PubMLST scheme 8 for chewBBACA..."
chewBBACA.py PrepExternalSchema \
  -g "${SOURCE}" \
  -o "${OUT}" \
  --ptf "${TRAIN}" \
  --bsr 0.6 \
  --t 11 \
  --st 0.2 \
  --cpu "${SLURM_CPUS_PER_TASK}"

python scripts/27_validate_chewbbaca_schema.py --schema "${OUT}"
