#!/bin/bash
#SBATCH --job-name=st10042_pubmlst
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=st10042_pubmlst-%j.out

# Do not enable nounset until after the site/login bash configuration has loaded.
# Puma's /etc/bashrc references BASHRCSOURCED before defining it, which fails under set -u.
set -eo pipefail
source ~/.bashrc
set -u

conda activate st10042

: "${PUBMLST_API_KEY:?Set PUBMLST_API_KEY before submitting the job}"

# Current global ST10042 population + all available PubMLST assemblies.
python scripts/01_query_pubmlst_st10042.py \
  --out data/pubmlst_st10042_all.tsv \
  --json-out data/pubmlst_st10042_all.json \
  --download-contigs data/pubmlst_st10042_assemblies

# Focal Peru metadata subset (assemblies are already present above).
python scripts/01_query_pubmlst_st10042.py \
  --country Peru \
  --out data/peru_st10042_current.tsv \
  --json-out data/peru_st10042_current.json

echo
echo "Global ST10042 records:"
awk 'END{print NR-1}' data/pubmlst_st10042_all.tsv
echo "Downloaded PubMLST assemblies:"
find data/pubmlst_st10042_assemblies -name '*.fasta' | wc -l
echo "Peru ST10042 records:"
awk 'END{print NR-1}' data/peru_st10042_current.tsv
