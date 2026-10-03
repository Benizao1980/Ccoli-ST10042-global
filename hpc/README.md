# HPC quick start

```bash
git clone git@github.com:Benizao1980/Ccoli-ST10042-global.git
cd Ccoli-ST10042-global

module load anaconda 2>/dev/null || true
conda env create -f environment.yml
conda activate st10042
```

## PubMLST authentication
PubMLST restricts access to records added after 31 Dec 2024 unless authenticated.
Create a personal API key from your PubMLST/BIGSdb profile and export it in the shell:

```bash
export PUBMLST_API_KEY='your-key-here'
```

Do **not** commit the key. The query script sends it as `X-API-Key` when present.

## Query all current ST10042

```bash
python scripts/01_query_pubmlst_st10042.py \
  --out data/pubmlst_st10042_all.tsv \
  --json-out data/pubmlst_st10042_all.json
```

## Query Peru and download PubMLST assemblies

```bash
python scripts/01_query_pubmlst_st10042.py \
  --country Peru \
  --out data/peru_st10042_current.tsv \
  --json-out data/peru_st10042_current.json \
  --download-contigs data/peru_assemblies
```

Expected current checkpoint from the exploratory PubMLST query: 23 Peru ST10042 records
(16 Lima, 7 Iquitos). Treat that as a check, not a hard-coded filter.
