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

Do **not** commit the key. The query script sends it as `X-API-Key` on the search,
record and FASTA requests.

## Recommended first job: retrieve the current global ST10042 population

```bash
sbatch hpc/01_pubmlst.sbatch
```

This creates:

- `data/pubmlst_st10042_all.tsv` — all current ST10042 metadata returned by PubMLST
- `data/pubmlst_st10042_all.json` — full selected PubMLST records
- `data/pubmlst_st10042_assemblies/` — all retrievable ST10042 PubMLST contig FASTAs
- `data/peru_st10042_current.tsv` — focal Peru subset
- `data/peru_st10042_current.json` — full Peru records

The exploratory query used to start the project contained 23 confirmed Peru ST10042
records (16 Lima, 7 Iquitos). Treat that as a checkpoint rather than a hard-coded filter.

## Resolve the Azevedo European reads

```bash
sbatch hpc/02_resolve_ena.sbatch
```

That job only resolves the 217 published read accessions. We can download/assemble
them after checking that the manifest reproduces the paper before moving to cgMLST.
