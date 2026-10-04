# HPC quick start

```bash
git clone git@github.com:Benizao1980/Ccoli-ST10042-global.git
cd Ccoli-ST10042-global

module load anaconda 2>/dev/null || true
conda env create -f environment.yml
conda activate st10042
```

For an existing environment after repository updates:

```bash
conda env update -f environment.yml
```

## PubMLST authentication: OAuth, not a personal API key

The ST10042 isolate search uses the BIGSdb `/isolates/search` endpoint, which is an
HTTP POST request. PubMLST does not permit personal `X-API-Key` credentials for POST
requests. The full contemporary dataset therefore uses OAuth via the official
`bigsdb-downloader` utility.

### One-time setup on the login node

In your PubMLST/BIGSdb profile, obtain an **OAuth client key and client secret**
(they are distinct from the simple personal API key). Then run:

```bash
mkdir -p ~/.bigsdb_tokens

bigsdb-downloader \
  --key_name PubMLST \
  --site PubMLST \
  --db pubmlst_campylobacter_isolates \
  --token_dir ~/.bigsdb_tokens \
  --setup
```

The program will ask for the client key and secret, give you an authorization URL,
and then ask for the verifier code shown after you authorize the client in PubMLST.
It stores the resulting access token under `~/.bigsdb_tokens`. Session tokens are
renewed automatically.

Do not place keys, secrets or token files in this repository.

## Recommended first job

After OAuth setup succeeds:

```bash
bash hpc/submit_pubmlst.sh
```

The wrapper uses the known-good Puma configuration for this account:
`--account=cooperma --partition=standard`, with memory requested as
`--mem-per-cpu=4G`. On Puma, this is more reliable than `--mem` for this
allocation and avoids falling onto the windfall partition.

This creates:

- `data/pubmlst_st10042_all.tsv` — all current authenticated ST10042 metadata
- `data/pubmlst_st10042_all.json` — full selected PubMLST records
- `data/pubmlst_st10042_assemblies/` — retrievable PubMLST ST10042 contig FASTAs
- `data/peru_st10042_current.tsv` — focal Peru subset
- `data/peru_st10042_current.json` — full Peru records

Anonymous mode remains available for diagnostics:

```bash
python scripts/01_query_pubmlst_st10042.py \
  --auth anonymous \
  --out data/pubmlst_st10042_public.tsv \
  --json-out data/pubmlst_st10042_public.json
```

Do not use the anonymous result as the final global dataset because contemporary
records may be hidden by PubMLST's access policy.

## Resolve the Azevedo European reads

```bash
sbatch hpc/02_resolve_ena.sbatch
```

That job resolves the 217 published read accessions. We can download/assemble them
after checking that the manifest reproduces the paper before moving to cgMLST.

## Check for inherited Slurm overrides

If a job reports an unexpected `ReqMem`, inspect the submitting shell:

```bash
env | grep '^SBATCH_' || true
```

On Puma, one standard CPU corresponds to 5 GB memory. The safe submission wrapper above explicitly requests 5 GB and removes inherited memory overrides before calling `sbatch`.


## Known-good Puma settings

For the cooperma allocation, use:

```bash
#SBATCH --account=cooperma
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
```

Submit with:

```bash
bash hpc/submit_pubmlst.sh
```


## Download and assemble the accessible Azevedo reads

Build the exact paired-read manifest:

```bash
python scripts/14_prepare_azevedo_spades_manifest.py
```

Create the separate assembly environment once:

```bash
conda env create -f environment-assembly.yml
```

Download all paired FASTQs with MD5 verification:

```bash
conda activate st10042
bash hpc/submit_azevedo_download.sh
```

After the download array completes cleanly, assemble with conservative fastp
preprocessing followed by SPAdes `--isolate`:

```bash
conda activate st10042-assembly
bash hpc/submit_spades.sh
```

The submission wrappers calculate the array size from
`results/azevedo_spades_manifest.tsv`; they do not hard-code 110 samples.


## Coverage-filter the full Azevedo set

The nine-isolate validation pilot showed that the INNUca-like post-assembly
coverage filter removes duplicate low-depth assembly fragments from usable genomes
while leaving the genuinely poor UK-2/UK-4/UK-6 samples clearly poor.

Scale the validated filter to every paired Azevedo assembly:

```bash
conda activate st10042-assembly
bash hpc/submit_azevedo_filter.sh
```

The wrapper derives the array size from
`results/azevedo_spades_manifest.tsv` and runs at most 10 tasks concurrently.

After the array completes:

```bash
conda activate st10042
python scripts/23_summarise_azevedo_innuca_like.py
```

Filtered FASTAs are written under
`data/europe_assemblies_innuca_like_trimmed/<Strain_ID>/contigs.innuca_like.fasta`.
Each sample also has `qc.json` and `contig_coverage.tsv`. Mapping BAM/index files
are deleted after successful filtering to avoid retaining unnecessary intermediate
data.

Do **not** put PubMLST REST calls inside the Slurm array. Scheme-8 typing/QC is the
next, separate step and should be performed at controlled request concurrency.
