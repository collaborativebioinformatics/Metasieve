# Metasieve

Mining environmental metagenomes for deeply novel biological sequence and structure.

This is a metagenomic bioprospecting pipeline designed to identify
environmental DNA and protein sequences that remain unexplained after
taxonomic, functional, and structural database searches.

The orchestrator is a standalone Python package. External tools
(metaSPAdes, Kraken2, SeqScreen) are executed through `subprocess.run`
and may run on the host (`local`) or inside Docker / Singularity with
identity bind-mounts for databases. Structure prediction uses
**Hugging Face ESMFold** (`facebook/esmfold_v1`) in-process via
`transformers` and PyTorch.

## Team

Dr. Todd Treangen\
Dr. Jennifer Lu\
Hiba Ben Aribi\
Francesco Picchi\
Felix Quintana\
Natalie Kokroko\
Jingyue Wu

## Workflow overview

![](metasieve_overview.png)

Raw paired-end short reads are assembled, classified, and stripped down to
Kraken2-unclassified contigs. SeqScreen (fast mode) screens those contigs;
translated ORFs are then folded with Hugging Face ESMFold. Every PDB is
named and indexed so it traces back to the originating contig and ORF.

```
paired-end FASTQ
        │
        ▼
   metaSPAdes          Step 1  assembly
        │ contigs.fasta
        ▼
     Kraken2           Step 2  taxonomy
        │
        ▼
 unclassified FASTA    Step 3  Kraken2 'U' / taxid 0
        │
        ▼
    SeqScreen --mode fast   Step 4  functional screening
        │
        ▼
   ORF / translation FASTA
        │
        ▼
   Hugging Face ESMFold    Step 5  structure prediction
        │  facebook/esmfold_v1
        ▼
   05_structures/ + structures_manifest.csv
```

## Software requirements

| Stage | Tool | Notes |
| --- | --- | --- |
| Assembly | metaSPAdes (SPAdes) | Paired-end short reads |
| Taxonomy | Kraken2 | Local database via `--kraken_db` |
| Functions | SeqScreen | Fast mode; local `--seqscreen_db` |
| Structure | Hugging Face ESMFold | `transformers` + PyTorch; default `facebook/esmfold_v1` |
| Orchestration | Python ≥ 3.10 | Docker, Singularity/Apptainer, or host binaries |

Foldseek is a planned downstream search step against the predicted structures
and is not run by this pipeline yet.

Install a CUDA build of PyTorch for GPU folding, for example:

```bash
pip install torch --index-url https://download.pytorch.org/whl/cu124
pip install -e .
```

The first ESMFold call downloads `facebook/esmfold_v1` into the Hugging Face
cache (override with `--esmfold_cache`, or pass a local snapshot as
`--esmfold_model /path/to/snapshot`).

## Project layout

```
Metasieve/
├── main.py                          # CLI entry (adds src/ to sys.path)
├── pyproject.toml
├── requirements.txt
├── environment.yml
├── assets/
│   └── params.example.yml
├── src/
│   └── metasieve/
│       ├── __init__.py
│       ├── __main__.py              # python -m metasieve
│       ├── cli.py                   # argparse
│       ├── config.py                # dataclass + YAML overlay
│       ├── pipeline.py              # 5-step orchestrator
│       ├── wrappers.py              # subprocess CLI wrappers
│       ├── folding.py               # Hugging Face ESMFold
│       ├── parsers.py               # FASTA / Kraken / SeqScreen / manifests
│       ├── containers.py            # Docker / Singularity / local
│       ├── samples.py               # paired-end FASTQ discovery
│       ├── logging_setup.py
│       └── exceptions.py
├── tests/
├── LICENSE
└── README.md
```

## Install

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

Or with conda:

```bash
conda env create -f environment.yml
conda activate metasieve
```

You can also run without installing:

```bash
pip install -r requirements.txt
python main.py --help
```

## Quick start

```bash
python main.py \
    --container_engine singularity \
    --reads 'data/*_{1,2}.fastq.gz' \
    --kraken_db /data/dbs/kraken2 \
    --seqscreen_db /data/dbs/seqscreen_databases \
    --esmfold_model facebook/esmfold_v1 \
    --esmfold_cache /data/models/huggingface \
    --esmfold_device auto \
    --bind /data/dbs \
    --threads 16 \
    --outdir results
```

Or copy `assets/params.example.yml` and run:

```bash
python main.py --config params.yml --container_engine singularity
```

`--container_engine` choices: `docker`, `singularity`, `local`.
Apptainer is used automatically if `singularity` is not on PATH.
That flag applies to metaSPAdes, Kraken2, and SeqScreen. ESMFold runs in
the Python process on CUDA when available.

## Parameters

| Parameter | Role |
| --- | --- |
| `--reads` | Paired-end glob, e.g. `data/*_{1,2}.fastq.gz` |
| `--r1` / `--r2` | Explicit mates for a single sample |
| `--outdir` | Results directory |
| `--kraken_db` | Kraken2 database directory |
| `--seqscreen_db` | SeqScreen databases directory |
| `--esmfold_model` | Hugging Face id or local snapshot (`facebook/esmfold_v1`) |
| `--esmfold_cache` | Hugging Face cache directory for downloaded weights |
| `--esmfold_device` | `auto`, `cuda`, or `cpu` |
| `--esmfold_num_recycles` | ESMFold recycle iterations (default: 4) |
| `--esmfold_chunk_size` | Trunk chunk size to limit GPU memory (default: 128) |
| `--container_engine` | `docker`, `singularity`, or `local` |
| `--bind` | Extra host paths to identity-mount, comma-separated |
| `--threads` | CPU threads for assembly / Kraken2 / SeqScreen |
| `--memory_gb` | metaSPAdes `-m` (GB) |
| `--min_contig_len`, `--min_orf_aa`, `--max_orf_aa` | Length filters |
| `--require_orf_start` / `--no_require_orf_start` | Start-codon filter (default: require) |
| `--skip_esmfold` | Stop after ORF extraction |

Large Kraken2 / SeqScreen databases are bind-mounted at the same absolute
path inside the container (identity mounts). They are never copied into
the work directory.

## Outputs

```
results/
├── 01_assembly/<sample>/
├── 02_classification/<sample>/          # *.kraken2.out, *.kraken2.report
├── 03_filtering/<sample>/               # unclassified FASTA + contig manifest
├── 04_seqscreen/<sample>/               # SeqScreen working dir, ORFs, orf_manifest.csv
├── 05_structures/
│   ├── esmfold/<sample>/
│   ├── tracked_structures/              # stamped PDB with contig metadata
│   └── <sample>.structures_manifest.csv # contig → ORF → PDB map
└── pipeline_info/
    ├── metasieve.log
    └── config.json
```

### Tracking (contig → structure)

ORF FASTA IDs are filesystem-safe and encode the source contig:

```
{sample}__{contig_id}__ORF_0001
```

The FASTA header keeps the original contig ID, coordinates, and SeqScreen
query. Predicted files are published as:

```
{sample}__{contig_id}__ORF_0001.esmfold.pdb
```

`structures_manifest.csv` joins each PDB to `sample_id`, `contig_id`,
`orf_id`, coordinates, SeqScreen query, mean pLDDT, and a peptide checksum.
PDB files also carry `REMARK 0 METASIEVE TRACKING` records with the same
fields.

## Container volume binding

Docker and Singularity do not see host databases unless those directories
are mounted. Metasieve constructs identity binds so a host path such as
`/data/dbs/kraken2` is available at the same path inside the container:

```text
docker run --rm \
    -u $(id -u):$(id -g) \
    -v /data/dbs/kraken2:/data/dbs/kraken2:ro \
    -v /path/to/results/sample:/path/to/results/sample \
    -w /path/to/results/sample \
    IMAGE kraken2 --db /data/dbs/kraken2 ...

singularity exec \
    -B /data/dbs/kraken2:/data/dbs/kraken2 \
    --pwd /path/to/results/sample \
    docker://IMAGE kraken2 --db /data/dbs/kraken2 ...
```

`--bind` adds extra host roots (for example a shared `/data` filesystem).
Use `--container_engine local` only when every binary is already on PATH.

ESMFold loads once per sample, then folds each ORF. A failed fold is
retried, then skipped so a single GPU OOM does not abort the sample.

## License

MIT. See [LICENSE](LICENSE).
