# Metasieve

Mining environmental metagenomes for deeply novel biological sequence and structure.

This is a metagenomic bioprospecting pipeline designed to identify
environmental DNA and protein sequences that remain unexplained after
taxonomic, functional, and structural database searches.

The orchestrator is a standalone Python package. External tools
(GGCAT, Kraken2, SeqScreen) are executed on the host through
`subprocess.run`. Kraken2 writes unclassified unitigs with
`--unclassified-out`; that FASTA is the SeqScreen input.

## Team

Dr. Todd Treangen\
Dr. Jennifer Lu\
Hiba Ben Aribi\
Felix Quintana\
Natalie Kokroko\
Jingyue Wu

## Workflow overview

![](metasieve_overview.png)

Raw paired-end short reads are assembled into unitigs and classified with
Kraken2. Unclassified sequences come from Kraken2 `--unclassified-out` and
are screened with SeqScreen (fast mode). Only queries that SeqScreen did
not assign a taxid **and** did not hit UniRef/UniProt are six-frame
translated. Those protein sequences are folded with Hugging Face ESMFold.
Every PDB is named and indexed so it traces back to the originating unitig
and ORF.

```
paired-end FASTQ
        │
        ▼
   GGCAT               Step 1  unitig assembly
        │ unitigs.fasta
        ▼
     Kraken2           Step 2  taxonomy
        │  --unclassified-out
        ▼
 unclassified FASTA    Step 3  optional min-length filter
        │
        ▼
    SeqScreen --mode fast   Step 4  functional screening
        │  keep queries with no taxid and no UniRef/UniProt hit
        ▼
   six-frame ORFs on unexplained unitigs
        │
        ▼
   ORF_manifest.csv
        │
        ▼
   Report generation
   
   
   (WIP: Future work) Hugging Face ESMFold    Step 5  structure prediction
        │  facebook/esmfold_v1
        ▼
   05_structures/ + structures_manifest.csv
```

## Software requirements

| Stage           | Tool | Notes |
|-----------------| --- | --- |
| Assembly        | GGCAT | Paired-end short reads → maximal unitigs (`unitigs.fasta`) |
| Taxonomy        | Kraken2 | `--unclassified-out` FASTA is fed to SeqScreen |
| Functions       | SeqScreen | Fast mode; keep queries with no taxid and no UniRef/UniProt hit |
| (WIP) Structure | Hugging Face ESMFold | `transformers` + PyTorch; default `facebook/esmfold_v1` |
| Orchestration   | Python ≥ 3.10 | Host binaries on PATH (`ggcat`, `kraken2`, `seqscreen`) |

Structural module is a planned downstream search step
and is not run by this pipeline yet.

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
│       ├── folding.py               # Hugging Face ESMFold (WIP)
│       ├── parsers.py               # FASTA / SeqScreen / manifests
│       ├── samples.py               # paired-end FASTQ discovery
│       ├── logging_setup.py
│       └── exceptions.py
├── tests/
├── LICENSE
└── README.md
```

## Install

```bash
conda env create -f environment.yml
conda activate metasieve
```

## Quick start

```bash
python main.py \
    --reads 'data/*_{1,2}.fastq.gz' \
    --kraken_db /data/dbs/kraken2 \
    --seqscreen_db /data/dbs/seqscreen_databases \
    --threads 16 \
    --outdir results
```

Or copy `assets/params.example.yml` and run:

```bash
python main.py --config params.yml
```

To skip earlier steps, pass the intermediate FASTA you already have:

```bash
# Start at SeqScreen (Kraken2 unclassified contigs)
python main.py \
    --unclassified sample.unclassified.fasta \
    --seqscreen_db /data/dbs/seqscreen_databases \
    --sample_id sample \
    --outdir results

# Start at Kraken2 (assembled unitigs)
python main.py \
    --unitigs sample.unitigs.fasta \
    --kraken_db /data/dbs/kraken2 \
    --seqscreen_db /data/dbs/seqscreen_databases \
    --outdir results

```

`--start_from` can force the entry point (`assembly`, `classification`,
`unclassified`, `folding`); by default it is inferred from the input flags.

GGCAT, Kraken2, and SeqScreen must be on `PATH` for the steps you run.
If you start from `--unclassified`, only SeqScreen is required. ESMFold runs in
the Python process on CUDA when available.

## Parameters

| Parameter | Role |
| --- | --- |
| `--reads` | Paired-end glob, e.g. `data/*_{1,2}.fastq.gz` |
| `--r1` / `--r2` | Explicit mates for a single sample |
| `--unitigs` | Assembled unitigs FASTA; skip GGCAT |
| `--unclassified` | Kraken2 unclassified FASTA; skip GGCAT and Kraken2 |
| `--orfs` | ORF protein FASTA or split directory; skip to ESMFold |
| `--start_from` | `auto` (default), `assembly`, `classification`, `unclassified`, `folding` |
| `--sample_id` | Sample name (useful with a single restart FASTA) |
| `--outdir` | Results directory |
| `--kraken_db` | Kraken2 database directory |
| `--seqscreen_db` | SeqScreen databases directory |
| `--threads` | CPU threads for GGCAT / Kraken2 / SeqScreen |
| `--memory_gb` | GGCAT memory cap in GB (passed to `ggcat build -m`) |
| `--ggcat_kmer` | GGCAT k-mer length (default: 31) |
| `--ggcat_min_multiplicity` | GGCAT minimum k-mer multiplicity `-s` (default: 2) |
| `--min_contig_len`, `--min_orf_aa`, `--max_orf_aa` | Length filters (unitigs after GGCAT, then unclassified FASTA) |
| `--require_orf_start` / `--no_require_orf_start` | Start-codon filter (default: require) |

Kraken2 `--unclassified-out` is the SeqScreen input. `--min_contig_len`
can drop short unclassified unitigs before screening. After SeqScreen,
only queries with no taxid **and** no UniRef/UniProt protein hit are
six-frame translated for folding. After that step, `{sample}.contig_report.csv`
lists every screened contig and `{sample}.seqscreen_step_report.txt` lists
kept contigs plus found ORFs (protein sequences are also in
`{sample}.orf_manifest.csv` and `{sample}.orfs.faa`).

## Outputs

```
results/
├── 01_assembly/<sample>/
├── 02_classification/<sample>/          # kraken2.out, report, unclassified FASTA
├── 03_filtering/<sample>/               # length-filtered unclassified FASTA
├── 04_seqscreen/<sample>/               # SeqScreen, unexplained FASTA, ORFs
│   ├── <sample>.contig_report.csv       # every screened contig (kept vs skipped)
│   ├── <sample>.orf_manifest.csv        # found ORFs (coords + protein sequence)
│   ├── <sample>.seqscreen_step_report.txt
│   └── orfs_split/<contig_id>/          # one directory per contig for split ORFs
└── pipeline_info/
    ├── metasieve.log
    └── config.json
```



## License

MIT. See [LICENSE](LICENSE).
