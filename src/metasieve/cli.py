"""Argparse CLI for the Metasieve Python pipeline."""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

from metasieve import __version__
from metasieve.config import (
    PipelineConfig,
    load_yaml_config,
    mapping_to_config,
    overlay_cli,
)
from metasieve.exceptions import PipelineError
from metasieve.logging_setup import setup_logging
from metasieve.pipeline import Pipeline

LOGGER = logging.getLogger(__name__)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="metasieve",
        description=(
            "Metasieve — assemble paired-end reads, keep Kraken2-unclassified "
            "contigs, screen with SeqScreen (fast), and fold ORFs with Hugging Face ESMFold."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="YAML parameter file (CLI flags override file values).",
    )

    io = parser.add_argument_group("inputs / outputs")
    io.add_argument(
        "--reads",
        default=None,
        help="Paired-end FASTQ glob, e.g. 'data/*_{1,2}.fastq.gz'.",
    )
    io.add_argument("--r1", type=Path, default=None, help="Forward reads (single sample).")
    io.add_argument("--r2", type=Path, default=None, help="Reverse reads (single sample).")
    io.add_argument("--sample_id", "--sample-id", dest="sample_id", default=None)
    io.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Pipeline results directory.",
    )

    dbs = parser.add_argument_group("databases and model weights")
    dbs.add_argument("--kraken_db", "--kraken-db", dest="kraken_db", type=Path, default=None)
    dbs.add_argument(
        "--seqscreen_db", "--seqscreen-db", dest="seqscreen_db", type=Path, default=None
    )
    dbs.add_argument(
        "--esmfold_model",
        "--esmfold-model",
        dest="esmfold_model",
        default=None,
        help="Hugging Face model id or local snapshot (default: facebook/esmfold_v1).",
    )
    dbs.add_argument(
        "--esmfold_cache",
        "--esmfold-cache",
        dest="esmfold_cache",
        type=Path,
        default=None,
        help="Hugging Face cache directory for ESMFold weights.",
    )
    dbs.add_argument(
        "--esmfold_device",
        "--esmfold-device",
        dest="esmfold_device",
        choices=["auto", "cuda", "cpu"],
        default=None,
        help="Device for ESMFold (default: auto).",
    )

    tools = parser.add_argument_group("tool behaviour")
    tools.add_argument("--seqscreen_mode", dest="seqscreen_mode", default=None)
    tools.add_argument("--seqscreen_extra_args", dest="seqscreen_extra_args", default=None)
    tools.add_argument("--kraken2_confidence", dest="kraken2_confidence", type=float, default=None)
    tools.add_argument("--min_contig_len", dest="min_contig_len", type=int, default=None)
    tools.add_argument("--min_orf_aa", dest="min_orf_aa", type=int, default=None)
    tools.add_argument("--max_orf_aa", dest="max_orf_aa", type=int, default=None)
    orf_start = tools.add_mutually_exclusive_group()
    orf_start.add_argument(
        "--require_orf_start",
        dest="require_orf_start",
        action="store_true",
        help="Keep only ORFs that begin with ATG/GTG/TTG.",
    )
    orf_start.add_argument(
        "--no_require_orf_start",
        dest="require_orf_start",
        action="store_false",
        help="Keep ORFs regardless of start codon.",
    )
    parser.set_defaults(require_orf_start=None)
    tools.add_argument("--genetic_code", dest="genetic_code", type=int, default=None)
    tools.add_argument("--esmfold_num_recycles", dest="esmfold_num_recycles", type=int, default=None)
    tools.add_argument("--esmfold_chunk_size", dest="esmfold_chunk_size", type=int, default=None)
    tools.add_argument("--skip_esmfold", dest="skip_esmfold", action="store_true", default=None)

    runtime = parser.add_argument_group("runtime")
    runtime.add_argument("--threads", type=int, default=None, help="CPU threads for assembly/taxonomy/SeqScreen.")
    runtime.add_argument("--memory_gb", "--memory-gb", dest="memory_gb", type=int, default=None)
    runtime.add_argument(
        "--log_level",
        "--log-level",
        dest="log_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default=None,
    )
    runtime.add_argument(
        "--force",
        action="store_true",
        default=None,
        help="Remove existing per-sample step directories before running.",
    )
    return parser


def _cli_mapping(args: argparse.Namespace) -> dict:
    mapping = vars(args).copy()
    mapping.pop("config", None)
    return mapping


def config_from_args(args: argparse.Namespace) -> PipelineConfig:
    cfg = PipelineConfig()
    if args.config is not None:
        yaml_path = Path(args.config).expanduser().resolve()
        if not yaml_path.is_file():
            raise PipelineError(f"Config file not found: {yaml_path}")
        cfg = mapping_to_config(load_yaml_config(yaml_path), base=cfg)
    cfg = overlay_cli(cfg, _cli_mapping(args))
    cfg.validate()
    return cfg


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    # Bootstrap logging before validation so errors are formatted.
    setup_logging(args.log_level or "INFO")
    try:
        cfg = config_from_args(args)
    except PipelineError as exc:
        LOGGER.error("%s", exc)
        parser.print_usage(sys.stderr)
        return 2

    log_file = cfg.resolved_outdir() / "pipeline_info" / "metasieve.log"
    setup_logging(cfg.log_level, log_file=log_file)
    LOGGER.info("Metasieve v%s", __version__)
    LOGGER.info("Log file: %s", log_file)

    try:
        return Pipeline(cfg).run()
    except PipelineError as exc:
        LOGGER.error("%s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
