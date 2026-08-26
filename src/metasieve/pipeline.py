"""Five-step Metasieve orchestrator.

1. Assembly (metaSPAdes)
2. Taxonomic classification (Kraken2, with --unclassified-out)
3. Optional length filter on Kraken2 unclassified FASTA
4. Functional screening (SeqScreen fast) + ORF extraction
5. Structure prediction (Hugging Face ESMFold)
"""

from __future__ import annotations

import json
import logging
import shutil
from pathlib import Path

from metasieve.config import PipelineConfig
from metasieve.exceptions import PipelineError, ToolError
from metasieve.folding import ESMFoldPredictor
from metasieve.parsers import (
    build_structure_manifest,
    extract_seqscreen_orfs,
    filter_fasta_min_length,
)
from metasieve.samples import SampleReads, discover_samples
from metasieve.wrappers import ToolRunner, run_with_retries

LOGGER = logging.getLogger(__name__)


class Pipeline:
    def __init__(self, config: PipelineConfig) -> None:
        self.config = config
        self.tools = ToolRunner()
        self.outdir = config.resolved_outdir()

    def run(self) -> int:
        samples = discover_samples(
            reads=self.config.reads,
            r1=self.config.r1,
            r2=self.config.r2,
            sample_id=self.config.sample_id,
        )
        LOGGER.info("Samples: %s", ", ".join(s.sample_id for s in samples))
        LOGGER.info("Output directory: %s", self.outdir)

        info_dir = self.outdir / "pipeline_info"
        info_dir.mkdir(parents=True, exist_ok=True)
        (info_dir / "config.json").write_text(
            json.dumps(self.config.to_dict(), indent=2, default=str) + "\n",
            encoding="utf-8",
        )

        failures = 0
        for sample in samples:
            try:
                self.run_sample(sample)
            except (PipelineError, ToolError, OSError) as exc:
                failures += 1
                LOGGER.error("Sample %s failed: %s", sample.sample_id, exc)

        if failures:
            LOGGER.error("%d / %d samples failed", failures, len(samples))
            return 1
        LOGGER.info("Pipeline complete. Results: %s", self.outdir)
        LOGGER.info("Structures: %s", self.outdir / "05_structures")
        return 0

    def run_sample(self, sample: SampleReads) -> None:
        LOGGER.info("======== Sample %s ========", sample.sample_id)
        prefix = sample.sample_id

        assembly_dir = self._step_dir("01_assembly", prefix)
        class_dir = self._step_dir("02_classification", prefix)
        filter_dir = self._step_dir("03_filtering", prefix)
        seq_dir = self._step_dir("04_seqscreen", prefix)
        struct_root = self.outdir / "05_structures"

        # Step 1 — Assembly
        LOGGER.info("[1/5] Assembly (metaSPAdes)")
        contigs = self.tools.metaspades(
            r1=sample.r1,
            r2=sample.r2,
            outdir=assembly_dir,
            threads=self.config.threads,
            memory_gb=self.config.memory_gb,
        )
        published_contigs = assembly_dir / f"{prefix}.contigs.fasta"
        if contigs != published_contigs:
            shutil.copy2(contigs, published_contigs)
            contigs = published_contigs

        # Step 2 — Taxonomy (unclassified FASTA comes from Kraken2)
        LOGGER.info("[2/5] Taxonomic classification (Kraken2 --unclassified-out)")
        _assignments, _report, unclassified_raw = self.tools.kraken2(
            contigs=contigs,
            db=self.config.kraken_db,  # type: ignore[arg-type]
            outdir=class_dir,
            threads=self.config.threads,
            confidence=self.config.kraken2_confidence,
            prefix=prefix,
        )

        # Step 3 — Optional length filter, then hand off to SeqScreen
        LOGGER.info("[3/5] Preparing unclassified contigs for SeqScreen")
        unclassified_fasta = filter_dir / f"{prefix}.unclassified.fasta"
        n_unclassified = filter_fasta_min_length(
            unclassified_raw,
            unclassified_fasta,
            min_length=self.config.min_contig_len,
        )
        if n_unclassified == 0:
            LOGGER.warning("Sample %s: no unclassified contigs; skipping remaining steps", prefix)
            return

        # Step 4 — SeqScreen + ORF extraction
        LOGGER.info("[4/5] Functional screening (SeqScreen --mode %s)", self.config.seqscreen_mode)
        seqscreen_work = self.tools.seqscreen(
            fasta=unclassified_fasta,
            db=self.config.seqscreen_db,  # type: ignore[arg-type]
            workdir=seq_dir,
            threads=self.config.threads,
            mode=self.config.seqscreen_mode,
            extra_args=self.config.seqscreen_extra_args,
        )
        LOGGER.info("[4/5] Extracting translated ORFs")
        orfs = extract_seqscreen_orfs(
            sample_id=prefix,
            unclassified_fasta=unclassified_fasta,
            seqscreen_dir=seqscreen_work,
            output_fasta=seq_dir / f"{prefix}.orfs.faa",
            output_manifest=seq_dir / f"{prefix}.orf_manifest.csv",
            output_split_dir=seq_dir / "orfs_split",
            output_stats=seq_dir / f"{prefix}.orf_stats.json",
            min_aa=self.config.min_orf_aa,
            max_aa=self.config.max_orf_aa,
            require_start=self.config.require_orf_start,
            genetic_code=self.config.genetic_code,
        )
        if orfs.n_orfs == 0:
            LOGGER.warning("Sample %s: no ORFs; skipping structure prediction", prefix)
            return

        # Step 5 — Hugging Face ESMFold
        LOGGER.info("[5/5] Structure prediction (Hugging Face ESMFold, %d ORFs)", orfs.n_orfs)
        sidecars, structures = self._fold_orfs(prefix, orfs.split_fastas, struct_root)

        tracked_dir = struct_root / "tracked_structures"
        manifest = struct_root / f"{prefix}.structures_manifest.csv"
        build_structure_manifest(
            orf_manifest=orfs.manifest,
            output_manifest=manifest,
            output_dir=tracked_dir,
            sidecar_paths=sidecars,
            structure_paths=structures,
        )
        LOGGER.info("Tracking CSV: %s", manifest)

    def _fold_orfs(
        self,
        sample_id: str,
        orf_fastas: list[Path],
        struct_root: Path,
    ) -> tuple[list[Path], list[Path]]:
        cfg = self.config
        sidecars: list[Path] = []
        structures: list[Path] = []
        if cfg.skip_esmfold:
            LOGGER.info("Skipping ESMFold (--skip_esmfold)")
            return sidecars, structures

        esm_dir = struct_root / "esmfold" / sample_id
        predictor = ESMFoldPredictor.from_config(cfg)
        try:
            predictor.load()
        except PipelineError:
            raise
        except Exception as exc:
            raise PipelineError(
                f"Failed to load Hugging Face ESMFold ({cfg.esmfold_model}): {exc}"
            ) from exc

        for fasta in orf_fastas:
            record = fasta.stem
            LOGGER.info("Folding %s", record)
            result = self._try_fold(
                label=f"ESMFold {record}",
                func=lambda f=fasta: predictor.predict_fasta(f, esm_dir, sample_id),
            )
            if result:
                pdb, sidecar = result
                structures.append(pdb)
                sidecars.append(sidecar)

        return sidecars, structures

    def _try_fold(self, *, label: str, func):
        try:
            return run_with_retries(func, retries=self.config.fold_retries, label=label)
        except ToolError as exc:
            LOGGER.error("%s skipped after failures: %s", label, exc)
            return None

    def _step_dir(self, step: str, sample_id: str) -> Path:
        path = self.outdir / step / sample_id
        if path.exists() and self.config.force:
            shutil.rmtree(path)
        path.mkdir(parents=True, exist_ok=True)
        return path
