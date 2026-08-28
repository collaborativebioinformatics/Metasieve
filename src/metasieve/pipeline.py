"""Five-step Metasieve orchestrator.

1. Assembly (GGCAT unitigs)
2. Taxonomic classification (Kraken2, with --unclassified-out)
3. Optional length filter on Kraken2 unclassified FASTA
4. Functional screening (SeqScreen fast) + ORF extraction from unexplained unitigs
5. Structure prediction (Hugging Face ESMFold)

Earlier steps can be skipped by supplying an intermediate FASTA
(``--unitigs``, ``--unclassified``, or ``--orfs``).
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
    prepare_orf_fastas,
    write_orf_manifest_from_fastas,
)
from metasieve.samples import (
    START_CLASSIFICATION,
    START_FOLDING,
    START_UNCLASSIFIED,
    SampleReads,
    discover_pipeline_samples,
)
from metasieve.wrappers import ToolRunner, run_with_retries

LOGGER = logging.getLogger(__name__)


class Pipeline:
    def __init__(self, config: PipelineConfig) -> None:
        self.config = config
        self.tools = ToolRunner()
        self.outdir = config.resolved_outdir()

    def run(self) -> int:
        samples = discover_pipeline_samples(
            start_from=self.config.start_from,
            reads=self.config.reads,
            r1=self.config.r1,
            r2=self.config.r2,
            sample_id=self.config.sample_id,
            unitigs=self.config.unitigs,
            unclassified=self.config.unclassified,
            orfs=self.config.orfs,
        )
        LOGGER.info("Start from: %s", self.config.start_from)
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
        start = self.config.start_from
        struct_root = self.outdir / "05_structures"

        if start == START_FOLDING:
            self._run_folding_from_input(sample, struct_root)
            return

        unclassified_fasta = self._unclassified_fasta(sample, prefix, start)
        if unclassified_fasta is None:
            return

        seq_dir = self._step_dir("04_seqscreen", prefix)
        LOGGER.info("[4/5] Functional screening (SeqScreen --mode %s)", self.config.seqscreen_mode)
        seqscreen_work = self.tools.seqscreen(
            fasta=unclassified_fasta,
            db=self.config.seqscreen_db,  # type: ignore[arg-type]
            workdir=seq_dir,
            threads=self.config.threads,
            mode=self.config.seqscreen_mode,
            extra_args=self.config.seqscreen_extra_args,
        )
        LOGGER.info("[4/4] Extracting ORFs from SeqScreen-unexplained unitigs")
        orfs = extract_seqscreen_orfs(
            sample_id=prefix,
            unclassified_fasta=unclassified_fasta,
            seqscreen_dir=seqscreen_work,
            output_fasta=seq_dir / f"{prefix}.orfs.faa",
            output_manifest=seq_dir / f"{prefix}.orf_manifest.csv",
            output_split_dir=seq_dir / "orfs_split",
            output_stats=seq_dir / f"{prefix}.orf_stats.json",
            output_unexplained_fasta=seq_dir / f"{prefix}.unexplained.fasta",
            output_contig_report=seq_dir / f"{prefix}.contig_report.csv",
            output_step_report=seq_dir / f"{prefix}.seqscreen_step_report.txt",
            min_aa=self.config.min_orf_aa,
            max_aa=self.config.max_orf_aa,
            require_start=self.config.require_orf_start,
            genetic_code=self.config.genetic_code,
        )
        LOGGER.info("Step 4 contig report: %s", orfs.contig_report)
        LOGGER.info("Step 4 ORF/contig summary: %s", orfs.step_report)
        if orfs.n_unexplained_contigs == 0:
            LOGGER.warning(
                "Sample %s: no SeqScreen-unexplained unitigs (no taxid and no UniRef hit); skipping remaining steps",
                prefix,
            )
            return
        if orfs.n_orfs == 0:
            LOGGER.warning("Sample %s: no ORFs; skipping structure prediction", prefix)
            return

        #self._fold_and_track(prefix, orfs.split_fastas, orfs.manifest, struct_root)

    def _unclassified_fasta(
        self, sample: SampleReads, prefix: str, start: str
    ) -> Path | None:
        if start == START_UNCLASSIFIED:
            source = sample.unclassified
            if source is None:
                raise PipelineError(f"Sample {prefix} has no --unclassified FASTA")
            LOGGER.info("[1/5] Assembly skipped (starting from unclassified FASTA)")
            LOGGER.info("[2/5] Kraken2 skipped (starting from unclassified FASTA)")
            return self._filter_unclassified(source, prefix)

        unitigs = self._unitigs_fasta(sample, prefix, start)
        LOGGER.info("[2/5] Taxonomic classification (Kraken2 --unclassified-out)")
        class_dir = self._step_dir("02_classification", prefix)
        _assignments, _report, unclassified_raw = self.tools.kraken2(
            contigs=unitigs,
            db=self.config.kraken_db,  # type: ignore[arg-type]
            outdir=class_dir,
            threads=self.config.threads,
            confidence=self.config.kraken2_confidence,
            prefix=prefix,
        )
        return self._filter_unclassified(unclassified_raw, prefix)

    def _unitigs_fasta(self, sample: SampleReads, prefix: str, start: str) -> Path:
        assembly_dir = self._step_dir("01_assembly", prefix)
        published = assembly_dir / f"{prefix}.unitigs.fasta"
        if start == START_CLASSIFICATION:
            source = sample.unitigs
            if source is None:
                raise PipelineError(f"Sample {prefix} has no --unitigs FASTA")
            LOGGER.info("[1/5] Assembly skipped (starting from unitigs FASTA)")
            if source.resolve() != published.resolve():
                shutil.copy2(source, published)
            return published

        LOGGER.info("[1/5] Assembly (GGCAT unitigs)")
        if sample.r1 is None or sample.r2 is None:
            raise PipelineError(f"Sample {prefix} is missing paired reads")
        unitigs = self.tools.ggcat(
            r1=sample.r1,
            r2=sample.r2,
            outdir=assembly_dir,
            threads=self.config.threads,
            memory_gb=self.config.memory_gb,
            kmer=self.config.ggcat_kmer,
            min_multiplicity=self.config.ggcat_min_multiplicity,
            min_unitig_len=self.config.min_contig_len,
            force=self.config.force,
        )
        if unitigs != published:
            shutil.copy2(unitigs, published)
        return published

    def _filter_unclassified(self, source: Path, prefix: str) -> Path | None:
        LOGGER.info("[3/5] Preparing unclassified unitigs for SeqScreen")
        filter_dir = self._step_dir("03_filtering", prefix)
        unclassified_fasta = filter_dir / f"{prefix}.unclassified.fasta"
        n_unclassified = filter_fasta_min_length(
            source,
            unclassified_fasta,
            min_length=self.config.min_contig_len,
        )
        if n_unclassified == 0:
            LOGGER.warning(
                "Sample %s: no unclassified unitigs; skipping remaining steps", prefix
            )
            return None
        return unclassified_fasta

    def _run_folding_from_input(self, sample: SampleReads, struct_root: Path) -> None:
        prefix = sample.sample_id
        source = sample.orfs
        if source is None:
            raise PipelineError(f"Sample {prefix} has no --orfs input")
        LOGGER.info("[1/5] Assembly skipped (starting from ORFs)")
        LOGGER.info("[2/5] Kraken2 skipped (starting from ORFs)")
        LOGGER.info("[3/5] Length filter skipped (starting from ORFs)")
        LOGGER.info("[4/5] SeqScreen skipped (starting from ORFs)")
        seq_dir = self._step_dir("04_seqscreen", prefix)
        split_fastas = prepare_orf_fastas(source, seq_dir / "orfs_split")
        manifest = seq_dir / f"{prefix}.orf_manifest.csv"
        sibling = source if source.is_file() else None
        if sibling is not None:
            candidate = sibling.with_name(f"{prefix}.orf_manifest.csv")
            if candidate.is_file() and candidate.resolve() != manifest.resolve():
                shutil.copy2(candidate, manifest)
        if not manifest.is_file():
            write_orf_manifest_from_fastas(
                sample_id=prefix,
                split_fastas=split_fastas,
                output_manifest=manifest,
            )
        if not split_fastas:
            LOGGER.warning("Sample %s: no ORF FASTA files; skipping structure prediction", prefix)
            return
        self._fold_and_track(prefix, split_fastas, manifest, struct_root)

    def _fold_and_track(
        self,
        prefix: str,
        split_fastas: list[Path],
        orf_manifest: Path,
        struct_root: Path,
    ) -> None:
        LOGGER.info("[5/5] Structure prediction (Hugging Face ESMFold, %d ORFs)", len(split_fastas))
        sidecars, structures = self._fold_orfs(prefix, split_fastas, struct_root)
        tracked_dir = struct_root / "tracked_structures"
        manifest = struct_root / f"{prefix}.structures_manifest.csv"
        build_structure_manifest(
            orf_manifest=orf_manifest,
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
