"""CLI wrappers: ``subprocess.run`` around assembly, taxonomy, and SeqScreen.

Tools are executed directly on the host (no Docker / Singularity).
ESMFold runs in-process via Hugging Face (:mod:`metasieve.folding`).
"""

from __future__ import annotations

import logging
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from metasieve.exceptions import ToolError

LOGGER = logging.getLogger(__name__)
_TAIL = 20000


@dataclass
class CommandResult:
    argv: list[str]
    returncode: int
    stdout: str
    stderr: str
    stdout_log: Path | None = None
    stderr_log: Path | None = None


def run_command(
    argv: Sequence[str],
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
    log_dir: Path | None = None,
    log_prefix: str = "cmd",
) -> CommandResult:
    """Run *argv* with ``subprocess.run(check=True)`` and capture stdout/stderr."""
    argv = [str(part) for part in argv]
    rendered = shlex.join(argv)
    LOGGER.info("Running: %s", rendered)
    if cwd:
        LOGGER.debug("cwd=%s", cwd)

    try:
        completed = subprocess.run(
            argv,
            cwd=str(cwd) if cwd else None,
            env=env,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise ToolError(
            f"Executable not found: {argv[0]}",
            cmd=argv,
        ) from exc
    except subprocess.CalledProcessError as exc:
        stdout = exc.stdout or ""
        stderr = exc.stderr or ""
        stdout_log, stderr_log = _write_logs(log_dir, log_prefix, stdout, stderr)
        LOGGER.error("Command failed (exit %s): %s", exc.returncode, rendered)
        if stdout:
            LOGGER.error("stdout (tail):\n%s", stdout[-_TAIL:])
        if stderr:
            LOGGER.error("stderr (tail):\n%s", stderr[-_TAIL:])
        raise ToolError(
            f"Command failed with exit code {exc.returncode}: {argv[0]}",
            cmd=argv,
            returncode=exc.returncode,
            stdout=stdout,
            stderr=stderr,
        ) from exc

    stdout_log, stderr_log = _write_logs(
        log_dir, log_prefix, completed.stdout or "", completed.stderr or ""
    )
    if completed.stdout:
        LOGGER.debug("stdout (tail):\n%s", completed.stdout[-_TAIL:])
    if completed.stderr:
        LOGGER.debug("stderr (tail):\n%s", completed.stderr[-_TAIL:])
    return CommandResult(
        argv=argv,
        returncode=completed.returncode,
        stdout=completed.stdout or "",
        stderr=completed.stderr or "",
        stdout_log=stdout_log,
        stderr_log=stderr_log,
    )


def _write_logs(
    log_dir: Path | None,
    prefix: str,
    stdout: str,
    stderr: str,
) -> tuple[Path | None, Path | None]:
    if log_dir is None:
        return None, None
    log_dir.mkdir(parents=True, exist_ok=True)
    stdout_log = log_dir / f"{prefix}.stdout.log"
    stderr_log = log_dir / f"{prefix}.stderr.log"
    stdout_log.write_text(stdout, encoding="utf-8")
    stderr_log.write_text(stderr, encoding="utf-8")
    return stdout_log, stderr_log


def run_with_retries(func, *, retries: int, label: str):
    """Call *func* up to *retries* times; re-raise the last ToolError."""
    last: ToolError | None = None
    attempts = max(1, retries)
    for attempt in range(1, attempts + 1):
        try:
            return func()
        except ToolError as exc:
            last = exc
            LOGGER.warning("%s attempt %d/%d failed: %s", label, attempt, attempts, exc)
    assert last is not None
    LOGGER.error("%s failed after %d attempts", label, attempts)
    raise last


class ToolRunner:
    """Host-binary wrappers for each pipeline CLI tool."""

    def _run(
        self,
        argv: Sequence[str],
        *,
        workdir: Path,
        log_prefix: str = "cmd",
    ) -> CommandResult:
        return run_command(
            argv,
            cwd=workdir,
            log_dir=workdir,
            log_prefix=log_prefix,
        )

    # ── Step 1: metaSPAdes ────────────────────────────────────────────────

    def metaspades(
        self,
        *,
        r1: Path,
        r2: Path,
        outdir: Path,
        threads: int,
        memory_gb: int,
    ) -> Path:
        r1 = r1.resolve()
        r2 = r2.resolve()
        outdir = outdir.resolve()
        outdir.mkdir(parents=True, exist_ok=True)
        work = outdir / "spades_work"

        argv = [
            "metaspades.py",
            "-1",
            str(r1),
            "-2",
            str(r2),
            "-o",
            str(work),
            "-t",
            str(threads),
            "-m",
            str(memory_gb),
        ]
        self._run(
            argv,
            workdir=outdir,
            log_prefix="metaspades",
        )

        src = work / "contigs.fasta"
        if not src.is_file() or src.stat().st_size == 0:
            raise ToolError(f"metaSPAdes produced an empty contigs file: {src}")
        dest = outdir / "contigs.fasta"
        shutil.copy2(src, dest)
        scaffolds = work / "scaffolds.fasta"
        if scaffolds.is_file() and scaffolds.stat().st_size > 0:
            shutil.copy2(scaffolds, outdir / "scaffolds.fasta")
        LOGGER.info("Assembly contigs: %s", dest)
        return dest

    # ── Step 2: Kraken2 ───────────────────────────────────────────────────

    def kraken2(
        self,
        *,
        contigs: Path,
        db: Path,
        outdir: Path,
        threads: int,
        confidence: float,
        prefix: str,
    ) -> tuple[Path, Path, Path]:
        contigs = contigs.resolve()
        db = db.resolve()
        outdir = outdir.resolve()
        outdir.mkdir(parents=True, exist_ok=True)
        assignments = outdir / f"{prefix}.kraken2.out"
        report = outdir / f"{prefix}.kraken2.report"
        unclassified = outdir / f"{prefix}.unclassified.fasta"

        argv = [
            "kraken2",
            "--db",
            str(db),
            "--threads",
            str(threads),
            "--confidence",
            str(confidence),
            "--output",
            str(assignments),
            "--report",
            str(report),
            "--unclassified-out",
            str(unclassified),
            "--use-names",
            str(contigs),
        ]
        self._run(
            argv,
            workdir=outdir,
            log_prefix="kraken2",
        )
        if not assignments.is_file():
            raise ToolError(f"Kraken2 did not write assignments: {assignments}")
        if not report.is_file():
            raise ToolError(f"Kraken2 did not write report: {report}")
        if not unclassified.is_file():
            unclassified.write_text("", encoding="utf-8")
            LOGGER.warning("Kraken2 wrote no --unclassified-out file; treating as empty")
        LOGGER.info("Kraken2 unclassified contigs: %s", unclassified)
        return assignments, report, unclassified

    # ── Step 4: SeqScreen ─────────────────────────────────────────────────

    def seqscreen(
        self,
        *,
        fasta: Path,
        db: Path,
        workdir: Path,
        threads: int,
        mode: str = "fast",
        extra_args: str = "",
    ) -> Path:
        fasta = fasta.resolve()
        db = db.resolve()
        workdir = workdir.resolve()
        workdir.mkdir(parents=True, exist_ok=True)
        working = workdir / "seqscreen_out"
        working.mkdir(parents=True, exist_ok=True)

        argv = [
            "seqscreen",
            "--fasta",
            str(fasta),
            "--databases",
            str(db),
            "--working",
            str(working),
            "--threads",
            str(threads),
            "--mode",
            mode,
        ]
        if extra_args:
            argv.extend(shlex.split(extra_args))

        self._run(
            argv,
            workdir=workdir,
            log_prefix="seqscreen",
        )
        _normalise_seqscreen_report(working)
        LOGGER.info("SeqScreen working dir: %s", working)
        return working


def _normalise_seqscreen_report(working: Path) -> None:
    canonical = working / "report_generation" / "seqscreen_report.tsv"
    if canonical.is_file():
        return
    canonical.parent.mkdir(parents=True, exist_ok=True)
    found = sorted(working.rglob("*seqscreen_report*.tsv"))
    if found:
        shutil.copy2(found[0], canonical)
        return
    LOGGER.warning("SeqScreen report TSV not found; creating an empty stub.")
    canonical.write_text("query\n", encoding="utf-8")
