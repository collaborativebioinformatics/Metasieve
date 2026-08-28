"""CLI wrappers: ``subprocess.run`` around assembly, taxonomy, and SeqScreen.

Tools are executed directly on the host (no Docker / Singularity).
Assembly uses GGCAT unitigs. ESMFold runs in-process via Hugging Face
(:mod:`metasieve.folding`).
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

    # ── Step 1: GGCAT unitigs ─────────────────────────────────────────────

    def ggcat(
        self,
        *,
        r1: Path,
        r2: Path,
        outdir: Path,
        threads: int,
        memory_gb: int,
        kmer: int = 31,
        min_multiplicity: int = 2,
        min_unitig_len: int = 200,
        force: bool = False,
    ) -> Path:
        """Build maximal unitigs from paired-end reads with GGCAT."""
        from metasieve.parsers import filter_fasta_min_length

        r1 = r1.resolve()
        r2 = r2.resolve()
        outdir = outdir.resolve()
        outdir.mkdir(parents=True, exist_ok=True)
        dest = outdir / "unitigs.fasta"
        if dest.is_file() and dest.stat().st_size > 0 and not force:
            LOGGER.info("Reusing existing unitigs: %s", dest)
            return dest

        work = outdir / "ggcat_work"
        tmp = work / "tmp"
        if work.exists() and force:
            shutil.rmtree(work)
        work.mkdir(parents=True, exist_ok=True)
        tmp.mkdir(parents=True, exist_ok=True)

        raw = work / "unitigs.raw.fasta"
        for stale in (raw, Path(str(raw) + ".lz4")):
            if stale.is_file():
                stale.unlink()
        argv = [
            "ggcat",
            "build",
            "-k",
            str(kmer),
            "-j",
            str(threads),
            "-m",
            str(memory_gb),
            "-s",
            str(min_multiplicity),
            "-t",
            str(tmp),
            "-o",
            str(raw),
            str(r1),
            str(r2),
        ]
        if not raw.exists():
            self._run(
                argv,
                workdir=outdir,
                log_prefix="ggcat",
            )

        produced = _ggcat_output_fasta(raw)
        if produced is None:
            raise ToolError(f"GGCAT produced no unitig FASTA at {raw}")
        n_kept = filter_fasta_min_length(produced, dest, min_length=min_unitig_len)
        if n_kept == 0 or not dest.is_file() or dest.stat().st_size == 0:
            raise ToolError(f"GGCAT produced no unitigs >= {min_unitig_len} bp: {dest}")
        LOGGER.info("Assembly unitigs: %s (%d sequences)", dest, n_kept)
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
        if not assignments.exists():
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
        LOGGER.info("Kraken2 unclassified sequences: %s", unclassified)
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
        canonical = working / "report_generation" / "seqscreen_report.tsv"

        if not canonical.exists():
            self._run(
                argv,
                workdir=workdir,
                log_prefix="seqscreen",
            )
        _normalise_seqscreen_report(working)
        LOGGER.info("SeqScreen working dir: %s", working)
        return working


def _ggcat_output_fasta(raw: Path) -> Path | None:
    """Return the GGCAT FASTA path, including a possible ``.lz4`` suffix."""
    candidates = [
        raw,
        Path(str(raw) + ".lz4"),
        raw.with_suffix(raw.suffix + ".lz4"),
        raw.with_suffix(".fasta"),
        raw.with_suffix(".fa"),
    ]
    for path in candidates:
        if path.is_file() and path.stat().st_size > 0:
            if path.suffix == ".lz4":
                return _decompress_lz4(path, raw if raw.suffix != ".lz4" else path.with_suffix(""))
            return path
    parent = raw.parent
    matches = sorted(parent.glob("*.fasta")) + sorted(parent.glob("*.fa")) + sorted(parent.glob("*.fasta.lz4"))
    for path in matches:
        if path.is_file() and path.stat().st_size > 0:
            if str(path).endswith(".lz4"):
                uncompressed = path.with_name(path.name[: -len(".lz4")])
                return _decompress_lz4(path, uncompressed)
            return path
    return None


def _decompress_lz4(src: Path, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    run_command(
        ["lz4", "-d", "-f", str(src), str(dest)],
        cwd=src.parent,
        log_dir=src.parent,
        log_prefix="lz4",
    )
    if not dest.is_file() or dest.stat().st_size == 0:
        raise ToolError(f"Failed to decompress GGCAT output {src}")
    return dest


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
