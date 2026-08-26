"""Discover paired-end FASTQ samples from a glob or explicit R1/R2 paths."""

from __future__ import annotations

import glob
import logging
import re
from dataclasses import dataclass
from pathlib import Path

from metasieve.exceptions import PipelineError

LOGGER = logging.getLogger(__name__)

_MATE_PATTERNS = (
    re.compile(r"^(?P<sample>.+?)_R?(?P<mate>[12])(?P<rest>\.(?:fastq|fq)(?:\.gz)?)$", re.I),
    re.compile(r"^(?P<sample>.+?)\.(?P<mate>[12])(?P<rest>\.(?:fastq|fq)(?:\.gz)?)$", re.I),
)


@dataclass(frozen=True)
class SampleReads:
    sample_id: str
    r1: Path
    r2: Path


def expand_brace_glob(pattern: str) -> list[Path]:
    """Expand ``{1,2}`` style globs then run POSIX globbing."""
    expanded = _expand_braces(pattern)
    files: list[Path] = []
    for item in expanded:
        matches = [Path(p) for p in glob.glob(item)]
        if not matches and Path(item).exists():
            matches = [Path(item)]
        files.extend(matches)
    unique = sorted({path.resolve() for path in files if path.is_file()})
    return unique


def _expand_braces(pattern: str) -> list[str]:
    match = re.search(r"\{([^{}]+)\}", pattern)
    if not match:
        return [pattern]
    options = match.group(1).split(",")
    prefix, suffix = pattern[: match.start()], pattern[match.end() :]
    out: list[str] = []
    for option in options:
        out.extend(_expand_braces(prefix + option + suffix))
    return out


def _mate_key(path: Path) -> tuple[str, str] | None:
    name = path.name
    for regex in _MATE_PATTERNS:
        match = regex.match(name)
        if match:
            return match.group("sample"), match.group("mate")
    return None


def pair_read_files(files: list[Path]) -> list[SampleReads]:
    buckets: dict[str, dict[str, Path]] = {}
    unmatched: list[Path] = []
    for path in files:
        parsed = _mate_key(path)
        if parsed is None:
            unmatched.append(path)
            continue
        sample, mate = parsed
        buckets.setdefault(sample, {})[mate] = path

    if unmatched:
        LOGGER.warning(
            "FASTQ files not recognised as paired mates and ignored: %s",
            ", ".join(str(p) for p in unmatched),
        )

    samples: list[SampleReads] = []
    for sample, mates in sorted(buckets.items()):
        if "1" not in mates or "2" not in mates:
            raise PipelineError(
                f"Sample '{sample}' is missing a mate. Found: {sorted(mates)}"
            )
        samples.append(SampleReads(sample_id=sample, r1=mates["1"], r2=mates["2"]))
    if not samples:
        raise PipelineError("No paired FASTQ files could be grouped into samples.")
    return samples


def discover_samples(
    *,
    reads: str | None,
    r1: Path | None,
    r2: Path | None,
    sample_id: str | None,
) -> list[SampleReads]:
    if r1 is not None and r2 is not None:
        r1 = r1.expanduser().resolve()
        r2 = r2.expanduser().resolve()
        if not r1.is_file() or not r2.is_file():
            raise PipelineError(f"Read files not found: {r1} / {r2}")
        name = sample_id or _infer_sample_id(r1)
        return [SampleReads(sample_id=name, r1=r1, r2=r2)]

    if not reads:
        raise PipelineError("Provide --reads or both --r1 and --r2")

    files = expand_brace_glob(reads)
    if not files:
        raise PipelineError(
            f"No paired reads matched --reads '{reads}'. "
            "Expected a glob such as 'data/*_{1,2}.fastq.gz'."
        )
    LOGGER.info("Matched %d FASTQ files from --reads", len(files))
    return pair_read_files(files)


def _infer_sample_id(r1: Path) -> str:
    parsed = _mate_key(r1)
    if parsed:
        return parsed[0]
    return r1.stem.split(".")[0]
