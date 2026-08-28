"""Discover samples from paired FASTQ or restart FASTA inputs."""

from __future__ import annotations

import glob
import logging
import re
from dataclasses import dataclass
from pathlib import Path

from metasieve.exceptions import PipelineError

LOGGER = logging.getLogger(__name__)

START_ASSEMBLY = "assembly"
START_CLASSIFICATION = "classification"
START_UNCLASSIFIED = "unclassified"
START_FOLDING = "folding"
START_AUTO = "auto"
START_STAGES = (
    START_ASSEMBLY,
    START_CLASSIFICATION,
    START_UNCLASSIFIED,
    START_FOLDING,
)
_START_ALIASES = {
    "assembly": START_ASSEMBLY,
    "ggcat": START_ASSEMBLY,
    "1": START_ASSEMBLY,
    "classification": START_CLASSIFICATION,
    "kraken": START_CLASSIFICATION,
    "kraken2": START_CLASSIFICATION,
    "unitigs": START_CLASSIFICATION,
    "2": START_CLASSIFICATION,
    "unclassified": START_UNCLASSIFIED,
    "filtering": START_UNCLASSIFIED,
    "filter": START_UNCLASSIFIED,
    "seqscreen": START_UNCLASSIFIED,
    "3": START_UNCLASSIFIED,
    "4": START_UNCLASSIFIED,
    "folding": START_FOLDING,
    "esmfold": START_FOLDING,
    "orfs": START_FOLDING,
    "5": START_FOLDING,
}
_NT_SUFFIXES = {".fa", ".fasta", ".fna", ".fsa", ".fas"}
_MATE_PATTERNS = (
    re.compile(r"^(?P<sample>.+?)_R?(?P<mate>[12])(?P<rest>\.(?:fastq|fq)(?:\.gz)?)$", re.I),
    re.compile(r"^(?P<sample>.+?)\.(?P<mate>[12])(?P<rest>\.(?:fastq|fq)(?:\.gz)?)$", re.I),
)


@dataclass(frozen=True)
class SampleReads:
    sample_id: str
    r1: Path | None = None
    r2: Path | None = None
    unitigs: Path | None = None
    unclassified: Path | None = None
    orfs: Path | None = None


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


def normalize_start_from(value: str | None) -> str:
    text = (value or START_AUTO).strip().lower()
    if text in {START_AUTO, ""}:
        return START_AUTO
    if text in START_STAGES:
        return text
    mapped = _START_ALIASES.get(text)
    if mapped is None:
        raise PipelineError(
            f"Invalid --start_from '{value}'. Choose auto, assembly, "
            "classification, unclassified, or folding."
        )
    return mapped


def infer_start_from(
    *,
    start_from: str | None,
    reads: str | None,
    r1: Path | None,
    r2: Path | None,
    unitigs: str | None,
    unclassified: str | None,
    orfs: str | None,
) -> str:
    """Pick the first pipeline step from --start_from or provided inputs."""
    explicit = normalize_start_from(start_from)
    if explicit != START_AUTO:
        return explicit
    if orfs:
        return START_FOLDING
    if unclassified:
        return START_UNCLASSIFIED
    if unitigs:
        return START_CLASSIFICATION
    if reads or (r1 and r2):
        return START_ASSEMBLY
    return START_ASSEMBLY


def discover_pipeline_samples(
    *,
    start_from: str,
    reads: str | None,
    r1: Path | None,
    r2: Path | None,
    sample_id: str | None,
    unitigs: str | None,
    unclassified: str | None,
    orfs: str | None,
) -> list[SampleReads]:
    if start_from == START_FOLDING:
        return discover_fasta_samples(orfs, sample_id=sample_id, kind="orfs")
    if start_from == START_UNCLASSIFIED:
        return discover_fasta_samples(
            unclassified, sample_id=sample_id, kind="unclassified"
        )
    if start_from == START_CLASSIFICATION:
        return discover_fasta_samples(unitigs, sample_id=sample_id, kind="unitigs")
    return discover_samples(reads=reads, r1=r1, r2=r2, sample_id=sample_id)


def discover_fasta_samples(
    pattern: str | Path | None,
    *,
    sample_id: str | None,
    kind: str,
) -> list[SampleReads]:
    if not pattern:
        raise PipelineError(f"Provide --{kind} to start from that step")
    raw = str(pattern).strip()
    path = Path(raw).expanduser()
    entries: list[Path] = []
    if path.is_dir():
        if kind == "orfs":
            entries = [path.resolve()]
        else:
            entries = sorted(
                p.resolve()
                for p in path.iterdir()
                if p.is_file() and p.suffix.lower() in _NT_SUFFIXES
            )
            if not entries:
                raise PipelineError(f"--{kind} directory has no FASTA files: {path}")
    elif path.is_file():
        entries = [path.resolve()]
    else:
        entries = expand_brace_glob(raw)
        if not entries and Path(raw).exists():
            entries = [Path(raw).expanduser().resolve()]

    if not entries:
        raise PipelineError(f"No files matched --{kind} '{pattern}'")

    if sample_id and len(entries) > 1:
        raise PipelineError(
            f"--sample_id cannot be used with multiple --{kind} inputs ({len(entries)} found)"
        )

    samples: list[SampleReads] = []
    seen: set[str] = set()
    for entry in entries:
        sid = sample_id or _infer_fasta_sample_id(entry, kind)
        if sid in seen:
            raise PipelineError(
                f"Duplicate sample_id '{sid}' from --{kind}. "
                "Rename files or pass one input with --sample_id."
            )
        seen.add(sid)
        kwargs: dict[str, Path] = {kind: entry}
        samples.append(SampleReads(sample_id=sid, **kwargs))
        LOGGER.info("Restart input %s: %s -> sample %s", kind, entry, sid)
    return samples


def _infer_fasta_sample_id(path: Path, kind: str) -> str:
    if path.is_dir():
        return path.name
    name = path.name
    suffixes = [
        f".{kind}.fasta",
        f".{kind}.fa",
        f".{kind}.fna",
        f".{kind}.faa",
        ".unclassified.fasta",
        ".unclassified.fa",
        ".unitigs.fasta",
        ".unitigs.fa",
        ".orfs.faa",
        ".orfs.fa",
        ".contigs.fasta",
        ".fasta",
        ".fa",
        ".fna",
        ".faa",
        ".fsa",
        ".fas",
    ]
    lower = name
    for suffix in suffixes:
        if lower.endswith(suffix) or lower.lower().endswith(suffix.lower()):
            stem = name[: -len(suffix)]
            if stem:
                return stem
    return path.stem.split(".")[0]


def _infer_sample_id(r1: Path) -> str:
    parsed = _mate_key(r1)
    if parsed:
        return parsed[0]
    return r1.stem.split(".")[0]
