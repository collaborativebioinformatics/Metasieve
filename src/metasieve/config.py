"""Pipeline configuration: defaults, YAML files, and CLI overlays."""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, fields
from pathlib import Path
from typing import Any, Mapping

import yaml

from metasieve.exceptions import PipelineError

LOGGER = logging.getLogger(__name__)

DEFAULT_ESMFOLD_MODEL = "facebook/esmfold_v1"


@dataclass
class PipelineConfig:
    """All runtime settings for a Metasieve run."""

    reads: str | None = None
    r1: Path | None = None
    r2: Path | None = None
    sample_id: str | None = None
    outdir: Path = Path("results")

    kraken_db: Path | None = None
    seqscreen_db: Path | None = None

    esmfold_model: str = DEFAULT_ESMFOLD_MODEL
    esmfold_cache: Path | None = None
    esmfold_device: str = "auto"
    esmfold_num_recycles: int = 4
    esmfold_chunk_size: int = 128
    skip_esmfold: bool = False

    seqscreen_mode: str = "fast"
    seqscreen_extra_args: str = ""
    kraken2_confidence: float = 0.01
    min_contig_len: int = 200
    min_orf_aa: int = 50
    max_orf_aa: int = 1024
    require_orf_start: bool = True
    genetic_code: int = 11

    threads: int = 16
    memory_gb: int = 64
    fold_retries: int = 2

    log_level: str = "INFO"
    force: bool = False

    def resolved_outdir(self) -> Path:
        return self.outdir.expanduser().resolve()

    def validate(self) -> None:
        missing: list[str] = []
        if not self.reads and not (self.r1 and self.r2):
            missing.append("--reads (or --r1 and --r2)")
        if self.kraken_db is None:
            missing.append("--kraken_db")
        if self.seqscreen_db is None:
            missing.append("--seqscreen_db")
        if missing:
            raise PipelineError("Missing required parameters: " + ", ".join(missing))

        device = (self.esmfold_device or "auto").lower().strip()
        if device not in {"auto", "cuda", "cpu"}:
            raise PipelineError(
                f"Invalid --esmfold_device '{self.esmfold_device}'. Choose auto, cuda, or cpu."
            )
        self.esmfold_device = device

        _require_dir(self.kraken_db, "--kraken_db")
        _require_dir(self.seqscreen_db, "--seqscreen_db")
        if self.esmfold_cache is not None:
            self.esmfold_cache = self.esmfold_cache.expanduser()
            self.esmfold_cache.mkdir(parents=True, exist_ok=True)

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        for key, value in list(payload.items()):
            if isinstance(value, Path):
                payload[key] = str(value)
            elif isinstance(value, list):
                payload[key] = [str(v) if isinstance(v, Path) else v for v in value]
        return payload


def _require_dir(path: Path | None, flag: str) -> None:
    if path is None:
        return
    if not path.is_dir():
        raise PipelineError(f"{flag} is not a directory: {path}")


def _coerce_path(value: Any) -> Path | None:
    if value is None or value == "":
        return None
    return Path(str(value)).expanduser()


def _coerce_bool(value: Any, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


PATH_FIELDS = {
    "r1",
    "r2",
    "outdir",
    "kraken_db",
    "seqscreen_db",
    "esmfold_cache",
}

BOOL_FIELDS = {
    "require_orf_start",
    "skip_esmfold",
    "force",
}


def mapping_to_config(data: Mapping[str, Any], base: PipelineConfig | None = None) -> PipelineConfig:
    """Build a config from a YAML/CLI mapping, overlaying *base* defaults."""
    cfg = PipelineConfig() if base is None else PipelineConfig(**asdict(base))
    known = {item.name for item in fields(PipelineConfig)}
    unknown = [key for key in data if key not in known]
    if unknown:
        LOGGER.warning("Ignoring unknown config keys: %s", ", ".join(sorted(unknown)))

    updates: dict[str, Any] = {}
    for key in known:
        if key not in data or data[key] is None:
            continue
        value = data[key]
        if key in PATH_FIELDS:
            updates[key] = _coerce_path(value)
        elif key in BOOL_FIELDS:
            updates[key] = _coerce_bool(value, getattr(cfg, key))
        elif key == "outdir":
            updates[key] = Path(str(value)).expanduser()
        else:
            updates[key] = value

    for key, value in updates.items():
        setattr(cfg, key, value)
    return cfg


def load_yaml_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise PipelineError(f"Config file must be a mapping: {path}")
    return data


def overlay_cli(cfg: PipelineConfig, cli_values: Mapping[str, Any]) -> PipelineConfig:
    """Apply CLI values that the user actually set (non-None)."""
    filtered = {key: value for key, value in cli_values.items() if value is not None}
    return mapping_to_config(filtered, base=cfg)
