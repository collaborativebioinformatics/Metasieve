"""Hugging Face ESMFold structure prediction.

Loads ``facebook/esmfold_v1`` (or a local snapshot) once via
``transformers.EsmForProteinFolding`` and reuses the weights for every ORF.
"""

from __future__ import annotations

import logging
from pathlib import Path

from metasieve.config import PipelineConfig
from metasieve.exceptions import PipelineError, ToolError
from metasieve.parsers import read_fasta, write_fold_sidecar

LOGGER = logging.getLogger(__name__)

DEFAULT_MODEL_ID = "facebook/esmfold_v1"


class ESMFoldPredictor:
    """In-process ESMFold backed by the Hugging Face Transformers checkpoint."""

    def __init__(
        self,
        *,
        model_id: str = DEFAULT_MODEL_ID,
        cache_dir: Path | None = None,
        device: str = "auto",
        num_recycles: int = 4,
        chunk_size: int = 128,
    ) -> None:
        self.model_id = model_id
        self.cache_dir = Path(cache_dir).expanduser() if cache_dir else None
        self.device_request = device
        self.num_recycles = num_recycles
        self.chunk_size = chunk_size
        self._model = None
        self._device: str | None = None

    @classmethod
    def from_config(cls, config: PipelineConfig) -> ESMFoldPredictor:
        return cls(
            model_id=config.esmfold_model,
            cache_dir=config.esmfold_cache,
            device=config.esmfold_device,
            num_recycles=config.esmfold_num_recycles,
            chunk_size=config.esmfold_chunk_size,
        )

    def load(self) -> None:
        if self._model is not None:
            return
        try:
            import torch
            from transformers import EsmForProteinFolding
        except ImportError as exc:
            raise PipelineError(
                "ESMFold requires torch and transformers. "
                "Install with: pip install 'torch' 'transformers>=4.36'"
            ) from exc

        self._device = _resolve_device(self.device_request, torch)
        kwargs: dict = {}
        if self.cache_dir is not None:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            kwargs["cache_dir"] = str(self.cache_dir)
        local = Path(self.model_id).expanduser()
        source = str(local.resolve()) if local.exists() else self.model_id
        LOGGER.info("Loading Hugging Face ESMFold from %s (device=%s)", source, self._device)
        model = EsmForProteinFolding.from_pretrained(source, **kwargs)
        if self.chunk_size:
            model.trunk.set_chunk_size(int(self.chunk_size))
        if self._device == "cuda":
            model = model.cuda()
            try:
                model.esm = model.esm.half()
            except Exception as exc:  # pragma: no cover - optional memory tweak
                LOGGER.warning("Could not cast ESM encoder to float16: %s", exc)
        else:
            model = model.to(self._device)
        model.eval()
        self._model = model
        LOGGER.info("ESMFold ready")

    def predict_fasta(self, fasta: Path, outdir: Path, sample_id: str) -> tuple[Path, Path]:
        fasta = Path(fasta).resolve()
        outdir = Path(outdir).resolve()
        outdir.mkdir(parents=True, exist_ok=True)
        records = list(read_fasta(fasta))
        if not records:
            raise ToolError(f"No sequences in {fasta}")
        if len(records) > 1:
            LOGGER.warning("%s contains %d records; folding the first", fasta, len(records))
        sequence = records[0].sequence.replace("*", "").replace(" ", "").upper()
        if not sequence:
            raise ToolError(f"Empty peptide sequence in {fasta}")

        pdb_text, plddt = self.predict_sequence(sequence)
        if not pdb_text.strip():
            raise ToolError(f"ESMFold produced an empty PDB for {fasta.stem}")

        pdb_dest = outdir / f"{fasta.stem}.esmfold.pdb"
        pdb_dest.write_text(pdb_text if pdb_text.endswith("\n") else pdb_text + "\n", encoding="utf-8")
        sidecar = write_fold_sidecar(
            fasta=fasta,
            structure=pdb_dest,
            method="esmfold",
            sample_id=sample_id,
            output=outdir / f"{fasta.stem}.esmfold.json",
            model_name=self.model_id,
            plddt=plddt,
        )
        LOGGER.info("ESMFold structure: %s (pLDDT=%s)", pdb_dest, plddt or "n/a")
        return pdb_dest, sidecar

    def predict_sequence(self, sequence: str) -> tuple[str, str]:
        self.load()
        assert self._model is not None
        try:
            import torch

            context = torch.no_grad()
        except ImportError:
            from contextlib import nullcontext

            context = nullcontext()

        try:
            with context:
                output = self._model.infer(sequence, num_recycles=self.num_recycles)
            pdb_text = self._model.output_to_pdb(output)[0]
        except Exception as exc:
            raise ToolError(f"Hugging Face ESMFold inference failed: {exc}") from exc
        return pdb_text, mean_plddt(output)


def mean_plddt(output: object) -> str:
    """Return mean pLDDT on a 0–100 scale, or an empty string if unavailable."""
    value: float | None = None
    mean_attr = getattr(output, "mean_plddt", None)
    plddt_attr = getattr(output, "plddt", None)
    try:
        if mean_attr is not None:
            tensor = mean_attr.detach().float().cpu().reshape(-1) if hasattr(mean_attr, "detach") else mean_attr
            value = float(tensor[0] if hasattr(tensor, "__getitem__") else tensor)
        elif plddt_attr is not None:
            tensor = plddt_attr.detach().float().cpu() if hasattr(plddt_attr, "detach") else plddt_attr
            value = float(tensor.mean() if hasattr(tensor, "mean") else sum(tensor) / len(tensor))
    except (TypeError, ValueError, IndexError):
        return ""
    if value is None:
        return ""
    if value <= 1.0:
        value *= 100.0
    return f"{value:.2f}"


def _resolve_device(request: str, torch_mod) -> str:
    name = (request or "auto").lower().strip()
    if name == "auto":
        return "cuda" if torch_mod.cuda.is_available() else "cpu"
    if name in {"cuda", "cpu"}:
        if name == "cuda" and not torch_mod.cuda.is_available():
            LOGGER.warning("CUDA requested but not available; using CPU")
            return "cpu"
        return name
    raise PipelineError(f"Unknown --esmfold_device '{request}'. Choose auto, cuda, or cpu.")
