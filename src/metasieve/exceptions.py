"""Pipeline exceptions."""

from __future__ import annotations

from typing import Sequence


class PipelineError(RuntimeError):
    """Unrecoverable pipeline error."""


class ToolError(PipelineError):
    """A wrapped CLI tool exited non-zero."""

    def __init__(
        self,
        message: str,
        *,
        cmd: Sequence[str] | None = None,
        returncode: int | None = None,
        stdout: str | None = None,
        stderr: str | None = None,
    ) -> None:
        super().__init__(message)
        self.cmd = list(cmd) if cmd is not None else []
        self.returncode = returncode
        self.stdout = stdout or ""
        self.stderr = stderr or ""
