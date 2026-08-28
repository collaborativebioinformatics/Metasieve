"""Allow ``python -m metasieve``."""

from __future__ import annotations

from metasieve.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
