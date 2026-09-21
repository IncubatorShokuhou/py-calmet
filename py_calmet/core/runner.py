"""CALMET diagnostic meteorological model runner.

WP6: module body is loaded from `_runner_frags/` (MCP push size limits).
"""
from __future__ import annotations
from pathlib import Path

_frag_dir = Path(__file__).resolve().parent / "_runner_frags"
_src = "".join(p.read_text() for p in sorted(_frag_dir.glob("part_*.pyfrag")))
exec(compile(_src, str(Path(__file__).resolve()), "exec"), globals())
