"""Diagnostic wind construction (interp, OA, slope flow, mass consistency).

WP6: module body is loaded from `_winds_frags/` (MCP push size limits).
"""
from __future__ import annotations
from pathlib import Path

_frag_dir = Path(__file__).resolve().parent / "_winds_frags"
_src = "".join(p.read_text() for p in sorted(_frag_dir.glob("part_*.pyfrag")))
exec(compile(_src, str(Path(__file__).resolve()), "exec"), globals())
