"""METLST list-file writer (human-readable run summary)."""
from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any


def write_metlst(
    path: str | Path,
    *,
    title: str = "py-calmet list file",
    inp_summary: dict[str, Any] | None = None,
    grid: dict[str, Any] | None = None,
    mode: str = "",
    nhrs: int = 0,
    notes: list[str] | None = None,
) -> None:
    """Write a Fortran-METLST-style text summary (not bit-identical)."""
    lines = [
        f"{title}",
        f"Generated: {datetime.now().astimezone().strftime('%Y-%m-%d %H:%M %Z')}",
        f"Mode: {mode}",
        f"Hours: {nhrs}",
        "",
        "--- Grid ---",
    ]
    for k, v in (grid or {}).items():
        lines.append(f"  {k} = {v}")
    lines.append("")
    lines.append("--- Key INP ---")
    for k, v in (inp_summary or {}).items():
        lines.append(f"  {k} = {v}")
    if notes:
        lines.append("")
        lines.append("--- Notes ---")
        lines.extend(f"  {n}" for n in notes)
    lines.append("")
    Path(path).write_text("\n".join(lines))
