"""METLST list-file writer (human-readable run summary + IPR*/LPRINT sections)."""
from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np


def write_metlst(
    path: str | Path,
    *,
    title: str = "py-calmet list file",
    inp_summary: dict[str, Any] | None = None,
    grid: dict[str, Any] | None = None,
    mode: str = "",
    nhrs: int = 0,
    notes: list[str] | None = None,
    lprint: bool = False,
    ipr_flags: dict[str, int] | None = None,
    print_fields: dict[str, bool] | None = None,
    field_samples: dict[str, Any] | None = None,
    iuvout: list[int] | None = None,
    iwout: list[int] | None = None,
    itout: list[int] | None = None,
    hour_dumps: list[dict[str, Any]] | None = None,
) -> None:
    """Write a Fortran-METLST-style text summary (not bit-identical).

    When ``lprint`` is True (or any IPR* / field print flag is set), appends
    verbosity sections mirroring CALMET list-file surface/layer dumps.
    """
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

    ipr = ipr_flags or {}
    pfields = print_fields or {}
    verbose = bool(lprint) or any(int(v) != 0 for v in ipr.values()) or any(pfields.values())

    if verbose:
        lines.append("")
        lines.append("--- Print / IPR controls ---")
        lines.append(f"  LPRINT = {bool(lprint)}")
        for k in sorted(ipr):
            lines.append(f"  {k} = {ipr[k]}")
        for k, v in pfields.items():
            lines.append(f"  {k} = {bool(v)}")
        if iuvout:
            lines.append(f"  IUVOUT layers = {iuvout}")
        if iwout:
            lines.append(f"  IWOUT layers = {iwout}")
        if itout:
            lines.append(f"  ITOUT layers = {itout}")

        if field_samples:
            lines.append("")
            lines.append("--- Field samples (domain mean / min / max) ---")
            for name, arr in field_samples.items():
                a = np.asarray(arr, dtype=np.float64)
                if a.size == 0:
                    continue
                lines.append(
                    f"  {name}: mean={float(np.nanmean(a)):.4g}  "
                    f"min={float(np.nanmin(a)):.4g}  max={float(np.nanmax(a)):.4g}"
                )

        if hour_dumps:
            lines.append("")
            lines.append("--- Hour dumps ---")
            for h in hour_dumps:
                lines.append(f"  hour {h.get('hour', '?')}: {h}")

    # IPR0..IPR8 echo (even when zero) when any provided
    if ipr:
        lines.append("")
        lines.append("--- IPR0–IPR8 ---")
        for i in range(9):
            key = f"IPR{i}"
            lines.append(f"  {key} = {int(ipr.get(key, 0))}")

    if notes:
        lines.append("")
        lines.append("--- Notes ---")
        lines.extend(f"  {n}" for n in notes)
    lines.append("")
    Path(path).write_text("\n".join(lines))
