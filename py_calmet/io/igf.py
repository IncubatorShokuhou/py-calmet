"""IGF-CALMET.DAT reader stub (prior CALMET.DAT as initial-guess field).

MM4/MM5.DAT are out of scope — use wrfout → 3D.DAT. IGF is a prior
CALMET.DAT (or header-compatible) file used when IGFMET≠0.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class IgfHeader:
    path: Path
    title: str = ""
    nx: int = 0
    ny: int = 0
    nz: int = 0
    nhrs: int = 0
    notes: list[str] = field(default_factory=list)
    ok: bool = False


@dataclass
class IgfData:
    header: IgfHeader
    U: np.ndarray | None = None  # (nhrs, nz, ny, nx)
    V: np.ndarray | None = None
    meta: dict[str, Any] = field(default_factory=dict)


def read_igf_header(path: str | Path) -> IgfHeader:
    """Parse what we can from an IGF / prior CALMET.DAT path.

    Tries the binary CALMET.DAT reader first; on failure returns a
    header-only stub so NIGF/IGFDAT switches are still functional.
    """
    p = Path(path)
    hdr = IgfHeader(path=p, title=p.name)
    if not p.is_file():
        hdr.notes.append(f"missing file: {p}")
        return hdr
    try:
        from .calmet_dat import read_calmet_dat

        ds = read_calmet_dat(p)
        hdr.nx = int(ds.nx)
        hdr.ny = int(ds.ny)
        hdr.nz = int(ds.nz)
        hdr.nhrs = int(ds.nt)
        if "U" in ds.fields_3d:
            u = ds.fields_3d["U"]
            hdr.nhrs, hdr.nz, hdr.ny, hdr.nx = (int(x) for x in u.shape[:4])
        hdr.ok = True
        hdr.notes.append("parsed via read_calmet_dat")
        return hdr
    except Exception as exc:  # noqa: BLE001
        hdr.notes.append(f"header stub (calmet_dat failed: {exc})")
        try:
            raw = p.read_bytes()[:256]
            hdr.title = raw.split(b"\n")[0][:80].decode("latin1", errors="replace")
        except Exception:
            pass
        return hdr


def read_igf(path: str | Path) -> IgfData:
    """Load IGF winds when the file is a readable CALMET.DAT; else header-only."""
    hdr = read_igf_header(path)
    data = IgfData(header=hdr)
    if not hdr.ok:
        return data
    try:
        from .calmet_dat import read_calmet_dat

        ds = read_calmet_dat(path)
        data.U = ds.fields_3d.get("U")
        data.V = ds.fields_3d.get("V")
        data.meta = {
            "nx": ds.nx, "ny": ds.ny, "nz": ds.nz, "nt": ds.nt,
            "path": str(path),
        }
    except Exception as exc:  # noqa: BLE001
        hdr.notes.append(f"fields unavailable: {exc}")
        hdr.ok = False
    return data
