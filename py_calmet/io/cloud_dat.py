"""CLOUD.DAT reader / writer (gridded cloud fraction, IFORMC=2 formatted)."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List
import numpy as np


@dataclass
class CloudRecord:
    year: int
    jday: int
    hour: int
    ccgrid: np.ndarray  # (ny, nx) fraction 0–1
    label: str = "CLOUDFRA"


@dataclass
class CloudData:
    nx: int
    ny: int
    records: List[CloudRecord] = field(default_factory=list)


def _code_to_yjh(code: int) -> tuple[int, int, int]:
    year = code // 100000
    rem = code % 100000
    jday = rem // 100
    hour = rem % 100
    return year, jday, hour


def read_cloud(path: str | Path, nx: int, ny: int) -> CloudData:
    """Read formatted CLOUD.DAT (label + YYYYJJJHH + nx*ny values)."""
    text = Path(path).read_text(errors="replace")
    tokens = text.replace("\n", " ").split()
    records: List[CloudRecord] = []
    i = 0
    ncell = nx * ny
    while i < len(tokens):
        lab = tokens[i]
        if lab.upper().startswith("CLOUD"):
            i += 1
            if i >= len(tokens):
                break
            try:
                code = int(float(tokens[i])); i += 1
            except ValueError:
                continue
            vals = []
            while len(vals) < ncell and i < len(tokens):
                try:
                    vals.append(float(tokens[i])); i += 1
                except ValueError:
                    break
            if len(vals) < ncell:
                break
            arr = np.array(vals[:ncell], dtype=np.float64).reshape((ny, nx), order="C")
            # Fortran writes ((cc(i,j),i=1,nx),j=1,ny) → row-major j,i if we reshape (ny,nx) with
            # vals in i-fastest: reshape (ny,nx) after filling i then j → use order F from flat
            arr = np.array(vals[:ncell], dtype=np.float64).reshape((nx, ny), order="F").T
            y, jd, h = _code_to_yjh(code)
            records.append(CloudRecord(year=y, jday=jd, hour=h, ccgrid=arr, label=lab))
        else:
            i += 1
    return CloudData(nx=nx, ny=ny, records=records)


def write_cloud(path: str | Path, data: CloudData, *, iformc: int = 2) -> None:
    """Write formatted CLOUD.DAT (IFORMC=2)."""
    lines = []
    for r in data.records:
        code = r.year * 100000 + r.jday * 100 + r.hour
        # ((cc(i,j), i=1,nx), j=1,ny)
        flat = np.asarray(r.ccgrid, dtype=np.float64).T.ravel(order="F")
        # Actually: Fortran column i varies fastest for fixed j → values[j*nx + i]
        flat = np.asarray(r.ccgrid, dtype=np.float64).ravel(order="C")  # j-major rows
        # CALMET: ((ccgrid(i,j),i=1,nx),j=1,ny) → for j in rows, for i in cols
        flat = []
        cc = np.asarray(r.ccgrid, dtype=np.float64)
        for j in range(data.ny):
            for i in range(data.nx):
                flat.append(float(cc[j, i]))
        body = " ".join(f"{v:.4f}" for v in flat)
        lines.append(f"{r.label} {code} {body}")
    Path(path).write_text("\n".join(lines) + "\n")


def cloud_for_hour(
    data: CloudData,
    year: int,
    jday: int,
    hour: int,
) -> np.ndarray | None:
    code = year * 100000 + jday * 100 + hour
    best = None
    for r in data.records:
        rc = r.year * 100000 + r.jday * 100 + r.hour
        if rc == code:
            return r.ccgrid
        if rc <= code and (best is None or rc > best[0]):
            best = (rc, r)
    return best[1].ccgrid if best else (data.records[0].ccgrid if data.records else None)
