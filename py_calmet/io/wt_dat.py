"""WT.DAT reader — overwater / sea-surface temperature (project soft-spot).

Project contract (README / PROGRESS soft spot): ``WTDAT`` supplies water
temperature used with SEA.DAT / ITWPROG precedence. Official CALMET User's
Guide §8.10 also documents WT.DAT as *terrain weighting factors*; that
layout is not implemented here (OutOfScope for this WP). This module
implements the documented soft-spot path: legacy / alternate SST input.

Precedence (wired in runner / overwater.resolve_water_temp)
-----------------------------------------------------------
1. ``ITWPROG ≠ 0`` + 3D.DAT → prognostic T over water LU cells
2. SEA.DAT stations present → nearest-station SST
3. ``WTDAT`` file present → this reader
4. else air temperature

Accepted formats
----------------
1. Domain-mean hourly SST::

       YYYY JJJ HH TSEA

2. Keyword::

       YYYY JJJ HH TSEA=285.0

3. Gridded (CLOUD-like)::

       WTSEA YYYYJJJHH v11 v12 ... (ny*nx values, j-major)

Header ``*`` / ``WT`` lines skipped. Temperatures in K; ``9999`` = missing.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List

import numpy as np

MISS = 9999.0


@dataclass
class WtRecord:
    year: int
    jday: int
    hour: int
    tsea: float = MISS  # domain-mean K
    grid: np.ndarray | None = None  # optional (ny, nx)

    @property
    def code(self) -> int:
        return self.year * 100000 + self.jday * 100 + self.hour

    def has_tsea(self) -> bool:
        if self.grid is not None:
            return True
        return self.tsea < MISS - 1.0 and self.tsea == self.tsea


@dataclass
class WtData:
    records: List[WtRecord] = field(default_factory=list)
    nx: int | None = None
    ny: int | None = None


def _code_to_yjh(code: int) -> tuple[int, int, int]:
    year = code // 100000
    rem = code % 100000
    return year, rem // 100, rem % 100


def read_wt(path: str | Path, nx: int | None = None, ny: int | None = None) -> WtData:
    """Read WT.DAT water-temperature file."""
    text = Path(path).read_text(errors="replace")
    lines = text.splitlines()
    records: List[WtRecord] = []
    tokens = text.replace("\n", " ").split()
    i = 0
    ncell = (nx * ny) if (nx and ny) else 0
    while i < len(tokens):
        lab = tokens[i]
        if lab.upper().startswith("WTSEA") or lab.upper() == "WTGRID":
            i += 1
            if i >= len(tokens):
                break
            try:
                code = int(float(tokens[i]))
                i += 1
            except ValueError:
                continue
            if ncell <= 0:
                break
            vals: list[float] = []
            while len(vals) < ncell and i < len(tokens):
                try:
                    vals.append(float(tokens[i]))
                    i += 1
                except ValueError:
                    break
            if len(vals) < ncell:
                break
            arr = np.array(vals[:ncell], dtype=np.float64).reshape((ny, nx))
            y, jd, h = _code_to_yjh(code)
            records.append(WtRecord(year=y, jday=jd, hour=h, grid=arr, tsea=float(np.nanmean(arr))))
        else:
            i += 1
    if records:
        return WtData(records=records, nx=nx, ny=ny)

    for line in lines:
        s = line.strip()
        if not s or s.startswith("*") or s.upper().startswith("WT"):
            continue
        if s.upper().startswith("WTSEA"):
            continue
        p = s.replace(",", " ").split()
        if len(p) < 4:
            continue
        try:
            year, jday, hour = int(float(p[0])), int(float(p[1])), int(float(p[2]))
        except ValueError:
            continue
        if hour == 24:
            hour = 0
            jday += 1
        tsea = MISS
        if "=" in p[3]:
            for tok in p[3:]:
                if tok.upper().startswith("TSEA="):
                    try:
                        tsea = float(tok.split("=", 1)[1])
                    except ValueError:
                        pass
        else:
            try:
                tsea = float(p[3])
            except ValueError:
                continue
        records.append(WtRecord(year=year, jday=jday, hour=hour, tsea=tsea))
    return WtData(records=records, nx=nx, ny=ny)


def wt_for_hour(
    data: WtData,
    year: int,
    jday: int,
    hour: int,
) -> WtRecord | None:
    if not data.records:
        return None
    code = year * 100000 + jday * 100 + hour
    best = None
    for r in data.records:
        if r.code == code:
            return r
        if r.code <= code and (best is None or r.code > best.code):
            best = r
    return best if best is not None else data.records[0]


def wt_sst_grid(
    rec: WtRecord,
    nx: int,
    ny: int,
    landuse: np.ndarray,
    iwat1: int,
    iwat2: int,
    tempk_fallback: np.ndarray,
) -> np.ndarray:
    """Apply WT.DAT SST on water LU cells; land keeps ``tempk_fallback``."""
    t = np.asarray(tempk_fallback, dtype=np.float64).copy()
    water = (landuse >= iwat1) & (landuse <= iwat2)
    if rec.grid is not None and rec.grid.shape == (ny, nx):
        t = np.where(water, rec.grid, t)
    elif rec.has_tsea():
        t = np.where(water, float(rec.tsea), t)
    return t
