"""SURF.DAT reader (dataset 2.1)."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import List
import numpy as np

# CALMET SURF.DAT missing sentinel (same family as winds 9999).
MISS = 9999.0
DEFAULT_TEMPK = 288.15  # K — ISA; avoid phantom ~10k K in PBL/flux
DEFAULT_RH = 70  # % — matches 3D.DAT / noobs fallback
DEFAULT_PRES = 1012.0  # mb — matches noobs air_density default
DEFAULT_SKY = 0  # tenths — clear; missing must not explode ELUSTR theta1


def is_surf_missing(v: float | int) -> bool:
    """True for non-finite or Fortran-style ≥9000 missing codes."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return True
    return (not np.isfinite(x)) or x >= 9000.0


def surf_tempk(rec: "SurfRecord", default: float = DEFAULT_TEMPK) -> float:
    """Surface temperature (K); 9999 / non-finite → ``default``."""
    t = float(rec.tempk)
    return default if is_surf_missing(t) else t


def surf_rh(rec: "SurfRecord", default: int = DEFAULT_RH) -> int:
    """Relative humidity (%); 9999 / non-finite → ``default``."""
    r = float(rec.rh)
    return int(default) if is_surf_missing(r) else int(r)


def surf_pres(rec: "SurfRecord", default: float = DEFAULT_PRES) -> float:
    """Station pressure (mb); 9999 / non-finite → ``default``."""
    p = float(rec.pres)
    return default if is_surf_missing(p) else p


def surf_sky(rec: "SurfRecord", default: int = DEFAULT_SKY) -> int:
    """Sky cover (tenths 0–10); 9999 / out-of-range → ``default``.

    Missing sky would otherwise poison nighttime ELUSTR (``jcc**2``) and
    inflate cloud fraction before the [0,1] clip.
    """
    s = float(rec.sky)
    if is_surf_missing(s) or s < 0.0 or s > 10.0:
        return int(default)
    return int(s)


@dataclass
class SurfRecord:
    year: int
    jday: int
    hour: int
    ws: float
    wd: float
    ceil: int
    sky: int
    tempk: float
    rh: int
    pres: float
    ppcode: int


@dataclass
class SurfData:
    station_ids: List[int]
    records: List[SurfRecord] = field(default_factory=list)


def read_surf(path: str | Path) -> SurfData:
    lines = Path(path).read_text().splitlines()
    i = 0
    i += 1  # header
    ncom = int(lines[i].split()[0]); i += 1
    i += ncom
    i += 1  # NONE
    i += 1  # timezone
    # begin/end + nsta
    parts = lines[i].split(); i += 1
    nsta = int(parts[-1])
    ids = [int(x) for x in lines[i].split()]; i += 1
    records = []
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        # time line
        if i + 1 >= len(lines):
            break
        tparts = lines[i].split(); i += 1
        year, jday, hour = int(tparts[0]), int(tparts[1]), int(tparts[2])
        # one line per station for this hour
        for _ in range(nsta):
            p = lines[i].split(); i += 1
            records.append(SurfRecord(
                year=year, jday=jday, hour=hour,
                ws=float(p[0]), wd=float(p[1]),
                ceil=int(float(p[2])), sky=int(float(p[3])),
                tempk=float(p[4]), rh=int(float(p[5])),
                pres=float(p[6]), ppcode=int(float(p[7])),
            ))
    return SurfData(station_ids=ids, records=records)
