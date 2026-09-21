"""UP.DAT reader (dataset 2.1, comma-delimited levels)."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import List
import numpy as np

# CALMET UP.DAT missing ≈ 999 (same family as VERTAV / obs_profile wind gate).
MISS = 999.0


def is_up_missing(v: float | int) -> bool:
    """True for non-finite or Fortran-style ≥998 missing codes (T/ws/wd)."""
    try:
        x = float(v)
    except (TypeError, ValueError):
        return True
    return (not np.isfinite(x)) or abs(x) >= 998.0


def up_tempk(lev: "UpLevel") -> float:
    """UP level temperature → Kelvin; missing → NaN (MIXDT/Holzworth filter).

    Raw ``temp_c == 999`` must not become 1272 K and silently poison
    Holzworth Zi or lapse rates that forget a miss gate.
    """
    if is_up_missing(lev.temp_c):
        return float("nan")
    return float(lev.temp_c) + 273.15


@dataclass
class UpLevel:
    pres: float
    height: float  # m MSL
    temp_c: float
    wd: float
    ws: float


@dataclass
class UpSounding:
    station_id: int
    year: int
    month: int
    day: int
    hour: int
    levels: List[UpLevel] = field(default_factory=list)


@dataclass
class UpData:
    soundings: List[UpSounding] = field(default_factory=list)


def read_up(path: str | Path) -> UpData:
    lines = Path(path).read_text().splitlines()
    i = 0
    i += 1
    ncom = int(lines[i].split()[0]); i += 1
    i += ncom
    i += 1  # NONE
    i += 1  # tz
    i += 1  # period header
    i += 1  # F F F F slash flags
    soundings = []
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        # sounding header: id year month day hour ...
        if ',' not in line:
            parts = line.split()
            if len(parts) < 6:
                i += 1
                continue
            sid = int(parts[0])
            year, month, day, hour = int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])
            nlev = int(parts[-1])
            i += 1
            # collect comma data until nlev levels
            vals: List[float] = []
            while len(vals) < nlev * 5 and i < len(lines):
                chunk = lines[i].replace(',', ' ').split()
                vals.extend(float(x) for x in chunk)
                i += 1
            levels = []
            for k in range(nlev):
                base = k * 5
                levels.append(UpLevel(
                    pres=vals[base], height=vals[base + 1],
                    temp_c=vals[base + 2], wd=vals[base + 3], ws=vals[base + 4],
                ))
            soundings.append(UpSounding(
                station_id=sid, year=year, month=month, day=day, hour=hour, levels=levels
            ))
        else:
            i += 1
    return UpData(soundings=soundings)
