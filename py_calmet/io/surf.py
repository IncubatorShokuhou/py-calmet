"""SURF.DAT reader (dataset 2.1)."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import List
import numpy as np


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
