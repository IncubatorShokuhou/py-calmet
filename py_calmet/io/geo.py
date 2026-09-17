"""GEO.DAT reader (dataset 2.0)."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np


@dataclass
class GeoData:
    nx: int
    ny: int
    xorigkm: float
    yorigkm: float
    dgridkm: float
    pmap: str
    utmzn: str
    datum: str
    landuse: np.ndarray  # [ny, nx] int
    elev: np.ndarray  # [ny, nx] float meters
    htfac: float = 1.0


def read_geo(path: str | Path) -> GeoData:
    lines = Path(path).read_text().splitlines()
    # header: dataset, ncom, comment, pmap, utmzn, datum+date, nx ny xorig yorig dx dy, units
    i = 0
    i += 1  # dataset
    ncom = int(lines[i].split()[0]); i += 1
    i += ncom  # comments
    pmap = lines[i].strip(); i += 1
    utmzn = lines[i].strip(); i += 1
    datum = lines[i][:8].strip(); i += 1
    parts = lines[i].split()
    nx, ny = int(parts[0]), int(parts[1])
    xorigkm, yorigkm = float(parts[2]), float(parts[3])
    dgridkm = float(parts[4]); i += 1
    i += 1  # units
    # IOPT1 line
    i += 1
    # LU: ny rows of nx
    lu = []
    for _ in range(ny):
        lu.append([int(x) for x in lines[i].split()])
        i += 1
    landuse = np.asarray(lu, dtype=np.int32)
    # HTFAC
    htfac = float(lines[i].split()[0]); i += 1
    elev = []
    for _ in range(ny):
        elev.append([float(x) * htfac for x in lines[i].split()])
        i += 1
    return GeoData(
        nx=nx, ny=ny, xorigkm=xorigkm, yorigkm=yorigkm, dgridkm=dgridkm,
        pmap=pmap, utmzn=utmzn, datum=datum,
        landuse=landuse, elev=np.asarray(elev, dtype=np.float64), htfac=htfac,
    )
