"""PRECIP.DAT reader (formatted IFORMP=2; free-form rates mm/h)."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List


@dataclass
class PrecipRecord:
    year: int
    jday: int
    hour: int
    rates: List[float]  # mm/h per station; 9999 = missing


@dataclass
class PrecipData:
    nsta: int
    records: List[PrecipRecord] = field(default_factory=list)
    station_ids: List[int] = field(default_factory=list)


def read_precip(path: str | Path, npsta: int | None = None) -> PrecipData:
    """Read formatted PRECIP.DAT (IFORMP=2 style).

    Formats accepted:
      * Header comments (``*`` / dataset lines) then ``YYYY JJJ HH r1 r2 ...``
      * Optional station-id line after timezone / begin-end block (SURF-like)
    """
    lines = Path(path).read_text(errors="replace").splitlines()
    i = 0
    # Skip dataset / comment header if present
    if i < len(lines) and "PRECIP" in lines[i].upper():
        i += 1
        if i < len(lines) and lines[i].strip() and lines[i].split()[0].isdigit():
            ncom = int(lines[i].split()[0]); i += 1
            i += ncom
        # optional NONE / timezone / begin-end nsta
        while i < len(lines) and (
            lines[i].strip().upper() in ("NONE",) or lines[i].strip().startswith("UTC")
        ):
            i += 1
        if i < len(lines):
            parts = lines[i].split()
            if len(parts) >= 7 and parts[0].isdigit():
                # begin/end + nsta
                nsta_hdr = int(parts[-1])
                i += 1
                if npsta is None:
                    npsta = nsta_hdr
                if i < len(lines) and not _looks_like_time(lines[i]):
                    ids = [int(float(x)) for x in lines[i].split()]
                    i += 1
                else:
                    ids = list(range(1, (npsta or nsta_hdr) + 1))
            else:
                ids = []
                nsta_hdr = npsta or 1
        else:
            ids = []
            nsta_hdr = npsta or 1
    else:
        # Skip leading comment lines
        while i < len(lines) and (not lines[i].strip() or lines[i].lstrip().startswith("*")):
            i += 1
        ids = []
        nsta_hdr = npsta or 0

    records: List[PrecipRecord] = []
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line or line.startswith("*"):
            continue
        p = line.replace(",", " ").split()
        if len(p) < 4:
            continue
        try:
            year, jday, hour = int(p[0]), int(p[1]), int(p[2])
            rates = [float(x) for x in p[3:]]
        except ValueError:
            continue
        if hour == 24:
            hour = 0
            jday += 1
        if npsta is not None and len(rates) < npsta:
            rates.extend([9999.0] * (npsta - len(rates)))
        records.append(PrecipRecord(year=year, jday=jday, hour=hour, rates=rates))
        if nsta_hdr == 0:
            nsta_hdr = len(rates)
    nsta = npsta or nsta_hdr or (len(records[0].rates) if records else 0)
    if not ids:
        ids = list(range(1, nsta + 1))
    return PrecipData(nsta=nsta, records=records, station_ids=ids)


def _looks_like_time(line: str) -> bool:
    p = line.split()
    if len(p) < 3:
        return False
    try:
        int(p[0]); int(p[1]); int(p[2])
        return True
    except ValueError:
        return False


def rates_for_hour(
    data: PrecipData,
    year: int,
    jday: int,
    hour: int,
    missing: float = 0.0,
) -> list[float]:
    """Return station rates (mm/h) for matching hour; missing → ``missing``."""
    for r in data.records:
        if r.year == year and r.jday == jday and r.hour == hour:
            return [missing if v > 9000.0 else float(v) for v in r.rates]
    # nearest earlier
    best = None
    code = year * 100000 + jday * 100 + hour
    for r in data.records:
        rc = r.year * 100000 + r.jday * 100 + r.hour
        if rc <= code and (best is None or rc > best[0]):
            best = (rc, r)
    if best is None:
        if not data.records:
            return [missing] * data.nsta
        r = data.records[0]
    else:
        r = best[1]
    return [missing if v > 9000.0 else float(v) for v in r.rates]
