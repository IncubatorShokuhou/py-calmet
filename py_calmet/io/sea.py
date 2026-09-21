"""SEA.DAT reader (overwater stations, dataset 2.0 / 2.1 / 2.11)."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List


@dataclass
class SeaRecord:
    x_km: float
    y_km: float
    z_anem: float
    year1: int
    jday1: int
    hour1: int
    year2: int
    jday2: int
    hour2: int
    dt_air_sea: float  # air − sea (K); missing 9999
    t_air: float
    rh: float
    zi: float
    tgradb: float
    tgrada: float
    ws: float
    wd: float
    twave: float = -999.0
    hwave: float = -999.0
    z_tair: float = 10.0
    z_sst: float = 0.6
    lon_east: float = 0.0

    @property
    def t_sea(self) -> float:
        """SST from air temp and air−sea ΔT (dtow = Tair − Tsea).

        Missing ``t_air`` (≥9000 / non-finite) → NaN so OA falls back rather
        than poisoning overwater PBL with a phantom ~10k K SST. Missing ΔT
        alone keeps ``t_air`` as SST (Fortran-ish soft fallback).
        """
        import math
        ta = float(self.t_air)
        if (not math.isfinite(ta)) or ta >= 9000.0:
            return float("nan")
        dt = float(self.dt_air_sea)
        if (not math.isfinite(dt)) or dt >= 9000.0:
            return ta
        return ta - dt

    @property
    def begin_code(self) -> int:
        return self.year1 * 100000 + self.jday1 * 100 + self.hour1

    @property
    def end_code(self) -> int:
        return self.year2 * 100000 + self.jday2 * 100 + self.hour2


@dataclass
class SeaStation:
    station_id: int
    name: str
    version: float
    pmap: str = "UTM"
    records: List[SeaRecord] = field(default_factory=list)


@dataclass
class SeaData:
    stations: List[SeaStation] = field(default_factory=list)

    @property
    def nowsta(self) -> int:
        return len(self.stations)


def _parse_version(ver: str) -> float:
    s = (ver or "2.0").strip()
    try:
        return float(s[:4])
    except Exception:
        try:
            return float(s.split()[0])
        except Exception:
            return 2.0


def read_sea(path: str | Path) -> SeaStation:
    """Read one SEA.DAT file (one overwater station)."""
    lines = Path(path).read_text(errors="replace").splitlines()
    i = 0
    # dataset dataver datamod
    hdr = lines[i]
    i += 1
    parts = hdr.split()
    dataset = parts[0] if parts else ""
    dataver = parts[1] if len(parts) > 1 else "2.0"
    if "SEA" not in dataset.upper():
        # tolerate missing label
        pass
    ver = _parse_version(dataver)
    ncom = int(lines[i].split()[0]); i += 1
    i += ncom
    pmap = lines[i].strip(); i += 1
    # projection params
    if pmap.upper().startswith("UTM"):
        i += 1  # zone hem
    elif pmap.upper().startswith("LCC"):
        i += 2  # 4 lats + feast/fnorth
    elif pmap.upper().startswith("PS"):
        i += 1
    elif pmap.upper().startswith("EM"):
        i += 1
    else:
        i += 2  # LAZA/TTM
    i += 1  # datum date
    i += 1  # xyunit
    if ver >= 2.11:
        i += 1  # timezone
        i += 1  # start/end
    # station id + name
    id_parts = lines[i].split(None, 1); i += 1
    sid = int(float(id_parts[0]))
    name = id_parts[1].strip() if len(id_parts) > 1 else f"SEA{sid}"
    records: List[SeaRecord] = []
    while i < len(lines):
        line = lines[i].strip()
        i += 1
        if not line or line.startswith("*"):
            continue
        p = line.replace(",", " ").split()
        try:
            if ver >= 2.11:
                # x y z ztair zsst y1 j1 h1 y2 j2 h2 dt tair rh zi tgb tga ws wd tw hw
                rec = SeaRecord(
                    x_km=float(p[0]), y_km=float(p[1]), z_anem=float(p[2]),
                    z_tair=float(p[3]), z_sst=float(p[4]),
                    year1=int(p[5]), jday1=int(p[6]), hour1=int(p[7]),
                    year2=int(p[8]), jday2=int(p[9]), hour2=int(p[10]),
                    dt_air_sea=float(p[11]), t_air=float(p[12]), rh=float(p[13]),
                    zi=float(p[14]), tgradb=float(p[15]), tgrada=float(p[16]),
                    ws=float(p[17]), wd=float(p[18]),
                    twave=float(p[19]) if len(p) > 19 else -999.0,
                    hwave=float(p[20]) if len(p) > 20 else -999.0,
                )
            elif ver >= 2.1:
                rec = SeaRecord(
                    x_km=float(p[0]), y_km=float(p[1]), z_anem=float(p[2]),
                    year1=int(p[3]), jday1=int(p[4]), hour1=int(p[5]),
                    year2=int(p[6]), jday2=int(p[7]), hour2=int(p[8]),
                    dt_air_sea=float(p[9]), t_air=float(p[10]), rh=float(p[11]),
                    zi=float(p[12]), tgradb=float(p[13]), tgrada=float(p[14]),
                    ws=float(p[15]), wd=float(p[16]),
                    twave=float(p[17]) if len(p) > 17 else -999.0,
                    hwave=float(p[18]) if len(p) > 18 else -999.0,
                )
            else:
                # pre-2.1: x y lon z ...
                lon_w = float(p[2])
                rec = SeaRecord(
                    x_km=float(p[0]), y_km=float(p[1]), z_anem=float(p[3]),
                    year1=int(p[4]), jday1=int(p[5]), hour1=int(p[6]),
                    year2=int(p[7]), jday2=int(p[8]), hour2=int(p[9]),
                    dt_air_sea=float(p[10]), t_air=float(p[11]), rh=float(p[12]),
                    zi=float(p[13]), tgradb=float(p[14]), tgrada=float(p[15]),
                    ws=float(p[16]), wd=float(p[17]),
                    lon_east=-lon_w,
                )
            if rec.twave > 9998:
                rec.twave = -999.0
            if rec.hwave > 9998:
                rec.hwave = -999.0
            records.append(rec)
        except (ValueError, IndexError):
            continue
    return SeaStation(station_id=sid, name=name, version=ver, pmap=pmap, records=records)


def read_sea_files(paths: List[str | Path]) -> SeaData:
    return SeaData(stations=[read_sea(p) for p in paths])


def pick_sea_record(station: SeaStation, year: int, jday: int, hour: int) -> SeaRecord | None:
    """Pick record covering YYYYJJJHH (hour-ending style, CALMET)."""
    code = year * 100000 + jday * 100 + hour
    best = None
    for r in station.records:
        if r.begin_code <= code + 1 <= r.end_code + 1 or r.begin_code <= code <= r.end_code:
            best = r
            if r.begin_code <= code < r.end_code or r.end_code == code:
                return r
    return best or (station.records[0] if station.records else None)
