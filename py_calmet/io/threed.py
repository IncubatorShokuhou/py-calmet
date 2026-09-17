"""3D.DAT reader (dataset 2.1, ioutmm5 format 92 uncompressed)."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple
import numpy as np


@dataclass
class ThreeDData:
    ni: int
    nj: int
    nk: int
    x0_km: float
    y0_km: float
    dx_km: float
    sigma: np.ndarray
    elev: np.ndarray  # [nj, ni] from metadata
    # time-major arrays
    hours: List[str]  # YYYYMMDDHH
    # surface
    wd10: np.ndarray  # [nt, nj, ni]
    ws10: np.ndarray
    t2: np.ndarray
    # upper
    pres: np.ndarray  # [nt, nj, ni, nk]
    height_msl: np.ndarray
    tempk: np.ndarray
    wd: np.ndarray
    ws: np.ndarray
    w: np.ndarray
    rh: np.ndarray


def read_3d(path: str | Path) -> ThreeDData:
    lines = Path(path).read_text().splitlines()
    i = 0
    i += 1  # dataset
    ncom = int(lines[i].split()[0]); i += 1
    i += ncom
    i += 1  # iout flags
    # map line: LCC lat0 lon0 ... x0 y0 dx ni nj nk  OR similar
    map_parts = lines[i].split(); i += 1
    # find trailing numbers: x0 y0 dx ni nj nk near end
    # format from make: LCC lat lon lat1 lat2 x0 y0 dx ni nj nk
    x0 = float(map_parts[-6]); y0 = float(map_parts[-5]); dx = float(map_parts[-4])
    ni = int(map_parts[-3]); nj = int(map_parts[-2]); nk = int(map_parts[-1])
    i += 1  # surface var options
    # grid data header
    hdr = lines[i]; i += 1
    # YYYYMMDDHH + nhrs + ni + nj + nk
    ymdh = hdr[:10]
    rest = hdr[10:].split()
    nhrs = int(rest[0])
    i += 1  # extraction subdomain
    sigma = np.array([float(lines[i + k].split()[0]) for k in range(nk)], dtype=np.float64)
    i += nk
    elev = np.zeros((nj, ni), dtype=np.float64)
    for jj in range(nj):
        for ii_ in range(ni):
            parts = lines[i].split(); i += 1
            # i j lat lon elev lu ...
            elev[jj, ii_] = float(parts[4])

    wd10 = np.zeros((nhrs, nj, ni)); ws10 = np.zeros_like(wd10); t2 = np.zeros_like(wd10)
    pres = np.zeros((nhrs, nj, ni, nk)); height_msl = np.zeros_like(pres)
    tempk = np.zeros_like(pres); wd = np.zeros_like(pres); ws = np.zeros_like(pres)
    w = np.zeros_like(pres); rh = np.zeros_like(pres)
    hours: List[str] = []

    for t in range(nhrs):
        for jj in range(nj):
            for ii_ in range(ni):
                sline = lines[i]; i += 1
                stamp = sline[:10]
                if jj == 0 and ii_ == 0:
                    hours.append(stamp)
                # stamp(10) i(3) j(3) spres(7.1) rain(5.2) sc(2) radsw(8.1) radlw(8.1) t2(8.1) q2(8.2) wd10(8.1) ws10(8.1) sst(8.1)
                body = sline[10:]
                # fixed-ish parse via split after removing stamp
                # use known widths from format
                ii_idx = int(sline[10:13]); jj_idx = int(sline[13:16])
                # remaining floats
                rest = sline[16:]
                # spres rain sc radsw radlw t2 q2 wd10 ws10 sst
                # sc is int in middle — split carefully
                vals = rest.replace('  ', ' ').split()
                # vals: spres rain sc radsw radlw t2 q2 wd10 ws10 sst
                t2[t, jj, ii_] = float(vals[5])
                wd10[t, jj, ii_] = float(vals[7])
                ws10[t, jj, ii_] = float(vals[8])
                for k in range(nk):
                    uline = lines[i]; i += 1
                    # pmb(4) z(6) temp(6.1) wd(4) ws(5.1) w(6.2) rh(3) vapmr(5.2)
                    pmb = int(uline[0:4]); z = int(uline[4:10])
                    temp = float(uline[10:16]); wdi = int(uline[16:20])
                    wsi = float(uline[20:25]); wi = float(uline[25:31])
                    rhi = int(uline[31:34])
                    pres[t, jj, ii_, k] = pmb
                    height_msl[t, jj, ii_, k] = z
                    tempk[t, jj, ii_, k] = temp
                    wd[t, jj, ii_, k] = wdi
                    ws[t, jj, ii_, k] = wsi
                    w[t, jj, ii_, k] = wi
                    rh[t, jj, ii_, k] = rhi

    return ThreeDData(
        ni=ni, nj=nj, nk=nk, x0_km=x0, y0_km=y0, dx_km=dx, sigma=sigma, elev=elev,
        hours=hours, wd10=wd10, ws10=ws10, t2=t2,
        pres=pres, height_msl=height_msl, tempk=tempk, wd=wd, ws=ws, w=w, rh=rh,
    )
