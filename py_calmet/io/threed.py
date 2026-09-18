"""3D.DAT reader (dataset 2.x, IOUTMM5 formats 81–95 uncompressed).

Format codes (CALMET RDMM5):
  81  P Z T WD WS
  82  + RH Q
  83  + QC QR
  84  + QI QS
  85  + QG
  91  P Z T WD WS W
  92  + RH Q          (default / golden path)
  93–95  compressed moisture variants (read uncompressed layout when present)
"""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import List
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
    elev: np.ndarray  # [nj, ni]
    hours: List[str]
    wd10: np.ndarray
    ws10: np.ndarray
    t2: np.ndarray
    rain: np.ndarray
    pres: np.ndarray
    height_msl: np.ndarray
    tempk: np.ndarray
    wd: np.ndarray
    ws: np.ndarray
    w: np.ndarray
    rh: np.ndarray
    ioutmm5: int = 92
    qc: np.ndarray | None = None
    qr: np.ndarray | None = None
    qi: np.ndarray | None = None
    qs: np.ndarray | None = None
    qg: np.ndarray | None = None
    qv: np.ndarray | None = None  # vapor mixing ratio g/kg


def _parse_iout_flags(line: str) -> int:
    """Infer IOUTMM5 from the ioutw/ioutq/ioutc/iouti/ioutg flag line."""
    parts = line.split()
    # Common layouts: "IOUTW IOUTQ IOUTC IOUTI IOUTG" or a single integer
    try:
        ints = [int(float(x)) for x in parts[:5]]
    except Exception:
        return 92
    if len(ints) == 1 and ints[0] >= 81:
        return ints[0]
    while len(ints) < 5:
        ints.append(0)
    ioutw, ioutq, ioutc, iouti, ioutg = ints[:5]
    # CALMET: ioutmm5 = 81 + 10*ioutw + ioutq + ioutc + iouti + ioutg
    # (with cloud/ice/graupel bits as 0/1 flags in practice)
    return 81 + 10 * ioutw + ioutq + ioutc + iouti + ioutg


def _parse_upper_line(uline: str, ioutmm5: int) -> dict:
    """Parse one upper-air record according to IOUTMM5."""
    # Fixed-width core always starts: i4, i6, f6.1, i4, f5.1
    # Be tolerant of free-format whitespace as well.
    out = {
        "pmb": 0.0, "z": 0.0, "temp": 0.0, "wd": 0.0, "ws": 0.0,
        "w": 0.0, "rh": 70.0, "q": 0.0,
        "qc": 0.0, "qr": 0.0, "qi": 0.0, "qs": 0.0, "qg": 0.0,
    }
    # Prefer fixed-width when line is long enough
    try:
        if len(uline) >= 25 and uline[0:4].strip().lstrip("+-").isdigit():
            out["pmb"] = float(uline[0:4])
            out["z"] = float(uline[4:10])
            out["temp"] = float(uline[10:16])
            out["wd"] = float(uline[16:20])
            out["ws"] = float(uline[20:25])
            pos = 25
            code = int(ioutmm5)
            if code in (91, 92, 93, 94, 95) or code >= 91:
                if len(uline) >= pos + 6:
                    out["w"] = float(uline[pos:pos + 6]); pos += 6
            if code in (82, 83, 84, 85, 92, 93, 94, 95):
                if len(uline) >= pos + 3:
                    out["rh"] = float(uline[pos:pos + 3]); pos += 3
                if len(uline) >= pos + 5:
                    out["q"] = float(uline[pos:pos + 5]); pos += 5
            # Moisture extras — consume remaining floats (f5.2 or f6.3)
            rest = uline[pos:].strip()
            if rest and code in (83, 84, 85, 93, 94, 95):
                toks = rest.replace("  ", " ").split()
                names = ["qc", "qr", "qi", "qs", "qg"]
                for name, tok in zip(names, toks):
                    # Skip compression flag if negative sentinel
                    try:
                        v = float(tok)
                    except ValueError:
                        break
                    if name == "qc" and v < -0.0001 and code >= 93:
                        # compressed → zeros already set
                        break
                    out[name] = v
            return out
    except Exception:
        pass
    # Free-format fallback
    toks = uline.replace("  ", " ").split()
    try:
        out["pmb"] = float(toks[0]); out["z"] = float(toks[1])
        out["temp"] = float(toks[2]); out["wd"] = float(toks[3])
        out["ws"] = float(toks[4])
        idx = 5
        code = int(ioutmm5)
        if code >= 91 and idx < len(toks):
            out["w"] = float(toks[idx]); idx += 1
        if code in (82, 83, 84, 85, 92, 93, 94, 95) and idx < len(toks):
            out["rh"] = float(toks[idx]); idx += 1
            if idx < len(toks):
                out["q"] = float(toks[idx]); idx += 1
        for name in ("qc", "qr", "qi", "qs", "qg"):
            if idx < len(toks):
                out[name] = float(toks[idx]); idx += 1
    except Exception:
        pass
    return out


def read_3d(path: str | Path) -> ThreeDData:
    lines = Path(path).read_text().splitlines()
    i = 0
    i += 1  # dataset
    ncom = int(lines[i].split()[0]); i += 1
    i += ncom
    iout_line = lines[i]; i += 1
    ioutmm5 = _parse_iout_flags(iout_line)
    map_parts = lines[i].split(); i += 1
    x0 = float(map_parts[-6]); y0 = float(map_parts[-5]); dx = float(map_parts[-4])
    ni = int(map_parts[-3]); nj = int(map_parts[-2]); nk = int(map_parts[-1])
    i += 1  # surface var options
    hdr = lines[i]; i += 1
    rest = hdr[10:].split()
    nhrs = int(rest[0])
    ext = lines[i]
    i += 1
    i0 = int(ext[0:4])
    j0 = int(ext[4:8])
    sigma = np.array([float(lines[i + k].split()[0]) for k in range(nk)], dtype=np.float64)
    i += nk
    elev = np.zeros((nj, ni), dtype=np.float64)
    for _jj in range(nj):
        for _ii in range(ni):
            parts = lines[i].split()
            i += 1
            ii = int(np.clip(int(parts[0]) - i0, 0, ni - 1))
            jj = int(np.clip(int(parts[1]) - j0, 0, nj - 1))
            elev[jj, ii] = float(parts[4])

    wd10 = np.zeros((nhrs, nj, ni)); ws10 = np.zeros_like(wd10); t2 = np.zeros_like(wd10)
    rain = np.zeros_like(wd10)
    pres = np.zeros((nhrs, nj, ni, nk)); height_msl = np.zeros_like(pres)
    tempk = np.zeros_like(pres); wd = np.zeros_like(pres); ws = np.zeros_like(pres)
    w = np.zeros_like(pres); rh = np.zeros_like(pres)
    qv = np.zeros_like(pres)
    qc = np.zeros_like(pres); qr = np.zeros_like(pres)
    qi = np.zeros_like(pres); qs = np.zeros_like(pres); qg = np.zeros_like(pres)
    hours: List[str] = []

    for t in range(nhrs):
        for _row in range(nj):
            for _col in range(ni):
                sline = lines[i]
                i += 1
                stamp = sline[:10]
                if _row == 0 and _col == 0:
                    hours.append(stamp)
                ii_idx = int(sline[10:13])
                jj_idx = int(sline[13:16])
                ii = int(np.clip(ii_idx - i0, 0, ni - 1))
                jj = int(np.clip(jj_idx - j0, 0, nj - 1))
                rest = sline[16:]
                vals = rest.replace("  ", " ").split()
                if len(vals) > 1:
                    try:
                        rain[t, jj, ii] = float(vals[1])
                    except Exception:
                        rain[t, jj, ii] = 0.0
                if len(vals) > 5:
                    t2[t, jj, ii] = float(vals[5])
                if len(vals) > 8:
                    wd10[t, jj, ii] = float(vals[7])
                    ws10[t, jj, ii] = float(vals[8])
                for k in range(nk):
                    uline = lines[i]
                    i += 1
                    parsed = _parse_upper_line(uline, ioutmm5)
                    pres[t, jj, ii, k] = parsed["pmb"]
                    height_msl[t, jj, ii, k] = parsed["z"]
                    tempk[t, jj, ii, k] = parsed["temp"]
                    wd[t, jj, ii, k] = parsed["wd"]
                    ws[t, jj, ii, k] = parsed["ws"]
                    w[t, jj, ii, k] = parsed["w"]
                    rh[t, jj, ii, k] = parsed["rh"]
                    qv[t, jj, ii, k] = parsed["q"]
                    qc[t, jj, ii, k] = parsed["qc"]
                    qr[t, jj, ii, k] = parsed["qr"]
                    qi[t, jj, ii, k] = parsed["qi"]
                    qs[t, jj, ii, k] = parsed["qs"]
                    qg[t, jj, ii, k] = parsed["qg"]

    # If RH missing (formats 81/91), leave default 70
    has_cloud = int(ioutmm5) in (83, 84, 85, 93, 94, 95)
    return ThreeDData(
        ni=ni, nj=nj, nk=nk, x0_km=x0, y0_km=y0, dx_km=dx, sigma=sigma, elev=elev,
        hours=hours, wd10=wd10, ws10=ws10, t2=t2, rain=rain,
        pres=pres, height_msl=height_msl, tempk=tempk, wd=wd, ws=ws, w=w, rh=rh,
        ioutmm5=ioutmm5,
        qc=qc if has_cloud else None,
        qr=qr if has_cloud else None,
        qi=qi if int(ioutmm5) in (84, 85, 94, 95) else None,
        qs=qs if int(ioutmm5) in (84, 85, 94, 95) else None,
        qg=qg if int(ioutmm5) in (85, 95) else None,
        qv=qv,
    )
