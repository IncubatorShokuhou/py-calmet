"""INP option helpers: legacy time, grid QA, file casefold, FCORIOL, IRAD, etc."""
from __future__ import annotations

from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import numpy as np


def resolve_data_file(
    inputs_dir: Path,
    case_dir: Path,
    name: str | None,
    fallback: str,
    *,
    lcfiles: bool = True,
) -> Path | None:
    """Honor INP filename; optionally case-fold match when LCFILES=T."""
    candidates: list[Path] = []
    if name:
        n = str(name).strip().strip("'\"")
        if n:
            candidates.extend([case_dir / n, inputs_dir / n, Path(n)])
    candidates.extend([case_dir / fallback, inputs_dir / fallback])
    for c in candidates:
        if c.is_file():
            return c
    if lcfiles:
        # Case-insensitive search in case_dir then inputs_dir
        for folder in (case_dir, inputs_dir):
            if not folder.is_dir():
                continue
            targets = set()
            if name:
                targets.add(str(name).strip().strip("'\"").lower())
            targets.add(fallback.lower())
            for p in folder.iterdir():
                if p.is_file() and p.name.lower() in targets:
                    return p
    return None


def run_window_with_legacy(inp) -> tuple[datetime, datetime, int, int]:
    """(start, end, nhrs, nsecdt) honoring IBTZ/IRLG legacy when needed."""
    ibyr = inp.get_int("IBYR", 2020)
    ibmo = inp.get_int("IBMO", 6)
    ibdy = inp.get_int("IBDY", 15)
    ibhr = inp.get_int("IBHR", 0)
    ibsec = inp.get_int("IBSEC", 0)
    ieyr = inp.get_int("IEYR", ibyr)
    iemo = inp.get_int("IEMO", ibmo)
    iedy = inp.get_int("IEDY", ibdy)
    iehr = inp.get_int("IEHR", ibhr + 3)
    iesec = inp.get_int("IESEC", 0)
    nsecdt = inp.get_int("NSECDT", 3600)
    start = datetime(ibyr, ibmo, ibdy, ibhr) + timedelta(seconds=ibsec)
    end = datetime(ieyr, iemo, iedy, iehr) + timedelta(seconds=iesec)

    # IRLG: legacy run length in hours when end == start or IRLG > 0 and end unset
    irlg = inp.get_int("IRLG", 0)
    if irlg > 0 and end <= start:
        end = start + timedelta(hours=irlg)

    span_sec = max(0, int((end - start).total_seconds()))
    nhrs = max(1, span_sec // max(nsecdt, 1))
    return start, end, nhrs, nsecdt


def ibtz_hours(inp) -> int:
    """Timezone hours west of GMT: prefer ABTZ, else legacy IBTZ."""
    abtz = inp.get("ABTZ")
    if abtz:
        s = str(abtz).upper().replace(" ", "")
        if "UTC" in s:
            rest = s.split("UTC", 1)[-1]
            sign = 1
            if rest.startswith("-"):
                sign = 1
                rest = rest[1:]
            elif rest.startswith("+"):
                sign = -1
                rest = rest[1:]
            try:
                return sign * int(rest[:2])
            except Exception:
                pass
    return int(inp.get_int("IBTZ", 0))


def effective_fcoriol(inp, lat0: float) -> float:
    """Use FCORIOL from INP when not the 999 sentinel; else 2Ωsinφ."""
    f = float(inp.get_float("FCORIOL", 999.0))
    if abs(f - 999.0) < 1e-6:
        from .met_utils import coriolis
        return float(coriolis(float(lat0)))
    return abs(f)


def qa_grid_vs_geo(inp, geo) -> list[str]:
    """Compare NX/NY/DGRIDKM/XORIGKM/YORIGKM to GEO.DAT; return notes."""
    notes = []
    nx = inp.get_int("NX", 0)
    ny = inp.get_int("NY", 0)
    if nx > 0 and nx != geo.nx:
        notes.append(f"NX INP={nx} vs GEO={geo.nx} (using GEO)")
    if ny > 0 and ny != geo.ny:
        notes.append(f"NY INP={ny} vs GEO={geo.ny} (using GEO)")
    dg = inp.get_float("DGRIDKM", 0.0)
    if dg > 0 and abs(dg - geo.dgridkm) > 1e-6:
        notes.append(f"DGRIDKM INP={dg} vs GEO={geo.dgridkm} (using GEO)")
    xo = inp.get_float("XORIGKM", 0.0)
    yo = inp.get_float("YORIGKM", 0.0)
    # Only flag if INP explicitly set nonzero and differs
    if abs(xo) + abs(yo) > 0 and (
        abs(xo - geo.xorigkm) > 1e-3 or abs(yo - geo.yorigkm) > 1e-3
    ):
        notes.append(
            f"X/YORIGKM INP=({xo},{yo}) vs GEO=({geo.xorigkm},{geo.yorigkm}) (using GEO)"
        )
    return notes


def qa_isteppgs(inp, nsecdt: int) -> list[str]:
    notes = []
    isteppgs = inp.get_int("ISTEPPGS", nsecdt)
    isteppg = inp.get_int("ISTEPPG", 1)
    if isteppgs > 0 and nsecdt > 0 and isteppgs % nsecdt != 0:
        notes.append(
            f"ISTEPPGS={isteppgs} not a multiple of NSECDT={nsecdt} (QA warning)"
        )
    if isteppg > 0 and isteppgs == nsecdt:
        # legacy hours consistent when isteppg * 3600 ≈ isteppgs
        pass
    return notes


def apply_nflagp(rates: np.ndarray, nflagp: int, cutp: float = 0.01) -> np.ndarray:
    """Precip QC: NFLAGP=0 pass-through; 1 zero missing; 2 also floor at CUTP.

    Missing includes negative rates and Fortran ≥9000 sentinels (PRECIP.DAT
    9999). Without the ≥9000 gate a raw 9999 mm/h would survive into Barnes OA.
    """
    r = np.asarray(rates, dtype=np.float64).copy()
    n = int(nflagp)
    if n >= 1:
        r = np.where((r < 0.0) | (r >= 9000.0) | ~np.isfinite(r), 0.0, r)
    if n >= 2:
        r = np.where(r < cutp, 0.0, r)
    return r


def average_temperature(
    temp2d: np.ndarray,
    *,
    iavet: int = 1,
    tradkm: float = 500.0,
    dgridkm: float = 1.0,
    numts: int = 5,
) -> np.ndarray:
    """Light spatial smoother for 2-D temperature when IAVET≠0.

    Uses a box of half-width min(NUMTS, tradkm/dgrid) cells. NUMTS=1 or
    IAVET=0 → identity (golden-safe).
    """
    if int(iavet) == 0 or int(numts) <= 1:
        return temp2d
    t = np.asarray(temp2d, dtype=np.float64)
    half = int(numts)
    if tradkm > 0 and dgridkm > 0:
        half = min(half, max(1, int(round(tradkm / dgridkm))))
    if half <= 1:
        return t
    # Separable box via cumulative sum
    pad = np.pad(t, half, mode="edge")
    c = np.cumsum(np.cumsum(pad, axis=0), axis=1)
    ny, nx = t.shape
    out = np.empty_like(t)
    w = (2 * half + 1) ** 2
    for j in range(ny):
        j0, j1 = j, j + 2 * half
        for i in range(nx):
            i0, i1 = i, i + 2 * half
            s = c[j1, i1] - (c[j0 - 1, i1] if j0 > 0 else 0)
            s -= c[j1, i0 - 1] if i0 > 0 else 0
            s += c[j0 - 1, i0 - 1] if (j0 > 0 and i0 > 0) else 0
            out[j, i] = s / w
    return out


def reject_mm4_mm5(inp) -> None:
    """MM4DAT / legacy MM5 paths are out of scope (wrfout→3D.DAT only)."""
    mm4 = inp.get("MM4DAT")
    if not mm4:
        return
    name = str(mm4).strip().strip("'\"").lower()
    if not name or name in {"mm4.dat", "mm5.dat"}:
        return  # default placeholder — ignored
    if "mm4" in name or "mm5" in name:
        raise NotImplementedError(
            f"MM4DAT={mm4!r} is out of scope; use wrfout → 3D.DAT (M3DDAT)."
        )


def multi_file_list(cfg_val, inp_key_prefix: str, inp, n: int) -> list[str]:
    """Collect up to n filenames from config list or INP UPDAT/M3DDAT/IGFDAT."""
    names: list[str] = []
    if cfg_val:
        raw = cfg_val if isinstance(cfg_val, (list, tuple)) else [cfg_val]
        names.extend(str(x).strip().strip("'\"") for x in raw if str(x).strip())
    # Also scan raw INP for UPDAT= / M3DDAT= style (single) — already in cfg
    return names[: max(n, 1)] if n > 0 else names[:1]


def parse_us1_coords(inp) -> tuple[float, float] | None:
    """US1 freeform → (x_km, y_km) or None."""
    val = inp.get("US1")
    if not val:
        # scan US*
        for key, v in sorted(getattr(inp, "raw", {}).items()):
            if key.startswith("US") and key[2:].isdigit():
                val = v
                break
    if not val:
        return None
    parts = str(val).replace("'", " ").split()
    try:
        return float(parts[2]), float(parts[3])
    except Exception:
        return None


def write_test_stubs(case_dir: Path, inp, *, enabled: bool) -> list[str]:
    """Touch DIAG/PROG/TST* / DCSTGD stub files when LDB or names set + enabled."""
    written = []
    if not enabled:
        return written
    mapping = {
        "DIADAT": "diag.dat",
        "PRGDAT": "prog.dat",
        "TSTPRT": "test.dat",
        "TSTOUT": "test.out",
        "TSTKIN": "test.kin",
        "TSTFRD": "test.frd",
        "TSTSLP": "test.slp",
        "DCSTGD": "dcst.grd",
    }
    for key, fallback in mapping.items():
        name = inp.get(key) or fallback
        p = case_dir / str(name).strip().strip("'\"")
        if not p.exists():
            p.write_text(f"py-calmet stub for {key}\n")
            written.append(str(p))
    return written
