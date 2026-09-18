"""Precipitation fields for CALMET.DAT RMM (mm/h).

NPSTA meanings (CALMET):
  -1 — use prognostic precip from 3D.DAT / MM5 (gridded)
   0 — no precip
  >0 — station precip objective analysis
"""
from __future__ import annotations

import numpy as np


def prognostic_precip_grid(
    rain_mm: np.ndarray,
    nx: int,
    ny: int,
    xorig_km: float,
    yorig_km: float,
    dgrid_km: float,
    threed,
) -> np.ndarray:
    """Map 3D.DAT surface rain (mm per step) onto the CALMET grid.

    ``rain_mm`` is ``(nj, ni)`` on the prognostic grid. Nearest-neighbor
    (same mapping as wind interp).
    """
    from .winds import _nearest_3d_index

    out = np.zeros((ny, nx), dtype=np.float64)
    rain = np.asarray(rain_mm, dtype=np.float64)
    for j in range(ny):
        for i in range(nx):
            xc = xorig_km + (i + 0.5) * dgrid_km
            yc = yorig_km + (j + 0.5) * dgrid_km
            ii, jj = _nearest_3d_index(xc, yc, threed, dgrid_km)
            out[j, i] = max(float(rain[jj, ii]), 0.0)
    return out


def barnes_precip(
    stn_x_m: np.ndarray,
    stn_y_m: np.ndarray,
    stn_rmm: np.ndarray,
    nx: int,
    ny: int,
    xorig_m: float,
    yorig_m: float,
    dgrid_m: float,
    sigmap_km: float = 100.0,
    cutp: float = 0.01,
) -> np.ndarray:
    """Barnes-like OA of station precip rates (mm/h) → gridded RMM.

    ``sigmap_km`` is CALMET SIGMAP (influence radius). Values below ``cutp``
    are zeroed (CALMET CUTP).
    """
    xs = np.atleast_1d(np.asarray(stn_x_m, dtype=np.float64))
    ys = np.atleast_1d(np.asarray(stn_y_m, dtype=np.float64))
    rr = np.atleast_1d(np.asarray(stn_rmm, dtype=np.float64))
    nstn = int(xs.size)
    out = np.zeros((ny, nx), dtype=np.float64)
    rk = max(float(sigmap_km) * 1000.0, 1.0)
    for j in range(ny):
        for i in range(nx):
            xc = xorig_m + (i + 0.5) * dgrid_m
            yc = yorig_m + (j + 0.5) * dgrid_m
            dist2 = (xc - xs) ** 2 + (yc - ys) ** 2
            w = np.exp(-dist2 / rk**2)
            den = float(w.sum())
            if den > 1e-12:
                out[j, i] = float(np.dot(w, rr)) / den
    out = np.where(out < cutp, 0.0, out)
    return out


def resolve_precip(
    *,
    npsta: int,
    nx: int,
    ny: int,
    rain_prog: np.ndarray | None = None,
    threed=None,
    xorig_km: float = 0.0,
    yorig_km: float = 0.0,
    dgrid_km: float = 1.0,
    stn_x_m: np.ndarray | None = None,
    stn_y_m: np.ndarray | None = None,
    stn_rmm: np.ndarray | None = None,
    sigmap_km: float = 100.0,
    cutp: float = 0.01,
) -> np.ndarray:
    """Build RMM (ny, nx) from NPSTA switch."""
    n = int(npsta)
    if n == 0:
        return np.zeros((ny, nx), dtype=np.float64)
    if n < 0:
        if rain_prog is None or threed is None:
            return np.zeros((ny, nx), dtype=np.float64)
        return prognostic_precip_grid(
            rain_prog, nx, ny, xorig_km, yorig_km, dgrid_km, threed
        )
    if stn_x_m is None or stn_y_m is None or stn_rmm is None:
        return np.zeros((ny, nx), dtype=np.float64)
    return barnes_precip(
        stn_x_m,
        stn_y_m,
        stn_rmm,
        nx,
        ny,
        xorig_km * 1000.0,
        yorig_km * 1000.0,
        dgrid_km * 1000.0,
        sigmap_km=sigmap_km,
        cutp=cutp,
    )
