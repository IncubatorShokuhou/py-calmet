"""Cloud-cover diagnostics (CALMET CLOUD3 / CLOUD4-style).

MCLOUD / ICLOUD method mapping (py-calmet):
  0 / default — use SURF sky tenths when available (caller)
  3 — Teixeira RH@~850 mb (CLOUD3)
  4 — MM5toGrads layered RH → total cloud fraction (CLOUD4-lite)
"""
from __future__ import annotations

import numpy as np


def cloud3_from_rh(rh_pct: np.ndarray) -> np.ndarray:
    """Teixeira (2001) cloud fraction from RH (%) at ~850 mb (CALMET CLOUD3).

    A = 0.02 * (-1 + sqrt(1 + 100*(1-rh))) / (1-rh)  for RH < 99%;
    A = 1 when RH >= 99%.
    """
    rh = np.asarray(rh_pct, dtype=np.float64)
    cc = np.empty_like(rh, dtype=np.float64)
    sat = rh >= 99.0
    cc[sat] = 1.0
    rhd = np.clip(rh[~sat] / 100.0, 0.0, 0.989)
    denom = np.maximum(1.0 - rhd, 1e-6)
    cc[~sat] = 0.02 * (-1.0 + np.sqrt(1.0 + 100.0 * (1.0 - rhd))) / denom
    return np.clip(cc, 0.0, 1.0)


def cloud4_from_rh_profile(
    rh_3d: np.ndarray,
    pres_mb: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """CLOUD4-lite: layered RH → total cloud fraction + ceiling height (m AGL proxy).

    ``rh_3d`` / ``pres_mb`` shaped ``(nk, ny, nx)`` or ``(ny, nx, nk)``.
    Returns ``(ccp, ceil_m)`` on ``(ny, nx)``.
    """
    rh = np.asarray(rh_3d, dtype=np.float64)
    pr = np.asarray(pres_mb, dtype=np.float64)
    if rh.ndim != 3:
        raise ValueError(f"rh_3d must be 3-D, got {rh.shape}")
    # Normalize to (nk, ny, nx)
    if rh.shape[0] < rh.shape[-1] and rh.shape[-1] > 4:
        # likely (ny, nx, nk)
        rh = np.moveaxis(rh, -1, 0)
        pr = np.moveaxis(pr, -1, 0)
    nk, ny, nx = rh.shape
    clfrlo = np.zeros((ny, nx), dtype=np.float64)
    clfrmi = np.zeros((ny, nx), dtype=np.float64)
    clfrhi = np.zeros((ny, nx), dtype=np.float64)
    # Pressure decreases with k in CALMET (surface → aloft)
    for k in range(nk):
        p = pr[k]
        r = rh[k]
        lo = (p < 970.0) & (p >= 800.0)
        mi = (p < 800.0) & (p >= 450.0)
        hi = p < 450.0
        # Also treat near-surface high-P as low cloud layer candidate
        lo = lo | ((p >= 970.0) & (k == 0))
        clfrlo = np.where(lo, np.maximum(clfrlo, r), clfrlo)
        clfrmi = np.where(mi, np.maximum(clfrmi, r), clfrmi)
        clfrhi = np.where(hi, np.maximum(clfrhi, r), clfrhi)

    clo = np.clip(4.0 * clfrlo / 100.0 - 3.0, 0.0, 1.0)
    cmi = np.clip(4.0 * clfrmi / 100.0 - 3.0, 0.0, 1.0)
    chi = np.clip(2.5 * clfrhi / 100.0 - 1.5, 0.0, 1.0)
    ccp = np.maximum(np.maximum(clo, cmi), chi)
    # Ceiling proxy: higher cloud → lower ceiling index
    ceil = np.where(clo >= 0.5, 1000.0, np.where(cmi >= 0.5, 3000.0, np.where(chi >= 0.5, 6000.0, 0.0)))
    ceil = np.where(ccp < 1e-9, 0.0, ceil)
    return ccp, ceil


def rh850_from_3d(
    rh: np.ndarray,
    pres_mb: np.ndarray,
    target_mb: float = 850.0,
) -> np.ndarray:
    """Nearest-level RH (%) to ``target_mb`` on each column.

    ``rh`` / ``pres_mb``: ``(nk, ny, nx)`` or time-slice ``(nj, ni, nk)`` from 3D.DAT.
    """
    rh_a = np.asarray(rh, dtype=np.float64)
    pr_a = np.asarray(pres_mb, dtype=np.float64)
    if rh_a.ndim == 3 and rh_a.shape[0] > rh_a.shape[-1]:
        # (nj, ni, nk) → (nk, nj, ni)
        rh_a = np.moveaxis(rh_a, -1, 0)
        pr_a = np.moveaxis(pr_a, -1, 0)
    nk = rh_a.shape[0]
    # Pick level closest to target using domain-mean pressure
    pmean = np.array([float(np.nanmean(pr_a[k])) for k in range(nk)])
    k850 = int(np.argmin(np.abs(pmean - target_mb)))
    return rh_a[k850]


def resolve_cloud_fraction(
    *,
    mcloud: int = 0,
    icloud: int = 0,
    sky_tenths: float | np.ndarray | None = None,
    rh_pct_2d: np.ndarray | None = None,
    rh_3d: np.ndarray | None = None,
    pres_3d: np.ndarray | None = None,
    shape: tuple[int, int] | None = None,
) -> np.ndarray:
    """Select cloud fraction [0,1] from MCLOUD/ICLOUD and available fields.

    Prefer MCLOUD when set (INP 2.2); else fall back to ICLOUD; else SURF sky.
    """
    method = int(mcloud) if int(mcloud) not in (0, 999) else int(icloud)
    if method in (3,) and rh_pct_2d is not None:
        return cloud3_from_rh(rh_pct_2d)
    if method in (3,) and rh_3d is not None and pres_3d is not None:
        return cloud3_from_rh(rh850_from_3d(rh_3d, pres_3d))
    if method in (4,) and rh_3d is not None and pres_3d is not None:
        rh = np.asarray(rh_3d, dtype=np.float64)
        pr = np.asarray(pres_3d, dtype=np.float64)
        if rh.ndim == 3 and rh.shape[-1] < rh.shape[0]:
            # (nj,ni,nk)
            ccp, _ = cloud4_from_rh_profile(np.moveaxis(rh, -1, 0), np.moveaxis(pr, -1, 0))
        else:
            ccp, _ = cloud4_from_rh_profile(rh, pr)
        return ccp
    # Default / ICLOUD 0,1: surface sky tenths
    if sky_tenths is not None:
        cc = np.asarray(sky_tenths, dtype=np.float64) * 0.1
        if cc.ndim == 0 and shape is not None:
            return np.full(shape, float(cc), dtype=np.float64)
        return np.clip(cc, 0.0, 1.0)
    if shape is not None:
        return np.zeros(shape, dtype=np.float64)
    return np.array(0.0, dtype=np.float64)
