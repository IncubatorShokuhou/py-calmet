"""MIXDT / MIXDT2 — potential-temperature lapse above Zi (CALMET).

Used by daytime Maul–Carson growth: gamma = dθ/dz in a DZZI-deep layer above
the previous-hour convective mixing height, floored at DPTMIN.
"""
from __future__ import annotations

import numpy as np


def _interp_temp_at(z: np.ndarray, t: np.ndarray, z_tgt: float, miss: float = 999.0) -> float | None:
    """Linear T at ``z_tgt`` from valid sounding levels; None if undersampled."""
    z = np.asarray(z, dtype=np.float64)
    t = np.asarray(t, dtype=np.float64)
    valid = np.isfinite(z) & np.isfinite(t) & (t < miss - 0.01)
    if int(valid.sum()) < 2:
        return None
    zv, tv = z[valid], t[valid]
    order = np.argsort(zv)
    zv, tv = zv[order], tv[order]
    if z_tgt <= zv[0]:
        return float(tv[0])
    if z_tgt >= zv[-1]:
        return float(tv[-1])
    return float(np.interp(z_tgt, zv, tv))


def mixdt_sounding(
    zl: np.ndarray,
    tz: np.ndarray,
    htold: float,
    *,
    dptmin: float = 0.001,
    dzzi: float = 200.0,
    miss: float = 999.0,
) -> tuple[float, float, float]:
    """CALMET MIXDT: (tht, thtp, dtheta) from a 1-D sounding.

    ``dtheta`` is potential-temperature lapse (K/m) in ``[htold, htold+dzzi]``,
    with dry-adiabatic conversion ``+0.0098`` and floor ``dptmin``.
    """
    htold = max(float(htold), 0.0)
    dzzi = max(float(dzzi), 1.0)
    htpdz = htold + dzzi
    tht = _interp_temp_at(zl, tz, htold, miss=miss)
    thtp = _interp_temp_at(zl, tz, htpdz, miss=miss)
    if tht is None or thtp is None:
        dtheta = max(float(dptmin), 1e-4)
        return float(htold), float(htpdz), dtheta
    dtheta = (thtp - tht) / dzzi + 0.0098
    dtheta = max(float(dtheta), float(dptmin))
    return float(tht), float(thtp), float(dtheta)


def mixdt2_column(
    zl: np.ndarray,
    tz: np.ndarray,
    htold: float,
    *,
    dptmin: float = 0.001,
    dzzi: float = 200.0,
    miss: float = 999.0,
) -> tuple[float, float, float, float]:
    """CALMET MIXDT2 (prognostic column): (tsf, tht, thtp, dtheta)."""
    zl = np.asarray(zl, dtype=np.float64)
    tz = np.asarray(tz, dtype=np.float64)
    valid = np.isfinite(tz) & (tz < miss - 0.01)
    if np.any(valid):
        tsf = float(tz[valid][0])
    elif tz.size:
        tsf = float(tz[0])
    else:
        tsf = 288.0
    tht, thtp, dtheta = mixdt_sounding(zl, tz, htold, dptmin=dptmin, dzzi=dzzi, miss=miss)
    return tsf, tht, thtp, dtheta


def gamma_field_from_sounding(
    zl: np.ndarray,
    tz: np.ndarray,
    ziconv: np.ndarray,
    *,
    dptmin: float = 0.001,
    dzzi: float = 200.0,
) -> np.ndarray:
    """Scalar sounding → 2-D gamma field matching ``ziconv`` shape."""
    out = np.empty(np.shape(ziconv), dtype=np.float64)
    for idx in np.ndindex(out.shape):
        _, _, g = mixdt_sounding(
            zl, tz, float(ziconv[idx]), dptmin=dptmin, dzzi=dzzi
        )
        out[idx] = g
    return out


def _as_nk_ny_nx(arr: np.ndarray, ny: int, nx: int) -> np.ndarray:
    a = np.asarray(arr, dtype=np.float64)
    if a.ndim != 3:
        raise ValueError(f"expected 3-D array, got {a.shape}")
    if a.shape == (a.shape[0], ny, nx):
        return a
    if a.shape == (ny, nx, a.shape[-1]):
        return np.moveaxis(a, -1, 0)
    # Heuristic: trailing dims match grid
    if a.shape[-2:] == (ny, nx):
        return a
    if a.shape[:2] == (ny, nx):
        return np.moveaxis(a, -1, 0)
    raise ValueError(f"cannot map shape {a.shape} onto grid ({ny},{nx})")


def gamma_field_from_3d(
    height_agl: np.ndarray,
    tempk: np.ndarray,
    ziconv: np.ndarray,
    *,
    dptmin: float = 0.001,
    dzzi: float = 200.0,
) -> np.ndarray:
    """Per-column MIXDT2 from 3D.DAT-like arrays → (ny, nx) gamma."""
    zi = np.asarray(ziconv, dtype=np.float64)
    ny, nx = zi.shape
    z = _as_nk_ny_nx(height_agl, ny, nx)
    t = _as_nk_ny_nx(tempk, ny, nx)
    out = np.empty((ny, nx), dtype=np.float64)
    for j in range(ny):
        for i in range(nx):
            _, _, _, g = mixdt2_column(
                z[:, j, i], t[:, j, i], float(zi[j, i]), dptmin=dptmin, dzzi=dzzi
            )
            out[j, i] = g
    return out


def resolve_gamma(
    *,
    itprog: int,
    ziconv: np.ndarray,
    dptmin: float = 0.001,
    dzzi: float = 200.0,
    sounding_z: np.ndarray | None = None,
    sounding_t: np.ndarray | None = None,
    height_agl_3d: np.ndarray | None = None,
    tempk_3d: np.ndarray | None = None,
) -> np.ndarray | float:
    """Pick MIXDT (obs) / MIXDT2 (prog) / constant DPTMIN from ITPROG.

    ITPROG: 0 = obs sounding, 1 = hybrid (prefer prog if present else obs),
    2 = prognostic columns only. Falls back to scalar ``dptmin``.
    """
    it = int(itprog)
    if it >= 1 and height_agl_3d is not None and tempk_3d is not None:
        return gamma_field_from_3d(
            height_agl_3d, tempk_3d, ziconv, dptmin=dptmin, dzzi=dzzi
        )
    if it in (0, 1) and sounding_z is not None and sounding_t is not None:
        return gamma_field_from_sounding(
            sounding_z, sounding_t, ziconv, dptmin=dptmin, dzzi=dzzi
        )
    return max(float(dptmin), 1e-4)
