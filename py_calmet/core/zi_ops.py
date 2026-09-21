"""Mixing-height post-processing: spatial average (IAVEZI) and relaxation (IZICRLX)."""
from __future__ import annotations

import numpy as np


def average_zi_upwind(
    zi: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    *,
    iavezi: int = 1,
    mnmdav: int = 1,
    hafang: float = 30.0,
    ilevzi: int = 1,
) -> np.ndarray:
    """Spatial average of Zi along upwind direction (CALMET AVEMHQ-ish).

    When ``iavezi==0`` or ``mnmdav<=1``, returns ``zi`` unchanged (identity —
    golden-safe for MNMDAV=1).
    """
    if int(iavezi) == 0 or int(mnmdav) <= 1:
        return zi
    z = np.asarray(zi, dtype=np.float64)
    ny, nx = z.shape
    n = int(mnmdav)
    half = np.deg2rad(float(hafang))
    # Layer for wind direction (1-based ilevzi → 0-based); u,v are (nz,ny,nx) or (ny,nx)
    if u.ndim == 3:
        k = int(np.clip(int(ilevzi) - 1, 0, u.shape[0] - 1))
        uu, vv = u[k], v[k]
    else:
        uu, vv = u, v
    wd = np.arctan2(-uu, -vv)  # met from-direction, radians
    out = z.copy()
    for j in range(ny):
        for i in range(nx):
            ang0 = float(wd[j, i])
            acc = float(z[j, i])
            cnt = 1
            # Walk upwind up to n cells within HAFANG cone
            for step in range(1, n + 1):
                # upwind vector (from where wind comes)
                di = int(np.round(np.sin(ang0) * step))
                dj = int(np.round(np.cos(ang0) * step))
                ii, jj = i + di, j + dj
                if not (0 <= ii < nx and 0 <= jj < ny):
                    break
                # angle check vs upwind
                dang = abs(((np.arctan2(di, dj) - ang0 + np.pi) % (2 * np.pi)) - np.pi)
                if dang > half:
                    continue
                acc += float(z[jj, ii])
                cnt += 1
            out[j, i] = acc / cnt
    return out


def relax_zi(
    zi: np.ndarray,
    zi_prev: np.ndarray,
    *,
    izicrlx: int = 1,
    tzicrlx: float = 800.0,
    dt_sec: float = 3600.0,
    daytime: np.ndarray | bool = True,
) -> np.ndarray:
    """Exponential relaxation of convective Zi toward previous hour.

    ``zi_new = zi_prev + (zi - zi_prev) * (1 - exp(-dt/tau))`` when
    ``izicrlx!=0`` and daytime; night / off → pass-through.
    """
    if int(izicrlx) == 0 or float(tzicrlx) <= 0:
        return zi
    z = np.asarray(zi, dtype=np.float64)
    zp = np.asarray(zi_prev, dtype=np.float64)
    if zp.shape != z.shape:
        return z
    alpha = 1.0 - np.exp(-float(dt_sec) / float(tzicrlx))
    blended = zp + (z - zp) * alpha
    day = np.asarray(daytime, dtype=bool)
    if day.shape != z.shape:
        day = np.full(z.shape, bool(daytime))
    # Only relax where previous Zi was already grown (avoid first-hour pull-to-zero)
    active = day & (zp > 1.0)
    return np.where(active, blended, z)


def mixht_holzworth(
    tempk_sfc: np.ndarray,
    sounding_z: np.ndarray,
    sounding_t: np.ndarray,
    *,
    zimin: float = 50.0,
    zimax: float = 3000.0,
) -> np.ndarray:
    """Simple Holzworth dry-adiabatic intercept mixing height (IMIXH=±3 style).

    Walks the sounding from the surface potential-temperature upward until
    the environmental pot-temp exceeds surface pot-temp.

    Levels with non-finite or unphysical Kelvin temps (UP.DAT 999 → NaN via
    ``up_tempk``, or raw ≥900 K) are skipped so missing T cannot yank Zi.
    """
    ts = np.asarray(tempk_sfc, dtype=np.float64)
    z = np.asarray(sounding_z, dtype=np.float64)
    t = np.asarray(sounding_t, dtype=np.float64)
    # Keep meteorological Kelvin only (filters UP missing 999°C→1272 K etc.)
    ok = np.isfinite(z) & np.isfinite(t) & (t > 150.0) & (t < 400.0)
    z, t = z[ok], t[ok]
    if z.size < 2:
        return np.full(ts.shape, zimin, dtype=np.float64)
    order = np.argsort(z)
    z, t = z[order], t[order]
    # pot temp along sounding (approx, p~const lapse)
    th_s = t + 0.0098 * z
    th_sfc = float(np.nanmean(ts)) if ts.ndim else float(ts)
    zi_col = zimin
    for i in range(1, z.size):
        if th_s[i] >= th_sfc:
            # linear intercept
            dth = th_s[i] - th_s[i - 1]
            if abs(dth) < 1e-6:
                zi_col = float(z[i])
            else:
                frac = (th_sfc - th_s[i - 1]) / dth
                zi_col = float(z[i - 1] + frac * (z[i] - z[i - 1]))
            break
    else:
        zi_col = float(z[-1])
    zi_col = float(np.clip(zi_col, zimin, zimax))
    return np.full(ts.shape, zi_col, dtype=np.float64)
