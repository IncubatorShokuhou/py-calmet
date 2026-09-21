"""Diagnostic wind construction (interp, OA, slope flow, mass consistency)."""
from __future__ import annotations
import numpy as np
from .met_utils import wind_uv, ZO_EXTRAP, layer_mids, G


def _nearest_3d_index(
    x_km: float,
    y_km: float,
    threed,
    dgrid_km: float,
) -> tuple[int, int]:
    """Nearest 3D.DAT mass point for a CALMET cell-center (km)."""
    dx = float(getattr(threed, "dx_km", 0.0) or dgrid_km)
    dx = max(dx, 1e-6)
    x0 = float(getattr(threed, "x0_km", x_km - dx))
    y0 = float(getattr(threed, "y0_km", y_km - dx))
    ii = int(np.round((x_km - x0) / dx - 0.5))
    jj = int(np.round((y_km - y0) / dx - 0.5))
    ii = int(np.clip(ii, 0, threed.ni - 1))
    jj = int(np.clip(jj, 0, threed.nj - 1))
    return ii, jj


def interp_3d_to_calmet(
    threed,
    zface: np.ndarray,
    nx: int,
    ny: int,
    xorig_km: float,
    yorig_km: float,
    dgrid_km: float,
    hour_index: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Map 3D.DAT winds onto CALMET layers (RDMM5-style log below first level).

    Horizontal mapping uses geographic cell centers vs 3D.DAT origin/spacing
    (the common 1-cell MM5 halo is a special case of this, not a hardcoded +1).
    """
    zmid = layer_mids(zface)
    nz = len(zmid)
    U = np.zeros((nz, ny, nx), dtype=np.float64)
    V = np.zeros_like(U)
    t = hour_index
    for j in range(ny):
        for i in range(nx):
            xc = xorig_km + (i + 0.5) * dgrid_km
            yc = yorig_km + (j + 0.5) * dgrid_km
            ii, jj = _nearest_3d_index(xc, yc, threed, dgrid_km)
            elev = float(threed.elev[jj, ii])
            zs = threed.height_msl[t, jj, ii, :] - elev
            zs = np.maximum.accumulate(np.maximum(zs, 1.0))
            us, vs = wind_uv(threed.wd[t, jj, ii, :], threed.ws[t, jj, ii, :])
            for L, zm in enumerate(zmid):
                if zm < zs[0]:
                    ratio = (np.log(zm) - np.log(ZO_EXTRAP)) / (
                        np.log(zs[0]) - np.log(ZO_EXTRAP)
                    )
                    U[L, j, i] = ratio * us[0]
                    V[L, j, i] = ratio * vs[0]
                else:
                    U[L, j, i] = np.interp(zm, zs, us)
                    V[L, j, i] = np.interp(zm, zs, vs)
    return U, V


def obs_surface_uv(ws: float, wd: float, nx: int, ny: int) -> tuple[np.ndarray, np.ndarray]:
    u, v = wind_uv(wd, ws)
    return np.full((ny, nx), u, dtype=np.float64), np.full((ny, nx), v, dtype=np.float64)


def obs_profile_similt(
    u_sfc: float,
    v_sfc: float,
    z_anem: float,
    z0: float,
    el: float,
    zi: float,
    zface: np.ndarray,
    sounding_levels,
    stn_elev: float,
    zimin: float,
    nx: int,
    ny: int,
    p_exp: float = 0.17,
    iextrp: int = -4,
    fextr2: list[float] | np.ndarray | None = None,
    bias: list[float] | np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Obs vertical profile controlled by IEXTRP.

    * ``IEXTRP = ±1``: surface wind in layer 0 only; aloft from UA sounding.
    * ``IEXTRP = ±2`` or ``-4`` (golden-safe): power-law speed + UA direction blend.
    * ``IEXTRP = ±3``: apply FEXTR2 layer factors to surface wind.
    * ``IEXTRP = +4``: Van Ulden–Holtslag SIMILT below Zi; UA above.

    Negative IEXTRP applies optional layer ``BIAS`` (additive m/s on speed).
    """
    from . import similt as _similt

    zmid = layer_mids(zface)
    nz = len(zmid)
    z_agl = np.array([lev.height - stn_elev for lev in sounding_levels], dtype=np.float64)
    wd = np.array([lev.wd for lev in sounding_levels], dtype=np.float64)
    ws = np.array([lev.ws for lev in sounding_levels], dtype=np.float64)
    order = np.argsort(z_agl)
    z_agl, wd, ws = z_agl[order], wd[order], ws[order]
    mask = z_agl > 0
    z_agl, wd, ws = z_agl[mask], wd[mask], ws[mask]
    uu, vv = wind_uv(wd, ws)
    ws1 = float(np.hypot(u_sfc, v_sfc))
    u_s = float(u_sfc / max(ws1, 1e-6))
    v_s = float(v_sfc / max(ws1, 1e-6))
    U = np.zeros((nz, ny, nx))
 