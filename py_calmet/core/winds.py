"""Diagnostic wind construction for the tiny-domain configurations."""
from __future__ import annotations
import numpy as np
from .met_utils import wind_uv, ZO_EXTRAP, layer_mids
from .similt import similt_profile


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
    """Map 3D.DAT winds onto CALMET layers (RDMM5-style log below first level)."""
    zmid = layer_mids(zface)
    nz = len(zmid)
    U = np.zeros((nz, ny, nx), dtype=np.float64)
    V = np.zeros_like(U)
    t = hour_index
    for j in range(ny):
        for i in range(nx):
            ii = min(i + 1, threed.ni - 1)
            jj = min(j + 1, threed.nj - 1)
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
) -> tuple[np.ndarray, np.ndarray]:
    """Obs vertical profile: power-law speed + UA direction blend.

    Layer 1 matches the anemometer wind. Aloft, speed follows a stable-layer
    power law gently blended toward upper-air speeds; direction blends from
    the surface toward the sounding (IEXTRP=-4 style approximation).
    """
    zmid = layer_mids(zface)
    z_agl = np.array([lev.height - stn_elev for lev in sounding_levels], dtype=np.float64)
    wd = np.array([lev.wd for lev in sounding_levels], dtype=np.float64)
    ws = np.array([lev.ws for lev in sounding_levels], dtype=np.float64)
    order = np.argsort(z_agl)
    z_agl, wd, ws = z_agl[order], wd[order], ws[order]
    mask = z_agl > 0
    z_agl, wd, ws = z_agl[mask], wd[mask], ws[mask]
    uu, vv = wind_uv(wd, ws)
    ws1 = float(np.hypot(u_sfc, v_sfc))
    wd_sfc = float(np.rad2deg(np.arctan2(-u_sfc, -v_sfc)) % 360.0)
    U = np.zeros((len(zmid), ny, nx))
    V = np.zeros_like(U)
    for L, zm in enumerate(zmid):
        if L == 0:
            U[L], V[L] = u_sfc, v_sfc
            continue
        spd = ws1 * (zm / z_anem) ** p_exp
        u_ua = float(np.interp(zm, z_agl, uu))
        v_ua = float(np.interp(zm, z_agl, vv))
        spd_ua = float(np.hypot(u_ua, v_ua))
        w = min(1.0, np.log(zm / z_anem) / np.log(80.0))
        spd = (1.0 - 0.4 * w) * spd + 0.4 * w * spd_ua
        wd_a = float(np.interp(zm, z_agl, wd))
        wdir = (1.0 - w) * wd_sfc + w * wd_a
        u, v = wind_uv(wdir, spd)
        U[L] = u
        V[L] = v
    return U, V


def objective_analyze(
    ug: np.ndarray,
    vg: np.ndarray,
    u_obs: np.ndarray,
    v_obs: np.ndarray,
    xs_m: float,
    ys_m: float,
    xorig_m: float,
    yorig_m: float,
    dgrid_m: float,
    r1_m: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Single-station Barnes-like OA of obs onto IGF (IPROG=14 style)."""
    nz, ny, nx = ug.shape
    U = ug.copy()
    V = vg.copy()
    for j in range(ny):
        for i in range(nx):
            xc = xorig_m + (i + 0.5) * dgrid_m
            yc = yorig_m + (j + 0.5) * dgrid_m
            r2 = (xc - xs_m) ** 2 + (yc - ys_m) ** 2
            w = np.exp(-r2 / max(r1_m, 1.0) ** 2)
            for k in range(nz):
                U[k, j, i] = (1 - w) * ug[k, j, i] + w * u_obs[k, j, i]
                V[k, j, i] = (1 - w) * vg[k, j, i] + w * v_obs[k, j, i]
    return U, V


def light_terrain_adjust(
    U: np.ndarray,
    V: np.ndarray,
    elev: np.ndarray,
    dgrid_m: float,
    alpha: float = 0.1,
) -> tuple[np.ndarray, np.ndarray]:
    """Light kinematic tilt of near-surface wind along terrain gradient."""
    Uo, Vo = U.copy(), V.copy()
    dzdx = np.gradient(elev, dgrid_m, axis=1)
    dzdy = np.gradient(elev, dgrid_m, axis=0)
    for k in range(min(2, U.shape[0])):
        scale = alpha * (1.0 - k / 3.0)
        speed = np.hypot(Uo[k], Vo[k])
        Uo[k] = Uo[k] - scale * dzdx * speed
        Vo[k] = Vo[k] - scale * dzdy * speed
    return Uo, Vo
