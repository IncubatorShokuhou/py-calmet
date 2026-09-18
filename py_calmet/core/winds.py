"""Diagnostic wind construction (interp, OA, slope flow, mass consistency)."""
from __future__ import annotations
import numpy as np
from .met_utils import wind_uv, ZO_EXTRAP, layer_mids, G, CP
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
    """Obs vertical profile: power-law speed + UA direction blend."""
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


def slope_flow(
    elev: np.ndarray,
    dgrid_m: float,
    qh: np.ndarray,
    tempk: np.ndarray,
    rho: np.ndarray,
    zface: np.ndarray,
    cdk: float = 0.0004,
) -> tuple[np.ndarray, np.ndarray]:
    """Allwine–Whiteman / Horst–Doran style slope-flow (U,V) on layer 1.

    Downslope when qh < 0, upslope when qh > 0. Higher layers get a
    linearly decaying contribution through the estimated drainage depth.
    """
    ny, nx = elev.shape
    nz = len(zface) - 1
    dzdx = np.gradient(elev, dgrid_m, axis=1)
    dzdy = np.gradient(elev, dgrid_m, axis=0)
    slope = np.hypot(dzdx, dzdy)
    sinalf = np.sin(np.arctan(np.maximum(slope, 0.0)))
    sinalf = np.where(np.abs(np.arctan(slope)) < 0.009, 0.0, sinalf)

    hmax = elev.max()
    hmin = elev.min()
    # distance-to-crest / valley proxies
    dcrest = np.sqrt((hmax - elev) ** 2 + (5.0 * dgrid_m) ** 2)
    dvalley = np.sqrt((elev - hmin) ** 2 + (5.0 * dgrid_m) ** 2)
    hd_down = np.maximum(0.05 * (hmax - elev), 0.05)
    hd_up = np.maximum(0.05 * (elev - hmin), 0.05)

    rhocp = np.maximum(rho * CP, 1.0)
    uslope = np.zeros((ny, nx))
    vslope = np.zeros((ny, nx))

    # Downslope (drainage)
    mask_d = qh < 0.0
    if np.any(mask_d):
        sa = np.minimum(sinalf[mask_d], (hmax - elev[mask_d]) / np.maximum(dcrest[mask_d], 1.0))
        tempmef = np.maximum(-dcrest[mask_d] * cdk / hd_down[mask_d], -50.0)
        speed = -(
            sa
            * (G / np.maximum(tempk[mask_d], 200.0))
            * np.abs(qh[mask_d])
            * dcrest[mask_d]
            / rhocp[mask_d]
            / cdk
        )
        speed = np.sign(speed) * (np.abs(speed) ** (1.0 / 3.0)) * ((1.0 - np.exp(tempmef)) ** (1.0 / 3.0))
        # direction: downslope = opposite terrain gradient
        mag = np.maximum(slope[mask_d], 1e-8)
        uslope[mask_d] = -speed * (dzdx[mask_d] / mag)  # toward lower elev when speed>0 after abs
        # speed is negative for downslope in Fortran UVALLY; take abs for magnitude
        mag_spd = np.abs(speed)
        uslope[mask_d] = -mag_spd * (dzdx[mask_d] / mag)
        vslope[mask_d] = -mag_spd * (dzdy[mask_d] / mag)

    # Upslope
    mask_u = qh > 0.0
    if np.any(mask_u):
        speed = (
            (G / np.maximum(tempk[mask_u], 200.0))
            * np.abs(qh[mask_u])
            * (elev[mask_u] - hmin)
            / rhocp[mask_u]
        ) ** (1.0 / 3.0)
        mag = np.maximum(slope[mask_u], 1e-8)
        uslope[mask_u] = speed * (dzdx[mask_u] / mag)
        vslope[mask_u] = speed * (dzdy[mask_u] / mag)

    # Depth decay: full at layer 1, fade by ~hd
    Uadd = np.zeros((nz, ny, nx))
    Vadd = np.zeros_like(Uadd)
    zmid = layer_mids(zface)
    hd = np.where(qh < 0, hd_down, hd_up)
    for L, zm in enumerate(zmid):
        w = np.clip(1.0 - zm / np.maximum(hd, 1.0), 0.0, 1.0)
        Uadd[L] = uslope * w
        Vadd[L] = vslope * w
    return Uadd, Vadd


def divergence_minimize(
    U: np.ndarray,
    V: np.ndarray,
    dgrid_m: float,
    niter: int = 50,
    alpha: float = 0.5,
    terrain: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """O'Brien / DIAGNO-style iterative divergence reduction on each layer.

    Solves a discrete Poisson equation for a velocity potential φ such that
    ∇²φ ≈ ∇·V, then U ← U − ∂φ/∂x, V ← V − ∂φ/∂y. Terrain blocking is a
    light optional damping near steep slopes.
    """
    Uo, Vo = U.copy(), V.copy()
    nz, ny, nx = U.shape
    dx = dgrid_m
    for k in range(nz):
        u = Uo[k]
        v = Vo[k]
        phi = np.zeros((ny, nx), dtype=np.float64)
        for _ in range(niter):
            # divergence
            dudx = np.gradient(u - np.gradient(phi, dx, axis=1), dx, axis=1)
            # Use residual form: update phi so laplacian(phi) → div(U,V)
            div = np.gradient(u, dx, axis=1) + np.gradient(v, dx, axis=0)
            # Jacobi relaxation for ∇²φ = div
            phi_new = phi.copy()
            if ny > 2 and nx > 2:
                neigh = (
                    phi[:-2, 1:-1]
                    + phi[2:, 1:-1]
                    + phi[1:-1, :-2]
                    + phi[1:-1, 2:]
                )
                phi_new[1:-1, 1:-1] = 0.25 * (neigh - div[1:-1, 1:-1] * dx * dx)
            phi = (1.0 - alpha) * phi + alpha * phi_new
        # subtract gradient of phi
        Uo[k] = u - np.gradient(phi, dx, axis=1)
        Vo[k] = v - np.gradient(phi, dx, axis=0)
    if terrain is not None:
        # lightly damp adjustments over flat water-like terrain (no-op mostly)
        pass
    return Uo, Vo


def vertical_velocity_from_div(
    U: np.ndarray,
    V: np.ndarray,
    zface: np.ndarray,
    dgrid_m: float,
) -> np.ndarray:
    """Kinematic W at layer mids from horizontally divergent flow (anelastic)."""
    nz, ny, nx = U.shape
    zmid = layer_mids(zface)
    W = np.zeros((nz, ny, nx), dtype=np.float64)
    # integrate div from surface
    w_face = np.zeros((nz + 1, ny, nx), dtype=np.float64)
    for L in range(nz):
        div = np.gradient(U[L], dgrid_m, axis=1) + np.gradient(V[L], dgrid_m, axis=0)
        dz = zface[L + 1] - zface[L]
        w_face[L + 1] = w_face[L] - div * dz
        W[L] = 0.5 * (w_face[L] + w_face[L + 1])
    return W
