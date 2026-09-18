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
    r2_m: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Single-station Barnes-like OA of obs onto IGF (IPROG=14 style).

    Surface (layer 0) uses R1; aloft layers use R2 (defaults to R1).
    """
    nz, ny, nx = ug.shape
    U = ug.copy()
    V = vg.copy()
    if r2_m is None:
        r2_m = r1_m
    for j in range(ny):
        for i in range(nx):
            xc = xorig_m + (i + 0.5) * dgrid_m
            yc = yorig_m + (j + 0.5) * dgrid_m
            r2 = (xc - xs_m) ** 2 + (yc - ys_m) ** 2
            for k in range(nz):
                rk = r1_m if k == 0 else r2_m
                w = np.exp(-r2 / max(rk, 1.0) ** 2)
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
    """Light kinematic tilt of near-surface wind along terrain gradient (IKINE path)."""
    Uo, Vo = U.copy(), V.copy()
    dzdx = np.gradient(elev, dgrid_m, axis=1)
    dzdy = np.gradient(elev, dgrid_m, axis=0)
    for k in range(min(2, U.shape[0])):
        scale = alpha * (1.0 - k / 3.0)
        speed = np.hypot(Uo[k], Vo[k])
        Uo[k] = Uo[k] - scale * dzdx * speed
        Vo[k] = Vo[k] - scale * dzdy * speed
    return Uo, Vo


def terrain_radius_stats(
    elev: np.ndarray,
    dgrid_m: float,
    terrad_km: float = 5.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """TERSET-like HMAX/HMIN and birdwise distances within TERRAD (km)."""
    ny, nx = elev.shape
    dx_km = dgrid_m * 0.001
    irange = max(int(round(terrad_km / max(dx_km, 1e-6))), 0)
    hmax = elev.copy()
    hmin = elev.copy()
    disthmax = np.zeros((ny, nx), dtype=np.float64)
    disthmin = np.zeros((ny, nx), dtype=np.float64)
    for j in range(ny):
        l1 = max(0, j - irange)
        l2 = min(ny, j + irange + 1)
        for i in range(nx):
            k1 = max(0, i - irange)
            k2 = min(nx, i + irange + 1)
            for jj in range(l1, l2):
                for ii in range(k1, k2):
                    dist_km = np.hypot((ii - i) * dx_km, (jj - j) * dx_km)
                    if dist_km > terrad_km:
                        continue
                    if elev[jj, ii] > hmax[j, i]:
                        hmax[j, i] = elev[jj, ii]
                        disthmax[j, i] = dist_km * 1000.0
                    if elev[jj, ii] < hmin[j, i]:
                        hmin[j, i] = elev[jj, ii]
                        disthmin[j, i] = dist_km * 1000.0
    return hmax, hmin, disthmax, disthmin


def _drainage_uv(delhi: float, delhj: float, uvally: float) -> tuple[float, float]:
    """Resolve Mahrt slope-flow speed into U,V (Fortran SLOPE angle convention)."""
    if delhi == 0.0 and delhj == 0.0:
        return 0.0, 0.0
    if delhi == 0.0:
        thet = 270.0 if delhj < 0.0 else 90.0
    else:
        thetp = np.degrees(np.arctan(delhj / delhi))
        if delhi < 0.0:
            thet = thetp + 180.0
        elif delhj > 0.0:
            thet = thetp
        else:
            thet = thetp + 360.0
    if 0.0 <= thet <= 90.0:
        thetd = 90.0 - thet
    else:
        thetd = 450.0 - thet
    ang = np.radians(270.0 - thetd)
    return float(-np.cos(ang) * uvally), float(-np.sin(ang) * uvally)


def slope_flow(
    elev: np.ndarray,
    dgrid_m: float,
    qh: np.ndarray,
    tempk: np.ndarray,
    rho: np.ndarray,
    zface: np.ndarray,
    landuse: np.ndarray | None = None,
    terrad_km: float = 5.0,
    iwat1: int = 55,
    iwat2: int = 55,
    cdk: float = 0.08,
) -> tuple[np.ndarray, np.ndarray]:
    """Mahrt / Horst–Doran slope flow aligned with CALMET SLOPE (cdk=0.08)."""
    ny, nx = elev.shape
    nz = len(zface) - 1
    hmax, hmin, disthmax, disthmin = terrain_radius_stats(elev, dgrid_m, terrad_km)
    rhocp = 1229.9
    Uadd = np.zeros((nz, ny, nx), dtype=np.float64)
    Vadd = np.zeros_like(Uadd)

    for j in range(ny):
        for i in range(nx):
            if landuse is not None and iwat1 <= int(landuse[j, i]) <= iwat2:
                continue
            im1 = elev[j, i - 1] if i > 0 else elev[j, i]
            ip1 = elev[j, i + 1] if i < nx - 1 else elev[j, i]
            jm1 = elev[j - 1, i] if j > 0 else elev[j, i]
            jp1 = elev[j + 1, i] if j < ny - 1 else elev[j, i]
            delhi = (ip1 - im1) * (0.5 / dgrid_m)
            delhj = (jp1 - jm1) * (0.5 / dgrid_m)
            aalpha = np.arctan(np.hypot(delhi, delhj))
            if abs(aalpha) < 0.009:
                continue
            sinalf = float(np.sin(abs(aalpha)))
            temp = float(max(tempk[j, i], 200.0))
            qhij = float(qh[j, i])
            if qhij < 0.0:
                dcrest = np.sqrt((hmax[j, i] - elev[j, i]) ** 2 + disthmax[j, i] ** 2)
                hd = 0.05 * (hmax[j, i] - elev[j, i])
                if hd <= 0.05:
                    continue
                sinalf = min(sinalf, (hmax[j, i] - elev[j, i]) / max(dcrest, 1e-6))
                tempmef = max(-dcrest * cdk / hd, -50.0)
                uvally = -(
                    (sinalf * (9.81 / temp) * abs(qhij) * dcrest / rhocp / cdk) ** (1.0 / 3.0)
                    * (1.0 - np.exp(tempmef)) ** (1.0 / 3.0)
                )
            else:
                hd = 0.05 * (elev[j, i] - hmin[j, i])
                if hd <= 0.05:
                    continue
                uvally = (
                    (9.81 / temp) * abs(qhij) * (elev[j, i] - hmin[j, i]) / rhocp
                ) ** (1.0 / 3.0)
            delu, delv = _drainage_uv(delhi, delhj, uvally)
            for k in range(nz):
                if hd > zface[k + 1]:
                    Uadd[k, j, i] = delu
                    Vadd[k, j, i] = delv
                elif hd > zface[k]:
                    ratio = (hd - zface[k]) / (zface[k + 1] - zface[k])
                    Uadd[k, j, i] = delu * ratio
                    Vadd[k, j, i] = delv * ratio
    return Uadd, Vadd


def froude_adjust(
    U: np.ndarray,
    V: np.ndarray,
    elev: np.ndarray,
    zface: np.ndarray,
    tempk: np.ndarray,
    dgrid_m: float,
    gamma: np.ndarray | float = 0.01,
    critfn: float = 1.0,
    terrad_km: float = 5.0,
) -> tuple[np.ndarray, np.ndarray]:
    """MELSAR Froude-number blocking adjustment (CALMET FRADJ)."""
    Uo, Vo = U.copy(), V.copy()
    nz, ny, nx = U.shape
    hmax, _, _, _ = terrain_radius_stats(elev, dgrid_m, terrad_km)
    zmid = layer_mids(zface)
    if np.isscalar(gamma):
        gam2d = np.full((ny, nx), float(gamma), dtype=np.float64)
    else:
        gam2d = np.asarray(gamma, dtype=np.float64)
    tau = -0.01
    for j in range(ny):
        for i in range(nx):
            im1 = elev[j, i - 1] if i > 0 else elev[j, i]
            ip1 = elev[j, i + 1] if i < nx - 1 else elev[j, i]
            jm1 = elev[j - 1, i] if j > 0 else elev[j, i]
            jp1 = elev[j + 1, i] if j < ny - 1 else elev[j, i]
            delhi = (ip1 - im1) * (0.5 / dgrid_m)
            delhj = (jp1 - jm1) * (0.5 / dgrid_m)
            if abs(delhi) < 1.5e-4 and abs(delhj) < 1.5e-4:
                delhi = 0.0
                delhj = 0.0
            temp = float(max(tempk[j, i], 200.0))
            gamma2 = float(gam2d[j, i]) - tau
            if gamma2 <= 0.0:
                continue
            for k in range(nz):
                speed = float(np.hypot(Uo[k, j, i], Vo[k, j, i]))
                if speed <= 0.0:
                    continue
                obshgt = hmax[j, i] - (elev[j, i] + zmid[k])
                if obshgt <= 0.0:
                    continue
                froude = speed / (np.sqrt(9.8 * gamma2 / temp) * obshgt)
                if froude > critfn:
                    continue
                # drainage direction of terrain (toward lower elev)
                if delhi == 0.0 and delhj == 0.0:
                    continue
                if delhi == 0.0:
                    thet = 270.0 if delhj < 0.0 else 90.0
                else:
                    thetp = np.degrees(np.arctan(delhj / delhi))
                    if delhi < 0.0:
                        thet = thetp + 180.0
                    elif delhj > 0.0:
                        thet = thetp
                    else:
                        thet = thetp + 360.0
                thetd = 90.0 - thet if 0.0 <= thet <= 90.0 else 450.0 - thet
                ang = np.radians(270.0 - thetd)
                drx = -np.cos(ang)
                dry = -np.sin(ang)
                un = Uo[k, j, i] / speed
                vn = Vo[k, j, i] / speed
                aa = drx * un + dry * vn
                if aa <= 0.0:
                    continue
                tt = -dry * un + drx * vn
                if tt > 0.0:
                    tax, tay = -dry, drx
                else:
                    tax, tay = dry, -drx
                Uo[k, j, i] = tax * speed
                Vo[k, j, i] = tay * speed
    return Uo, Vo


def smooth_winds(
    U: np.ndarray,
    V: np.ndarray,
    nsmth: list[int] | np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """CALMET SMOOTH: 5-point filter with NSMTH passes per layer."""
    Uo, Vo = U.copy(), V.copy()
    nz, ny, nx = U.shape
    if nsmth is None:
        nsmth = [2] + [4] * (nz - 1)
    nsmth = list(nsmth) + [nsmth[-1]] * max(0, nz - len(nsmth))
    for k in range(nz):
        for _ in range(int(nsmth[k])):
            Un = Uo[k].copy()
            Vn = Vo[k].copy()
            for j in range(ny):
                for i in range(nx):
                    uim1 = Uo[k, j, i - 1] if i > 0 else Uo[k, j, i]
                    uip1 = Uo[k, j, i + 1] if i < nx - 1 else Uo[k, j, i]
                    ujm1 = Uo[k, j - 1, i] if j > 0 else Uo[k, j, i]
                    ujp1 = Uo[k, j + 1, i] if j < ny - 1 else Uo[k, j, i]
                    vim1 = Vo[k, j, i - 1] if i > 0 else Vo[k, j, i]
                    vip1 = Vo[k, j, i + 1] if i < nx - 1 else Vo[k, j, i]
                    vjm1 = Vo[k, j - 1, i] if j > 0 else Vo[k, j, i]
                    vjp1 = Vo[k, j + 1, i] if j < ny - 1 else Vo[k, j, i]
                    Un[j, i] = 0.5 * Uo[k, j, i] + 0.125 * (uip1 + uim1 + ujp1 + ujm1)
                    Vn[j, i] = 0.5 * Vo[k, j, i] + 0.125 * (vip1 + vim1 + vjp1 + vjm1)
            Uo[k] = Un
            Vo[k] = Vn
    return Uo, Vo


def divergence_minimize(
    U: np.ndarray,
    V: np.ndarray,
    dgrid_m: float,
    niter: int = 50,
    alpha: float = 0.5,
    terrain: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """O'Brien / DIAGNO-style iterative divergence reduction on each layer.

    Only used when IOBR=1 (or IKINE=1 via TOPOF2+MINIM). Default case configs
    set IOBR=0, so this is typically skipped by the runner.
    """
    Uo, Vo = U.copy(), V.copy()
    nz, ny, nx = U.shape
    dx = dgrid_m
    for k in range(nz):
        u = Uo[k]
        v = Vo[k]
        phi = np.zeros((ny, nx), dtype=np.float64)
        for _ in range(niter):
            div = np.gradient(u, dx, axis=1) + np.gradient(v, dx, axis=0)
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
        Uo[k] = u - np.gradient(phi, dx, axis=1)
        Vo[k] = v - np.gradient(phi, dx, axis=0)
    return Uo, Vo


def vertical_velocity_from_div(
    U: np.ndarray,
    V: np.ndarray,
    zface: np.ndarray,
    dgrid_m: float,
) -> np.ndarray:
    """Kinematic W at layer mids from horizontally divergent flow (anelastic)."""
    nz, ny, nx = U.shape
    W = np.zeros((nz, ny, nx), dtype=np.float64)
    w_face = np.zeros((nz + 1, ny, nx), dtype=np.float64)
    for L in range(nz):
        div = np.gradient(U[L], dgrid_m, axis=1) + np.gradient(V[L], dgrid_m, axis=0)
        dz = zface[L + 1] - zface[L]
        w_face[L + 1] = w_face[L] - div * dz
        W[L] = 0.5 * (w_face[L] + w_face[L + 1])
    return W
