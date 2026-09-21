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
    """Surface U/V from SURF.DAT ws/wd.

    CALMET missing sentinel ``9999`` (and non-finite) → calm zeros, not a
    phantom ~10 km/s wind that would poison OA / profiles.
    """
    ws_f, wd_f = float(ws), float(wd)
    if (
        not (np.isfinite(ws_f) and np.isfinite(wd_f))
        or ws_f >= 9000.0
        or wd_f >= 9000.0
    ):
        z = np.zeros((ny, nx), dtype=np.float64)
        return z, z.copy()
    u, v = wind_uv(wd_f, ws_f)
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
    # UP.DAT missing ≈ 999 (same gate as domain_avg_wind_from_sounding)
    mask = (
        (z_agl > 0)
        & np.isfinite(ws)
        & np.isfinite(wd)
        & (ws < 998.0)
        & (wd < 998.0)
    )
    z_agl, wd, ws = z_agl[mask], wd[mask], ws[mask]
    have_ua = z_agl.size >= 1
    if have_ua:
        uu, vv = wind_uv(wd, ws)
    else:
        uu = vv = np.zeros(0, dtype=np.float64)
    ws1 = float(np.hypot(u_sfc, v_sfc))
    u_s = float(u_sfc / max(ws1, 1e-6))
    v_s = float(v_sfc / max(ws1, 1e-6))
    U = np.zeros((nz, ny, nx))
    V = np.zeros_like(U)
    mode = abs(int(iextrp))

    def _ua_at(zm: float) -> tuple[float, float]:
        if not have_ua:
            return float(u_sfc), float(v_sfc)
        return float(np.interp(zm, z_agl, uu)), float(np.interp(zm, z_agl, vv))

    if mode == 4 and int(iextrp) > 0:
        # True SIMILT (positive IEXTRP=4 only; -4 keeps golden power-law)
        us, vs = _similt.similt_profile(
            u_sfc, v_sfc, z_anem, max(z0, 1e-4), el, zi, zmid, zimin=zimin
        )
        for L, zm in enumerate(zmid):
            if np.isnan(us[L]):
                U[L], V[L] = _ua_at(zm)
            else:
                U[L], V[L] = float(us[L]), float(vs[L])
    elif mode == 1:
        for L, zm in enumerate(zmid):
            if L == 0:
                U[L], V[L] = u_sfc, v_sfc
            else:
                U[L], V[L] = _ua_at(zm)
    elif mode == 3:
        fx = list(fextr2) if fextr2 is not None else [1.0] * nz
        fx = fx + [fx[-1]] * max(0, nz - len(fx))
        for L in range(nz):
            U[L] = u_sfc * float(fx[L])
            V[L] = v_sfc * float(fx[L])
    else:
        # ±2 and -4 (and unknown): golden-tuned power-law blend
        for L, zm in enumerate(zmid):
            if L == 0:
                U[L], V[L] = u_sfc, v_sfc
                continue
            spd = ws1 * (zm / max(z_anem, 1e-3)) ** p_exp
            u_ua, v_ua = _ua_at(zm)
            spd_ua = float(np.hypot(u_ua, v_ua))
            w = min(1.0, np.log(max(zm, z_anem) / max(z_anem, 1e-3)) / np.log(80.0))
            if have_ua:
                spd = (1.0 - 0.4 * w) * spd + 0.4 * w * spd_ua
                u_a = u_ua / max(spd_ua, 1e-6)
                v_a = v_ua / max(spd_ua, 1e-6)
                bu = (1.0 - w) * u_s + w * u_a
                bv = (1.0 - w) * v_s + w * v_a
            else:
                bu, bv = u_s, v_s
            wdir = float(np.rad2deg(np.arctan2(-bu, -bv)) % 360.0)
            u, v = wind_uv(wdir, spd)
            U[L] = u
            V[L] = v

    if int(iextrp) < 0 and bias is not None and len(bias) > 0:
        b = list(bias) + [0.0] * max(0, nz - len(bias))
        for L in range(nz):
            spd = float(np.hypot(U[L, 0, 0], V[L, 0, 0]))
            if spd < 1e-9:
                continue
            scale = (spd + float(b[L])) / spd
            U[L] *= scale
            V[L] *= scale
    return U, V


def objective_analyze(
    ug: np.ndarray,
    vg: np.ndarray,
    u_obs: np.ndarray,
    v_obs: np.ndarray,
    xs_m: float | np.ndarray,
    ys_m: float | np.ndarray,
    xorig_m: float,
    yorig_m: float,
    dgrid_m: float,
    r1_m: float,
    r2_m: float | None = None,
    rprog_m: float = 0.0,
    rmax1_m: float | None = None,
    rmax2_m: float | None = None,
    rmax3_m: float | None = None,
    rmin_m: float = 0.0,
    lvary: bool = False,
    icalm: int = 0,
    is_water: np.ndarray | None = None,
    nintr2: list[int] | np.ndarray | None = None,
    barriers=None,
    xorig_km: float | None = None,
    yorig_km: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Barnes-like OA of obs onto IGF (IPROG=14 style).

    * Single-station: ``xs_m``/``ys_m`` scalars (legacy path; golden-compatible).
    * Multi-station: 1-D ``xs_m``/``ys_m`` with matching station U/V in
      ``u_obs``/``v_obs`` shaped ``(nstn, nz)`` **or** still ``(nz, ny, nx)``
      gridded obs (then stations are sampled at ``xs_m``/``ys_m``).
    * ``rprog_m``: diagnostic influence radius for the IGF/prog field used as
      an extra "station" at every grid point (CALMET RPROG). ``0`` disables.
    Surface (layer 0) uses R1; aloft layers use R2 (defaults to R1).
    """
    nz, ny, nx = ug.shape
    U = ug.copy()
    V = vg.copy()
    if r2_m is None:
        r2_m = r1_m

    xs = np.atleast_1d(np.asarray(xs_m, dtype=np.float64))
    ys = np.atleast_1d(np.asarray(ys_m, dtype=np.float64))
    nstn = int(xs.size)

    # Normalize obs to per-station profiles (nstn, nz)
    u_obs_a = np.asarray(u_obs, dtype=np.float64)
    v_obs_a = np.asarray(v_obs, dtype=np.float64)
    if u_obs_a.ndim == 3 and u_obs_a.shape == (nz, ny, nx) and nstn == 1:
        # legacy: obs already gridded — keep previous formula exactly
        for j in range(ny):
            for i in range(nx):
                xc = xorig_m + (i + 0.5) * dgrid_m
                yc = yorig_m + (j + 0.5) * dgrid_m
                dist2 = (xc - float(xs[0])) ** 2 + (yc - float(ys[0])) ** 2
                for k in range(nz):
                    rk = r1_m if k == 0 else r2_m
                    w = np.exp(-dist2 / max(rk, 1.0) ** 2)
                    if rprog_m and rprog_m > 0.0:
                        # Blend toward IGF with prog weight at every cell
                        wp = np.exp(0.0)  # self-distance 0 → weight 1 * scaled
                        # Effective: obs weight w, prog weight wp_scaled
                        wp = (r1_m / max(rprog_m, 1.0)) ** 2 if k == 0 else (r2_m / max(rprog_m, 1.0)) ** 2
                        # Standard CALMET-ish: final = (w*obs + wp*ug) / (w+wp)
                        # with residual IGF when both small — collapse to legacy when rprog=0
                        denom = w + wp
                        U[k, j, i] = (w * u_obs_a[k, j, i] + wp * ug[k, j, i]) / denom
                        V[k, j, i] = (w * v_obs_a[k, j, i] + wp * vg[k, j, i]) / denom
                    else:
                        U[k, j, i] = (1.0 - w) * ug[k, j, i] + w * u_obs_a[k, j, i]
                        V[k, j, i] = (1.0 - w) * vg[k, j, i] + w * v_obs_a[k, j, i]
        return U, V

    # Multi-station (or station profiles)
    if u_obs_a.ndim == 3 and u_obs_a.shape == (nz, ny, nx):
        # Sample gridded obs at station locations
        u_stn = np.zeros((nstn, nz), dtype=np.float64)
        v_stn = np.zeros((nstn, nz), dtype=np.float64)
        for s in range(nstn):
            ii = int(np.clip(np.floor((xs[s] - xorig_m) / dgrid_m), 0, nx - 1))
            jj = int(np.clip(np.floor((ys[s] - yorig_m) / dgrid_m), 0, ny - 1))
            u_stn[s] = u_obs_a[:, jj, ii]
            v_stn[s] = v_obs_a[:, jj, ii]
    elif u_obs_a.ndim == 2 and u_obs_a.shape[0] == nstn:
        u_stn, v_stn = u_obs_a, v_obs_a
    elif u_obs_a.ndim == 2 and u_obs_a.shape[1] == nstn:
        u_stn, v_stn = u_obs_a.T, v_obs_a.T
    else:
        raise ValueError(
            f"u_obs shape {u_obs_a.shape} incompatible with nstn={nstn}, nz={nz}"
        )

    for j in range(ny):
        for i in range(nx):
            xc = xorig_m + (i + 0.5) * dgrid_m
            yc = yorig_m + (j + 0.5) * dgrid_m
            dist2 = (xc - xs) ** 2 + (yc - ys) ** 2
            for k in range(nz):
                rk = r1_m if k == 0 else r2_m
                rmax = rmax1_m if k == 0 else rmax2_m
                # RMAX3: over-water cells use alternate cutoff when provided
                if (
                    is_water is not None
                    and rmax3_m is not None
                    and rmax3_m > 0
                    and bool(is_water[j, i])
                ):
                    rmax = rmax3_m
                dist = np.sqrt(dist2)
                # RMIN: floor distance to avoid singularity at station
                dist_w = np.maximum(dist, float(rmin_m) if rmin_m and rmin_m > 0 else 0.0)
                w_stn = np.exp(-(dist_w ** 2) / max(rk, 1.0) ** 2)
                if rmax is not None and rmax > 0:
                    w_stn = np.where(dist <= rmax, w_stn, 0.0)
                    # LVARY: expand radius until at least one station weighs in
                    if lvary and not np.any(w_stn > 0) and dist.size:
                        order = np.argsort(dist)
                        keep = order[: max(1, min(3, dist.size))]
                        w_stn = np.zeros_like(w_stn)
                        w_stn[keep] = np.exp(-(dist_w[keep] ** 2) / max(rk, 1.0) ** 2)
                if barriers is not None:
                    from .barriers import station_clear_mask
                    x0k = (xorig_m / 1000.0) if xorig_km is None else xorig_km
                    y0k = (yorig_m / 1000.0) if yorig_km is None else yorig_km
                    clear = station_clear_mask(
                        xc / 1000.0, yc / 1000.0,
                        xs / 1000.0, ys / 1000.0,
                        barriers, layer_k=k,
                    )
                    w_stn = np.where(clear, w_stn, 0.0)
                if nintr2 is not None:
                    nmax = int(nintr2[k]) if k < len(nintr2) else int(nintr2[-1])
                    if nmax > 0 and nmax < w_stn.size:
                        # Keep only the nmax nearest stations
                        order = np.argsort(dist2)
                        mask = np.zeros_like(w_stn, dtype=bool)
                        mask[order[:nmax]] = True
                        w_stn = np.where(mask, w_stn, 0.0)
                num_u = float(np.dot(w_stn, u_stn[:, k]))
                num_v = float(np.dot(w_stn, v_stn[:, k]))
                den = float(w_stn.sum())
                if rprog_m and rprog_m > 0.0:
                    wp = (rk / max(rprog_m, 1.0)) ** 2
                    num_u += wp * ug[k, j, i]
                    num_v += wp * vg[k, j, i]
                    den += wp
                if den > 1e-12:
                    u_oa = num_u / den
                    v_oa = num_v / den
                    if int(icalm) != 0:
                        # ICALM≠0: discard calm obs contributions (ws < 0.5 m/s)
                        spd_oa = (u_oa ** 2 + v_oa ** 2) ** 0.5
                        if spd_oa < 0.5:
                            continue  # keep IGF
                    U[k, j, i] = u_oa
                    V[k, j, i] = v_oa
                # else leave IGF
    return U, V


def topographic_kinematic_w(
    U: np.ndarray,
    V: np.ndarray,
    elev: np.ndarray,
    zface: np.ndarray,
    tempk: np.ndarray,
    dgrid_m: float,
    alpha: float = 0.1,
    gamma: float | np.ndarray = 0.01,
) -> np.ndarray:
    """TOPOF2-style terrain-induced vertical velocity (IKINE=1).

    Computes stability-dependent exponential decay of slope-forced W
    (Yocke 1979 / CALMET TOPOF2), returned at layer midpoints ``(nz,ny,nx)``.
    Does **not** modify U,V — pair with ``divergence_minimize`` / O'Brien.
    """
    nz, ny, nx = U.shape
    if np.isscalar(gamma):
        gam2d = np.full((ny, nx), float(gamma), dtype=np.float64)
    else:
        gam2d = np.asarray(gamma, dtype=np.float64)
    tau = -0.01
    hinv = 500.0
    W_face = np.zeros((nz + 1, ny, nx), dtype=np.float64)
    dzi = 0.5 / dgrid_m
    for j in range(ny):
        for i in range(nx):
            temp = float(max(tempk[j, i], 200.0))
            gamma2 = float(gam2d[j, i]) - tau
            if gamma2 < 0.0:
                s = -1.0
            else:
                s = float(np.sqrt(G * gamma2 / temp))
            xws = float(np.hypot(U[nz - 1, j, i], V[nz - 1, j, i]))
            xws = max(xws, 1e-6)
            bk = 2.0 / hinv if s <= 0.0 else s / xws
            im1 = elev[j, i - 1] if i > 0 else elev[j, i]
            ip1 = elev[j, i + 1] if i < nx - 1 else elev[j, i]
            jm1 = elev[j - 1, i] if j > 0 else elev[j, i]
            jp1 = elev[j + 1, i] if j < ny - 1 else elev[j, i]
            delhi = (ip1 - im1) * dzi
            delhj = (jp1 - jm1) * dzi
            wtopo = float(U[0, j, i] * delhi + V[0, j, i] * delhj)
            w1 = np.zeros(nz + 1, dtype=np.float64)
            w1[0] = wtopo
            w_tf = np.zeros(nz + 1, dtype=np.float64)
            w_tf[0] = wtopo
            for k in range(nz):
                bkz = min(bk * float(zface[k + 1]), 50.0)
                w1[k + 1] = wtopo * np.exp(-bkz)
                dz = float(zface[k + 1] - zface[k])
                dwdz1 = (w1[k + 1] - w1[k]) / max(dz, 1e-6)
                dwdz = alpha * dwdz1
                w_tf[k + 1] = dwdz * dz + w_tf[k]
            # Terrain-following transform
            w_tf[0] = w_tf[0] - U[0, j, i] * delhi - V[0, j, i] * delhj
            for k in range(nz):
                w_tf[k + 1] = w_tf[k + 1] - U[k, j, i] * delhi - V[k, j, i] * delhj
            W_face[:, j, i] = w_tf
    W = np.zeros((nz, ny, nx), dtype=np.float64)
    for k in range(nz):
        W[k] = 0.5 * (W_face[k] + W_face[k + 1])
    return W


def light_terrain_adjust(
    U: np.ndarray,
    V: np.ndarray,
    elev: np.ndarray,
    dgrid_m: float,
    alpha: float = 0.1,
) -> tuple[np.ndarray, np.ndarray]:
    """Legacy near-surface tilt (pre-TOPOF2). Prefer topographic_kinematic_w."""
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
                    (sinalf * (G / temp) * abs(qhij) * dcrest / rhocp / cdk) ** (1.0 / 3.0)
                    * (1.0 - np.exp(tempmef)) ** (1.0 / 3.0)
                )
            else:
                hd = 0.05 * (elev[j, i] - hmin[j, i])
                if hd <= 0.05:
                    continue
                uvally = (
                    (G / temp) * abs(qhij) * (elev[j, i] - hmin[j, i]) / rhocp
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
                froude = speed / (np.sqrt(G * gamma2 / temp) * obshgt)
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
    divlim: float = 5e-6,
    W: np.ndarray | None = None,
    zface: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """O'Brien / MINIM-style iterative divergence reduction (IOBR=1).

    Per-layer Poisson adjustment of U,V until |div| ≤ ``divlim`` or ``niter``.
    When ``W`` and ``zface`` are supplied, vertical derivative ∂W/∂z is included
    in the divergence (anelastic continuity) — the O'Brien coupling path.
    Gated off when IOBR=0 (runner); goldens keep IOBR=0.
    """
    _ = terrain
    Uo, Vo = U.copy(), V.copy()
    nz, ny, nx = U.shape
    dx = dgrid_m
    for k in range(nz):
        u = Uo[k]
        v = Vo[k]
        for it in range(niter):
            div = np.gradient(u, dx, axis=1) + np.gradient(v, dx, axis=0)
            if W is not None and zface is not None:
                dz = float(zface[k + 1] - zface[k])
                if k == 0:
                    dwdz = (W[k] - 0.0) / max(dz, 1e-6)
                else:
                    dwdz = (W[k] - W[k - 1]) / max(dz, 1e-6)
                div = div + dwdz
            divmax = float(np.max(np.abs(div)))
            if divmax <= divlim:
                break
            # Four-pass directional adjustment (CALMET MINIM-inspired)
            for _idir in range(4):
                for j in range(ny):
                    for i in range(nx):
                        d = div[j, i]
                        if abs(d) < 1e-20:
                            continue
                        ut = -0.5 * d * dx  # alpha1..4 = 0.5 → AL=2 → factor 0.5
                        vt = -0.5 * d * dx
                        if i + 1 < nx:
                            u[j, i + 1] += 0.5 * ut
                        if i - 1 >= 0:
                            u[j, i - 1] -= 0.5 * ut
                        if j + 1 < ny:
                            v[j + 1, i] += 0.5 * vt
                        if j - 1 >= 0:
                            v[j - 1, i] -= 0.5 * vt
                div = np.gradient(u, dx, axis=1) + np.gradient(v, dx, axis=0)
                if W is not None and zface is not None:
                    dz = float(zface[k + 1] - zface[k])
                    if k == 0:
                        dwdz = W[k] / max(dz, 1e-6)
                    else:
                        dwdz = (W[k] - W[k - 1]) / max(dz, 1e-6)
                    div = div + dwdz
        Uo[k] = u
        Vo[k] = v
    return Uo, Vo


def obrien_adjust(
    U: np.ndarray,
    V: np.ndarray,
    W: np.ndarray,
    zface: np.ndarray,
    dgrid_m: float,
    niter: int = 50,
    divlim: float = 5e-6,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Full O'Brien path: adjust U,V for continuity including kinematic W, then rebuild W."""
    U2, V2 = divergence_minimize(
        U, V, dgrid_m, niter=niter, divlim=divlim, W=W, zface=zface
    )
    W2 = vertical_velocity_from_div(U2, V2, zface, dgrid_m)
    # Blend top boundary: force W_top → 0 (classic O'Brien)
    nz = U2.shape[0]
    if nz >= 2:
        w_top = W2[-1]
        for k in range(nz):
            fac = float(k) / float(nz - 1) if nz > 1 else 1.0
            W2[k] = W2[k] - fac * w_top
    return U2, V2, W2


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
