"""Overwater bulk-flux approximation (COARE-lite + warm-layer / cool-skin hooks).

Full Fairall COARE 3.0 is large; this module provides a documented SciPy/NumPy
bulk-flux scheme bound to ICOARE / ITWPROG / IWARM / ICOOL / SEA.DAT. It is
intentionally simpler than COARE but produces stable ustar / sensible heat / Zi.
"""
from __future__ import annotations

import numpy as np

from .met_utils import VK, CP, G


def _psi_m_stable(zL: np.ndarray) -> np.ndarray:
    zL = np.asarray(zL, dtype=np.float64)
    c = np.minimum(50.0, 0.35 * zL)
    return -((1.0 + zL) + 0.6667 * (zL - 14.28) / np.exp(c) + 8.525)


def _psi_m_unstable(zL: np.ndarray) -> np.ndarray:
    zL = np.asarray(zL, dtype=np.float64)
    x = (1.0 - 16.0 * zL) ** 0.25
    return (
        2.0 * np.log((1.0 + x) / 2.0)
        + np.log((1.0 + x * x) / 2.0)
        - 2.0 * np.arctan(x)
        + 0.5 * np.pi
    )


def cool_skin_delta(
    qnet: np.ndarray | float,
    ustar: np.ndarray,
    *,
    rho: float | np.ndarray = 1025.0,
    cp_w: float = 4000.0,
) -> np.ndarray:
    """Fairall-ish cool-skin ΔT (K, skin colder → positive ΔT = Tbulk−Tskin).

    Lightweight: ΔT ≈ min(0.3, 0.02 * max(Qnet,0) / (ρ cp u*)).
    """
    q = np.maximum(np.asarray(qnet, dtype=np.float64), 0.0)
    us = np.maximum(np.asarray(ustar, dtype=np.float64), 1e-3)
    dT = 0.02 * q / (np.asarray(rho, dtype=np.float64) * cp_w * us)
    return np.clip(dT, 0.0, 0.5)


def warm_layer_delta(
    qsw: np.ndarray | float,
    ws: np.ndarray,
    *,
    iwarm: int = 1,
) -> np.ndarray:
    """Diurnal warm-layer ΔT (K, surface warmer than bulk).

    Lightweight: ΔT ≈ 0.5 * (QSW/1000) * exp(−WS/8) when IWARM≠0.
    """
    if int(iwarm) == 0:
        return np.zeros_like(np.asarray(ws, dtype=np.float64))
    q = np.maximum(np.asarray(qsw, dtype=np.float64), 0.0)
    w = np.maximum(np.asarray(ws, dtype=np.float64), 0.1)
    return np.clip(0.5 * (q / 1000.0) * np.exp(-w / 8.0), 0.0, 2.0)


def adjust_sst_skin(
    t_sea_bulk: np.ndarray,
    *,
    qsw: np.ndarray | float = 0.0,
    qnet: np.ndarray | float | None = None,
    ustar: np.ndarray | None = None,
    ws: np.ndarray | None = None,
    iwarm: int = 0,
    icool: int = 0,
) -> np.ndarray:
    """Apply IWARM warm-layer and ICOOL cool-skin to bulk SST → skin temperature."""
    t = np.asarray(t_sea_bulk, dtype=np.float64).copy()
    if ws is None:
        ws = np.full(t.shape, 5.0)
    if int(iwarm) != 0:
        t = t + warm_layer_delta(qsw, ws, iwarm=iwarm)
    if int(icool) != 0:
        us = ustar if ustar is not None else 0.2 * np.ones_like(t)
        qn = qnet if qnet is not None else np.maximum(np.asarray(qsw, dtype=np.float64) - 50.0, 0.0)
        t = t - cool_skin_delta(qn, us)
    return t


def coare_lite_fluxes(
    u10: np.ndarray,
    v10: np.ndarray,
    t_air: np.ndarray,
    t_sea: np.ndarray,
    z0_water: float | np.ndarray = 0.0001,
    z_ref: float = 10.0,
    rho: np.ndarray | float = 1.2,
    *,
    dshelf_km: float = 0.0,
    twave: float | np.ndarray | None = None,
    hwave: float | np.ndarray | None = None,
    iwarm: int = 0,
    icool: int = 0,
    qsw: np.ndarray | float = 0.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bulk overwater fluxes: (ustar, qh, el).

    Optional ``twave``/``hwave`` from SEA.DAT lightly adjust Charnock z0
    (ICOARE wave methods 1/2 lite). IWARM/ICOOL adjust skin SST before flux.
    """
    ws = np.maximum(np.hypot(u10, v10), 0.5)
    t_air = np.asarray(t_air, dtype=np.float64)
    t_sea = np.asarray(t_sea, dtype=np.float64)
    rho_a = np.asarray(rho, dtype=np.float64)
    z0 = np.maximum(np.asarray(z0_water, dtype=np.float64), 1e-5)
    ustar = 0.035 * ws
    # Pre-adjust SST with warm/cool using seed ustar
    t_skin = adjust_sst_skin(
        t_sea, qsw=qsw, ustar=ustar, ws=ws, iwarm=iwarm, icool=icool
    )
    el = np.full(ws.shape, -1000.0, dtype=np.float64)
    for _ in range(5):
        zL = z_ref / np.where(np.abs(el) < 1.0, np.sign(el) * 1.0, el)
        psi = np.where(zL < 0.0, _psi_m_unstable(zL), _psi_m_stable(np.maximum(zL, 0.0)))
        cdn = (VK / (np.log(z_ref / z0) - psi)) ** 2
        if dshelf_km and dshelf_km > 0.0:
            cdn = cdn * (1.0 + 0.15 * np.tanh(10.0 / max(dshelf_km, 0.1)))
        ustar = np.sqrt(cdn) * ws
        # Charnock + optional wave enhancement
        z0 = 0.011 * ustar**2 / G + 0.11 * 1.5e-5 / np.maximum(ustar, 1e-3)
        if twave is not None and hwave is not None:
            tw = np.asarray(twave, dtype=np.float64)
            hw = np.asarray(hwave, dtype=np.float64)
            ok = (tw > 0.0) & (hw > 0.0) & (tw < 9000.0) & (hw < 9000.0)
            # Taylor–Yelland-ish: z0 ~ 1200 * H * (H/Lp)^4.5  (lite scale)
            lp = 1.56 * tw**2
            z0_w = 1200.0 * hw * (hw / np.maximum(lp, 1.0)) ** 4.5
            z0 = np.where(ok, np.maximum(z0, np.minimum(z0_w, 0.05)), z0)
        z0 = np.maximum(z0, 1e-5)
        t_skin = adjust_sst_skin(
            t_sea, qsw=qsw, ustar=ustar, ws=ws, iwarm=iwarm, icool=icool
        )
        ch = 1.1e-3
        qh = rho_a * CP * ch * ws * (t_skin - t_air)
        qh_safe = np.where(np.abs(qh) < 1e-6, np.where(qh >= 0, 1e-6, -1e-6), qh)
        el = -rho_a * CP * t_air * ustar**3 / (VK * G * qh_safe)
        el = np.clip(el, -1.0e5, 1.0e5)
    ustar = np.maximum(ustar, 0.01)
    return ustar, qh, el


def mixht_overwater(
    ustar: np.ndarray,
    el: np.ndarray,
    fcori: float | np.ndarray,
    constw: float = 0.16,
    ziminw: float = 50.0,
    zimaxw: float = 3000.0,
    *,
    qh: np.ndarray | None = None,
    threshw: float = 0.05,
) -> np.ndarray:
    """Overwater mixing height (Venkatram-like with CONSTW + optional THRESHW)."""
    f = np.maximum(np.asarray(fcori, dtype=np.float64), 1e-5)
    zi = constw * ustar / f
    el_pos = np.maximum(el, 0.0)
    zi2 = np.where(el_pos > 0.0, 0.4 * np.sqrt(ustar * el_pos / f), zi)
    zi = np.minimum(zi, zi2)
    if qh is not None and threshw and threshw > 0:
        # Mild convective boost when upward buoyancy exceeds threshw (W/m²/m * Zi proxy)
        conv = np.asarray(qh, dtype=np.float64) > 0.0
        zi = np.where(conv, np.maximum(zi, constw * ustar / f * (1.0 + 0.1)), zi)
    return np.clip(zi, ziminw, zimaxw)


def apply_overwater_pbl(
    landuse: np.ndarray,
    iwat1: int,
    iwat2: int,
    u10: np.ndarray,
    v10: np.ndarray,
    tempk: np.ndarray,
    rho: np.ndarray,
    ustar: np.ndarray,
    el: np.ndarray,
    qh: np.ndarray,
    zi: np.ndarray,
    fcori: float,
    *,
    icoare: int = 0,
    t_sea: np.ndarray | None = None,
    constw: float = 0.16,
    ziminw: float = 50.0,
    zimaxw: float = 3000.0,
    dshelf: float = 0.0,
    iwarm: int = 0,
    icool: int = 0,
    qsw: np.ndarray | float = 0.0,
    twave: float | np.ndarray | None = None,
    hwave: float | np.ndarray | None = None,
    threshw: float = 0.05,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Replace land PBL params with COARE-lite over water LU cells when ICOARE≠0."""
    if int(icoare) == 0:
        return ustar, el, qh, zi
    water = (landuse >= iwat1) & (landuse <= iwat2)
    if not np.any(water):
        return ustar, el, qh, zi
    sst = tempk if t_sea is None else t_sea
    u_w, qh_w, el_w = coare_lite_fluxes(
        u10, v10, tempk, sst, rho=rho, dshelf_km=dshelf,
        twave=twave, hwave=hwave, iwarm=iwarm, icool=icool, qsw=qsw,
    )
    zi_w = mixht_overwater(
        u_w, el_w, fcori, constw, ziminw, zimaxw, qh=qh_w, threshw=threshw
    )
    ustar = np.where(water, u_w, ustar)
    el = np.where(water, el_w, el)
    qh = np.where(water, qh_w, qh)
    zi = np.where(water, zi_w, zi)
    return ustar, el, qh, zi


def sea_sst_grid(
    sea_records: list,
    nx: int,
    ny: int,
    xorig_km: float,
    yorig_km: float,
    dgrid_km: float,
    tempk_fallback: np.ndarray,
) -> tuple[np.ndarray, float | None, float | None]:
    """Nearest-station SST / wave fields from SEA.DAT records onto the grid.

    Returns (t_sea, twave_or_None, hwave_or_None).
    """
    t_sea = np.asarray(tempk_fallback, dtype=np.float64).copy()
    if not sea_records:
        return t_sea, None, None
    xs = np.array([r.x_km for r in sea_records], dtype=np.float64)
    ys = np.array([r.y_km for r in sea_records], dtype=np.float64)
    sst = np.array([r.t_sea for r in sea_records], dtype=np.float64)
    tws = np.array([getattr(r, "twave", -999.0) for r in sea_records], dtype=np.float64)
    hws = np.array([getattr(r, "hwave", -999.0) for r in sea_records], dtype=np.float64)
    tw_grid = np.full((ny, nx), -999.0)
    hw_grid = np.full((ny, nx), -999.0)
    for j in range(ny):
        for i in range(nx):
            xc = xorig_km + (i + 0.5) * dgrid_km
            yc = yorig_km + (j + 0.5) * dgrid_km
            d2 = (xs - xc) ** 2 + (ys - yc) ** 2
            k = int(np.argmin(d2))
            t_sea[j, i] = sst[k]
            tw_grid[j, i] = tws[k]
            hw_grid[j, i] = hws[k]
    tw_out = tw_grid if np.any(tws > 0) else None
    hw_out = hw_grid if np.any(hws > 0) else None
    return t_sea, tw_out, hw_out
