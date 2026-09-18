"""Overwater bulk-flux approximation (COARE-lite).

Full Fairall COARE 3.0 is large; this module provides a documented SciPy/NumPy
bulk-flux scheme bound to ICOARE / ITWPROG / water LU flags. It is intentionally
simpler than COARE but produces stable ustar / sensible heat / Zi over water.
"""
from __future__ import annotations

import numpy as np

from .met_utils import VK, CP, G


def _psi_m_stable(zL: np.ndarray) -> np.ndarray:
    """Stable momentum stability function (Beljaars–Holtslag-ish)."""
    zL = np.asarray(zL, dtype=np.float64)
    c = np.minimum(50.0, 0.35 * zL)
    return -((1.0 + zL) + 0.6667 * (zL - 14.28) / np.exp(c) + 8.525)


def _psi_m_unstable(zL: np.ndarray) -> np.ndarray:
    """Unstable Dyer (16) PSI_m."""
    zL = np.asarray(zL, dtype=np.float64)
    x = (1.0 - 16.0 * zL) ** 0.25
    return (
        2.0 * np.log((1.0 + x) / 2.0)
        + np.log((1.0 + x * x) / 2.0)
        - 2.0 * np.arctan(x)
        + 0.5 * np.pi
    )


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
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bulk overwater fluxes: (ustar, qh_buoyancy, el).

    Uses Charnock-ish z0 and iterative Monin–Obukhov. ``dshelf_km`` lightly
    scales Cd toward coastal values when > 0 (ICOARE shelf hook).
    """
    ws = np.maximum(np.hypot(u10, v10), 0.5)
    t_air = np.asarray(t_air, dtype=np.float64)
    t_sea = np.asarray(t_sea, dtype=np.float64)
    rho_a = np.asarray(rho, dtype=np.float64)
    z0 = np.maximum(np.asarray(z0_water, dtype=np.float64), 1e-5)
    # Charnock update seed
    ustar = 0.035 * ws
    el = np.full(ws.shape, -1000.0, dtype=np.float64)
    for _ in range(5):
        zL = z_ref / np.where(np.abs(el) < 1.0, np.sign(el) * 1.0, el)
        psi = np.where(zL < 0.0, _psi_m_unstable(zL), _psi_m_stable(np.maximum(zL, 0.0)))
        cdn = (VK / (np.log(z_ref / z0) - psi)) ** 2
        if dshelf_km and dshelf_km > 0.0:
            # Mild coastal enhancement (documented approximation)
            cdn = cdn * (1.0 + 0.15 * np.tanh(10.0 / max(dshelf_km, 0.1)))
        ustar = np.sqrt(cdn) * ws
        # Roughness Charnock
        z0 = np.maximum(0.011 * ustar**2 / G + 0.11 * 1.5e-5 / np.maximum(ustar, 1e-3), 1e-5)
        dt = t_air - t_sea
        # Bulk sensible (positive upward when sea warmer)
        ch = 1.1e-3
        qh = rho_a * CP * ch * ws * (t_sea - t_air)  # sea→air when t_sea > t_air
        # MO length from buoyancy flux ≈ qh
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
) -> np.ndarray:
    """Overwater mixing height (Venkatram-like with CONSTW)."""
    f = np.maximum(np.asarray(fcori, dtype=np.float64), 1e-5)
    zi = constw * ustar / f
    # Stable enhancement
    el_pos = np.maximum(el, 0.0)
    zi2 = np.where(el_pos > 0.0, 0.4 * np.sqrt(ustar * el_pos / f), zi)
    zi = np.minimum(zi, zi2)
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
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Replace land PBL params with COARE-lite over water LU cells when ICOARE≠0."""
    if int(icoare) == 0:
        return ustar, el, qh, zi
    water = (landuse >= iwat1) & (landuse <= iwat2)
    if not np.any(water):
        return ustar, el, qh, zi
    sst = tempk if t_sea is None else t_sea
    u_w, qh_w, el_w = coare_lite_fluxes(
        u10, v10, tempk, sst, rho=rho, dshelf_km=dshelf
    )
    zi_w = mixht_overwater(u_w, el_w, fcori, constw, ziminw, zimaxw)
    ustar = np.where(water, u_w, ustar)
    el = np.where(water, el_w, el)
    qh = np.where(water, qh_w, qh)
    zi = np.where(water, zi_w, zi)
    return ustar, el, qh, zi
