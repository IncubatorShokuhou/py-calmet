"""PBL / surface parameterizations (ELUSTR stable branch + MIXHT night)."""
from __future__ import annotations
import numpy as np
from .met_utils import VK, CP, G, coriolis


def air_density(tempk: np.ndarray, pres_mb: float = 1012.0) -> np.ndarray:
    """Ideal-gas density (kg/m^3) from T and surface pressure."""
    # rho = P / (R_specific * T); R_d ≈ 287.05
    return (pres_mb * 100.0) / (287.05 * tempk)


def elustr_stable(
    u1: np.ndarray,
    v1: np.ndarray,
    z0: np.ndarray,
    zm: float,
    tempk: np.ndarray,
    rho: np.ndarray,
    sky_tenths: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Nighttime ELUSTR: returns ustar, el, qh."""
    ws = np.maximum(np.hypot(u1, v1), 0.5)
    xlnzz0 = np.log(zm / np.maximum(z0, 1e-4))
    cdn = VK / xlnzz0
    jcc = np.asarray(sky_tenths, dtype=np.float64)
    theta1 = 0.09 * (1.0 - 5.0e-3 * jcc ** 2)
    theta2 = tempk * cdn * ws ** 2 / (184.428 * zm)
    thetas = np.minimum(theta1, theta2)
    thetas = np.maximum(thetas, 1.0e-9)
    u02 = 46.107 * zm * thetas / tempk
    cu2 = np.maximum(0.0, ws * ws - 4.0 * u02 / cdn)
    ustar = 0.5 * cdn * (ws + np.sqrt(cu2))
    ustar = np.maximum(ustar, 0.05)
    xcrit = 0.05
    for _ in range(3):
        mask = ustar * thetas > xcrit
        if not np.any(mask):
            break
        thetas = np.where(mask, xcrit / ustar, thetas)
        u02 = 46.107 * zm * thetas / tempk
        cu2 = np.maximum(0.0, ws * ws - 4.0 * u02 / cdn)
        ustar = np.maximum(0.5 * cdn * (ws + np.sqrt(cu2)), 0.05)
    qh = -CP * rho * ustar * thetas
    el = -253.8226 * rho * tempk * ustar ** 3 / qh
    return ustar, el, qh


def mixht_night(
    ustar: np.ndarray,
    el: np.ndarray,
    fcori: float | np.ndarray,
    constn: float = 2400.0,
    zimin: float = 50.0,
    zimax: float = 3000.0,
) -> np.ndarray:
    """Stable mechanical mixing height (Venkatram + Zilitinkevich)."""
    zi1 = constn * ustar ** 1.5
    zi2 = 0.4 * np.sqrt(ustar * el / fcori)
    zi = np.minimum(np.minimum(zi1, zi2), zimax)
    zi = np.maximum(zi, zimin)
    return zi


def ipgt_from_el(el: np.ndarray, zimin: float = 50.0) -> np.ndarray:
    """Rough PGT class from MO length (nighttime-leaning)."""
    # Simple mapping used for tiny-domain parity (classes 5–6)
    ipgt = np.full(el.shape, 5, dtype=np.int32)
    ipgt = np.where(el > 200.0, 6, ipgt)
    ipgt = np.where(el < 50.0, 4, ipgt)
    return ipgt


def wstar_field(ziconv: np.ndarray, qh: np.ndarray, tempk: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Convective velocity scale; zero when qh <= 0."""
    out = np.zeros_like(qh)
    mask = qh > 0
    if np.any(mask):
        wt = qh[mask] / (rho[mask] * CP)
        out[mask] = (G / tempk[mask] * wt * np.maximum(ziconv[mask], 1.0)) ** (1.0 / 3.0)
    return out
