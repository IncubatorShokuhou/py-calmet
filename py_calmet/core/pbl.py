"""PBL / surface parameterizations (ELUSTR, MIXHT night + daytime convective)."""
from __future__ import annotations
import numpy as np
from .met_utils import VK, CP, G, coriolis


def air_density(tempk: np.ndarray, pres_mb: float | np.ndarray = 1012.0) -> np.ndarray:
    """Ideal-gas density (kg/m^3) from T and surface pressure."""
    return (np.asarray(pres_mb, dtype=np.float64) * 100.0) / (287.05 * tempk)


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


def elustr_unstable(
    u1: np.ndarray,
    v1: np.ndarray,
    z0: np.ndarray,
    zm: float,
    tempk: np.ndarray,
    rho: np.ndarray,
    qsw: np.ndarray,
    albedo: float = 0.2,
    bowen: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Daytime ELUSTR-ish: positive sensible heat from net shortwave.

    Simplified energy-budget qh = (1-albedo)*qsw / (1+bowen); ustar from
    neutral/unstable CD; L from Monin–Obukhov.
    """
    ws = np.maximum(np.hypot(u1, v1), 0.5)
    xlnzz0 = np.log(zm / np.maximum(z0, 1e-4))
    cdn = VK / xlnzz0
    ustar = np.maximum(cdn * ws, 0.05)
    qh = np.maximum((1.0 - albedo) * qsw / (1.0 + bowen), 1.0)
    # MO length (negative for unstable)
    el = -rho * CP * tempk * ustar ** 3 / (VK * G * qh)
    el = np.minimum(el, -1.0)
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


def mixht_day_carson(
    qh: np.ndarray,
    rho: np.ndarray,
    tempk: np.ndarray,
    ustar: np.ndarray,
    fcori: float | np.ndarray,
    dt_sec: float = 3600.0,
    ziconv_prev: np.ndarray | None = None,
    dtheta: float = 0.01,
    threshl: float = 0.05,
    constb: float = 1.41,
    zimin: float = 50.0,
    zimax: float = 3000.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Daytime mixing height: Carson convective + mechanical (Venkatram 1980b).

    Returns (zi, ziconv). Energy-balance Carson growth:
        d(zi)/dt ≈ (2 * wt) / dtheta   with wt = qh/(rho*cp)
    gated by THRESHL on buoyancy flux per meter.
    """
    wt = qh / (np.maximum(rho, 0.5) * CP)
    # buoyancy flux proxy
    buoy = G / np.maximum(tempk, 200.0) * wt
    htold = np.zeros_like(qh) if ziconv_prev is None else ziconv_prev.copy()
    # Carson: zi_new^2 ≈ zi_old^2 + 2*wt*dt / dtheta  (potential-temp jump)
    dth = max(dtheta, 1e-3)
    grow = np.where(buoy > threshl, 2.0 * wt * dt_sec / dth, 0.0)
    ziconv = np.sqrt(np.maximum(htold ** 2 + grow, 0.0))
    ziconv = np.minimum(ziconv, zimax)

    # Mechanical daytime (neutral): hmech = cmech * ustar / N^0.5 approx
    # Use BVF proxy from dtheta/dz ~ dth/200m
    tave = tempk
    bvf = (G * dth / np.maximum(tave * 200.0, 1.0)) ** 0.25
    cmech = constb / np.sqrt(np.maximum(np.asarray(fcori, dtype=np.float64), 1e-5))
    hmech = cmech * ustar / np.maximum(bvf, 1e-3)

    zi = np.maximum(np.maximum(zimin, hmech), ziconv)
    zi = np.minimum(zi, zimax)
    return zi, ziconv


def ipgt_from_el(el: np.ndarray, zimin: float = 50.0) -> np.ndarray:
    """PGT class from MO length (Golder-like)."""
    ipgt = np.full(el.shape, 4, dtype=np.int32)  # neutral default
    ipgt = np.where(el > 200.0, 6, ipgt)       # F stable
    ipgt = np.where((el > 50.0) & (el <= 200.0), 5, ipgt)  # E
    ipgt = np.where((el < -50.0) & (el >= -200.0), 3, ipgt)  # C
    ipgt = np.where((el < -200.0) & (el >= -1000.0), 2, ipgt)  # B
    ipgt = np.where(el < -1000.0, 1, ipgt)  # A
    return ipgt


def wstar_field(ziconv: np.ndarray, qh: np.ndarray, tempk: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """Convective velocity scale; zero when qh <= 0."""
    out = np.zeros_like(qh)
    mask = qh > 0
    if np.any(mask):
        wt = qh[mask] / (rho[mask] * CP)
        out[mask] = (G / tempk[mask] * wt * np.maximum(ziconv[mask], 1.0)) ** (1.0 / 3.0)
    return out


def relative_humidity_2d(q2_gkg: np.ndarray, tempk: np.ndarray, pres_mb: float | np.ndarray) -> np.ndarray:
    """RH (%) from mixing ratio (g/kg), T, P."""
    qv = np.maximum(q2_gkg / 1000.0, 1e-12)
    p_pa = np.asarray(pres_mb, dtype=np.float64) * 100.0
    e = qv * p_pa / (0.622 + qv)
    tc = tempk - 273.15
    a = np.array(
        [6.107799961, 4.436518521e-1, 1.428945805e-2, 2.650648471e-4,
         3.031240396e-6, 2.034080948e-8, 6.136820929e-11]
    )
    es_mb = a[0] + tc * (a[1] + tc * (a[2] + tc * (a[3] + tc * (a[4] + tc * (a[5] + tc * a[6])))))
    es = np.maximum(es_mb, 0.01) * 100.0
    return np.clip(100.0 * e / es, 1.0, 100.0).astype(np.float64)
