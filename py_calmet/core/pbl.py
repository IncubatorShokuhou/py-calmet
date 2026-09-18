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
    # Zilitinkevich only for stable (el > 0); else mechanical-only
    el_pos = np.maximum(el, 0.0)
    zi2 = np.where(el_pos > 0.0, 0.4 * np.sqrt(ustar * el_pos / fcori), zi1)
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
    dtheta: float = 0.001,
    threshl: float = 0.05,
    constb: float = 1.41,
    conste: float = 0.15,
    zimin: float = 50.0,
    zimax: float = 3000.0,
    dptt_prev: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Daytime mixing height: Maul–Carson convective (MIXHMC) + mechanical.

    Returns (zi, ziconv, dptt) so the inversion jump ``dptt`` can be carried
    into the next hour (CALMET MIXHMC). Callers that only need heights can
    ignore the third array.
    """
    qh = np.asarray(qh, dtype=np.float64)
    rho = np.maximum(np.asarray(rho, dtype=np.float64), 0.5)
    wt = qh / (rho * CP)  # <w'Theta'> K m/s
    htold = np.zeros_like(qh) if ziconv_prev is None else np.asarray(ziconv_prev, dtype=np.float64).copy()
    dptt = np.zeros_like(qh) if dptt_prev is None else np.asarray(dptt_prev, dtype=np.float64).copy()

    # pot-temp lapse above zi (K/m); scalar or 2-D from MIXDT/MIXDT2
    gamma = np.asarray(dtheta, dtype=np.float64)
    if gamma.ndim == 0:
        gamma = np.full(qh.shape, max(float(gamma), 1e-4))
    else:
        gamma = np.maximum(gamma, 1e-4)
    onedte = dt_sec * (1.0 + conste)
    twodte = 2.0 * dt_sec * conste

    # Threshold buoyancy flux (K m/s): wto = thresh * h / (rho*cp)
    # thresh is W/m^2/m → divide by rho*cp gives K/s / m * h = K m/s
    wto = threshl * htold / (rho * CP)

    ziconv = np.zeros_like(qh)
    net = wt - wto
    grow = net > 0.0
    # Weakly convective: relax toward equilibrium zi = rho*cp*wt/thresh
    weak = (wt > 0.0) & ~grow
    if np.any(weak) and threshl > 0:
        ziceq = rho * CP * wt / max(threshl, 1e-6)
        ziconv = np.where(weak, ziceq + (htold - ziceq) * np.exp(-dt_sec / 800.0), ziconv)

    if np.any(grow):
        dpttp1 = np.sqrt(np.maximum(gamma * twodte * net, 0.0))
        unsqrt = htold ** 2 + 2.0 * (net * onedte - dptt * htold) / gamma
        unsqrt = np.maximum(unsqrt, 0.0)
        zic = np.sqrt(unsqrt) + dpttp1 / gamma
        ziconv = np.where(grow, np.minimum(np.maximum(zic, 0.0), zimax), ziconv)
        dptt = np.where(grow, dpttp1, dptt)

    # Mechanical daytime (Venkatram): CMECH * ustar
    cmech = constb / np.sqrt(np.maximum(np.asarray(fcori, dtype=np.float64), 1e-5))
    hmech = np.minimum(cmech * ustar, zimax)

    zi = np.maximum(np.maximum(zimin, hmech), ziconv)
    zi = np.minimum(zi, zimax)
    return zi, ziconv, dptt


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


# Holtslag & van Ulden (1983) short-wave / energy-budget constants (CALMET defaults)
HA1 = 990.0
HA2 = -30.0
HB1 = -0.75
HB2 = 3.4
HC1 = 5.31e-13
HC2 = 60.0
HC3 = 0.12


def sine_solar_elevation(
    lat_deg: float | np.ndarray,
    lon_deg_east: float | np.ndarray,
    jday: int,
    hour_utc: float,
    ibtz: int = 0,
) -> np.ndarray:
    """Sine of solar elevation (CALMET SOLAR, half-hour centered).

    ``ibtz`` is the base time zone offset such that local = UTC - ibtz
    (CALMET: 5=EST …). For ABTZ=UTC+0000 use ibtz=0.
    lon is **east** longitude (CALMET ≥ 050328 convention).
    """
    lat = np.asarray(lat_deg, dtype=np.float64)
    lon = np.asarray(lon_deg_east, dtype=np.float64)
    d = (float(jday) - 1.0) * 0.9856479
    radd = np.deg2rad(d)
    xsind, xcosd = np.sin(radd), np.cos(radd)
    rad2d = 2.0 * radd
    sin2d, cos2d = np.sin(rad2d), np.cos(rad2d)
    em = (
        12.0
        + 0.12357 * xsind
        - 0.004289 * xcosd
        + 0.153809 * sin2d
        + 0.060783 * cos2d
    )
    sigma = (
        279.9348
        + d
        + 1.914827 * xsind
        - 0.079525 * xcosd
        + 0.019938 * sin2d
        - 0.00162 * cos2d
    )
    sincd = 0.39784989 * np.sin(np.deg2rad(sigma))
    capd = np.arcsin(sincd)
    coscd = np.cos(capd)
    radlat = np.deg2rad(lat)
    sinlat, coslat = np.sin(radlat), np.cos(radlat)
    # half-hour after clock hour, in GMT = hour+0.5 + ibtz  (ihr-1 - 0.5 + ibtz with ihr=hour+1)
    gmt = float(hour_utc) + 0.5 + float(ibtz)
    solha = 15.0 * (gmt - em) + lon
    return sinlat * sincd + coslat * coscd * np.cos(np.deg2rad(solha))


def shortwave_radiation(
    sinalp: np.ndarray,
    ccfrac: float | np.ndarray = 0.0,
    *,
    ha1: float = HA1,
    ha2: float = HA2,
    hb1: float = HB1,
    hb2: float = HB2,
) -> np.ndarray:
    """QSW (W/m^2) from sine solar elevation and cloud fraction."""
    cc = np.asarray(ccfrac, dtype=np.float64)
    qsw = (ha1 * np.asarray(sinalp, dtype=np.float64) + ha2) * (1.0 + hb1 * cc ** hb2)
    return np.maximum(qsw, 0.0)


def heat_flux_energy_budget(
    qsw: np.ndarray,
    tempk: np.ndarray,
    sinalp: np.ndarray,
    ccfrac: float | np.ndarray = 0.0,
    albedo: float | np.ndarray = 0.2,
    bowen: float | np.ndarray = 1.0,
    hcg: float | np.ndarray = 0.15,
    qf: float | np.ndarray = 0.0,
    landuse: np.ndarray | None = None,
    iwat1: int = 55,
    iwat2: int = 55,
    *,
    hc1: float = HC1,
    hc2: float = HC2,
    hc3: float = HC3,
) -> np.ndarray:
    """Daytime sensible heat flux (Holtslag–van Ulden); night → -0.1 over land."""
    qsw = np.asarray(qsw, dtype=np.float64)
    tempk = np.asarray(tempk, dtype=np.float64)
    sinalp = np.asarray(sinalp, dtype=np.float64)
    cc = np.asarray(ccfrac, dtype=np.float64)
    alb = np.asarray(albedo, dtype=np.float64)
    bo = np.asarray(bowen, dtype=np.float64)
    hg = np.asarray(hcg, dtype=np.float64)
    qa = np.asarray(qf, dtype=np.float64)
    qh = np.full(tempk.shape, -0.1, dtype=np.float64)
    day = sinalp > 0.0
    if landuse is not None:
        water = (landuse >= iwat1) & (landuse <= iwat2)
        day = day & ~water
        qh = np.where(water, 0.0, qh)
    if np.any(day):
        qstar = (
            (1.0 - alb) * qsw
            + hc1 * tempk ** 6
            - 5.67e-8 * tempk ** 4
            + hc2 * cc
        ) / (hc3 + 1.0)
        qh_day = bo * (qstar * (1.0 - hg) + qa) / (1.0 + bo)
        qh = np.where(day, qh_day, qh)
    return qh
