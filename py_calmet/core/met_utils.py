"""Shared meteorological helpers."""
from __future__ import annotations
import numpy as np

VK = 0.4
G = 9.81
CP = 996.0
ZO_EXTRAP = 0.5  # RDMM5 log-profile roughness for levels below first MM5 level


def wind_uv(wd_deg: np.ndarray | float, ws: np.ndarray | float):
    """Meteorological wind direction (from) + speed -> U, V (to the east/north)."""
    rad = np.deg2rad(wd_deg)
    u = -ws * np.sin(rad)
    v = -ws * np.cos(rad)
    return u, v


def coriolis(lat_deg: float) -> float:
    return max(abs(2.0 * 7.2921e-5 * np.sin(np.deg2rad(lat_deg))), 0.25e-4)


def layer_mids(zface: np.ndarray) -> np.ndarray:
    return 0.5 * (zface[:-1] + zface[1:])


def utm_to_latlon(
    easting_m: float,
    northing_m: float,
    zone: int,
    northern: bool = True,
) -> tuple[float, float]:
    """WGS84 UTM → (latitude_deg, longitude_deg_east).

    Pure-NumPy inverse so solar / Coriolis do not depend on pyproj, and so a
    missing optional extra cannot silently fall back to a Maine default.
    """
    a = 6378137.0
    f = 1.0 / 298.257223563
    k0 = 0.9996
    e2 = f * (2.0 - f)
    ep2 = e2 / (1.0 - e2)
    x = float(easting_m) - 500000.0
    y = float(northing_m) if northern else float(northing_m) - 1.0e7
    m = y / k0
    mu = m / (a * (1.0 - e2 / 4.0 - 3.0 * e2**2 / 64.0 - 5.0 * e2**3 / 256.0))
    e1 = (1.0 - np.sqrt(1.0 - e2)) / (1.0 + np.sqrt(1.0 - e2))
    fp = (
        mu
        + (3.0 * e1 / 2.0 - 27.0 * e1**3 / 32.0) * np.sin(2.0 * mu)
        + (21.0 * e1**2 / 16.0 - 55.0 * e1**4 / 32.0) * np.sin(4.0 * mu)
        + (151.0 * e1**3 / 96.0) * np.sin(6.0 * mu)
        + (1097.0 * e1**4 / 512.0) * np.sin(8.0 * mu)
    )
    sinf = np.sin(fp)
    cosf = np.cos(fp)
    tanf = np.tan(fp)
    c1 = ep2 * cosf**2
    t1 = tanf**2
    n1 = a / np.sqrt(1.0 - e2 * sinf**2)
    r1 = a * (1.0 - e2) / (1.0 - e2 * sinf**2) ** 1.5
    d = x / (n1 * k0)
    lat = fp - (n1 * tanf / r1) * (
        d**2 / 2.0
        - (5.0 + 3.0 * t1 + 10.0 * c1 - 4.0 * c1**2 - 9.0 * ep2) * d**4 / 24.0
        + (61.0 + 90.0 * t1 + 298.0 * c1 + 45.0 * t1**2 - 252.0 * ep2 - 3.0 * c1**2)
        * d**6
        / 720.0
    )
    lon = (
        d
        - (1.0 + 2.0 * t1 + c1) * d**3 / 6.0
        + (5.0 - 2.0 * c1 + 28.0 * t1 - 3.0 * c1**2 + 8.0 * ep2 + 24.0 * t1**2)
        * d**5
        / 120.0
    ) / cosf
    lon0 = (int(zone) - 1) * 6 - 180 + 3
    return float(np.degrees(lat)), float(np.degrees(lon) + lon0)


def relative_rmse(
    pred: np.ndarray,
    ref: np.ndarray,
    eps: float = 1e-6,
    floor: float | None = None,
) -> float:
    """Relative RMSE. Optional ``floor`` clamps |ref| in the denominator
    (useful on coarse complex-terrain fields where |ref|≈0 inflates the metric).
    """
    denom = np.abs(ref) + eps
    if floor is not None:
        denom = np.maximum(denom, floor)
    return float(np.sqrt(np.mean(((pred - ref) / denom) ** 2)))


def max_rel_err(pred: np.ndarray, ref: np.ndarray, eps: float = 1e-6) -> float:
    return float(np.max(np.abs(pred - ref) / (np.abs(ref) + eps)))
