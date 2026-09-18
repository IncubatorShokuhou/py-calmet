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
