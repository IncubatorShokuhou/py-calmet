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
    """Map 3D.DAT winds onto CALMET 