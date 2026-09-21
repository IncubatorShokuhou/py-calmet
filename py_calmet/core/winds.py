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
    x0 = float(getattr(threed, "x0_km", x_km - d