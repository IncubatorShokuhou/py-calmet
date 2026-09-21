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
                        np.log(zs[0]) - np.log(ZO_EXTRA