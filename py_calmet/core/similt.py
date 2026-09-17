"""Van Ulden & Holtslag (1985) surface-wind extrapolation (CALMET SIMILT)."""
from __future__ import annotations
import numpy as np


def psiud(zL: float) -> float:
    """COARE/Dyer PSI_m with factor 16 (CALMET PSIUD)."""
    zL = float(zL)
    if zL < 0:
        x = (1.0 - 16.0 * zL) ** 0.25
        psik = (
            2.0 * np.log((1.0 + x) / 2.0)
            + np.log((1.0 + x * x) / 2.0)
            - 2.0 * np.arctan(x)
            + 2.0 * np.arctan(1.0)
        )
        y = (1.0 - 10.15 * zL) ** (1.0 / 3.0)
        psic = (
            1.5 * np.log((1.0 + y + y * y) / 3.0)
            - np.sqrt(3.0) * np.arctan((1.0 + 2.0 * y) / np.sqrt(3.0))
            + 4.0 * np.arctan(1.0) / np.sqrt(3.0)
        )
        f = zL * zL / (1.0 + zL * zL)
        return float((1.0 - f) * psik + f * psic)
    c = min(50.0, 0.35 * zL)
    return float(-((1.0 + zL) + 0.6667 * (zL - 14.28) / np.exp(c) + 8.525))


# Table 2 from Van Ulden & Holtslag (1985): 1/L bins and turning angle at 200 m
_ELTBLI = np.array(
    [-0.03333333, -0.01, -0.0027027, 0.0, 0.002857, 0.007692, 0.0166667, 0.05, 0.1111111]
)
_DHTBL2 = np.array([12.0, 10.0, 9.0, 12.0, 18.0, 28.0, 35.0, 38.0, 39.0])
_D1, _D2, _ZTBL2 = 1.58, 1.0, 200.0


def _dh200(elinv: float) -> float:
    if abs(elinv) < 5e-5:
        return float(_DHTBL2[3])
    for i in range(9):
        if _ELTBLI[i] >= elinv:
            if i == 0:
                return float(_DHTBL2[0])
            return float(
                _DHTBL2[i - 1]
                + (_DHTBL2[i] - _DHTBL2[i - 1])
                * (elinv - _ELTBLI[i - 1])
                / (_ELTBLI[i] - _ELTBLI[i - 1])
            )
    return float(_DHTBL2[8])


def similt_profile(
    u_anem: float,
    v_anem: float,
    z_anem: float,
    z0: float,
    el: float,
    zi: float,
    zmid: np.ndarray,
    zimin: float = 50.0,
    northern: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Extrapolate anemometer wind through CALMET layer midpoints (SIMILT).

    Levels at or below max(zi, zimin) use similarity; above are left as
    NaN for the caller to fill from upper-air / prognostic data.
    """
    nz = len(zmid)
    U = np.full(nz, np.nan)
    V = np.full(nz, np.nan)

    eladj = el
    if abs(eladj) < 5.0 * z0:
        eladj = -5.0 * z0 if eladj < 0 else 5.0 * z0
    elinv = 0.0 if abs(eladj) >= 10000.0 else 1.0 / eladj

    ws1 = float(np.hypot(u_anem, v_anem))
    # Fortran: angle = 270 - atan2(v,u)*factor; wd = amod(angle,360)
    angle = 270.0 - np.rad2deg(np.arctan2(v_anem, u_anem))
    wd1 = angle % 360.0
    if wd1 == 0.0:
        wd1 = 360.0

    xlnz1z0 = np.log(z_anem / z0)
    psim1 = 0.0 if abs(elinv) < 5e-5 else psiud(z_anem / eladj)

    dh200 = _dh200(elinv)
    dzan = dh200 * _D1 * (1.0 - np.exp(-_D2 * z_anem / _ZTBL2))

    z_limit = max(zi + 0.001, zimin)
    for k, zm in enumerate(zmid):
        # SIMILT uses zface(k) <= z_limit for level k-1; approximate with zm
        if zm > z_limit:
            break
        xlnzz0 = np.log(zm / z0)
        psim = 0.0 if abs(elinv) < 5e-5 else psiud(zm / eladj)
        ws = ws1 * (xlnzz0 - psim) / (xlnz1z0 - psim1)
        dwd = dh200 * _D1 * (1.0 - np.exp(-_D2 * zm / _ZTBL2))
        if northern:
            wd = wd1 + dwd - dzan
            if wd > 360.0:
                wd -= 360.0
        else:
            wd = wd1 - dwd + dzan
            if wd <= 0.0:
                wd = 360.0 + wd
        wdrad = np.deg2rad(wd)
        U[k] = -ws * np.sin(wdrad)
        V[k] = -ws * np.cos(wdrad)
    return U, V
