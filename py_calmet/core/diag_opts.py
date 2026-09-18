"""Diagnostic wind-module helpers: CGAMMA (ZUPT), domain-avg UA wind (IUPWND/ZUPWND).

Deepens IDIOPT1–5 / ISURFT / IUPT / ZUPT / IUPWND / ZUPWND beyond METLST echo:
when IDIOPTn=0 (default), compute internally from observations / 3D; when =1,
preprocessed diagnostic-file input is not supported (QA note only).
"""
from __future__ import annotations

from typing import Any

import numpy as np

from .met_utils import wind_uv


def cgamma_from_sounding(
    zl: np.ndarray,
    tz: np.ndarray,
    *,
    zupt: float = 200.0,
    ziconv: float = 0.0,
    ts_sfc: float | None = None,
    miss: float = 999.0,
) -> float:
    """CALMET CGAMMA-style average T lapse (K/m) in the layer up to ``ZUPT``.

    If ``ziconv >= zupt``, return dry-adiabatic ``-0.0098``. Otherwise
    ``(T(zupt) - Ts) / zupt`` with optional surface temperature override.
    """
    zupt = max(float(zupt), 1.0)
    if float(ziconv) >= zupt:
        return -0.0098
    zl = np.asarray(zl, dtype=np.float64)
    tz = np.asarray(tz, dtype=np.float64)
    valid = np.isfinite(zl) & np.isfinite(tz) & (tz < miss - 0.01)
    if int(valid.sum()) < 2:
        return 0.0
    zv, tv = zl[valid], tz[valid]
    order = np.argsort(zv)
    zv, tv = zv[order], tv[order]
    ts = float(ts_sfc) if ts_sfc is not None else float(tv[0])
    tt = float(np.interp(zupt, zv, tv))
    return (tt - ts) / zupt


def domain_avg_wind_from_sounding(
    levels,
    *,
    zlo: float = 1.0,
    zhi: float = 1000.0,
    stn_elev: float = 0.0,
) -> tuple[float, float]:
    """Layer-average U/V from one UP sounding between ``ZUPWND`` heights AGL.

    Mirrors CALMET ``VERTAV`` used when ``IDIOPT3=0`` and ``IUPWND`` selects a
    station (or -1 → caller averages several soundings).
    """
    zlo = float(zlo)
    zhi = max(float(zhi), zlo + 1.0)
    zs: list[float] = []
    us: list[float] = []
    vs: list[float] = []
    for lev in levels:
        hag = float(lev.height) - float(stn_elev)
        if hag < zlo or hag > zhi:
            continue
        if not np.isfinite(lev.ws) or not np.isfinite(lev.wd):
            continue
        if float(lev.ws) >= 998.0 or float(lev.wd) >= 998.0:
            continue
        u, v = wind_uv(float(lev.wd), float(lev.ws))
        zs.append(hag)
        us.append(float(u))
        vs.append(float(v))
    if not us:
        return 0.0, 0.0
    return float(np.mean(us)), float(np.mean(vs))


def resolve_diag_gamma(
    *,
    idiopt2: int,
    zupt: float,
    ziconv_mean: float,
    sounding_z: np.ndarray | None,
    sounding_t: np.ndarray | None,
    temp_sfc: float | None = None,
    daytime: bool = True,
) -> tuple[float, str]:
    """Pick Froude/TOPOF2 lapse: CGAMMA when IDIOPT2=0 + sounding, else proxy.

    Returns ``(gamma_K_per_m, note)``.
    """
    if int(idiopt2) != 0:
        g = 0.01 if not daytime else 0.005
        return g, "IDIOPT2=1: preprocessed lapse not supported; using day/night proxy"
    if sounding_z is not None and sounding_t is not None:
        g = cgamma_from_sounding(
            sounding_z, sounding_t, zupt=zupt, ziconv=ziconv_mean, ts_sfc=temp_sfc
        )
        # Froude / kinematic W expect a stable (positive) pot-temp proxy when
        # blocking; convert environmental T-lapse toward pot-temp-ish magnitude.
        g_use = abs(g + 0.0098) if g < 0 else max(g, 1e-4)
        if not daytime:
            g_use = max(g_use, 0.005)
        return float(g_use), f"CGAMMA(ZUPT={zupt:g}) → {g:.5f} K/m (diag {g_use:.5f})"
    g = 0.01 if not daytime else 0.005
    return g, "no sounding for CGAMMA; day/night proxy"


def qa_idiopt(idiopts: list[int], *, irtype: int = 1) -> list[str]:
    """QA notes for IDIOPT1–5 (preprocessed paths are soft spots)."""
    notes: list[str] = []
    names = ("IDIOPT1", "IDIOPT2", "IDIOPT3", "IDIOPT4", "IDIOPT5")
    for name, val in zip(names, idiopts):
        if int(val) == 0:
            continue
        if name in ("IDIOPT4", "IDIOPT5") and int(irtype) != 0:
            notes.append(
                f"{name}=1 with IRTYPE={irtype}: Fortran forbids preprocessed "
                "winds when computing full met fields; ignored"
            )
        else:
            notes.append(
                f"{name}=1: preprocessed diagnostic-file input not supported "
                "(computed internally / ignored)"
            )
    return notes


def diag_meta_dict(
    *,
    idiopt1: int,
    idiopt2: int,
    idiopt3: int,
    idiopt4: int,
    idiopt5: int,
    zupt: float,
    iupwnd: int,
    zupwnd: list[float],
    um: float = 0.0,
    vm: float = 0.0,
    gamma_diag: float | None = None,
) -> dict[str, Any]:
    return {
        "idiopt1": int(idiopt1),
        "idiopt2": int(idiopt2),
        "idiopt3": int(idiopt3),
        "idiopt4": int(idiopt4),
        "idiopt5": int(idiopt5),
        "zupt": float(zupt),
        "iupwnd": int(iupwnd),
        "zupwnd": [float(x) for x in zupwnd],
        "um_domain": float(um),
        "vm_domain": float(vm),
        "gamma_diag": None if gamma_diag is None else float(gamma_diag),
    }
