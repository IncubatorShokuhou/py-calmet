"""Daytime convective ZI golden (Maul–Carson path vs Fortran)."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pytest

from py_calmet import run_calmet, read_calmet_dat
from py_calmet.core.met_utils import relative_rmse
from thresholds import DAYTIME_ZI_THRESH

ROOT = Path(__file__).resolve().parents[1]
GOLDENS = ROOT / "cases" / "daytime_zi" / "goldens"
INPUTS = GOLDENS / "inputs"


def _corr(a, b) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    if a.std() < 1e-9 or b.std() < 1e-9:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


@pytest.mark.skipif(not (GOLDENS / "noobs" / "CALMET.DAT").exists(), reason="daytime golden missing")
def test_daytime_zi_carson_path():
    res = run_calmet(GOLDENS / "noobs", mode="noobs", inputs_dir=INPUTS)
    gold = read_calmet_dat(GOLDENS / "noobs" / "CALMET.DAT")
    Zi = gold.get_2d_field("ZI")
    Qsw = gold.get_2d_field("QSW")
    thr = DAYTIME_ZI_THRESH

    day = Qsw.mean(axis=(1, 2)) > 50.0
    night = Qsw.mean(axis=(1, 2)) < 1.0
    assert day.any() and night.any()

    # Solar short-wave tightly matched (Holtslag constants)
    assert relative_rmse(res.QSW[day], Qsw[day]) <= thr["QSW"]

    # Night mechanical ZI still close
    assert relative_rmse(res.ZI[night], Zi[night]) <= thr["ZI_night"]

    # Daytime convective path engaged
    assert (res.EL[day] < 0).mean() > 0.5
    assert res.WSTAR[day].mean() > 0.1
    assert res.ZI[day].mean() > res.ZI[night].mean()
    # Growing through the afternoon
    assert res.ZI[day][-1].mean() >= res.ZI[day][0].mean()

    # Shape vs Fortran (absolute rate differs: constant DPTMIN gamma)
    assert _corr(res.ZI[day], Zi[day]) >= thr["ZI_corr_day"]

    Ug, Vg = gold.get_3d_field("U"), gold.get_3d_field("V")
    assert relative_rmse(res.U, Ug) <= thr["U"]
    assert relative_rmse(res.V, Vg) <= thr["V"]


def test_mixht_day_grows_unit():
    from py_calmet.core import pbl

    qh = np.full((4, 4), 120.0)
    rho = np.full((4, 4), 1.2)
    tempk = np.full((4, 4), 300.0)
    ustar = np.full((4, 4), 0.35)
    zi0, zc0 = pbl.mixht_day_carson(qh, rho, tempk, ustar, 1e-4, dt_sec=3600.0, dtheta=0.001)
    zi1, zc1 = pbl.mixht_day_carson(
        qh, rho, tempk, ustar, 1e-4, dt_sec=3600.0, ziconv_prev=zc0, dtheta=0.001
    )
    assert zc1.mean() > zc0.mean()
    assert zi1.mean() >= 50.0
