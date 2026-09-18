"""Parity / smoke tests for the wrf_demo case (NCAR Katrina wrfout)."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pytest

from py_calmet import run_calmet, read_calmet_dat
from py_calmet.core.met_utils import relative_rmse
from thresholds import WRF_DEMO_THRESH

ROOT = Path(__file__).resolve().parents[1]
GOLDENS = ROOT / "cases" / "wrf_demo" / "goldens"
INPUTS = GOLDENS / "inputs"


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    if a.std() < 1e-9 or b.std() < 1e-9:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


@pytest.mark.parametrize("mode", ["obs", "obs_model", "noobs"])
def test_wrf_demo_parity(mode):
    if not (GOLDENS / mode / "CALMET.DAT").exists():
        pytest.skip("wrf_demo golden missing")
    res = run_calmet(GOLDENS / mode, mode=mode, inputs_dir=INPUTS)
    gold = read_calmet_dat(GOLDENS / mode / "CALMET.DAT")
    Ug = gold.get_3d_field("U")
    Vg = gold.get_3d_field("V")
    thr = WRF_DEMO_THRESH[mode]
    stats = {
        "U": relative_rmse(res.U, Ug),
        "V": relative_rmse(res.V, Vg),
        "U_corr": _corr(res.U, Ug),
        "V_corr": _corr(res.V, Vg),
        "ZI": relative_rmse(res.ZI, gold.get_2d_field("ZI")),
        "USTAR": relative_rmse(res.USTAR, gold.get_2d_field("USTAR")),
        "SPD": relative_rmse(np.hypot(res.U, res.V), np.hypot(Ug, Vg)),
    }
    assert res.U.shape == Ug.shape
    for key in ("U", "V", "ZI", "USTAR", "SPD"):
        assert stats[key] <= thr[key], f"{mode} {key}={stats[key]:.4e} > {thr[key]}"
    assert stats["U_corr"] >= thr["U_corr"], f"{mode} U_corr={stats['U_corr']:.3f}"
    assert stats["V_corr"] >= thr["V_corr"], f"{mode} V_corr={stats['V_corr']:.3f}"
    # finite output fields
    for name in ("W", "T", "IPGT", "RHO", "QSW", "IRH", "WSTAR", "EL"):
        arr = getattr(res, name)
        assert np.isfinite(arr).all(), name
