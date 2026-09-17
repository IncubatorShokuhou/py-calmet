"""Parity tests: pure-NumPy py-calmet vs Fortran CALMET goldens."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pytest

from py_calmet import run_calmet, read_calmet_dat
from py_calmet.core.met_utils import relative_rmse
from thresholds import THRESH

ROOT = Path(__file__).resolve().parents[1]
GOLDENS = ROOT / "cases" / "small_domain" / "goldens"
INPUTS = GOLDENS / "inputs"


def _run_and_compare(mode: str):
    case = GOLDENS / mode
    res = run_calmet(case, mode=mode, inputs_dir=INPUTS)
    gold = read_calmet_dat(case / "CALMET.DAT")
    Ug = gold.get_3d_field("U")
    Vg = gold.get_3d_field("V")
    zi = gold.get_2d_field("ZI")
    ust = gold.get_2d_field("USTAR")
    thr = THRESH[mode]
    return {
        "U": relative_rmse(res.U, Ug),
        "V": relative_rmse(res.V, Vg),
        "U_lev1": relative_rmse(res.U[:, 0], Ug[:, 0]),
        "V_lev1": relative_rmse(res.V[:, 0], Vg[:, 0]),
        "ZI": relative_rmse(res.ZI, zi),
        "USTAR": relative_rmse(res.USTAR, ust),
        "SPD": relative_rmse(np.hypot(res.U, res.V), np.hypot(Ug, Vg)),
        "thr": thr,
        "res": res,
        "gold_U": Ug,
    }


@pytest.mark.parametrize("mode", ["obs", "obs_model", "noobs"])
def test_core_parity(mode):
    stats = _run_and_compare(mode)
    thr = stats["thr"]
    for key in ("U", "V", "U_lev1", "V_lev1", "ZI", "USTAR", "SPD"):
        assert stats[key] <= thr[key], (
            f"{mode} {key} relative RMSE {stats[key]:.4e} exceeds {thr[key]}"
        )


@pytest.mark.parametrize("mode", ["obs", "obs_model", "noobs"])
def test_shapes(mode):
    stats = _run_and_compare(mode)
    assert stats["res"].U.shape == stats["gold_U"].shape
