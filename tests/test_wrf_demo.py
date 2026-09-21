"""Primary (and only) WRF/Fortran validation gate: real Katrina wrf_demo.

Goldens under ``cases/wrf_demo/goldens/`` are Fortran CALMET 6.5.0 outputs
from the NCAR Katrina tutorial wrfout mountain window with real terrain
``GEO.DAT``. Synthetic ``small_domain`` / ``daytime_zi`` cases were removed;
this module is the suite's sole Fortran-compare gate.

Reference preference: live ``calmet.x`` when vendored; otherwise archived
Fortran-from-WRF-case ``CALMET.DAT`` under goldens/ (still real-WRF, never
synthetic). Run ``scripts/compare_wrf_fortran.py`` for a CLI report.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from py_calmet import run_calmet, read_calmet_dat
from py_calmet.core.met_utils import relative_rmse
from thresholds import WRF_DEMO_THRESH, WRF_DEMO_UV_FLOOR

ROOT = Path(__file__).resolve().parents[1]
CASE = ROOT / "cases" / "wrf_demo"
GOLDENS = CASE / "goldens"
INPUTS = GOLDENS / "inputs"
SHARED = CASE / "shared"

pytestmark = pytest.mark.wrf_primary


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    if a.std() < 1e-9 or b.std() < 1e-9:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def _parse_geo_elev(path: Path) -> np.ndarray:
    lines = path.read_text().splitlines()
    idx = next(i for i, ln in enumerate(lines) if "HTFAC" in ln)
    rows = []
    for ln in lines[idx + 1 :]:
        parts = ln.split()
        if not parts:
            break
        try:
            rows.append([float(x) for x in parts])
        except ValueError:
            break
        if len(rows) >= 12:
            break
    arr = np.asarray(rows, dtype=np.float64)
    assert arr.shape == (12, 12), f"geo elev shape {arr.shape}"
    return arr


def test_wrf_demo_geo_dat_in_path():
    """GEO.DAT (real Katrina mountain terrain) must be on the validation path."""
    for path in (SHARED / "geo.dat", INPUTS / "geo.dat"):
        assert path.is_file(), f"missing {path}"
        text = path.read_text()
        assert "Katrina" in text or "wrfout" in text.lower()
        assert "14N" in text or "UTM" in text
        elev = _parse_geo_elev(path)
        assert float(elev.max()) > 1000.0, "expected mountain terrain in geo.dat"
        assert float(elev.min()) >= 0.0
    assert (SHARED / "geo.dat").read_bytes() == (INPUTS / "geo.dat").read_bytes()
    assert (SHARED / "3d.dat").read_bytes() == (INPUTS / "3d.dat").read_bytes()


def test_wrf_demo_3d_is_real_wrf_case():
    """3D.DAT header must identify the Katrina WRF hourly-interp case."""
    head = (SHARED / "3d.dat").read_text(errors="replace")[:800]
    assert "3D.DAT" in head
    assert "2005" in head
    assert "Katrina" in head or "wrfout" in head.lower() or "hourly" in head.lower()


@pytest.mark.parametrize("mode", ["noobs", "obs_model", "obs"])
def test_wrf_demo_parity(mode):
    gold_path = GOLDENS / mode / "CALMET.DAT"
    assert gold_path.is_file(), (
        f"wrf_demo Fortran golden missing for mode={mode}; "
        "primary gate cannot skip — restore cases/wrf_demo/goldens/"
    )
    assert (INPUTS / "geo.dat").is_file(), "geo.dat required on wrf_demo input path"
    res = run_calmet(GOLDENS / mode, mode=mode, inputs_dir=INPUTS)
    gold = read_calmet_dat(gold_path)
    Ug = gold.get_3d_field("U")
    Vg = gold.get_3d_field("V")
    thr = WRF_DEMO_THRESH[mode]
    fl = WRF_DEMO_UV_FLOOR
    stats = {
        "U": relative_rmse(res.U, Ug, floor=fl),
        "V": relative_rmse(res.V, Vg, floor=fl),
        "U_corr": _corr(res.U, Ug),
        "V_corr": _corr(res.V, Vg),
        "ZI": relative_rmse(res.ZI, gold.get_2d_field("ZI")),
        "USTAR": relative_rmse(res.USTAR, gold.get_2d_field("USTAR")),
        "SPD": relative_rmse(np.hypot(res.U, res.V), np.hypot(Ug, Vg), floor=fl),
    }
    assert res.U.shape == Ug.shape
    for key in ("U", "V", "ZI", "USTAR", "SPD"):
        assert stats[key] <= thr[key], (
            f"PRIMARY wrf_demo {mode} {key}={stats[key]:.4e} > {thr[key]} "
            f"(Fortran-from-WRF-case golden)"
        )
    assert stats["U_corr"] >= thr["U_corr"], (
        f"PRIMARY wrf_demo {mode} U_corr={stats['U_corr']:.3f} < {thr['U_corr']}"
    )
    assert stats["V_corr"] >= thr["V_corr"], (
        f"PRIMARY wrf_demo {mode} V_corr={stats['V_corr']:.3f} < {thr['V_corr']}"
    )
    for name in ("W", "T", "IPGT", "RHO", "QSW", "IRH", "WSTAR", "EL"):
        arr = getattr(res, name)
        assert np.isfinite(arr).all(), name
