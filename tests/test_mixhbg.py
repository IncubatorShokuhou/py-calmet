"""Batchvarova–Gryning (IMIXH=±2) mixing-height unit + smoke tests."""
from __future__ import annotations

import numpy as np
import pytest

from py_calmet.config import CalmetConfig
from py_calmet.core import pbl, diag_opts


def test_imixh2_no_longer_raises():
    CalmetConfig(imixh=2).check_unsupported()
    CalmetConfig(imixh=-2).check_unsupported()


def test_fbg_positive_args():
    # FBG(zeta, zetap, ...) with positive args must be finite
    v = pbl._fbg(100.0, 80.0, -10.0, 5.0, 2.0)
    assert np.isfinite(v)


def test_mixhbg_grows_under_strong_flux():
    qh = np.full((4, 4), 150.0)
    rho = np.full((4, 4), 1.2)
    tempk = np.full((4, 4), 300.0)
    ustar = np.full((4, 4), 0.35)
    el = np.full((4, 4), -50.0)  # unstable
    zi0, zc0, dp0 = pbl.mixht_day_bg(
        qh, rho, tempk, ustar, el, 1e-4, dt_sec=3600.0, dtheta=0.005, threshl=0.05
    )
    zi1, zc1, _ = pbl.mixht_day_bg(
        qh,
        rho,
        tempk,
        ustar,
        el,
        1e-4,
        dt_sec=3600.0,
        ziconv_prev=zc0,
        dtheta=0.005,
        threshl=0.05,
    )
    assert zc0.mean() >= 50.0
    assert zc1.mean() >= zc0.mean() - 1e-6
    assert zi1.mean() >= zi0.mean() - 1e-6
    assert np.all(dp0 == 0.0)
    assert np.all(zi1 <= 3000.0)


def test_mixhbg_stable_collapses():
    qh = np.full((3, 3), 20.0)
    rho = np.full((3, 3), 1.2)
    tempk = np.full((3, 3), 290.0)
    ustar = np.full((3, 3), 0.2)
    el = np.full((3, 3), 100.0)  # stable → no convective growth
    _, zc, _ = pbl.mixht_day_bg(
        qh, rho, tempk, ustar, el, 1e-4, dt_sec=3600.0, dtheta=0.01
    )
    assert float(zc.mean()) == pytest.approx(0.0)


def test_mixhbg_weak_flux_relaxes():
    # Positive but sub-threshold buoyancy → equilibrium relaxation
    qh = np.full((2, 2), 5.0)  # small
    rho = np.full((2, 2), 1.2)
    tempk = np.full((2, 2), 295.0)
    ustar = np.full((2, 2), 0.25)
    el = np.full((2, 2), -80.0)
    prev = np.full((2, 2), 400.0)
    _, zc, _ = pbl.mixht_day_bg(
        qh,
        rho,
        tempk,
        ustar,
        el,
        1e-4,
        dt_sec=3600.0,
        ziconv_prev=prev,
        dtheta=0.01,
        threshl=0.05,  # wto = 0.05*400/(1.2*1005) ≈ 0.0166; wt=5/(1.2*1005)≈0.004
        izicrlx=1,
        tzicrlx=800.0,
    )
    # Should move toward small equilibrium, not stay at 400
    assert float(zc.mean()) < 400.0
    assert float(zc.mean()) > 0.0


def test_cgamma_adiabatic_when_zi_above_zupt():
    z = np.array([0.0, 100.0, 300.0, 600.0])
    t = np.array([290.0, 289.0, 287.0, 285.0])
    g = diag_opts.cgamma_from_sounding(z, t, zupt=200.0, ziconv=250.0)
    assert g == pytest.approx(-0.0098)


def test_cgamma_from_sounding_finite():
    z = np.array([0.0, 100.0, 200.0, 500.0])
    t = np.array([290.0, 289.0, 288.0, 285.0])
    g = diag_opts.cgamma_from_sounding(z, t, zupt=200.0, ziconv=50.0)
    assert np.isfinite(g)
    # (288-290)/200 = -0.01
    assert g == pytest.approx(-0.01, abs=1e-4)


def test_domain_avg_wind():
    class _L:
        def __init__(self, h, wd, ws):
            self.height = h
            self.wd = wd
            self.ws = ws

    levels = [_L(100.0, 270.0, 10.0), _L(500.0, 270.0, 10.0), _L(2000.0, 90.0, 20.0)]
    um, vm = diag_opts.domain_avg_wind_from_sounding(levels, zlo=1.0, zhi=1000.0)
    # WD=270 → u=+ws, v≈0
    assert um == pytest.approx(10.0, abs=0.1)
    assert abs(vm) < 0.5


def test_qa_idiopt4_with_irtype():
    notes = diag_opts.qa_idiopt([0, 0, 0, 1, 0], irtype=1)
    assert any("IDIOPT4" in n for n in notes)
