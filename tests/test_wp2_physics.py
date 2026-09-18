"""WP2: DIAGNO OA / IKINE / IOBR / clouds / precip / COARE-lite unit tests."""
from __future__ import annotations

import numpy as np
import pytest

from py_calmet.core import clouds, precip, overwater, winds
from py_calmet.core.similt import similt_profile


def test_cloud3_teixeira_limits():
    rh = np.array([0.0, 50.0, 99.0, 100.0])
    cc = clouds.cloud3_from_rh(rh)
    assert cc[0] == pytest.approx(0.18, abs=0.02)  # Teixera floor ~0.18 at RH=0
    assert 0.0 < cc[1] < 1.0
    assert cc[2] == 1.0
    assert cc[3] == 1.0


def test_cloud4_layered_rh():
    ny, nx, nk = 3, 3, 5
    # pressures: surface→aloft
    pres = np.zeros((nk, ny, nx))
    rh = np.zeros((nk, ny, nx))
    for k, p in enumerate([1000.0, 900.0, 700.0, 500.0, 300.0]):
        pres[k] = p
        rh[k] = 40.0
    rh[1] = 95.0  # low cloud layer
    ccp, ceil = clouds.cloud4_from_rh_profile(rh, pres)
    assert ccp.shape == (ny, nx)
    assert float(ccp.mean()) > 0.5
    assert float(ceil.mean()) > 0.0


def test_resolve_mcloud3_uses_rh():
    rh2 = np.full((4, 4), 90.0)
    cc = clouds.resolve_cloud_fraction(mcloud=3, rh_pct_2d=rh2)
    assert cc.shape == (4, 4)
    assert float(cc.mean()) > 0.4


def test_barnes_precip_and_prognostic():
    rmm = precip.barnes_precip(
        stn_x_m=np.array([500.0]),
        stn_y_m=np.array([500.0]),
        stn_rmm=np.array([2.5]),
        nx=3,
        ny=3,
        xorig_m=0.0,
        yorig_m=0.0,
        dgrid_m=500.0,
        sigmap_km=1.0,
        cutp=0.01,
    )
    assert rmm.shape == (3, 3)
    assert rmm.max() > 1.0
    z = precip.resolve_precip(npsta=0, nx=2, ny=2)
    assert np.all(z == 0.0)


def test_coare_lite_ustar_positive():
    u = np.full((5, 5), 5.0)
    v = np.zeros((5, 5))
    t_air = np.full((5, 5), 290.0)
    t_sea = np.full((5, 5), 292.0)
    ust, qh, el = overwater.coare_lite_fluxes(u, v, t_air, t_sea)
    assert np.all(ust > 0.01)
    assert np.all(qh > 0.0)  # sea warmer → upward buoyancy
    zi = overwater.mixht_overwater(ust, el, fcori=1e-4)
    assert np.all(zi >= 50.0)


def test_apply_overwater_only_on_water_lu():
    lu = np.array([[20, 55], [20, 55]])
    u = np.ones((2, 2)) * 4.0
    v = np.zeros((2, 2))
    t = np.full((2, 2), 288.0)
    rho = np.full((2, 2), 1.2)
    ust0 = np.full((2, 2), 0.2)
    el0 = np.full((2, 2), -100.0)
    qh0 = np.full((2, 2), 10.0)
    zi0 = np.full((2, 2), 500.0)
    ust, el, qh, zi = overwater.apply_overwater_pbl(
        lu, 55, 55, u, v, t, rho, ust0, el0, qh0, zi0, 1e-4, icoare=10
    )
    assert ust[0, 0] == pytest.approx(0.2)  # land unchanged
    assert ust[0, 1] != pytest.approx(0.2)  # water replaced


def test_topof2_kinematic_w_nonzero_on_slope():
    nz, ny, nx = 4, 5, 5
    zface = np.array([0.0, 20.0, 40.0, 80.0, 160.0])
    U = np.ones((nz, ny, nx))
    V = np.zeros_like(U)
    elev = np.tile(np.linspace(0, 200, nx), (ny, 1))
    tempk = np.full((ny, nx), 280.0)
    W = winds.topographic_kinematic_w(U, V, elev, zface, tempk, dgrid_m=1000.0, alpha=0.1)
    assert W.shape == (nz, ny, nx)
    assert np.max(np.abs(W)) > 0.0


def test_obrien_reduces_divergence():
    nz, ny, nx = 3, 6, 6
    zface = np.array([0.0, 50.0, 100.0, 200.0])
    rng = np.random.default_rng(0)
    U = rng.normal(size=(nz, ny, nx))
    V = rng.normal(size=(nz, ny, nx))
    dgrid = 1000.0
    div0 = np.gradient(U[0], dgrid, axis=1) + np.gradient(V[0], dgrid, axis=0)
    U2, V2 = winds.divergence_minimize(U, V, dgrid, niter=20, divlim=1e-8)
    div1 = np.gradient(U2[0], dgrid, axis=1) + np.gradient(V2[0], dgrid, axis=0)
    assert float(np.max(np.abs(div1))) <= float(np.max(np.abs(div0))) + 1e-9


def test_iextrp_similt_positive4():
    zface = np.array([0.0, 20.0, 40.0, 80.0, 160.0])
    # Fake sounding levels
    class Lev:
        def __init__(self, h, wd, ws):
            self.height = h
            self.wd = wd
            self.ws = ws
    levels = [Lev(100, 270, 5), Lev(500, 280, 8), Lev(1500, 290, 10)]
    U4, V4 = winds.obs_profile_similt(
        3.0, 0.0, 10.0, 0.1, -50.0, 200.0, zface, levels, 0.0, 50.0, 2, 2, iextrp=4
    )
    Um, Vm = winds.obs_profile_similt(
        3.0, 0.0, 10.0, 0.1, -50.0, 200.0, zface, levels, 0.0, 50.0, 2, 2, iextrp=-4
    )
    # +4 (SIMILT) differs from golden -4 power-law path
    assert not np.allclose(U4, Um)


def test_multi_station_oa_with_rprog():
    nz, ny, nx = 2, 5, 5
    Ug = np.ones((nz, ny, nx))
    Vg = np.zeros_like(Ug)
    u_stn = np.array([[5.0, 6.0], [0.0, 0.0]])  # (nstn=2, nz) — wait shape (nstn, nz)
    u_stn = np.array([[5.0, 6.0], [4.0, 5.0]])
    v_stn = np.zeros((2, 2))
    xs = np.array([500.0, 4500.0])
    ys = np.array([2500.0, 2500.0])
    U, V = winds.objective_analyze(
        Ug, Vg, u_stn, v_stn,
        xs_m=xs, ys_m=ys,
        xorig_m=0.0, yorig_m=0.0, dgrid_m=1000.0,
        r1_m=2000.0, r2_m=3000.0, rprog_m=5000.0,
        rmax1_m=10000.0, rmax2_m=10000.0, nintr2=[2, 2],
    )
    assert U.shape == (nz, ny, nx)
    # Near first station, U closer to 5 than to IGF=1
    assert U[0, 2, 0] > 2.0


def test_similt_profile_smoke():
    zmid = np.array([10.0, 30.0, 60.0, 120.0])
    u, v = similt_profile(3.0, 1.0, 10.0, 0.1, -100.0, 400.0, zmid)
    assert u.shape == zmid.shape
    assert np.isfinite(u[0])


def test_ikine_iobr_gated_run_smoke(tmp_path):
    """IKINE=1 / IOBR=1 alter W and stay finite on a tiny synthetic domain."""
    from py_calmet.core.winds import topographic_kinematic_w, obrien_adjust, vertical_velocity_from_div

    nz, ny, nx = 3, 4, 4
    zface = np.array([0.0, 20.0, 50.0, 100.0])
    U = np.ones((nz, ny, nx)) * 2.0
    V = np.ones((nz, ny, nx)) * 0.5
    elev = np.linspace(0, 100, ny * nx).reshape(ny, nx)
    tempk = np.full((ny, nx), 285.0)
    Wk = topographic_kinematic_w(U, V, elev, zface, tempk, 500.0, alpha=0.1)
    W0 = vertical_velocity_from_div(U, V, zface, 500.0)
    U2, V2, W2 = obrien_adjust(U, V, W0 + Wk, zface, 500.0, niter=10, divlim=1e-5)
    assert np.all(np.isfinite(U2)) and np.all(np.isfinite(W2))
    assert Wk.shape == W2.shape


def test_npsta_minus1_uses_rain_field():
    class T:
        ni = 3
        nj = 3
        x0_km = 0.0
        y0_km = 0.0
        dx_km = 1.0

    rain = np.ones((3, 3)) * 1.5
    rmm = precip.resolve_precip(
        npsta=-1, nx=2, ny=2, rain_prog=rain, threed=T(),
        xorig_km=0.0, yorig_km=0.0, dgrid_km=1.0,
    )
    assert rmm.shape == (2, 2)
    assert float(rmm.mean()) == pytest.approx(1.5)
