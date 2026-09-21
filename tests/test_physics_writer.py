"""Unit tests for expanded physics helpers and CALMET.DAT writer round-trip."""
from __future__ import annotations
from datetime import datetime
from pathlib import Path
import numpy as np
import pytest

from py_calmet.core import winds, pbl
from py_calmet.core.runner import run_calmet, CalmetResult
from py_calmet.io.calmet_dat import (
    RunControl,
    read_calmet_dat,
    write_calmet_dat,
    write_calmet_netcdf,
)

ROOT = Path(__file__).resolve().parents[1]
TINY = ROOT / "cases" / "wrf_demo" / "goldens"  # real WRF case


def test_slope_flow_nonzero_on_slope():
    elev = np.linspace(0, 500, 12).reshape(3, 4)
    elev = np.tile(elev, (2, 1))[:4, :4]  # 4x4 ramp
    qh = np.full((4, 4), -20.0)
    tempk = np.full((4, 4), 290.0)
    rho = np.full((4, 4), 1.2)
    zface = np.array([0.0, 20.0, 40.0, 80.0])
    U, V = winds.slope_flow(elev, 1000.0, qh, tempk, rho, zface)
    assert U.shape[0] == 3
    assert np.hypot(U[0], V[0]).mean() > 0


def test_divergence_minimize_reduces_div():
    ny = nx = 8
    U = np.zeros((2, ny, nx))
    V = np.zeros_like(U)
    # divergent field
    x = np.linspace(-1, 1, nx)
    y = np.linspace(-1, 1, ny)
    xx, yy = np.meshgrid(x, y)
    U[0] = xx
    V[0] = yy
    dgrid = 1000.0
    div0 = np.gradient(U[0], dgrid, axis=1) + np.gradient(V[0], dgrid, axis=0)
    U2, V2 = winds.divergence_minimize(U, V, dgrid, niter=40)
    div1 = np.gradient(U2[0], dgrid, axis=1) + np.gradient(V2[0], dgrid, axis=0)
    assert np.abs(div1[2:-2, 2:-2]).mean() < np.abs(div0[2:-2, 2:-2]).mean()


def test_mixht_day_carson_grows():
    qh = np.full((5, 5), 100.0)
    rho = np.full((5, 5), 1.2)
    tempk = np.full((5, 5), 300.0)
    ustar = np.full((5, 5), 0.3)
    zi0, zc0, dp0 = pbl.mixht_day_carson(qh, rho, tempk, ustar, 1e-4, dt_sec=3600.0, dtheta=0.001)
    zi1, zc1, dp1 = pbl.mixht_day_carson(
        qh, rho, tempk, ustar, 1e-4, dt_sec=3600.0, ziconv_prev=zc0, dptt_prev=dp0, dtheta=0.001
    )
    assert zc1.mean() >= zc0.mean()
    assert zi1.mean() >= 50.0
    assert dp1.mean() > 0.0
    # Carrying the inversion jump forward must not explode growth vs a zero dptt restart.
    zi2, zc2, _ = pbl.mixht_day_carson(
        qh, rho, tempk, ustar, 1e-4, dt_sec=3600.0, ziconv_prev=zc0, dptt_prev=None, dtheta=0.001
    )
    assert zc1.mean() <= zc2.mean() + 1e-6


def test_ipgt_wstar_rho():
    el = np.array([[-500.0, 10.0, 300.0]])
    ip = pbl.ipgt_from_el(el)
    assert ip[0, 0] <= 3 and ip[0, 2] >= 5
    qh = np.array([[50.0, -10.0]])
    ws = pbl.wstar_field(np.array([[1000.0, 0.0]]), qh, np.array([[300.0, 300.0]]), np.array([[1.2, 1.2]]))
    assert ws[0, 0] > 0 and ws[0, 1] == 0
    rho = pbl.air_density(np.array([290.0]), 1013.0)
    assert 1.0 < float(np.asarray(rho).ravel()[0]) < 1.5


def test_calmet_dat_roundtrip(tmp_path):
    res = run_calmet(TINY / "obs", mode="obs", inputs_dir=TINY / "inputs")
    gold = read_calmet_dat(TINY / "obs" / "CALMET.DAT")
    rc = gold.run_control
    out = tmp_path / "out.dat"
    write_calmet_dat(
        out,
        result=res,
        run_control=rc,
        z0=res.z0,
        landuse=np.full(res.elev.shape, 20, dtype=np.int32),
        elev=res.elev,
        xssta=gold.static_fields.get("XSSTA"),
        yssta=gold.static_fields.get("YSSTA"),
        xusta=gold.static_fields.get("XUSTA"),
        yusta=gold.static_fields.get("YUSTA"),
        nears=gold.static_fields.get("NEARS"),
    )
    back = read_calmet_dat(out)
    assert back.nx == rc.nx and back.ny == rc.ny and back.nz == rc.nz
    assert back.nt == res.U.shape[0]
    assert np.allclose(back.get_3d_field("U"), res.U, rtol=1e-4, atol=1e-3)
    assert np.allclose(back.get_2d_field("ZI"), res.ZI, rtol=1e-4, atol=1e-2)


def test_netcdf_writer(tmp_path):
    res = run_calmet(TINY / "noobs", mode="noobs", inputs_dir=TINY / "inputs")
    gold = read_calmet_dat(TINY / "noobs" / "CALMET.DAT")
    out = tmp_path / "out.nc"
    write_calmet_netcdf(out, res, gold.run_control)
    assert out.exists() and out.stat().st_size > 1000


def test_utm_to_latlon_statue_of_liberty():
    from py_calmet.core.met_utils import utm_to_latlon

    # WGS84 UTM 18N for 40.6892 N, 74.0445 W
    lat, lon = utm_to_latlon(580735.87, 4504695.17, 18, northern=True)
    assert abs(lat - 40.6892) < 1e-4
    assert abs(lon - (-74.0445)) < 1e-4


def test_run_window_crosses_midnight():
    from py_calmet.core.runner import _run_window
    from py_calmet.io.inp import CalmetInp

    inp = CalmetInp(
        raw={
            "IBYR": "2005",
            "IBMO": "8",
            "IBDY": "28",
            "IBHR": "22",
            "IEYR": "2005",
            "IEMO": "8",
            "IEDY": "29",
            "IEHR": "2",
            "NSECDT": "3600",
        }
    )
    start, end, nhrs, nsecdt = _run_window(inp)
    assert nhrs == 4
    assert nsecdt == 3600
    assert start.hour == 22 and end.hour == 2


def test_write_calmet_dat_npsta(tmp_path):
    from dataclasses import replace

    res = run_calmet(TINY / "obs", mode="obs", inputs_dir=TINY / "inputs")
    gold = read_calmet_dat(TINY / "obs" / "CALMET.DAT")
    rc = replace(gold.run_control, npsta=1)
    out = tmp_path / "precip.dat"
    write_calmet_dat(
        out,
        result=res,
        run_control=rc,
        z0=res.z0,
        landuse=np.full(res.elev.shape, 20, dtype=np.int32),
        elev=res.elev,
        xssta=gold.static_fields.get("XSSTA"),
        yssta=gold.static_fields.get("YSSTA"),
        xusta=gold.static_fields.get("XUSTA"),
        yusta=gold.static_fields.get("YUSTA"),
        nears=gold.static_fields.get("NEARS"),
        xpsta=np.array([1.0], dtype=np.float32),
        ypsta=np.array([2.0], dtype=np.float32),
        rmm=np.zeros_like(res.ZI),
    )
    back = read_calmet_dat(out)
    assert "XPSTA" in back.static_fields
    assert "RMM" in back.fields_2d
    assert np.allclose(back.static_fields["XPSTA"], [1.0])


def test_obs_profile_direction_wraps_0_360():
    from types import SimpleNamespace

    zface = np.array([0.0, 20.0, 40.0, 80.0, 160.0])
    # Surface from 350°, sounding from 10° — linear deg interp would go the long way.
    levels = [
        SimpleNamespace(height=200.0, wd=10.0, ws=5.0),
        SimpleNamespace(height=500.0, wd=10.0, ws=6.0),
    ]
    from py_calmet.core.met_utils import wind_uv

    u_sfc, v_sfc = wind_uv(350.0, 3.0)
    U, V = winds.obs_profile_similt(
        u_sfc=float(u_sfc),
        v_sfc=float(v_sfc),
        z_anem=10.0,
        z0=0.25,
        el=-50.0,
        zi=200.0,
        zface=zface,
        sounding_levels=levels,
        stn_elev=0.0,
        zimin=50.0,
        nx=2,
        ny=2,
    )
    wd_mid = float(np.rad2deg(np.arctan2(-U[1, 0, 0], -V[1, 0, 0])) % 360.0)
    # Blended direction should stay near 350–10, not near 180.
    assert wd_mid > 300.0 or wd_mid < 40.0


def test_3d_interp_uses_origin_not_plus_one_halo():
    from types import SimpleNamespace

    # CALMET and 3D share origin/spacing (no halo): i=0 must map to 3D i=0.
    nx = ny = 3
    zface = np.array([0.0, 20.0, 40.0])
    wd = np.zeros((1, 3, 3, 2))
    ws = np.ones((1, 3, 3, 2))
    wd[0, 0, 0, :] = 90.0  # easterly at 3D (0,0) → U negative
    wd[0, 1, 1, :] = 180.0
    height = np.full((1, 3, 3, 2), 50.0)
    height[..., 1] = 150.0
    threed = SimpleNamespace(
        ni=3,
        nj=3,
        x0_km=0.0,
        y0_km=0.0,
        dx_km=1.0,
        elev=np.zeros((3, 3)),
        height_msl=height,
        wd=wd,
        ws=ws,
    )
    U, V = winds.interp_3d_to_calmet(threed, zface, nx, ny, 0.0, 0.0, 1.0, 0)
    # Cell (0,0) center is 0.5 km → still nearest 3D (0,0) → WD 90° → U<0, V~0
    assert U[0, 0, 0] < -0.5
    assert abs(V[0, 0, 0]) < 0.2
