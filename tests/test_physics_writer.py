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
TINY = ROOT / "cases" / "small_domain" / "goldens"


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
    zi0, zc0 = pbl.mixht_day_carson(qh, rho, tempk, ustar, 1e-4, dt_sec=3600.0)
    zi1, zc1 = pbl.mixht_day_carson(
        qh, rho, tempk, ustar, 1e-4, dt_sec=3600.0, ziconv_prev=zc0
    )
    assert zc1.mean() >= zc0.mean()
    assert zi1.mean() >= 50.0


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
