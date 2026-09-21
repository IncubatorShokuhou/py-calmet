"""Smoke tests for input and CALMET.DAT readers (wrf_demo / real WRF)."""
from pathlib import Path
import numpy as np
from py_calmet.io import read_geo, read_surf, read_up, read_3d, read_inp, read_calmet_dat

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "cases" / "wrf_demo" / "goldens" / "inputs"
GOLD = ROOT / "cases" / "wrf_demo" / "goldens"


def test_read_geo():
    g = read_geo(INPUTS / "geo.dat")
    assert g.nx == 12 and g.ny == 12
    assert g.elev.shape == (12, 12)
    assert g.landuse[0, 0] == 20
    # Katrina mountain window — real terrain, not flat synthetic
    assert float(g.elev.max()) > 1000.0


def test_read_surf():
    s = read_surf(INPUTS / "surf.dat")
    assert len(s.records) == 3
    assert s.records[0].ws > 0.0


def test_read_up():
    u = read_up(INPUTS / "up.dat")
    assert len(u.soundings) >= 3
    assert len(u.soundings[0].levels) == 10


def test_read_3d():
    d = read_3d(INPUTS / "3d.dat")
    assert d.ni == 14 and d.nj == 14 and d.nk == 10
    assert d.ws.shape[0] >= 3


def test_read_calmet_dat_obs():
    m = read_calmet_dat(GOLD / "obs" / "CALMET.DAT")
    assert m.nx == 12 and m.ny == 12 and m.nz == 8 and m.nt == 3
    assert m.get_3d_field("U").shape == (3, 8, 12, 12)
    assert "ZI" in m.fields_2d
    assert m.get_time_bounds()[0][0].year == 2005


def test_read_inp_modes():
    assert read_inp(GOLD / "obs" / "calmet.inp").mode == "obs"
    assert read_inp(GOLD / "obs_model" / "calmet.inp").mode == "obs_model"
    assert read_inp(GOLD / "noobs" / "calmet.inp").mode == "noobs"


def test_yyyyjjjhh_years_not_divisible_by_10():
    from datetime import datetime
    from py_calmet.io.calmet_dat import _parse_yyyyjjjhh, _yyyyjjjhh

    for dt in (
        datetime(2005, 8, 28, 0),
        datetime(2005, 8, 28, 3),
        datetime(1999, 12, 31, 23),
        datetime(2020, 6, 15, 0),
    ):
        code, sec = _yyyyjjjhh(dt)
        parsed = _parse_yyyyjjjhh(code, sec)
        assert parsed.replace(minute=0, second=0, microsecond=0) == dt.replace(
            minute=0, second=0, microsecond=0
        )


def test_wrf_demo_calmet_dat_times():
    path = GOLD / "noobs" / "CALMET.DAT"
    m = read_calmet_dat(path)
    t0 = m.get_time_bounds()[0][0]
    assert t0.year == 2005
    assert t0.month == 8
    assert t0.day == 28
    # Cell centers, not SW corners.
    assert m.x[0] == m.run_control.xorigr + 0.5 * m.run_control.dgrid
    assert m.y[0] == m.run_control.yorigr + 0.5 * m.run_control.dgrid


def test_package_version_matches_pyproject():
    from py_calmet import __version__

    assert __version__ == "1.0.0"
