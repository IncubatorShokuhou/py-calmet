"""Light Fortran/golden alignment: missing-obs sentinels must not poison winds."""
from __future__ import annotations

import numpy as np
import pytest

from py_calmet.core.winds import obs_profile_similt, obs_surface_uv
from py_calmet.io.surf import SurfRecord, surf_pres, surf_rh, surf_sky, surf_tempk
from py_calmet.io.up import UpLevel
from py_calmet.core import pbl


def test_surf_9999_is_calm_not_phantom():
    u, v = obs_surface_uv(9999.0, 9999.0, 4, 4)
    assert float(np.hypot(u, v).max()) == pytest.approx(0.0)
    u2, v2 = obs_surface_uv(3.5, 220.0, 2, 2)
    assert float(np.hypot(u2, v2).mean()) == pytest.approx(3.5, abs=1e-6)


def test_up_missing_levels_do_not_inflate_profile():
    levels = [
        UpLevel(1000.0, 50.0, 15.0, 220.0, 5.0),
        UpLevel(850.0, 1500.0, 5.0, 999.0, 999.0),  # missing mid-level
        UpLevel(700.0, 3000.0, -5.0, 250.0, 15.0),
    ]
    zface = np.array([0.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 2000.0, 3000.0])
    U, V = obs_profile_similt(
        3.0, 0.0, 10.0, 0.1, -100.0, 500.0, zface, levels, 0.0, 50.0, 3, 3, iextrp=-4
    )
    spd = np.hypot(U[:, 0, 0], V[:, 0, 0])
    # With the 999 sentinel filtered, speeds stay meteorological (not hundreds of m/s)
    assert float(spd.max()) < 40.0
    assert float(spd[0]) == pytest.approx(3.0, abs=1e-6)


def test_up_all_missing_falls_back_to_surface_powerlaw():
    levels = [
        UpLevel(1000.0, 50.0, 15.0, 999.0, 999.0),
        UpLevel(850.0, 1500.0, 5.0, 999.0, 999.0),
    ]
    zface = np.array([0.0, 20.0, 50.0, 100.0, 200.0])
    U, V = obs_profile_similt(
        4.0, 0.0, 10.0, 0.1, -100.0, 500.0, zface, levels, 0.0, 50.0, 2, 2, iextrp=-4
    )
    assert float(np.hypot(U[0, 0, 0], V[0, 0, 0])) == pytest.approx(4.0, abs=1e-6)
    assert float(np.hypot(U[:, 0, 0], V[:, 0, 0]).max()) < 20.0


def _surf_rec(**kw):
    base = dict(
        year=2020, jday=1, hour=0, ws=3.0, wd=180.0, ceil=50, sky=3,
        tempk=291.0, rh=65, pres=1012.0, ppcode=0,
    )
    base.update(kw)
    return SurfRecord(**base)


def test_surf_missing_temp_rh_pres_sky_use_defaults():
    miss = _surf_rec(tempk=9999.0, rh=9999, pres=9999.0, sky=9999)
    assert surf_tempk(miss) == pytest.approx(288.15)
    assert surf_rh(miss) == 70
    assert surf_pres(miss) == pytest.approx(1012.0)
    assert surf_sky(miss) == 0
    ok = _surf_rec()
    assert surf_tempk(ok) == pytest.approx(291.0)
    assert surf_rh(ok) == 65
    assert surf_pres(ok) == pytest.approx(1012.0)
    assert surf_sky(ok) == 3


def test_surf_missing_thermo_does_not_poison_air_density_or_elustr():
    """9999 T/P/sky must not yield nonsense rho / ustar (phantom ~10k K)."""
    miss = _surf_rec(tempk=9999.0, rh=9999, pres=9999.0, sky=9999)
    t = surf_tempk(miss)
    p = surf_pres(miss)
    sky = surf_sky(miss)
    rho = float(pbl.air_density(np.array([[t]]), p)[0, 0])
    assert 0.9 < rho < 1.5
    u = np.full((2, 2), 3.0)
    v = np.zeros((2, 2))
    z0 = np.full((2, 2), 0.1)
    temp2d = np.full((2, 2), t)
    rho2d = pbl.air_density(temp2d, p)
    ust, el, qh = pbl.elustr_stable(u, v, z0, 10.0, temp2d, rho2d, sky)
    assert float(np.nanmax(np.abs(ust))) < 5.0
    assert float(np.nanmax(np.abs(qh))) < 500.0
    # Contrast: raw 9999 would make rho tiny and heat flux explode
    rho_bad = float(pbl.air_density(np.array([[9999.0]]), 9999.0)[0, 0])
    assert rho_bad < 0.5  # poisoned path still detectable
