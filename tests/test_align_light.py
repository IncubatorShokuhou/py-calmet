"""Light Fortran/golden alignment: missing-obs sentinels must not poison winds."""
from __future__ import annotations

import numpy as np
import pytest

from py_calmet.core.winds import obs_profile_similt, obs_surface_uv
from py_calmet.io.up import UpLevel


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
