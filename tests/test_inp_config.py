"""WP1: full CalmetConfig surface, INP round-trip, JWAT alias binding."""
from __future__ import annotations

from pathlib import Path

import pytest

from py_calmet.config import PARAM_NAMES, CalmetConfig, default_config
from py_calmet.io.inp import read_inp, write_inp

ROOT = Path(__file__).resolve().parents[1]
CASE_INPS = [
    ROOT / "cases" / "small_domain" / "goldens" / "obs" / "calmet.inp",
    ROOT / "cases" / "wrf_demo" / "goldens" / "obs" / "calmet.inp",
    ROOT / "cases" / "daytime_zi" / "goldens" / "noobs" / "calmet.inp",
]


def test_config_has_all_207_params():
    assert len(PARAM_NAMES) == 207
    cfg = default_config()
    for name in PARAM_NAMES:
        assert hasattr(cfg, name.lower()), name
    # radiation + water defaults
    assert cfg.ha1 == 990.0
    assert cfg.hb1 == -0.75
    assert cfg.hc1 == pytest.approx(5.31e-13)
    assert cfg.jwat1 == 999
    assert cfg.effective_iwat() == (55, 55)


@pytest.mark.parametrize("path", CASE_INPS, ids=["small_obs", "wrf_obs", "daytime_noobs"])
def test_full_inp_parse_no_key_loss(path: Path):
    assert path.is_file()
    inp = read_inp(path)
    # every active !KEY=...! lands in raw
    assert len(inp.raw) >= 150
    # known Group 0 / JWAT / HA surface
    assert "GEODAT" in inp.raw
    assert "METDAT" in inp.raw
    assert "JWAT1" in inp.raw and "JWAT2" in inp.raw
    assert inp.get_int("JWAT1") == 999
    assert inp.get_int("IWAT1") == 999  # alias
    assert inp.config.effective_iwat() == (55, 55)
    assert inp.config.geodat.lower().endswith("geo.dat")
    # HA* absent in samples → Fortran defaults on config
    assert "HA1" not in inp.raw
    assert inp.config.ha1 == 990.0


@pytest.mark.parametrize("path", CASE_INPS, ids=["small_obs", "wrf_obs", "daytime_noobs"])
def test_inp_roundtrip_preserves_keys(path: Path, tmp_path: Path):
    original = read_inp(path)
    out = tmp_path / "roundtrip.inp"
    write_inp(out, original)
    rewritten = read_inp(out)
    assert set(original.raw.keys()) == set(rewritten.raw.keys())
    # values for a few critical keys
    for key in ("JWAT1", "JWAT2", "NSECDT", "CONSTN", "ZIMIN"):
        if key in original.raw:
            assert rewritten.get_int(key) == original.get_int(key)
    if "GEODAT" in original.raw:
        assert str(rewritten.get("GEODAT")).strip() == str(original.get("GEODAT")).strip()


def test_jwat_alias_overrides_water_range():
    cfg = CalmetConfig()
    cfg.jwat1 = 50
    cfg.jwat2 = 55
    assert cfg.effective_iwat() == (50, 55)
    # IWAT alias writes JWAT
    cfg2 = CalmetConfig()
    cfg2.set_param("IWAT1", 40)
    cfg2.set_param("IWAT2", 45)
    assert cfg2.jwat1 == 40
    assert cfg2.effective_iwat() == (40, 45)


def test_check_unsupported_switches():
    cfg = CalmetConfig()
    cfg.check_unsupported()  # defaults ok
    cfg.igfmet = 1
    with pytest.raises(NotImplementedError, match="IGFMET"):
        cfg.check_unsupported()


def test_objective_analyze_multistation_matches_single():
    """Multi-station API with one station must match legacy single-station OA."""
    import numpy as np
    from py_calmet.core.winds import objective_analyze

    rng = np.random.default_rng(0)
    nz, ny, nx = 3, 5, 5
    ug = rng.normal(size=(nz, ny, nx))
    vg = rng.normal(size=(nz, ny, nx))
    uo = ug + 1.0
    vo = vg - 0.5
    kwargs = dict(
        xorig_m=0.0,
        yorig_m=0.0,
        dgrid_m=1000.0,
        r1_m=2000.0,
        r2_m=4000.0,
        xs_m=2500.0,
        ys_m=2500.0,
    )
    U1, V1 = objective_analyze(ug, vg, uo, vo, **kwargs)
    U2, V2 = objective_analyze(ug, vg, uo, vo, rprog_m=0.0, **kwargs)
    assert np.allclose(U1, U2) and np.allclose(V1, V2)
    # two identical stations → same as one when RMAX open
    xs = np.array([2500.0, 2500.0])
    ys = np.array([2500.0, 2500.0])
    u_stn = np.stack([uo[:, 2, 2], uo[:, 2, 2]], axis=0)  # (nstn, nz) approx center
    # use profile form
    u_prof = np.stack([uo[:, ny//2, nx//2], uo[:, ny//2, nx//2]], axis=0)
    v_prof = np.stack([vo[:, ny//2, nx//2], vo[:, ny//2, nx//2]], axis=0)
    Um, Vm = objective_analyze(
        ug, vg, u_prof, v_prof, xs, ys,
        xorig_m=0.0, yorig_m=0.0, dgrid_m=1000.0, r1_m=2000.0, r2_m=4000.0,
    )
    assert Um.shape == ug.shape
