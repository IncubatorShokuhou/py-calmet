"""WP6: IDIOPT1/4/5 DIAG.DAT ingestion + WTDAT water-T precedence + Part B."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from py_calmet.core import diag_opts, overwater, winds
from py_calmet.core.met_utils import G
from py_calmet.io.diag_dat import read_diag, diag_for_hour
from py_calmet.io.wt_dat import read_wt, wt_for_hour, wt_sst_grid


def test_read_diag_positional_and_keyword(tmp_path: Path):
    p = tmp_path / "diag.dat"
    p.write_text(
        "DIAG.DAT 1.0\n"
        "* comment\n"
        "2020 167 0 288.15 0.005 2.0 0.5 3.0 0.0 5.0 1.0\n"
        "2020 167 1 TSFC=290.0 USFC=4.0 VSFC=-1.0 UUP=6.0 VUP=2.0 GAMMA=0.01\n"
    )
    data = read_diag(p)
    assert len(data.records) == 2
    r0 = diag_for_hour(data, 2020, 167, 0)
    assert r0 is not None
    assert r0.has_tsfc() and r0.tsfc == pytest.approx(288.15)
    assert r0.has_sfc_uv() and r0.usfc == pytest.approx(3.0)
    assert r0.has_up_uv() and r0.uup == pytest.approx(5.0)
    r1 = diag_for_hour(data, 2020, 167, 1)
    assert r1 is not None and r1.tsfc == pytest.approx(290.0)
    assert r1.gamma == pytest.approx(0.01)


def test_idiopt1_applies_diag_tsfc():
    t = np.full((3, 3), 280.0)
    out = diag_opts.apply_diag_sfc_temp(t, idiopt1=1, tsfc=295.0)
    assert float(out.mean()) == pytest.approx(295.0)
    out0 = diag_opts.apply_diag_sfc_temp(t, idiopt1=0, tsfc=295.0)
    assert float(out0.mean()) == pytest.approx(280.0)


def test_idiopt4_5_gated_by_irtype():
    u = np.ones((4, 4))
    v = np.zeros((4, 4))
    u2, v2 = diag_opts.apply_diag_sfc_uv(
        u, v, idiopt4=1, irtype=1, usfc=7.0, vsfc=2.0
    )
    assert float(u2.mean()) == pytest.approx(1.0)
    u3, v3 = diag_opts.apply_diag_sfc_uv(
        u, v, idiopt4=1, irtype=0, usfc=7.0, vsfc=2.0
    )
    assert float(u3.mean()) == pytest.approx(7.0)
    assert float(v3.mean()) == pytest.approx(2.0)

    U = np.zeros((3, 2, 2))
    V = np.zeros((3, 2, 2))
    U[0] = 1.0
    U2, V2 = diag_opts.apply_diag_upper_uv(
        U, V, idiopt5=1, irtype=0, uup=9.0, vup=-3.0
    )
    assert float(U2[1].mean()) == pytest.approx(9.0)
    assert float(U2[0].mean()) == pytest.approx(1.0)


def test_idiopt2_uses_diag_gamma():
    g, note = diag_opts.resolve_diag_gamma(
        idiopt2=1,
        zupt=200.0,
        ziconv_mean=0.0,
        sounding_z=None,
        sounding_t=None,
        daytime=True,
        diag_gamma=0.008,
    )
    assert "DIAG.DAT" in note
    assert g == pytest.approx(0.008)
    g2, note2 = diag_opts.resolve_diag_gamma(
        idiopt2=1,
        zupt=200.0,
        ziconv_mean=0.0,
        sounding_z=None,
        sounding_t=None,
        daytime=True,
        diag_gamma=None,
    )
    assert "missing" in note2
    assert g2 == pytest.approx(0.005)


def test_qa_idiopt_reports_ingestion():
    notes = diag_opts.qa_idiopt([1, 0, 0, 1, 0], irtype=0, diag_loaded=True)
    assert any("ingesting" in n and "IDIOPT1" in n for n in notes)
    notes2 = diag_opts.qa_idiopt([0, 0, 0, 1, 0], irtype=1, diag_loaded=True)
    assert any("forbids" in n for n in notes2)


def test_read_wt_and_sst_grid(tmp_path: Path):
    p = tmp_path / "wt.dat"
    p.write_text("WT.DAT 1.0\n* sst\n2020 167 0 285.0\n2020 167 1 TSEA=286.5\n")
    data = read_wt(p)
    assert len(data.records) == 2
    r = wt_for_hour(data, 2020, 167, 0)
    assert r is not None and r.tsea == pytest.approx(285.0)
    lu = np.array([[10, 55], [55, 10]])
    fb = np.full((2, 2), 290.0)
    sst = wt_sst_grid(r, 2, 2, lu, 55, 55, fb)
    assert sst[0, 1] == pytest.approx(285.0)
    assert sst[0, 0] == pytest.approx(290.0)


def test_water_t_precedence_wt_over_air():
    class _Rec:
        def has_tsea(self):
            return True

        tsea = 283.0
        grid = None

    lu = np.full((2, 2), 55)
    tempk = np.full((2, 2), 290.0)
    t_sea, src, _, _ = overwater.resolve_water_temp(
        itwprog=0,
        landuse=lu,
        iwat1=55,
        iwat2=55,
        tempk=tempk,
        sea_records=None,
        wt_record=_Rec(),
    )
    assert src == "WT.DAT"
    assert t_sea is not None
    assert float(t_sea.mean()) == pytest.approx(283.0)


def test_water_t_sea_beats_wt():
    class _Sea:
        x_km = 0.5
        y_km = 0.5
        twave = -999.0
        hwave = -999.0

        @property
        def t_sea(self):
            return 280.0

    class _Wt:
        def has_tsea(self):
            return True

        tsea = 283.0
        grid = None

    lu = np.full((2, 2), 55)
    tempk = np.full((2, 2), 290.0)
    t_sea, src, _, _ = overwater.resolve_water_temp(
        itwprog=0,
        landuse=lu,
        iwat1=55,
        iwat2=55,
        tempk=tempk,
        xorig_km=0.0,
        yorig_km=0.0,
        dgrid_km=1.0,
        sea_records=[_Sea()],
        wt_record=_Wt(),
    )
    assert src == "SEA.DAT"
    assert float(t_sea.mean()) == pytest.approx(280.0)


def test_part_b_g_constant_in_froude():
    assert G == pytest.approx(9.81)
    src = Path(winds.__file__).read_text()
    assert "sqrt(G *" in src
    assert "(G / temp)" in src


def test_part_b_irtype0_zeros_like_pattern():
    runner = Path(__file__).resolve().parents[1] / "py_calmet" / "core" / "runner.py"
    text = runner.read_text()
    assert "wstar = np.zeros_like(zi)" in text
    assert "wstar = np.zeros_like(wstar)" not in text
