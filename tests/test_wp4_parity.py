"""WP4: Partial→Implemented wiring (listfile, OA extras, Zi ops, IGF, OutOfScope MM4)."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from py_calmet.config import CalmetConfig
from py_calmet.core import run_options, zi_ops, winds
from py_calmet.io.igf import read_igf_header
from py_calmet.io.metlst import write_metlst


def test_mm4dat_out_of_scope():
    CalmetConfig().check_unsupported()
    with pytest.raises(NotImplementedError, match="MM4DAT"):
        CalmetConfig(mm4dat="legacy_mm5.dat").check_unsupported()


def test_igfmet_no_longer_hard_stop():
    CalmetConfig(igfmet=1).check_unsupported()


def test_lcfiles_casefold(tmp_path: Path):
    (tmp_path / "GEO.DAT").write_text("x")
    p = run_options.resolve_data_file(
        tmp_path, tmp_path, "geo.dat", "geo.dat", lcfiles=True
    )
    assert p is not None and p.name == "GEO.DAT"


def test_fcoriol_honors_inp():
    class _I:
        def get_float(self, k, d=0.0):
            return 1.0e-4 if k == "FCORIOL" else d

    assert run_options.effective_fcoriol(_I(), 40.0) == pytest.approx(1e-4)


def test_zi_average_noop_and_relax():
    zi = np.full((4, 4), 500.0)
    u = np.ones((2, 4, 4))
    v = np.zeros((2, 4, 4))
    assert np.allclose(
        zi_ops.average_zi_upwind(zi, u, v, iavezi=1, mnmdav=1), zi
    )
    z2 = zi_ops.average_zi_upwind(zi, u, v, iavezi=1, mnmdav=3, hafang=90.0)
    assert z2.shape == zi.shape
    prev = np.full_like(zi, 100.0)
    out = zi_ops.relax_zi(
        zi, prev, izicrlx=1, tzicrlx=800.0, dt_sec=3600.0, daytime=True
    )
    assert float(out.mean()) < 500.0
    assert float(out.mean()) > 100.0


def test_holzworth_mixht():
    z = np.array([0.0, 200.0, 500.0, 1000.0, 2000.0])
    t = np.array([290.0, 288.0, 286.0, 285.0, 284.0])
    zi = zi_ops.mixht_holzworth(
        np.full((2, 2), 290.0), z, t, zimin=50.0, zimax=3000.0
    )
    assert zi.shape == (2, 2)
    assert float(zi.mean()) > 50.0


def test_oa_rmin_lvary_icalm():
    ug = np.ones((2, 3, 3))
    vg = np.zeros((2, 3, 3))
    uo = np.full((2, 3, 3), 5.0)
    vo = np.zeros((2, 3, 3))
    U, V = winds.objective_analyze(
        ug,
        vg,
        uo,
        vo,
        xs_m=1500.0,
        ys_m=1500.0,
        xorig_m=0.0,
        yorig_m=0.0,
        dgrid_m=1000.0,
        r1_m=5000.0,
        r2_m=5000.0,
        rmin_m=100.0,
        lvary=True,
        icalm=0,
        rmax1_m=10000.0,
        rmax2_m=10000.0,
    )
    assert U.shape == ug.shape
    assert float(np.max(U)) > 1.0


def test_nflagp_and_temp_average():
    rates = np.array([-1.0, 0.005, 2.0])
    out = run_options.apply_nflagp(rates, nflagp=2, cutp=0.01)
    assert out.tolist() == [0.0, 0.0, 2.0]
    t = np.arange(16.0).reshape(4, 4)
    assert np.allclose(
        run_options.average_temperature(t, iavet=1, numts=1), t
    )
    sm = run_options.average_temperature(
        t, iavet=1, numts=2, tradkm=5.0, dgridkm=1.0
    )
    assert sm.shape == t.shape
    assert float(sm.std()) < float(t.std())


def test_metlst_lprint_ipr(tmp_path: Path):
    write_metlst(
        tmp_path / "a.lst",
        mode="obs",
        nhrs=2,
        lprint=True,
        ipr_flags={"IPR0": 1, "IPR3": 1},
        print_fields={"USTAR": True, "MIXHT": True},
        field_samples={"ZI": np.array([[100.0, 200.0]])},
        iuvout=[1, 2],
    )
    text = (tmp_path / "a.lst").read_text()
    assert "LPRINT = True" in text
    assert "IPR0 = 1" in text
    assert "ZI:" in text


def test_igf_header_missing(tmp_path: Path):
    hdr = read_igf_header(tmp_path / "missing.dat")
    assert hdr.ok is False
    assert "missing" in hdr.notes[0].lower()


def test_ibtz_and_irlg_window():
    class _I:
        def __init__(self):
            self._d = {
                "IBYR": 2020,
                "IBMO": 1,
                "IBDY": 1,
                "IBHR": 0,
                "IBSEC": 0,
                "IEYR": 2020,
                "IEMO": 1,
                "IEDY": 1,
                "IEHR": 0,
                "IESEC": 0,
                "NSECDT": 3600,
                "IRLG": 5,
                "IBTZ": 5,
                "ABTZ": None,
            }

        def get_int(self, k, d=0):
            v = self._d.get(k, d)
            return int(d if v is None else v)

        def get(self, k, d=None):
            return self._d.get(k, d)

    start, end, nhrs, _ = run_options.run_window_with_legacy(_I())
    assert nhrs == 5
    assert run_options.ibtz_hours(_I()) == 5


def test_grid_qa_notes():
    class _G:
        nx = 10
        ny = 10
        dgridkm = 1.0
        xorigkm = 0.0
        yorigkm = 0.0

    class _I:
        def get_int(self, k, d=0):
            return 12 if k in ("NX", "NY") else d

        def get_float(self, k, d=0.0):
            return d

    notes = run_options.qa_grid_vs_geo(_I(), _G())
    assert any("NX" in n for n in notes)


def test_test_stubs(tmp_path: Path):
    class _I:
        def get(self, k, d=None):
            return None

    written = run_options.write_test_stubs(tmp_path, _I(), enabled=True)
    assert len(written) >= 5
    assert any(Path(p).is_file() for p in written)
