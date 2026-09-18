"""WP3: MIXDT, SEA/PRECIP/CLOUD I/O, barriers, lake breeze, coord, IOUTMM5, outputs."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from py_calmet.core import mixdt, barriers, overwater, coord
from py_calmet.core.pbl import mixht_day_carson
from py_calmet.io.sea import read_sea, SeaRecord, pick_sea_record
from py_calmet.io.precip_dat import read_precip, rates_for_hour
from py_calmet.io.cloud_dat import read_cloud, write_cloud, CloudData, CloudRecord, cloud_for_hour
from py_calmet.io.metlst import write_metlst
from py_calmet.io.pacout import write_pacout, mixed_layer_uv
from py_calmet.io.threed import _parse_upper_line, _parse_iout_flags
from py_calmet.config import CalmetConfig


def test_mixdt_stable_lapse_from_inversion():
    # Strong inversion above 500 m: T rises with height → large gamma
    zl = np.array([0.0, 200.0, 500.0, 700.0, 1000.0])
    tz = np.array([290.0, 288.0, 286.0, 288.0, 290.0])  # inversion above 500
    tht, thtp, g = mixdt.mixdt_sounding(zl, tz, htold=500.0, dptmin=0.001, dzzi=200.0)
    assert g > 0.01  # pot-temp lapse stronger than dry adiabatic alone
    assert thtp > tht


def test_mixdt_floors_at_dptmin():
    zl = np.array([0.0, 500.0, 1000.0, 1500.0])
    tz = np.array([300.0, 295.0, 290.0, 285.0])  # dry-adiabatic-ish
    _, _, g = mixdt.mixdt_sounding(zl, tz, htold=400.0, dptmin=0.005, dzzi=200.0)
    assert g >= 0.005


def test_mixdt2_and_gamma_field():
    z = np.linspace(0, 2000, 8)
    t = 290.0 - 0.006 * z
    tsf, tht, thtp, g = mixdt.mixdt2_column(z, t, 300.0, dptmin=0.001, dzzi=200.0)
    assert tsf == pytest.approx(290.0, abs=0.1)
    zi = np.full((3, 3), 300.0)
    gf = mixdt.gamma_field_from_sounding(z, t, zi, dptmin=0.001, dzzi=200.0)
    assert gf.shape == (3, 3)
    assert np.all(gf >= 0.001)


def test_mixht_day_accepts_2d_gamma():
    qh = np.full((4, 4), 100.0)
    rho = np.full((4, 4), 1.2)
    t = np.full((4, 4), 290.0)
    ust = np.full((4, 4), 0.3)
    gamma = np.full((4, 4), 0.01)
    zi, zic, dptt = mixht_day_carson(
        qh, rho, t, ust, 1e-4, dt_sec=3600.0,
        ziconv_prev=np.full((4, 4), 200.0),
        dtheta=gamma, threshl=0.05,
    )
    assert zi.shape == (4, 4)
    assert float(zi.mean()) > 200.0


def test_barrier_blocks_opposite_side():
    bar = barriers.BarrierSet(
        xbbar=np.array([0.0]),
        ybbar=np.array([0.0]),
        xebar=np.array([0.0]),
        yebar=np.array([10.0]),  # N-S barrier at x=0
        kbar=5,
    )
    # Station at x=-1, cell at x=+1 → blocked
    assert barriers.same_side(1.0, 5.0, -1.0, 5.0, bar) is False
    # Same side → clear
    assert barriers.same_side(1.0, 5.0, 2.0, 5.0, bar) is True
    mask = barriers.station_clear_mask(1.0, 5.0, np.array([-1.0, 2.0]), np.array([5.0, 5.0]), bar, 0)
    assert mask.tolist() == [False, True]


def test_lake_breeze_blends_inside_box():
    U = np.zeros((2, 5, 5))  # calm → onshore blend should inject wind
    V = np.zeros_like(U)
    box = barriers.LakeBreezeBox(
        xg1=0.0, xg2=5.0, yg1=0.0, yg2=5.0,
        xbcst=0.0, ybcst=2.5, xecst=0.0, yecst=2.5,
        metbxid=[],
    )
    cfg = barriers.LakeBreezeConfig(boxes=[box])
    U2, V2 = barriers.apply_lake_breeze(
        U, V, xorig_km=0.0, yorig_km=0.0, dgrid_km=1.0, cfg=cfg, blend=0.5
    )
    assert U2.shape == U.shape
    # Surface changed somewhere (onshore injection)
    assert float(np.max(np.abs(U2[0]))) > 0.1 or float(np.max(np.abs(V2[0]))) > 0.1


def test_coord_utm_roundtrip_and_lcc():
    proj = coord.MapProjection(pmap="UTM", iutmzn=19, utmhem="N")
    lat, lon = coord.project_xy_to_ll(500000.0, 5000000.0, proj)
    assert 40.0 < lat < 50.0
    lcc = coord.MapProjection(
        pmap="LCC", rlat0=40.0, rlon0=-100.0, xlat1=30.0, xlat2=50.0
    )
    x, y = coord.project_ll_to_xy(40.0, -100.0, lcc)
    lat2, lon2 = coord.project_xy_to_ll(x, y, lcc)
    assert lat2 == pytest.approx(40.0, abs=0.05)
    assert lon2 == pytest.approx(-100.0, abs=0.05)


def test_coare_iwarm_icool_change_flux():
    u = np.full((3, 3), 6.0)
    v = np.zeros((3, 3))
    t_air = np.full((3, 3), 290.0)
    t_sea = np.full((3, 3), 290.0)
    ust0, qh0, _ = overwater.coare_lite_fluxes(u, v, t_air, t_sea, iwarm=0, icool=0, qsw=800.0)
    ust1, qh1, _ = overwater.coare_lite_fluxes(u, v, t_air, t_sea, iwarm=1, icool=1, qsw=800.0)
    # Warm layer heats skin → more upward flux when air≈sea bulk
    assert float(qh1.mean()) > float(qh0.mean())


def test_sea_precip_cloud_io(tmp_path: Path):
    # SEA.DAT mini
    sea = tmp_path / "sea.dat"
    sea.write_text(
        "SEA.DAT         2.1             test\n"
        "1\n"
        "comment\n"
        "UTM\n"
        "  19N\n"
        "WGS-84  10-04-2020\n"
        "KM  \n"
        "101 SEA1\n"
        "100.0 200.0 10.0 2020 167 0 2020 167 1 -2.0 290.0 80.0 500.0 0.01 0.005 5.0 180.0 8.0 1.5\n"
    )
    st = read_sea(sea)
    assert st.version == pytest.approx(2.1)
    assert len(st.records) == 1
    assert st.records[0].t_sea == pytest.approx(292.0)  # 290 - (-2)
    r = pick_sea_record(st, 2020, 167, 0)
    assert r is not None

    # PRECIP.DAT
    prc = tmp_path / "precip.dat"
    prc.write_text("2020 167 0 1.5 0.0\n2020 167 1 0.5 0.2\n")
    pd = read_precip(prc, npsta=2)
    assert pd.nsta == 2
    rates = rates_for_hour(pd, 2020, 167, 0)
    assert rates[0] == pytest.approx(1.5)

    # CLOUD.DAT round-trip
    cc = np.linspace(0, 1, 12).reshape(3, 4)
    cd = CloudData(nx=4, ny=3, records=[CloudRecord(2020, 167, 0, cc)])
    cpath = tmp_path / "cloud.dat"
    write_cloud(cpath, cd)
    cd2 = read_cloud(cpath, 4, 3)
    assert len(cd2.records) == 1
    assert cd2.records[0].ccgrid.shape == (3, 4)
    assert cloud_for_hour(cd2, 2020, 167, 0) is not None


def test_metlst_pacout(tmp_path: Path):
    write_metlst(tmp_path / "list.txt", mode="noobs", nhrs=3, grid={"nx": 2})
    assert (tmp_path / "list.txt").is_file()
    U = np.ones((2, 3, 4, 4))
    V = np.zeros_like(U)
    zi = np.full((2, 4, 4), 500.0)
    zface = np.array([0.0, 20.0, 50.0, 100.0])
    um, vm = mixed_layer_uv(U, V, zi, zface)
    assert um.shape == (2, 4, 4)
    write_pacout(
        tmp_path / "pac.npz",
        U_mix=um, V_mix=vm, zi=zi,
        ustar=np.ones_like(zi) * 0.2,
        wstar=np.zeros_like(zi),
        el=np.full_like(zi, -100.0),
        ipgt=np.full(zi.shape, 4, dtype=np.int32),
        rmm=np.zeros_like(zi),
        rho=np.ones_like(zi),
        tempk=np.full_like(zi, 290.0),
        qsw=np.zeros_like(zi),
        irh=np.full(zi.shape, 70, dtype=np.int32),
    )
    assert list(tmp_path.glob("pac*.npz")) or list(tmp_path.glob("pac*.npz.npz")) or (tmp_path / "pac.npz").exists() or (tmp_path / "pac.npz.npz").exists()


def test_ioutmm5_parse_variants():
    assert _parse_iout_flags("1 1 0 0 0") == 92  # 81+10+1
    line81 = f"{1000:4d}{100:6d}{280.0:6.1f}{270:4d}{5.0:5.1f}"
    p = _parse_upper_line(line81, 81)
    assert p["temp"] == pytest.approx(280.0)
    assert p["ws"] == pytest.approx(5.0)
    # format 92 with w + rh + q
    line92 = f"{1000:4d}{100:6d}{280.0:6.1f}{270:4d}{5.0:5.1f}{0.10:6.2f}{70:3d}{5.00:5.2f}"
    p92 = _parse_upper_line(line92, 92)
    assert p92["w"] == pytest.approx(0.10, abs=0.01)
    assert p92["rh"] == pytest.approx(70)
    # format 83: P Z T WD WS RH Q QC QR (no W)
    line83 = (
        f"{1000:4d}{100:6d}{280.0:6.1f}{270:4d}{5.0:5.1f}"
        f"{70:3d}{5.00:5.2f}{0.100:6.3f}{0.050:6.3f}"
    )
    p83 = _parse_upper_line(line83, 83)
    assert p83["rh"] == pytest.approx(70)
    assert p83["qc"] == pytest.approx(0.1, abs=0.01)


def test_config_allows_llbreze_and_nbar():
    cfg = CalmetConfig(llbreze=True, nbar=1, igfmet=0, imixh=1)
    cfg.check_unsupported()  # should not raise
    CalmetConfig(igfmet=1).check_unsupported()  # IGFMET now implemented
    with pytest.raises(NotImplementedError):
        CalmetConfig(mm4dat="foo_mm4.dat").check_unsupported()
