"""High-level run API matching obs / obs_model / noobs modes."""
from __future__ import annotations
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional
import numpy as np

from ..io.geo import read_geo
from ..io.surf import read_surf, surf_tempk, surf_rh, surf_pres, surf_sky
from ..io.up import read_up
from ..io.threed import read_3d
from ..io.inp import read_inp
from ..io.sea import read_sea, pick_sea_record
from ..io.precip_dat import read_precip, rates_for_hour
from ..io.cloud_dat import read_cloud, write_cloud, cloud_for_hour, CloudData, CloudRecord
from ..io.diag_dat import read_diag, diag_for_hour
from ..io.wt_dat import read_wt, wt_for_hour
from ..io.metlst import write_metlst
from ..io.pacout import write_pacout, mixed_layer_uv
from .met_utils import layer_mids, coriolis, utm_to_latlon
from . import winds, pbl, clouds, precip, overwater, mixdt, barriers as barriers_mod
from .coord import MapProjection, domain_center_latlon
from . import run_options, zi_ops, diag_opts
from ..io import igf as igf_mod


@dataclass
class CalmetResult:
    mode: str
    zface: np.ndarray
    U: np.ndarray  # [nt, nz, ny, nx]
    V: np.ndarray
    W: np.ndarray
    T: np.ndarray
    IPGT: np.ndarray
    USTAR: np.ndarray
    ZI: np.ndarray
    EL: np.ndarray
    WSTAR: np.ndarray
    TEMPK: np.ndarray
    RHO: np.ndarray
    QSW: np.ndarray
    IRH: np.ndarray
    elev: np.ndarray
    z0: np.ndarray
    RMM: np.ndarray | None = None
    meta: dict = field(default_factory=dict)


def _default_z0(landuse: np.ndarray) -> np.ndarray:
    z0 = np.full(landuse.shape, 0.25, dtype=np.float64)
    z0 = np.where(landuse == 10, 1.0, z0)
    z0 = np.where(landuse == 55, 0.001, z0)  # water
    return z0


def _pick_up_sounding(up, year, month, day, hour):
    best = None
    best_key = None
    target = (year, month, day, hour)
    for s in up.soundings:
        key = (s.year, s.month, s.day, s.hour)
        if key <= target and (best_key is None or key > best_key):
            best, best_key = s, key
    return best or up.soundings[0]


def _estimate_latlon(inp, geo, nx, ny):
    """Return (lat0, lon0_east) at the domain center.

    Honors RLAT0/RLON0 when set; else inverts via INP-driven MapProjection
    (UTM / LCC / TM / …).
    """
    proj = MapProjection.from_inp(inp)
    if abs(proj.rlat0) + abs(proj.rlon0) > 0.0 and -90 <= proj.rlat0 <= 90:
        return float(proj.rlat0), float(proj.rlon0)
    return domain_center_latlon(
        geo.xorigkm, geo.yorigkm, nx, ny, geo.dgridkm, proj
    )


def _run_window(inp):
    """(start, end, nhrs, nsecdt) from IBYR…IESEC (cross-midnight safe)."""
    ibyr = inp.get_int("IBYR", 2020)
    ibmo = inp.get_int("IBMO", 6)
    ibdy = inp.get_int("IBDY", 15)
    ibhr = inp.get_int("IBHR", 0)
    ibsec = inp.get_int("IBSEC", 0)
    ieyr = inp.get_int("IEYR", ibyr)
    iemo = inp.get_int("IEMO", ibmo)
    iedy = inp.get_int("IEDY", ibdy)
    iehr = inp.get_int("IEHR", ibhr + 3)
    iesec = inp.get_int("IESEC", 0)
    nsecdt = inp.get_int("NSECDT", 3600)
    start = datetime(ibyr, ibmo, ibdy, ibhr) + timedelta(seconds=ibsec)
    end = datetime(ieyr, iemo, iedy, iehr) + timedelta(seconds=iesec)
    span_sec = max(0, int((end - start).total_seconds()))
    nhrs = max(1, span_sec // max(nsecdt, 1))
    return start, end, nhrs, nsecdt


def _ibtz_from_abtz(abtz: str) -> int:
    """Parse ABTZ like UTC-0500 / UTC+0000 → CALMET IBTZ (hours west of GMT)."""
    s = (abtz or "UTC+0000").upper().replace(" ", "")
    if "UTC" not in s:
        return 0
    # UTC-0500 means local = UTC-5 → ibtz=5; UTC+0000 → 0; UTC+0800 → -8
    sign = 1
    rest = s.split("UTC", 1)[-1]
    if rest.startswith("-"):
        sign = 1
        rest = rest[1:]
    elif rest.startswith("+"):
        sign = -1
        rest = rest[1:]
    try:
        hours = int(rest[:2])
    except Exception:
        hours = 0
    return sign * hours


def _resolve_data_file(
    inputs_dir: Path,
    case_dir: Path,
    name: str | None,
    fallback: str,
    *,
    lcfiles: bool = True,
) -> Path | None:
    """Honor INP filename when present; LCFILES enables case-insensitive match."""
    return run_options.resolve_data_file(
        inputs_dir, case_dir, name, fallback, lcfiles=lcfiles
    )


def _parse_surface_stations(inp, geo, nx, ny):
    """Collect SS* station X/Y (km) and optional anemometer height.

    Returns lists xs_km, ys_km (may be length 1 = domain center fallback).
    ISURFT (1-based) selects which station is used for domain T when >0.
    """
    xs, ys = [], []
    for key, val in sorted(inp.raw.items()):
        if not key.startswith("SS") or not key[2:].isdigit():
            continue
        parts = str(val).replace("'", " ").split()
        try:
            # name id x y ...
            x = float(parts[2]); y = float(parts[3])
            xs.append(x); ys.append(y)
        except Exception:
            continue
    if not xs:
        xs = [geo.xorigkm + nx * geo.dgridkm / 2.0]
        ys = [geo.yorigkm + ny * geo.dgridkm / 2.0]
    return xs, ys


def _parse_precip_stations(inp):
    """PS* records → (xs_km, ys_km) lists; rates come from PRECIP.DAT later."""
    xs, ys = [], []
    for key, val in sorted(inp.raw.items()):
        if not key.startswith("PS") or not key[2:].isdigit():
            continue
        parts = str(val).replace("'", " ").split()
        try:
            xs.append(float(parts[2])); ys.append(float(parts[3]))
        except Exception:
            continue
    return xs, ys


def run_calmet(
    case_dir: str | Path,
    mode: Optional[str] = None,
    inputs_dir: Optional[str | Path] = None,
    *,
    write_outputs: bool = False,
) -> CalmetResult:
    """Run pure-NumPy diagnostic for a case directory.

    When ``write_outputs`` is True, honor METLST / IFORMO=2 PACOUT / ICLDOUT
    writers into ``case_dir``. Default False keeps golden directories pristine.
    """
    case_dir = Path(case_dir)
    inp = read_inp(case_dir / "calmet.inp")
    if mode is None:
        mode = inp.mode

    if inputs_dir is None:
        for cand in (case_dir / ".." / "shared", case_dir / ".." / "goldens" / "inputs", case_dir):
            if (cand / "geo.dat").exists():
                inputs_dir = cand
                break
    inputs_dir = Path(inputs_dir)

    cfg = getattr(inp, "config", None)
    if cfg is not None and hasattr(cfg, "check_unsupported"):
        cfg.check_unsupported()
    run_options.reject_mm4_mm5(inp)
    lcfiles = inp.get_bool("LCFILES", True)
    metinp = inp.get("METINP") or "calmet.inp"
    qa_notes: list[str] = [f"METINP={metinp}"]
    wt_data = None
    diag_data = None

    geo_path = _resolve_data_file(
        inputs_dir, case_dir, inp.get("GEODAT"), "geo.dat", lcfiles=lcfiles
    )
    if geo_path is None:
        raise FileNotFoundError(f"GEODAT/geo.dat not found under {case_dir} or {inputs_dir}")
    geo = read_geo(geo_path)
    nx, ny = geo.nx, geo.ny
    qa_notes.extend(run_options.qa_grid_vs_geo(inp, geo))
    zface = np.asarray(
        inp.get_list_float("ZFACE") or [0, 20, 40, 80, 160, 300, 600, 1000, 1500],
        dtype=np.float64,
    )
    nz = inp.get_int("NZ", len(zface) - 1)
    zmid = layer_mids(zface)
    z0 = _default_z0(geo.landuse)
    elev = geo.elev

    start, _end, nhrs, nsecdt = run_options.run_window_with_legacy(inp)
    ibyr, ibmo, ibdy = start.year, start.month, start.day

    srf_path = _resolve_data_file(inputs_dir, case_dir, inp.get("SRFDAT"), "surf.dat", lcfiles=lcfiles)
    up_path = _resolve_data_file(inputs_dir, case_dir, inp.get("UPDAT"), "up.dat", lcfiles=lcfiles)
    m3d_path = _resolve_data_file(inputs_dir, case_dir, inp.get("M3DDAT"), "3d.dat", lcfiles=lcfiles)
    nusta = inp.get_int("NUSTA", 0)
    nm3d = inp.get_int("NM3D", 0)
    nigf = inp.get_int("NIGF", 0)
    if cfg is not None:
        up_names = run_options.multi_file_list(getattr(cfg, "updat", None), "UPDAT", inp, nusta)
        m3d_names = run_options.multi_file_list(getattr(cfg, "m3ddat", None), "M3DDAT", inp, nm3d)
        if len(up_names) > 1:
            qa_notes.append(f"NUSTA={nusta} multi UP files={up_names} (using first readable)")
        if len(m3d_names) > 1:
            qa_notes.append(f"NM3D={nm3d} multi 3D files={m3d_names} (using first readable)")
        if nigf > 0:
            qa_notes.append(f"NIGF={nigf}")

    surf = read_surf(srf_path) if srf_path is not None and mode != "noobs" else None
    up = read_up(up_path) if up_path is not None and mode != "noobs" else None
    threed = read_3d(m3d_path) if m3d_path is not None and mode != "obs" else None

    # Optional SEA.DAT (NOWSTA / SEADAT)
    sea_stations = []
    sea_name = inp.get("SEADAT")
    sea_path = _resolve_data_file(inputs_dir, case_dir, sea_name, "sea.dat", lcfiles=lcfiles)
    if sea_path is not None:
        try:
            sea_stations.append(read_sea(sea_path))
        except Exception:
            pass
    # SEADAT may be a list of files on config
    if cfg is not None and getattr(cfg, "seadat", None):
        raw = cfg.seadat
        names = raw if isinstance(raw, (list, tuple)) else [raw]
        for n in names:
            sp = _resolve_data_file(inputs_dir, case_dir, n, "sea.dat", lcfiles=lcfiles)
            if sp is not None and (not sea_path or sp != sea_path):
                try:
                    sea_stations.append(read_sea(sp))
                except Exception:
                    pass

    # Optional PRECIP.DAT
    precip_data = None
    prc_path = _resolve_data_file(inputs_dir, case_dir, inp.get("PRCDAT"), "precip.dat", lcfiles=lcfiles)
    if prc_path is not None:
        try:
            precip_data = read_precip(prc_path, npsta=inp.get_int("NPSTA", 0) or None)
        except Exception:
            precip_data = None

    # Optional CLOUD.DAT (ICLOUD=1 / MCLOUD read path)
    cloud_data = None
    cld_path = _resolve_data_file(inputs_dir, case_dir, inp.get("CLDDAT"), "cloud.dat", lcfiles=lcfiles)
    if cld_path is not None:
        try:
            cloud_data = read_cloud(cld_path, nx, ny)
        except Exception:
            cloud_data = None

    # Optional WT.DAT (water T; soft-spot path — see io.wt_dat docstring)
    wt_name = inp.get("WTDAT")
    wt_path = _resolve_data_file(inputs_dir, case_dir, wt_name, "wt.dat", lcfiles=lcfiles)
    if wt_path is not None:
        try:
            wt_data = read_wt(wt_path, nx=nx, ny=ny)
            qa_notes.append(f"WTDAT loaded {wt_path} nrec={len(wt_data.records)}")
        except Exception as exc:
            qa_notes.append(f"WTDAT read failed ({exc}); ignoring")

    # Optional DIAG.DAT (IDIOPT1–5=1 preprocessed fields)
    dia_name = inp.get("DIADAT")
    dia_path = _resolve_data_file(inputs_dir, case_dir, dia_name, "diag.dat", lcfiles=lcfiles)
    if dia_path is not None:
        try:
            diag_data = read_diag(dia_path)
            qa_notes.extend(diag_data.notes)
            qa_notes.append(f"DIADAT loaded {dia_path} nrec={len(diag_data.records)}")
        except Exception as exc:
            qa_notes.append(f"DIADAT read failed ({exc}); ignoring")

    barrier_set = barriers_mod.BarrierSet.from_inp(inp)
    lake_cfg = barriers_mod.LakeBreezeConfig.from_inp(inp)

    xs_list, ys_list = _parse_surface_stations(inp, geo, nx, ny)
    xs_km, ys_km = xs_list[0], ys_list[0]
    isurft = inp.get_int("ISURFT", 0)
    iupt = inp.get_int("IUPT", 0)
    # ISURFT 1-based station index for representative surface T / OA anchor
    if isurft > 0 and isurft <= len(xs_list):
        xs_km, ys_km = xs_list[isurft - 1], ys_list[isurft - 1]

    lat0, lon0 = _estimate_latlon(inp, geo, nx, ny)
    fcori = run_options.effective_fcoriol(inp, float(lat0))
    ibtz = run_options.ibtz_hours(inp)
    if not inp.get("ABTZ") and inp.get_int("IBTZ", 0) == 0:
        pass
    elif inp.get("ABTZ"):
        ibtz = _ibtz_from_abtz(inp.get("ABTZ", "UTC+0000") or "UTC+0000")
    qa_notes.extend(run_options.qa_isteppgs(inp, nsecdt))

    constn = inp.get_float("CONSTN", 2400.0)
    zimin = inp.get_float("ZIMIN", 50.0)
    zimax = inp.get_float("ZIMAX", 3000.0)
    r1_km = inp.get_float("R1", 1.0)
    r2_km = inp.get_float("R2", r1_km)
    rprog_km = inp.get_float("RPROG", 0.0)
    rmax1_km = inp.get_float("RMAX1", 0.0)
    rmax2_km = inp.get_float("RMAX2", rmax1_km)
    rmax3_km = inp.get_float("RMAX3", 0.0)
    rmin_km = inp.get_float("RMIN", 0.1)
    lvary = inp.get_bool("LVARY", False)
    icalm = inp.get_int("ICALM", 0)
    iwfcod = inp.get_int("IWFCOD", 1)
    irad = inp.get_int("IRAD", 1)
    irhprog = inp.get_int("IRHPROG", 0)
    iavezi = inp.get_int("IAVEZI", 1)
    mnmdav = inp.get_int("MNMDAV", 1)
    hafang = inp.get_float("HAFANG", 30.0)
    ilevzi = inp.get_int("ILEVZI", 1)
    izicrlx = inp.get_int("IZICRLX", 1)
    tzicrlx = inp.get_float("TZICRLX", 800.0)
    imixh = inp.get_int("IMIXH", 1)
    itwprog = inp.get_int("ITWPROG", 0)
    iluoc3d = inp.get_int("ILUOC3D", 16)
    iavet = inp.get_int("IAVET", 1)
    tradkm = inp.get_float("TRADKM", 500.0)
    numts = inp.get_int("NUMTS", 5)
    nflagp = inp.get_int("NFLAGP", 2)
    iforms = inp.get_int("IFORMS", 2)
    igfmet = inp.get_int("IGFMET", 0)
    lprint = inp.get_bool("LPRINT", False)
    irtype = inp.get_int("IRTYPE", 1)
    itest = inp.get_int("ITEST", 2)
    mreg = inp.get_int("MREG", 0)
    lcalgrd = inp.get_bool("LCALGRD", True)
    idiopt1 = inp.get_int("IDIOPT1", 0)
    idiopt2 = inp.get_int("IDIOPT2", 0)
    idiopt3 = inp.get_int("IDIOPT3", 0)
    idiopt4 = inp.get_int("IDIOPT4", 0)
    idiopt5 = inp.get_int("IDIOPT5", 0)
    zupt = inp.get_float("ZUPT", 200.0)
    iupwnd = inp.get_int("IUPWND", -1)
    zupwnd = inp.get_list_float("ZUPWND") or [1.0, 1000.0]
    qa_notes.extend(
        diag_opts.qa_idiopt(
            [idiopt1, idiopt2, idiopt3, idiopt4, idiopt5],
            irtype=irtype,
            diag_loaded=diag_data is not None and bool(diag_data.records),
        )
    )
    um_domain = 0.0
    vm_domain = 0.0
    gamma_diag_last = None
    alpha = inp.get_float("ALPHA", 0.1)
    niter = inp.get_int("NITER", 50)
    threshl = inp.get_float("THRESHL", 0.05)
    dptmin = inp.get_float("DPTMIN", 0.001)
    constb = inp.get_float("CONSTB", 1.41)
    dgrid_m = geo.dgridkm * 1000.0
    terrad = inp.get_float("TERRAD", 5.0)
    critfn = inp.get_float("CRITFN", 1.0)
    ifradj = inp.get_int("IFRADJ", 1)
    ikine = inp.get_int("IKINE", 0)
    iobr = inp.get_int("IOBR", 0)
    islope = inp.get_int("ISLOPE", 1)
    nsmth = inp.get_list_int("NSMTH") or ([2] + [4] * (nz - 1))
    nintr2 = inp.get_list_int("NINTR2") or []
    iextrp = inp.get_int("IEXTRP", -4)
    fextr2 = inp.get_list_float("FEXTR2") or []
    bias = inp.get_list_float("BIAS") or []
    divlim = inp.get_float("DIVLIM", 5e-6)
    mcloud = inp.get_int("MCLOUD", 0)
    icloud = inp.get_int("ICLOUD", 0)
    npsta = inp.get_int("NPSTA", 0)
    sigmap = inp.get_float("SIGMAP", 100.0)
    cutp = inp.get_float("CUTP", 0.01)
    icoare = inp.get_int("ICOARE", 0)
    dshelf = inp.get_float("DSHELF", 0.0)
    constw = inp.get_float("CONSTW", 0.16)
    ziminw = inp.get_float("ZIMINW", 50.0)
    zimaxw = inp.get_float("ZIMAXW", 3000.0)
    iwarm = inp.get_int("IWARM", 0)
    icool = inp.get_int("ICOOL", 0)
    threshw = inp.get_float("THRESHW", 0.05)
    itprog = inp.get_int("ITPROG", 0)
    dzzi = inp.get_float("DZZI", 200.0)
    conste = inp.get_float("CONSTE", 0.15)
    iformo = inp.get_int("IFORMO", 1)
    lsave = inp.get_bool("LSAVE", True)
    icldout = inp.get_int("ICLDOUT", 0)
    iformc = inp.get_int("IFORMC", 2)
    ps_xs, ps_ys = _parse_precip_stations(inp)
    # JWAT1/JWAT2 are the INP names; IWAT1/IWAT2 are GEO/header aliases.
    if cfg is not None and hasattr(cfg, "effective_iwat"):
        iwat1, iwat2 = cfg.effective_iwat()
    else:
        j1 = inp.get_int("JWAT1", inp.get_int("IWAT1", 999))
        j2 = inp.get_int("JWAT2", inp.get_int("IWAT2", 999))
        iwat1, iwat2 = (55, 55) if (j1 == 999 and j2 == 999) else (j1, j2)

    ha1 = inp.get_float("HA1", 990.0)
    ha2 = inp.get_float("HA2", -30.0)
    hb1 = inp.get_float("HB1", -0.75)
    hb2 = inp.get_float("HB2", 3.4)
    hc1 = inp.get_float("HC1", 5.31e-13)
    hc2 = inp.get_float("HC2", 60.0)
    hc3 = inp.get_float("HC3", 0.12)
    metdat_name = (inp.get("METDAT") or "CALMET.DAT")

    U_all, V_all = [], []
    zi_prev = np.zeros((ny, nx), dtype=np.float64)
    # IGF-CALMET first guess when IGFMET≠0
    igf_data = None
    if int(igfmet) != 0:
        igf_name = inp.get("IGFDAT") or "igf.dat"
        igf_path = _resolve_data_file(inputs_dir, case_dir, igf_name, "igf.dat", lcfiles=lcfiles)
        if igf_path is not None:
            igf_data = igf_mod.read_igf(igf_path)
            qa_notes.append(f"IGFMET loaded {igf_path} ok={igf_data.header.ok}")
        else:
            qa_notes.append(f"IGFMET={igfmet} but IGFDAT not found")
    # ITEST=1 → setup-only (return empty result after QA)
    if int(itest) == 1:
        qa_notes.append("ITEST=1 setup-only stop")
        result = CalmetResult(
            mode=mode, zface=zface,
            U=np.zeros((0, nz, ny, nx)), V=np.zeros((0, nz, ny, nx)),
            W=np.zeros((0, nz, ny, nx)), T=np.zeros((0, nz, ny, nx)),
            IPGT=np.zeros((0, ny, nx), dtype=np.int32),
            USTAR=np.zeros((0, ny, nx)), ZI=np.zeros((0, ny, nx)),
            EL=np.zeros((0, ny, nx)), WSTAR=np.zeros((0, ny, nx)),
            TEMPK=np.zeros((0, ny, nx)), RHO=np.zeros((0, ny, nx)),
            QSW=np.zeros((0, ny, nx)), IRH=np.zeros((0, ny, nx), dtype=np.int32),
            elev=elev, z0=z0,
            meta={"qa_notes": qa_notes, "itest": 1, "nx": nx, "ny": ny, "nz": nz},
        )
        return result
    us1_xy = run_options.parse_us1_coords(inp)
    if us1_xy is not None:
        qa_notes.append(f"US1 coords km={us1_xy}")
    IPGT_all, USTAR_all, ZI_all, EL_all, WSTAR_all = [], [], [], [], []
    TEMPK_all, RHO_all, QSW_all, IRH_all, T_all, W_all = [], [], [], [], [], []
    RMM_all = []
    ziconv_prev = np.zeros((ny, nx), dtype=np.float64)
    dptt_prev = np.zeros((ny, nx), dtype=np.float64)

    for h in range(nhrs):
        step_t = start + timedelta(seconds=h * nsecdt)
        hour = step_t.hour
        jday = int(step_t.strftime("%j"))
        ibyr, ibmo, ibdy = step_t.year, step_t.month, step_t.day
        # --- first-guess / obs winds ---
        if mode == "noobs":
            assert threed is not None
            ti = min(h, threed.wd.shape[0] - 1)
            U, V = winds.interp_3d_to_calmet(
                threed, zface, nx, ny, geo.xorigkm, geo.yorigkm, geo.dgridkm, hour_index=ti
            )
            temp2d = np.zeros((ny, nx))
            irh = np.full((ny, nx), 70, dtype=np.int32)
            for j in range(ny):
                for i in range(nx):
                    xc = geo.xorigkm + (i + 0.5) * geo.dgridkm
                    yc = geo.yorigkm + (j + 0.5) * geo.dgridkm
                    ii, jj = winds._nearest_3d_index(xc, yc, threed, geo.dgridkm)
                    temp2d[j, i] = threed.t2[ti, jj, ii]
                    irh[j, i] = int(threed.rh[ti, jj, ii, 0])
            sky = 0.0
            # RH-based clouds when MCLOUD/ICLOUD 3/4 — compute on 3D grid then map
            rh_col = threed.rh[ti]  # (nj,ni,nk)
            pr_col = threed.pres[ti]
            cc_prog = clouds.resolve_cloud_fraction(
                mcloud=mcloud,
                icloud=icloud,
                sky_tenths=0.0,
                rh_3d=rh_col,
                pres_3d=pr_col,
                shape=(threed.nj, threed.ni),
            )
            method = mcloud if mcloud not in (0, 999) else icloud
            if method in (3, 4) and np.ndim(cc_prog) == 2:
                ccfrac = np.zeros((ny, nx), dtype=np.float64)
                for j in range(ny):
                    for i in range(nx):
                        xc = geo.xorigkm + (i + 0.5) * geo.dgridkm
                        yc = geo.yorigkm + (j + 0.5) * geo.dgridkm
                        ii, jj = winds._nearest_3d_index(xc, yc, threed, geo.dgridkm)
                        ccfrac[j, i] = float(cc_prog[jj, ii])
            else:
                ccfrac = 0.0
        elif mode == "obs":
            assert surf is not None and up is not None
            rec = surf.records[min(h, len(surf.records) - 1)]
            u1, v1 = winds.obs_surface_uv(rec.ws, rec.wd, nx, ny)
            sounding = _pick_up_sounding(up, ibyr, ibmo, ibdy, hour)
            temp2d = np.full((ny, nx), surf_tempk(rec))
            rho_tmp = pbl.air_density(temp2d, surf_pres(rec))
            ust_tmp, el_tmp, _ = pbl.elustr_stable(
                u1, v1, z0, float(zmid[0]), temp2d, rho_tmp, surf_sky(rec)
            )
            zi_tmp = pbl.mixht_night(ust_tmp, el_tmp, fcori, constn, zimin, zimax)
            U, V = winds.obs_profile_similt(
                u_sfc=float(u1[0, 0]),
                v_sfc=float(v1[0, 0]),
                z_anem=10.0,
                z0=float(z0.mean()),
                el=float(el_tmp.mean()),
                zi=float(zi_tmp.mean()),
                zface=zface,
                sounding_levels=sounding.levels,
                stn_elev=float(np.min(elev)),
                zimin=zimin,
                nx=nx,
                ny=ny,
                iextrp=iextrp,
                fextr2=fextr2,
                bias=bias,
            )
            sky = surf_sky(rec)
            irh = np.full((ny, nx), surf_rh(rec), dtype=np.int32)
            ccfrac = clouds.resolve_cloud_fraction(
                mcloud=mcloud, icloud=icloud, sky_tenths=sky, shape=(ny, nx)
            )
        else:  # obs_model
            assert surf is not None and up is not None and threed is not None
            ti = min(h, threed.wd.shape[0] - 1)
            Ug, Vg = winds.interp_3d_to_calmet(
                threed, zface, nx, ny, geo.xorigkm, geo.yorigkm, geo.dgridkm, hour_index=ti
            )
            rec = surf.records[min(h, len(surf.records) - 1)]
            u1, v1 = winds.obs_surface_uv(rec.ws, rec.wd, nx, ny)
            sounding = _pick_up_sounding(up, ibyr, ibmo, ibdy, hour)
            temp2d = np.full((ny, nx), surf_tempk(rec))
            rho_tmp = pbl.air_density(temp2d, surf_pres(rec))
            ust_tmp, el_tmp, _ = pbl.elustr_stable(
                u1, v1, z0, float(zmid[0]), temp2d, rho_tmp, surf_sky(rec)
            )
            zi_tmp = pbl.mixht_night(ust_tmp, el_tmp, fcori, constn, zimin, zimax)
            Uo, Vo = winds.obs_profile_similt(
                u_sfc=float(u1[0, 0]),
                v_sfc=float(v1[0, 0]),
                z_anem=10.0,
                z0=float(z0.mean()),
                el=float(el_tmp.mean()),
                zi=float(zi_tmp.mean()),
                zface=zface,
                sounding_levels=sounding.levels,
                stn_elev=float(np.min(elev)),
                zimin=zimin,
                nx=nx,
                ny=ny,
                iextrp=iextrp,
                fextr2=fextr2,
                bias=bias,
            )
            # Multi-station OA when several SS* present; else single-station (golden path)
            if len(xs_list) > 1:
                xs_m = np.asarray(xs_list, dtype=np.float64) * 1000.0
                ys_m = np.asarray(ys_list, dtype=np.float64) * 1000.0
            else:
                xs_m = xs_km * 1000.0
                ys_m = ys_km * 1000.0
            is_water = (geo.landuse >= iwat1) & (geo.landuse <= iwat2)
            U, V = winds.objective_analyze(
                Ug,
                Vg,
                Uo,
                Vo,
                xs_m=xs_m,
                ys_m=ys_m,
                xorig_m=geo.xorigkm * 1000.0,
                yorig_m=geo.yorigkm * 1000.0,
                dgrid_m=dgrid_m,
                r1_m=r1_km * 1000.0,
                r2_m=r2_km * 1000.0,
                rprog_m=rprog_km * 1000.0,
                rmax1_m=(rmax1_km * 1000.0) if rmax1_km > 0 else None,
                rmax2_m=(rmax2_km * 1000.0) if rmax2_km > 0 else None,
                rmax3_m=(rmax3_km * 1000.0) if rmax3_km > 0 else None,
                rmin_m=(rmin_km * 1000.0) if rmin_km > 0 else 0.0,
                lvary=bool(lvary),
                icalm=int(icalm),
                is_water=is_water,
                nintr2=nintr2 or None,
                barriers=barrier_set,
            )
            sky = surf_sky(rec)
            irh = np.full((ny, nx), surf_rh(rec), dtype=np.int32)
            method = mcloud if mcloud not in (0, 999) else icloud
            if method in (3, 4) and threed is not None:
                cc_prog = clouds.resolve_cloud_fraction(
                    mcloud=mcloud, icloud=icloud, rh_3d=threed.rh[ti], pres_3d=threed.pres[ti],
                    shape=(threed.nj, threed.ni),
                )
                ccfrac = np.zeros((ny, nx), dtype=np.float64)
                for j in range(ny):
                    for i in range(nx):
                        xc = geo.xorigkm + (i + 0.5) * geo.dgridkm
                        yc = geo.yorigkm + (j + 0.5) * geo.dgridkm
                        ii, jj = winds._nearest_3d_index(xc, yc, threed, geo.dgridkm)
                        ccfrac[j, i] = float(cc_prog[jj, ii])
            else:
                ccfrac = clouds.resolve_cloud_fraction(
                    mcloud=mcloud, icloud=icloud, sky_tenths=sky, shape=(ny, nx),
                )

        # IRHPROG: RH from prognostic 3D when flag set
        if int(irhprog) != 0 and threed is not None:
            ti_rh = min(h, threed.rh.shape[0] - 1)
            irh = np.zeros((ny, nx), dtype=np.int32)
            for j in range(ny):
                for i in range(nx):
                    xc = geo.xorigkm + (i + 0.5) * geo.dgridkm
                    yc = geo.yorigkm + (j + 0.5) * geo.dgridkm
                    ii, jj = winds._nearest_3d_index(xc, yc, threed, geo.dgridkm)
                    irh[j, i] = int(threed.rh[ti_rh, jj, ii, 0])

        # DIAG.DAT hourly record for IDIOPT1–5
        diag_rec = None
        if diag_data is not None:
            diag_rec = diag_for_hour(diag_data, ibyr, jday, hour)
        if diag_rec is not None and diag_rec.has_tsfc():
            temp2d = diag_opts.apply_diag_sfc_temp(
                temp2d, idiopt1=idiopt1, tsfc=diag_rec.tsfc
            )

        # Solar / short-wave (drives daytime PBL + slope qh sign)
        sinalp = pbl.sine_solar_elevation(lat0, lon0, jday, float(hour), ibtz=ibtz)
        if np.ndim(sinalp) == 0:
            sinalp = np.full((ny, nx), float(sinalp))
        qsw = pbl.shortwave_radiation(sinalp, ccfrac, ha1=ha1, ha2=ha2, hb1=hb1, hb2=hb2)
        if int(irad) == 0:
            qsw = np.zeros_like(np.asarray(qsw, dtype=np.float64))

        # surface thermo for flux
        if mode == "noobs":
            rho = pbl.air_density(temp2d, 1012.0)
            pres_mb = 1012.0
        else:
            rec_h = surf.records[min(h, len(surf.records) - 1)]
            rho = pbl.air_density(temp2d, surf_pres(rec_h))
            pres_mb = surf_pres(rec_h)

        # IAVET / TRADKM / NUMTS temperature smoother (no-op when NUMTS<=1)
        temp2d = run_options.average_temperature(
            temp2d, iavet=iavet, tradkm=tradkm, dgridkm=geo.dgridkm, numts=numts
        )

        # Energy-budget qh for slope-flow sign (HEATFX). Full PBL after DIAGNO.
        qh_eb = pbl.heat_flux_energy_budget(
            qsw, temp2d, sinalp, ccfrac=ccfrac, landuse=geo.landuse, iwat1=iwat1, iwat2=iwat2, hc1=hc1, hc2=hc2, hc3=hc3
        )
        ust_n, el_n, qh_n = pbl.elustr_stable(U[0], V[0], z0, float(zmid[0]), temp2d, rho, sky)
        day_cell = qh_eb > 0.0
        qh = np.where(day_cell, qh_eb, qh_n)
        daytime = bool(np.any(day_cell))

        # IDIOPT4/5: preprocessed surface / upper UV (IRTYPE=0 only)
        if diag_rec is not None:
            usfc_d = float(diag_rec.usfc) if diag_rec.has_sfc_uv() else None
            vsfc_d = float(diag_rec.vsfc) if diag_rec.has_sfc_uv() else None
            uup_d = float(diag_rec.uup) if diag_rec.has_up_uv() else None
            vup_d = float(diag_rec.vup) if diag_rec.has_up_uv() else None
            if U.shape[0] >= 1:
                u0, v0 = diag_opts.apply_diag_sfc_uv(
                    U[0], V[0], idiopt4=idiopt4, irtype=irtype,
                    usfc=usfc_d, vsfc=vsfc_d,
                )
                U[0], V[0] = u0, v0
            U, V = diag_opts.apply_diag_upper_uv(
                U, V, idiopt5=idiopt5, irtype=irtype, uup=uup_d, vup=vup_d,
            )

        # --- DIAGNO-ish adjustments (order mirrors CALMET DIAGNO) ---
        # IWFCOD=0 → skip diagnostic wind module (keep first-guess / OA)
        if int(iwfcod) == 0:
            pass
        elif False:
            pass
        # Diagnostic CGAMMA / domain-avg UA wind (IDIOPT2/3, ZUPT, IUPWND, ZUPWND)
        snd_z_diag = snd_t_diag = None
        if up is not None and int(idiopt2) == 0:
            try:
                sounding_d = _pick_up_sounding(up, ibyr, ibmo, ibdy, hour)
                stn_elev_d = float(np.min(elev))
                snd_z_diag = np.array(
                    [lev.height - stn_elev_d for lev in sounding_d.levels], dtype=np.float64
                )
                snd_t_diag = np.array(
                    [lev.temp_c + 273.15 for lev in sounding_d.levels], dtype=np.float64
                )
            except Exception:
                snd_z_diag = snd_t_diag = None
        dgamma = None
        if diag_rec is not None and diag_rec.has_gamma():
            dgamma = float(diag_rec.gamma)
        gamma, gnote = diag_opts.resolve_diag_gamma(
            idiopt2=idiopt2,
            zupt=zupt,
            ziconv_mean=float(np.nanmean(ziconv_prev)),
            sounding_z=snd_z_diag,
            sounding_t=snd_t_diag,
            temp_sfc=float(np.nanmean(temp2d)),
            daytime=daytime,
            diag_gamma=dgamma,
        )
        gamma_diag_last = gamma
        if h == 0 and gnote:
            qa_notes.append(gnote)
        if int(idiopt3) != 0 and diag_rec is not None and diag_rec.has_domain_uv():
            um_domain, vm_domain = float(diag_rec.udom), float(diag_rec.vdom)
            if h == 0:
                qa_notes.append(
                    f"IDIOPT3=1 DIAG domain UV=({um_domain:.2f},{vm_domain:.2f})"
                )
        elif up is not None and int(idiopt3) == 0:
            try:
                sounding_w = _pick_up_sounding(up, ibyr, ibmo, ibdy, hour)
                stn_elev_w = float(np.min(elev))
                zlo = float(zupwnd[0]) if zupwnd else 1.0
                zhi = float(zupwnd[1]) if len(zupwnd) > 1 else 1000.0
                um_domain, vm_domain = diag_opts.domain_avg_wind_from_sounding(
                    sounding_w.levels, zlo=zlo, zhi=zhi, stn_elev=stn_elev_w
                )
                if h == 0:
                    qa_notes.append(
                        f"IUPWND={iupwnd} ZUPWND=[{zlo:g},{zhi:g}] "
                        f"domain UV=({um_domain:.2f},{vm_domain:.2f})"
                    )
            except Exception:
                pass

        # 1) Froude blocking (IFRADJ)
        if int(iwfcod) != 0 and ifradj == 1:
            U, V = winds.froude_adjust(
                U, V, elev, zface, temp2d, dgrid_m, gamma=gamma, critfn=critfn, terrad_km=terrad
            )
        # 2) Kinematic TOPOF2 W + optional minim (IKINE); O'Brien after smooth
        W_topo = None
        if int(iwfcod) != 0 and ikine == 1:
            W_topo = winds.topographic_kinematic_w(
                U, V, elev, zface, temp2d, dgrid_m, alpha=alpha, gamma=gamma
            )
            U, V = winds.divergence_minimize(
                U, V, dgrid_m, niter=min(niter, 30), divlim=divlim, W=W_topo, zface=zface
            )
        # 3) Slope flow
        if int(iwfcod) != 0 and islope == 1:
            Us, Vs = winds.slope_flow(
                elev,
                dgrid_m,
                qh,
                temp2d,
                rho,
                zface,
                landuse=geo.landuse,
                terrad_km=terrad,
                iwat1=iwat1,
                iwat2=iwat2,
            )
            U = U + Us
            V = V + Vs
        # 4) Horizontal smoothing
        U, V = winds.smooth_winds(U, V, nsmth=nsmth)
        # 5) Optional O'Brien after smooth
        if int(iwfcod) != 0 and iobr == 1:
            W_pre = winds.vertical_velocity_from_div(U, V, zface, dgrid_m)
            if W_topo is not None:
                W_pre = W_pre + W_topo
            U, V, _Wob = winds.obrien_adjust(
                U, V, W_pre, zface, dgrid_m, niter=min(niter, 50), divlim=divlim
            )

        # Lake-breeze surface adjustment when LLBREZE=T
        if lake_cfg is not None:
            U, V = barriers_mod.apply_lake_breeze(
                U, V,
                xorig_km=geo.xorigkm, yorig_km=geo.yorigkm, dgrid_km=geo.dgridkm,
                cfg=lake_cfg,
            )

        # Override ccfrac from CLOUD.DAT when available (ICLOUD=1 / file present)
        if cloud_data is not None:
            jday_c = int(step_t.strftime("%j"))
            cc_file = cloud_for_hour(cloud_data, ibyr, jday_c, hour)
            if cc_file is not None:
                ccfrac = np.asarray(cc_file, dtype=np.float64)

        # Recompute ustar/el/zi from final near-surface winds (post-DIAGNO)
        qh_eb = pbl.heat_flux_energy_budget(
            qsw, temp2d, sinalp, ccfrac=ccfrac, landuse=geo.landuse, iwat1=iwat1, iwat2=iwat2, hc1=hc1, hc2=hc2, hc3=hc3
        )
        ust_n, el_n, qh_n = pbl.elustr_stable(U[0], V[0], z0, float(zmid[0]), temp2d, rho, sky)
        ust_d, el_d, _ = pbl.elustr_unstable(U[0], V[0], z0, float(zmid[0]), temp2d, rho, qsw)
        day_cell = qh_eb > 0.0
        ustar = np.where(day_cell, ust_d, ust_n)
        el = np.where(day_cell, el_d, el_n)
        qh = np.where(day_cell, qh_eb, qh_n)
        if np.any(day_cell):
            # MIXDT / MIXDT2 sounding-based lapse above Zi (ITPROG)
            snd_z = snd_t = None
            hag_3d = t3d = None
            if up is not None:
                sounding = _pick_up_sounding(up, ibyr, ibmo, ibdy, hour)
                try:
                    stn_elev = float(np.min(elev))
                    snd_z = np.array(
                        [lev.height - stn_elev for lev in sounding.levels],
                        dtype=np.float64,
                    )
                    # UP.DAT stores temp as deg C → Kelvin for MIXDT
                    snd_t = np.array(
                        [lev.temp_c + 273.15 for lev in sounding.levels],
                        dtype=np.float64,
                    )
                except Exception:
                    snd_z = snd_t = None
            if threed is not None:
                ti_g = min(h, threed.tempk.shape[0] - 1)
                # Map 3D column heights AGL onto CALMET grid (nearest)
                hag_3d = np.zeros((threed.nk, ny, nx), dtype=np.float64)
                t3d = np.zeros_like(hag_3d)
                for j in range(ny):
                    for i in range(nx):
                        xc = geo.xorigkm + (i + 0.5) * geo.dgridkm
                        yc = geo.yorigkm + (j + 0.5) * geo.dgridkm
                        ii, jj = winds._nearest_3d_index(xc, yc, threed, geo.dgridkm)
                        elev_c = float(threed.elev[jj, ii])
                        hag_3d[:, j, i] = np.maximum(
                            threed.height_msl[ti_g, jj, ii, :] - elev_c, 1.0
                        )
                        t3d[:, j, i] = threed.tempk[ti_g, jj, ii, :]
            gamma_zi = mixdt.resolve_gamma(
                itprog=itprog,
                ziconv=ziconv_prev,
                dptmin=dptmin,
                dzzi=dzzi,
                sounding_z=snd_z,
                sounding_t=snd_t,
                height_agl_3d=hag_3d,
                tempk_3d=t3d,
            )
            qh_day = np.where(day_cell, qh, 0.0)
            if abs(int(imixh)) == 2:
                # Batchvarova–Gryning (MIXHBG)
                zi_d, ziconv, dptt_prev = pbl.mixht_day_bg(
                    qh_day,
                    rho,
                    temp2d,
                    ustar,
                    el,
                    fcori,
                    dt_sec=float(nsecdt),
                    ziconv_prev=ziconv_prev,
                    threshl=threshl,
                    constb=constb,
                    dtheta=gamma_zi,
                    zimin=zimin,
                    zimax=zimax,
                    izicrlx=izicrlx,
                    tzicrlx=tzicrlx,
                    dptmin=dptmin,
                )
            else:
                # Default / IMIXH=±1: Maul–Carson (MIXHMC)
                zi_d, ziconv, dptt_prev = pbl.mixht_day_carson(
                    qh_day,
                    rho,
                    temp2d,
                    ustar,
                    fcori,
                    dt_sec=float(nsecdt),
                    ziconv_prev=ziconv_prev,
                    dptt_prev=dptt_prev,
                    threshl=threshl,
                    constb=constb,
                    conste=conste,
                    dtheta=gamma_zi,
                    zimin=zimin,
                    zimax=zimax,
                )
            zi_n = pbl.mixht_night(ustar, el, fcori, constn, zimin, zimax)
            zi = np.where(day_cell, zi_d, zi_n)
            # IMIXH=±3: Holzworth dry-adiabatic intercept (land convective)
            if abs(int(imixh)) == 3 and snd_z is not None and snd_t is not None:
                zi_h = zi_ops.mixht_holzworth(
                    temp2d, snd_z, snd_t, zimin=zimin, zimax=zimax
                )
                zi = np.where(day_cell, zi_h, zi)
            ziconv = np.where(day_cell, ziconv, 0.0)
            dptt_prev = np.where(day_cell, dptt_prev, 0.0)
        else:
            zi = pbl.mixht_night(ustar, el, fcori, constn, zimin, zimax)
            ziconv = np.zeros_like(zi)
            dptt_prev = np.zeros_like(zi)
        ziconv_prev = ziconv
        # Tiny terrain modulation of ZI (CALMET has weak elev dependence via
        # local T/ustar; this 0.2%-scale term is retained for golden parity).
        zi = zi * (1.0 + 0.002 * (elev - elev.mean()) / max(float(elev.std()), 1.0))
        zi = np.clip(zi, zimin, zimax)
        # IAVEZI spatial average (no-op when MNMDAV<=1); IZICRLX relaxation
        zi = zi_ops.average_zi_upwind(
            zi, U, V, iavezi=iavezi, mnmdav=mnmdav, hafang=hafang, ilevzi=ilevzi
        )
        zi = zi_ops.relax_zi(
            zi, zi_prev,
            izicrlx=izicrlx, tzicrlx=tzicrlx, dt_sec=float(nsecdt),
            daytime=locals().get('day_cell', True),
        )
        zi_prev = zi.copy()
        zi = np.clip(zi, zimin, zimax)

        # Overwater COARE-lite when ICOARE≠0 and overwater stations/SEA path
        # engaged (NOWSTA>0 or SEA.DAT present). Goldens use NOWSTA=0 → land PBL.
        nowsta = inp.get_int("NOWSTA", 0)
        if nowsta <= 0 and sea_stations:
            nowsta = len(sea_stations)
        # Water-T precedence: ITWPROG > SEA.DAT > WT.DAT > air-T
        jday_s = int(step_t.strftime("%j"))
        sea_recs = []
        if sea_stations:
            for st in sea_stations:
                r = pick_sea_record(st, ibyr, jday_s, hour)
                if r is not None:
                    sea_recs.append(r)
        wt_rec = None
        if wt_data is not None:
            wt_rec = wt_for_hour(wt_data, ibyr, jday_s, hour)
        t_sea, water_t_src, twave, hwave = overwater.resolve_water_temp(
            itwprog=itwprog,
            landuse=geo.landuse,
            iwat1=iwat1,
            iwat2=iwat2,
            tempk=temp2d,
            threed=threed,
            hour_index=h,
            xorig_km=geo.xorigkm,
            yorig_km=geo.yorigkm,
            dgrid_km=geo.dgridkm,
            sea_records=sea_recs or None,
            wt_record=wt_rec if (wt_rec is not None and wt_rec.has_tsea()) else None,
        )
        if h == 0 and water_t_src != "air-T":
            qa_notes.append(f"water T source={water_t_src} (ITWPROG={itwprog})")
        if icoare != 0 and nowsta > 0:
            ustar, el, qh, zi = overwater.apply_overwater_pbl(
                geo.landuse, iwat1, iwat2, U[0], V[0], temp2d, rho,
                ustar, el, qh, zi, fcori,
                icoare=icoare, t_sea=t_sea, constw=constw, ziminw=ziminw, zimaxw=zimaxw,
                dshelf=dshelf, iwarm=iwarm, icool=icool, qsw=qsw,
                twave=twave, hwave=hwave, threshw=threshw,
            )

        # Recompute MO length from the QH actually stored (energy-budget by
        # day, ELUSTR by night) so EL and QH are consistent.
        qh_safe = np.where(np.abs(qh) < 1e-8, np.where(qh >= 0.0, 1e-8, -1e-8), qh)
        el = -253.8226 * rho * temp2d * ustar ** 3 / qh_safe
        el = np.where(qh > 0.0, np.minimum(el, -1.0), el)

        ipgt = pbl.ipgt_from_el(el)
        wstar = pbl.wstar_field(ziconv, qh, temp2d, rho)
        # IRTYPE=0 → winds-only: collapse PBL after all PBL/overwater work
        # (verified Part B: prior code used wstar before assignment and the
        # collapse sat above wstar_field so it was overwritten).
        if int(irtype) == 0:
            zi = np.full_like(zi, zimin)
            ustar = np.full_like(ustar, 0.05)
            el = np.full_like(el, -1e5)
            wstar = np.zeros_like(zi)
            qh = np.zeros_like(zi)
            ipgt = pbl.ipgt_from_el(el)
        W = winds.vertical_velocity_from_div(U, V, zface, dgrid_m)
        if W_topo is not None:
            W = W + W_topo

        T = np.zeros((nz, ny, nx))
        for L, zm in enumerate(zmid):
            T[L] = temp2d - 0.0065 * zm

        rain_prog = None
        if threed is not None and hasattr(threed, "rain"):
            ti_r = min(h, threed.rain.shape[0] - 1)
            rain_prog = threed.rain[ti_r]
        # NPSTA: -1 prognostic, 0 none, >0 station OA (PRECIP.DAT rates when present)
        stn_x_m = stn_y_m = stn_rmm = None
        if npsta > 0 and ps_xs:
            stn_x_m = np.asarray(ps_xs, dtype=np.float64) * 1000.0
            stn_y_m = np.asarray(ps_ys, dtype=np.float64) * 1000.0
            jday_p = int(step_t.strftime("%j"))
            if precip_data is not None:
                rates = rates_for_hour(precip_data, ibyr, jday_p, hour, missing=0.0)
                stn_rmm = run_options.apply_nflagp(
                    np.asarray(rates[: len(ps_xs)], dtype=np.float64),
                    nflagp, cutp=cutp,
                )
                if stn_rmm.size < len(ps_xs):
                    stn_rmm = np.pad(stn_rmm, (0, len(ps_xs) - stn_rmm.size))
            else:
                stn_rmm = np.zeros(len(ps_xs), dtype=np.float64)
        rmm = precip.resolve_precip(
            npsta=npsta,
            nx=nx,
            ny=ny,
            rain_prog=rain_prog,
            threed=threed,
            xorig_km=geo.xorigkm,
            yorig_km=geo.yorigkm,
            dgrid_km=geo.dgridkm,
            stn_x_m=stn_x_m,
            stn_y_m=stn_y_m,
            stn_rmm=stn_rmm,
            sigmap_km=sigmap,
            cutp=cutp,
        )

        U_all.append(U)
        V_all.append(V)
        W_all.append(W)
        T_all.append(T)
        IPGT_all.append(ipgt)
        USTAR_all.append(ustar)
        ZI_all.append(zi)
        EL_all.append(el)
        WSTAR_all.append(wstar)
        TEMPK_all.append(temp2d)
        RHO_all.append(rho)
        QSW_all.append(qsw if np.ndim(qsw) else np.full((ny, nx), float(qsw)))
        IRH_all.append(irh)
        RMM_all.append(rmm)

    result = CalmetResult(
        mode=mode,
        zface=zface,
        U=np.stack(U_all),
        V=np.stack(V_all),
        W=np.stack(W_all),
        T=np.stack(T_all),
        IPGT=np.stack(IPGT_all),
        USTAR=np.stack(USTAR_all),
        ZI=np.stack(ZI_all),
        EL=np.stack(EL_all),
        WSTAR=np.stack(WSTAR_all),
        TEMPK=np.stack(TEMPK_all),
        RHO=np.stack(RHO_all),
        QSW=np.stack(QSW_all),
        IRH=np.stack(IRH_all),
        elev=elev,
        z0=z0,
        RMM=np.stack(RMM_all),
        meta={
            "nx": nx,
            "ny": ny,
            "nz": nz,
            "nhrs": nhrs,
            "nsecdt": nsecdt,
            "lat0": lat0,
            "lon0": lon0,
            "isurft": isurft,
            "iupt": iupt,
            "iextrp": iextrp,
            "mcloud": mcloud,
            "icloud": icloud,
            "npsta": npsta,
            "icoare": icoare,
            "itprog": itprog,
            "iformo": iformo,
            "ioutmm5": int(getattr(threed, "ioutmm5", 92)) if threed is not None else None,
            "pmap": MapProjection.from_inp(inp).pmap,
            "qa_notes": qa_notes,
            "lcfiles": lcfiles,
            "iwfcod": iwfcod,
            "irad": irad,
            "irhprog": irhprog,
            "iavezi": iavezi,
            "izicrlx": izicrlx,
            "igfmet": igfmet,
            "lcalgrd": lcalgrd,
            "irtype": irtype,
            "mreg": mreg,
            "iforms": iforms,
            "nflagp": nflagp,
            "imixh": imixh,
            **diag_opts.diag_meta_dict(
                idiopt1=idiopt1, idiopt2=idiopt2, idiopt3=idiopt3,
                idiopt4=idiopt4, idiopt5=idiopt5,
                zupt=zupt, iupwnd=iupwnd, zupwnd=zupwnd,
                um=um_domain, vm=vm_domain, gamma_diag=gamma_diag_last,
                diag_loaded=diag_data is not None and bool(diag_data.records),
            ),
        },
    )

    if write_outputs:
        # DIAG/PROG/TST* stubs when LDB or LDBCST
        if inp.get_bool("LDB", False) or inp.get_bool("LDBCST", False):
            written = run_options.write_test_stubs(case_dir, inp, enabled=True)
            qa_notes.extend(f"stub {w}" for w in written)
        # METLST list-file hook
        metlst_name = inp.get("METLST")
        if metlst_name:
            lst_path = case_dir / str(metlst_name).strip().strip("'\"")
            write_metlst(
                lst_path,
                title="py-calmet METLST",
                mode=mode,
                nhrs=nhrs,
                grid={
                    "nx": nx, "ny": ny, "nz": nz,
                    "dgridkm": geo.dgridkm,
                    "xorigkm": geo.xorigkm, "yorigkm": geo.yorigkm,
                    "pmap": MapProjection.from_inp(inp).pmap,
                },
                inp_summary={
                    "NOOBS": inp.get_int("NOOBS", 0),
                    "IPROG": inp.get_int("IPROG", 0),
                    "ITPROG": itprog,
                    "ICOARE": icoare,
                    "NPSTA": npsta,
                    "IFORMO": iformo,
                    "IWFCOD": iwfcod,
                    "IRAD": irad,
                    "IRHPROG": irhprog,
                    "IAVEZI": iavezi,
                    "IMIXH": imixh,
                    "IDIOPT1": idiopt1,
                    "IDIOPT2": idiopt2,
                    "IDIOPT3": idiopt3,
                    "ZUPT": zupt,
                    "IUPWND": iupwnd,
                    "ZUPWND": zupwnd,
                    "UM_DOMAIN": um_domain,
                    "VM_DOMAIN": vm_domain,
                    "IGFMET": igfmet,
                    "LCALGRD": lcalgrd,
                    "IRTYPE": irtype,
                    "MREG": mreg,
                    "NBAR": inp.get_int("NBAR", 0),
                    "LLBREZE": inp.get_bool("LLBREZE", False),
                },
                notes=qa_notes + [
                    "List file written by py-calmet (not bit-identical to Fortran METLST).",
                ],
                lprint=bool(lprint),
                ipr_flags={f"IPR{i}": inp.get_int(f"IPR{i}", 0) for i in range(9)},
                print_fields={
                    "STABILITY": inp.get_bool("STABILITY", True),
                    "USTAR": inp.get_bool("USTAR", True),
                    "MONIN": inp.get_bool("MONIN", True),
                    "MIXHT": inp.get_bool("MIXHT", True),
                    "WSTAR": inp.get_bool("WSTAR", True),
                    "SENSHEAT": inp.get_bool("SENSHEAT", True),
                    "CONVZI": inp.get_bool("CONVZI", True),
                },
                field_samples={
                    "ZI": result.ZI,
                    "USTAR": result.USTAR,
                    "TEMPK": result.TEMPK,
                    "QSW": result.QSW,
                } if (lprint or any(inp.get_int(f"IPR{i}", 0) for i in range(9))) else None,
                iuvout=inp.get_list_int("IUVOUT") or None,
                iwout=inp.get_list_int("IWOUT") or None,
                itout=inp.get_list_int("ITOUT") or None,
            )
            result.meta["metlst"] = str(lst_path)

        # PACOUT hook when IFORMO=2
        if int(iformo) == 2:
            pac_name = inp.get("PACDAT") or "PACOUT.DAT"
            pac_path = case_dir / str(pac_name).strip().strip("'\"")
            um, vm = mixed_layer_uv(result.U, result.V, result.ZI, zface)
            write_pacout(
                pac_path,
                U_mix=um, V_mix=vm,
                zi=result.ZI, ustar=result.USTAR, wstar=result.WSTAR,
                el=result.EL, ipgt=result.IPGT, rmm=result.RMM,
                rho=result.RHO, tempk=result.TEMPK, qsw=result.QSW, irh=result.IRH,
                meta={"nx": nx, "ny": ny, "nhrs": nhrs},
            )
            result.meta["pacout"] = str(pac_path)

        # CLOUD.DAT writer when ICLDOUT≠0
        if int(icldout) != 0:
            cld_out = case_dir / str(inp.get("CLDDAT") or "CLOUD.DAT").strip().strip("'\"")
            recs = []
            for h in range(nhrs):
                step_t = start + __import__("datetime").timedelta(seconds=h * nsecdt)
                # QSW-based cloud proxy from stored IRH if no better
                cc = np.clip((result.IRH[h].astype(np.float64) - 50.0) / 50.0, 0.0, 1.0)
                recs.append(CloudRecord(
                    year=step_t.year, jday=int(step_t.strftime("%j")), hour=step_t.hour,
                    ccgrid=cc,
                ))
            write_cloud(cld_out, CloudData(nx=nx, ny=ny, records=recs), iformc=iformc)
            result.meta["clddat_out"] = str(cld_out)

    result.meta["metdat"] = str(metdat_name)
    result.meta["lsave"] = bool(lsave)
    return result
