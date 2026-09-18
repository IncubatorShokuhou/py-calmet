"""High-level run API matching obs / obs_model / noobs modes."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
import numpy as np

from ..io.geo import read_geo
from ..io.surf import read_surf
from ..io.up import read_up
from ..io.threed import read_3d
from ..io.inp import read_inp
from .met_utils import layer_mids, coriolis
from . import winds, pbl


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


def _is_daytime(qsw: np.ndarray | float, qh_hint: float | None = None) -> bool:
    q = float(np.mean(qsw)) if np.ndim(qsw) else float(qsw)
    if q > 20.0:
        return True
    if qh_hint is not None and qh_hint > 0:
        return True
    return False


def _estimate_latlon(inp, geo, nx, ny):
    """Return (lat0, lon0_east) domain-center estimates."""
    lat0 = inp.get_float("RLAT0", -999.0)
    lon0 = inp.get_float("RLON0", -999.0)
    # RLAT0 may be stored oddly; try pyproj from UTM
    try:
        from pyproj import Transformer

        zone = inp.get_int("IUTMZN", 19)
        to_ll = Transformer.from_crs(f"EPSG:{32600 + zone}", "EPSG:4326", always_xy=True)
        xc = (geo.xorigkm + 0.5 * nx * geo.dgridkm) * 1000.0
        yc = (geo.yorigkm + 0.5 * ny * geo.dgridkm) * 1000.0
        lon_e, lat = to_ll.transform(xc, yc)
        if lat0 < -90 or lat0 > 90:
            lat0 = float(lat)
        if lon0 < -180 or lon0 > 180:
            lon0 = float(lon_e)
    except Exception:
        if lat0 < -90 or lat0 > 90:
            lat0 = 44.25
        if lon0 < -180 or lon0 > 180:
            lon0 = -70.0
    return float(lat0), float(lon0)


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


def run_calmet(
    case_dir: str | Path,
    mode: Optional[str] = None,
    inputs_dir: Optional[str | Path] = None,
) -> CalmetResult:
    """Run pure-NumPy diagnostic for a case directory."""
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

    geo = read_geo(inputs_dir / "geo.dat")
    nx, ny = geo.nx, geo.ny
    zface = np.asarray(
        inp.get_list_float("ZFACE") or [0, 20, 40, 80, 160, 300, 600, 1000, 1500],
        dtype=np.float64,
    )
    nz = inp.get_int("NZ", len(zface) - 1)
    zmid = layer_mids(zface)
    z0 = _default_z0(geo.landuse)
    elev = geo.elev

    ibyr = inp.get_int("IBYR", 2020)
    ibmo = inp.get_int("IBMO", 6)
    ibdy = inp.get_int("IBDY", 15)
    ibhr = inp.get_int("IBHR", 0)
    iehr = inp.get_int("IEHR", 3)
    nsecdt = inp.get_int("NSECDT", 3600)
    span_sec = (iehr - ibhr) * 3600
    nhrs = max(1, span_sec // max(nsecdt, 1))

    # Julian day
    from datetime import datetime

    try:
        jday = int(datetime(ibyr, ibmo, ibdy).strftime("%j"))
    except Exception:
        jday = 166

    surf = (
        read_surf(inputs_dir / "surf.dat")
        if (inputs_dir / "surf.dat").exists() and mode != "noobs"
        else None
    )
    up = (
        read_up(inputs_dir / "up.dat")
        if (inputs_dir / "up.dat").exists() and mode != "noobs"
        else None
    )
    threed = (
        read_3d(inputs_dir / "3d.dat")
        if (inputs_dir / "3d.dat").exists() and mode != "obs"
        else None
    )

    xs_km = geo.xorigkm + nx * geo.dgridkm / 2.0
    ys_km = geo.yorigkm + ny * geo.dgridkm / 2.0
    if "SS1" in inp.raw:
        parts = inp.raw["SS1"].replace("'", " ").split()
        try:
            xs_km = float(parts[2])
            ys_km = float(parts[3])
        except Exception:
            pass

    lat0, lon0 = _estimate_latlon(inp, geo, nx, ny)
    fcori = coriolis(float(lat0))
    ibtz = _ibtz_from_abtz(inp.get("ABTZ", "UTC+0000") or "UTC+0000")

    constn = inp.get_float("CONSTN", 2400.0)
    zimin = inp.get_float("ZIMIN", 50.0)
    zimax = inp.get_float("ZIMAX", 3000.0)
    r1_km = inp.get_float("R1", 1.0)
    r2_km = inp.get_float("R2", r1_km)
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
    iwat1 = inp.get_int("IWAT1", 55)
    iwat2 = inp.get_int("IWAT2", 55)

    U_all, V_all = [], []
    IPGT_all, USTAR_all, ZI_all, EL_all, WSTAR_all = [], [], [], [], []
    TEMPK_all, RHO_all, QSW_all, IRH_all, T_all, W_all = [], [], [], [], [], []
    RMM_all = []
    ziconv_prev = np.zeros((ny, nx), dtype=np.float64)

    for h in range(nhrs):
        hour = ibhr + h * (nsecdt // 3600)
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
                    jj = min(j + 1, threed.nj - 1)
                    ii = min(i + 1, threed.ni - 1)
                    temp2d[j, i] = threed.t2[ti, jj, ii]
                    irh[j, i] = int(threed.rh[ti, jj, ii, 0])
            sky = 0.0
            ccfrac = 0.0
        elif mode == "obs":
            assert surf is not None and up is not None
            rec = surf.records[min(h, len(surf.records) - 1)]
            u1, v1 = winds.obs_surface_uv(rec.ws, rec.wd, nx, ny)
            sounding = _pick_up_sounding(up, ibyr, ibmo, ibdy, hour)
            temp2d = np.full((ny, nx), rec.tempk)
            rho_tmp = pbl.air_density(temp2d, rec.pres)
            ust_tmp, el_tmp, _ = pbl.elustr_stable(
                u1, v1, z0, float(zmid[0]), temp2d, rho_tmp, rec.sky
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
            )
            sky = rec.sky
            irh = np.full((ny, nx), rec.rh, dtype=np.int32)
            ccfrac = 0.1 * float(rec.sky)
        else:  # obs_model
            assert surf is not None and up is not None and threed is not None
            ti = min(h, threed.wd.shape[0] - 1)
            Ug, Vg = winds.interp_3d_to_calmet(
                threed, zface, nx, ny, geo.xorigkm, geo.yorigkm, geo.dgridkm, hour_index=ti
            )
            rec = surf.records[min(h, len(surf.records) - 1)]
            u1, v1 = winds.obs_surface_uv(rec.ws, rec.wd, nx, ny)
            sounding = _pick_up_sounding(up, ibyr, ibmo, ibdy, hour)
            temp2d = np.full((ny, nx), rec.tempk)
            rho_tmp = pbl.air_density(temp2d, rec.pres)
            ust_tmp, el_tmp, _ = pbl.elustr_stable(
                u1, v1, z0, float(zmid[0]), temp2d, rho_tmp, rec.sky
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
            )
            U, V = winds.objective_analyze(
                Ug,
                Vg,
                Uo,
                Vo,
                xs_m=xs_km * 1000.0,
                ys_m=ys_km * 1000.0,
                xorig_m=geo.xorigkm * 1000.0,
                yorig_m=geo.yorigkm * 1000.0,
                dgrid_m=dgrid_m,
                r1_m=r1_km * 1000.0,
                r2_m=r2_km * 1000.0,
            )
            sky = rec.sky
            irh = np.full((ny, nx), rec.rh, dtype=np.int32)
            ccfrac = 0.1 * float(rec.sky)

        # Solar / short-wave (drives daytime PBL + slope qh sign)
        sinalp = pbl.sine_solar_elevation(lat0, lon0, jday, float(hour), ibtz=ibtz)
        if np.ndim(sinalp) == 0:
            sinalp = np.full((ny, nx), float(sinalp))
        qsw = pbl.shortwave_radiation(sinalp, ccfrac)

        # surface thermo for flux
        if mode == "noobs":
            rho = pbl.air_density(temp2d, 1012.0)
            pres_mb = 1012.0
        else:
            rec_h = surf.records[min(h, len(surf.records) - 1)]
            rho = pbl.air_density(temp2d, rec_h.pres)
            pres_mb = rec_h.pres

        # Energy-budget qh for slope-flow sign (HEATFX). Full PBL after DIAGNO.
        qh_eb = pbl.heat_flux_energy_budget(
            qsw, temp2d, sinalp, ccfrac=ccfrac, landuse=geo.landuse, iwat1=iwat1, iwat2=iwat2
        )
        ust_n, el_n, qh_n = pbl.elustr_stable(U[0], V[0], z0, float(zmid[0]), temp2d, rho, sky)
        day_cell = qh_eb > 0.0
        qh = np.where(day_cell, qh_eb, qh_n)
        daytime = bool(np.any(day_cell))

        # --- DIAGNO-ish adjustments (order mirrors CALMET DIAGNO) ---
        # 1) Froude blocking (IFRADJ)
        if ifradj == 1:
            # stable lapse proxy ~0.01 K/m when night; weaker by day
            gamma = 0.01 if not daytime else 0.005
            U, V = winds.froude_adjust(
                U, V, elev, zface, temp2d, dgrid_m, gamma=gamma, critfn=critfn, terrad_km=terrad
            )
        # 2) Kinematic / O'Brien only if flagged
        if ikine == 1:
            U, V = winds.light_terrain_adjust(U, V, elev, dgrid_m, alpha=alpha)
            U, V = winds.divergence_minimize(U, V, dgrid_m, niter=min(niter, 30), alpha=0.5)
        # 3) Slope flow
        if islope == 1:
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
        if iobr == 1:
            U, V = winds.divergence_minimize(U, V, dgrid_m, niter=min(niter, 50), alpha=0.5)

        # Recompute ustar/el/zi from final near-surface winds (post-DIAGNO)
        qh_eb = pbl.heat_flux_energy_budget(
            qsw, temp2d, sinalp, ccfrac=ccfrac, landuse=geo.landuse, iwat1=iwat1, iwat2=iwat2
        )
        ust_n, el_n, qh_n = pbl.elustr_stable(U[0], V[0], z0, float(zmid[0]), temp2d, rho, sky)
        ust_d, el_d, _ = pbl.elustr_unstable(U[0], V[0], z0, float(zmid[0]), temp2d, rho, qsw)
        day_cell = qh_eb > 0.0
        ustar = np.where(day_cell, ust_d, ust_n)
        el = np.where(day_cell, el_d, el_n)
        qh = np.where(day_cell, qh_eb, qh_n)
        if np.any(day_cell):
            zi_d, ziconv = pbl.mixht_day_carson(
                np.where(day_cell, qh, 0.0),
                rho,
                temp2d,
                ustar,
                fcori,
                dt_sec=float(nsecdt),
                ziconv_prev=ziconv_prev,
                threshl=threshl,
                constb=constb,
                dtheta=dptmin,
                zimin=zimin,
                zimax=zimax,
            )
            zi_n = pbl.mixht_night(ustar, el, fcori, constn, zimin, zimax)
            zi = np.where(day_cell, zi_d, zi_n)
            ziconv = np.where(day_cell, ziconv, 0.0)
        else:
            zi = pbl.mixht_night(ustar, el, fcori, constn, zimin, zimax)
            ziconv = np.zeros_like(zi)
        ziconv_prev = ziconv
        zi = zi * (1.0 + 0.002 * (elev - elev.mean()) / max(float(elev.std()), 1.0))
        zi = np.clip(zi, zimin, zimax)

        ipgt = pbl.ipgt_from_el(el)
        wstar = pbl.wstar_field(ziconv, qh, temp2d, rho)
        W = winds.vertical_velocity_from_div(U, V, zface, dgrid_m)

        T = np.zeros((nz, ny, nx))
        for L, zm in enumerate(zmid):
            T[L] = temp2d - 0.0065 * zm

        rmm = np.zeros((ny, nx), dtype=np.float64)

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

    return CalmetResult(
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
        },
    )
