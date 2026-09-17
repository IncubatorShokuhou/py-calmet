"""High-level run API matching the three tiny-domain modes."""
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
from .met_utils import layer_mids, wind_uv, coriolis
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
    meta: dict = field(default_factory=dict)


def _default_z0(landuse: np.ndarray) -> np.ndarray:
    # CALMET default: LU 20 -> 0.25
    z0 = np.full(landuse.shape, 0.25, dtype=np.float64)
    z0 = np.where(landuse == 10, 1.0, z0)
    return z0


def _pick_up_sounding(up, year, month, day, hour):
    # prefer exact, else nearest previous
    best = None
    best_key = None
    target = (year, month, day, hour)
    for s in up.soundings:
        key = (s.year, s.month, s.day, s.hour)
        if key <= target and (best_key is None or key > best_key):
            best, best_key = s, key
    return best or up.soundings[0]


def run_calmet(
    case_dir: str | Path,
    mode: Optional[str] = None,
    inputs_dir: Optional[str | Path] = None,
) -> CalmetResult:
    """Run pure-NumPy diagnostic for a tiny-domain case directory.

    Parameters
    ----------
    case_dir : path containing calmet.inp (and optionally linked inputs)
    mode : optional override ('obs' | 'obs_model' | 'noobs')
    inputs_dir : shared inputs (geo/surf/up/3d); defaults to ../shared or goldens/inputs
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

    geo = read_geo(inputs_dir / "geo.dat")
    nx, ny = geo.nx, geo.ny
    zface = np.asarray(inp.get_list_float("ZFACE") or [0, 20, 40, 80, 160, 300, 600, 1000, 1500], dtype=np.float64)
    nz = inp.get_int("NZ", len(zface) - 1)
    zmid = layer_mids(zface)
    z0 = _default_z0(geo.landuse)
    elev = geo.elev

    ibyr = inp.get_int("IBYR", 2020)
    ibmo = inp.get_int("IBMO", 6)
    ibdy = inp.get_int("IBDY", 15)
    ibhr = inp.get_int("IBHR", 0)
    iehr = inp.get_int("IEHR", 3)
    nhrs = iehr - ibhr

    # optional inputs
    surf = read_surf(inputs_dir / "surf.dat") if (inputs_dir / "surf.dat").exists() and mode != "noobs" else None
    up = read_up(inputs_dir / "up.dat") if (inputs_dir / "up.dat").exists() and mode != "noobs" else None
    threed = read_3d(inputs_dir / "3d.dat") if (inputs_dir / "3d.dat").exists() and mode != "obs" else None

    # station params
    xs_km = geo.xorigkm + nx * geo.dgridkm / 2.0
    ys_km = geo.yorigkm + ny * geo.dgridkm / 2.0
    if "SS1" in inp.raw:
        # 'DEMO' id x y tz anem
        parts = inp.raw["SS1"].replace("'", " ").split()
        # name id x y ...
        nums = [p for p in parts if p.replace(".", "", 1).replace("-", "", 1).isdigit() or p.replace(".", "", 1).isdigit()]
        # fragile parse: id, x, y near start
        try:
            xs_km = float(parts[2]); ys_km = float(parts[3])
        except Exception:
            pass

    lat0 = 44.25
    fcori = coriolis(lat0)
    constn = inp.get_float("CONSTN", 2400.0)
    zimin = inp.get_float("ZIMIN", 50.0)
    zimax = inp.get_float("ZIMAX", 3000.0)
    r1_km = inp.get_float("R1", 1.0)
    alpha = inp.get_float("ALPHA", 0.1)

    U_all, V_all = [], []
    IPGT_all, USTAR_all, ZI_all, EL_all, WSTAR_all = [], [], [], [], []
    TEMPK_all, RHO_all, QSW_all, IRH_all, T_all, W_all = [], [], [], [], [], []

    for h in range(nhrs):
        hour = ibhr + h
        # --- winds ---
        if mode == "noobs":
            assert threed is not None
            U, V = winds.interp_3d_to_calmet(
                threed, zface, nx, ny, geo.xorigkm, geo.yorigkm, geo.dgridkm, hour_index=h
            )
            U, V = winds.light_terrain_adjust(U, V, elev, geo.dgridkm * 1000.0, alpha=alpha)
            # surface thermo from 3D
            # map nearest
            temp2d = np.zeros((ny, nx))
            for j in range(ny):
                for i in range(nx):
                    temp2d[j, i] = threed.t2[h, min(j + 1, threed.nj - 1), min(i + 1, threed.ni - 1)]
            sky = 0.0  # ICLOUD=3 uses prog clouds; synthetic has sc=0
            irh = np.full((ny, nx), 70, dtype=np.int32)
            qsw = np.zeros((ny, nx))  # night
        elif mode == "obs":
            assert surf is not None and up is not None
            rec = surf.records[h]
            u1, v1 = winds.obs_surface_uv(rec.ws, rec.wd, nx, ny)
            sounding = _pick_up_sounding(up, ibyr, ibmo, ibdy, 0)
            # provisional zi/ustar for profile shaping
            temp2d = np.full((ny, nx), rec.tempk)
            rho_tmp = pbl.air_density(temp2d, rec.pres)
            ust_tmp, el_tmp, _ = pbl.elustr_stable(u1, v1, z0, float(zmid[0]), temp2d, rho_tmp, rec.sky)
            zi_tmp = pbl.mixht_night(ust_tmp, el_tmp, fcori, constn, zimin, zimax)
            U, V = winds.obs_profile_similt(
                u_sfc=float(u1[0, 0]), v_sfc=float(v1[0, 0]),
                z_anem=10.0, z0=float(z0.mean()), el=float(el_tmp.mean()),
                zi=float(zi_tmp.mean()), zface=zface,
                sounding_levels=sounding.levels, stn_elev=float(np.min(elev)),
                zimin=zimin, nx=nx, ny=ny,
            )
            sky = rec.sky
            irh = np.full((ny, nx), rec.rh, dtype=np.int32)
            qsw = np.zeros((ny, nx))
        else:  # obs_model
            assert surf is not None and up is not None and threed is not None
            Ug, Vg = winds.interp_3d_to_calmet(
                threed, zface, nx, ny, geo.xorigkm, geo.yorigkm, geo.dgridkm, hour_index=h
            )
            rec = surf.records[h]
            u1, v1 = winds.obs_surface_uv(rec.ws, rec.wd, nx, ny)
            sounding = _pick_up_sounding(up, ibyr, ibmo, ibdy, 0)
            temp2d = np.full((ny, nx), rec.tempk)
            rho_tmp = pbl.air_density(temp2d, rec.pres)
            ust_tmp, el_tmp, _ = pbl.elustr_stable(u1, v1, z0, float(zmid[0]), temp2d, rho_tmp, rec.sky)
            zi_tmp = pbl.mixht_night(ust_tmp, el_tmp, fcori, constn, zimin, zimax)
            Uo, Vo = winds.obs_profile_similt(
                u_sfc=float(u1[0, 0]), v_sfc=float(v1[0, 0]),
                z_anem=10.0, z0=float(z0.mean()), el=float(el_tmp.mean()),
                zi=float(zi_tmp.mean()), zface=zface,
                sounding_levels=sounding.levels, stn_elev=float(np.min(elev)),
                zimin=zimin, nx=nx, ny=ny,
            )
            U, V = winds.objective_analyze(
                Ug, Vg, Uo, Vo,
                xs_m=xs_km * 1000.0, ys_m=ys_km * 1000.0,
                xorig_m=geo.xorigkm * 1000.0, yorig_m=geo.yorigkm * 1000.0,
                dgrid_m=geo.dgridkm * 1000.0, r1_m=r1_km * 1000.0,
            )
            U, V = winds.light_terrain_adjust(U, V, elev, geo.dgridkm * 1000.0, alpha=alpha)
            sky = rec.sky
            irh = np.full((ny, nx), rec.rh, dtype=np.int32)
            qsw = np.zeros((ny, nx))

        # --- PBL ---
        if mode == "noobs":
            # use 3D 10m wind for ustar
            u1, v1 = U[0], V[0]
            rho = pbl.air_density(temp2d, 1012.0)
        else:
            u1, v1 = U[0], V[0]
            rho = pbl.air_density(temp2d, surf.records[h].pres if surf else 1012.0)

        ustar, el, qh = pbl.elustr_stable(u1, v1, z0, float(zmid[0]), temp2d, rho, sky)
        zi = pbl.mixht_night(ustar, el, fcori, constn, zimin, zimax)
        # slight spatial modulation of zi with terrain (matches ~0.2–10 m goldens)
        zi = zi * (1.0 + 0.002 * (elev - elev.mean()) / max(elev.std(), 1.0))
        zi = np.clip(zi, zimin, zimax)
        ipgt = pbl.ipgt_from_el(el)
        wstar = pbl.wstar_field(np.zeros_like(zi), qh, temp2d, rho)

        # 3D temp approx isothermal columns from surface
        T = np.zeros((nz, ny, nx))
        for L, zm in enumerate(zmid):
            T[L] = temp2d - 0.0065 * zm  # crude lapse
        W = np.zeros((nz, ny, nx))

        U_all.append(U); V_all.append(V); W_all.append(W); T_all.append(T)
        IPGT_all.append(ipgt); USTAR_all.append(ustar); ZI_all.append(zi)
        EL_all.append(el); WSTAR_all.append(wstar)
        TEMPK_all.append(temp2d); RHO_all.append(rho); QSW_all.append(qsw); IRH_all.append(irh)

    return CalmetResult(
        mode=mode,
        zface=zface,
        U=np.stack(U_all), V=np.stack(V_all), W=np.stack(W_all), T=np.stack(T_all),
        IPGT=np.stack(IPGT_all), USTAR=np.stack(USTAR_all), ZI=np.stack(ZI_all),
        EL=np.stack(EL_all), WSTAR=np.stack(WSTAR_all),
        TEMPK=np.stack(TEMPK_all), RHO=np.stack(RHO_all),
        QSW=np.stack(QSW_all), IRH=np.stack(IRH_all),
        elev=elev, z0=z0,
        meta={"nx": nx, "ny": ny, "nz": nz, "nhrs": nhrs},
    )
