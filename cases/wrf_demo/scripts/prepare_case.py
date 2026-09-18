#!/usr/bin/env python3
"""Build wrf_demo CALMET inputs (GEO/SURF/UP/3D/INP) from NCAR tutorial wrfout.

Synthesizes a few SURF/UP stations from near-surface WRF fields when real
obs are unavailable (documented in cases/wrf_demo/README.md).
"""
from __future__ import annotations

import math
import shutil
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]  # cases/wrf_demo
REPO = ROOT.parents[1]
SHARED = ROOT / "shared"
GOLD_IN = ROOT / "goldens" / "inputs"
RAW = ROOT / "raw"

sys.path.insert(0, str(REPO / "scripts"))
from wrfout_to_3d import convert_wrfout  # noqa: E402

# CALMET grid: interior of 14x14 MM5 subset
NX = NY = 12
DGRIDKM = 30.0  # match WRF DX
XORIGKM = 30.0  # one MM5 cell inset (MM5 x0=0)
YORIGKM = 30.0
NZ = 8
ZFACE = [0.0, 20.0, 40.0, 80.0, 160.0, 300.0, 600.0, 1000.0, 1500.0]
UTMZN = 16
UTMHEM = "N"
DATUM = "WGS-84"
PMAP = "UTM"
SFC_NAME, SFC_ID = "WRF1", 88801
UP_NAME, UP_ID = "WRF1", 88801
ANEM_HT = 10.0
ABTZ = "UTC+0000"
STN_TZ = 0

# WRF subset window (0-based mass indices into 90x73 domain)
I0, J0, NI, NJ = 4, 0, 14, 14


def _julian(dt: datetime) -> int:
    return int(dt.strftime("%j"))


def write_geo(path: Path, elev: np.ndarray, lu: np.ndarray) -> None:
    """Write GEO.DAT for CALMET 12x12 interior (j=1..12 ↔ MM5 j=2..13)."""
    # elev/lu are [14,14] MM5; take [1:13, 1:13]
    e = elev[1 : 1 + NY, 1 : 1 + NX]
    l = lu[1 : 1 + NY, 1 : 1 + NX]
    lines = [
        "GEO.DAT         2.0             Header structure with coordinate parameters",
        "   1",
        "wrf_demo GEO from NCAR wrfout_d01_2005-08-28 (Katrina tutorial)",
        PMAP,
        f"  {UTMZN}{UTMHEM}",
        f"{DATUM:<8s}08-28-2005",
        f"{NX:8d}{NY:8d}{XORIGKM:12.3f}{YORIGKM:12.3f}{DGRIDKM:12.3f}{DGRIDKM:12.3f}",
        "KM  M",
        "0                 -  LAND USE DATA  - IOPT1:  0=DEFAULT CATEGORIES  1=NEW CATEGORIES",
    ]
    for j in range(NY):
        lines.append(" ".join(f"{int(l[j, i]):3d}" for i in range(NX)))
    lines.append("1.0               -  TERRAIN HEIGHTS - HTFAC (units conversion to meters)")
    for j in range(NY):
        lines.append("".join(f"{e[j, i]:8.1f}" for i in range(NX)))
    lines += [
        "0 - IOPT2 (z0)   -- (0=default z0-lu table,    1=new z0-lu table,    2=gridded)",
        "0 - IOPT3 (alb) -- (0=default albedo-lu table,1=new albedo-lu table,2=gridded)",
        "0 - IOPT4 (Bo)  -- (0=default Bowen-lu table, 1=new Bowen-lu table, 2=gridded)",
        "0 - IOPT5 (HCG) -- (0=default HCG-lu table,   1=new HCG-lu table,   2=gridded)",
        "0 - IOPT6 (QF)  -- (0=default QF-lu table,    1=new QF-lu table,    2=gridded)",
        "0 - IOPT7 (LAI) -- (0=default XLAI-lu table,  1=new XLAI-lu table,  2=gridded)",
    ]
    path.write_text("\n".join(lines) + "\n")


def write_surf(path: Path, meta: dict) -> None:
    """Write SURF for 3 timesteps at WRF hours 00/03/06 (NSECDT=10800)."""
    times: list[datetime] = meta["times"]
    run_times = times[:3]
    t0, t_end = run_times[0], times[3]  # end exclusive = 09Z
    year = t0.year
    jul0 = _julian(t0)
    cy, cj = NJ // 2, NI // 2
    lines = [
        "SURF.DAT        2.1             Hour Start and End Times with Seconds",
        "   1",
        "Synthetic SURF from wrfout near-surface (center cell); not real observations",
        "NONE",
        ABTZ,
        f"{year:5d}{jul0:5d}{t0.hour:4d}{0:6d}{t_end.year:6d}{_julian(t_end):5d}{t_end.hour:4d}{0:6d}{1:5d}",
        f"{SFC_ID:8d}",
    ]
    for it, t in enumerate(run_times):
        t_next = times[it + 1]
        ws = float(meta["ws10"][it, cy, cj])
        wd = float(meta["wd10"][it, cy, cj])
        temp = float(meta["t2"][it, cy, cj])
        rh = int(np.clip(100.0 * float(meta["q2_gkg"][it, cy, cj]) / 12.0, 30, 95))
        # Prefer RH from Q2 roughly; clamp
        pres = float(meta["psfc_mb"][it, cy, cj])
        ceil, sky, pp = 50, 4, 0
        lines.append(
            f"{t.year:4d}{_julian(t):5d}{t.hour:4d}{0:5d}"
            f"{t_next.year:6d}{_julian(t_next):5d}{t_next.hour:4d}{0:5d}"
        )
        lines.append(
            f"{ws:8.3f}{wd:9.3f}{ceil:5d}{sky:5d}{temp:9.3f}{rh:5d}{pres:9.3f}{pp:5d}"
        )
    path.write_text("\n".join(lines) + "\n")


def write_up(path: Path, meta: dict) -> None:
    times: list[datetime] = meta["times"]
    t0 = times[0]
    year = t0.year
    # Cover day before through day after
    d0 = t0 - timedelta(days=1)
    d1 = t0 + timedelta(days=1)
    lines = [
        "UP.DAT          2.1             Hour Start and End Times with Seconds",
        "   1",
        "Synthetic UP from wrfout column at domain center; not real radiosonde",
        "NONE",
        ABTZ,
        f"{d0.year:5d}{_julian(d0):6d}{0:5d}{0:5d}{d1.year:5d}{_julian(d1):5d}{0:5d}{0:5d}{500.:6.0f}{2:5d}{2:5d}",
        "     F    F    F    F",
    ]
    cy, cj = NJ // 2, NI // 2
    nk = meta["pres_mb"].shape[1]
    # Soundings at 00Z and 12Z around the case day
    sound_hours = [
        (t0.year, t0.month, t0.day - 1 if t0.day > 1 else 27, 0),
        (t0.year, t0.month, t0.day - 1 if t0.day > 1 else 27, 12),
        (t0.year, t0.month, t0.day, 0),
        (t0.year, t0.month, t0.day, 12),
        (t0.year, t0.month, t0.day + 1, 0),
    ]
    # Use WRF time index 0 for 00Z profile; time 1 (~03Z) for slight variation at 12Z
    for yi, mo, dy, hr in sound_hours:
        ti = 0 if hr == 0 else min(1, len(times) - 1)
        nlev = nk
        hdr = (
            f"{'':9s}{UP_ID:8d}"
            f"{'':4s}{yi:4d}{mo:4d}{dy:3d}{hr:3d}{0:5d}"
            f"{'':4s}{yi:4d}{mo:4d}{dy:3d}{hr:3d}{0:5d}"
            f"{'':8s}{nlev:5d}"
        )
        lines.append(hdr)
        parts = []
        for k in range(nk):
            p = float(meta["pres_mb"][ti, k, cy, cj])
            z = float(meta["z_msl"][ti, k, cy, cj])
            tc = float(meta["tempk"][ti, k, cy, cj]) - 273.15
            wd = float(meta["wd"][ti, k, cy, cj])
            ws = float(meta["ws"][ti, k, cy, cj])
            parts.append(f"{p:.1f},{z:.0f},{tc:.1f},{wd:.0f},{ws:.1f}")
        for i in range(0, nlev, 4):
            lines.append(",".join(parts[i : i + 4]))
    path.write_text("\n".join(lines) + "\n")


def make_inp(mode: str, path: Path, year: int, month: int, day: int, ibhr: int, iehr: int = 9) -> None:
    if mode == "obs":
        noobs, iprog, nm3d, nssta, nusta = 0, 0, 0, 1, 1
        itprog, irhprog, icloud, iextrp = 0, 0, 0, -4
        m3d_block = ""
        up_block = "UP1.DAT       input     1  ! UPDAT=up.dat!    !END!\n"
        srf_line = "! SRFDAT=surf.dat      !"
    elif mode == "obs_model":
        noobs, iprog, nm3d, nssta, nusta = 0, 14, 1, 1, 1
        itprog, irhprog, icloud, iextrp = 0, 0, 0, -4
        m3d_block = "MM51.DAT       input     1  ! M3DDAT=3d.dat!    !END!\n"
        up_block = "UP1.DAT       input     1  ! UPDAT=up.dat!    !END!\n"
        srf_line = "! SRFDAT=surf.dat      !"
    else:
        noobs, iprog, nm3d, nssta, nusta = 2, 14, 1, 0, 0
        itprog, irhprog, icloud, iextrp = 2, 1, 3, 1
        m3d_block = "MM51.DAT       input     1  ! M3DDAT=3d.dat!    !END!\n"
        up_block = ""
        srf_line = "* SRFDAT=             *"

    sfc_x = XORIGKM + (NX / 2.0) * DGRIDKM
    sfc_y = YORIGKM + (NY / 2.0) * DGRIDKM
    bias = ", ".join(["0"] * NZ)
    nsmth = ", ".join(["2"] + ["4"] * (NZ - 1))
    fextr2 = ", ".join(["0"] * NZ)
    zface = ", ".join(str(z) for z in ZFACE)

    grp7 = ""
    grp8 = ""
    if nssta > 0:
        grp7 = f"""
INPUT GROUP: 7 -- Surface meteorological station parameters
--------------
! SS1  ='{SFC_NAME:<4s}'  {SFC_ID:6d}    {sfc_x:10.3f}  {sfc_y:10.3f}  {STN_TZ:3d}  {ANEM_HT:5.0f}  !
!END!
"""
    if nusta > 0:
        grp8 = f"""
INPUT GROUP: 8 -- Upper air meteorological station parameters
--------------
! US1  ='{UP_NAME:<4s}'  {UP_ID:6d}  {sfc_x:10.3f}  {sfc_y:10.3f}  {STN_TZ:3d}  !
!END!
"""

    text = f"""CALMET.INP      2.1             Hour Start and End Times with Seconds
wrf_demo mode={mode}
12x12 @ 30km, UTM {UTMZN}{UTMHEM}, from NCAR Katrina wrfout
---------------- Run title (3 lines) ------------------------------------------

                    CALMET MODEL CONTROL FILE
                    --------------------------

-------------------------------------------------------------------------------

INPUT GROUP: 0 -- Input and Output File Names


Subgroup (a)
------------
Default Name  Type          File Name
------------  ----          ---------
GEO.DAT       input    ! GEODAT=geo.dat       !
SURF.DAT      input    {srf_line}
CLOUD.DAT     input    * CLDDAT=            *
PRECIP.DAT    input    * PRCDAT=            *
WT.DAT        input    * WTDAT=             *

CALMET.LST    output   ! METLST=CALMET.LST     !
CALMET.DAT    output   ! METDAT=CALMET.DAT    !
PACOUT.DAT    output   * PACDAT=            *

         T = lower case      ! LCFILES = T !

    Number of upper air stations (NUSTA)  No default     ! NUSTA =  {nusta}  !
    Number of overwater met stations
                                 (NOWSTA) No default     ! NOWSTA =  0  !

    Number of MM4/MM5/3D.DAT files
                                 (NM3D) No default       ! NM3D =  {nm3d}  !

    Number of IGF-CALMET.DAT files
                                 (NIGF)   No default     ! NIGF =  0  !

                       !END!
--------------------------------------------------------------------------------
Subgroup (b)
---------------------------------
Upper air files (one per station)
---------------------------------
{up_block}--------------------------------------------------------------------------------
Subgroup (c)
-----------------------------------------
Overwater station files (one per station)
-----------------------------------------
--------------------------------------------------------------------------------
Subgroup (d)
------------------------------------------------
MM4/MM5/3D.DAT files (consecutive or overlapping)
------------------------------------------------
{m3d_block}--------------------------------------------------------------------------------
Subgroup (e)
-------------------------------------------------
IGF-CALMET.DAT files (consecutive or overlapping)
-------------------------------------------------
--------------------------------------------------------------------------------
Subgroup (f)
----------------
Other file names
----------------

Default Name  Type       File Name
------------  ----       ---------
DIAG.DAT      input      * DIADAT=                  *
PROG.DAT      input      * PRGDAT=                  *

TEST.PRT      output     * TSTPRT=                  *
TEST.OUT      output     * TSTOUT=                  *
TEST.KIN      output     * TSTKIN=                  *
TEST.FRD      output     * TSTFRD=                  *
TEST.SLP      output     * TSTSLP=                  *
DCST.GRD      output     * DCSTGD=                  *

                         !END!


-------------------------------------------------------------------------------

INPUT GROUP: 1 -- General run control parameters
--------------

     Starting date:    Year   (IBYR)  --    No default   ! IBYR  =  {year}  !
                       Month  (IBMO)  --    No default   ! IBMO  =  {month}  !
                       Day    (IBDY)  --    No default   ! IBDY  =  {day}  !
     Starting time:    Hour   (IBHR)  --    No default   ! IBHR  =  {ibhr}  !
                       Second (IBSEC) --    No default   ! IBSEC =  0  !

     Ending date:      Year   (IEYR)  --    No default   ! IEYR  =  {year}  !
                       Month  (IEMO)  --    No default   ! IEMO  =  {month}  !
                       Day    (IEDY)  --    No default   ! IEDY  =  {day}  !
     Ending time:      Hour   (IEHR)  --    No default   ! IEHR  =  {iehr}  !
                       Second (IESEC) --    No default   ! IESEC =  0  !

      UTC time zone         (ABTZ) -- No default       ! ABTZ= UTC+0000 !

     Length of modeling time-step (seconds)
     (NSECDT)                        Default:3600     ! NSECDT =  10800  !

     Run type            (IRTYPE) -- Default: 1       ! IRTYPE=  1  !

     Compute special data fields required
     by CALGRID (LCALGRD)            Default: T    ! LCALGRD = T !

      Flag to stop run after
      SETUP phase (ITEST)             Default: 2       ! ITEST=  2   !

     Test options specified to see if
     they conform to regulatory
     values? (MREG)                   No Default       ! MREG =   0   !

!END!

-------------------------------------------------------------------------------

INPUT GROUP: 2 -- Map Projection and Grid control parameters
--------------

     Map projection
     (PMAP)                     Default: UTM    ! PMAP = UTM  !

     False Easting and Northing (km) at the projection origin
     (FEAST)                    Default=0.0   ! FEAST  = 0.0 !
     (FNORTH)                   Default=0.0   ! FNORTH = 0.0 !

     UTM zone (1 to 60)
     (IUTMZN)                   No Default      ! IUTMZN =  {UTMZN}   !

     Hemisphere for UTM projection?
     (UTMHEM)                   Default: N    ! UTMHEM = {UTMHEM}  !

     Datum-region
     (DATUM)                    Default: WGS-84    ! DATUM = WGS-84  !

     NX  NY
     (NX)                       No Default      ! NX =  {NX}  !
     (NY)                       No Default      ! NY =  {NY}  !

     Grid spacing (km)
     (DGRIDKM)                  No Default      ! DGRIDKM =  {DGRIDKM}  !

     Reference coordinates of the SW corner of grid cell (1,1) (km)
     (XORIGKM)                  No Default      ! XORIGKM = {XORIGKM}  !
     (YORIGKM)                  No Default      ! YORIGKM = {YORIGKM}  !

     Number of vertical layers
     (NZ)                       No Default      ! NZ =  {NZ}  !

     Cell face heights (m) in the
     vertical grid (ZFACE(NZ+1))    No defaults
     ! ZFACE = {zface}  !

!END!

-------------------------------------------------------------------------------

INPUT GROUP: 3 -- Output options
--------------

        ! LSAVE = T !
        ! IOUTU = 1 !
        ! LDB = F !
        ! NN1 = 1 !
        ! NN2 = 1 !
        ! LDBCST = F !
        ! IOUTD = 0 !
        ! IOUTE = 0 !
        ! IOUTQ = 0 !
        ! IOUTC = 0 !

!END!

-------------------------------------------------------------------------------

INPUT GROUP: 4 -- Meteorological data options
--------------

    NO OBSERVATION MODE             (NOOBS)  Default: 0     ! NOOBS =  {noobs}   !

    Number of surface stations
    (NSSTA)  No default     ! NSSTA =  {nssta}  !

    Number of precipitation stations
    (NPSTA)  Default: 0     ! NPSTA =  0  !

       ! ICLOUD  =  {icloud}  !
       ! IRTYPE  =  1  !

       ! LFARAD  = F !

!END!

-------------------------------------------------------------------------------

INPUT GROUP: 5 -- Wind field options and parameters
--------------

       ! IWFCOD = 1 !
       ! IFRADJ = 1 !
       ! IKINE  = 0 !
       ! IOBRIEN = 0 !
       ! KBAR   = 999 !
       ! RMIN   = 0.1 !
       ! RMAX1  = 1.0 !
       ! RMAX2  = 1.0 !
       ! RMAX3  = 1.0 !
       ! R1     = 1.0 !
       ! R2     = 1.0 !
       ! RPROG  = 0.0 !

       ! IEXTRP = {iextrp} !
       ! IPROG  = {iprog}  !
       ! LVARY  = F !

       ! BIAS = {bias} !
       ! RMIN1 = 0.0 !
       ! RMIN2 = 0.0 !

       ! IPROG is used with 3D.DAT when IPROG=14 !

       ! ALPHA  = 0.1 !
       ! NITER  = 50 !
       ! NSMTH = {nsmth} !
       ! FEXTR2 = {fextr2} !

       ! IDIOPT1 = 0 !
       ! ISURFT  = 1 !
       ! IDIOPT2 = 0 !
       ! IUPT    = 1 !
       ! ZUPT    = 200. !
       ! IDIOPT3 = 0 !
       ! IUPWND  = -1 !
       ! ZUPWND  = 1., 1000. !
       ! IDIOPT4 = 0 !
       ! IDIOPT5 = 0 !

       ! LLBREZE = F !
       ! NBOX = 0 !

!END!

-------------------------------------------------------------------------------

INPUT GROUP: 6 -- Mixing Height, Temperature and Precipitation Parameters
--------------

       ! CONSTB  = 1.41 !
       ! CONSTE  = 0.15 !
       ! CONSTN  = 2400. !
       ! CONSTW  = 0.16 !
       ! FCORIOL = 1.0E-4 !

       ! IAVEZI = 1 !
       ! MNMDAV = 1 !
       ! HAFANG = 30. !
       ! ILEVZI = 1 !

       ! IMIXH   = 1 !
       ! THRESHL = 0.05 !
       ! THRESHW = 0.05 !
       ! ITWPROG = 0 !
       ! ILUOC3D = 16 !

       ! DPTMIN = 0.001 !
       ! DZZI   = 200. !
       ! ZIMIN  = 50. !
       ! ZIMAX  = 3000. !
       ! ZIMINW = 50. !
       ! ZIMAXW = 3000. !

       ! ICOARE = 10 !
       ! DSHELF = 0.0 !
       ! IWARM  = 0 !
       ! ICOOL  = 0 !

       ! IRHPROG = {irhprog} !
       ! ITPROG  = {itprog} !
       ! IRAD   = 1 !
       ! TRADKM = 500. !
       ! NUMTS  = 5 !
       ! IAVET  = 1 !
       ! TGDEFB = -0.0098 !
       ! TGDEFA = -0.0045 !
       ! JWAT1  = 999 !
       ! JWAT2  = 999 !

       ! NFLAGP = 2 !
       ! SIGMAP = 100. !
       ! CUTP   = 0.01 !

!END!
{grp7}{grp8}
"""
    path.write_text(text)


def main() -> None:
    SHARED.mkdir(parents=True, exist_ok=True)
    GOLD_IN.mkdir(parents=True, exist_ok=True)
    wrfout = RAW / "wrfout_d01_2005-08-28_00_00_00"
    if not wrfout.exists():
        raise SystemExit(
            f"Missing {wrfout}. Download NCAR wrf_tutorial_data "
            "(see cases/wrf_demo/README.md)."
        )

    meta = convert_wrfout(
        wrfout,
        SHARED / "3d.dat",
        i0=I0,
        j0=J0,
        ni=NI,
        nj=NJ,
        nk_out=10,
        time_indices=[0, 1, 2, 3],  # 00,03,06,09 UTC
    )
    global DGRIDKM, XORIGKM, YORIGKM, UTMZN
    DGRIDKM = round(float(meta["dx_km"]), 3)
    XORIGKM = round(float(meta["x0_km"] + meta["dx_km"]), 3)  # one MM5 cell inset
    YORIGKM = round(float(meta["y0_km"] + meta["dx_km"]), 3)
    UTMZN = int(meta["utm_zone"])
    # Rewrite 3D.DAT map header is already UTM-local via converter x0/y0
    write_geo(SHARED / "geo.dat", meta["elev"], meta["lu"])
    write_surf(SHARED / "surf.dat", meta)
    write_up(SHARED / "up.dat", meta)

    t0 = meta["times"][0]
    calmet_x = REPO / "vendor" / "calmet-fortran" / "src" / "calmet.x"
    for mode in ("obs", "obs_model", "noobs"):
        d = ROOT / mode
        d.mkdir(parents=True, exist_ok=True)
        make_inp(mode, d / "calmet.inp", t0.year, t0.month, t0.day, t0.hour, iehr=9)
        for f in ("geo.dat", "surf.dat", "up.dat", "3d.dat"):
            target = d / f
            if target.exists() or target.is_symlink():
                target.unlink()
            target.symlink_to(f"../shared/{f}")
        cx = d / "calmet.x"
        if cx.exists() or cx.is_symlink():
            cx.unlink()
        if calmet_x.exists():
            cx.symlink_to(calmet_x)

    # Archive inputs into goldens/inputs
    for f in ("geo.dat", "surf.dat", "up.dat", "3d.dat"):
        shutil.copy2(SHARED / f, GOLD_IN / f)
        # also copy inp stubs into goldens later after Fortran runs

    print("Prepared wrf_demo shared inputs:")
    print(f"  times={meta['times']}")
    print(f"  elev range {meta['elev'].min():.1f}-{meta['elev'].max():.1f} m")
    print(f"  3d.dat size {(SHARED / '3d.dat').stat().st_size / 1e6:.2f} MB")


if __name__ == "__main__":
    main()
