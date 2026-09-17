#!/usr/bin/env python3
"""Generate tiny synthetic CALMET demo domain (12x12, 3 hours) and three mode configs."""
from __future__ import annotations
import math
import os
from pathlib import Path
from pyproj import Transformer

ROOT = Path("/workspace/py-calmet/cases/small_domain")
SHARED = ROOT / "shared"

# Domain: UTM 19N, inland Maine-ish
UTMZN = 19
UTMHEM = "N"
DATUM = "WGS-84"
PMAP = "UTM"
NX = NY = 12
DGRIDKM = 1.0
XORIGKM = 400.0
YORIGKM = 4900.0
NZ = 8
ZFACE = [0.0, 20.0, 40.0, 80.0, 160.0, 300.0, 600.0, 1000.0, 1500.0]

# Time window: 2020-06-15 00:00–03:00 UTC (3 hours)
YEAR, MONTH, DAY = 2020, 6, 15
JULIAN = 167  # June 15, 2020
IBHR, IEHR = 0, 3
NHRS = IEHR - IBHR  # 3
M3D_NHRS = NHRS + 1  # one extra hour for prog interpolation bracket
ABTZ = "UTC+0000"
STN_TZ = 0  # relative to UTC for station params

# Stations near domain center
SFC_NAME, SFC_ID = "DEMO", 99901
UP_NAME, UP_ID = "DEMO", 99901
SFC_X = XORIGKM + (NX / 2.0) * DGRIDKM
SFC_Y = YORIGKM + (NY / 2.0) * DGRIDKM
ANEM_HT = 10.0
ELEV_M = 150

# 3D.DAT grid: slightly larger so CALMET is interior
M3D_I0, M3D_J0 = 1, 1
M3D_NI = M3D_NJ = 14  # covers 14 km
M3D_DX_KM = 1.0
M3D_X0 = XORIGKM - M3D_DX_KM  # one cell west/south of CALMET
M3D_Y0 = YORIGKM - M3D_DX_KM
# sigma (approx pressure/1013) descending from near-surface? Actually MM5 half-sigma
# typically decreases upward: ~0.995 near surface to smaller aloft
SIGMA = [0.995, 0.980, 0.950, 0.900, 0.850, 0.800, 0.700, 0.600, 0.500, 0.400]
PRES_MB = [int(round(s * 1013)) for s in SIGMA]
NZP = len(SIGMA)

to_ll = Transformer.from_crs("EPSG:32619", "EPSG:4326", always_xy=True)


def utm_to_ll(x_km: float, y_km: float) -> tuple[float, float]:
    lon, lat = to_ll.transform(x_km * 1000.0, y_km * 1000.0)
    return lat, lon


def write_geo(path: Path) -> None:
    lines = []
    lines.append("GEO.DAT         2.0             Header structure with coordinate parameters")
    lines.append("   1")
    lines.append("Tiny synthetic demo domain for py-calmet goldens")
    lines.append("UTM")
    lines.append(f"  {UTMZN}{UTMHEM}")
    lines.append(f"{DATUM:<8s}06-15-2020")
    lines.append(
        f"{NX:8d}{NY:8d}{XORIGKM:12.3f}{YORIGKM:12.3f}{DGRIDKM:12.3f}{DGRIDKM:12.3f}"
    )
    lines.append("KM  M")
    # IOPT1=0 default LU categories
    lines.append("0                 -  LAND USE DATA  - IOPT1:  0=DEFAULT CATEGORIES  1=NEW CATEGORIES")
    # LU rows: for each Y (north to south? MAKEGEO writes j outer as i rows of nx)
    # From make_GEO: for i in 1..nx: for j in 1..ny: print LU — so nx rows of ny? Wait:
    # for i in seq 1 nx; for j in seq 1 ny; printf — each i is a row with ny values.
    # Sample GEO has 99 values per line (nx=99), and 99 such lines (ny).
    # FILLGEO typically reads NY rows of NX values (row=j, col=i).
    # Sample: each line has 99 values = NX, and there are 99 lines = NY.
    for j in range(NY):
        row = "".join(f"{20:6d}" for _ in range(NX))  # agricultural
        lines.append(row)
    lines.append("1.0               -  TERRAIN HEIGHTS - HTFAC")
    # gentle slope + bump
    elevs = []
    for j in range(NY):
        for i in range(NX):
            e = ELEV_M + 5.0 * i + 3.0 * j + 20.0 * math.exp(
                -((i - 5) ** 2 + (j - 6) ** 2) / 8.0
            )
            elevs.append(e)
    # write NY rows of NX
    idx = 0
    for j in range(NY):
        row = "".join(f"{elevs[idx + i]:8.1f}" for i in range(NX))
        lines.append(row)
        idx += NX
    lines.append("0 - IOPT2 (z0)   -- (0=default z0-lu table,    1=new z0-lu table,    2=gridded)")
    lines.append("0 - IOPT3 (alb) -- (0=default albedo-lu table,1=new albedo-lu table,2=gridded)")
    lines.append("0 - IOPT4 (Bo)  -- (0=default Bowen-lu table, 1=new Bowen-lu table, 2=gridded)")
    lines.append("0 - IOPT5 (HCG) -- (0=default HCG-lu table,   1=new HCG-lu table,   2=gridded)")
    lines.append("0 - IOPT6 (QF)  -- (0=default QF-lu table,    1=new QF-lu table,    2=gridded)")
    lines.append("0 - IOPT7 (LAI) -- (0=default XLAI-lu table,  1=new XLAI-lu table,  2=gridded)")
    path.write_text("\n".join(lines) + "\n")


def write_surf(path: Path) -> None:
    lines = []
    lines.append("SURF.DAT        2.1             Hour Start and End Times with Seconds")
    lines.append("   1")
    lines.append("Synthetic SURF.DAT for py-calmet tiny domain")
    lines.append("NONE")
    lines.append(ABTZ)
    # begin y j h s  end y j h s  nsta
    lines.append(
        f"{YEAR:5d}{JULIAN:5d}{IBHR:4d}{0:6d}{YEAR:6d}{JULIAN:5d}{IEHR:4d}{0:6d}{1:5d}"
    )
    lines.append(f"{SFC_ID:8d}")
    # hourly records: begin/end then one station line
    # ws wd ceil sky tempK rh% pres_mb ppcode
    for h in range(IBHR, IEHR):
        lines.append(
            f"{YEAR:4d}{JULIAN:5d}{h:4d}{0:5d}{YEAR:6d}{JULIAN:5d}{h+1:4d}{0:5d}"
        )
        ws = 3.5 + 0.2 * (h - IBHR)
        wd = 220.0 + 5.0 * (h - IBHR)
        ceil = 50
        sky = 3
        temp = 291.0 + 0.5 * (h - IBHR)
        rh = 65
        pres = 1012.0
        pp = 0
        lines.append(
            f"{ws:8.3f}{wd:9.3f}{ceil:5d}{sky:5d}{temp:9.3f}{rh:5d}{pres:9.3f}{pp:5d}"
        )
    path.write_text("\n".join(lines) + "\n")


def write_up(path: Path) -> None:
    """UP.DAT 2.1, comma-delimited levels, UTC+0000, station ID matches US1."""
    lines = []
    lines.append("UP.DAT          2.1             Hour Start and End Times with Seconds")
    lines.append("   1")
    lines.append("Synthetic UP.DAT for py-calmet tiny domain")
    lines.append("NONE")
    lines.append(ABTZ)
    # begin y j h s end y j h s ptop jdat ifmt
    # Cover from 00Z day before through 00Z next day so run is bracketed
    # Julian 166 = June 14, 167 = June 15, 168 = June 16
    lines.append(
        f"{YEAR:5d}{166:6d}{0:5d}{0:5d}{YEAR:5d}{168:5d}{0:5d}{0:5d}{500.:6.0f}{2:5d}{2:5d}"
    )
    lines.append("     F    F    F    F")

    # Soundings at 00Z and 12Z on Jun 14, 15, 16 (snapshots: begin==end)
    sound_times = [
        (YEAR, 6, 14, 0),
        (YEAR, 6, 14, 12),
        (YEAR, 6, 15, 0),
        (YEAR, 6, 15, 12),
        (YEAR, 6, 16, 0),
    ]
    # pressure(mb), height(m MSL), temp(C), wd, ws
    base_levels = [
        (1000.0, 150.0, 18.0, 220.0, 4.0),
        (975.0, 370.0, 16.5, 225.0, 5.0),
        (950.0, 600.0, 15.0, 230.0, 6.0),
        (925.0, 840.0, 13.5, 235.0, 7.0),
        (900.0, 1090.0, 12.0, 240.0, 8.0),
        (850.0, 1550.0, 9.0, 245.0, 10.0),
        (800.0, 2050.0, 6.0, 250.0, 12.0),
        (700.0, 3100.0, 0.0, 255.0, 15.0),
        (600.0, 4300.0, -7.0, 260.0, 18.0),
        (500.0, 5700.0, -15.0, 265.0, 22.0),
    ]
    nlev = len(base_levels)
    for yi, mo, dy, hr in sound_times:
        # format(9x,i8,2(4x,i4,i4,i3,i3,i5),8x,i5)
        # 9 spaces + i8 ID + (4x i4 i4 i3 i3 i5)*2 + 8x + i5 nlev
        # Put junk in first 9 cols then ID
        hdr = (
            f"{'':9s}{UP_ID:8d}"
            f"{'':4s}{yi:4d}{mo:4d}{dy:3d}{hr:3d}{0:5d}"
            f"{'':4s}{yi:4d}{mo:4d}{dy:3d}{hr:3d}{0:5d}"
            f"{'':8s}{nlev:5d}"
        )
        lines.append(hdr)
        # comma-delimited all levels on one or more lines
        parts = []
        for k, (p, z, t, wd, ws) in enumerate(base_levels):
            # slight time variation
            t2 = t + 0.3 * math.sin(hr / 12.0 * math.pi)
            wd2 = wd + hr * 0.5
            ws2 = ws + 0.1 * hr
            parts.append(f"{p:.1f},{z:.0f},{t2:.1f},{wd2:.0f},{ws2:.1f}")
        # wrap ~4 levels per line
        for i in range(0, nlev, 4):
            lines.append(",".join(parts[i : i + 4]))
    path.write_text("\n".join(lines) + "\n")


def write_3d(path: Path) -> None:
    """3D.DAT dataset 2.1 following hrrr2calmet.py conventions."""
    lines = []
    lines.append(
        f"{'3D.DAT':<16s}{'2.1':<16s}{'Header Structure with Comment Lines':<64s}"
    )
    lines.append("   1")
    lines.append(f"{'Synthetic 3D.DAT for py-calmet tiny domain':<132s}")
    # ioutw ioutq ioutc iouti ioutg
    lines.append("  1  1  0  0  0")  # ioutw,ioutq,ioutc,iouti,ioutg -> ioutmm5=92
    # map projection (skipped by reader but write UTM-ish LCC-compatible line)
    # hrrr uses: LCC lat0 lon0 lat1 lat2 x0 y0 dx ni0 nj0 nk0
    # Use LCC that approx matches region for header completeness
    lat0, lon0 = 44.3, -70.2
    lines.append(
        f"{'LCC':<4s}{lat0:9.4f}{lon0:10.4f}{30.0:7.2f}{60.0:7.2f}"
        f"{M3D_X0:10.3f}{M3D_Y0:10.3f}{M3D_DX_KM:8.3f}"
        f"{M3D_NI:4d}{M3D_NJ:3d}{NZP:3d}"
    )
    # surface variable options (21 ints)
    lines.append("  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0")
    # grid data: YYYYMMDDHH + nhrs + ni + nj + nk  format i4,3i2,i5,3i4
    ymdh = f"{YEAR:04d}{MONTH:02d}{DAY:02d}{IBHR:02d}"
    lines.append(f"{ymdh}{M3D_NHRS:5d}{M3D_NI:4d}{M3D_NJ:4d}{NZP:4d}")
    # extraction subdomain
    # compute lat/lon range
    lats, lons = [], []
    for jj in range(M3D_NJ):
        for ii in range(M3D_NI):
            x = M3D_X0 + (ii + 0.5) * M3D_DX_KM
            y = M3D_Y0 + (jj + 0.5) * M3D_DX_KM
            lat, lon = utm_to_ll(x, y)
            lats.append(lat)
            lons.append(lon)
    lines.append(
        f"{M3D_I0:4d}{M3D_J0:4d}{M3D_I0+M3D_NI-1:4d}{M3D_J0+M3D_NJ-1:4d}"
        f"{1:4d}{NZP:4d}{min(lons):10.4f}{max(lons):10.4f}"
        f"{min(lats):9.4f}{max(lats):9.4f}"
    )
    for s in SIGMA:
        lines.append(f"{s:6.3f}")
    # surface gridpoint metadata: i j lat lon elev lu + dummies
    for jj in range(M3D_NJ):
        for ii in range(M3D_NI):
            i_idx = M3D_I0 + ii
            j_idx = M3D_J0 + jj
            x = M3D_X0 + (ii + 0.5) * M3D_DX_KM
            y = M3D_Y0 + (jj + 0.5) * M3D_DX_KM
            lat, lon = utm_to_ll(x, y)
            elev = int(
                ELEV_M
                + 5.0 * ii
                + 3.0 * jj
                + 20.0 * math.exp(-((ii - 6) ** 2 + (jj - 7) ** 2) / 8.0)
            )
            lu = 20
            lines.append(
                f"{i_idx:4d}{j_idx:4d}{lat:9.4f}{lon:10.4f}{elev:5d}{lu:3d}"
                f" {-999:9.4f}{-999:10.4f}{-999:5d}"
            )

    # data hours: for each hour, for each gridpoint: surface then upper levels
    # surface template from hrrr2calmet:
    # '%10s%03d%03d%7.1f%5.2f%2d%8.1f%8.1f%8.1f%8.2f%8.1f%8.1f%8.1f'
    # upper: '%4d%6d%6.1f%4d%5.1f%6.2f%3d%5.2f%6.3f'
    for h_off in range(M3D_NHRS):
        hr = IBHR + h_off
        stamp = f"{YEAR:04d}{MONTH:02d}{DAY:02d}{hr:02d}"
        for jj in range(M3D_NJ):
            for ii in range(M3D_NI):
                i_idx = M3D_I0 + ii
                j_idx = M3D_J0 + jj
                elev = ELEV_M + 5.0 * ii + 3.0 * jj
                spres = 1013.0 - elev / 8.5
                rain = 0.00
                sc = 0
                radsw = 600.0 + 50.0 * h_off
                radlw = 300.0
                t2 = 291.0 + 0.5 * h_off + 0.01 * (ii + jj)
                q2 = 0.008  # g/kg? hrrr uses g/kg after /1000 of kg/kg — actually sh2/1000 gives tiny; sample uses ~8 g/kg style as f8.2
                # Looking at hrrr: q2 = ds[1].sh2 / 1000 — if sh2 is kg/kg (~0.008), then q2=8e-6 which is wrong.
                # Probably sh2 already in g/kg or bug. Use ~8 g/kg.
                q2 = 8.0
                wd10 = 220.0 + 5.0 * h_off
                ws10 = 3.5 + 0.2 * h_off
                sst = t2 + 1.0
                lines.append(
                    f"{stamp}{i_idx:03d}{j_idx:03d}{spres:7.1f}{rain:5.2f}{sc:2d}"
                    f"{radsw:8.1f}{radlw:8.1f}{t2:8.1f}{q2:8.2f}"
                    f"{wd10:8.1f}{ws10:8.1f}{sst:8.1f}"
                )
                for k, (sig, pmb) in enumerate(zip(SIGMA, PRES_MB)):
                    # height approx hydrostatics from elev
                    z = int(elev + (1.0 - sig) * 9000.0)
                    tempk = 291.0 - (1.0 - sig) * 60.0 + 0.2 * h_off
                    wd = int(220 + 10 * (1.0 - sig) + h_off)
                    ws = 4.0 + 15.0 * (1.0 - sig)
                    w = 0.01
                    rh = int(70 - 20 * (1.0 - sig))
                    vapmr = 6.0 * sig
                    cldmr = 0.001
                    # format 92: i4,i6,f6.1,i4,f5.1,f6.2,i3,f5.2
                    lines.append(
                        f"{pmb:4d}{z:6d}{tempk:6.1f}{wd:4d}{ws:5.1f}"
                        f"{w:6.2f}{rh:3d}{vapmr:5.2f}"
                    )
    path.write_text("\n".join(lines) + "\n")


def make_inp(mode: str, path: Path) -> None:
    """Write CALMET.INP for obs / obs_model / noobs."""
    if mode == "obs":
        noobs, iprog, nm3d, nssta, nusta, npsta = 0, 0, 0, 1, 1, 0
        itprog, irhprog, icloud = 0, 0, 0
        iextrp = -4
        nowsta = 0
        m3d_block = ""
        up_block = "UP1.DAT       input     1  ! UPDAT=up.dat!    !END!\n"
        srf_line = "! SRFDAT=surf.dat      !"
        geodat = "geo.dat"
    elif mode == "obs_model":
        noobs, iprog, nm3d, nssta, nusta, npsta = 0, 14, 1, 1, 1, 0
        itprog, irhprog, icloud = 0, 0, 0
        iextrp = -4
        nowsta = 0
        m3d_block = "MM51.DAT       input     1  ! M3DDAT=3d.dat!    !END!\n"
        up_block = "UP1.DAT       input     1  ! UPDAT=up.dat!    !END!\n"
        srf_line = "! SRFDAT=surf.dat      !"
        geodat = "geo.dat"
    elif mode == "noobs":
        noobs, iprog, nm3d, nssta, nusta, npsta = 2, 14, 1, 0, 0, 0
        itprog, irhprog, icloud = 2, 1, 3
        iextrp = 1  # required when NOOBS=2
        nowsta = 0
        m3d_block = "MM51.DAT       input     1  ! M3DDAT=3d.dat!    !END!\n"
        up_block = ""
        srf_line = "* SRFDAT=             *"
        geodat = "geo.dat"
    else:
        raise ValueError(mode)

    # BIAS/NSMTH length NZ
    bias = ", ".join(["0"] * NZ)
    nsmth = ", ".join(["2"] + ["4"] * (NZ - 1))
    fextr2 = ", ".join(["0"] * NZ)
    iuv = ", ".join(["1"] + ["0"] * (NZ - 1))

    # Groups 7-9 only if stations present
    grp7 = ""
    grp8 = ""
    grp9 = ""
    if nssta > 0:
        grp7 = f"""
INPUT GROUP: 7 -- Surface meteorological station parameters
--------------
! SS1  ='{SFC_NAME:<4s}'  {SFC_ID:6d}    {SFC_X:10.3f}  {SFC_Y:10.3f}  {STN_TZ:3d}  {ANEM_HT:5.0f}  !
!END!
"""
    if nusta > 0:
        grp8 = f"""
INPUT GROUP: 8 -- Upper air meteorological station parameters
--------------
! US1  ='{UP_NAME:<4s}'  {UP_ID:6d}  {SFC_X:10.3f}  {SFC_Y:10.3f}  {STN_TZ:3d}  !
!END!
"""
    # NPSTA=0: no group 9

    text = f"""CALMET.INP      2.1             Hour Start and End Times with Seconds
Tiny synthetic domain -- mode={mode}
12x12 @ 1km, UTM 19N, 2020-06-15 00-03 UTC
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
         F = UPPER CASE

    Number of upper air stations (NUSTA)  No default     ! NUSTA =  {nusta}  !
    Number of overwater met stations
                                 (NOWSTA) No default     ! NOWSTA =  {nowsta}  !

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

     Starting date:    Year   (IBYR)  --    No default   ! IBYR  =  {YEAR}  !
                       Month  (IBMO)  --    No default   ! IBMO  =  {MONTH}  !
                       Day    (IBDY)  --    No default   ! IBDY  =  {DAY}  !
     Starting time:    Hour   (IBHR)  --    No default   ! IBHR  =  {IBHR}  !
                       Second (IBSEC) --    No default   ! IBSEC =  0  !

     Ending date:      Year   (IEYR)  --    No default   ! IEYR  =  {YEAR}  !
                       Month  (IEMO)  --    No default   ! IEMO  =  {MONTH}  !
                       Day    (IEDY)  --    No default   ! IEDY  =  {DAY}  !
     Ending time:      Hour   (IEHR)  --    No default   ! IEHR  =  {IEHR}  !
                       Second (IESEC) --    No default   ! IESEC =  0  !

      UTC time zone         (ABTZ) -- No default       ! ABTZ= {ABTZ} !
         (character*8)

     Length of modeling time-step (seconds)
     (NSECDT)                        Default:3600     ! NSECDT =  3600  !

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
     (PMAP)                     Default: UTM    ! PMAP = {PMAP}  !

     False Easting and Northing (km) at the projection origin
     (Used only if PMAP= TTM, LCC, or LAZA)
     (FEAST)                    Default=0.0   ! FEAST  = 0.0 !
     (FNORTH)                   Default=0.0   ! FNORTH = 0.0 !

     UTM zone (1 to 60)
     (Used only if PMAP=UTM)
     (IUTMZN)                   No Default      ! IUTMZN =  {UTMZN}   !

     Hemisphere for UTM projection?
     (Used only if PMAP=UTM)
     (UTMHEM)                   Default: N    ! UTMHEM = {UTMHEM}  !
         N   - Northern hemisphere
         S   - Southern hemisphere

     Latitude and Longitude of N.pole or projection origin
     (Used only if PMAP= LCC, PS, EM, TTM, or LAZA)
     (RLAT0)                    No Default      * RLAT0= 0N  *
     (RLON0)                    No Default      * RLON0= 0E  *

     Matching parallels
     (Used only if PMAP= LCC or PS)
     (XLAT1)                    No Default      * XLAT1= 30N *
     (XLAT2)                    No Default      * XLAT2= 60N *

     DATUM for output coordinates
     (DATUM)                    Default: WGS-84    ! DATUM = {DATUM}  !

            No. X grid cells (NX)      No default     ! NX =   {NX}  !
            No. Y grid cells (NY)      No default     ! NY =   {NY}  !
     Grid spacing (DGRIDKM)            No default     ! DGRIDKM = {DGRIDKM} !
                                             Units: km
        X coordinate of SW corner
        of grid cell (1,1)
        (lower left corner of grid)
        X coordinate (XORIGKM)         No default     ! XORIGKM = {XORIGKM:.3f} !
        Y coordinate (YORIGKM)         No default     ! YORIGKM = {YORIGKM:.3f} !
                                             Units: km

     Number of vertical layers (NZ)    No default     ! NZ = {NZ} !

     Cell face heights
        vertical grid (ZFACE(NZ+1))    No defaults
                                             Units: m
        ! ZFACE = {", ".join(str(z) for z in ZFACE)}  !

!END!

-------------------------------------------------------------------------------

INPUT GROUP: 3 -- Output Options
--------------

     Disk Output
     -----------
     Save met. fields in unformatted
     output file ? (LSAVE)               Default: T     ! LSAVE = T !
     Type of unformatted output file:
     (IFORMO)                            Default: 1     ! IFORMO = 1 !
        1 = CALPUFF/CALGRID type file (CALMET.DAT)
        2 = MESOPUFF-II type file (PACOUT.DAT)

     Line printer output flag
     (LPRINT)                            Default: F     ! LPRINT = F !
     Print interval (hours)
     (IPRINF)                            Default: 1     ! IPRINF = 1 !

     ! IUVOUT = {iuv} !
     ! IWOUT  = {iuv} !
     ! ITOUT  = {iuv} !

     ! STABILITY = F !
     ! USTAR     = F !
     ! MONIN     = F !
     ! MIXHT     = F !
     ! WSTAR     = F !
     ! PRECIP    = F !
     ! SENSHEAT  = F !
     ! CONVZI    = F !

     Testing and debug print options
     ! LDB    = F   !
     ! NN1    = 1   !
     ! NN2    = 1   !
     ! LDBCST = F   !
     ! IOUTD  = 0   !
     ! NZPRN2 = 1   !
     ! IPR0   = 0   !
     ! IPR1   = 0   !
     ! IPR2   = 0   !
     ! IPR3   = 0   !
     ! IPR4   = 0   !
     ! IPR5   = 0   !
     ! IPR6   = 0   !
     ! IPR7   = 0   !
     ! IPR8   = 0   !

!END!

-------------------------------------------------------------------------------

INPUT GROUP: 4 -- Meteorological data options
--------------

    NO OBSERVATION MODE             (NOOBS)  Default: 0     ! NOOBS =  {noobs}   !

       Number of surface stations   (NSSTA)  No default     ! NSSTA =  {nssta}  !

       Number of precipitation stations
                                    (NPSTA)  No default     ! NPSTA =  {npsta}  !

       Gridded cloud fields:
                                   (ICLOUD)  Default: 0     ! ICLOUD =  {icloud}  !

       Surface meteorological data file format
                                   (IFORMS)  Default: 2     ! IFORMS =  2  !
       Precipitation data file format
                                   (IFORMP)  Default: 2     ! IFORMP =  2  !
       Cloud data file format
                                   (IFORMC)  Default: 2     ! IFORMC =  2  !

!END!

-------------------------------------------------------------------------------

INPUT GROUP: 5 -- Wind Field Options and Parameters
--------------

       Model (IWFCOD)                    Default: 1      ! IWFCOD =  1  !
       Compute Froude number adjustment
       effects? (IFRADJ)                 Default: 1      ! IFRADJ =  1  !
       Compute kinematic effects? (IKINE) Default: 0      ! IKINE  =  0  !
       Use O'Brien procedure for adjustment
       of the vertical velocity? (IOBR)  Default: 0      ! IOBR   =  0  !
       Compute slope flow effects? (ISLOPE) Default: 1    ! ISLOPE =  1  !

       Extrapolate surface wind observations
       to upper layers (IEXTRP)          Default: -4     ! IEXTRP = {iextrp}  !
       Extrapolate surface winds even if
       calm? (ICALM)                     Default: 0      ! ICALM  =  0  !

       Layer-dependent biases (BIAS(NZ))
                               ! BIAS =  {bias}  !

       Minimum distance from nearest upper air station
       to surface station (RMIN2)
                                             Default: 4.     ! RMIN2 = -1.0 !

       Use gridded prognostic wind field model
       output fields (IPROG)              Default: 0      ! IPROG =  {iprog}  !

       Timestep (seconds) of the prognostic
       model input data   (ISTEPPGS)         Default: 3600   ! ISTEPPGS =  3600   !

       Use coarse CALMET fields as initial guess fields (IGFMET)
                                             Default: 0      ! IGFMET =  0  !

       Use varying radius of influence       Default: F      ! LVARY =  T!

       Maximum radius of influence over land
       in the surface layer (RMAX1)          No default      ! RMAX1 = 50. !
       Maximum radius of influence over land
       aloft (RMAX2)                         No default      ! RMAX2 = 50. !
       Maximum radius of influence over water
       (RMAX3)                               No default      ! RMAX3 = 50. !

       Minimum radius of influence used in
       the wind field interpolation (RMIN)   Default: 0.1    ! RMIN = 0.1 !
       Radius of influence of terrain
       features (TERRAD)                     No default      ! TERRAD = 5. !

       Relative weighting surface (R1)       No default      ! R1 = 1. !
       Relative weighting aloft (R2)         No default      ! R2 = 1. !
       Relative weighting of the prognostic
       wind field data (RPROG)               No default      ! RPROG = 0.0 !

       Maximum acceptable divergence
       (DIVLIM)                              Default: 5.E-6  ! DIVLIM = 5.0E-6 !
       Maximum number of iterations in the
       divergence minimization procedure
       (NITER)                               Default: 50     ! NITER = 50 !
       Number of stations used in each layer
       for the interpolation of data to a
       grid point (NINTR2(NZ))
                               ! NINTR2 = {", ".join(["99"] * NZ)} !

       Critical Froude number (CRITFN)       Default: 1.0    ! CRITFN = 1.0 !
       Empirical factor controlling the
       influence of kinematic effects
       (ALPHA)                               Default: 0.1    ! ALPHA = 0.1 !

       Number of barriers to interpolation
       of the wind fields (NBAR)             Default: 0      ! NBAR = 0 !
       ! KBAR = {NZ} !
       ! XBBAR = 0.0 !
       ! YBBAR = 0.0 !
       ! XEBAR = 0.0 !
       ! YEBAR = 0.0 !

       Method used to compute surface
       temperatures (IDIOPT1)                Default: 0      ! IDIOPT1 = 0 !
       Surface met station to use for the
       surface temperature (ISURFT)
                                             No default      ! ISURFT = {"1" if nssta > 0 else "-2"} !

       Method used to compute domain-averaged
       temperature lapse rate (IDIOPT2)      Default: 0      ! IDIOPT2 = 0 !
       Upper air station to use for the
       domain-scale lapse rate (IUPT)
                                             No default      ! IUPT = {"1" if nusta > 0 else "-2"} !
       Depth (m) through which the domain-scale
       lapse rate is computed (ZUPT)         Default: 200.   ! ZUPT = 200. !

       Method used to compute domain-averaged
       wind components (IDIOPT3)             Default: 0      ! IDIOPT3 = 0 !
       Upper air station to use for the
       domain-scale winds (IUPWND)
                                             Default: -1     ! IUPWND = -1 !
       Bottom and top of layer (m) through
       which domain-scale winds are computed
       (ZUPWND(1), ZUPWND(2))
                                             Default: 1., 1000. ! ZUPWND = 1., 1000. !

       Selection of observed surface wind
       components (IDIOPT4)                  Default: 0      ! IDIOPT4 = 0 !
       Selection of observed upper air wind
       components (IDIOPT5)                  Default: 0      ! IDIOPT5 = 0 !

       ! LLBREZE = F !
       ! NBOX = 0 !
       ! XG1 = 0.0 !
       ! XG2 = 0.0 !
       ! YG1 = 0.0 !
       ! YG2 = 0.0 !
       ! XBCST = 0.0 !
       ! YBCST = 0.0 !
       ! XECST = 0.0 !
       ! YECST = 0.0 !
       ! NLB = 0 !
       ! METBXID = 0 !

       ! NSMTH = {nsmth} !
       ! FEXTR2 = {fextr2} !

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

{grp7}{grp8}{grp9}
"""
    path.write_text(text)


def main() -> None:
    SHARED.mkdir(parents=True, exist_ok=True)
    write_geo(SHARED / "geo.dat")
    write_surf(SHARED / "surf.dat")
    write_up(SHARED / "up.dat")
    write_3d(SHARED / "3d.dat")
    for mode in ("obs", "obs_model", "noobs"):
        d = ROOT / mode
        d.mkdir(parents=True, exist_ok=True)
        # symlink/copy shared inputs
        for f in ("geo.dat", "surf.dat", "up.dat", "3d.dat"):
            src = SHARED / f
            dst = d / f
            if dst.exists() or dst.is_symlink():
                dst.unlink()
            os.symlink(os.path.relpath(src, d), dst)
        make_inp(mode, d / "calmet.inp")
        # binary symlink
        bin_src = Path("/workspace/py-calmet/vendor/calmet-fortran/src/calmet.x")
        bin_dst = d / "calmet.x"
        if bin_dst.exists() or bin_dst.is_symlink():
            bin_dst.unlink()
        os.symlink(bin_src, bin_dst)
        print(f"Prepared {d}")
    print("Done. Shared files:", list(SHARED.iterdir()))


if __name__ == "__main__":
    main()
