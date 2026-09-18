#!/usr/bin/env python3
"""Convert a WRF wrfout NetCDF file to CALMET 3D.DAT (dataset 2.1, IOUTMM5=92).

Reads U/V/T/P/PB/PH/PHB/QVAPOR/HGT/LU_INDEX/U10/V10/T2/Q2/PSFC/SWDOWN/GLW
with netCDF4, destaggers winds, builds hydrostatic heights, and writes the
text 3D.DAT format used by hrrr2calmet / tiny-domain goldens.
"""
from __future__ import annotations

import argparse
import math
from datetime import datetime
from pathlib import Path

import numpy as np

try:
    from netCDF4 import Dataset
except ImportError as exc:  # pragma: no cover
    raise SystemExit("netCDF4 is required: pip install netCDF4") from exc

R_D = 287.05
CP = 1004.0
P0 = 1.0e5
G = 9.81


def uv_to_wdws(u, v):
    ws = np.hypot(u, v)
    wd = (np.degrees(np.arctan2(-u, -v)) + 360.0) % 360.0
    return wd, ws


def rh_from_qv_t_p(qv, tk, p_pa):
    """Approximate RH (%) from water-vapor mixing ratio, T, P."""
    qv = np.maximum(qv, 1e-12)
    e = qv * p_pa / (0.622 + qv)
    tc = tk - 273.15
    # Lowe (1977) saturation vapor pressure (mb) → Pa
    a = np.array(
        [
            6.107799961,
            4.436518521e-1,
            1.428945805e-2,
            2.650648471e-4,
            3.031240396e-6,
            2.034080948e-8,
            6.136820929e-11,
        ]
    )
    es_mb = a[0] + tc * (
        a[1]
        + tc
        * (a[2] + tc * (a[3] + tc * (a[4] + tc * (a[5] + tc * a[6]))))
    )
    es = np.maximum(es_mb, 0.01) * 100.0
    return np.clip(100.0 * e / es, 1.0, 100.0)


def destagger_u(u):
    """U staggered on west_east_stag → mass points."""
    return 0.5 * (u[..., :, :-1] + u[..., :, 1:])


def destagger_v(v):
    """V staggered on south_north_stag → mass points."""
    return 0.5 * (v[..., :-1, :] + v[..., 1:, :])


def destagger_w_full(w):
    """W on bottom_top_stag → mass levels."""
    return 0.5 * (w[..., :-1, :, :] + w[..., 1:, :, :])


def mm5_lu_from_wrf(lu_index: np.ndarray) -> np.ndarray:
    """Map USGS 24-category LU_INDEX toward CALMET/MM5-ish codes (best-effort)."""
    # Keep urban=1-ish as 10, water=16 → 55 (CALMET water), else 20 (default land)
    out = np.full(lu_index.shape, 20, dtype=np.int32)
    out = np.where(np.isin(lu_index, [1]), 10, out)  # urban
    out = np.where(np.isin(lu_index, [16]), 55, out)  # water
    return out


def write_3d_dat(
    path: Path,
    *,
    times: list[datetime],
    x0_km: float,
    y0_km: float,
    dx_km: float,
    lat: np.ndarray,  # [nj, ni]
    lon: np.ndarray,
    elev: np.ndarray,
    lu: np.ndarray,
    sigma: np.ndarray,
    # surface [nt, nj, ni]
    psfc_mb,
    rain_cm,
    sc,
    radsw,
    radlw,
    t2,
    q2_gkg,
    wd10,
    ws10,
    sst,
    # upper [nt, nk, nj, ni]
    pres_mb,
    z_msl,
    tempk,
    wd,
    ws,
    w,
    rh,
    vapmr,
    lat0: float,
    lon0: float,
    lat1: float,
    lat2: float,
    comment: str = "WRF wrfout → 3D.DAT via scripts/wrfout_to_3d.py",
    i0: int = 1,
    j0: int = 1,
):
    nt, nj, ni = t2.shape
    nk = len(sigma)
    lines = []
    lines.append(f"{'3D.DAT':<16s}{'2.1':<16s}{'Header Structure with Comment Lines':<64s}")
    lines.append("   1")
    lines.append(f"{comment:<132s}")
    lines.append("  1  1  0  0  0")  # ioutw ioutq → format 92
    lines.append(
        f"{'LCC':<4s}{lat0:9.4f}{lon0:10.4f}{lat1:7.2f}{lat2:7.2f}"
        f"{x0_km:10.3f}{y0_km:10.3f}{dx_km:8.3f}{ni:4d}{nj:3d}{nk:3d}"
    )
    lines.append("  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0")
    ymdh0 = times[0].strftime("%Y%m%d%H")
    lines.append(f"{ymdh0}{nt:5d}{ni:4d}{nj:4d}{nk:4d}")
    lines.append(
        f"{i0:4d}{j0:4d}{i0 + ni - 1:4d}{j0 + nj - 1:4d}"
        f"{1:4d}{nk:4d}{float(lon.min()):10.4f}{float(lon.max()):10.4f}"
        f"{float(lat.min()):9.4f}{float(lat.max()):9.4f}"
    )
    for s in sigma:
        lines.append(f"{s:6.3f}")
    for jj in range(nj):
        for ii in range(ni):
            lines.append(
                f"{i0 + ii:4d}{j0 + jj:4d}{lat[jj, ii]:9.4f}{lon[jj, ii]:10.4f}"
                f"{int(round(elev[jj, ii])):5d}{int(lu[jj, ii]):3d}"
                f" {-999:9.4f}{-999:10.4f}{-999:5d}"
            )

    for t, stamp_dt in enumerate(times):
        stamp = stamp_dt.strftime("%Y%m%d%H")
        for jj in range(nj):
            for ii in range(ni):
                lines.append(
                    f"{stamp}{i0 + ii:03d}{j0 + jj:03d}"
                    f"{psfc_mb[t, jj, ii]:7.1f}{rain_cm[t, jj, ii]:5.2f}{int(sc[t, jj, ii]):2d}"
                    f"{radsw[t, jj, ii]:8.1f}{radlw[t, jj, ii]:8.1f}"
                    f"{t2[t, jj, ii]:8.1f}{q2_gkg[t, jj, ii]:8.2f}"
                    f"{wd10[t, jj, ii]:8.1f}{ws10[t, jj, ii]:8.1f}{sst[t, jj, ii]:8.1f}"
                )
                for k in range(nk):
                    lines.append(
                        f"{int(round(pres_mb[t, k, jj, ii])):4d}"
                        f"{int(round(z_msl[t, k, jj, ii])):6d}"
                        f"{tempk[t, k, jj, ii]:6.1f}"
                        f"{int(round(wd[t, k, jj, ii])) % 360:4d}"
                        f"{ws[t, k, jj, ii]:5.1f}"
                        f"{w[t, k, jj, ii]:6.2f}"
                        f"{int(round(rh[t, k, jj, ii])):3d}"
                        f"{vapmr[t, k, jj, ii]:5.2f}"
                    )
    path.write_text("\n".join(lines) + "\n")


def convert_wrfout(
    wrfout: Path,
    out_3d: Path,
    *,
    i0: int = 30,
    j0: int = 25,
    ni: int = 14,
    nj: int = 14,
    nk_out: int = 10,
    time_indices: list[int] | None = None,
) -> dict:
    ds = Dataset(wrfout)
    ntime = len(ds.dimensions["Time"])
    if time_indices is None:
        time_indices = list(range(min(ntime, 4)))
    times = []
    for ti in time_indices:
        raw = b"".join(ds.variables["Times"][ti]).decode().strip()
        times.append(datetime.strptime(raw, "%Y-%m-%d_%H:%M:%S"))

    # subset mass-point window (0-based inclusive start)
    is0, js0 = i0, j0
    is1, js1 = i0 + ni, j0 + nj

    lat = np.asarray(ds.variables["XLAT"][0, js0:js1, is0:is1], dtype=np.float64)
    lon = np.asarray(ds.variables["XLONG"][0, js0:js1, is0:is1], dtype=np.float64)
    elev = np.asarray(ds.variables["HGT"][0, js0:js1, is0:is1], dtype=np.float64)
    lu = mm5_lu_from_wrf(np.asarray(ds.variables["LU_INDEX"][0, js0:js1, is0:is1]))

    # Project mass-point lat/lon to UTM so CALMET GEO↔3D lat/lon match.
    try:
        from pyproj import Transformer
    except ImportError as exc:  # pragma: no cover
        raise SystemExit("pyproj required for UTM mapping: pip install pyproj") from exc

    # Choose UTM zone from mean longitude
    lon_m = float(lon.mean())
    zone = int((lon_m + 180) // 6) + 1
    epsg = 32600 + zone if float(lat.mean()) >= 0 else 32700 + zone
    to_utm = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    x_m, y_m = to_utm.transform(lon, lat)
    # Mean grid spacing from adjacent cell centers
    dx_km = float(np.mean(np.diff(x_m, axis=1)) / 1000.0)
    dy_km = float(np.mean(np.diff(y_m, axis=0)) / 1000.0)
    dx_km = 0.5 * (abs(dx_km) + abs(dy_km))
    # SW corner of cell (0,0): center minus half spacing
    x0_km = float(x_m[0, 0] / 1000.0 - 0.5 * dx_km)
    y0_km = float(y_m[0, 0] / 1000.0 - 0.5 * dx_km)

    lat0 = float(getattr(ds, "CEN_LAT", float(lat.mean())))
    lon0 = float(getattr(ds, "CEN_LON", float(lon.mean())))
    lat1 = float(getattr(ds, "TRUELAT1", lat0) or lat0)
    lat2 = float(getattr(ds, "TRUELAT2", lat0) or lat0)
    if abs(lat1) < 1e-6 and abs(lat2) < 1e-6:
        lat1 = lat2 = lat0

    nt = len(time_indices)
    nk_full = len(ds.dimensions["bottom_top"])
    # pick lowest nk_out levels (near surface first in WRF)
    k_sel = list(range(min(nk_out, nk_full)))
    nk = len(k_sel)

    # sigma proxy from mean pressure ratio
    p_full = np.asarray(ds.variables["P"][time_indices[0], :, js0:js1, is0:is1]) + np.asarray(
        ds.variables["PB"][time_indices[0], :, js0:js1, is0:is1]
    )
    psfc0 = np.asarray(ds.variables["PSFC"][time_indices[0], js0:js1, is0:is1])
    sigma = np.array([(p_full[k] / psfc0).mean() for k in k_sel], dtype=np.float64)

    psfc_mb = np.zeros((nt, nj, ni))
    rain_cm = np.zeros((nt, nj, ni))
    sc = np.zeros((nt, nj, ni), dtype=np.int32)
    radsw = np.zeros((nt, nj, ni))
    radlw = np.zeros((nt, nj, ni))
    t2 = np.zeros((nt, nj, ni))
    q2_gkg = np.zeros((nt, nj, ni))
    wd10 = np.zeros((nt, nj, ni))
    ws10 = np.zeros((nt, nj, ni))
    sst = np.zeros((nt, nj, ni))

    pres_mb = np.zeros((nt, nk, nj, ni))
    z_msl = np.zeros((nt, nk, nj, ni))
    tempk = np.zeros((nt, nk, nj, ni))
    wd = np.zeros((nt, nk, nj, ni))
    ws = np.zeros((nt, nk, nj, ni))
    w = np.zeros((nt, nk, nj, ni))
    rh = np.zeros((nt, nk, nj, ni))
    vapmr = np.zeros((nt, nk, nj, ni))

    has_sw = "SWDOWN" in ds.variables
    has_glw = "GLW" in ds.variables
    has_rainc = "RAINC" in ds.variables
    has_rainnc = "RAINNC" in ds.variables

    for it, ti in enumerate(time_indices):
        u10 = np.asarray(ds.variables["U10"][ti, js0:js1, is0:is1])
        v10 = np.asarray(ds.variables["V10"][ti, js0:js1, is0:is1])
        wd10[it], ws10[it] = uv_to_wdws(u10, v10)
        t2[it] = np.asarray(ds.variables["T2"][ti, js0:js1, is0:is1])
        q2 = np.asarray(ds.variables["Q2"][ti, js0:js1, is0:is1])
        q2_gkg[it] = q2 * 1000.0  # kg/kg → g/kg
        psfc_mb[it] = np.asarray(ds.variables["PSFC"][ti, js0:js1, is0:is1]) / 100.0
        sst[it] = t2[it]
        if has_sw:
            radsw[it] = np.asarray(ds.variables["SWDOWN"][ti, js0:js1, is0:is1])
        if has_glw:
            radlw[it] = np.asarray(ds.variables["GLW"][ti, js0:js1, is0:is1])
        if has_rainc or has_rainnc:
            rain = np.zeros((nj, ni))
            if has_rainc:
                rain += np.asarray(ds.variables["RAINC"][ti, js0:js1, is0:is1])
            if has_rainnc:
                rain += np.asarray(ds.variables["RAINNC"][ti, js0:js1, is0:is1])
            rain_cm[it] = rain / 10.0  # mm → cm

        # 3D fields
        U = destagger_u(np.asarray(ds.variables["U"][ti, :, js0:js1, is0 : is1 + 1]))
        V = destagger_v(np.asarray(ds.variables["V"][ti, :, js0 : js1 + 1, is0:is1]))
        Wf = destagger_w_full(np.asarray(ds.variables["W"][ti, :, js0:js1, is0:is1]))
        P = np.asarray(ds.variables["P"][ti]) + np.asarray(ds.variables["PB"][ti])
        P = P[:, js0:js1, is0:is1]
        PH = (np.asarray(ds.variables["PH"][ti]) + np.asarray(ds.variables["PHB"][ti])) / G
        # mass-level height
        z_full = 0.5 * (PH[:-1, js0:js1, is0:is1] + PH[1:, js0:js1, is0:is1])
        theta = np.asarray(ds.variables["T"][ti, :, js0:js1, is0:is1]) + 300.0
        tk = theta * (P / P0) ** (R_D / CP)
        qv = np.asarray(ds.variables["QVAPOR"][ti, :, js0:js1, is0:is1])

        for kk, k in enumerate(k_sel):
            uu, vv = U[k], V[k]
            wd_k, ws_k = uv_to_wdws(uu, vv)
            wd[it, kk] = wd_k
            ws[it, kk] = ws_k
            w[it, kk] = Wf[k]
            pres_mb[it, kk] = P[k] / 100.0
            z_msl[it, kk] = z_full[k]
            tempk[it, kk] = tk[k]
            rh[it, kk] = rh_from_qv_t_p(qv[k], tk[k], P[k])
            vapmr[it, kk] = qv[k] * 1000.0  # g/kg

    ds.close()

    write_3d_dat(
        out_3d,
        times=times,
        x0_km=x0_km,
        y0_km=y0_km,
        dx_km=dx_km,
        lat=lat,
        lon=lon,
        elev=elev,
        lu=lu,
        sigma=sigma,
        psfc_mb=psfc_mb,
        rain_cm=rain_cm,
        sc=sc,
        radsw=radsw,
        radlw=radlw,
        t2=t2,
        q2_gkg=q2_gkg,
        wd10=wd10,
        ws10=ws10,
        sst=sst,
        pres_mb=pres_mb,
        z_msl=z_msl,
        tempk=tempk,
        wd=wd,
        ws=ws,
        w=w,
        rh=rh,
        vapmr=vapmr,
        lat0=lat0,
        lon0=lon0,
        lat1=lat1,
        lat2=lat2,
        i0=1,
        j0=1,
    )
    return {
        "times": times,
        "ni": ni,
        "nj": nj,
        "nk": nk,
        "dx_km": dx_km,
        "lat": lat,
        "lon": lon,
        "elev": elev,
        "lu": lu,
        "x0_km": x0_km,
        "y0_km": y0_km,
        "utm_zone": zone,
        "utm_epsg": epsg,
        "lat0": lat0,
        "lon0": lon0,
        "t2": t2,
        "wd10": wd10,
        "ws10": ws10,
        "psfc_mb": psfc_mb,
        "q2_gkg": q2_gkg,
        "radsw": radsw,
        "pres_mb": pres_mb,
        "z_msl": z_msl,
        "tempk": tempk,
        "wd": wd,
        "ws": ws,
        "rh": rh,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("wrfout", type=Path)
    p.add_argument("-o", "--output", type=Path, required=True)
    p.add_argument("--i0", type=int, default=30)
    p.add_argument("--j0", type=int, default=25)
    p.add_argument("--ni", type=int, default=14)
    p.add_argument("--nj", type=int, default=14)
    p.add_argument("--nk", type=int, default=10)
    args = p.parse_args()
    meta = convert_wrfout(
        args.wrfout,
        args.output,
        i0=args.i0,
        j0=args.j0,
        ni=args.ni,
        nj=args.nj,
        nk_out=args.nk,
    )
    print(
        f"Wrote {args.output}  nt={len(meta['times'])} "
        f"grid={meta['ni']}x{meta['nj']}x{meta['nk']} dx={meta['dx_km']} km"
    )


if __name__ == "__main__":
    main()
