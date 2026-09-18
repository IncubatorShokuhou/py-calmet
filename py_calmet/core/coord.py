"""Map projections used by CALMET (UTM / LCC / TM / PS / EM / LAZA).

INP-driven: PMAP, RLAT0/RLON0, XLAT1/XLAT2, FEAST/FNORTH, IUTMZN/UTMHEM, DATUM.
Pure NumPy; WGS84 ellipsoid. Good enough for solar/Coriolis and station remap —
not a bit-identical COORDLIB clone.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .met_utils import utm_to_latlon

A_WGS84 = 6378137.0
F_WGS84 = 1.0 / 298.257223563
E2 = F_WGS84 * (2.0 - F_WGS84)


@dataclass
class MapProjection:
    pmap: str = "UTM"
    iutmzn: int = 19
    utmhem: str = "N"
    rlat0: float = 0.0
    rlon0: float = 0.0
    xlat1: float = 0.0
    xlat2: float = 0.0
    feast: float = 0.0
    fnorth: float = 0.0
    datum: str = "WGS-84"

    @classmethod
    def from_inp(cls, inp) -> "MapProjection":
        def _f(key, default=0.0):
            try:
                return float(inp.get_float(key, default))
            except Exception:
                return float(default)

        def _parse_ll(raw, default=0.0):
            if raw is None:
                return float(default)
            s = str(raw).strip().upper().replace("'", "")
            if not s or s.startswith("*"):
                return float(default)
            # Accept 40.0N / 70.0W / -70.0 / 40N
            hemi = 1.0
            if s[-1] in "NS":
                hemi = 1.0 if s[-1] == "N" else -1.0
                s = s[:-1]
            elif s[-1] in "EW":
                hemi = 1.0 if s[-1] == "E" else -1.0
                s = s[:-1]
            try:
                return hemi * float(s)
            except Exception:
                return float(default)

        pmap = str(inp.get("PMAP", "UTM") or "UTM").strip().upper()
        return cls(
            pmap=pmap,
            iutmzn=int(inp.get_int("IUTMZN", 19)),
            utmhem=str(inp.get("UTMHEM", "N") or "N").strip().upper()[:1] or "N",
            rlat0=_parse_ll(inp.get("RLAT0"), _f("RLAT0", 0.0)),
            rlon0=_parse_ll(inp.get("RLON0"), _f("RLON0", 0.0)),
            xlat1=_parse_ll(inp.get("XLAT1"), _f("XLAT1", 0.0)),
            xlat2=_parse_ll(inp.get("XLAT2"), _f("XLAT2", 0.0)),
            feast=_f("FEAST", 0.0),
            fnorth=_f("FNORTH", 0.0),
            datum=str(inp.get("DATUM", "WGS-84") or "WGS-84").strip(),
        )


def _lcc_forward(lat, lon, lat0, lon0, lat1, lat2, feast=0.0, fnorth=0.0):
    """Lambert Conformal Conic → (x_m, y_m)."""
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    lat0 = np.deg2rad(lat0)
    lon0 = np.deg2rad(lon0)
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    if abs(lat1 - lat2) < 1e-12:
        n = np.sin(lat1)
    else:
        n = np.log(np.cos(lat1) / np.cos(lat2)) / np.log(
            np.tan(np.pi / 4 + lat2 / 2) / np.tan(np.pi / 4 + lat1 / 2)
        )
    f = (np.cos(lat1) * (np.tan(np.pi / 4 + lat1 / 2) ** n)) / n
    rho = A_WGS84 * f * (np.tan(np.pi / 4 + lat / 2) ** (-n))
    rho0 = A_WGS84 * f * (np.tan(np.pi / 4 + lat0 / 2) ** (-n))
    theta = n * (lon - lon0)
    x = feast + rho * np.sin(theta)
    y = fnorth + rho0 - rho * np.cos(theta)
    return float(x), float(y)


def _lcc_inverse(x, y, lat0, lon0, lat1, lat2, feast=0.0, fnorth=0.0):
    """Lambert Conformal Conic inverse → (lat, lon_east)."""
    x = float(x) - feast
    y = float(y) - fnorth
    lat0 = np.deg2rad(lat0)
    lon0 = np.deg2rad(lon0)
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
    if abs(lat1 - lat2) < 1e-12:
        n = np.sin(lat1)
    else:
        n = np.log(np.cos(lat1) / np.cos(lat2)) / np.log(
            np.tan(np.pi / 4 + lat2 / 2) / np.tan(np.pi / 4 + lat1 / 2)
        )
    f = (np.cos(lat1) * (np.tan(np.pi / 4 + lat1 / 2) ** n)) / n
    rho0 = A_WGS84 * f * (np.tan(np.pi / 4 + lat0 / 2) ** (-n))
    rho = np.sign(n) * np.sqrt(x * x + (rho0 - y) ** 2)
    theta = np.arctan2(x, rho0 - y)
    lat = 2.0 * np.arctan((A_WGS84 * f / rho) ** (1.0 / n)) - np.pi / 2
    lon = theta / n + lon0
    return float(np.rad2deg(lat)), float(np.rad2deg(lon))


def _tm_forward(lat, lon, lat0, lon0, feast=0.0, fnorth=0.0, k0=0.9996):
    """Transverse Mercator (approx) → (x_m, y_m)."""
    lat = np.deg2rad(lat)
    lon = np.deg2rad(lon)
    lon0 = np.deg2rad(lon0)
    lat0 = np.deg2rad(lat0)
    ep2 = E2 / (1.0 - E2)
    n = A_WGS84 / np.sqrt(1.0 - E2 * np.sin(lat) ** 2)
    t = np.tan(lat) ** 2
    c = ep2 * np.cos(lat) ** 2
    a = (lon - lon0) * np.cos(lat)
    m = A_WGS84 * (
        (1 - E2 / 4 - 3 * E2**2 / 64 - 5 * E2**3 / 256) * lat
        - (3 * E2 / 8 + 3 * E2**2 / 32 + 45 * E2**3 / 1024) * np.sin(2 * lat)
        + (15 * E2**2 / 256 + 45 * E2**3 / 1024) * np.sin(4 * lat)
        - (35 * E2**3 / 3072) * np.sin(6 * lat)
    )
    m0 = A_WGS84 * (
        (1 - E2 / 4 - 3 * E2**2 / 64 - 5 * E2**3 / 256) * lat0
        - (3 * E2 / 8 + 3 * E2**2 / 32 + 45 * E2**3 / 1024) * np.sin(2 * lat0)
        + (15 * E2**2 / 256 + 45 * E2**3 / 1024) * np.sin(4 * lat0)
        - (35 * E2**3 / 3072) * np.sin(6 * lat0)
    )
    x = feast + k0 * n * (
        a + (1 - t + c) * a**3 / 6 + (5 - 18 * t + t**2 + 72 * c - 58 * ep2) * a**5 / 120
    )
    y = fnorth + k0 * (
        m
        - m0
        + n
        * np.tan(lat)
        * (
            a**2 / 2
            + (5 - t + 9 * c + 4 * c**2) * a**4 / 24
            + (61 - 58 * t + t**2 + 600 * c - 330 * ep2) * a**6 / 720
        )
    )
    return float(x), float(y)


def project_ll_to_xy(lat: float, lon_east: float, proj: MapProjection) -> tuple[float, float]:
    """Geographic → projected meters."""
    pmap = proj.pmap.upper().strip()
    if pmap in ("UTM",):
        # Rough forward via zone central meridian TM
        lon0 = (proj.iutmzn - 1) * 6 - 180 + 3
        x, y = _tm_forward(lat, lon_east, 0.0, lon0, feast=500000.0, fnorth=0.0 if proj.utmhem != "S" else 1.0e7)
        return x, y
    if pmap in ("LCC",):
        lat1 = proj.xlat1 if abs(proj.xlat1) > 1e-9 else proj.rlat0 - 5.0
        lat2 = proj.xlat2 if abs(proj.xlat2) > 1e-9 else proj.rlat0 + 5.0
        return _lcc_forward(lat, lon_east, proj.rlat0, proj.rlon0, lat1, lat2, proj.feast, proj.fnorth)
    if pmap in ("TTM", "TM", "LAZA", "EM", "PS"):
        # Approximate with TM about (rlat0, rlon0)
        return _tm_forward(lat, lon_east, proj.rlat0, proj.rlon0, proj.feast, proj.fnorth)
    return _tm_forward(lat, lon_east, proj.rlat0, proj.rlon0, proj.feast, proj.fnorth)


def project_xy_to_ll(x_m: float, y_m: float, proj: MapProjection) -> tuple[float, float]:
    """Projected meters → (lat, lon_east)."""
    pmap = proj.pmap.upper().strip()
    if pmap in ("UTM",):
        return utm_to_latlon(x_m, y_m, proj.iutmzn, northern=proj.utmhem != "S")
    if pmap in ("LCC",):
        lat1 = proj.xlat1 if abs(proj.xlat1) > 1e-9 else proj.rlat0 - 5.0
        lat2 = proj.xlat2 if abs(proj.xlat2) > 1e-9 else proj.rlat0 + 5.0
        return _lcc_inverse(x_m, y_m, proj.rlat0, proj.rlon0, lat1, lat2, proj.feast, proj.fnorth)
    # Fallback: treat as TM about origin — use finite-difference inverse via UTM-like
    # For PS/EM/LAZA/TTM: approximate with LCC single-parallel (lat1=lat2=rlat0)
    return _lcc_inverse(
        x_m, y_m, proj.rlat0, proj.rlon0, proj.rlat0, proj.rlat0, proj.feast, proj.fnorth
    )


def domain_center_latlon(
    xorig_km: float,
    yorig_km: float,
    nx: int,
    ny: int,
    dgrid_km: float,
    proj: MapProjection,
) -> tuple[float, float]:
    """Lat/lon at domain center from INP-driven projection."""
    xc = (xorig_km + 0.5 * nx * dgrid_km) * 1000.0
    yc = (yorig_km + 0.5 * ny * dgrid_km) * 1000.0
    return project_xy_to_ll(xc, yc, proj)
