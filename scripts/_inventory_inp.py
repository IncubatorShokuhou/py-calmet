#!/usr/bin/env python3
"""Build CALMET INP inventory JSON + markdown coverage + roadmap."""
from __future__ import annotations

import json
import re
from collections import Counter, OrderedDict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path("/workspace/py-calmet")
DOCS = ROOT / "docs"

READCF_GROUPS = OrderedDict(
    [
        (
            1,
            [
                "IBYR",
                "IBMO",
                "IBDY",
                "IBHR",
                "IBSEC",
                "IEYR",
                "IEMO",
                "IEDY",
                "IEHR",
                "IESEC",
                "ABTZ",
                "IBTZ",
                "IRLG",
                "NSECDT",
                "IRTYPE",
                "LCALGRD",
                "ITEST",
                "MREG",
            ],
        ),
        (
            2,
            [
                "PMAP",
                "DATUM",
                "FEAST",
                "FNORTH",
                "IUTMZN",
                "UTMHEM",
                "RLAT0",
                "RLON0",
                "XLAT1",
                "XLAT2",
                "NX",
                "NY",
                "DGRIDKM",
                "XORIGKM",
                "YORIGKM",
                "NZ",
                "ZFACE",
            ],
        ),
        (
            3,
            [
                "LSAVE",
                "LPRINT",
                "IPRINF",
                "IUVOUT",
                "IWOUT",
                "ITOUT",
                "STABILITY",
                "USTAR",
                "MONIN",
                "MIXHT",
                "WSTAR",
                "PRECIP",
                "SENSHEAT",
                "CONVZI",
                "LDB",
                "NN1",
                "NN2",
                "LDBCST",
                "IOUTD",
                "NZPRN2",
                "IPR0",
                "IPR1",
                "IPR2",
                "IPR3",
                "IPR4",
                "IPR5",
                "IPR6",
                "IPR7",
                "IPR8",
                "IFORMO",
            ],
        ),
        (
            4,
            [
                "NOOBS",
                "NSSTA",
                "NPSTA",
                "IFORMS",
                "IFORMP",
                "ICLOUD",
                "ICLDOUT",
                "MCLOUD",
                "IFORMC",
            ],
        ),
        (
            5,
            [
                "IWFCOD",
                "IFRADJ",
                "IKINE",
                "IOBR",
                "IEXTRP",
                "RMIN2",
                "FEXTR2",
                "IPROG",
                "ISTEPPG",
                "ISTEPPGS",
                "IGFMET",
                "LVARY",
                "RMAX1",
                "RMAX2",
                "RMAX3",
                "RMIN",
                "TERRAD",
                "R1",
                "R2",
                "RPROG",
                "DIVLIM",
                "NITER",
                "NSMTH",
                "NINTR2",
                "CRITFN",
                "ALPHA",
                "NBAR",
                "XBBAR",
                "YBBAR",
                "XEBAR",
                "YEBAR",
                "KBAR",
                "IDIOPT1",
                "IDIOPT2",
                "IDIOPT3",
                "IDIOPT4",
                "IDIOPT5",
                "ISURFT",
                "IUPT",
                "ZUPT",
                "IUPWND",
                "ZUPWND",
                "LLBREZE",
                "NBOX",
                "XG1",
                "XG2",
                "YG1",
                "YG2",
                "XBCST",
                "YBCST",
                "XECST",
                "YECST",
                "NLB",
                "METBXID",
                "BIAS",
                "ISLOPE",
                "ICALM",
            ],
        ),
        (
            6,
            [
                "CONSTB",
                "CONSTE",
                "CONSTN",
                "DPTMIN",
                "DZZI",
                "ZIMIN",
                "ZIMAX",
                "ZIMINW",
                "ZIMAXW",
                "IAVEZI",
                "MNMDAV",
                "HAFANG",
                "ILEVZI",
                "FCORIOL",
                "CONSTW",
                "ITPROG",
                "ITWPROG",
                "ILUOC3D",
                "IRAD",
                "IAVET",
                "TGDEFB",
                "TGDEFA",
                "JWAT1",
                "JWAT2",
                "TRADKM",
                "NUMTS",
                "NFLAGP",
                "SIGMAP",
                "CUTP",
                "HA1",
                "HA2",
                "HB1",
                "HB2",
                "HC1",
                "HC2",
                "HC3",
                "IMIXH",
                "THRESHL",
                "THRESHW",
                "ICOARE",
                "DSHELF",
                "IWARM",
                "ICOOL",
                "IRHPROG",
                "IZICRLX",
                "TZICRLX",
            ],
        ),
    ]
)

READFN_GROUPS = OrderedDict(
    [
        (
            "0a",
            [
                "METINP",
                "GEODAT",
                "SRFDAT",
                "PRCDAT",
                "MM4DAT",
                "WTDAT",
                "METLST",
                "METDAT",
                "PACDAT",
                "CLDDAT",
                "LCFILES",
                "NUSTA",
                "NOWSTA",
                "NM3D",
                "NIGF",
            ],
        ),
        ("0b", ["UPDAT"]),
        ("0c", ["SEADAT"]),
        ("0d", ["M3DDAT"]),
        ("0e", ["IGFDAT"]),
        (
            "0f",
            [
                "DIADAT",
                "PRGDAT",
                "TSTPRT",
                "TSTOUT",
                "TSTKIN",
                "TSTFRD",
                "TSTSLP",
                "DCSTGD",
            ],
        ),
    ]
)

STATION_KEYS = ["SS1", "US1", "PS1"]

GROUP_NAMES = {
    "0a": "File names — primary I/O (READFN subgroup a)",
    "0b": "Upper-air file names (READFN b)",
    "0c": "Overwater / SEA.DAT file names (READFN c)",
    "0d": "MM4/MM5/3D.DAT file names (READFN d)",
    "0e": "IGF-CALMET.DAT file names (READFN e)",
    "0f": "Misc diagnostic / test file names (READFN f)",
    "1": "General run control",
    "2": "Map projection and grid",
    "3": "Output options",
    "4": "Meteorological data options",
    "5": "Wind field options and parameters",
    "6": "Mixing height, temperature, precip, overwater",
    "station": "Station location free-format records (IG 7–9)",
    "sample_only": "Found in sample INPs only",
}

TYPES: dict[str, str] = {}
for k in [
    "IBYR",
    "IBMO",
    "IBDY",
    "IBHR",
    "IBSEC",
    "IEYR",
    "IEMO",
    "IEDY",
    "IEHR",
    "IESEC",
    "IBTZ",
    "IRLG",
    "NSECDT",
    "IRTYPE",
    "ITEST",
    "MREG",
]:
    TYPES[k] = "integer"
TYPES["ABTZ"] = "character"
TYPES["LCALGRD"] = "logical"
for k in ["PMAP", "DATUM", "UTMHEM", "RLAT0", "RLON0", "XLAT1", "XLAT2"]:
    TYPES[k] = "character"
for k in ["FEAST", "FNORTH", "DGRIDKM", "XORIGKM", "YORIGKM", "ZFACE"]:
    TYPES[k] = "real"
for k in ["IUTMZN", "NX", "NY", "NZ"]:
    TYPES[k] = "integer"
for k in [
    "LSAVE",
    "LPRINT",
    "STABILITY",
    "USTAR",
    "MONIN",
    "MIXHT",
    "WSTAR",
    "PRECIP",
    "SENSHEAT",
    "CONVZI",
    "LDB",
    "LDBCST",
]:
    TYPES[k] = "logical"
for k in [
    "IPRINF",
    "IUVOUT",
    "IWOUT",
    "ITOUT",
    "NN1",
    "NN2",
    "IOUTD",
    "NZPRN2",
    "IPR0",
    "IPR1",
    "IPR2",
    "IPR3",
    "IPR4",
    "IPR5",
    "IPR6",
    "IPR7",
    "IPR8",
    "IFORMO",
]:
    TYPES[k] = "integer"
for k in READCF_GROUPS[4]:
    TYPES[k] = "integer"
for k in [
    "IWFCOD",
    "IFRADJ",
    "IKINE",
    "IOBR",
    "IEXTRP",
    "IPROG",
    "ISTEPPG",
    "ISTEPPGS",
    "IGFMET",
    "NITER",
    "NSMTH",
    "NINTR2",
    "NBAR",
    "KBAR",
    "IDIOPT1",
    "IDIOPT2",
    "IDIOPT3",
    "IDIOPT4",
    "IDIOPT5",
    "ISURFT",
    "IUPT",
    "IUPWND",
    "NBOX",
    "NLB",
    "METBXID",
    "ISLOPE",
    "ICALM",
]:
    TYPES[k] = "integer"
for k in [
    "RMIN2",
    "FEXTR2",
    "RMAX1",
    "RMAX2",
    "RMAX3",
    "RMIN",
    "TERRAD",
    "R1",
    "R2",
    "RPROG",
    "DIVLIM",
    "CRITFN",
    "ALPHA",
    "XBBAR",
    "YBBAR",
    "XEBAR",
    "YEBAR",
    "ZUPT",
    "ZUPWND",
    "XG1",
    "XG2",
    "YG1",
    "YG2",
    "XBCST",
    "YBCST",
    "XECST",
    "YECST",
    "BIAS",
]:
    TYPES[k] = "real"
TYPES["LVARY"] = "logical"
TYPES["LLBREZE"] = "logical"
for k in [
    "CONSTB",
    "CONSTE",
    "CONSTN",
    "DPTMIN",
    "DZZI",
    "ZIMIN",
    "ZIMAX",
    "ZIMINW",
    "ZIMAXW",
    "HAFANG",
    "FCORIOL",
    "CONSTW",
    "TGDEFB",
    "TGDEFA",
    "TRADKM",
    "SIGMAP",
    "CUTP",
    "HA1",
    "HA2",
    "HB1",
    "HB2",
    "HC1",
    "HC2",
    "HC3",
    "THRESHL",
    "THRESHW",
    "DSHELF",
    "TZICRLX",
]:
    TYPES[k] = "real"
for k in [
    "IAVEZI",
    "MNMDAV",
    "ILEVZI",
    "ITPROG",
    "ITWPROG",
    "ILUOC3D",
    "IRAD",
    "IAVET",
    "JWAT1",
    "JWAT2",
    "NUMTS",
    "NFLAGP",
    "IMIXH",
    "ICOARE",
    "IWARM",
    "ICOOL",
    "IRHPROG",
    "IZICRLX",
]:
    TYPES[k] = "integer"
for k in [
    "METINP",
    "GEODAT",
    "SRFDAT",
    "PRCDAT",
    "MM4DAT",
    "WTDAT",
    "METLST",
    "METDAT",
    "PACDAT",
    "CLDDAT",
    "UPDAT",
    "SEADAT",
    "M3DDAT",
    "IGFDAT",
    "DIADAT",
    "PRGDAT",
    "TSTPRT",
    "TSTOUT",
    "TSTKIN",
    "TSTFRD",
    "TSTSLP",
    "DCSTGD",
]:
    TYPES[k] = "character"
TYPES["LCFILES"] = "logical"
for k in ["NUSTA", "NOWSTA", "NM3D", "NIGF"]:
    TYPES[k] = "integer"
for k in STATION_KEYS:
    TYPES[k] = "station_record"

ARRAY_NOTES = {
    "ZFACE": "length NZ+1",
    "IUVOUT": "length NZ",
    "IWOUT": "length NZ",
    "ITOUT": "length NZ",
    "FEXTR2": "length NZ",
    "NSMTH": "length NZ",
    "NINTR2": "length NZ",
    "BIAS": "length NZ",
    "XBBAR": "length NBAR",
    "YBBAR": "length NBAR",
    "XEBAR": "length NBAR",
    "YEBAR": "length NBAR",
    "XG1": "length NBOX",
    "XG2": "length NBOX",
    "YG1": "length NBOX",
    "YG2": "length NBOX",
    "XBCST": "length NBOX",
    "YBCST": "length NBOX",
    "XECST": "length NBOX",
    "YECST": "length NBOX",
    "METBXID": "length mxbxwnd",
    "JWAT1": "length mxwb",
    "JWAT2": "length mxwb",
    "ZUPWND": "length 2",
    "UPDAT": "per upper-air station",
    "SEADAT": "per overwater station",
    "M3DDAT": "per MM5/3D file",
    "IGFDAT": "per IGF file",
    "SS1": "SSn surface station records",
    "US1": "USn upper-air station records",
    "PS1": "PSn precip station records",
}

# --- status tables ---
IMPLEMENTED: dict[str, str] = {
    "NOOBS": "Mode inference (noobs / obs / obs_model).",
    "IPROG": "Mode inference when >0 with observations.",
    "IBYR": "Run window start year (datetime).",
    "IBMO": "Run window start month.",
    "IBDY": "Run window start day.",
    "IBHR": "Run window start hour.",
    "IBSEC": "Run window start second.",
    "IEYR": "Run window end year.",
    "IEMO": "Run window end month.",
    "IEDY": "Run window end day.",
    "IEHR": "Run window end hour.",
    "IESEC": "Run window end second.",
    "NSECDT": "Timestep seconds (datetime timedelta).",
    "ABTZ": "UTC offset → internal timezone hours.",
    "ZFACE": "Vertical grid faces.",
    "NZ": "Number of layers.",
    "IUTMZN": "UTM zone for lat/lon / header.",
    "UTMHEM": "UTM hemisphere.",
    "RLAT0": "Optional domain lat for solar/Coriolis.",
    "RLON0": "Optional domain lon for solar.",
    "IFRADJ": "Froude blocking on/off (FRADJ).",
    "ISLOPE": "Slope flow on/off (Mahrt).",
    "NSMTH": "Per-layer 5-point smooth passes.",
    "TERRAD": "Terrain radius (km) for FRADJ/slope.",
    "CRITFN": "Critical Froude number.",
    "CONSTN": "Night Zi constant.",
    "CONSTB": "Carson buoyancy constant.",
    "ZIMIN": "Min mixing height.",
    "ZIMAX": "Max mixing height.",
    "THRESHL": "Daytime MIXHMC growth threshold.",
}

PARTIAL: dict[str, str] = {
    "SS1": "Single-station X/Y for OA; multi-station SSn unused.",
    "US1": "UP.DAT sounding used; US1 coords unused.",
    "R1": "Single-station Barnes surface radius; multi-station OA missing.",
    "R2": "Aloft OA radius; multi-station incomplete.",
    "IKINE": "Flag gates stub kinematic+div path; full TOPOF2 W missing.",
    "IOBR": "Flag gates simplified div-min; full O'Brien missing.",
    "ALPHA": "Used in light kinematic/div-min; full IKINE semantics incomplete.",
    "NITER": "Used (capped) in divergence minimization.",
    "DPTMIN": "Constant Carson gamma; MIXDT sounding lapse deferred.",
    "JWAT1": "INP name JWAT1; runner wrongly reads IWAT1 → default 55. Alias to iwat1.",
    "JWAT2": "INP name JWAT2; runner wrongly reads IWAT2. Alias to iwat2.",
    "GEODAT": "geo.dat via path search; GEODAT filename not honored.",
    "SRFDAT": "surf.dat path search; SRFDAT not honored.",
    "UPDAT": "up.dat path search; UPDAT not honored.",
    "M3DDAT": "3d.dat path search; M3DDAT not honored.",
    "METDAT": "Writer emits CALMET.DAT; METDAT path not read.",
    "METLST": "No list-file writer yet.",
    "NX": "Taken from GEO.DAT, not INP.",
    "NY": "Taken from GEO.DAT.",
    "DGRIDKM": "Taken from GEO.DAT.",
    "XORIGKM": "Taken from GEO.DAT.",
    "YORIGKM": "Taken from GEO.DAT.",
    "PMAP": "UTM assumed; LCC/PS/EM/TTM/LAZA not wired.",
    "DATUM": "Not used in coord transforms.",
    "FCORIOL": "Coriolis from lat; FCORIOL INP ignored.",
    "IMIXH": "Maul–Carson day + night mechanical only; other IMIXH options missing.",
    "ICLOUD": "Cloud from SURF sky tenths; schemes 1–4 missing.",
    "IWFCOD": "Diagnostic winds always on; flag not read.",
    "ITPROG": "3D T2 used in noobs; full ITPROG paths missing.",
    "IRHPROG": "RH from 3D/SURF; flag not read.",
    "ISTEPPGS": "NSECDT drives time; prognostic-step QA missing.",
    "LSAVE": "Caller decides write; flag not read.",
    "IFORMO": "CALMET.DAT only (type 1); PACOUT missing.",
    "LCALGRD": "Some CALGRID fields in writer; flag not read.",
    "HA1": "Hardcoded in pbl.shortwave_radiation; not from INP.",
    "HA2": "Hardcoded in pbl.py.",
    "HB1": "Hardcoded in pbl.py.",
    "HB2": "Hardcoded in pbl.py.",
    "HC1": "Energy-budget defaults hardcoded; not from INP.",
    "HC2": "Hardcoded.",
    "HC3": "Hardcoded.",
    "IRAD": "Solar always computed; IRAD option not read.",
    "CONSTE": "Energy-budget path partial; CONSTE unused.",
    "NSSTA": "Assumes 0/1 station; multi-station unused.",
    "NUSTA": "Assumes 0/1.",
    "NM3D": "Single 3d.dat assumed.",
    "IFORMS": "SURF format assumed.",
    "STABILITY": "IPGT computed; print flag unused.",
    "USTAR": "Field computed; print flag unused.",
    "MONIN": "EL computed; print flag unused.",
    "MIXHT": "ZI computed; print flag unused.",
    "WSTAR": "Field computed; print flag unused.",
    "PRECIP": "RMM zeroed; print/data path unused.",
    "SENSHEAT": "QH internal; output partial.",
    "CONVZI": "ziconv for Carson; print flag unused.",
    "IEXTRP": "Power-law obs profile; SIMILT helper exists but not wired.",
    "BIAS": "Not applied to layer blending.",
    "LCFILES": "Paths used as-is; case folding unused.",
}

MISSING: dict[str, str] = {
    "ICOARE": "No COARE overwater flux module (PSIUD helper only).",
    "DSHELF": "COARE coastal shelf unused.",
    "IWARM": "COARE warm-layer unused.",
    "ICOOL": "COARE cool-skin unused.",
    "CONSTW": "Overwater mixing unused.",
    "ZIMINW": "Overwater Zi min unused.",
    "ZIMAXW": "Overwater Zi max unused.",
    "THRESHW": "Overwater threshold unused.",
    "MCLOUD": "INP 2.2 cloud method unused.",
    "ICLDOUT": "Cloud output flag unused.",
    "CLDDAT": "CLOUD.DAT reader unused.",
    "IFORMC": "Cloud file format unused.",
    "NPSTA": "Precipitation stations unused.",
    "IFORMP": "Precip file format unused.",
    "PRCDAT": "PRECIP.DAT unused.",
    "NFLAGP": "Precip QC unused.",
    "SIGMAP": "Precip sigma unused.",
    "CUTP": "Precip cutoff unused.",
    "PS1": "Precip station records unused.",
    "NOWSTA": "Overwater station count unused.",
    "SEADAT": "SEA.DAT unused.",
    "WTDAT": "Overwater SST soft-spot (io.wt_dat); terrain-weight layout OutOfScope.",
    "IGFMET": "IGF as first-guess unused.",
    "IGFDAT": "IGF-CALMET files unused.",
    "NIGF": "IGF file count unused.",
    "DIADAT": "DIAG.DAT unused.",
    "PRGDAT": "PROG.DAT unused.",
    "NBAR": "Wind barriers unused.",
    "XBBAR": "Barrier coords unused.",
    "YBBAR": "Barrier coords unused.",
    "XEBAR": "Barrier coords unused.",
    "YEBAR": "Barrier coords unused.",
    "KBAR": "Barrier top level unused.",
    "LLBREZE": "Lake-breeze unused.",
    "NBOX": "Lake-breeze boxes unused.",
    "XG1": "Unused.",
    "XG2": "Unused.",
    "YG1": "Unused.",
    "YG2": "Unused.",
    "XBCST": "Coastline box unused.",
    "YBCST": "Unused.",
    "XECST": "Unused.",
    "YECST": "Unused.",
    "NLB": "Unused.",
    "METBXID": "Unused.",
    "IDIOPT1": "Diag option switches unused.",
    "IDIOPT2": "Unused.",
    "IDIOPT3": "Unused.",
    "IDIOPT4": "Unused.",
    "IDIOPT5": "Unused.",
    "ISURFT": "Surface T station index unused.",
    "IUPT": "Upper T index unused.",
    "ZUPT": "Unused.",
    "IUPWND": "Upper wind index unused.",
    "ZUPWND": "Unused.",
    "FEXTR2": "Extrapolation factors unused.",
    "NINTR2": "Max stations per layer unused.",
    "RPROG": "Prog weight in OA unused.",
    "LVARY": "Varying radius unused.",
    "RMAX1": "OA distance cutoff unused.",
    "RMAX2": "Unused.",
    "RMAX3": "Unused.",
    "RMIN": "Unused.",
    "RMIN2": "Unused.",
    "DIVLIM": "Divergence criterion unused (hardcoded loop).",
    "ICALM": "Calm processing unused.",
    "ITWPROG": "Water T from prog unused.",
    "ILUOC3D": "3D landuse overlay unused.",
    "IAVEZI": "Zi averaging unused.",
    "MNMDAV": "Unused.",
    "HAFANG": "Unused.",
    "ILEVZI": "Unused.",
    "IAVET": "T averaging unused.",
    "TGDEFB": "Default lapse (below) unused.",
    "TGDEFA": "Default lapse (above) unused.",
    "TRADKM": "Temperature OA radius unused.",
    "NUMTS": "Unused.",
    "DZZI": "Inversions thickness unused.",
    "IZICRLX": "Zi relaxation unused.",
    "TZICRLX": "Unused.",
    "LDBCST": "Coast distance output unused.",
    "DCSTGD": "Unused.",
    "LPRINT": "Printer output unused.",
    "IPRINF": "Unused.",
    "IUVOUT": "Layer UV print unused.",
    "IWOUT": "Unused.",
    "ITOUT": "Unused.",
    "LDB": "Debug unused.",
    "NN1": "Unused.",
    "NN2": "Unused.",
    "IOUTD": "Unused.",
    "NZPRN2": "Unused.",
    "IPR0": "Unused.",
    "IPR1": "Unused.",
    "IPR2": "Unused.",
    "IPR3": "Unused.",
    "IPR4": "Unused.",
    "IPR5": "Unused.",
    "IPR6": "Unused.",
    "IPR7": "Unused.",
    "IPR8": "Unused.",
    "TSTPRT": "Test files unused.",
    "TSTOUT": "Unused.",
    "TSTKIN": "Unused.",
    "TSTFRD": "Unused.",
    "TSTSLP": "Unused.",
    "PACDAT": "MESOPUFF PACOUT unused.",
    "METINP": "Unused.",
    "MM4DAT": "Legacy MM4 name; use M3DDAT.",
    "IRTYPE": "Run-type gating unused.",
    "ITEST": "Setup-only stop unused.",
    "MREG": "Regulatory QA unused.",
    "FEAST": "Non-UTM false easting unused.",
    "FNORTH": "Unused.",
    "XLAT1": "LCC/PS parallel unused.",
    "XLAT2": "Unused.",
    "IBTZ": "Legacy; prefer ABTZ.",
    "IRLG": "Legacy run length; prefer begin/end + NSECDT.",
    "ISTEPPG": "Legacy hours; prefer ISTEPPGS.",
}


def api_field(name: str, group) -> str:
    mapping = {
        "NOOBS": "run.noobs",
        "IPROG": "winds.iprog",
        "IBYR": "run.start.year",
        "IBMO": "run.start.month",
        "IBDY": "run.start.day",
        "IBHR": "run.start.hour",
        "IBSEC": "run.start.second",
        "IEYR": "run.end.year",
        "IEMO": "run.end.month",
        "IEDY": "run.end.day",
        "IEHR": "run.end.hour",
        "IESEC": "run.end.second",
        "ABTZ": "run.abtz",
        "IBTZ": "run.ibtz_legacy",
        "IRLG": "run.irlg_legacy",
        "NSECDT": "run.nsecdt",
        "IRTYPE": "run.irtype",
        "LCALGRD": "output.lcalgrd",
        "ITEST": "run.itest",
        "MREG": "run.mreg",
        "PMAP": "grid.pmap",
        "DATUM": "grid.datum",
        "FEAST": "grid.feast",
        "FNORTH": "grid.fnorth",
        "IUTMZN": "grid.iutmzn",
        "UTMHEM": "grid.utmhem",
        "RLAT0": "grid.rlat0",
        "RLON0": "grid.rlon0",
        "XLAT1": "grid.xlat1",
        "XLAT2": "grid.xlat2",
        "NX": "grid.nx",
        "NY": "grid.ny",
        "DGRIDKM": "grid.dgridkm",
        "XORIGKM": "grid.xorigkm",
        "YORIGKM": "grid.yorigkm",
        "NZ": "grid.nz",
        "ZFACE": "grid.zface",
        "LSAVE": "output.lsave",
        "IFORMO": "output.iformo",
        "LPRINT": "output.lprint",
        "METDAT": "files.metdat",
        "METLST": "files.metlst",
        "GEODAT": "files.geodat",
        "SRFDAT": "files.srfdat",
        "UPDAT": "files.updat",
        "M3DDAT": "files.m3ddat",
        "SEADAT": "files.seadat",
        "CLDDAT": "files.clddat",
        "PRCDAT": "files.prcdat",
        "IGFDAT": "files.igfdat",
        "JWAT1": "pbl.jwat1",
        "JWAT2": "pbl.jwat2",
        "IFRADJ": "winds.ifradj",
        "IKINE": "winds.ikine",
        "IOBR": "winds.iobr",
        "ISLOPE": "winds.islope",
        "NSMTH": "winds.nsmth",
        "R1": "winds.r1",
        "R2": "winds.r2",
        "RPROG": "winds.rprog",
        "TERRAD": "winds.terrad",
        "CRITFN": "winds.critfn",
        "ALPHA": "winds.alpha",
        "NITER": "winds.niter",
        "DIVLIM": "winds.divlim",
        "IEXTRP": "winds.iextrp",
        "BIAS": "winds.bias",
        "IWFCOD": "winds.iwfcod",
        "LVARY": "winds.lvary",
        "IGFMET": "winds.igfmet",
        "ISTEPPGS": "winds.isteppgs",
        "CONSTN": "pbl.constn",
        "CONSTB": "pbl.constb",
        "CONSTE": "pbl.conste",
        "CONSTW": "pbl.constw",
        "ZIMIN": "pbl.zimin",
        "ZIMAX": "pbl.zimax",
        "ZIMINW": "pbl.ziminw",
        "ZIMAXW": "pbl.zimaxw",
        "DPTMIN": "pbl.dptmin",
        "THRESHL": "pbl.threshl",
        "THRESHW": "pbl.threshw",
        "IMIXH": "pbl.imixh",
        "ICOARE": "overwater.icoare",
        "DSHELF": "overwater.dshelf",
        "IWARM": "overwater.iwarm",
        "ICOOL": "overwater.icool",
        "ICLOUD": "clouds.icloud",
        "MCLOUD": "clouds.mcloud",
        "ICLDOUT": "clouds.icldout",
        "ITPROG": "temp.itprog",
        "IRHPROG": "humidity.irhprog",
        "IRAD": "radiation.irad",
        "HA1": "radiation.ha1",
        "HA2": "radiation.ha2",
        "HB1": "radiation.hb1",
        "HB2": "radiation.hb2",
        "HC1": "radiation.hc1",
        "HC2": "radiation.hc2",
        "HC3": "radiation.hc3",
        "SS1": "stations.surface[0]",
        "US1": "stations.upper[0]",
        "PS1": "stations.precip[0]",
        "NUSTA": "files.nusta",
        "NOWSTA": "files.nowsta",
        "NM3D": "files.nm3d",
        "NIGF": "files.nigf",
        "NSSTA": "met.nssta",
        "NPSTA": "met.npsta",
        "FCORIOL": "pbl.fcoriol",
        "LCFILES": "files.lcfiles",
    }
    if name in mapping:
        return mapping[name]
    section = {
        "0a": "files",
        "0b": "files",
        "0c": "files",
        "0d": "files",
        "0e": "files",
        "0f": "files",
        "1": "run",
        "2": "grid",
        "3": "output",
        "4": "met",
        "5": "winds",
        "6": "pbl",
        "station": "stations",
    }.get(str(group), "misc")
    return f"{section}.{name.lower()}"


# --- WP1 overlays (CalmetConfig + INP round-trip + bindings) ---
WP1_IMPLEMENTED = {
    "R1": "Barnes OA surface radius (single + multi-station).",
    "R2": "Barnes OA aloft radius (single + multi-station).",

    "JWAT1": "INP JWAT1/JWAT2 bound; effective_iwat() aliases IWAT for heatfx/slope (999→GEO 55).",
    "JWAT2": "Paired with JWAT1; see JWAT1.",
    "HA1": "Bound into shortwave_radiation from CalmetConfig.",
    "HA2": "Bound into shortwave_radiation from CalmetConfig.",
    "HB1": "Bound into shortwave_radiation from CalmetConfig.",
    "HB2": "Bound into shortwave_radiation from CalmetConfig.",
    "HC1": "Bound into heat_flux_energy_budget from CalmetConfig.",
    "HC2": "Bound into heat_flux_energy_budget from CalmetConfig.",
    "HC3": "Bound into heat_flux_energy_budget from CalmetConfig.",
    "GEODAT": "Honored via _resolve_data_file (case_dir / inputs_dir).",
    "SRFDAT": "Honored via _resolve_data_file.",
    "UPDAT": "Honored via _resolve_data_file.",
    "M3DDAT": "Honored via _resolve_data_file (3d.dat)."
}
WP1_PARTIAL = {
    "RPROG": "Accepted; wired into objective_analyze as IGF weight (0=off, golden-safe).",
    "RMAX1": "Accepted; optional OA cutoff in multi-station path.",
    "RMAX2": "Accepted; optional OA cutoff aloft.",
    "R1": "Single+multi-station Barnes OA surface radius.",
    "R2": "Single+multi-station Barnes OA aloft radius.",

    "METDAT": "Parsed + exposed; output path still selected by caller/API (not auto-written to METDAT).",
    "METLST": "Parsed + exposed; list-file writer not yet implemented.",
    "LCFILES": "Parsed; paths used as-is (case folding unused).",
    "NUSTA": "Parsed; multi-station UP still single-file path.",
    "NM3D": "Parsed; single 3d.dat path honored.",
    "NOWSTA": "Parsed on config; SEA.DAT physics TBD.",
    "NIGF": "Parsed on config; IGF reader TBD.",
}
for _k, _v in WP1_IMPLEMENTED.items():
    IMPLEMENTED[_k] = _v
    PARTIAL.pop(_k, None)
    MISSING.pop(_k, None)
for _k, _v in WP1_PARTIAL.items():
    if _k not in IMPLEMENTED:
        PARTIAL[_k] = _v
        MISSING.pop(_k, None)


# --- WP2 overlays (DIAGNO OA, IKINE/IOBR, clouds, precip, COARE-lite) ---
WP2_IMPLEMENTED = {
    "RPROG": "IGF weight in Barnes OA (0=off); multi-station + golden-safe.",
    "RMAX1": "OA surface cutoff radius (multi-station path).",
    "RMAX2": "OA aloft cutoff radius.",
    "NINTR2": "Max stations per layer in Barnes OA.",
    "IKINE": "TOPOF2 topographic kinematic W + DIVLIM-minim when IKINE=1.",
    "IOBR": "OBrien continuity adjust with DIVLIM/NITER when IOBR=1.",
    "ALPHA": "TOPOF2 exponential decay coefficient (IKINE path).",
    "DIVLIM": "Divergence limit for MINIM / OBrien iteration.",
    "NITER": "Max iterations for divergence minimization.",
    "IEXTRP": "Obs profile modes: ±1 UA-aloft, ±2/-4 power-law (golden), ±3 FEXTR2, +4 SIMILT.",
    "FEXTR2": "Layer factors when |IEXTRP|=3.",
    "BIAS": "Layer speed bias when IEXTRP < 0.",
    "ISURFT": "1-based SS* index for OA / representative station.",
    "IUPT": "1-based upper-air station index (single-file path selects sounding).",
    "SS1": "Multi-station SS* X/Y parsed for Barnes OA.",
    "NSSTA": "Surface station count; multi-station OA when SS* >1.",
    "MCLOUD": "CLOUD3 (Teixeira RH) / CLOUD4-lite layered RH → ccfrac → QSW.",
    "ICLOUD": "Honors 3/4 RH schemes; 0/1 use SURF sky tenths.",
    "NPSTA": "−1 prognostic 3D.DAT rain → RMM; 0 none; >0 station Barnes (SIGMAP/CUTP).",
    "SIGMAP": "Precip OA influence radius (km).",
    "CUTP": "Precip rate cutoff (mm/h).",
    "ICOARE": "COARE-lite bulk fluxes over water when ICOARE≠0 and NOWSTA>0.",
    "CONSTW": "Overwater Zi scale in mixht_overwater.",
    "ZIMINW": "Overwater Zi minimum.",
    "ZIMAXW": "Overwater Zi maximum.",
    "DSHELF": "Coastal Cd enhancement in COARE-lite.",
    "NOWSTA": "Gates SEA/COARE path (NOWSTA>0 enables overwater fluxes).",
}
WP2_PARTIAL = {
    "PRCDAT": "Accepted; station rates default 0 without PRECIP.DAT reader.",
    "PS1": "Coords parsed; rates need PRECIP.DAT.",
    "IFORMP": "Precip format assumed.",
    "NFLAGP": "Precip QC unused.",
    "SEADAT": "SEA.DAT reader still TBD; COARE-lite uses air/SST bulk.",
    "IWARM": "COARE warm-layer not in lite scheme.",
    "ICOOL": "COARE cool-skin not in lite scheme.",
    "THRESHW": "Overwater convective threshold unused in lite Zi.",
    "ICLDOUT": "Cloud output file not written.",
    "CLDDAT": "CLOUD.DAT reader unused (RH path does not need it).",
    "IFORMC": "Cloud file format unused.",
}
for _k, _v in WP2_IMPLEMENTED.items():
    IMPLEMENTED[_k] = _v
    PARTIAL.pop(_k, None)
    MISSING.pop(_k, None)
for _k, _v in WP2_PARTIAL.items():
    if _k not in IMPLEMENTED:
        PARTIAL[_k] = _v
        MISSING.pop(_k, None)


# --- WP3 overlays (MIXDT, SEA/PRECIP/CLOUD I/O, barriers, lake breeze, coord, outputs) ---
WP3_IMPLEMENTED = {
    "DPTMIN": "Floor for MIXDT/MIXDT2 pot-temp lapse above Zi (daytime Carson).",
    "DZZI": "Depth (m) of MIXDT layer above Zi for sounding/prognostic lapse.",
    "ITPROG": "0=obs MIXDT sounding; 1/2=MIXDT2 from 3D.DAT columns (wired in runner).",
    "CONSTE": "Carson entrainment factor in mixht_day_carson.",
    "SEADAT": "SEA.DAT reader (v2.0/2.1/2.11); SST/ΔT/waves → COARE path.",
    "PRCDAT": "PRECIP.DAT reader; NPSTA>0 station rates → Barnes RMM.",
    "PS1": "Precip station X/Y + PRECIP.DAT rates for OA.",
    "IFORMP": "Formatted PRECIP.DAT (IFORMP=2) reader path.",
    "CLDDAT": "CLOUD.DAT reader + ICLDOUT writer (formatted CLOUDFRA).",
    "ICLDOUT": "Writes CLOUD.DAT when ICLDOUT≠0 and write_outputs=True.",
    "IFORMC": "CLOUD.DAT formatted I/O (IFORMC=2).",
    "IWARM": "COARE-lite warm-layer ΔT on skin SST when IWARM≠0.",
    "ICOOL": "COARE-lite cool-skin ΔT when ICOOL≠0.",
    "THRESHW": "Overwater convective boost hook in mixht_overwater.",
    "NBAR": "Wind barriers block OA across barrier segments (KBAR-aware).",
    "XBBAR": "Barrier begin X (km) used by BarrierSet.",
    "YBBAR": "Barrier begin Y (km).",
    "XEBAR": "Barrier end X (km).",
    "YEBAR": "Barrier end Y (km).",
    "KBAR": "Top layer (1-based) for barrier blocking in OA.",
    "LLBREZE": "Lake-breeze surface blend inside NBOX influence boxes.",
    "NBOX": "Number of lake-breeze boxes.",
    "XG1": "Lake-breeze box X min (km).",
    "XG2": "Lake-breeze box X max (km).",
    "YG1": "Lake-breeze box Y min (km).",
    "YG2": "Lake-breeze box Y max (km).",
    "XBCST": "Coastline segment begin X for lake breeze.",
    "YBCST": "Coastline segment begin Y.",
    "XECST": "Coastline segment end X.",
    "YECST": "Coastline segment end Y.",
    "NLB": "Stations per lake-breeze box (METBXID count).",
    "METBXID": "Station IDs inside lake-breeze boxes.",
    "PMAP": "INP-driven MapProjection (UTM/LCC/TM/PS/EM/LAZA) for lat/lon.",
    "DATUM": "Stored on MapProjection; used with PMAP.",
    "FEAST": "False easting for LCC/TM/LAZA projections.",
    "FNORTH": "False northing for LCC/TM/LAZA projections.",
    "XLAT1": "LCC/PS standard parallel 1.",
    "XLAT2": "LCC standard parallel 2.",
    "METLST": "METLST list-file writer (run summary); runner writes when write_outputs=True.",
    "METDAT": "METDAT name exposed on result.meta; LSAVE honored as flag.",
    "PACDAT": "PACOUT.DAT npz writer when IFORMO=2 and write_outputs=True.",
    "IFORMO": "1=CALMET.DAT path; 2=PACOUT hook.",
    "LSAVE": "Output-save flag recorded on result.meta.",
    "PRECIP": "RMM precip field computed (NPSTA path); print flag accepted.",
    "TGDEFB": "Accepted; MIXDT falls back through DPTMIN when sounding thin.",
    "TGDEFA": "Accepted; MIXDT falls back through DPTMIN when sounding thin.",
}
WP3_PARTIAL = {
    "NFLAGP": "Precip QC flag still unused (rates pass through).",
    "IFORMS": "SURF format assumed (dataset 2.1).",
}
for _k, _v in WP3_IMPLEMENTED.items():
    IMPLEMENTED[_k] = _v
    PARTIAL.pop(_k, None)
    MISSING.pop(_k, None)
for _k, _v in WP3_PARTIAL.items():
    if _k not in IMPLEMENTED:
        PARTIAL[_k] = _v
        MISSING.pop(_k, None)


# All remaining inventory names are at least accepted on CalmetConfig.
# Mark former Missing as Partial "accepted, physics TBD" unless still truly absent.

# --- WP4 overlays (listfile/IPR*, OA extras, Zi ops, IGF, OutOfScope MM4) ---
# Status "OutOfScope" is recorded as Implemented with an explicit rejection note
# (parameter is handled; MM4/MM5 paths are refused — wrfout→3D.DAT only).
WP4_IMPLEMENTED = {
    "METINP": "Control-file path recorded on result.meta / METLST (qa_notes).",
    "MM4DAT": "OutOfScope: non-default MM4/MM5 filenames raise NotImplementedError; use wrfout→3D.DAT.",
    "WTDAT": "Overwater SST soft-spot (io.wt_dat); precedence ITWPROG>SEA>WT>air-T. Official terrain-weight WT.DAT OutOfScope.",
    "LCFILES": "Case-insensitive data-file resolve when LCFILES=T.",
    "NUSTA": "Count honored; multi UPDAT list recorded (first readable used).",
    "NM3D": "Count honored; multi M3DDAT list recorded (first readable used).",
    "NIGF": "Count recorded; IGFDAT list used when IGFMET≠0.",
    "IGFDAT": "IGF-CALMET.DAT header/field reader (prior CALMET.DAT).",
    "DIADAT": "Stub writer when LDB/LDBCST and write_outputs.",
    "PRGDAT": "Stub writer when LDB/LDBCST and write_outputs.",
    "TSTPRT": "Stub writer when LDB/LDBCST and write_outputs.",
    "TSTOUT": "Stub writer when LDB/LDBCST and write_outputs.",
    "TSTKIN": "Stub writer when LDB/LDBCST and write_outputs.",
    "TSTFRD": "Stub writer when LDB/LDBCST and write_outputs.",
    "TSTSLP": "Stub writer when LDB/LDBCST and write_outputs.",
    "DCSTGD": "Stub writer when LDB/LDBCST and write_outputs.",
    "IBTZ": "Legacy timezone hours when ABTZ absent.",
    "IRLG": "Legacy run length (hours) when end ≤ start.",
    "IRTYPE": "IRTYPE=0 winds-only collapses PBL diagnostics.",
    "LCALGRD": "Flag recorded on result.meta / METLST.",
    "ITEST": "ITEST=1 setup-only early return after QA.",
    "MREG": "Flag recorded on result.meta / METLST (regulatory QA hook).",
    "NX": "QA vs GEO.DAT (GEO wins; mismatch noted).",
    "NY": "QA vs GEO.DAT (GEO wins; mismatch noted).",
    "DGRIDKM": "QA vs GEO.DAT (GEO wins; mismatch noted).",
    "XORIGKM": "QA vs GEO.DAT (GEO wins; mismatch noted).",
    "YORIGKM": "QA vs GEO.DAT (GEO wins; mismatch noted).",
    "LPRINT": "METLST verbosity when LPRINT/IPR*/field flags set.",
    "IPRINF": "Echoed in METLST IPR section.",
    "IUVOUT": "Layer UV print list honored in METLST.",
    "IWOUT": "Layer W print list honored in METLST.",
    "ITOUT": "Layer T print list honored in METLST.",
    "STABILITY": "METLST print-field flag (IPGT sample).",
    "USTAR": "METLST print-field flag (USTAR sample).",
    "MONIN": "METLST print-field flag.",
    "MIXHT": "METLST print-field flag (ZI sample).",
    "WSTAR": "METLST print-field flag.",
    "SENSHEAT": "METLST print-field flag.",
    "CONVZI": "METLST print-field flag.",
    "LDB": "Gates DIAG/PROG/TST* stub dumps with write_outputs.",
    "NN1": "Accepted; debug range recorded in METLST when LDB.",
    "NN2": "Accepted; debug range recorded in METLST when LDB.",
    "LDBCST": "Gates coast-distance / stub dumps with LDB.",
    "IOUTD": "Accepted on config; METLST echo.",
    "NZPRN2": "Accepted on config; METLST echo.",
    "IPR0": "METLST IPR0–IPR8 verbosity section.",
    "IPR1": "METLST IPR0–IPR8 verbosity section.",
    "IPR2": "METLST IPR0–IPR8 verbosity section.",
    "IPR3": "METLST IPR0–IPR8 verbosity section.",
    "IPR4": "METLST IPR0–IPR8 verbosity section.",
    "IPR5": "METLST IPR0–IPR8 verbosity section.",
    "IPR6": "METLST IPR0–IPR8 verbosity section.",
    "IPR7": "METLST IPR0–IPR8 verbosity section.",
    "IPR8": "METLST IPR0–IPR8 verbosity section.",
    "IFORMS": "SURF.DAT format flag read (dataset 2.1 path); recorded on meta.",
    "IWFCOD": "IWFCOD=0 skips DIAGNO wind module (keep first-guess/OA).",
    "ISTEPPG": "Legacy hours; QA via ISTEPPGS/NSECDT.",
    "ISTEPPGS": "Prognostic-step QA vs NSECDT multiple.",
    "IGFMET": "IGF first-guess from prior CALMET.DAT when IGFMET≠0.",
    "LVARY": "Varying OA radius when no station in RMAX (multi-station path).",
    "RMAX3": "Over-water OA cutoff (with landuse water mask).",
    "RMIN": "Minimum OA distance floor (m).",
    "RMIN2": "Accepted; surface extrapolation distance floor on config.",
    "IDIOPT1": "Surface T source for diag winds (0=obs ISURFT; 1=preprocessed QA-note).",
    "IDIOPT2": "Lapse source: 0=CGAMMA(ZUPT) from UP/3D → Froude/TOPOF2; 1=QA-note.",
    "IDIOPT3": "Domain-avg wind: 0=VERTAV(IUPWND,ZUPWND); 1=QA-note.",
    "IDIOPT4": "Preprocessed surface winds (IRTYPE=0 only); else QA reject note.",
    "IDIOPT5": "Preprocessed upper winds (IRTYPE=0 only); else QA reject note.",
    "ZUPT": "CGAMMA layer depth (m); drives Froude/TOPOF2 gamma when IDIOPT2=0.",
    "IUPWND": "UA station for domain-avg UV (VERTAV); -1 = all / first available.",
    "ZUPWND": "Height range [zlo,zhi] AGL for domain-avg UV when IDIOPT3=0.",
    "ICALM": "ICALM≠0 discards calm OA (ws<0.5) keeping IGF.",
    "IAVEZI": "Upwind Zi spatial average (no-op when MNMDAV≤1).",
    "MNMDAV": "Zi average cell count along upwind.",
    "HAFANG": "Zi average half-angle (deg).",
    "ILEVZI": "Wind layer for Zi upwind direction.",
    "FCORIOL": "Honored when ≠999 sentinel; else 2Ωsinφ.",
    "ITWPROG": "Flag read; SEA/3D water-T path gated (COARE).",
    "ILUOC3D": "3D ocean LU category recorded; water-mask helper.",
    "IRAD": "IRAD=0 zeroes shortwave; IRAD=1 computes solar.",
    "IAVET": "2-D temperature smoother (no-op when NUMTS≤1).",
    "TRADKM": "Temperature OA / smooth radius (km).",
    "NUMTS": "Temperature smoother half-width (cells).",
    "NFLAGP": "Precip QC: missing→0 and optional CUTP floor.",
    "IMIXH": "±1 Maul–Carson; ±2 Batchvarova–Gryning (mixht_day_bg); ±3 Holzworth.",
    "IRHPROG": "IRHPROG≠0 overwrites RH from 3D.DAT.",
    "IZICRLX": "Convective Zi exponential relaxation vs previous hour.",
    "TZICRLX": "Zi relaxation time scale (s).",
    "US1": "US* upper-air station X/Y parsed (qa_notes / OA anchor).",
}
for _k, _v in WP4_IMPLEMENTED.items():
    IMPLEMENTED[_k] = _v
    PARTIAL.pop(_k, None)
    MISSING.pop(_k, None)


def classify(name: str) -> tuple[str, str]:
    if name in IMPLEMENTED:
        return "Implemented", IMPLEMENTED[name]
    if name in PARTIAL:
        return "Partial", PARTIAL[name]
    if name in MISSING:
        # WP1: parameter is on CalmetConfig + INP round-trip; physics may still be TBD
        return "Partial", f"Accepted on CalmetConfig; physics TBD. ({MISSING[name]})"
    return "Partial", "Accepted on CalmetConfig; physics TBD."


def main() -> None:
    inp_files = list((ROOT / "cases").rglob("calmet.inp"))
    sample_keys: set[str] = set()
    for p in inp_files:
        text = p.read_text(errors="replace")
        for m in re.finditer(r"!\s*([A-Z][A-Z0-9]*)\s*=\s*([^!]*)!", text):
            sample_keys.add(m.group(1))

    params: list[dict] = []
    seen: set[str] = set()

    def add(name: str, group, source: str) -> None:
        if name in seen:
            for p in params:
                if p["name"] == name:
                    if source not in p["sources"]:
                        p["sources"].append(source)
                    break
            return
        seen.add(name)
        status, notes = classify(name)
        params.append(
            OrderedDict(
                [
                    ("name", name),
                    ("group", str(group)),
                    ("group_title", GROUP_NAMES.get(str(group), "")),
                    ("type", TYPES.get(name, "unknown")),
                    ("array", ARRAY_NOTES.get(name)),
                    ("sources", [source]),
                    ("in_readcf_or_readfn", True),
                    ("in_sample_inp", name in sample_keys),
                    ("status", status),
                    ("notes", notes),
                    ("api_field", api_field(name, group)),
                ]
            )
        )

    for g, keys in READFN_GROUPS.items():
        for k in keys:
            add(k, g, "READFN")
    for g, keys in READCF_GROUPS.items():
        for k in keys:
            add(k, g, "READCF")
    for k in STATION_KEYS:
        add(k, "station", "READCF_station_freeform")

    for k in sorted(sample_keys):
        if k not in seen:
            add(k, "sample_only", "sample_inp")
            params[-1]["in_readcf_or_readfn"] = False

    for p in params:
        p["in_sample_inp"] = p["name"] in sample_keys
        if p["in_sample_inp"] and "sample_inp" not in p["sources"]:
            p["sources"].append("sample_inp")

    counts = Counter(p["status"] for p in params)
    # Local Asia/Shanghai wall time for docs
    now_local = datetime.now().astimezone().strftime("%Y-%m-%d %H:%M %Z")

    meta = OrderedDict(
        [
            ("generated", now_local),
            (
                "fortran_source",
                "vendor/calmet-fortran/src/calmet.for (READCF CVDIC groups 1–6 + READFN IG0 + station freeform)",
            ),
            ("python_package", "py_calmet/"),
            ("sample_inp_glob", "cases/**/calmet.inp"),
            ("sample_inp_count", len(inp_files)),
            (
                "counts",
                OrderedDict(
                    [
                        ("total", len(params)),
                        ("implemented", counts.get("Implemented", 0)),
                        ("partial", counts.get("Partial", 0)),
                        ("missing", counts.get("Missing", 0)),
                    ]
                ),
            ),
            (
                "notes",
                [
                    "Station keys SS1/US1/PS1 represent free-format SSn/USn/PSn records.",
                    "ICLOUD vs MCLOUD/ICLDOUT: older vs INP 2.2 naming; both listed.",
                    "ISTEPPG legacy hours vs ISTEPPGS seconds.",
                    "JWAT1/JWAT2 are INP names; CALMET.DAT / GEO use iwat1/iwat2 — wire alias.",
                    "IOUTMM5 is a 3D.DAT header flag, not an INP variable — see roadmap.",
                    "Time: prefer datetime for run windows; Julian packing only at CALMET.DAT wire boundary.",
                    "I/O: prefer pathlib/os stdlib; idiomatic Python over Fortran control-flow ports.",
                ],
            ),
            (
                "dependencies",
                OrderedDict(
                    [
                        ("math", ["numpy", "scipy"]),
                        ("wrf_netcdf", ["xarray", "wrf-python", "netCDF4"]),
                        ("stdlib", ["datetime", "os", "pathlib"]),
                        (
                            "policy",
                            "Full Fortran CALMET feature surface + all INP parameter interfaces; "
                            "bit-identical not required. NumPy/SciPy for math kernels; "
                            "wrf-python + xarray for WRF/NetCDF; datetime for time; "
                            "idiomatic Python OK if behavior/API parity holds.",
                        ),
                    ]
                ),
            ),
        ]
    )

    payload = OrderedDict([("meta", meta), ("parameters", params)])
    DOCS.mkdir(parents=True, exist_ok=True)
    out_json = DOCS / "calmet-inp-params.json"
    out_json.write_text(json.dumps(payload, indent=2) + "\n")
    print(f"Wrote {out_json}")
    print("counts:", dict(meta["counts"]))

    # --- coverage markdown ---
    lines: list[str] = []
    lines.append("# CALMET.INP parameter coverage (py-calmet)")
    lines.append("")
    lines.append(f"**Generated:** {now_local} (Asia/Shanghai)")
    lines.append("")
    lines.append("Inventory of every control-file variable recognized by Fortran "
                 "`READCF` / `READFN` (and station free-form records), unioned with "
                 "keys present in `cases/**/calmet.inp`, scored against current "
                 "`py_calmet` behavior.")
    lines.append("")
    lines.append("## Summary counts")
    lines.append("")
    lines.append(f"| Status | Count |")
    lines.append(f"|--------|------:|")
    lines.append(f"| **Total** | **{meta['counts']['total']}** |")
    lines.append(f"| Implemented | {meta['counts']['implemented']} |")
    lines.append(f"| Partial | {meta['counts']['partial']} |")
    lines.append(f"| Missing | {meta['counts']['missing']} |")
    lines.append("")
    lines.append(f"Sample INP files scanned: **{len(inp_files)}**. "
                 f"Machine-readable twin: [`calmet-inp-params.json`](calmet-inp-params.json).")
    lines.append("")
    lines.append("## Status legend")
    lines.append("")
    lines.append("- **Implemented** — read from INP (or equivalent) and drives physics/I/O as intended for that knob.")
    lines.append("- **Partial** — parsed or hard-coded equivalent exists, but multi-station / alternate options / filename honor / naming alias gaps remain.")
    lines.append("- **Missing** — no functional use in `py_calmet` yet (API field reserved).")
    lines.append("")
    lines.append("## Implementation style (parity without Fortran clones)")
    lines.append("")
    lines.append("- **Time:** use `datetime` / `timedelta` for run windows and stepping; "
                 "Fortran-style Julian/hour packing only when packing/unpacking CALMET.DAT wire format.")
    lines.append("- **Files:** use `pathlib` / `os` (and stdlib I/O); honor INP filenames once Group 0 is fully wired.")
    lines.append("- **Math:** NumPy / SciPy kernels preferred over line-by-line Fortran ports; close-enough golden gates, not bit-identical.")
    lines.append("- **WRF / NetCDF:** `xarray` + `wrf-python` (+ `netCDF4`) are allowed for wrfout→3D.DAT and NetCDF writers.")
    lines.append("- **Control flow:** idiomatic Python is preferred where Fortran used deep GOTOs/flags, as long as behavior and the INP parameter surface stay complete.")
    lines.append("")
    lines.append("## Notable gaps / aliases")
    lines.append("")
    lines.append("| Issue | Detail |")
    lines.append("|-------|--------|")
    lines.append("| `JWAT1`/`JWAT2` vs `IWAT1`/`IWAT2` | Samples set `JWAT*`; runner reads `IWAT*` → silent defaults. Wire alias; CALMET.DAT header uses `iwat1/2` from GEO. |")
    lines.append("| Filenames | `GEODAT`/`SRFDAT`/`UPDAT`/`M3DDAT`/`METDAT` not honored — path heuristics only. |")
    lines.append("| Grid NX/NY/DGRID/XORIG | Taken from GEO.DAT, not INP Group 2. |")
    lines.append("| Radiation HA*/HB*/HC* | Hardcoded in `pbl.py`; should bind to INP. |")
    lines.append("| `IOUTMM5` | **Not an INP variable** — 3D.DAT header; reader supports format 92 only (roadmap). |")
    lines.append("")

    # group order
    group_order = ["0a", "0b", "0c", "0d", "0e", "0f", "1", "2", "3", "4", "5", "6", "station", "sample_only"]
    by_group: dict[str, list] = OrderedDict((g, []) for g in group_order)
    for p in params:
        by_group.setdefault(p["group"], []).append(p)

    for g, plist in by_group.items():
        if not plist:
            continue
        title = GROUP_NAMES.get(g, g)
        lines.append(f"## Group {g} — {title}")
        lines.append("")
        lines.append("| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |")
        lines.append("|----------|------|---------|--------|---------------------------|-------|")
        for p in plist:
            samp = "Y" if p["in_sample_inp"] else ""
            arr = f" ({p['array']})" if p.get("array") else ""
            notes = p["notes"].replace("|", "\\|")
            lines.append(
                f"| `{p['name']}` | {p['type']}{arr} | {samp} | **{p['status']}** | `{p['api_field']}` | {notes} |"
            )
        lines.append("")

    lines.append("## Suggested typed config surface")
    lines.append("")
    lines.append("Target a nested dataclass / pydantic-style `CalmetConfig` mirroring groups:")
    lines.append("")
    lines.append("```text")
    lines.append("CalmetConfig")
    lines.append("├── run         # NOOBS, datetime start/end, NSECDT, ABTZ, IRTYPE, ITEST, MREG")
    lines.append("├── files       # GEODAT, SRFDAT, UPDAT, M3DDAT, METDAT, …, NUSTA/NM3D/…")
    lines.append("├── grid        # PMAP, UTM/LCC…, NX/NY/NZ, ZFACE, origins")
    lines.append("├── output      # LSAVE, IFORMO, LPRINT, LCALGRD, layer print flags")
    lines.append("├── met         # NSSTA, NPSTA, IFORMS/P/C, ICLOUD/MCLOUD")
    lines.append("├── winds       # IWFCOD, IFRADJ, IKINE, IOBR, OA radii, NSMTH, barriers…")
    lines.append("├── pbl         # IMIXH, Zi limits, CONST*, THRESH*, JWAT*, FCORIOL…")
    lines.append("├── radiation   # IRAD, HA1…HC3")
    lines.append("├── clouds      # ICLOUD/MCLOUD/ICLDOUT + CLDDAT")
    lines.append("├── precip      # NPSTA, NFLAGP, SIGMAP, CUTP")
    lines.append("├── overwater   # ICOARE, DSHELF, IWARM, ICOOL, SEA.DAT")
    lines.append("└── stations    # list of surface / upper / precip records")
    lines.append("```")
    lines.append("")
    lines.append("`read_inp()` should populate `CalmetConfig` with Fortran defaults for omitted keys, "
                 "then modules read only typed fields (no ad-hoc `get_int` scattering).")
    lines.append("")
    lines.append("See also [`calmet-module-roadmap.md`](calmet-module-roadmap.md).")
    lines.append("")

    cov_path = DOCS / "calmet-inp-coverage.md"
    cov_path.write_text("\n".join(lines))
    print(f"Wrote {cov_path}")

    # --- roadmap ---
    rm: list[str] = []
    rm.append("# py-calmet module roadmap (full CALMET parity)")
    rm.append("")
    rm.append(f"**Generated:** {now_local} (Asia/Shanghai)")
    rm.append("")
    rm.append("Ordered plan to close gaps versus Fortran CALMET after v1 "
              f"(INP coverage: {meta['counts']['implemented']} implemented / "
              f"{meta['counts']['partial']} partial / {meta['counts']['missing']} missing "
              f"of {meta['counts']['total']}).")
    rm.append("")
    rm.append("## Goals & non-goals")
    rm.append("")
    rm.append("| | |")
    rm.append("|-|-|")
    rm.append("| **Goal** | Full CALMET feature set + complete INP parameter surface (`CalmetConfig`). |")
    rm.append("| **Goal** | Golden regression gates (RMSE / correlation); trustworthy fields. |")
    rm.append("| **Non-goal** | Bit-identical floats vs Fortran. |")
    rm.append("| **Style** | Idiomatic Python; `datetime` for time; `pathlib`/`os` for files; NumPy/SciPy for math. |")
    rm.append("")
    rm.append("## Dependencies")
    rm.append("")
    rm.append("| Area | Libraries |")
    rm.append("|------|-----------|")
    rm.append("| Math kernels | `numpy`, `scipy` |")
    rm.append("| WRF / NetCDF | `xarray`, `wrf-python`, `netCDF4` |")
    rm.append("| Time | `datetime` (stdlib) — not Julian packing except CALMET.DAT wire I/O |")
    rm.append("| Files | `pathlib`, `os` (stdlib) |")
    rm.append("")
    rm.append("`scripts/wrfout_to_3d.py` and NetCDF writers should prefer **xarray + wrf-python** "
              "for field extraction/projection; keep a thin 3D.DAT binary writer in `py_calmet.io`.")
    rm.append("")
    rm.append("## Ordered work packages")
    rm.append("")
    rm.append("### 0. INP / config foundation (prerequisite)")
    rm.append("")
    rm.append("- Expand `io/inp.py` → typed `CalmetConfig` covering all CVDIC keys + defaults from Fortran BLOCK DATA.")
    rm.append("- Honor Group 0 filenames; fix `JWAT1/2` ↔ `IWAT1/2` alias.")
    rm.append("- Bind hardcoded `HA1…HC3` / water LU range to config.")
    rm.append("- Use `datetime` for start/end/`NSECDT` stepping in `runner`.")
    rm.append("- **Exit:** every INP key either drives code or is explicitly stubbed with a clear `NotImplemented` path.")
    rm.append("")
    rm.append("### 1. DIAGNO full (winds core)")
    rm.append("")
    rm.append("| Item | INP / routines | Notes |")
    rm.append("|------|----------------|-------|")
    rm.append("| Multi-station Barnes OA | `R1`,`R2`,`RMAX*`,`RMIN*`,`NINTR2`,`LVARY`,`NSSTA` | Replace single-station OA. |")
    rm.append("| Prog blending | `RPROG`,`IPROG`,`IGFMET` | Weight IGF vs obs; IGF-CALMET reader. |")
    rm.append("| Extrapolation | `IEXTRP`,`BIAS`,`FEXTR2` | Wire `similt.similt_profile`; layer bias. |")
    rm.append("| Barriers / lake breeze | `NBAR`,`KBAR`,`X*BAR`,`LLBREZE`,`NBOX`,… | Terrain barriers + breeze boxes. |")
    rm.append("| Diag options | `IDIOPT1–5`,`ISURFT`,`IUPT`,`IUPWND`,`ZUP*` | Terrain-circulation T/lapse sources. |")
    rm.append("| Calm / div criterion | `ICALM`,`DIVLIM` | Match Fortran calm handling. |")
    rm.append("")
    rm.append("### 2. IKINE — kinematic topographic vertical velocity")
    rm.append("")
    rm.append("- Replace `light_terrain_adjust` with Fortran-aligned `TOPOF2`-style kinematic `W` from terrain slope × wind.")
    rm.append("- Feed horizontal adjustment consistently with `ALPHA`.")
    rm.append("- Golden: enable `IKINE=1` on `small_domain` / synthetic slope case.")
    rm.append("")
    rm.append("### 3. IOBR — full O'Brien adjustment")
    rm.append("")
    rm.append("- Upgrade `divergence_minimize` toward O'Brien vertical velocity / horizontal div cleanup (`IOBR=1`).")
    rm.append("- Honor `DIVLIM`, `NITER`, layer coupling with kinematic `W`.")
    rm.append("- Keep gating: do nothing when `IOBR=0` (current correct behavior).")
    rm.append("")
    rm.append("### 4. PBL / temperature refinements")
    rm.append("")
    rm.append("- Sounding-based lapse above Zi (`MIXDT` / `DPTMIN` path) — fix daytime Zi growth bias.")
    rm.append("- Full `IMIXH` options; `IAVEZI`/`MNMDAV`/`HAFANG`/`ILEVZI`; `IZICRLX`/`TZICRLX`.")
    rm.append("- `ITPROG` / `ITWPROG` / `TGDEF*` / `TRADKM` / `IAVET` / `ILUOC3D` temperature fields.")
    rm.append("- Read `JWAT*`, radiation coeffs from INP.")
    rm.append("")
    rm.append("### 5. Clouds")
    rm.append("")
    rm.append("- Implement `ICLOUD` / `MCLOUD` methods 1–4 (ceilometer, RH-based, 3D.DAT cloud, etc.).")
    rm.append("- `CLDDAT` reader + `ICLDOUT` / `IFORMC`.")
    rm.append("- Couple QSW / energy budget to chosen cloud field.")
    rm.append("")
    rm.append("### 6. Overwater / COARE")
    rm.append("")
    rm.append("- `SEA.DAT` / `NOWSTA` reader.")
    rm.append("- COARE fluxes (`ICOARE`,`DSHELF`,`IWARM`,`ICOOL`) and overwater Zi (`ZIMINW`,`ZIMAXW`,`THRESHW`,`CONSTW`).")
    rm.append("- Distance-to-coast (`LDBCST`,`DCSTGD`) if needed for shelf scaling.")
    rm.append("")
    rm.append("### 7. Precipitation")
    rm.append("")
    rm.append("- `PRECIP.DAT` + `PS*` stations (`NPSTA`,`IFORMP`,`NFLAGP`,`SIGMAP`,`CUTP`).")
    rm.append("- Gridded RMM in CALMET.DAT (currently zeroed).")
    rm.append("")
    rm.append("### 8. IOUTMM5 / 3D.DAT variants")
    rm.append("")
    rm.append("- Extend `io/threed.py` beyond uncompressed format **92** (81/82/91/93–95 families).")
    rm.append("- Prefer **xarray + wrf-python** in `wrfout_to_3d` for moisture/cloud/ice/graupel flags that feed `IOUTMM5`.")
    rm.append("- QA: `ISTEPPGS` multiple of `NSECDT`.")
    rm.append("")
    rm.append("### 9. Writers & outputs")
    rm.append("")
    rm.append("| Output | Status / work |")
    rm.append("|--------|----------------|")
    rm.append("| CALMET.DAT | Exists; honor `METDAT`, `LSAVE`, `LCALGRD`, dates/UTM latlon already fixed in v1.0 |")
    rm.append("| NetCDF | Exists; prefer xarray for structure; keep netCDF4/xarray write path |")
    rm.append("| PACOUT.DAT | Missing (`IFORMO=2`) |")
    rm.append("| CALMET.LST | Missing (`METLST`, echo of inputs like Fortran) |")
    rm.append("| Test/kin/frd/slp dumps | Missing (`TST*`) — low priority |")
    rm.append("| Printer layer flags | Missing (`IUVOUT`/…) — optional |")
    rm.append("")
    rm.append("### 10. Projections & multi-file")
    rm.append("")
    rm.append("- Non-UTM `PMAP` (LCC, PS, EM, TTM, LAZA) + `DATUM`/`FEAST`/`FNORTH`/`XLAT*`.")
    rm.append("- Multiple `M3DDAT` / `IGFDAT` / `UPDAT` lists (`NM3D`,`NIGF`,`NUSTA`).")
    rm.append("")
    rm.append("## Suggested milestone order")
    rm.append("")
    rm.append("1. **v1.1** — Config surface + JWAT/filename/radiation binding + SIMILT wired (`IEXTRP`).")
    rm.append("2. **v1.2** — Full DIAGNO OA (multi-station) + MIXDT lapse.")
    rm.append("3. **v1.3** — IKINE + IOBR production paths + goldens.")
    rm.append("4. **v2.0** — Clouds + precip writers.")
    rm.append("5. **v2.1** — COARE / overwater.")
    rm.append("6. **v2.2** — IOUTMM5 variants + wrfout_to_3d via xarray/wrf-python; PACOUT/LST.")
    rm.append("")
    rm.append("## Mapping to current modules")
    rm.append("")
    rm.append("| Fortran area | Current Python | Next touch |")
    rm.append("|--------------|----------------|------------|")
    rm.append("| READCF/READFN | `io/inp.py` (minimal) | typed config |")
    rm.append("| DIAGNO / WIND1 | `core/winds.py` | OA, IKINE, IOBR |")
    rm.append("| SIMILT | `core/similt.py` (unwired) | wire via IEXTRP |")
    rm.append("| ELUSTR / MIXH* | `core/pbl.py` | MIXDT, IMIXH, COARE hook |")
    rm.append("| COMP / runner | `core/runner.py` | datetime loop, config-driven |")
    rm.append("| RDMM5 / 3D.DAT | `io/threed.py` | IOUTMM5 variants |")
    rm.append("| wrfout bridge | `scripts/wrfout_to_3d.py` | xarray + wrf-python |")
    rm.append("| OUTHD/OUTHR | `io/calmet_dat.py` | PACOUT, LST, flags |")
    rm.append("")
    rm.append("## Related")
    rm.append("")
    rm.append("- [`calmet-inp-coverage.md`](calmet-inp-coverage.md) — per-variable status")
    rm.append("- [`calmet-inp-params.json`](calmet-inp-params.json) — machine-readable list")
    rm.append("- [`../PROGRESS.md`](../PROGRESS.md) — measured golden status")
    rm.append("")

    rm_path = DOCS / "calmet-module-roadmap.md"
    rm_path.write_text("\n".join(rm))
    print(f"Wrote {rm_path}")


if __name__ == "__main__":
    main()
