"""Typed CALMET control-file configuration (full INP parameter surface).

All 207 READCF/READFN/station parameters from the inventory are present
with Fortran BLOCK DATA defaults where known. Physics not yet implemented
still accepts the parameter; see docs/calmet-inp-coverage.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Any


INP_ALIASES: dict[str, tuple[str, ...]] = {
    'JWAT1': ('IWAT1',),
    'JWAT2': ('IWAT2',),
}


def _list_field(default: list) -> Any:
    return field(default_factory=lambda d=list(default): list(d))


@dataclass
class CalmetConfig:
    """Complete CALMET.INP parameter object (API parity surface)."""

    # --- Group 0a: File names — primary I/O (READFN subgroup a) ---
    metinp: str = 'calmet.inp'
    geodat: str = 'geo.dat'
    srfdat: str = 'surf.dat'
    prcdat: str = 'precip.dat'
    mm4dat: str = 'mm4.dat'
    wtdat: str = 'wt.dat'
    metlst: str = 'calmet.lst'
    metdat: str = 'calmet.dat'
    pacdat: str = 'pacout.dat'
    clddat: str = 'cloud.dat'
    lcfiles: bool = True
    nusta: int = 0
    nowsta: int = 0
    nm3d: int = 0
    nigf: int = 0

    # --- Group 0b: Upper-air file names (READFN b) ---
    updat: str | list[str] = 'up.dat'

    # --- Group 0c: Overwater / SEA.DAT file names (READFN c) ---
    seadat: str | list[str] = 'sea.dat'

    # --- Group 0d: MM4/MM5/3D.DAT file names (READFN d) ---
    m3ddat: str | list[str] = '3d.dat'

    # --- Group 0e: IGF-CALMET.DAT file names (READFN e) ---
    igfdat: str | list[str] = 'igf.dat'

    # --- Group 0f: Misc diagnostic / test file names (READFN f) ---
    diadat: str = 'diag.dat'
    prgdat: str = 'prog.dat'
    tstprt: str = 'test.dat'
    tstout: str = 'test.out'
    tstkin: str = 'test.kin'
    tstfrd: str = 'test.frd'
    tstslp: str = 'test.slp'
    dcstgd: str = 'dcst.grd'

    # --- Group 1: General run control ---
    ibyr: int = 2000
    ibmo: int = 1
    ibdy: int = 1
    ibhr: int = 0
    ibsec: int = 0
    ieyr: int = 2000
    iemo: int = 1
    iedy: int = 1
    iehr: int = 1
    iesec: int = 0
    abtz: str = 'UTC+0000'
    ibtz: int = 0
    irlg: int = 0
    nsecdt: int = 3600
    irtype: int = 1
    lcalgrd: bool = True
    itest: int = 2
    mreg: int = 0

    # --- Group 2: Map projection and grid ---
    pmap: str = 'UTM'
    datum: str = 'WGS-84'
    feast: float = 0.0
    fnorth: float = 0.0
    iutmzn: int = 0
    utmhem: str = 'N'
    rlat0: str = '0N'
    rlon0: str = '0E'
    xlat1: str = '0N'
    xlat2: str = '0N'
    nx: int = 0
    ny: int = 0
    dgridkm: float = 1.0
    xorigkm: float = 0.0
    yorigkm: float = 0.0
    nz: int = 8
    zface: list[float] = _list_field([0.0, 20.0, 40.0, 80.0, 160.0, 300.0, 600.0, 1000.0, 1500.0])

    # --- Group 3: Output options ---
    lsave: bool = True
    lprint: bool = False
    iprinf: int = 1
    iuvout: list[int] = _list_field([])
    iwout: list[int] = _list_field([])
    itout: list[int] = _list_field([])
    stability: bool = True
    ustar: bool = True
    monin: bool = True
    mixht: bool = True
    wstar: bool = True
    precip: bool = True
    sensheat: bool = True
    convzi: bool = True
    ldb: bool = False
    nn1: int = 1
    nn2: int = 1
    ldbcst: bool = False
    ioutd: int = 0
    nzprn2: int = 1
    ipr0: int = 0
    ipr1: int = 0
    ipr2: int = 0
    ipr3: int = 0
    ipr4: int = 0
    ipr5: int = 0
    ipr6: int = 0
    ipr7: int = 0
    ipr8: int = 0
    iformo: int = 1

    # --- Group 4: Meteorological data options ---
    noobs: int = 0
    nssta: int = 0
    npsta: int = -1
    iforms: int = 2
    iformp: int = 2
    icloud: int = 0
    icldout: int = 0
    mcloud: int = 0
    iformc: int = 2

    # --- Group 5: Wind field options and parameters ---
    iwfcod: int = 1
    ifradj: int = 1
    ikine: int = 0
    iobr: int = 0
    iextrp: int = -4
    rmin2: float = 4.0
    fextr2: list[float] = _list_field([])
    iprog: int = 0
    isteppg: int = 1
    isteppgs: int = 3600
    igfmet: int = 0
    lvary: bool = False
    rmax1: float = 50.0
    rmax2: float = 50.0
    rmax3: float = 50.0
    rmin: float = 0.1
    terrad: float = 10.0
    r1: float = 1.0
    r2: float = 1.0
    rprog: float = 0.0
    divlim: float = 5e-06
    niter: int = 50
    nsmth: list[int] = _list_field([])
    nintr2: list[int] = _list_field([])
    critfn: float = 1.0
    alpha: float = 0.1
    nbar: int = 0
    xbbar: list[float] = _list_field([])
    ybbar: list[float] = _list_field([])
    xebar: list[float] = _list_field([])
    yebar: list[float] = _list_field([])
    kbar: int = 0
    idiopt1: int = 0
    idiopt2: int = 0
    idiopt3: int = 0
    idiopt4: int = 0
    idiopt5: int = 0
    isurft: int = 0
    iupt: int = 0
    zupt: float = 200.0
    iupwnd: int = -1
    zupwnd: list[float] = _list_field([1.0, 1000.0])
    llbreze: bool = False
    nbox: int = 0
    xg1: list[float] = _list_field([])
    xg2: list[float] = _list_field([])
    yg1: list[float] = _list_field([])
    yg2: list[float] = _list_field([])
    xbcst: list[float] = _list_field([])
    ybcst: list[float] = _list_field([])
    xecst: list[float] = _list_field([])
    yecst: list[float] = _list_field([])
    nlb: int = 0
    metbxid: list[int] = _list_field([])
    bias: list[float] = _list_field([])
    islope: int = 1
    icalm: int = 0

    # --- Group 6: Mixing height, temperature, precip, overwater ---
    constb: float = 1.41
    conste: float = 0.15
    constn: float = 2400.0
    dptmin: float = 0.001
    dzzi: float = 200.0
    zimin: float = 50.0
    zimax: float = 3000.0
    ziminw: float = 50.0
    zimaxw: float = 3000.0
    iavezi: int = 1
    mnmdav: int = 1
    hafang: float = 30.0
    ilevzi: int = 1
    fcoriol: float = 999.0
    constw: float = 0.16
    itprog: int = 0
    itwprog: int = 0
    iluoc3d: int = 16
    irad: int = 1
    iavet: int = 1
    tgdefb: float = -0.0098
    tgdefa: float = -0.0045
    jwat1: int = 999
    jwat2: int = 999
    tradkm: float = 500.0
    numts: int = 5
    nflagp: int = 2
    sigmap: float = 100.0
    cutp: float = 0.01
    ha1: float = 990.0
    ha2: float = -30.0
    hb1: float = -0.75
    hb2: float = 3.4
    hc1: float = 5.31e-13
    hc2: float = 60.0
    hc3: float = 0.12
    imixh: int = 1
    threshl: float = 0.05
    threshw: float = 0.05
    icoare: int = 0
    dshelf: float = 0.0
    iwarm: int = 0
    icool: int = 0
    irhprog: int = 0
    izicrlx: int = 1
    tzicrlx: float = 800.0

    # --- Group station: Station location free-format records (IG 7–9) ---
    ss1: str = ''
    us1: str = ''
    ps1: str = ''

    extra: dict[str, Any] = field(default_factory=dict)

    def get(self, name: str, default: Any = None) -> Any:
        """Case-insensitive get by INP name (supports IWAT↔JWAT aliases)."""
        key = name.upper()
        for canon, aliases in INP_ALIASES.items():
            if key == canon or key in aliases:
                key = canon
                break
        attr = key.lower()
        if hasattr(self, attr) and attr != "extra":
            return getattr(self, attr)
        if key in self.extra:
            return self.extra[key]
        return default

    def set_param(self, name: str, value: Any) -> None:
        key = name.upper()
        for canon, aliases in INP_ALIASES.items():
            if key == canon or key in aliases:
                key = canon
                break
        attr = key.lower()
        if hasattr(self, attr) and attr != "extra":
            setattr(self, attr, value)
        else:
            self.extra[key] = value

    def as_inp_dict(self) -> dict[str, Any]:
        out: dict[str, Any] = {}
        for f in fields(self):
            if f.name == "extra":
                continue
            out[f.name.upper()] = getattr(self, f.name)
        out.update(self.extra)
        return out

    @property
    def mode(self) -> str:
        if int(self.noobs) >= 2:
            return "noobs"
        if int(self.iprog) > 0:
            return "obs_model"
        return "obs"

    def _first_int(self, value: Any, default: int) -> int:
        if value is None:
            return default
        if isinstance(value, (list, tuple)):
            if not value:
                return default
            return int(value[0])
        return int(value)

    def effective_iwat(self) -> tuple[int, int]:
        """Water LU range for heat-flux / slope checks.

        INP stores JWAT1/JWAT2 (GEO/header use IWAT1/IWAT2). Fortran default
        JWAT=999 means "no special water body for T"; heatfx still uses GEO
        IWAT defaults of 55. When JWAT is a real LU range, honor it.
        """
        j1 = self._first_int(self.jwat1, 999)
        j2 = self._first_int(self.jwat2, 999)
        if j1 == 999 and j2 == 999:
            return 55, 55
        return j1, j2

    def check_unsupported(self) -> None:
        """Raise NotImplementedError for switches that would silently wrong-result."""
        if int(self.igfmet) != 0:
            raise NotImplementedError('IGFMET!=0 (IGF first-guess) is not implemented')
        if bool(self.llbreze):
            raise NotImplementedError('LLBREZE=T (lake breeze) is not implemented')
        if abs(int(self.imixh)) == 2:
            raise NotImplementedError('IMIXH=±2 (Batchvarova–Gryning) is not implemented')
        if int(self.nbar) > 0:
            raise NotImplementedError('NBAR>0 (wind barriers) is not implemented')


PARAM_NAMES: tuple[str, ...] = tuple(
    f.name.upper() for f in fields(CalmetConfig) if f.name != "extra"
)


def default_config() -> CalmetConfig:
    return CalmetConfig()

