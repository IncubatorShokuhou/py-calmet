"""DIAG.DAT reader — preprocessed diagnostic inputs for IDIOPT1–5=1.

Fortran CALMET reads Douglas–Kessler-compatible DIAG.DAT when any of
IDIOPT1/2/3/4/5 is 1 (surface T, lapse, domain UV, surface UV, upper UV).
The full binary/legacy layout is not vendored in this tree; py-calmet uses a
documented ASCII subset that carries the same physical fields.

Accepted formats
----------------
1. Positional hourly records::

       YYYY JJJ HH TSFC GAMMA UDOM VDOM USFC VSFC UUP VUP

   Missing sentinel ``9999`` (or omit trailing columns). Temperatures in K;
   winds in m/s (U positive east, V positive north); GAMMA in K/m.

2. Keyword records (any subset)::

       YYYY JJJ HH TSFC=288.15 USFC=3.0 VSFC=0.5 UUP=5.0 VUP=1.0

Header lines starting with ``*``, ``DIAG``, or blank are skipped.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List

MISS = 9999.0


def _finite(x: float) -> bool:
    return x == x and x not in (float("inf"), float("-inf"))


@dataclass
class DiagRecord:
    year: int
    jday: int
    hour: int
    tsfc: float = MISS  # K
    gamma: float = MISS  # K/m
    udom: float = MISS  # m/s domain-avg
    vdom: float = MISS
    usfc: float = MISS  # m/s surface
    vsfc: float = MISS
    uup: float = MISS  # m/s upper
    vup: float = MISS

    @property
    def code(self) -> int:
        return self.year * 100000 + self.jday * 100 + self.hour

    def has_tsfc(self) -> bool:
        return self.tsfc < MISS - 1.0 and _finite(self.tsfc)

    def has_gamma(self) -> bool:
        return self.gamma < MISS - 1.0 and _finite(self.gamma)

    def has_domain_uv(self) -> bool:
        return (
            self.udom < MISS - 1.0
            and self.vdom < MISS - 1.0
            and _finite(self.udom)
            and _finite(self.vdom)
        )

    def has_sfc_uv(self) -> bool:
        return (
            self.usfc < MISS - 1.0
            and self.vsfc < MISS - 1.0
            and _finite(self.usfc)
            and _finite(self.vsfc)
        )

    def has_up_uv(self) -> bool:
        return (
            self.uup < MISS - 1.0
            and self.vup < MISS - 1.0
            and _finite(self.uup)
            and _finite(self.vup)
        )


@dataclass
class DiagData:
    records: List[DiagRecord] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


_KEYS = ("TSFC", "GAMMA", "UDOM", "VDOM", "USFC", "VSFC", "UUP", "VUP")


def _parse_kv(parts: list[str]) -> dict[str, float]:
    out: dict[str, float] = {}
    for p in parts:
        if "=" not in p:
            continue
        k, _, v = p.partition("=")
        ku = k.strip().upper()
        if ku in _KEYS:
            try:
                out[ku] = float(v)
            except ValueError:
                continue
    return out


def read_diag(path: str | Path) -> DiagData:
    """Read ASCII DIAG.DAT (IDIOPT preprocessed fields)."""
    lines = Path(path).read_text(errors="replace").splitlines()
    records: List[DiagRecord] = []
    notes: List[str] = []
    for line in lines:
        s = line.strip()
        if not s or s.startswith("*") or s.upper().startswith("DIAG"):
            continue
        p = s.replace(",", " ").split()
        if len(p) < 3:
            continue
        try:
            year, jday, hour = int(float(p[0])), int(float(p[1])), int(float(p[2]))
        except ValueError:
            continue
        if hour == 24:
            hour = 0
            jday += 1
        rec = DiagRecord(year=year, jday=jday, hour=hour)
        kv = _parse_kv(p[3:])
        if kv:
            for attr, key in (
                ("tsfc", "TSFC"),
                ("gamma", "GAMMA"),
                ("udom", "UDOM"),
                ("vdom", "VDOM"),
                ("usfc", "USFC"),
                ("vsfc", "VSFC"),
                ("uup", "UUP"),
                ("vup", "VUP"),
            ):
                if key in kv:
                    setattr(rec, attr, float(kv[key]))
        else:
            vals: list[float] = []
            for tok in p[3:]:
                try:
                    vals.append(float(tok))
                except ValueError:
                    break
            names = ("tsfc", "gamma", "udom", "vdom", "usfc", "vsfc", "uup", "vup")
            for i, name in enumerate(names):
                if i < len(vals):
                    setattr(rec, name, vals[i])
        records.append(rec)
    if not records:
        notes.append(f"DIAG.DAT {path}: no hourly records parsed")
    return DiagData(records=records, notes=notes)


def diag_for_hour(
    data: DiagData,
    year: int,
    jday: int,
    hour: int,
) -> DiagRecord | None:
    """Exact hour match, else nearest earlier, else first record."""
    if not data.records:
        return None
    code = year * 100000 + jday * 100 + hour
    best = None
    for r in data.records:
        if r.code == code:
            return r
        if r.code <= code and (best is None or r.code > best.code):
            best = r
    return best if best is not None else data.records[0]
