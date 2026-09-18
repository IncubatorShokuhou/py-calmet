"""Minimal CALMET.INP parser for keys used by the tiny-domain configs."""
from __future__ import annotations
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List
import re


@dataclass
class CalmetInp:
    raw: Dict[str, str] = field(default_factory=dict)
    mode: str = "obs"  # inferred

    def get(self, key: str, default=None):
        return self.raw.get(key, default)

    def get_int(self, key: str, default: int = 0) -> int:
        v = self.raw.get(key)
        if v is None:
            return default
        return int(float(str(v).split(',')[0].strip()))

    def get_float(self, key: str, default: float = 0.0) -> float:
        v = self.raw.get(key)
        if v is None:
            return default
        return float(str(v).split(',')[0].strip())

    def get_bool(self, key: str, default: bool = False) -> bool:
        v = self.raw.get(key)
        if v is None:
            return default
        return str(v).strip().upper().startswith('T')

    def get_list_float(self, key: str) -> List[float]:
        v = self.raw.get(key, '')
        return [float(x.strip()) for x in str(v).split(',') if x.strip()]

    def get_list_int(self, key: str) -> List[int]:
        v = self.raw.get(key, '')
        if not v:
            return []
        return [int(float(x.strip())) for x in str(v).split(',') if x.strip()]


def read_inp(path: str | Path) -> CalmetInp:
    text = Path(path).read_text()
    raw: Dict[str, str] = {}
    for m in re.finditer(r'!\s*([A-Z0-9]+)\s*=\s*([^!]*)!', text):
        raw[m.group(1)] = m.group(2).strip()
    for m in re.finditer(r'!\s*(SS\d+|US\d+)\s*=\s*([^!]*)!', text):
        raw[m.group(1)] = m.group(2).strip()
    inp = CalmetInp(raw=raw)
    noobs = inp.get_int('NOOBS', 0)
    iprog = inp.get_int('IPROG', 0)
    if noobs >= 2:
        inp.mode = 'noobs'
    elif iprog > 0:
        inp.mode = 'obs_model'
    else:
        inp.mode = 'obs'
    return inp
