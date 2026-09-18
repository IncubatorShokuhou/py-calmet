"""CALMET.INP reader/writer with full parameter surface and round-trip support."""
from __future__ import annotations

import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from ..config import INP_ALIASES, PARAM_NAMES, CalmetConfig, default_config

_ASSIGN_RE = re.compile(r"!\s*([A-Z][A-Z0-9]*)\s*=\s*([^!]*)!", re.IGNORECASE)
_STATION_RE = re.compile(r"^(SS|US|PS)\d+$")


def _parse_logical(s: str) -> bool:
    t = s.strip().upper().lstrip(".")
    if t.startswith("T"):
        return True
    if t.startswith("F"):
        return False
    raise ValueError(f"not a CALMET logical: {s!r}")


def _alias_canon(name: str) -> str:
    key = name.upper()
    for canon, aliases in INP_ALIASES.items():
        if key == canon or key in aliases:
            return canon
    return key


def _parse_value(raw: str, name: str, current: Any) -> Any:
    s = raw.strip()
    if s == "":
        return current

    if _STATION_RE.match(name) or name in ("SS1", "US1", "PS1"):
        return s

    if (s.startswith("'") and s.endswith("'")) or (s.startswith('"') and s.endswith('"')):
        return s[1:-1].strip()

    if isinstance(current, bool) or s.upper().lstrip(".") in (
        "T",
        "F",
        "TRUE",
        "FALSE",
    ):
        try:
            return _parse_logical(s.split(",")[0])
        except ValueError:
            pass

    if "," in s:
        parts = [p.strip() for p in s.split(",") if p.strip() != ""]
        if isinstance(current, list):
            if current and isinstance(current[0], float):
                return [float(p) for p in parts]
            if current and isinstance(current[0], int):
                return [int(float(p)) for p in parts]
            try:
                if any(("." in p or "e" in p.lower()) for p in parts):
                    return [float(p) for p in parts]
                return [int(float(p)) for p in parts]
            except ValueError:
                return parts
        # scalar field given a list in INP — use first token
        s = parts[0]

    if isinstance(current, bool):
        return _parse_logical(s)
    if isinstance(current, int) and not isinstance(current, bool):
        return int(float(s))
    if isinstance(current, float):
        return float(s)
    if isinstance(current, list):
        try:
            if "." in s or "e" in s.lower():
                return [float(s)]
            return [int(float(s))]
        except ValueError:
            return [s]
    return s


def _format_value(value: Any) -> str:
    if isinstance(value, bool):
        return "T" if value else "F"
    if isinstance(value, list):
        return ", ".join(_format_value(v) for v in value)
    if isinstance(value, float):
        if value == 0.0:
            return "0.0"
        abs_v = abs(value)
        if abs_v >= 1e5 or (0.0 < abs_v < 1e-4):
            return f"{value:.6g}"
        return f"{value:.8g}"
    if isinstance(value, int):
        return str(value)
    text = str(value)
    if text == "" or any(c.isspace() for c in text) or text.upper() in ("T", "F"):
        return f"'{text}'"
    return text


@dataclass
class CalmetInp:
    """Parsed CALMET.INP: typed config + raw map + optional source text."""

    config: CalmetConfig = field(default_factory=default_config)
    raw: dict[str, str] = field(default_factory=dict)
    unknown: dict[str, str] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)
    source_text: str | None = None
    source_path: Path | None = None

    @property
    def mode(self) -> str:
        return self.config.mode


    def _raw_has(self, key: str) -> bool:
        k = key.upper()
        if k in self.raw:
            return True
        canon = _alias_canon(k)
        if canon != k and canon in self.raw:
            return True
        for canon_name, aliases in INP_ALIASES.items():
            if k == canon_name or k in aliases:
                if canon_name in self.raw or any(a in self.raw for a in aliases):
                    return True
        return False

    def _raw_text(self, key: str) -> str | None:
        k = key.upper()
        if k in self.raw:
            return self.raw[k]
        canon = _alias_canon(k)
        if canon in self.raw:
            return self.raw[canon]
        for canon_name, aliases in INP_ALIASES.items():
            if k == canon_name or k in aliases:
                if canon_name in self.raw:
                    return self.raw[canon_name]
                for a in aliases:
                    if a in self.raw:
                        return self.raw[a]
        return None

    def get(self, key: str, default: Any = None) -> Any:
        if self._raw_has(key):
            # Prefer typed config when populated; else raw text
            typed = self.config.get(key, None)
            raw = self._raw_text(key)
            # If config still holds the class default and raw differs, parse raw
            return typed if typed is not None else (raw if raw is not None else default)
        if default is not None:
            return default
        return self.config.get(key, default)

    def get_int(self, key: str, default: int = 0) -> int:
        raw = self._raw_text(key)
        if raw is not None:
            try:
                return int(float(str(raw).split(",")[0].strip()))
            except ValueError:
                pass
        if not self._raw_has(key):
            return default
        v = self.config.get(key, None)
        if v is None:
            return default
        if isinstance(v, (list, tuple)):
            return int(v[0]) if v else default
        return int(v)

    def get_float(self, key: str, default: float = 0.0) -> float:
        raw = self._raw_text(key)
        if raw is not None:
            try:
                return float(str(raw).split(",")[0].strip())
            except ValueError:
                return default
        if not self._raw_has(key):
            return default
        v = self.config.get(key, None)
        if v is None:
            return default
        if isinstance(v, (list, tuple)):
            return float(v[0]) if v else default
        try:
            return float(v)
        except (TypeError, ValueError):
            return default

    def get_bool(self, key: str, default: bool = False) -> bool:
        raw = self._raw_text(key)
        if raw is not None:
            try:
                return _parse_logical(str(raw).split(",")[0])
            except ValueError:
                return default
        if not self._raw_has(key):
            return default
        v = self.config.get(key, None)
        if v is None:
            return default
        return bool(v)

    def get_list_float(self, key: str) -> list[float]:
        raw = self._raw_text(key)
        if raw is not None:
            return [float(x.strip()) for x in str(raw).split(",") if x.strip()]
        if not self._raw_has(key):
            return []
        v = self.config.get(key, None)
        if v is None:
            return []
        if isinstance(v, list):
            return [float(x) for x in v]
        return [float(v)]

    def get_list_int(self, key: str) -> list[int]:
        raw = self._raw_text(key)
        if raw is not None:
            return [int(float(x.strip())) for x in str(raw).split(",") if x.strip()]
        if not self._raw_has(key):
            return []
        v = self.config.get(key, None)
        if v is None:
            return []
        if isinstance(v, list):
            return [int(x) for x in v]
        return [int(v)]



def read_inp(path: str | Path) -> CalmetInp:
    path = Path(path)
    text = path.read_text(errors="replace")
    return parse_inp_text(text, source_path=path)


def parse_inp_text(text: str, source_path: Path | None = None) -> CalmetInp:
    cfg = default_config()
    raw: dict[str, str] = {}
    unknown: dict[str, str] = {}
    warns: list[str] = []
    known = set(PARAM_NAMES)

    for m in _ASSIGN_RE.finditer(text):
        name = m.group(1).upper()
        val_text = m.group(2).strip()
        raw[name] = val_text
        canon = _alias_canon(name)

        if _STATION_RE.match(canon):
            prefix = canon[:2]
            if prefix == "SS":
                if canon == "SS1" or not cfg.ss1:
                    cfg.ss1 = val_text
                else:
                    cfg.extra[canon] = val_text
            elif prefix == "US":
                if canon == "US1" or not cfg.us1:
                    cfg.us1 = val_text
                else:
                    cfg.extra[canon] = val_text
            else:
                if canon == "PS1" or not cfg.ps1:
                    cfg.ps1 = val_text
                else:
                    cfg.extra[canon] = val_text
            continue

        if canon in known:
            current = getattr(cfg, canon.lower())
            try:
                coerced = _parse_value(val_text, canon, current)
            except Exception as exc:  # noqa: BLE001
                warns.append(f"{canon}: could not parse {val_text!r} ({exc}); keeping default")
                continue
            setattr(cfg, canon.lower(), coerced)
        else:
            unknown[name] = val_text
            cfg.extra[name] = val_text
            warns.append(f"unknown INP key preserved: {name}")

    for w in warns:
        if w.startswith("unknown"):
            warnings.warn(w, stacklevel=2)

    return CalmetInp(
        config=cfg,
        raw=raw,
        unknown=unknown,
        warnings=warns,
        source_text=text,
        source_path=source_path,
    )


def write_inp(
    path: str | Path,
    inp: CalmetInp | CalmetConfig,
    *,
    template: str | Path | None = None,
) -> None:
    """Write INP. Prefer updating a template/source in place for round-trips."""
    path = Path(path)
    if isinstance(inp, CalmetConfig):
        cfg = inp
        src: str | None = None
        raw: dict[str, str] = {}
    else:
        cfg = inp.config
        src = inp.source_text
        raw = dict(inp.raw)

    if template is not None:
        tpath = Path(template) if not isinstance(template, Path) else template
        if isinstance(template, str) and ("\n" in template or "!" in template) and not tpath.exists():
            src = template
        else:
            src = Path(template).read_text(errors="replace")

    values = cfg.as_inp_dict()

    if src:
        def repl(m: re.Match[str]) -> str:
            name = m.group(1).upper()
            canon = _alias_canon(name)
            if canon not in values and name not in values:
                return m.group(0)
            val = values.get(canon, values.get(name))
            return f"! {name} = {_format_value(val)} !"

        path.write_text(_ASSIGN_RE.sub(repl, src))
        return

    lines = ["CALMET.INP      2.1             Generated by py-calmet", ""]
    for name in PARAM_NAMES:
        val = values[name]
        if val == "" or val == []:
            continue
        lines.append(f"! {name} = {_format_value(val)} !")
    for name, val in cfg.extra.items():
        lines.append(f"! {name} = {_format_value(val)} !")
    lines.append("")
    path.write_text("\n".join(lines))
