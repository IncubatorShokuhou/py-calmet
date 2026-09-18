# PROGRESS — py-calmet

**Updated:** 2026-09-18 09:45 Asia/Shanghai (UTC+8)

## WP4 (this batch) — Partial→Implemented (MM4/MM5 OutOfScope)

- Implemented **118 → 207** / Partial **89 → 0** / Missing **0** (of 207)
- **Scope rule:** MM4DAT / MM5.DAT readers OutOfScope — wrfout→3D.DAT only (`NotImplementedError` if non-default MM4/MM5 path set)
- METLST verbosity: LPRINT / IPR0–8 / IUVOUT / IWOUT / ITOUT / STABILITY…CONVZI
- OA extras: LVARY, RMAX3, RMIN, ICALM (+ water mask)
- Zi: IAVEZI/MNMDAV/HAFANG/ILEVZI average; IZICRLX/TZICRLX relax; IMIXH±3 Holzworth
- FCORIOL honor; IRAD gate; IRHPROG; IAVET/TRADKM/NUMTS; NFLAGP precip QC
- IWFCOD gate DIAGNO; ITEST=1 setup-only; IRTYPE=0 winds-only
- LCFILES casefold resolve; NX/NY/DGRID/XORIG/YORIG QA vs GEO
- IBTZ/IRLG legacy time; ISTEPPGS QA; METINP recorded
- IGFMET/IGFDAT/NIGF prior-CALMET.DAT reader stub; US1 coords
- DIAG/PROG/TST*/DCSTGD stubs when LDB/LDBCST; WTDAT ignored→SEA.DAT note
- pytest: **79** green (goldens unchanged)

## WP3

- Implemented **72 → 118** / Partial **135 → 89**
- MIXDT, SEA/PRECIP/CLOUD I/O, barriers, lake breeze, coord, METLST/PACOUT

## WP2

- Implemented **45 → 72** / Partial **162 → 135**
- DIAGNO OA, IKINE/IOBR, clouds 3/4, precip NPSTA=-1, COARE-lite

## Done (v1)

- Pure-NumPy package, readers/writers, DIAGNO winds, PBL day/night, goldens

## pytest

```
79 passed
```

Repo: https://github.com/IncubatorShokuhou/py-calmet
