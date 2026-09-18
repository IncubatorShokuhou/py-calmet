# PROGRESS — py-calmet

**Updated:** 2026-09-18 09:30 Asia/Shanghai (UTC+8)

## WP3 (this batch) — Partial→Implemented

- Implemented **72 → 118** / Partial **135 → 89** / Missing **0** (of 207)
- MIXDT / MIXDT2 sounding + 3D.DAT lapse above Zi (`ITPROG`, `DPTMIN`, `DZZI`, `CONSTE`)
- SEA.DAT reader (v2.0/2.1/2.11) + COARE-lite IWARM/ICOOL/THRESHW + wave z0
- PRECIP.DAT reader + NPSTA>0 Barnes rates; CLOUD.DAT reader/writer (`ICLDOUT`)
- Barriers (`NBAR`/`KBAR`/XY*) in OA; lake breeze (`LLBREZE`/`NBOX`/…)
- INP-driven map projections (`PMAP` UTM/LCC/TM/… + FEAST/FNORTH/XLAT*)
- METLST writer; PACOUT npz hook (`IFORMO=2`); more IOUTMM5 parse (81–95)
- pytest: **66** green (goldens unchanged)

## WP2

- Implemented **45 → 72** / Partial **162 → 135**
- DIAGNO OA, IKINE/IOBR, clouds 3/4, precip NPSTA=-1, COARE-lite

## Done (v1)

- Pure-NumPy package, readers/writers, DIAGNO winds, PBL day/night, goldens

## pytest

```
66 passed
```

Repo: https://github.com/IncubatorShokuhou/py-calmet
