# py-calmet progress

## WP5 (this batch) — Batchvarova–Gryning + diag depth

- **IMIXH=±2** implemented: `pbl.mixht_day_bg` / `_mixhbg_scalar` (Fortran MIXHBG+FBG analytical false-images + secant); wired in runner; `check_unsupported` no longer raises for ±2.
- **IDIOPT2/ZUPT**: CGAMMA-style lapse from UP sounding feeds Froude + TOPOF2 gamma (replacing hard-coded 0.01/0.005 when sounding present).
- **IDIOPT3/IUPWND/ZUPWND**: domain-avg UA wind via VERTAV-style layer average → meta/METLST.
- **IDIOPT1/4/5=1**: QA notes (preprocessed diag files unsupported); IDIOPT4/5+IRTYPE≠0 flagged.
- **wrfout**: `py_calmet.io.wrfout` (xarray first-class); MM4/MM5 stay OutOfScope.
- Tests: `tests/test_mixhbg.py`; full pytest green.

## WP4 — Partial→Implemented (MM4/MM5 OutOfScope)

- Implemented **118 → 207** / Partial **89 → 0** / Missing **0** (of 207)
- **Scope rule:** MM4DAT / MM5.DAT readers OutOfScope — wrfout→3D.DAT only (`NotImplementedError` if non-default MM4/MM5 path set)
- Zi: IAVEZI/MNMDAV/HAFANG/ILEVZI average; IZICRLX/TZICRLX relax; IMIXH±3 Holzworth

## WP3 — MIXDT, SEA/PRECIP/CLOUD I/O, barriers, coord, outputs

- Implemented **72 → 118** / Partial **135 → 89**

## WP2 — TOPOF2/OBrien, RH clouds, precip, COARE-lite

- Implemented **45 → 72** / Partial **162 → 135**
