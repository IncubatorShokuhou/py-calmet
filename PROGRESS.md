# py-calmet progress

## Align-light (马尾) — missing-obs sentinels

- **UP.DAT** missing `ws`/`wd` (≥998) filtered in `obs_profile_similt` (same gate as VERTAV/domain-avg); empty UA falls back to surface power-law — no more phantom ~100 m/s layers.
- **SURF.DAT** `9999` → calm in `obs_surface_uv` (not a 9999 m/s wind).
- Coverage note for **WTDAT** corrected (SST soft-spot, not “ignored”).
- Skipped: COARE/DIAGNO rewrites, bit-identical floats, SURF temp/RH missing (winds-only), GEO N→S row order, slope `rho` (Fortran RHOCP=1229.9 constant).


## WP6 (this batch) — IDIOPT soft spots + WTDAT + Fortran align

- **IDIOPT1/4/5=1**: real DIAG.DAT ingestion (`io.diag_dat` ASCII subset); wired in
  runner for surface T, surface UV (IRTYPE=0), upper UV (IRTYPE=0). IDIOPT2/3=1
  also consume DIAG GAMMA / domain UV when present.
- **WTDAT**: `io.wt_dat` reader; water-T precedence **ITWPROG > SEA.DAT > WT.DAT >
  air-T** via `overwater.resolve_water_temp`; stop documenting as ignored.
- **Part B (verified)**: IRTYPE=0 PBL collapse moved after `wstar_field` (was
  NameError / overwritten); Froude/TOPOF2/slope use shared `met_utils.G` (9.81).
- **MM4/MM5** remain OutOfScope.
- Assumption: DIAG.DAT uses documented ASCII fields (full Douglas–Kessler binary
  not vendored); WTDAT here is SST soft-spot, not official terrain-weight WT.DAT.
- Tests: `tests/test_wp6_idiopt_wtdat.py`; CI will validate (no local pytest).

## WP5 — Batchvarova–Gryning + diag depth

- **IMIXH=±2** implemented: `pbl.mixht_day_bg` / `_mixhbg_scalar` (Fortran MIXHBG+FBG analytical false-images + secant); wired in runner; `check_unsupported` no longer raises for ±2.
- **IDIOPT2/ZUPT**: CGAMMA-style lapse from UP sounding feeds Froude + TOPOF2 gamma (replacing hard-coded 0.01/0.005 when sounding present).
- **IDIOPT3/IUPWND/ZUPWND**: domain-avg UA wind via VERTAV-style layer average → meta/METLST.
- **IDIOPT1/4/5=1**: (superseded by WP6 DIAG.DAT ingestion).
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
