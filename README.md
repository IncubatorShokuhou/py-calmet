# py-calmet

Pure-NumPy diagnostic meteorological downscaling inspired by CALMET (CALPUFF system).

## Status — v1 + WP2 + WP3 + WP4 + WP5

**Coverage:** Implemented **207** / Partial **0** / Missing **0** (of 207 INP params).

**v1.0** — diagnostic core + goldens. **WP2** — multi-station Barnes OA + RPROG,
TOPOF2 IKINE, O'Brien IOBR, CLOUD3/4 RH clouds, NPSTA precip, COARE-lite overwater.
**WP3** — MIXDT/MIXDT2 lapse, SEA.DAT + IWARM/ICOOL, PRECIP/CLOUD.DAT I/O,
barriers/lake breeze, non-UTM projections, METLST/PACOUT, IOUTMM5 81–95.
**WP4** — Partial→Implemented surface clearance; MM4/MM5 OutOfScope.
**WP5** — IMIXH=±2 Batchvarova–Gryning (`mixht_day_bg`); IDIOPT*/ZUPT/IUPWND
drive CGAMMA lapse + domain-avg UA wind for Froude/TOPOF2; wrfout via xarray.

| Case | Grid | Modes | Notes |
|------|------|-------|-------|
| `cases/small_domain` | 12×12 @ 1 km, synthetic | `obs` / `obs_model` / `noobs` | Tight relative-RMSE gates |
| `cases/wrf_demo` | 12×12 @ ~28 km, NCAR Katrina wrfout | `obs` / `obs_model` / `noobs` | Floored U/V RMSE + correlation |
| `cases/daytime_zi` | 12×12 @ 1 km, hours 00–17 UTC | `noobs` | Carson/MIXHMC daytime path |

### v1 complete

- Readers + writers: GEO / SURF / UP / 3D.DAT / INP / CALMET.DAT / NetCDF
- Winds: 3D→CALMET interp, IEXTRP profiles (power-law / SIMILT), multi-station Barnes OA (R1/R2/RPROG/RMAX*/NINTR2), **FRADJ**, **slope**, **NSMTH**; **IKINE TOPOF2** + **IOBR O'Brien** (gated)
- PBL: night ELUSTR+MIXHT; daytime energy-budget QH + Maul–Carson **or Batchvarova–Gryning** ZI with MIXDT lapse; RH clouds (MCLOUD 3/4); NPSTA precip + PRECIP.DAT; COARE-lite + SEA.DAT/IWARM/ICOOL
- Golden parity: tiny-domain tight; wrf_demo U/V corr ≳ 0.97 (noobs); daytime ZI corr ≳ 0.99 vs Fortran shape

### OutOfScope / known soft spots

- **MM4DAT / MM5.DAT readers** — OutOfScope by project rule; meteorology input is **wrfout → 3D.DAT** (`scripts/wrfout_to_3d.py`, `py_calmet.io.wrfout` with xarray/wrf-python). Setting a non-default `MM4DAT` raises `NotImplementedError`.
- **IDIOPT1/4/5 = 1** — preprocessed diagnostic-file winds/T are not ingested; QA notes only (defaults `=0` compute internally).
- **WTDAT** — ignored; use SEA.DAT / ITWPROG for overwater temperature.

### Optional extras

```bash
pip install 'py-calmet[wrf]'   # xarray + wrf-python + netCDF4
python scripts/wrfout_to_3d.py wrfout_d01_… out/3d.dat
```
