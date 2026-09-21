# py-calmet

Pure-NumPy diagnostic meteorological downscaling inspired by CALMET (CALPUFF system).

## Status — v1 + WP2 + WP3 + WP4 + WP5 + WP6

**Coverage:** Implemented **207** / Partial **0** / Missing **0** (of 207 INP params).

**v1.0** — diagnostic core + goldens. **WP2** — multi-station Barnes OA + RPROG,
TOPOF2 IKINE, O'Brien IOBR, CLOUD3/4 RH clouds, NPSTA precip, COARE-lite overwater.
**WP3** — MIXDT/MIXDT2 lapse, SEA.DAT + IWARM/ICOOL, PRECIP/CLOUD.DAT I/O,
barriers/lake breeze, non-UTM projections, METLST/PACOUT, IOUTMM5 81–95.
**WP4** — Partial→Implemented surface clearance; MM4/MM5 OutOfScope.
**WP5** — IMIXH=±2 Batchvarova–Gryning (`mixht_day_bg`); IDIOPT*/ZUPT/IUPWND
drive CGAMMA lapse + domain-avg UA wind for Froude/TOPOF2; wrfout via xarray.
**WP6** — IDIOPT1/4/5=1 DIAG.DAT ingestion; WTDAT water-T with
ITWPROG>SEA>WT precedence; IRTYPE=0 PBL collapse fix; shared `G=9.81`.

| Case | Grid | Modes | Notes |
|------|------|-------|-------|
| `cases/wrf_demo` | 12×12 @ ~28 km, NCAR Katrina wrfout + GEO.DAT | `obs` / `obs_model` / `noobs` | **Sole** Fortran-compare gate (floored U/V RMSE + correlation) |

Synthetic `small_domain` / `daytime_zi` cases were removed — validation is WRF-only.

### v1 complete

- Readers + writers: GEO / SURF / UP / 3D.DAT / INP / CALMET.DAT / NetCDF
- Winds: 3D→CALMET interp, IEXTRP profiles (power-law / SIMILT), multi-station Barnes OA (R1/R2/RPROG/RMAX*/NINTR2), **FRADJ**, **slope**, **NSMTH**; **IKINE TOPOF2** + **IOBR O'Brien** (gated)
- PBL: night ELUSTR+MIXHT; daytime energy-budget QH + Maul–Carson **or Batchvarova–Gryning** ZI with MIXDT lapse; RH clouds (MCLOUD 3/4); NPSTA precip + PRECIP.DAT; COARE-lite + SEA.DAT/IWARM/ICOOL
- Golden parity: **wrf_demo only** (real Katrina WRF + mountain GEO.DAT); noobs U/V corr ≳ 0.97 vs Fortran

### OutOfScope / known soft spots

- **MM4DAT / MM5.DAT readers** — OutOfScope by project rule; meteorology input is **wrfout → 3D.DAT** (`scripts/wrfout_to_3d.py`, `py_calmet.io.wrfout` with xarray/wrf-python). Setting a non-default `MM4DAT` raises `NotImplementedError`.
- **Official WT.DAT terrain-weight layout** (CALMET User's Guide §8.10) — not implemented; this package's `WTDAT` soft-spot path is **overwater SST** (see `io.wt_dat`).

## Validation (WRF-only)

Primary gate: `tests/test_wrf_demo.py` and `scripts/compare_wrf_fortran.py` against
Fortran CALMET outputs for the Katrina wrfout case (with `geo.dat` terrain).
Live `calmet.x` is used when present; otherwise archived Fortran-from-WRF-case
goldens under `cases/wrf_demo/goldens/`.

```bash
pytest tests/test_wrf_demo.py -v
python scripts/compare_wrf_fortran.py
```

### Optional extras


```bash
pip install 'py-calmet[wrf]'   # xarray + wrf-python + netCDF4
python scripts/wrfout_to_3d.py wrfout_d01_… out/3d.dat
```
