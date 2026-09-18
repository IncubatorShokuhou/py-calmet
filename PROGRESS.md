# PROGRESS — py-calmet

**Updated:** 2026-09-18 01:14 UTC

## WP2 (this batch) — Partial→Implemented

- Implemented **45 → 72** / Partial **162 → 135** / Missing **0** (of 207)
- DIAGNO: multi-station Barnes OA + RPROG/RMAX*/NINTR2; IEXTRP modes (+4 SIMILT); ISURFT/SS*
- IKINE: TOPOF2 kinematic W; IOBR: O'Brien + DIVLIM (gated; goldens IKINE=IOBR=0)
- Clouds: CLOUD3 Teixeira + CLOUD4-lite; MCLOUD/ICLOUD 3/4 → QSW
- Precip: NPSTA=-1 prognostic rain from 3D.DAT → RMM; station Barnes + SIGMAP/CUTP
- Overwater: COARE-lite bulk fluxes (NOWSTA>0 gate; unit-tested)
- pytest: **52+** green (goldens unchanged)

## Done (v1)

- [x] Pure-NumPy package (`py_calmet`)
- [x] Readers + writers (CALMET.DAT round-trip, NetCDF)
- [x] DIAGNO-aligned winds: FRADJ, Mahrt slope (cdk=0.08), NSMTH, IOBR/IKINE gating
- [x] PBL: night + daytime energy-budget QH + Maul–Carson (MIXHMC) ZI + solar/QSW
- [x] `wrf_demo` NCAR Katrina case; `daytime_zi` convective golden (00–17 UTC)
- [x] pytest green — **23 passed**
- [x] Pushed Python-only tree (no `vendor/` / Fortran / raw wrfout)

## pytest

```
23 passed
```

Repo: https://github.com/IncubatorShokuhou/py-calmet

### Measured — small_domain (unchanged tight gates)

| Mode | U | V | ZI | USTAR | SPD |
|------|------|------|------|-------|------|
| obs | ~0.072 | ~0.102 | ~0.012 | ~0.003 | ~0.024 |
| obs_model | ~0.081 | ~0.016 | ~0.069 | ~0.036 | ~0.042 |
| noobs | ~0.019 | ~0.016 | ~0.063 | ~0.019 | ~0.001 |

### Measured — wrf_demo **before → after** DIAGNO fix

| Mode | metric | before | after |
|------|--------|--------|-------|
| noobs | U corr | ~0.73 | **~0.97** |
| noobs | V corr | ~0.85 | **~0.98** |
| noobs | U abs RMSE | ~1.49 | **~0.50** |
| noobs | V abs RMSE | ~1.75 | **~1.27** |
| noobs | U floored rel (0.5) | — | **~0.44** |
| noobs | V floored rel (0.5) | — | **~1.28** |
| noobs | raw U/V rel RMSE | ~8.5 / ~41 | ~4.2 / ~20 (still inflated by \|ref\|≈0) |
| obs | U/V corr | ~0.38 / ~0.43 | **~0.49 / ~0.63** |
| obs_model | U/V corr | ~0.73 / ~0.85 | **~0.97 / ~0.98** |

Root cause of prior degradation: runner applied divergence minimization despite
`IOBR=0` / `IKINE=0`, with an incorrect slope `cdk`. Fix: gate on INP flags;
Fortran-aligned slope + NSMTH + FRADJ.

### Measured — daytime_zi (noobs, hours 00–17)

| Field | Result |
|-------|--------|
| QSW day rel RMSE | ~0.012 |
| ZI night rel RMSE | ~0.062 |
| ZI day correlation | **~0.995** |
| Day EL / WSTAR | unstable / >0 (Carson path active) |

Absolute daytime ZI growth is faster than Fortran (constant `DPTMIN` gamma vs
sounding lapse) — documented v2 item; shape/correlation gated.

## v1 complete / v2 deferred

See README. In short: core diagnostic + writers + goldens done; COARE, full
clouds, IOUTMM5 variants, MIXDT sounding lapse, multi-station OA → **v2**.

## wrfout citation

- https://github.com/NCAR/wrf_tutorial_data — `wrfout_d01_2005-08-28_00_00_00`
- Documented in `cases/wrf_demo/README.md`
