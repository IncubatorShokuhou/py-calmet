# PROGRESS — py-calmet

**Updated:** 2026-09-18 08:20 CST (Asia/Shanghai)

## Done

- [x] Pure-NumPy package (`py_calmet`) v0.1.0
- [x] Readers: GEO / SURF / UP / 3D.DAT / INP / CALMET.DAT
- [x] **Writers:** CALMET.DAT (round-trip tested) + NetCDF
- [x] Diagnostic core: obs / obs_model / noobs
- [x] **Slope flow** (Allwine–Whiteman / Horst–Doran style) + **divergence minimization**
- [x] **Kinematic W** from horizontal divergence
- [x] PBL: nighttime ELUSTR + MIXHT; **daytime Carson convective ZI** + unstable ELUSTR
- [x] Output fields exercised: U,V,W,T,IPGT,USTAR,ZI,EL,WSTAR,TEMPK,RHO,QSW,IRH,RMM
- [x] Public wrfout demo case (`cases/wrf_demo`) from NCAR Katrina tutorial
- [x] WRF → 3D.DAT converter (`scripts/wrfout_to_3d.py`)
- [x] Three Fortran goldens on wrf_demo (obs / obs_model / noobs)
- [x] pytest green — **21 passed**
- [x] Pushed Python-only tree (no `vendor/` / Fortran / raw wrfout)

## pytest

```
21 passed
```

Repo: https://github.com/IncubatorShokuhou/py-calmet

### Measured relative RMSE — small_domain (tight gates)

| Mode | U | V | ZI | USTAR | SPD |
|------|------|------|------|-------|------|
| obs | ~0.07 | ~0.10 | ~0.01 | ~0.001 | ~0.02 |
| obs_model | ~0.08 | ~0.02 | ~0.08 | ~0.04 | ~0.05 |
| noobs | ~0.02 | ~0.02 | ~0.07 | ~0.02 | ~0.003 |

### Measured — wrf_demo (looser gates + correlation)

| Mode | U rel | V rel | U corr | V corr | ZI rel | USTAR rel |
|------|-------|-------|--------|--------|--------|-----------|
| obs | ~1.2 | ~7.5 | ≥0.3 | ≥0.3 | ~0.34 | ~0.48 |
| obs_model | ~8.5 | ~37 | ≥0.4 | ≥0.5 | ~0.81 | ~1.4 |
| noobs | ~8.5 | ~41 | ~0.73 | ~0.85 | ~0.88 | ~1.4 |

Gates: `tests/thresholds.py` (`THRESH`, `WRF_DEMO_THRESH`).

## wrfout citation

- https://github.com/NCAR/wrf_tutorial_data — `wrfout_d01_2005-08-28_00_00_00` (~63 MB)
- Documented in `cases/wrf_demo/README.md`
- SURF/UP synthesized from wrfout near-surface/column (not real stations)

## Remaining true blockers / later work

1. Full DIAGNO OA / kinematic adjustment still approximate on coarse complex terrain (`wrf_demo` correlations good; relative RMSE large where |ref| is small).
2. Daytime convective path unit-tested; wrf_demo golden window is night/early-morning UTC over Mexico — limited daytime radiative forcing in that archive.
3. IOUTMM5 variants beyond 92, overwater COARE, cloud schemes — out of v1.
4. Raw wrfout not in git (size); rebuild instructions in case README.
