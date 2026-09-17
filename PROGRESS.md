# PROGRESS — py-calmet

**Updated:** 2026-09-17 23:35 CST

## Done

- [x] Pure-NumPy package (`py_calmet`)
- [x] Readers: GEO / SURF / UP / 3D.DAT / INP / CALMET.DAT
- [x] Diagnostic core for obs / obs_model / noobs
- [x] PBL: nighttime ELUSTR + MIXHT (Venkatram / Zilitinkevich)
- [x] pytest green — **12 passed**
- [x] Pushed Python-only tree to `main` (no `vendor/` / Fortran)

## pytest

```
12 passed in 0.20s
```

Repo: https://github.com/IncubatorShokuhou/py-calmet  
Commit: `60f334a` (feat: pure-NumPy diagnostic core with golden parity tests)

### Measured relative RMSE (vs goldens)

| Mode | U | V | ZI | USTAR | SPD |
|------|------|------|------|-------|------|
| obs | 0.072 | 0.102 | 0.008 | 0.001 | 0.024 |
| obs_model | 0.085 | 0.019 | 0.078 | 0.038 | 0.045 |
| noobs | 0.020 | 0.016 | 0.066 | 0.019 | 0.003 |

Gates live in `tests/thresholds.py`.

## Remaining gaps

1. Full DIAGNO OA / slope / mass-consistency — simplified for tiny domain
2. Daytime convective PBL — night path only in this golden window
3. W / T-LEV / IPGT / WSTAR / RHO / QSW / IRH — best-effort, not gated
4. CALMET.DAT writer — not yet
5. Broader domains / IOUTMM5 variants — unvalidated
