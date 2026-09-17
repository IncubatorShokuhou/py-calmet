# PROGRESS — py-calmet

**Updated:** 2026-09-17 23:30 CST

## Done

- [x] Pure-NumPy package at `/workspace/py-calmet`
- [x] Readers: GEO.DAT, SURF.DAT, UP.DAT, 3D.DAT, CALMET.INP, CALMET.DAT (Fortran sequential)
- [x] Diagnostic core for three modes (obs / obs_model / noobs)
- [x] PBL: nighttime ELUSTR (ustar, EL, QH) + MIXHT Venkatram/Zilitinkevich ZI
- [x] pytest green: **12 passed** (readers + core parity × 3 modes)
- [x] Fortran vendor tree excluded from git (`.gitignore`)

## pytest summary

```
12 passed in 0.23s
```

Parity metric: relative RMSE  
`sqrt(mean(((pred-ref)/(|ref|+eps))^2))` vs `cases/small_domain/goldens/{mode}/CALMET.DAT`.

| Mode | U | V | U L1 | V L1 | ZI | USTAR | SPD |
|------|------|------|------|------|------|-------|------|
| obs | ≤0.10 | ≤0.15 | ≤0.01 | ≤0.01 | ≤0.02 | ≤0.01 | ≤0.05 |
| noobs | ≤0.05 | ≤0.05 | ≤0.05 | ≤0.05 | ≤0.10 | ≤0.05 | ≤0.05 |
| obs_model | ≤0.20 | ≤0.10 | ≤0.25 | ≤0.15 | ≤0.12 | ≤0.08 | ≤0.15 |

## Remaining gaps

1. **Full DIAGNO / OA / slope / divergence minimization** — not fully ported; obs aloft winds use a power-law + UA direction blend (IEXTRP=-4 approximation). Layer-1 obs winds and PBL match tightly.
2. **Daytime heat flux / convective ZI** — night-only path exercised by this golden window (QSW=0).
3. **W / T-LEV / IPGT / WSTAR / RHO / QSW / IRH** — produced best-effort; not gated in pytest.
4. **CALMET.DAT writer** — reader only in v0.0.1.
5. **Arbitrary domains / IOUTMM5 variants** — validated on the synthetic 12×12 case only.

## Commands

```bash
cd /workspace/py-calmet
python -m venv .venv && .venv/bin/pip install -e ".[dev]"
.venv/bin/pytest -q
```
