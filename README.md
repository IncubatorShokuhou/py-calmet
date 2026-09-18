# py-calmet

Pure-NumPy diagnostic meteorological downscaling inspired by CALMET (CALPUFF system).

## Status — v1 complete

**v1.0.0** — diagnostic core with Fortran-aligned DIAGNO steps (Froude, slope,
NSMTH smooth), Maul–Carson daytime ZI, CALMET.DAT/NetCDF writers, and two golden
suites (`small_domain`, `wrf_demo`) plus a daytime convective ZI golden.

| Case | Grid | Modes | Notes |
|------|------|-------|-------|
| `cases/small_domain` | 12×12 @ 1 km, synthetic | `obs` / `obs_model` / `noobs` | Tight relative-RMSE gates |
| `cases/wrf_demo` | 12×12 @ ~28 km, NCAR Katrina wrfout | `obs` / `obs_model` / `noobs` | Floored U/V RMSE + correlation |
| `cases/daytime_zi` | 12×12 @ 1 km, hours 00–17 UTC | `noobs` | Carson/MIXHMC daytime path |

### v1 complete

- Readers + writers: GEO / SURF / UP / 3D.DAT / INP / CALMET.DAT / NetCDF
- Winds: 3D→CALMET geographic interp, power-law obs profiles (`similt_profile` helper exists, not wired), Barnes OA (R1/R2), **FRADJ**, **slope flow** (Mahrt/cdk=0.08), **NSMTH** smooth; **IOBR/IKINE-gated** divergence minimization
- PBL: night ELUSTR+MIXHT; **daytime energy-budget QH + Maul–Carson ZI**; solar/QSW
- Golden parity: tiny-domain tight; wrf_demo U/V corr ≳ 0.97 (noobs); daytime ZI corr ≳ 0.99 vs Fortran shape

### Explicitly deferred to v2

- Full DIAGNO OA multi-station / RPROG blending refinements beyond single-station Barnes
- Sounding-based lapse rates above ZI (MIXDT) — daytime growth uses DPTMIN gamma
- IOUTMM5 variants beyond format 92; full cloud schemes (MCLOUD 2/3/4 detail)
- Overwater **COARE** fluxes / OCD marine mixing heights
- IKINE topographic vertical velocity + full O'Brien (IOBR=1) production tuning
- Raw wrfout in git (rebuild from NCAR tutorial; see `cases/wrf_demo/README.md`)

## Install

```bash
pip install -e ".[dev]"
pytest -q
```

## Quick start

```python
from py_calmet import run_calmet, read_calmet_dat, write_calmet_dat, write_calmet_netcdf

res = run_calmet("cases/small_domain/goldens/obs", mode="obs",
                 inputs_dir="cases/small_domain/goldens/inputs")
print(res.U.shape)  # (ntime, nz, ny, nx)

gold = read_calmet_dat("cases/small_domain/goldens/obs/CALMET.DAT")
```

WRF → 3D.DAT:

```bash
python scripts/wrfout_to_3d.py path/to/wrfout_d01 -o 3d.dat --i0 4 --j0 0 --ni 14 --nj 14
```

## Package layout

- `py_calmet/io/` — GEO / SURF / UP / 3D.DAT / INP / CALMET.DAT **readers + writers**, NetCDF writer
- `py_calmet/core/` — diagnostic winds + PBL (night + daytime Carson)
- `scripts/wrfout_to_3d.py` — wrfout NetCDF → 3D.DAT
- `tests/` — reader smoke, golden parity, wrf_demo, daytime_zi, physics/writer unit tests

## Parity thresholds

Relative RMSE gates live in `tests/thresholds.py`. Tiny-domain core U/V/ZI/USTAR
remain tight (~1e-2–1e-1). `wrf_demo` uses **floored** relative RMSE (0.5 m/s floor
on U/V) plus correlation floors. `daytime_zi` gates QSW, night ZI, and daytime ZI correlation.

## License / attribution

CALMET algorithm references: Scire et al. / Exponent CALPUFF system documentation.
This repository does **not** redistribute Fortran CALMET sources.
WRF sample: NCAR `wrf_tutorial_data` (Katrina tutorial).
