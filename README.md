# py-calmet

Pure-NumPy diagnostic meteorological downscaling inspired by CALMET (CALPUFF system).

## Status

**v0.1.0** — expanded diagnostic core + CALMET.DAT/NetCDF writers + `wrf_demo` case
from a public NCAR Katrina tutorial wrfout.

| Case | Grid | Modes |
|------|------|-------|
| `cases/small_domain` | 12×12 @ 1 km, synthetic | `obs` / `obs_model` / `noobs` |
| `cases/wrf_demo` | 12×12 @ ~28 km, NCAR wrfout subset | `obs` / `obs_model` / `noobs` |

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
- `py_calmet/core/` — diagnostic winds (3D interp, SIMILT, OA, **slope flow**, **divergence minimization**, kinematic W) + PBL (**night + daytime Carson**)
- `scripts/wrfout_to_3d.py` — wrfout NetCDF → 3D.DAT
- `tests/` — reader smoke, synthetic golden parity, wrf_demo parity, physics/writer unit tests
- `cases/small_domain/goldens/` — synthetic Fortran goldens
- `cases/wrf_demo/` — public wrfout-based case (see `cases/wrf_demo/README.md`)

## Parity thresholds

Relative RMSE gates live in `tests/thresholds.py`. Tiny-domain core U/V/ZI/USTAR
remain tight (~1e-2). `wrf_demo` uses looser relative-RMSE gates plus U/V correlation
floors (complex terrain; DIAGNO still approximate).

## License / attribution

CALMET algorithm references: Scire et al. / Exponent CALPUFF system documentation.
This repository does **not** redistribute Fortran CALMET sources.
WRF sample: NCAR `wrf_tutorial_data` (Katrina tutorial).
