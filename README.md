# py-calmet

Pure-NumPy diagnostic meteorological downscaling inspired by CALMET (CALPUFF system).

## Status — v1 + WP2 physics

**Coverage:** Implemented **72** / Partial **135** / Missing **0** (of 207 INP params).

**v1.0** — diagnostic core + goldens. **WP2** — multi-station Barnes OA + RPROG,
TOPOF2 IKINE, O'Brien IOBR, CLOUD3/4 RH clouds, NPSTA precip, COARE-lite overwater.

| Case | Grid | Modes | Notes |
|------|------|-------|-------|
| `cases/small_domain` | 12×12 @ 1 km, synthetic | `obs` / `obs_model` / `noobs` | Tight relative-RMSE gates |
| `cases/wrf_demo` | 12×12 @ ~28 km, NCAR Katrina wrfout | `obs` / `obs_model` / `noobs` | Floored U/V RMSE + correlation |
| `cases/daytime_zi` | 12×12 @ 1 km, hours 00–17 UTC | `noobs` | Carson/MIXHMC daytime path |

### v1 complete

- Readers + writers: GEO / SURF / UP / 3D.DAT / INP / CALMET.DAT / NetCDF
- Winds: 3D→CALMET interp, IEXTRP profiles (power-law / SIMILT), multi-station Barnes OA (R1/R2/RPROG/RMAX*/NINTR2), **FRADJ**, **slope**, **NSMTH**; **IKINE TOPOF2** + **IOBR O'Brien** (gated)
- PBL: night ELUSTR+MIXHT; daytime energy-budget QH + Maul–Carson ZI; RH clouds (MCLOUD 3/4); NPSTA precip; COARE-lite overwater
- Golden parity: tiny-domain tight; wrf_demo U/V corr ≳ 0.97 (noobs); daytime ZI corr ≳ 0.99 vs Fortran shape

### Still Partial / next WPs

- Sounding-based lapse rates above ZI (MIXDT) — daytime growth uses DPTMIN gamma
- Full Fairall COARE 3.0 (warm-layer/cool-skin); SEA.DAT reader
- PRECIP.DAT station rates; CLOUD.DAT I/O; IOUTMM5 variants beyond format 92
- Barriers / lake breeze; non-UTM projections; PACOUT / METLST writers
- Raw wrfout in git (rebuild from NCAR tutorial; see `cases/wrf_demo/README.md`)


## Design principles

- **Full Fortran CALMET feature + parameter-API parity.** Every INP/control parameter is accepted on `CalmetConfig` (207 keys); physics may still be Partial/TBD, but silent wrong-results for selected switches raise `NotImplementedError`. Bit-identical floats are **not** required.
- **Correctness over bit-identical Fortran ports.** Match Fortran goldens closely enough that core fields are trustworthy; do not chase bit-identical outputs.
- **Prefer NumPy / SciPy** for math kernels; **datetime** for time; **pathlib/os** for I/O; idiomatic Python over mechanical Fortran translation.
- **WRF / NetCDF:** optional `xarray` + `wrf-python` (+ `netCDF4`) for wrfout→3D.DAT and NetCDF writers (`pip install 'py-calmet[wrf]'`).
- **Regression gates** (relative RMSE / correlation in `tests/thresholds.py`) catch real regressions without over-fitting every grid cell.

See [`docs/calmet-inp-coverage.md`](docs/calmet-inp-coverage.md) for the live Implemented / Partial / Missing matrix.

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
