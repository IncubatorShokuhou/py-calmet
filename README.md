# py-calmet

Pure-NumPy diagnostic meteorological downscaling inspired by CALMET (CALPUFF system).

## Status

v0.0.1 — core wind + PBL parity against a tiny Fortran-CALMET golden domain
(`cases/small_domain`, 12×12 @ 1 km, 2020-06-15 00–03 UTC) for three modes:

| Mode | Description |
|------|-------------|
| `obs` | Surface + upper-air observations |
| `obs_model` | Observations + 3D.DAT initial-guess field |
| `noobs` | Model-only (3D.DAT / NOOBS=2) |

## Install

```bash
pip install -e ".[dev]"
pytest -q
```

## Quick start

```python
from py_calmet import run_calmet, read_calmet_dat

res = run_calmet("cases/small_domain/goldens/obs", mode="obs",
                 inputs_dir="cases/small_domain/goldens/inputs")
print(res.U.shape)  # (ntime, nz, ny, nx)

gold = read_calmet_dat("cases/small_domain/goldens/obs/CALMET.DAT")
```

## Package layout

- `py_calmet/io/` — GEO / SURF / UP / 3D.DAT / INP / CALMET.DAT readers
- `py_calmet/core/` — diagnostic winds (3D interp, SIMILT helpers, OA blend) + PBL (ELUSTR/MIXHT night)
- `tests/` — reader smoke tests + golden parity (thresholds in `tests/thresholds.py`)
- `cases/small_domain/goldens/` — archived Fortran outputs + shared inputs

## Parity thresholds

Relative RMSE gates are documented in `tests/thresholds.py` and `PROGRESS.md`.
Core layer-1 U/V and ZI/USTAR for `obs` are ~1e-2 or better; full 3-D winds for
`obs`/`obs_model` are best-effort approximations of CALMET DIAGNO pending a fuller port.

## License / attribution

CALMET algorithm references: Scire et al. / Exponent CALPUFF system documentation.
This repository does **not** redistribute Fortran CALMET sources.
