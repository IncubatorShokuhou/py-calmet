# py-calmet Implementation Plan

> **For agentic workers:** implement task-by-task; verify Fortran before Python rewrite.

**Goal:** Fortran CALMET goldens for three modes, then pure-NumPy `py-calmet` with parity tests, published Python-only to GitHub.

**Architecture:** Temporary Fortran vendor tree produces golden NetCDF/binary extracts; `py_calmet` reimplements the diagnostic core; `tests/` compare fields.

**Tech Stack:** gfortran, NumPy, netCDF4/xarray, pytest; optional wrf-python/netCDF for wrfout.

## Global Constraints
- Order: Fortran compile+run success → Python rewrite → parity → publish Python-only
- Public target: IncubatorShokuhou/py-calmet
- Parity: relative error ~1e-3–1e-2 on core fields (not bit-identical)
- User chose small public demo domain picked by agent

## Tasks
- [x] Compile `vendor/calmet-fortran` to `calmet.x`
- [x] Obtain/create small-domain GEO + SURF/UPPER + short wrfout/3D.DAT
- [x] Run obs / obs+model / no-obs; archive goldens under `cases/small_domain/goldens/`
- [x] Scaffold `py_calmet` package + readers/writers
- [x] Implement wind diagnostic core; parity vs goldens
- [x] Implement PBL/surface params; parity vs goldens
- [x] Docs, CI tests, push public repo; remove Fortran from published tree
