# py-calmet Design

## Goal
Pure-NumPy library that reproduces CALMET's core diagnostic downscaling (winds + near-surface/PBL params) so users can downscale without running Fortran CALMET. Final public repo is Python-only; Fortran is a temporary gold-standard harness.

## Approach (chosen)
1. Compile unofficial CALMET Fortran (`vendor/calmet-fortran`).
2. Build a small public demo domain (terrain + short wrfout + few obs).
3. Run three configs to golden: obs-dominant, obs+model, no-obs.
4. Implement `py_calmet` core against those goldens with relative-error gates (~1e-3–1e-2 on core wind/PBL fields).
5. Publish to `IncubatorShokuhou/py-calmet` (public); strip Fortran before final.

## v1 scope
- GEO terrain/land-use ingest
- SURF/UPPER and/or 3D.DAT (WRF via CALWRF or equivalent writer)
- Diagnostic U/V (and W if in scope of parity)
- IPGT, USTAR, EL, ZI, WSTAR, TEMPK, RHO, QSW, IRH, RMM when available
- NetCDF (or array) output comparable to CALMET.DAT extracts

## Out of v1 (later modules)
Full cloud schemes, overwater bulk flux, all IOUTMM5 variants, GUI, every preprocessor.

## Acceptance
Same case inputs → Python vs Fortran relative error within configured thresholds on core fields for all three run modes.
