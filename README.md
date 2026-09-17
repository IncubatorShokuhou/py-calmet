# py-calmet

Pure-NumPy diagnostic meteorological downscaling inspired by CALMET (CALPUFF system).

**Status:** early development. Fortran CALMET is used only as a temporary gold-standard harness; the published library is Python-only.

## Plan
1. Compile and run Fortran CALMET for three modes (obs / obs+model / no-obs) on a small public domain.
2. Reimplement the diagnostic wind + PBL core in NumPy.
3. Enforce relative-error parity tests against Fortran goldens.
4. Publish to GitHub without Fortran sources.

See `docs/superpowers/specs/2026-09-17-py-calmet-design.md`.
