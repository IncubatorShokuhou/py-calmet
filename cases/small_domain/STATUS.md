# Small-domain CALMET demo — STATUS

**Updated:** 2026-09-17 23:02 CST

## Summary

| Mode | Directory | Result | CALMET.DAT size |
|------|-----------|--------|-----------------|
| Observation-dominant | `cases/small_domain/obs/` | **SUCCESS** (exit 0) | 142876 bytes |
| Obs + model | `cases/small_domain/obs_model/` | **SUCCESS** (exit 0) | 143016 bytes |
| No-obs / model-only | `cases/small_domain/noobs/` | **SUCCESS** (exit 0) | 138204 bytes |

Goldens archived under `cases/small_domain/goldens/{obs,obs_model,noobs}/` (CALMET.DAT + CALMET.LST + calmet.inp + README).

## Domain choice

- **Synthetic tiny domain** (not a downloaded wrfout): designed for fast Fortran goldens.
- Horizontal: **12x12** cells @ **1 km**, UTM zone **19N**, SW origin **(400.000, 4900.000) km** (~44.25N, 70.25W, inland Maine).
- Vertical: **NZ=8**, `ZFACE = 0, 20, 40, 80, 160, 300, 600, 1000, 1500` m.
- Time: **2020-06-15 00:00-03:00 UTC** (must start before 05:00 LST; CALMET enforces this).
- Datum: WGS-84. One surface + one upper-air station at domain center for obs modes.

## Data sources

| File | Source |
|------|--------|
| `shared/geo.dat` | Synthesized (GEO.DAT 2.0); LU=20 (ag), gentle slope + bump terrain |
| `shared/surf.dat` | Synthesized (SURF.DAT 2.1); 1 station, 3 hourly records |
| `shared/up.dat` | Synthesized (UP.DAT 2.1, comma-delimited); soundings 00/12Z Jun 14-16 |
| `shared/3d.dat` | Synthesized (3D.DAT **2.1**); 14x14x10, 4 hours; format follows `/workspace/calmet-docs/hrrr2calmet.py` + Fortran RDMM5 formats 62/92 |
| CALMET.INP | Adapted from EPA/demo structure in `aluislfh/py_calpuff_wrf` sample |

No public wrfout was downloaded (kept domain tiny). Generator: `cases/small_domain/scripts/make_tiny_domain.py`.

Format references:
- GEO/SURF/UP headers from unofficial CALMET sample + ramespada/calpuff make_* scripts
- 3D.DAT writer conventions from hrrr2calmet.py (adjusted iout flags to `1 1 0 0 0` -> ioutmm5=92 uncompressed)

## Compile note

- Fortran: `/workspace/py-calmet/vendor/calmet-fortran/` (ramespada/calmet unofficial)
- Binary: `/workspace/py-calmet/vendor/calmet-fortran/src/calmet.x`
- Built with gfortran + make (rank-mismatch warnings ignored, as directed)
- Case dirs symlink this binary as `./calmet.x`

## Exact commands

```bash
# Regenerate inputs (optional)
python3 cases/small_domain/scripts/make_tiny_domain.py

# Mode 1 -- obs
cd /workspace/py-calmet/cases/small_domain/obs
./calmet.x calmet.inp
# -> calmet.dat, calmet.lst

# Mode 2 -- obs + model
cd /workspace/py-calmet/cases/small_domain/obs_model
./calmet.x calmet.inp

# Mode 3 -- noobs
cd /workspace/py-calmet/cases/small_domain/noobs
./calmet.x calmet.inp
```

LCFILES=T -> all filenames lowercase on Linux (`geo.dat`, `surf.dat`, `up.dat`, `3d.dat`, `calmet.inp`).

## Key log snippets (success)

obs / obs_model / noobs each show:

```
 ENTERING SETUP PHASE
 ENTERING COMPUTATIONAL PHASE
+Processing Year, Day, Hour, Sec from: 2020 167  0    0 to:2020 167  0 3600
+Processing Year, Day, Hour, Sec from: 2020 167  1    0 to:2020 167  1 3600
+Processing Year, Day, Hour, Sec from: 2020 167  2    0 to:2020 167  2 3600
 ENTERING TERMINATION PHASE
```

LST files end with `End of run -- Clock time: ...`

## Blockers / notes resolved during setup

1. **Start time:** CALMET requires begin time before 05:00 LST -> used 00-03 UTC.
2. **Case sensitivity:** LCFILES=T needs lowercase paths on Linux.
3. **3D.DAT length:** prognostic interpolation needs NHRS+1 hours in 3D.DAT for an NHRS-hour CALMET window.
4. **iout flags:** ioutc=1 selects compressed upper format (ioutmm5>=93); used `1 1 0 0 0` (format 92) for simple synthetic profiles.
5. **NOOBS=2:** requires |IEXTRP|=1 (no surface extrapolation).

## Paths

- Executable: `/workspace/py-calmet/vendor/calmet-fortran/src/calmet.x`
- Cases: `/workspace/py-calmet/cases/small_domain/{obs,obs_model,noobs}/`
- Shared inputs: `/workspace/py-calmet/cases/small_domain/shared/`
- Goldens: `/workspace/py-calmet/cases/small_domain/goldens/`
- Status: `/workspace/py-calmet/cases/small_domain/STATUS.md`
