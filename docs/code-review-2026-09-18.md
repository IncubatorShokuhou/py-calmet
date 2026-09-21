> **Note (2026-09):** synthetic `small_domain` / `daytime_zi` cases were removed; Fortran validation is **wrf_demo-only**.

# Code review — py-calmet (main @ 722b850)

**Date:** 2026-09-18
**Scope:** `main` as published (pure-NumPy CALMET-style diagnostic). Fortran was used only to produce goldens; the tree is Python-only.
**Cases reviewed:** `small_domain`, `wrf_demo`, `daytime_zi`.
**This PR:** review document plus high-value correctness / I/O / packaging fixes. No v2 features (COARE, full clouds, IOUTMM5 variants, MIXDT sounding lapse, multi-station OA).

---

## Strengths

- Clear v1 product: diagnostic U/V + PBL fields, three run modes, CALMET.DAT / NetCDF writers, golden suites that actually run Fortran-produced binaries through the Python reader.
- DIAGNO order in `runner.run_calmet` is documented and gated on INP flags (`IFRADJ`, `ISLOPE`, `IKINE`, `IOBR`). The earlier “always minimize divergence” bug is already gone on `main`.
- Physics helpers are split into readable modules (`winds`, `pbl`, `similt`, `met_utils`) rather than one 2k-line script. Fortran sequential I/O (`fortran_bin.py` + labeled records) is careful about endian markers.
- `small_domain` relative-RMSE gates are tight enough to be a real regression net on U/V/ZI/USTAR. `daytime_zi` QSW vs Fortran (~0.012 rel RMSE) shows the Holtslag solar path is not decorative.
- No secrets, no vendored `calmet.x` / `.f90`, no raw wrfout. `.gitignore` covers `vendor/`, `*.x`, and `cases/wrf_demo/raw/`.
- Optional deps are honestly optional in `pyproject.toml` except where the runner previously *pretended* they were (see Critical lat/lon).

---

## Critical

### C1. CALMET.DAT `YYYYJJJHH` decoder drops the year ones-digit

**Where:** `py_calmet/io/calmet_dat.py` `_parse_yyyyjjjhh`

NDATHR is `year * 100000 + jday * 100 + hour` (remainder after the year is **five** digits, `JJJHH`). The decoder used `% 1_000_000`, which only works when the year ends in `0`.

| Timestamp | Encoded | Old decode | True date |
|-----------|---------|------------|-----------|
| 2020-06-15 00Z | 202016700 | 2020-06-15 | OK (year ends in 0) |
| 2005-08-28 00Z | 200524000 | **2019-05-07** | 2005-08-28 |

`wrf_demo` goldens are 2005-08-28. Reader tests only asserted `year == 2020` on `small_domain`. Time-bounds on Katrina files were silently wrong; round-trip of a 2005 writer file would not match.

**Fix in this PR:** parse `% 100_000`; reject impossible JJJ/HH; tests for 2005/1999/2020 and `wrf_demo` golden time-bounds.

### C2. Domain lat/lon silently defaulted to inland Maine

**Where:** `py_calmet/core/runner.py` `_estimate_latlon` (before this PR)

Case INPs comment out `RLAT0`/`RLON0` (`* RLAT0= 0N *`). The runner then called **optional** `pyproj`. On `ImportError` (or any transform failure) it used `(44.25°N, 70°W)`.

- `small_domain` / `daytime_zi` happen to *be* that point, so QSW golden tests would still pass without pyproj.
- `wrf_demo` is UTM 14N Mexico (~20°N). Night-only Katrina hours hide the solar error; a daytime run would get Maine insolation and Coriolis.

`pyproj` is not a runtime dependency (`pyproject.toml` lists only `numpy`).

**Fix in this PR:** WGS84 UTM inverse in `met_utils.utm_to_latlon` (no pyproj); honor `IUTMZN` / `UTMHEM`; no Maine fallback.

---

## Important

### I1. Run length from `IEHR - IBHR` only (no calendar)

**Where:** `run_calmet`

`nhrs = (IEHR - IBHR) * 3600 // NSECDT` ignores `IEYR`/`IEMO`/`IEDY`. A 22Z→02Z window becomes 1 hour (the `max(1, …)` floor) instead of 4. Julian day and sounding pick used the start date for every step.

Current goldens are same-calendar-day, so pytest did not see this.

**Fix in this PR:** `_run_window` uses full start/end datetimes; each step advances `jday` / YMD.

### I2. MIXHMC inversion jump `dptt` was computed and discarded

**Where:** `pbl.mixht_day_carson`

Maul–Carson energy balance uses previous-hour `dptt` in
`h² + 2((w'θ − wto)(1+CONSTE)Δt − dptt·h)/γ`. The new jump was calculated as `dpttp1` but never returned. The runner also never passed `dptt_prev`. That makes every hour start with `dptt = 0`, which **over-grows** ZI — consistent with the documented “faster than Fortran / DPTMIN gamma” note, but it was a second, unstated bug.

**Fix in this PR:** return `(zi, ziconv, dptt)` and persist in the runner. Absolute daytime ZI will still differ from Fortran (constant `DPTMIN` vs MIXDT sounding lapse, deferred to v2). Correlation gate remains the contract.

### I3. 3D→CALMET mapping assumed a 1-cell halo and ignored origin

**Where:** `winds.interp_3d_to_calmet` (and noobs T2/RH)

Signature takes `xorig_km, yorig_km, dgrid_km` but used `ii = min(i+1, ni-1)`. That is correct **only** when 3D.DAT is inset by one cell of the same spacing (true for both shipped goldens: `M3D_X0 = XORIG − dx`). A same-size 3D grid would shift every column and clamp the last.

**Fix in this PR:** nearest 3D mass point from cell-center vs `threed.x0_km` / `dx_km`. Equivalent on current goldens; tested with a no-halo grid.

### I4. Daytime `EL` and stored `QH` disagreed

**Where:** `run_calmet` post-DIAGNO PBL

Daytime `ustar`/`el` came from `elustr_unstable` (simple `(1−albedo)QSW/(1+Bowen)`), then `QH` was overwritten with Holtslag–van Ulden energy-budget flux. `WSTAR` used the energy-budget `QH`; `EL`/`IPGT` used the other. Night ELUSTR was internally consistent.

**Fix in this PR:** after combining `ustar` and `QH`, recompute `EL` with the CALMET constant `253.8226 = CP/(κg)` so output fields match.

### I5. `write_calmet_dat` omitted `XPSTA`/`YPSTA` when `NPSTA ≥ 1`

The reader always consumes those static records if `npsta ≥ 1`. Writer wrote `RMM`/`IPCODE` for precip but skipped station locations, so a precip-enabled file could not be read back. Shipped goldens have `NPSTA = 0`.

**Fix in this PR:** write `XPSTA`/`YPSTA`; round-trip test with `npsta=1`.

### I6. Declared cell-center coordinates were SW corners

`CalmetDataset.x` / `.y` docs say cell centers but returned `XORIG + i·Δx`. Objective analysis in the runner already uses `+ 0.5` cell. Animations and any user georeference were offset by half a cell.

**Fix in this PR:** `+ 0.5 * dgrid`.

### I7. Package version split-brain

`pyproject.toml` → `1.0.0`; `py_calmet.__version__` → `0.1.0`.

**Fix in this PR:** `__version__ = "1.0.0"`.

### I8. Docs advertised SIMILT; runner uses a power-law

`core.similt.similt_profile` is a real Van Ulden–Holtslag implementation. `winds.py` imported it and never called it. `obs_profile_similt` is power-law speed + UA direction blend; `z0`/`el`/`zi` were unused. README listed “SIMILT profiles”.

Wiring SIMILT would move `obs` / `obs_model` goldens (v1 was tuned on the power-law). **Not wired.** README corrected; unused import removed; API args kept with an explicit comment.

### I9. Meteorological-direction interpolation across 0°/360°

`np.interp` on `wd` in degrees. 350° and 10° average to 180°. Current goldens sit near 220° so gates did not catch it.

**Fix in this PR:** blend unit vectors (same construction as `wind_uv`).

### I10. Parity gates — meaningful vs loose

| Suite | Verdict |
|-------|---------|
| `small_domain` U/V/ZI/USTAR | Meaningful. Measured RMSE sits ~50–80% of the cap. Does **not** gate W, T, IPGT, EL, QSW, IRH. |
| `wrf_demo` noobs / obs_model | U/V **correlation** ≳ 0.90 is a real structure check after the DIAGNO gating fix (PROGRESS: ~0.97). Floored rel RMSE on V (cap 2.0 vs ~1.28) is loose; raw rel RMSE is still dominated by near-zero cells (documented). |
| `wrf_demo` **obs** | Weak. `U_corr ≥ 0.40` / `V_corr ≥ 0.55` with measured ~0.49 / ~0.63. Floored U/V caps 1.2 / 2.0. This mode can regress a lot before CI fails. |
| `daytime_zi` | QSW and night ZI are tight; **day ZI is correlation-only** (absolute growth is a known v2 MIXDT gap). U/V still use unfloored rel RMSE 0.05 on a night-started synthetic — OK here, would be the wrong metric on Katrina-like terrain. |

No threshold loosening in this PR. A follow-up could add EL/QSW smoke vs Fortran on `small_domain` and raise `wrf_demo` obs correlation now.

---

## Minor

- Undocumented ZI terrain fudge: `zi *= 1 + 0.002 (elev−mean)/std` (~0.2% scale). Retained for golden parity; commented in the runner.
- `slope_flow(..., rho)` ignores `rho` and uses `rhocp = 1229.9`. Froude uses `9.8` while `met_utils.G = 9.81`.
- `heat_flux_energy_budget` night QH is `−0.1` W m⁻² (placeholder) vs a real longwave budget.
- GEO.DAT reader assumes generator row order (j increasing with file line). Official MAKEGEO is often north→south; interoperability with third-party GEO files is untested.
- SURF.DAT does not map 9999 missing flags; single-station OA fills the whole grid.
- `relative_humidity_2d` is unused; noobs IRH is 3D level-1 RH, not 2 m q2 (q2 is on the 3D line but not parsed).
- `FortranSequentialWriter` is duplicated in `calmet_dat.py` instead of `fortran_bin.py`. `BinaryIO` was used unimported (annotations were deferred).
- `calmet_dat.py` shebang is `#!/bin/python3`. Duplicate `calmet.lst` / `CALMET.LST` in goldens.
- `pyproject.toml` has no license classifier, URLs, or console script for `wrfout_to_3d.py`.
- `inp.py` is a regex subset of CALMET.INP; fine for shipped files, brittle on wrapped `! KEY =` lines.
- `objective_analyze` is single-station Barnes `exp(−r²/R²)` without the two-pass γ factor (accepted v1 scope).
- `divergence_minimize(..., terrain=)` is unused. `ipgt_from_el(..., zimin=)` unused.
- Case scripts (`make_tiny_domain.py`) hardcode `/workspace/py-calmet/...` paths from the original agent workspace.

---

## Security / safety

- No API keys, tokens, or credentials in the tree.
- No Fortran sources or `calmet.x` in git. `.gitignore` includes `vendor/` and `*.x`.
- Golden `CALMET.DAT` files are Fortran unformatted output (expected). They are not executables.
- `decode_string(..., errors="ignore")` / `errors="replace"` on labels is appropriate for Fortran `CHARACTER` padding, not a sanitizer issue (local scientific files).
- `write_calmet_netcdf` / `wrfout_to_3d` require optional `netCDF4`; failure is an ImportError, not a catch-all.

---

## Maintainability

- Type hints are present on most public functions; `write_calmet_dat(result: "object")` and `threed` untyped blobs are the sore spots.
- Dead / misleading: unused `_is_daytime` (removed), unused `similt_profile` import (removed), `obs_profile_similt` name vs behavior (documented).
- Nested Python loops over `(j, i)` (and slope/Froude/smooth) are fine at 12×12; a 200×200 operational grid will want vectorization.
- Tests import `from thresholds import THRESH` via pytest’s testdir path — works, slightly unconventional.

---

## Fixes included in this PR

| ID | Change |
|----|--------|
| C1 | `YYYYJJJHH` remainder is `% 100_000` |
| C2 | Pure-NumPy UTM→lat/lon; no Maine fallback |
| I1 | Full start/end window + per-step date |
| I2 | MIXHMC returns and persists `dptt` |
| I3 | Geographic 3D nearest-neighbor; 3D.DAT uses on-disk i,j |
| I4 | EL recomputed from the QH that is stored |
| I5 | Writer emits XPSTA/YPSTA |
| I6 | Cell-center x/y |
| I7 | `__version__` = 1.0.0 |
| I8 | README + comments: power-law, not SIMILT |
| I9 | Unit-vector wind-direction blend |

Not changed (on purpose): ZI terrain fudge, SIMILT wiring, MIXDT lapse, COARE, cloud schemes, threshold tightening, vectorizing DIAGNO loops.

---

## Assessment

The library is a credible v1 diagnostic: goldens exist, the tiny-domain wind/PBL gates are real, and DIAGNO flag gating is in place. It was **not** ready to treat `wrf_demo` CALMET.DAT timestamps or non-Maine solar geometry as trustworthy, and several “CALMET-aligned” PBL/I/O pieces were incomplete in ways pytest could not see.

After the fixes in this PR, I/O time stamps, run windows, lat/lon, MIXHMC state, and 3D mapping match the documented contract more closely. Residual v1 limitations (power-law vs SIMILT, constant `DPTMIN` gamma, loose `wrf_demo` obs correlation, missing EL/QSW golden gates) should stay labeled as such rather than implied Fortran parity.

Post-fix pytest: **31 passed**. Golden RMSE on `small_domain` and `wrf_demo` is unchanged to ~3 digits (those cases are night). `daytime_zi` day-ZI correlation vs Fortran is **0.9997** (was ~0.995) after persisting `dptt`; absolute mean is still high (Python ~1582 m vs Fortran ~1074 m) because MIXDT sounding lapse remains v2.

**Recommendation:** merge the review + fixes; keep v2 items out of this tree until MIXDT / multi-station OA have their own goldens.
