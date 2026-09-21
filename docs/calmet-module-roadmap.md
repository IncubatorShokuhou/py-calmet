# py-calmet module roadmap (full CALMET parity)

**Generated:** 2026-09-18 01:32 UTC (Asia/Shanghai)

Ordered plan to close gaps versus Fortran CALMET after v1 (INP coverage: 207 implemented / 0 partial / 0 missing of 207).

## Goals & non-goals

| | |
|-|-|
| **Goal** | Full CALMET feature set + complete INP parameter surface (`CalmetConfig`). |
| **Goal** | Golden regression gates (RMSE / correlation); trustworthy fields. |
| **Non-goal** | Bit-identical floats vs Fortran. |
| **Style** | Idiomatic Python; `datetime` for time; `pathlib`/`os` for files; NumPy/SciPy for math. |

## Dependencies

| Area | Libraries |
|------|-----------|
| Math kernels | `numpy`, `scipy` |
| WRF / NetCDF | `xarray`, `wrf-python`, `netCDF4` |
| Time | `datetime` (stdlib) — not Julian packing except CALMET.DAT wire I/O |
| Files | `pathlib`, `os` (stdlib) |

`scripts/wrfout_to_3d.py` and NetCDF writers should prefer **xarray + wrf-python** for field extraction/projection; keep a thin 3D.DAT binary writer in `py_calmet.io`.

## Ordered work packages

### 0. INP / config foundation (prerequisite)

- Expand `io/inp.py` → typed `CalmetConfig` covering all CVDIC keys + defaults from Fortran BLOCK DATA.
- Honor Group 0 filenames; fix `JWAT1/2` ↔ `IWAT1/2` alias.
- Bind hardcoded `HA1…HC3` / water LU range to config.
- Use `datetime` for start/end/`NSECDT` stepping in `runner`.
- **Exit:** every INP key either drives code or is explicitly stubbed with a clear `NotImplemented` path.

### 1. DIAGNO full (winds core)

| Item | INP / routines | Notes |
|------|----------------|-------|
| Multi-station Barnes OA | `R1`,`R2`,`RMAX*`,`RMIN*`,`NINTR2`,`LVARY`,`NSSTA` | Replace single-station OA. |
| Prog blending | `RPROG`,`IPROG`,`IGFMET` | Weight IGF vs obs; IGF-CALMET reader. |
| Extrapolation | `IEXTRP`,`BIAS`,`FEXTR2` | Wire `similt.similt_profile`; layer bias. |
| Barriers / lake breeze | `NBAR`,`KBAR`,`X*BAR`,`LLBREZE`,`NBOX`,… | Terrain barriers + breeze boxes. |
| Diag options | `IDIOPT1–5`,`ISURFT`,`IUPT`,`IUPWND`,`ZUP*` | Terrain-circulation T/lapse sources. |
| Calm / div criterion | `ICALM`,`DIVLIM` | Match Fortran calm handling. |

### 2. IKINE — kinematic topographic vertical velocity

- Replace `light_terrain_adjust` with Fortran-aligned `TOPOF2`-style kinematic `W` from terrain slope × wind.
- Feed horizontal adjustment consistently with `ALPHA`.
- Golden: enable `IKINE=1` on `wrf_demo` slope/terrain case.

### 3. IOBR — full O'Brien adjustment

- Upgrade `divergence_minimize` toward O'Brien vertical velocity / horizontal div cleanup (`IOBR=1`).
- Honor `DIVLIM`, `NITER`, layer coupling with kinematic `W`.
- Keep gating: do nothing when `IOBR=0` (current correct behavior).

### 4. PBL / temperature refinements

- Sounding-based lapse above Zi (`MIXDT` / `DPTMIN` path) — fix daytime Zi growth bias.
- Full `IMIXH` options; `IAVEZI`/`MNMDAV`/`HAFANG`/`ILEVZI`; `IZICRLX`/`TZICRLX`.
- `ITPROG` / `ITWPROG` / `TGDEF*` / `TRADKM` / `IAVET` / `ILUOC3D` temperature fields.
- Read `JWAT*`, radiation coeffs from INP.

### 5. Clouds

- Implement `ICLOUD` / `MCLOUD` methods 1–4 (ceilometer, RH-based, 3D.DAT cloud, etc.).
- `CLDDAT` reader + `ICLDOUT` / `IFORMC`.
- Couple QSW / energy budget to chosen cloud field.

### 6. Overwater / COARE

- `SEA.DAT` / `NOWSTA` reader.
- COARE fluxes (`ICOARE`,`DSHELF`,`IWARM`,`ICOOL`) and overwater Zi (`ZIMINW`,`ZIMAXW`,`THRESHW`,`CONSTW`).
- Distance-to-coast (`LDBCST`,`DCSTGD`) if needed for shelf scaling.

### 7. Precipitation

- `PRECIP.DAT` + `PS*` stations (`NPSTA`,`IFORMP`,`NFLAGP`,`SIGMAP`,`CUTP`).
- Gridded RMM in CALMET.DAT (currently zeroed).

### 8. IOUTMM5 / 3D.DAT variants

- Extend `io/threed.py` beyond uncompressed format **92** (81/82/91/93–95 families).
- Prefer **xarray + wrf-python** in `wrfout_to_3d` for moisture/cloud/ice/graupel flags that feed `IOUTMM5`.
- QA: `ISTEPPGS` multiple of `NSECDT`.

### 9. Writers & outputs

| Output | Status / work |
|--------|----------------|
| CALMET.DAT | Exists; honor `METDAT`, `LSAVE`, `LCALGRD`, dates/UTM latlon already fixed in v1.0 |
| NetCDF | Exists; prefer xarray for structure; keep netCDF4/xarray write path |
| PACOUT.DAT | Missing (`IFORMO=2`) |
| CALMET.LST | Missing (`METLST`, echo of inputs like Fortran) |
| Test/kin/frd/slp dumps | Missing (`TST*`) — low priority |
| Printer layer flags | Missing (`IUVOUT`/…) — optional |

### 10. Projections & multi-file

- Non-UTM `PMAP` (LCC, PS, EM, TTM, LAZA) + `DATUM`/`FEAST`/`FNORTH`/`XLAT*`.
- Multiple `M3DDAT` / `IGFDAT` / `UPDAT` lists (`NM3D`,`NIGF`,`NUSTA`).

## Suggested milestone order

1. **v1.1** — Config surface + JWAT/filename/radiation binding + SIMILT wired (`IEXTRP`).
2. **v1.2** — Full DIAGNO OA (multi-station) + MIXDT lapse.
3. **v1.3** — IKINE + IOBR production paths + goldens.
4. **v2.0** — Clouds + precip writers.
5. **v2.1** — COARE / overwater.
6. **v2.2** — IOUTMM5 variants + wrfout_to_3d via xarray/wrf-python; PACOUT/LST.

## Mapping to current modules

| Fortran area | Current Python | Next touch |
|--------------|----------------|------------|
| READCF/READFN | `io/inp.py` (minimal) | typed config |
| DIAGNO / WIND1 | `core/winds.py` | OA, IKINE, IOBR |
| SIMILT | `core/similt.py` (unwired) | wire via IEXTRP |
| ELUSTR / MIXH* | `core/pbl.py` | MIXDT, IMIXH, COARE hook |
| COMP / runner | `core/runner.py` | datetime loop, config-driven |
| RDMM5 / 3D.DAT | `io/threed.py` | IOUTMM5 variants |
| wrfout bridge | `scripts/wrfout_to_3d.py` | xarray + wrf-python |
| OUTHD/OUTHR | `io/calmet_dat.py` | PACOUT, LST, flags |

## Related

- [`calmet-inp-coverage.md`](calmet-inp-coverage.md) — per-variable status
- [`calmet-inp-params.json`](calmet-inp-params.json) — machine-readable list
- [`../PROGRESS.md`](../PROGRESS.md) — measured golden status
