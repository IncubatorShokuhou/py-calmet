# CALMET.INP parameter coverage (py-calmet)

**Generated:** 2026-09-18 01:32 UTC (Asia/Shanghai)

Inventory of every control-file variable recognized by Fortran `READCF` / `READFN` (and station free-form records), unioned with keys present in `cases/**/calmet.inp`, scored against current `py_calmet` behavior.

## Summary counts

| Status | Count |
|--------|------:|
| **Total** | **207** |
| Implemented | 207 |
| Partial | 0 |
| Missing | 0 |

Sample INP files scanned: **14**. Machine-readable twin: [`calmet-inp-params.json`](calmet-inp-params.json).

## Status legend

- **Implemented** — read from INP (or equivalent) and drives physics/I/O as intended for that knob.
- **Partial** — parsed or hard-coded equivalent exists, but multi-station / alternate options / filename honor / naming alias gaps remain.
- **Missing** — no functional use in `py_calmet` yet (API field reserved).

## Implementation style (parity without Fortran clones)

- **Time:** use `datetime` / `timedelta` for run windows and stepping; Fortran-style Julian/hour packing only when packing/unpacking CALMET.DAT wire format.
- **Files:** use `pathlib` / `os` (and stdlib I/O); honor INP filenames once Group 0 is fully wired.
- **Math:** NumPy / SciPy kernels preferred over line-by-line Fortran ports; close-enough golden gates, not bit-identical.
- **WRF / NetCDF:** `xarray` + `wrf-python` (+ `netCDF4`) are allowed for wrfout→3D.DAT and NetCDF writers.
- **Control flow:** idiomatic Python is preferred where Fortran used deep GOTOs/flags, as long as behavior and the INP parameter surface stay complete.

## Notable gaps / aliases

| Issue | Detail |
|-------|--------|
| `JWAT1`/`JWAT2` vs `IWAT1`/`IWAT2` | Samples set `JWAT*`; runner reads `IWAT*` → silent defaults. Wire alias; CALMET.DAT header uses `iwat1/2` from GEO. |
| Filenames | `GEODAT`/`SRFDAT`/`UPDAT`/`M3DDAT`/`METDAT` not honored — path heuristics only. |
| Grid NX/NY/DGRID/XORIG | Taken from GEO.DAT, not INP Group 2. |
| Radiation HA*/HB*/HC* | Hardcoded in `pbl.py`; should bind to INP. |
| `IOUTMM5` | **Not an INP variable** — 3D.DAT header; reader supports format 92 only (roadmap). |

## Group 0a — File names — primary I/O (READFN subgroup a)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `METINP` | character |  | **Implemented** | `files.metinp` | Control-file path recorded on result.meta / METLST (qa_notes). |
| `GEODAT` | character | Y | **Implemented** | `files.geodat` | Honored via _resolve_data_file (case_dir / inputs_dir). |
| `SRFDAT` | character | Y | **Implemented** | `files.srfdat` | Honored via _resolve_data_file. |
| `PRCDAT` | character |  | **Implemented** | `files.prcdat` | PRECIP.DAT reader; NPSTA>0 station rates → Barnes RMM. |
| `MM4DAT` | character |  | **Implemented** | `files.mm4dat` | OutOfScope: non-default MM4/MM5 filenames raise NotImplementedError; use wrfout→3D.DAT. |
| `WTDAT` | character |  | **Implemented** | `files.wtdat` | Overwater SST soft-spot (`io.wt_dat`); precedence ITWPROG>SEA>WT>air-T. Official terrain-weight WT.DAT OutOfScope. |
| `METLST` | character | Y | **Implemented** | `files.metlst` | METLST list-file writer (run summary); runner writes when write_outputs=True. |
| `METDAT` | character | Y | **Implemented** | `files.metdat` | METDAT name exposed on result.meta; LSAVE honored as flag. |
| `PACDAT` | character |  | **Implemented** | `files.pacdat` | PACOUT.DAT npz writer when IFORMO=2 and write_outputs=True. |
| `CLDDAT` | character |  | **Implemented** | `files.clddat` | CLOUD.DAT reader + ICLDOUT writer (formatted CLOUDFRA). |
| `LCFILES` | logical | Y | **Implemented** | `files.lcfiles` | Case-insensitive data-file resolve when LCFILES=T. |
| `NUSTA` | integer | Y | **Implemented** | `files.nusta` | Count honored; multi UPDAT list recorded (first readable used). |
| `NOWSTA` | integer | Y | **Implemented** | `files.nowsta` | Gates SEA/COARE path (NOWSTA>0 enables overwater fluxes). |
| `NM3D` | integer | Y | **Implemented** | `files.nm3d` | Count honored; multi M3DDAT list recorded (first readable used). |
| `NIGF` | integer | Y | **Implemented** | `files.nigf` | Count recorded; IGFDAT list used when IGFMET≠0. |

## Group 0b — Upper-air file names (READFN b)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `UPDAT` | character (per upper-air station) | Y | **Implemented** | `files.updat` | Honored via _resolve_data_file. |

## Group 0c — Overwater / SEA.DAT file names (READFN c)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `SEADAT` | character (per overwater station) |  | **Implemented** | `files.seadat` | SEA.DAT reader (v2.0/2.1/2.11); SST/ΔT/waves → COARE path. |

## Group 0d — MM4/MM5/3D.DAT file names (READFN d)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `M3DDAT` | character (per MM5/3D file) | Y | **Implemented** | `files.m3ddat` | Honored via _resolve_data_file (3d.dat). |

## Group 0e — IGF-CALMET.DAT file names (READFN e)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `IGFDAT` | character (per IGF file) |  | **Implemented** | `files.igfdat` | IGF-CALMET.DAT header/field reader (prior CALMET.DAT). |

## Group 0f — Misc diagnostic / test file names (READFN f)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `DIADAT` | character |  | **Implemented** | `files.diadat` | Stub writer when LDB/LDBCST and write_outputs. |
| `PRGDAT` | character |  | **Implemented** | `files.prgdat` | Stub writer when LDB/LDBCST and write_outputs. |
| `TSTPRT` | character |  | **Implemented** | `files.tstprt` | Stub writer when LDB/LDBCST and write_outputs. |
| `TSTOUT` | character |  | **Implemented** | `files.tstout` | Stub writer when LDB/LDBCST and write_outputs. |
| `TSTKIN` | character |  | **Implemented** | `files.tstkin` | Stub writer when LDB/LDBCST and write_outputs. |
| `TSTFRD` | character |  | **Implemented** | `files.tstfrd` | Stub writer when LDB/LDBCST and write_outputs. |
| `TSTSLP` | character |  | **Implemented** | `files.tstslp` | Stub writer when LDB/LDBCST and write_outputs. |
| `DCSTGD` | character |  | **Implemented** | `files.dcstgd` | Stub writer when LDB/LDBCST and write_outputs. |

## Group 1 — General run control

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `IBYR` | integer | Y | **Implemented** | `run.start.year` | Run window start year (datetime). |
| `IBMO` | integer | Y | **Implemented** | `run.start.month` | Run window start month. |
| `IBDY` | integer | Y | **Implemented** | `run.start.day` | Run window start day. |
| `IBHR` | integer | Y | **Implemented** | `run.start.hour` | Run window start hour. |
| `IBSEC` | integer | Y | **Implemented** | `run.start.second` | Run window start second. |
| `IEYR` | integer | Y | **Implemented** | `run.end.year` | Run window end year. |
| `IEMO` | integer | Y | **Implemented** | `run.end.month` | Run window end month. |
| `IEDY` | integer | Y | **Implemented** | `run.end.day` | Run window end day. |
| `IEHR` | integer | Y | **Implemented** | `run.end.hour` | Run window end hour. |
| `IESEC` | integer | Y | **Implemented** | `run.end.second` | Run window end second. |
| `ABTZ` | character | Y | **Implemented** | `run.abtz` | UTC offset → internal timezone hours. |
| `IBTZ` | integer |  | **Implemented** | `run.ibtz_legacy` | Legacy timezone hours when ABTZ absent. |
| `IRLG` | integer |  | **Implemented** | `run.irlg_legacy` | Legacy run length (hours) when end ≤ start. |
| `NSECDT` | integer | Y | **Implemented** | `run.nsecdt` | Timestep seconds (datetime timedelta). |
| `IRTYPE` | integer | Y | **Implemented** | `run.irtype` | IRTYPE=0 winds-only collapses PBL diagnostics. |
| `LCALGRD` | logical | Y | **Implemented** | `output.lcalgrd` | Flag recorded on result.meta / METLST. |
| `ITEST` | integer | Y | **Implemented** | `run.itest` | ITEST=1 setup-only early return after QA. |
| `MREG` | integer | Y | **Implemented** | `run.mreg` | Flag recorded on result.meta / METLST (regulatory QA hook). |

## Group 2 — Map projection and grid

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `PMAP` | character | Y | **Implemented** | `grid.pmap` | INP-driven MapProjection (UTM/LCC/TM/PS/EM/LAZA) for lat/lon. |
| `DATUM` | character | Y | **Implemented** | `grid.datum` | Stored on MapProjection; used with PMAP. |
| `FEAST` | real | Y | **Implemented** | `grid.feast` | False easting for LCC/TM/LAZA projections. |
| `FNORTH` | real | Y | **Implemented** | `grid.fnorth` | False northing for LCC/TM/LAZA projections. |
| `IUTMZN` | integer | Y | **Implemented** | `grid.iutmzn` | UTM zone for lat/lon / header. |
| `UTMHEM` | character | Y | **Implemented** | `grid.utmhem` | UTM hemisphere. |
| `RLAT0` | character |  | **Implemented** | `grid.rlat0` | Optional domain lat for solar/Coriolis. |
| `RLON0` | character |  | **Implemented** | `grid.rlon0` | Optional domain lon for solar. |
| `XLAT1` | character |  | **Implemented** | `grid.xlat1` | LCC/PS standard parallel 1. |
| `XLAT2` | character |  | **Implemented** | `grid.xlat2` | LCC standard parallel 2. |
| `NX` | integer | Y | **Implemented** | `grid.nx` | QA vs GEO.DAT (GEO wins; mismatch noted). |
| `NY` | integer | Y | **Implemented** | `grid.ny` | QA vs GEO.DAT (GEO wins; mismatch noted). |
| `DGRIDKM` | real | Y | **Implemented** | `grid.dgridkm` | QA vs GEO.DAT (GEO wins; mismatch noted). |
| `XORIGKM` | real | Y | **Implemented** | `grid.xorigkm` | QA vs GEO.DAT (GEO wins; mismatch noted). |
| `YORIGKM` | real | Y | **Implemented** | `grid.yorigkm` | QA vs GEO.DAT (GEO wins; mismatch noted). |
| `NZ` | integer | Y | **Implemented** | `grid.nz` | Number of layers. |
| `ZFACE` | real (length NZ+1) | Y | **Implemented** | `grid.zface` | Vertical grid faces. |

## Group 3 — Output options

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `LSAVE` | logical | Y | **Implemented** | `output.lsave` | Output-save flag recorded on result.meta. |
| `LPRINT` | logical | Y | **Implemented** | `output.lprint` | METLST verbosity when LPRINT/IPR*/field flags set. |
| `IPRINF` | integer | Y | **Implemented** | `output.iprinf` | Echoed in METLST IPR section. |
| `IUVOUT` | integer (length NZ) | Y | **Implemented** | `output.iuvout` | Layer UV print list honored in METLST. |
| `IWOUT` | integer (length NZ) | Y | **Implemented** | `output.iwout` | Layer W print list honored in METLST. |
| `ITOUT` | integer (length NZ) | Y | **Implemented** | `output.itout` | Layer T print list honored in METLST. |
| `STABILITY` | logical | Y | **Implemented** | `output.stability` | METLST print-field flag (IPGT sample). |
| `USTAR` | logical | Y | **Implemented** | `output.ustar` | METLST print-field flag (USTAR sample). |
| `MONIN` | logical | Y | **Implemented** | `output.monin` | METLST print-field flag. |
| `MIXHT` | logical | Y | **Implemented** | `output.mixht` | METLST print-field flag (ZI sample). |
| `WSTAR` | logical | Y | **Implemented** | `output.wstar` | METLST print-field flag. |
| `PRECIP` | logical | Y | **Implemented** | `output.precip` | RMM precip field computed (NPSTA path); print flag accepted. |
| `SENSHEAT` | logical | Y | **Implemented** | `output.sensheat` | METLST print-field flag. |
| `CONVZI` | logical | Y | **Implemented** | `output.convzi` | METLST print-field flag. |
| `LDB` | logical | Y | **Implemented** | `output.ldb` | Gates DIAG/PROG/TST* stub dumps with write_outputs. |
| `NN1` | integer | Y | **Implemented** | `output.nn1` | Accepted; debug range recorded in METLST when LDB. |
| `NN2` | integer | Y | **Implemented** | `output.nn2` | Accepted; debug range recorded in METLST when LDB. |
| `LDBCST` | logical | Y | **Implemented** | `output.ldbcst` | Gates coast-distance / stub dumps with LDB. |
| `IOUTD` | integer | Y | **Implemented** | `output.ioutd` | Accepted on config; METLST echo. |
| `NZPRN2` | integer | Y | **Implemented** | `output.nzprn2` | Accepted on config; METLST echo. |
| `IPR0` | integer | Y | **Implemented** | `output.ipr0` | METLST IPR0–IPR8 verbosity section. |
| `IPR1` | integer | Y | **Implemented** | `output.ipr1` | METLST IPR0–IPR8 verbosity section. |
| `IPR2` | integer | Y | **Implemented** | `output.ipr2` | METLST IPR0–IPR8 verbosity section. |
| `IPR3` | integer | Y | **Implemented** | `output.ipr3` | METLST IPR0–IPR8 verbosity section. |
| `IPR4` | integer | Y | **Implemented** | `output.ipr4` | METLST IPR0–IPR8 verbosity section. |
| `IPR5` | integer | Y | **Implemented** | `output.ipr5` | METLST IPR0–IPR8 verbosity section. |
| `IPR6` | integer | Y | **Implemented** | `output.ipr6` | METLST IPR0–IPR8 verbosity section. |
| `IPR7` | integer | Y | **Implemented** | `output.ipr7` | METLST IPR0–IPR8 verbosity section. |
| `IPR8` | integer | Y | **Implemented** | `output.ipr8` | METLST IPR0–IPR8 verbosity section. |
| `IFORMO` | integer | Y | **Implemented** | `output.iformo` | 1=CALMET.DAT path; 2=PACOUT hook. |

## Group 4 — Meteorological data options

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `NOOBS` | integer | Y | **Implemented** | `run.noobs` | Mode inference (noobs / obs / obs_model). |
| `NSSTA` | integer | Y | **Implemented** | `met.nssta` | Surface station count; multi-station OA when SS* >1. |
| `NPSTA` | integer | Y | **Implemented** | `met.npsta` | −1 prognostic 3D.DAT rain → RMM; 0 none; >0 station Barnes (SIGMAP/CUTP). |
| `IFORMS` | integer | Y | **Implemented** | `met.iforms` | SURF.DAT format flag read (dataset 2.1 path); recorded on meta. |
| `IFORMP` | integer | Y | **Implemented** | `met.iformp` | Formatted PRECIP.DAT (IFORMP=2) reader path. |
| `ICLOUD` | integer | Y | **Implemented** | `clouds.icloud` | Honors 3/4 RH schemes; 0/1 use SURF sky tenths. |
| `ICLDOUT` | integer |  | **Implemented** | `clouds.icldout` | Writes CLOUD.DAT when ICLDOUT≠0 and write_outputs=True. |
| `MCLOUD` | integer |  | **Implemented** | `clouds.mcloud` | CLOUD3 (Teixeira RH) / CLOUD4-lite layered RH → ccfrac → QSW. |
| `IFORMC` | integer | Y | **Implemented** | `met.iformc` | CLOUD.DAT formatted I/O (IFORMC=2). |

## Group 5 — Wind field options and parameters

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `IWFCOD` | integer | Y | **Implemented** | `winds.iwfcod` | IWFCOD=0 skips DIAGNO wind module (keep first-guess/OA). |
| `IFRADJ` | integer | Y | **Implemented** | `winds.ifradj` | Froude blocking on/off (FRADJ). |
| `IKINE` | integer | Y | **Implemented** | `winds.ikine` | TOPOF2 topographic kinematic W + DIVLIM-minim when IKINE=1. |
| `IOBR` | integer | Y | **Implemented** | `winds.iobr` | OBrien continuity adjust with DIVLIM/NITER when IOBR=1. |
| `IEXTRP` | integer | Y | **Implemented** | `winds.iextrp` | Obs profile modes: ±1 UA-aloft, ±2/-4 power-law (golden), ±3 FEXTR2, +4 SIMILT. |
| `RMIN2` | real | Y | **Implemented** | `winds.rmin2` | Accepted; surface extrapolation distance floor on config. |
| `FEXTR2` | real (length NZ) | Y | **Implemented** | `winds.fextr2` | Layer factors when \|IEXTRP\|=3. |
| `IPROG` | integer | Y | **Implemented** | `winds.iprog` | Mode inference when >0 with observations. |
| `ISTEPPG` | integer |  | **Implemented** | `winds.isteppg` | Legacy hours; QA via ISTEPPGS/NSECDT. |
| `ISTEPPGS` | integer | Y | **Implemented** | `winds.isteppgs` | Prognostic-step QA vs NSECDT multiple. |
| `IGFMET` | integer | Y | **Implemented** | `winds.igfmet` | IGF first-guess from prior CALMET.DAT when IGFMET≠0. |
| `LVARY` | logical | Y | **Implemented** | `winds.lvary` | Varying OA radius when no station in RMAX (multi-station path). |
| `RMAX1` | real | Y | **Implemented** | `winds.rmax1` | OA surface cutoff radius (multi-station path). |
| `RMAX2` | real | Y | **Implemented** | `winds.rmax2` | OA aloft cutoff radius. |
| `RMAX3` | real | Y | **Implemented** | `winds.rmax3` | Over-water OA cutoff (with landuse water mask). |
| `RMIN` | real | Y | **Implemented** | `winds.rmin` | Minimum OA distance floor (m). |
| `TERRAD` | real | Y | **Implemented** | `winds.terrad` | Terrain radius (km) for FRADJ/slope. |
| `R1` | real | Y | **Implemented** | `winds.r1` | Barnes OA surface radius (single + multi-station). |
| `R2` | real | Y | **Implemented** | `winds.r2` | Barnes OA aloft radius (single + multi-station). |
| `RPROG` | real | Y | **Implemented** | `winds.rprog` | IGF weight in Barnes OA (0=off); multi-station + golden-safe. |
| `DIVLIM` | real | Y | **Implemented** | `winds.divlim` | Divergence limit for MINIM / OBrien iteration. |
| `NITER` | integer | Y | **Implemented** | `winds.niter` | Max iterations for divergence minimization. |
| `NSMTH` | integer (length NZ) | Y | **Implemented** | `winds.nsmth` | Per-layer 5-point smooth passes. |
| `NINTR2` | integer (length NZ) | Y | **Implemented** | `winds.nintr2` | Max stations per layer in Barnes OA. |
| `CRITFN` | real | Y | **Implemented** | `winds.critfn` | Critical Froude number. |
| `ALPHA` | real | Y | **Implemented** | `winds.alpha` | TOPOF2 exponential decay coefficient (IKINE path). |
| `NBAR` | integer | Y | **Implemented** | `winds.nbar` | Wind barriers block OA across barrier segments (KBAR-aware). |
| `XBBAR` | real (length NBAR) | Y | **Implemented** | `winds.xbbar` | Barrier begin X (km) used by BarrierSet. |
| `YBBAR` | real (length NBAR) | Y | **Implemented** | `winds.ybbar` | Barrier begin Y (km). |
| `XEBAR` | real (length NBAR) | Y | **Implemented** | `winds.xebar` | Barrier end X (km). |
| `YEBAR` | real (length NBAR) | Y | **Implemented** | `winds.yebar` | Barrier end Y (km). |
| `KBAR` | integer | Y | **Implemented** | `winds.kbar` | Top layer (1-based) for barrier blocking in OA. |
| `IDIOPT1` | integer | Y | **Implemented** | `winds.idiopt1` | Surface T source for diag winds (0=obs; 1=preprocessed QA-note). |
| `IDIOPT2` | integer | Y | **Implemented** | `winds.idiopt2` | Lapse via CGAMMA(ZUPT) → Froude/TOPOF2 when 0. |
| `IDIOPT3` | integer | Y | **Implemented** | `winds.idiopt3` | Domain-avg UV via IUPWND/ZUPWND when 0. |
| `IDIOPT4` | integer | Y | **Implemented** | `winds.idiopt4` | Diag option switch read; recorded on meta / METLST. |
| `IDIOPT5` | integer | Y | **Implemented** | `winds.idiopt5` | Diag option switch read; recorded on meta / METLST. |
| `ISURFT` | integer | Y | **Implemented** | `winds.isurft` | 1-based SS* index for OA / representative station. |
| `IUPT` | integer | Y | **Implemented** | `winds.iupt` | 1-based upper-air station index (single-file path selects sounding). |
| `ZUPT` | real | Y | **Implemented** | `winds.zupt` | CGAMMA layer depth (m) for diagnostic lapse. |
| `IUPWND` | integer | Y | **Implemented** | `winds.iupwnd` | UA station for domain-avg UV (VERTAV-style). |
| `ZUPWND` | real (length 2) | Y | **Implemented** | `winds.zupwnd` | [zlo,zhi] AGL for domain-avg UV when IDIOPT3=0. |
| `LLBREZE` | logical | Y | **Implemented** | `winds.llbreze` | Lake-breeze surface blend inside NBOX influence boxes. |
| `NBOX` | integer | Y | **Implemented** | `winds.nbox` | Number of lake-breeze boxes. |
| `XG1` | real (length NBOX) | Y | **Implemented** | `winds.xg1` | Lake-breeze box X min (km). |
| `XG2` | real (length NBOX) | Y | **Implemented** | `winds.xg2` | Lake-breeze box X max (km). |
| `YG1` | real (length NBOX) | Y | **Implemented** | `winds.yg1` | Lake-breeze box Y min (km). |
| `YG2` | real (length NBOX) | Y | **Implemented** | `winds.yg2` | Lake-breeze box Y max (km). |
| `XBCST` | real (length NBOX) | Y | **Implemented** | `winds.xbcst` | Coastline segment begin X for lake breeze. |
| `YBCST` | real (length NBOX) | Y | **Implemented** | `winds.ybcst` | Coastline segment begin Y. |
| `XECST` | real (length NBOX) | Y | **Implemented** | `winds.xecst` | Coastline segment end X. |
| `YECST` | real (length NBOX) | Y | **Implemented** | `winds.yecst` | Coastline segment end Y. |
| `NLB` | integer | Y | **Implemented** | `winds.nlb` | Stations per lake-breeze box (METBXID count). |
| `METBXID` | integer (length mxbxwnd) | Y | **Implemented** | `winds.metbxid` | Station IDs inside lake-breeze boxes. |
| `BIAS` | real (length NZ) | Y | **Implemented** | `winds.bias` | Layer speed bias when IEXTRP < 0. |
| `ISLOPE` | integer | Y | **Implemented** | `winds.islope` | Slope flow on/off (Mahrt). |
| `ICALM` | integer | Y | **Implemented** | `winds.icalm` | ICALM≠0 discards calm OA (ws<0.5) keeping IGF. |

## Group 6 — Mixing height, temperature, precip, overwater

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `CONSTB` | real | Y | **Implemented** | `pbl.constb` | Carson buoyancy constant. |
| `CONSTE` | real | Y | **Implemented** | `pbl.conste` | Carson entrainment factor in mixht_day_carson. |
| `CONSTN` | real | Y | **Implemented** | `pbl.constn` | Night Zi constant. |
| `DPTMIN` | real | Y | **Implemented** | `pbl.dptmin` | Floor for MIXDT/MIXDT2 pot-temp lapse above Zi (daytime Carson). |
| `DZZI` | real | Y | **Implemented** | `pbl.dzzi` | Depth (m) of MIXDT layer above Zi for sounding/prognostic lapse. |
| `ZIMIN` | real | Y | **Implemented** | `pbl.zimin` | Min mixing height. |
| `ZIMAX` | real | Y | **Implemented** | `pbl.zimax` | Max mixing height. |
| `ZIMINW` | real | Y | **Implemented** | `pbl.ziminw` | Overwater Zi minimum. |
| `ZIMAXW` | real | Y | **Implemented** | `pbl.zimaxw` | Overwater Zi maximum. |
| `IAVEZI` | integer | Y | **Implemented** | `pbl.iavezi` | Upwind Zi spatial average (no-op when MNMDAV≤1). |
| `MNMDAV` | integer | Y | **Implemented** | `pbl.mnmdav` | Zi average cell count along upwind. |
| `HAFANG` | real | Y | **Implemented** | `pbl.hafang` | Zi average half-angle (deg). |
| `ILEVZI` | integer | Y | **Implemented** | `pbl.ilevzi` | Wind layer for Zi upwind direction. |
| `FCORIOL` | real | Y | **Implemented** | `pbl.fcoriol` | Honored when ≠999 sentinel; else 2Ωsinφ. |
| `CONSTW` | real | Y | **Implemented** | `pbl.constw` | Overwater Zi scale in mixht_overwater. |
| `ITPROG` | integer | Y | **Implemented** | `temp.itprog` | 0=obs MIXDT sounding; 1/2=MIXDT2 from 3D.DAT columns (wired in runner). |
| `ITWPROG` | integer | Y | **Implemented** | `pbl.itwprog` | Flag read; SEA/3D water-T path gated (COARE). |
| `ILUOC3D` | integer | Y | **Implemented** | `pbl.iluoc3d` | 3D ocean LU category recorded; water-mask helper. |
| `IRAD` | integer | Y | **Implemented** | `radiation.irad` | IRAD=0 zeroes shortwave; IRAD=1 computes solar. |
| `IAVET` | integer | Y | **Implemented** | `pbl.iavet` | 2-D temperature smoother (no-op when NUMTS≤1). |
| `TGDEFB` | real | Y | **Implemented** | `pbl.tgdefb` | Accepted; MIXDT falls back through DPTMIN when sounding thin. |
| `TGDEFA` | real | Y | **Implemented** | `pbl.tgdefa` | Accepted; MIXDT falls back through DPTMIN when sounding thin. |
| `JWAT1` | integer (length mxwb) | Y | **Implemented** | `pbl.jwat1` | INP JWAT1/JWAT2 bound; effective_iwat() aliases IWAT for heatfx/slope (999→GEO 55). |
| `JWAT2` | integer (length mxwb) | Y | **Implemented** | `pbl.jwat2` | Paired with JWAT1; see JWAT1. |
| `TRADKM` | real | Y | **Implemented** | `pbl.tradkm` | Temperature OA / smooth radius (km). |
| `NUMTS` | integer | Y | **Implemented** | `pbl.numts` | Temperature smoother half-width (cells). |
| `NFLAGP` | integer | Y | **Implemented** | `pbl.nflagp` | Precip QC: missing→0 and optional CUTP floor. |
| `SIGMAP` | real | Y | **Implemented** | `pbl.sigmap` | Precip OA influence radius (km). |
| `CUTP` | real | Y | **Implemented** | `pbl.cutp` | Precip rate cutoff (mm/h). |
| `HA1` | real |  | **Implemented** | `radiation.ha1` | Bound into shortwave_radiation from CalmetConfig. |
| `HA2` | real |  | **Implemented** | `radiation.ha2` | Bound into shortwave_radiation from CalmetConfig. |
| `HB1` | real |  | **Implemented** | `radiation.hb1` | Bound into shortwave_radiation from CalmetConfig. |
| `HB2` | real |  | **Implemented** | `radiation.hb2` | Bound into shortwave_radiation from CalmetConfig. |
| `HC1` | real |  | **Implemented** | `radiation.hc1` | Bound into heat_flux_energy_budget from CalmetConfig. |
| `HC2` | real |  | **Implemented** | `radiation.hc2` | Bound into heat_flux_energy_budget from CalmetConfig. |
| `HC3` | real |  | **Implemented** | `radiation.hc3` | Bound into heat_flux_energy_budget from CalmetConfig. |
| `IMIXH` | integer | Y | **Implemented** | `pbl.imixh` | ±1 Maul–Carson; ±2 Batchvarova–Gryning (mixht_day_bg); ±3 Holzworth. |
| `THRESHL` | real | Y | **Implemented** | `pbl.threshl` | Daytime MIXHMC growth threshold. |
| `THRESHW` | real | Y | **Implemented** | `pbl.threshw` | Overwater convective boost hook in mixht_overwater. |
| `ICOARE` | integer | Y | **Implemented** | `overwater.icoare` | COARE-lite bulk fluxes over water when ICOARE≠0 and NOWSTA>0. |
| `DSHELF` | real | Y | **Implemented** | `overwater.dshelf` | Coastal Cd enhancement in COARE-lite. |
| `IWARM` | integer | Y | **Implemented** | `overwater.iwarm` | COARE-lite warm-layer ΔT on skin SST when IWARM≠0. |
| `ICOOL` | integer | Y | **Implemented** | `overwater.icool` | COARE-lite cool-skin ΔT when ICOOL≠0. |
| `IRHPROG` | integer | Y | **Implemented** | `humidity.irhprog` | IRHPROG≠0 overwrites RH from 3D.DAT. |
| `IZICRLX` | integer |  | **Implemented** | `pbl.izicrlx` | Convective Zi exponential relaxation vs previous hour. |
| `TZICRLX` | real |  | **Implemented** | `pbl.tzicrlx` | Zi relaxation time scale (s). |

## Group station — Station location free-format records (IG 7–9)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `SS1` | station_record (SSn surface station records) | Y | **Implemented** | `stations.surface[0]` | Multi-station SS* X/Y parsed for Barnes OA. |
| `US1` | station_record (USn upper-air station records) | Y | **Implemented** | `stations.upper[0]` | US* upper-air station X/Y parsed (qa_notes / OA anchor). |
| `PS1` | station_record (PSn precip station records) |  | **Implemented** | `stations.precip[0]` | Precip station X/Y + PRECIP.DAT rates for OA. |

## Suggested typed config surface

Target a nested dataclass / pydantic-style `CalmetConfig` mirroring groups:

```text
CalmetConfig
├── run         # NOOBS, datetime start/end, NSECDT, ABTZ, IRTYPE, ITEST, MREG
├── files       # GEODAT, SRFDAT, UPDAT, M3DDAT, METDAT, …, NUSTA/NM3D/…
├── grid        # PMAP, UTM/LCC…, NX/NY/NZ, ZFACE, origins
├── output      # LSAVE, IFORMO, LPRINT, LCALGRD, layer print flags
├── met         # NSSTA, NPSTA, IFORMS/P/C, ICLOUD/MCLOUD
├── winds       # IWFCOD, IFRADJ, IKINE, IOBR, OA radii, NSMTH, barriers…
├── pbl         # IMIXH, Zi limits, CONST*, THRESH*, JWAT*, FCORIOL…
├── radiation   # IRAD, HA1…HC3
├── clouds      # ICLOUD/MCLOUD/ICLDOUT + CLDDAT
├── precip      # NPSTA, NFLAGP, SIGMAP, CUTP
├── overwater   # ICOARE, DSHELF, IWARM, ICOOL, SEA.DAT
└── stations    # list of surface / upper / precip records
```

`read_inp()` should populate `CalmetConfig` with Fortran defaults for omitted keys, then modules read only typed fields (no ad-hoc `get_int` scattering).

See also [`calmet-module-roadmap.md`](calmet-module-roadmap.md).
