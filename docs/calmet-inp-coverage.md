# CALMET.INP parameter coverage (py-calmet)

**Generated:** 2026-09-18 01:07 UTC (Asia/Shanghai)

Inventory of every control-file variable recognized by Fortran `READCF` / `READFN` (and station free-form records), unioned with keys present in `cases/**/calmet.inp`, scored against current `py_calmet` behavior.

## Summary counts

| Status | Count |
|--------|------:|
| **Total** | **207** |
| Implemented | 43 |
| Partial | 164 |
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
| `METINP` | character |  | **Partial** | `files.metinp` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `GEODAT` | character | Y | **Implemented** | `files.geodat` | Honored via _resolve_data_file (case_dir / inputs_dir). |
| `SRFDAT` | character | Y | **Implemented** | `files.srfdat` | Honored via _resolve_data_file. |
| `PRCDAT` | character |  | **Partial** | `files.prcdat` | Accepted on CalmetConfig; physics TBD. (PRECIP.DAT unused.) |
| `MM4DAT` | character |  | **Partial** | `files.mm4dat` | Accepted on CalmetConfig; physics TBD. (Legacy MM4 name; use M3DDAT.) |
| `WTDAT` | character |  | **Partial** | `files.wtdat` | Accepted on CalmetConfig; physics TBD. (WT.DAT unused.) |
| `METLST` | character | Y | **Partial** | `files.metlst` | Parsed + exposed; list-file writer not yet implemented. |
| `METDAT` | character | Y | **Partial** | `files.metdat` | Parsed + exposed; output path still selected by caller/API (not auto-written to METDAT). |
| `PACDAT` | character |  | **Partial** | `files.pacdat` | Accepted on CalmetConfig; physics TBD. (MESOPUFF PACOUT unused.) |
| `CLDDAT` | character |  | **Partial** | `files.clddat` | Accepted on CalmetConfig; physics TBD. (CLOUD.DAT reader unused.) |
| `LCFILES` | logical | Y | **Partial** | `files.lcfiles` | Parsed; paths used as-is (case folding unused). |
| `NUSTA` | integer | Y | **Partial** | `files.nusta` | Parsed; multi-station UP still single-file path. |
| `NOWSTA` | integer | Y | **Partial** | `files.nowsta` | Parsed on config; SEA.DAT physics TBD. |
| `NM3D` | integer | Y | **Partial** | `files.nm3d` | Parsed; single 3d.dat path honored. |
| `NIGF` | integer | Y | **Partial** | `files.nigf` | Parsed on config; IGF reader TBD. |

## Group 0b — Upper-air file names (READFN b)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `UPDAT` | character (per upper-air station) | Y | **Implemented** | `files.updat` | Honored via _resolve_data_file. |

## Group 0c — Overwater / SEA.DAT file names (READFN c)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `SEADAT` | character (per overwater station) |  | **Partial** | `files.seadat` | Accepted on CalmetConfig; physics TBD. (SEA.DAT unused.) |

## Group 0d — MM4/MM5/3D.DAT file names (READFN d)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `M3DDAT` | character (per MM5/3D file) | Y | **Implemented** | `files.m3ddat` | Honored via _resolve_data_file (3d.dat). |

## Group 0e — IGF-CALMET.DAT file names (READFN e)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `IGFDAT` | character (per IGF file) |  | **Partial** | `files.igfdat` | Accepted on CalmetConfig; physics TBD. (IGF-CALMET files unused.) |

## Group 0f — Misc diagnostic / test file names (READFN f)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `DIADAT` | character |  | **Partial** | `files.diadat` | Accepted on CalmetConfig; physics TBD. (DIAG.DAT unused.) |
| `PRGDAT` | character |  | **Partial** | `files.prgdat` | Accepted on CalmetConfig; physics TBD. (PROG.DAT unused.) |
| `TSTPRT` | character |  | **Partial** | `files.tstprt` | Accepted on CalmetConfig; physics TBD. (Test files unused.) |
| `TSTOUT` | character |  | **Partial** | `files.tstout` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `TSTKIN` | character |  | **Partial** | `files.tstkin` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `TSTFRD` | character |  | **Partial** | `files.tstfrd` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `TSTSLP` | character |  | **Partial** | `files.tstslp` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `DCSTGD` | character |  | **Partial** | `files.dcstgd` | Accepted on CalmetConfig; physics TBD. (Unused.) |

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
| `IBTZ` | integer |  | **Partial** | `run.ibtz_legacy` | Accepted on CalmetConfig; physics TBD. (Legacy; prefer ABTZ.) |
| `IRLG` | integer |  | **Partial** | `run.irlg_legacy` | Accepted on CalmetConfig; physics TBD. (Legacy run length; prefer begin/end + NSECDT.) |
| `NSECDT` | integer | Y | **Implemented** | `run.nsecdt` | Timestep seconds (datetime timedelta). |
| `IRTYPE` | integer | Y | **Partial** | `run.irtype` | Accepted on CalmetConfig; physics TBD. (Run-type gating unused.) |
| `LCALGRD` | logical | Y | **Partial** | `output.lcalgrd` | Some CALGRID fields in writer; flag not read. |
| `ITEST` | integer | Y | **Partial** | `run.itest` | Accepted on CalmetConfig; physics TBD. (Setup-only stop unused.) |
| `MREG` | integer | Y | **Partial** | `run.mreg` | Accepted on CalmetConfig; physics TBD. (Regulatory QA unused.) |

## Group 2 — Map projection and grid

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `PMAP` | character | Y | **Partial** | `grid.pmap` | UTM assumed; LCC/PS/EM/TTM/LAZA not wired. |
| `DATUM` | character | Y | **Partial** | `grid.datum` | Not used in coord transforms. |
| `FEAST` | real | Y | **Partial** | `grid.feast` | Accepted on CalmetConfig; physics TBD. (Non-UTM false easting unused.) |
| `FNORTH` | real | Y | **Partial** | `grid.fnorth` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IUTMZN` | integer | Y | **Implemented** | `grid.iutmzn` | UTM zone for lat/lon / header. |
| `UTMHEM` | character | Y | **Implemented** | `grid.utmhem` | UTM hemisphere. |
| `RLAT0` | character |  | **Implemented** | `grid.rlat0` | Optional domain lat for solar/Coriolis. |
| `RLON0` | character |  | **Implemented** | `grid.rlon0` | Optional domain lon for solar. |
| `XLAT1` | character |  | **Partial** | `grid.xlat1` | Accepted on CalmetConfig; physics TBD. (LCC/PS parallel unused.) |
| `XLAT2` | character |  | **Partial** | `grid.xlat2` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `NX` | integer | Y | **Partial** | `grid.nx` | Taken from GEO.DAT, not INP. |
| `NY` | integer | Y | **Partial** | `grid.ny` | Taken from GEO.DAT. |
| `DGRIDKM` | real | Y | **Partial** | `grid.dgridkm` | Taken from GEO.DAT. |
| `XORIGKM` | real | Y | **Partial** | `grid.xorigkm` | Taken from GEO.DAT. |
| `YORIGKM` | real | Y | **Partial** | `grid.yorigkm` | Taken from GEO.DAT. |
| `NZ` | integer | Y | **Implemented** | `grid.nz` | Number of layers. |
| `ZFACE` | real (length NZ+1) | Y | **Implemented** | `grid.zface` | Vertical grid faces. |

## Group 3 — Output options

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `LSAVE` | logical | Y | **Partial** | `output.lsave` | Caller decides write; flag not read. |
| `LPRINT` | logical | Y | **Partial** | `output.lprint` | Accepted on CalmetConfig; physics TBD. (Printer output unused.) |
| `IPRINF` | integer | Y | **Partial** | `output.iprinf` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IUVOUT` | integer (length NZ) | Y | **Partial** | `output.iuvout` | Accepted on CalmetConfig; physics TBD. (Layer UV print unused.) |
| `IWOUT` | integer (length NZ) | Y | **Partial** | `output.iwout` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `ITOUT` | integer (length NZ) | Y | **Partial** | `output.itout` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `STABILITY` | logical | Y | **Partial** | `output.stability` | IPGT computed; print flag unused. |
| `USTAR` | logical | Y | **Partial** | `output.ustar` | Field computed; print flag unused. |
| `MONIN` | logical | Y | **Partial** | `output.monin` | EL computed; print flag unused. |
| `MIXHT` | logical | Y | **Partial** | `output.mixht` | ZI computed; print flag unused. |
| `WSTAR` | logical | Y | **Partial** | `output.wstar` | Field computed; print flag unused. |
| `PRECIP` | logical | Y | **Partial** | `output.precip` | RMM zeroed; print/data path unused. |
| `SENSHEAT` | logical | Y | **Partial** | `output.sensheat` | QH internal; output partial. |
| `CONVZI` | logical | Y | **Partial** | `output.convzi` | ziconv for Carson; print flag unused. |
| `LDB` | logical | Y | **Partial** | `output.ldb` | Accepted on CalmetConfig; physics TBD. (Debug unused.) |
| `NN1` | integer | Y | **Partial** | `output.nn1` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `NN2` | integer | Y | **Partial** | `output.nn2` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `LDBCST` | logical | Y | **Partial** | `output.ldbcst` | Accepted on CalmetConfig; physics TBD. (Coast distance output unused.) |
| `IOUTD` | integer | Y | **Partial** | `output.ioutd` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `NZPRN2` | integer | Y | **Partial** | `output.nzprn2` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IPR0` | integer | Y | **Partial** | `output.ipr0` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IPR1` | integer | Y | **Partial** | `output.ipr1` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IPR2` | integer | Y | **Partial** | `output.ipr2` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IPR3` | integer | Y | **Partial** | `output.ipr3` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IPR4` | integer | Y | **Partial** | `output.ipr4` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IPR5` | integer | Y | **Partial** | `output.ipr5` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IPR6` | integer | Y | **Partial** | `output.ipr6` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IPR7` | integer | Y | **Partial** | `output.ipr7` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IPR8` | integer | Y | **Partial** | `output.ipr8` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IFORMO` | integer | Y | **Partial** | `output.iformo` | CALMET.DAT only (type 1); PACOUT missing. |

## Group 4 — Meteorological data options

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `NOOBS` | integer | Y | **Implemented** | `run.noobs` | Mode inference (noobs / obs / obs_model). |
| `NSSTA` | integer | Y | **Partial** | `met.nssta` | Assumes 0/1 station; multi-station unused. |
| `NPSTA` | integer | Y | **Partial** | `met.npsta` | Accepted on CalmetConfig; physics TBD. (Precipitation stations unused.) |
| `IFORMS` | integer | Y | **Partial** | `met.iforms` | SURF format assumed. |
| `IFORMP` | integer | Y | **Partial** | `met.iformp` | Accepted on CalmetConfig; physics TBD. (Precip file format unused.) |
| `ICLOUD` | integer | Y | **Partial** | `clouds.icloud` | Cloud from SURF sky tenths; schemes 1–4 missing. |
| `ICLDOUT` | integer |  | **Partial** | `clouds.icldout` | Accepted on CalmetConfig; physics TBD. (Cloud output flag unused.) |
| `MCLOUD` | integer |  | **Partial** | `clouds.mcloud` | Accepted on CalmetConfig; physics TBD. (INP 2.2 cloud method unused.) |
| `IFORMC` | integer | Y | **Partial** | `met.iformc` | Accepted on CalmetConfig; physics TBD. (Cloud file format unused.) |

## Group 5 — Wind field options and parameters

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `IWFCOD` | integer | Y | **Partial** | `winds.iwfcod` | Diagnostic winds always on; flag not read. |
| `IFRADJ` | integer | Y | **Implemented** | `winds.ifradj` | Froude blocking on/off (FRADJ). |
| `IKINE` | integer | Y | **Partial** | `winds.ikine` | Flag gates stub kinematic+div path; full TOPOF2 W missing. |
| `IOBR` | integer | Y | **Partial** | `winds.iobr` | Flag gates simplified div-min; full O'Brien missing. |
| `IEXTRP` | integer | Y | **Partial** | `winds.iextrp` | Power-law obs profile; SIMILT helper exists but not wired. |
| `RMIN2` | real | Y | **Partial** | `winds.rmin2` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `FEXTR2` | real (length NZ) | Y | **Partial** | `winds.fextr2` | Accepted on CalmetConfig; physics TBD. (Extrapolation factors unused.) |
| `IPROG` | integer | Y | **Implemented** | `winds.iprog` | Mode inference when >0 with observations. |
| `ISTEPPG` | integer |  | **Partial** | `winds.isteppg` | Accepted on CalmetConfig; physics TBD. (Legacy hours; prefer ISTEPPGS.) |
| `ISTEPPGS` | integer | Y | **Partial** | `winds.isteppgs` | NSECDT drives time; prognostic-step QA missing. |
| `IGFMET` | integer | Y | **Partial** | `winds.igfmet` | Accepted on CalmetConfig; physics TBD. (IGF as first-guess unused.) |
| `LVARY` | logical | Y | **Partial** | `winds.lvary` | Accepted on CalmetConfig; physics TBD. (Varying radius unused.) |
| `RMAX1` | real | Y | **Partial** | `winds.rmax1` | Accepted on CalmetConfig; physics TBD. (OA distance cutoff unused.) |
| `RMAX2` | real | Y | **Partial** | `winds.rmax2` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `RMAX3` | real | Y | **Partial** | `winds.rmax3` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `RMIN` | real | Y | **Partial** | `winds.rmin` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `TERRAD` | real | Y | **Implemented** | `winds.terrad` | Terrain radius (km) for FRADJ/slope. |
| `R1` | real | Y | **Partial** | `winds.r1` | Single-station Barnes surface radius; multi-station OA missing. |
| `R2` | real | Y | **Partial** | `winds.r2` | Aloft OA radius; multi-station incomplete. |
| `RPROG` | real | Y | **Partial** | `winds.rprog` | Accepted on CalmetConfig; physics TBD. (Prog weight in OA unused.) |
| `DIVLIM` | real | Y | **Partial** | `winds.divlim` | Accepted on CalmetConfig; physics TBD. (Divergence criterion unused (hardcoded loop).) |
| `NITER` | integer | Y | **Partial** | `winds.niter` | Used (capped) in divergence minimization. |
| `NSMTH` | integer (length NZ) | Y | **Implemented** | `winds.nsmth` | Per-layer 5-point smooth passes. |
| `NINTR2` | integer (length NZ) | Y | **Partial** | `winds.nintr2` | Accepted on CalmetConfig; physics TBD. (Max stations per layer unused.) |
| `CRITFN` | real | Y | **Implemented** | `winds.critfn` | Critical Froude number. |
| `ALPHA` | real | Y | **Partial** | `winds.alpha` | Used in light kinematic/div-min; full IKINE semantics incomplete. |
| `NBAR` | integer | Y | **Partial** | `winds.nbar` | Accepted on CalmetConfig; physics TBD. (Wind barriers unused.) |
| `XBBAR` | real (length NBAR) | Y | **Partial** | `winds.xbbar` | Accepted on CalmetConfig; physics TBD. (Barrier coords unused.) |
| `YBBAR` | real (length NBAR) | Y | **Partial** | `winds.ybbar` | Accepted on CalmetConfig; physics TBD. (Barrier coords unused.) |
| `XEBAR` | real (length NBAR) | Y | **Partial** | `winds.xebar` | Accepted on CalmetConfig; physics TBD. (Barrier coords unused.) |
| `YEBAR` | real (length NBAR) | Y | **Partial** | `winds.yebar` | Accepted on CalmetConfig; physics TBD. (Barrier coords unused.) |
| `KBAR` | integer | Y | **Partial** | `winds.kbar` | Accepted on CalmetConfig; physics TBD. (Barrier top level unused.) |
| `IDIOPT1` | integer | Y | **Partial** | `winds.idiopt1` | Accepted on CalmetConfig; physics TBD. (Diag option switches unused.) |
| `IDIOPT2` | integer | Y | **Partial** | `winds.idiopt2` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IDIOPT3` | integer | Y | **Partial** | `winds.idiopt3` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IDIOPT4` | integer | Y | **Partial** | `winds.idiopt4` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IDIOPT5` | integer | Y | **Partial** | `winds.idiopt5` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `ISURFT` | integer | Y | **Partial** | `winds.isurft` | Accepted on CalmetConfig; physics TBD. (Surface T station index unused.) |
| `IUPT` | integer | Y | **Partial** | `winds.iupt` | Accepted on CalmetConfig; physics TBD. (Upper T index unused.) |
| `ZUPT` | real | Y | **Partial** | `winds.zupt` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `IUPWND` | integer | Y | **Partial** | `winds.iupwnd` | Accepted on CalmetConfig; physics TBD. (Upper wind index unused.) |
| `ZUPWND` | real (length 2) | Y | **Partial** | `winds.zupwnd` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `LLBREZE` | logical | Y | **Partial** | `winds.llbreze` | Accepted on CalmetConfig; physics TBD. (Lake-breeze unused.) |
| `NBOX` | integer | Y | **Partial** | `winds.nbox` | Accepted on CalmetConfig; physics TBD. (Lake-breeze boxes unused.) |
| `XG1` | real (length NBOX) | Y | **Partial** | `winds.xg1` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `XG2` | real (length NBOX) | Y | **Partial** | `winds.xg2` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `YG1` | real (length NBOX) | Y | **Partial** | `winds.yg1` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `YG2` | real (length NBOX) | Y | **Partial** | `winds.yg2` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `XBCST` | real (length NBOX) | Y | **Partial** | `winds.xbcst` | Accepted on CalmetConfig; physics TBD. (Coastline box unused.) |
| `YBCST` | real (length NBOX) | Y | **Partial** | `winds.ybcst` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `XECST` | real (length NBOX) | Y | **Partial** | `winds.xecst` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `YECST` | real (length NBOX) | Y | **Partial** | `winds.yecst` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `NLB` | integer | Y | **Partial** | `winds.nlb` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `METBXID` | integer (length mxbxwnd) | Y | **Partial** | `winds.metbxid` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `BIAS` | real (length NZ) | Y | **Partial** | `winds.bias` | Not applied to layer blending. |
| `ISLOPE` | integer | Y | **Implemented** | `winds.islope` | Slope flow on/off (Mahrt). |
| `ICALM` | integer | Y | **Partial** | `winds.icalm` | Accepted on CalmetConfig; physics TBD. (Calm processing unused.) |

## Group 6 — Mixing height, temperature, precip, overwater

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `CONSTB` | real | Y | **Implemented** | `pbl.constb` | Carson buoyancy constant. |
| `CONSTE` | real | Y | **Partial** | `pbl.conste` | Energy-budget path partial; CONSTE unused. |
| `CONSTN` | real | Y | **Implemented** | `pbl.constn` | Night Zi constant. |
| `DPTMIN` | real | Y | **Partial** | `pbl.dptmin` | Constant Carson gamma; MIXDT sounding lapse deferred. |
| `DZZI` | real | Y | **Partial** | `pbl.dzzi` | Accepted on CalmetConfig; physics TBD. (Inversions thickness unused.) |
| `ZIMIN` | real | Y | **Implemented** | `pbl.zimin` | Min mixing height. |
| `ZIMAX` | real | Y | **Implemented** | `pbl.zimax` | Max mixing height. |
| `ZIMINW` | real | Y | **Partial** | `pbl.ziminw` | Accepted on CalmetConfig; physics TBD. (Overwater Zi min unused.) |
| `ZIMAXW` | real | Y | **Partial** | `pbl.zimaxw` | Accepted on CalmetConfig; physics TBD. (Overwater Zi max unused.) |
| `IAVEZI` | integer | Y | **Partial** | `pbl.iavezi` | Accepted on CalmetConfig; physics TBD. (Zi averaging unused.) |
| `MNMDAV` | integer | Y | **Partial** | `pbl.mnmdav` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `HAFANG` | real | Y | **Partial** | `pbl.hafang` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `ILEVZI` | integer | Y | **Partial** | `pbl.ilevzi` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `FCORIOL` | real | Y | **Partial** | `pbl.fcoriol` | Coriolis from lat; FCORIOL INP ignored. |
| `CONSTW` | real | Y | **Partial** | `pbl.constw` | Accepted on CalmetConfig; physics TBD. (Overwater mixing unused.) |
| `ITPROG` | integer | Y | **Partial** | `temp.itprog` | 3D T2 used in noobs; full ITPROG paths missing. |
| `ITWPROG` | integer | Y | **Partial** | `pbl.itwprog` | Accepted on CalmetConfig; physics TBD. (Water T from prog unused.) |
| `ILUOC3D` | integer | Y | **Partial** | `pbl.iluoc3d` | Accepted on CalmetConfig; physics TBD. (3D landuse overlay unused.) |
| `IRAD` | integer | Y | **Partial** | `radiation.irad` | Solar always computed; IRAD option not read. |
| `IAVET` | integer | Y | **Partial** | `pbl.iavet` | Accepted on CalmetConfig; physics TBD. (T averaging unused.) |
| `TGDEFB` | real | Y | **Partial** | `pbl.tgdefb` | Accepted on CalmetConfig; physics TBD. (Default lapse (below) unused.) |
| `TGDEFA` | real | Y | **Partial** | `pbl.tgdefa` | Accepted on CalmetConfig; physics TBD. (Default lapse (above) unused.) |
| `JWAT1` | integer (length mxwb) | Y | **Implemented** | `pbl.jwat1` | INP JWAT1/JWAT2 bound; effective_iwat() aliases IWAT for heatfx/slope (999→GEO 55). |
| `JWAT2` | integer (length mxwb) | Y | **Implemented** | `pbl.jwat2` | Paired with JWAT1; see JWAT1. |
| `TRADKM` | real | Y | **Partial** | `pbl.tradkm` | Accepted on CalmetConfig; physics TBD. (Temperature OA radius unused.) |
| `NUMTS` | integer | Y | **Partial** | `pbl.numts` | Accepted on CalmetConfig; physics TBD. (Unused.) |
| `NFLAGP` | integer | Y | **Partial** | `pbl.nflagp` | Accepted on CalmetConfig; physics TBD. (Precip QC unused.) |
| `SIGMAP` | real | Y | **Partial** | `pbl.sigmap` | Accepted on CalmetConfig; physics TBD. (Precip sigma unused.) |
| `CUTP` | real | Y | **Partial** | `pbl.cutp` | Accepted on CalmetConfig; physics TBD. (Precip cutoff unused.) |
| `HA1` | real |  | **Implemented** | `radiation.ha1` | Bound into shortwave_radiation from CalmetConfig. |
| `HA2` | real |  | **Implemented** | `radiation.ha2` | Bound into shortwave_radiation from CalmetConfig. |
| `HB1` | real |  | **Implemented** | `radiation.hb1` | Bound into shortwave_radiation from CalmetConfig. |
| `HB2` | real |  | **Implemented** | `radiation.hb2` | Bound into shortwave_radiation from CalmetConfig. |
| `HC1` | real |  | **Implemented** | `radiation.hc1` | Bound into heat_flux_energy_budget from CalmetConfig. |
| `HC2` | real |  | **Implemented** | `radiation.hc2` | Bound into heat_flux_energy_budget from CalmetConfig. |
| `HC3` | real |  | **Implemented** | `radiation.hc3` | Bound into heat_flux_energy_budget from CalmetConfig. |
| `IMIXH` | integer | Y | **Partial** | `pbl.imixh` | Maul–Carson day + night mechanical only; other IMIXH options missing. |
| `THRESHL` | real | Y | **Implemented** | `pbl.threshl` | Daytime MIXHMC growth threshold. |
| `THRESHW` | real | Y | **Partial** | `pbl.threshw` | Accepted on CalmetConfig; physics TBD. (Overwater threshold unused.) |
| `ICOARE` | integer | Y | **Partial** | `overwater.icoare` | Accepted on CalmetConfig; physics TBD. (No COARE overwater flux module (PSIUD helper only).) |
| `DSHELF` | real | Y | **Partial** | `overwater.dshelf` | Accepted on CalmetConfig; physics TBD. (COARE coastal shelf unused.) |
| `IWARM` | integer | Y | **Partial** | `overwater.iwarm` | Accepted on CalmetConfig; physics TBD. (COARE warm-layer unused.) |
| `ICOOL` | integer | Y | **Partial** | `overwater.icool` | Accepted on CalmetConfig; physics TBD. (COARE cool-skin unused.) |
| `IRHPROG` | integer | Y | **Partial** | `humidity.irhprog` | RH from 3D/SURF; flag not read. |
| `IZICRLX` | integer |  | **Partial** | `pbl.izicrlx` | Accepted on CalmetConfig; physics TBD. (Zi relaxation unused.) |
| `TZICRLX` | real |  | **Partial** | `pbl.tzicrlx` | Accepted on CalmetConfig; physics TBD. (Unused.) |

## Group station — Station location free-format records (IG 7–9)

| Variable | Type | Sample? | Status | Suggested `py_calmet` API | Notes |
|----------|------|---------|--------|---------------------------|-------|
| `SS1` | station_record (SSn surface station records) | Y | **Partial** | `stations.surface[0]` | Single-station X/Y for OA; multi-station SSn unused. |
| `US1` | station_record (USn upper-air station records) | Y | **Partial** | `stations.upper[0]` | UP.DAT sounding used; US1 coords unused. |
| `PS1` | station_record (PSn precip station records) |  | **Partial** | `stations.precip[0]` | Accepted on CalmetConfig; physics TBD. (Precip station records unused.) |

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
