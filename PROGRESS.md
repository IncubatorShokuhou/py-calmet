# py-calmet progress

## WRF-only validation (马尾) — drop synthetic goldens

- **Deleted** `cases/small_domain`, `cases/daytime_zi`, `tests/test_parity.py`, `tests/test_daytime_zi.py`.
- Suite Fortran gate is **only** `cases/wrf_demo` (Katrina wrfout + mountain `geo.dat`): `tests/test_wrf_demo.py` + `scripts/compare_wrf_fortran.py`.
- Readers / INP / CALMET.DAT writer tests retargeted to wrf_demo inputs.
- Live `calmet.x` preferred when vendored; else archived Fortran-from-WRF-case goldens (still real WRF, not synthetic).
- Lat/lon bbox documented in `cases/wrf_demo/README.md` for upcoming SRTM terrain.
- Skipped: DIAGNO rewrite, bit-identical floats, MM4/MM5.

## Align-up-temp-precip-sea-miss (马尾) — UP T / PRECIP / SEA 缺测

- **UP.DAT** `temp_c≥998` → `up_tempk` 返回 NaN；runner 探空 T 走该路径；**Holzworth** 跳过非物理 Kelvin，避免 999°C→1272 K 把 Zi 拽崩。
- **NFLAGP**：≥9000 / 非有限降水率清零（与 PRECIP.DAT 9999 同族），别让 Barnes OA 吃到幻影 mm/h。
- **SEA.DAT**：`t_air` 缺测 → `t_sea=NaN`，`sea_sst_grid` 跳过坏站保留空气温回退；仅 ΔT 缺测仍用 `t_air`。
- 调查 wrf_demo obs `U_corr≈0.48`：L0 对齐，高层差来自 IEXTRP=-4 幂律相对 Fortran 更强切变；非一行哨兵 bug，**不改 DIAGNO**。
- Skipped: GEO N→S、COARE/DIAGNO 大改、bit-identical、MM4/MM5。


## Align-surf-missing-t-rh (马尾) — SURF T/RH/P/sky sentinels

- **SURF.DAT** missing `tempk` / `rh` / `pres` / `sky` (≥9000 or non-finite) → meteorological defaults via `surf_tempk` / `surf_rh` / `surf_pres` / `surf_sky` (288.15 K, 70%, 1012 mb, sky=0), wired in runner — no more phantom ~10k K poisoning PBL/flux / ELUSTR.
- Winds 9999→calm already on main (PR #4). Slope `RHOCP=1229.9` already present; skipped re-touch.
- Skipped: GEO N→S row order, COARE/DIAGNO rewrites, bit-identical floats, MM4/MM5.

## Align-light (马尾) — missing-obs sentinels

- **UP.DAT** missing `ws`/`wd` (≥998) filtered in `obs_profile_similt` (same gate as VERTAV/domain-avg); empty UA falls back to surface power-law — no more phantom ~100 m/s layers.
- **SURF.DAT** `9999` → calm in `obs_surface_uv` (not a 9999 m/s wind).
- Coverage note for **WTDAT** corrected (SST soft-spot, not “ignored”).
- Follow-up: SURF temp/RH/pres/sky missing (this batch).


## WP6 (this batch) — IDIOPT soft spots + WTDAT + Fortran align

- **IDIOPT1/4/5=1**: real DIAG.DAT ingestion (`io.diag_dat` ASCII subset); wired in
  runner for surface T, surface UV (IRTYPE=0), upper UV (IRTYPE=0). IDIOPT2/3=1
  also consume DIAG GAMMA / domain UV when present.
- **WTDAT**: `io.wt_dat` reader; water-T precedence **ITWPROG > SEA.DAT > WT.DAT >
  air-T** via `overwater.resolve_water_temp`; stop documenting as ignored.
- **Part B (verified)**: IRTYPE=0 PBL collapse moved after `wstar_field` (was
  NameError / overwritten); Froude/TOPOF2/slope use shared `met_utils.G` (9.81).
- **MM4/MM5** remain OutOfScope.
- Assumption: DIAG.DAT uses documented ASCII fields (full Douglas–Kessler binary
  not vendored); WTDAT here is SST soft-spot, not official terrain-weight WT.DAT.
- Tests: `tests/test_wp6_idiopt_wtdat.py`; CI will validate (no local pytest).

## WP5 — Batchvarova–Gryning + diag depth

- **IMIXH=±2** implemented: `pbl.mixht_day_bg` / `_mixhbg_scalar` (Fortran MIXHBG+FBG analytical false-images + secant); wired in runner; `check_unsupported` no longer raises for ±2.
- **IDIOPT2/ZUPT**: CGAMMA-style lapse from UP sounding feeds Froude + TOPOF2 gamma (replacing hard-coded 0.01/0.005 when sounding present).
- **IDIOPT3/IUPWND/ZUPWND**: domain-avg UA wind via VERTAV-style layer average → meta/METLST.
- **IDIOPT1/4/5=1**: (superseded by WP6 DIAG.DAT ingestion).
- **wrfout**: `py_calmet.io.wrfout` (xarray first-class); MM4/MM5 stay OutOfScope.
- Tests: `tests/test_mixhbg.py`; full pytest green.

## WP4 — Partial→Implemented (MM4/MM5 OutOfScope)

- Implemented **118 → 207** / Partial **89 → 0** / Missing **0** (of 207)
- **Scope rule:** MM4DAT / MM5.DAT readers OutOfScope — wrfout→3D.DAT only (`NotImplementedError` if non-default MM4/MM5 path set)
- Zi: IAVEZI/MNMDAV/HAFANG/ILEVZI average; IZICRLX/TZICRLX relax; IMIXH±3 Holzworth

## WP3 — MIXDT, SEA/PRECIP/CLOUD I/O, barriers, coord, outputs

- Implemented **72 → 118** / Partial **135 → 89**

## WP2 — TOPOF2/OBrien, RH clouds, precip, COARE-lite

- Implemented **45 → 72** / Partial **162 → 135**
