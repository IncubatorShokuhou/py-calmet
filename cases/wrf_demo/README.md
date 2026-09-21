# wrf_demo — NCAR Katrina tutorial wrfout case

Public WRF sample used to exercise py-calmet beyond the synthetic tiny domain.

## Source wrfout

| Item | Value |
|------|-------|
| File | `wrfout_d01_2005-08-28_00_00_00` |
| Origin | [NCAR/wrf_tutorial_data](https://github.com/NCAR/wrf_tutorial_data) (WRF-Python / NCAR tutorial Katrina case) |
| License | NCAR tutorial sample data (educational / open redistribution via GitHub) |
| Size | ~63 MB (4 times: 00/03/06/09 UTC) |
| Domain | Mercator 30 km, Gulf of Mexico / Mexico coast |

Download:

```bash
git clone --depth 1 https://github.com/NCAR/wrf_tutorial_data.git
cp wrf_tutorial_data/wrfout_d01_2005-08-28_00_00_00 cases/wrf_demo/raw/
```

(An earlier surface-only 4.4 MB Zenodo candidate,
https://zenodo.org/records/18960698 `wrfout_d01_2017-01-12_21_00_00`,
lacks 3-D U/V/T and was not used.)

## What we build

1. Subset a 14×14 mass-point window over mountainous Mexico (`i0=4,j0=0`) with real terrain.
2. Convert → `3d.dat` (CALMET 3D.DAT 2.1 / IOUTMM5≈92) via `scripts/wrfout_to_3d.py`.
3. Temporally interpolate WRF 00Z↔03Z to hourly 00–03 UTC (CALMET requires `NSECDT` ≤ 3600 s).
4. Build matching `geo.dat` (UTM zone 14N), and **synthesize** `surf.dat` / `up.dat` from the
   wrfout center-column near-surface / sounding fields (documented as synthetic, not real obs).
5. Run Fortran `calmet.x` for `obs` / `obs_model` / `noobs`; archive under `goldens/`.

## Rebuild

```bash
python cases/wrf_demo/scripts/prepare_case.py   # or /tmp rebuild hourly helper
# then for each mode:
cd cases/wrf_demo/noobs && ./calmet.x calmet.inp
```

Raw wrfout is gitignored (`cases/wrf_demo/raw/`); goldens + shared text inputs are committed.

## Lat/lon bounding box (for SRTM / higher-res terrain)

Authoritative corners from NCAR Katrina `wrfout` mass-point `XLAT`/`XLONG` on the
committed window (`i0=4, j0=0`, `ni=nj=14`), cross-checked against `3d.dat` header
line `1 1 14 14 … -99.9285 -96.4206 19.1075 22.3873`.

| Domain | South | North | West | East |
|--------|------:|------:|-----:|-----:|
| 3D.DAT / MM5 14×14 (cell centers) | 19.107475 | 22.387329 | -99.928474 | -96.420570 |
| GEO.DAT / CALMET 12×12 interior (cell centers) | 19.362259 | 22.137604 | -99.658630 | -96.690407 |
| SRTM recommend (12×12 ±½ cell) | 19.236107 | 22.263756 | -99.793549 | -96.555489 |
| SRTM recommend (14×14 ±½ cell) | 18.981327 | 22.513477 | -100.063393 | -96.285652 |

Grid spacing ≈ 0.25° / ~28 km. UTM zone **14N** in `geo.dat`. Terrain elev in GEO
matches wrfout `HGT` (max ≈ 2884 m).

