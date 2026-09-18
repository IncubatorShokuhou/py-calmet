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
