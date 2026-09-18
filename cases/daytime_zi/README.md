# daytime_zi — convective mixing-height golden

Synthetic 12×12 @ 1 km domain (same geo/3D meteorology as `small_domain`) run
**00–17 UTC on 2020-06-15** so CALMET’s Maul–Carson ZI can accumulate from
before 05:00 LST through peak afternoon insolation.

## Why this case

The NCAR Katrina `wrf_demo` archive is night/early-morning only (SWDOWN=0).
This case proves the **daytime Carson / MIXHMC path** against a Fortran
`CALMET.DAT` golden (`goldens/noobs/`).

## Rebuild

```bash
# inputs already under goldens/; to regenerate Fortran golden:
cd cases/daytime_zi/work   # or recreate from scripts
./calmet.x calmet.inp
cp calmet.dat ../goldens/noobs/CALMET.DAT
```

Requires `vendor/calmet-fortran/src/calmet.x` (not shipped in the public tree).
