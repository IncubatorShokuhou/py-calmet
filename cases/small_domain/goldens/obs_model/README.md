# Golden: observation + model mode

## Config
- `NOOBS=0`, `IPROG=14` (3D.DAT as initial guess), `NM3D=1`
- Inputs: `geo.dat`, `surf.dat`, `up.dat`, `3d.dat`
- Same grid/time window as obs mode

## Run
```bash
cd cases/small_domain/obs_model
./calmet.x calmet.inp
```

## Result
- Exit 0, `CALMET.DAT` non-empty (~143 KB)
- Note: `3d.dat` includes 4 hours (00-03) so prognostic interpolation can bracket the 3-hour CALMET window
