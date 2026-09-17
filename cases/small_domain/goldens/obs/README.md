# Golden: observation-dominant mode

## Config
- `NOOBS=0`, `IPROG=0`, `NM3D=0`
- Inputs: `geo.dat`, `surf.dat`, `up.dat` (no 3D.DAT)
- Domain: 12x12 @ 1 km, UTM 19N, origin (400, 4900) km
- Time: 2020-06-15 00:00-03:00 UTC (julian 167), `ABTZ=UTC+0000`

## Run
```bash
cd cases/small_domain/obs
./calmet.x calmet.inp
```

## Result
- Exit 0, `CALMET.DAT` non-empty (~143 KB)
- List file: `CALMET.LST` ends with "End of run"
