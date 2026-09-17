# Golden: no-observation / model-only mode

## Config
- `NOOBS=2`, `IPROG=14`, `NM3D=1`, `NSSTA=0`, `NUSTA=0`
- `ITPROG=2`, `IRHPROG=1`, `ICLOUD=3`, `IEXTRP=1`
- Inputs: `geo.dat`, `3d.dat` only (surf/up unused)

## Run
```bash
cd cases/small_domain/noobs
./calmet.x calmet.inp
```

## Result
- Exit 0, `CALMET.DAT` non-empty (~138 KB)
