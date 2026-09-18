"""First-class wrfout → fields bridge via xarray (+ optional wrf-python).

MM4/MM5.DAT readers stay OutOfScope. Prefer this module (or
``scripts/wrfout_to_3d.py``) to ingest WRF meteorology into 3D.DAT / arrays.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any


def open_wrfout(path: str | Path):
    """Open a wrfout NetCDF with xarray (requires optional ``wrf`` extra)."""
    try:
        import xarray as xr
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "xarray is required for wrfout I/O; pip install 'py-calmet[wrf]'"
        ) from exc
    return xr.open_dataset(Path(path))


def wrfout_summary(path: str | Path) -> dict[str, Any]:
    """Lightweight metadata dict for QA / METLST (no MM5 dependency)."""
    ds = open_wrfout(path)
    try:
        dims = {k: int(v) for k, v in ds.sizes.items()}
        vars_ = sorted(ds.data_vars)
        times = None
        for key in ("Times", "time", "XTIME"):
            if key in ds.variables or key in ds.coords:
                try:
                    times = [str(x) for x in ds[key].values[:3]]
                except Exception:
                    times = None
                break
        return {
            "path": str(path),
            "dims": dims,
            "n_vars": len(vars_),
            "has_U": "U" in ds,
            "has_V": "V" in ds,
            "has_T": "T" in ds,
            "has_PH": "PH" in ds,
            "sample_times": times,
            "engine": "xarray",
        }
    finally:
        ds.close()


def try_wrf_python_destagger(field, stagger: str = "Z"):
    """Optional destagger via wrf-python when installed; else return field."""
    try:
        import wrf  # type: ignore
    except ImportError:
        return field
    try:
        return wrf.destagger(field, stagger=stagger if stagger else None)
    except Exception:
        return field
