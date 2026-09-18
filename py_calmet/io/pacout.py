"""PACOUT / MESOPAC-II style output hooks (IFORMO=2).

Writes a NumPy ``.npz`` archive with MESOPAC-like fields so downstream tools
can consume an IFORMO=2 path without the legacy binary packing. Not
bit-identical to Fortran PACOUT.DAT.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np


def write_pacout(
    path: str | Path,
    *,
    U_mix: np.ndarray,
    V_mix: np.ndarray,
    U_up: np.ndarray | None = None,
    V_up: np.ndarray | None = None,
    zi: np.ndarray,
    ustar: np.ndarray,
    wstar: np.ndarray,
    el: np.ndarray,
    ipgt: np.ndarray,
    rmm: np.ndarray,
    rho: np.ndarray,
    tempk: np.ndarray,
    qsw: np.ndarray,
    irh: np.ndarray,
    meta: dict[str, Any] | None = None,
) -> None:
    """Persist PACOUT-equivalent fields to ``path`` (``.npz`` or ``.npz.npz``)."""
    path = Path(path)
    if path.suffix.lower() not in (".npz", ".npz.npz"):
        path = path.with_suffix(path.suffix + ".npz") if path.suffix else path.with_suffix(".npz")
    payload = {
        "U_mix": np.asarray(U_mix),
        "V_mix": np.asarray(V_mix),
        "ZI": np.asarray(zi),
        "USTAR": np.asarray(ustar),
        "WSTAR": np.asarray(wstar),
        "EL": np.asarray(el),
        "IPGT": np.asarray(ipgt),
        "RMM": np.asarray(rmm),
        "RHO": np.asarray(rho),
        "TEMPK": np.asarray(tempk),
        "QSW": np.asarray(qsw),
        "IRH": np.asarray(irh),
    }
    if U_up is not None:
        payload["U_up"] = np.asarray(U_up)
    if V_up is not None:
        payload["V_up"] = np.asarray(V_up)
    if meta:
        for k, v in meta.items():
            payload[f"meta_{k}"] = np.asarray(v) if not np.isscalar(v) else v
    np.savez(path, **payload)


def mixed_layer_uv(U: np.ndarray, V: np.ndarray, zi: np.ndarray, zface: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Layer-average U/V within Zi (MESOPAC lower wind).

    ``U,V``: (nt,nz,ny,nx); ``zi``: (nt,ny,nx).
    """
    zmid = 0.5 * (zface[:-1] + zface[1:])
    nt, nz, ny, nx = U.shape
    um = np.zeros((nt, ny, nx), dtype=np.float64)
    vm = np.zeros_like(um)
    for t in range(nt):
        for j in range(ny):
            for i in range(nx):
                mask = zmid <= max(float(zi[t, j, i]), zmid[0])
                if not np.any(mask):
                    mask[0] = True
                um[t, j, i] = float(U[t, mask, j, i].mean())
                vm[t, j, i] = float(V[t, mask, j, i].mean())
    return um, vm
