"""Wind barriers (NBAR) and lake-breeze (LLBREZE) adjustments.

Barriers: stations on the opposite side of a barrier segment from a grid cell
are blocked for OA (CALMET IJOUT / BARRI). Lake breeze: within NBOX influence
boxes, surface wind is blended toward an onshore vector from METBXID stations.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np


@dataclass
class BarrierSet:
    """Barrier endpoints in km relative to domain origin (CALMET BARXY)."""

    xbbar: np.ndarray  # begin X (km)
    ybbar: np.ndarray
    xebar: np.ndarray  # end X (km)
    yebar: np.ndarray
    kbar: int = 999  # top layer index (1-based); barriers apply for k < kbar

    @classmethod
    def from_inp(cls, inp) -> "BarrierSet | None":
        nbar = int(inp.get_int("NBAR", 0))
        if nbar <= 0:
            return None
        xb = np.atleast_1d(np.asarray(inp.get_list_float("XBBAR") or [], dtype=np.float64))
        yb = np.atleast_1d(np.asarray(inp.get_list_float("YBBAR") or [], dtype=np.float64))
        xe = np.atleast_1d(np.asarray(inp.get_list_float("XEBAR") or [], dtype=np.float64))
        ye = np.atleast_1d(np.asarray(inp.get_list_float("YEBAR") or [], dtype=np.float64))
        # Pad / trim to nbar
        def _pad(a):
            a = np.asarray(a, dtype=np.float64)
            if a.size >= nbar:
                return a[:nbar]
            return np.pad(a, (0, nbar - a.size))

        kbar = int(inp.get_int("KBAR", 999))
        return cls(_pad(xb), _pad(yb), _pad(xe), _pad(ye), kbar=kbar)

    @property
    def nbar(self) -> int:
        return int(self.xbbar.size)


def _barrier_slopes(bar: BarrierSet) -> tuple[np.ndarray, np.ndarray]:
    dx = bar.xbbar - bar.xebar
    dy = bar.ybbar - bar.yebar
    slope = np.where(np.abs(dx) < 1e-12, 9.9e9, dy / np.where(np.abs(dx) < 1e-12, 1.0, dx))
    intercept = bar.ybbar - slope * bar.xbbar
    return slope, intercept


def same_side(
    x: float,
    y: float,
    xs: float,
    ys: float,
    bar: BarrierSet,
) -> bool:
    """True if (x,y) and (xs,ys) are on the same side of every barrier (clear)."""
    if bar.nbar <= 0:
        return True
    slope, intercept = _barrier_slopes(bar)
    for i in range(bar.nbar):
        d1 = ys - slope[i] * xs - intercept[i]
        d2 = y - slope[i] * x - intercept[i]
        if d1 * d2 >= 0.0:
            continue
        # Opposite sides — check whether segment actually blocks (angle test)
        ax = bar.xbbar[i] - x
        ay = bar.ybbar[i] - y
        bx = bar.xebar[i] - x
        by = bar.yebar[i] - y
        cx = xs - x
        cy = ys - y

        def _unidot(u, v):
            nu = np.hypot(u[0], u[1])
            nv = np.hypot(v[0], v[1])
            if nu < 1e-12 or nv < 1e-12:
                return 1.0
            return (u[0] * v[0] + u[1] * v[1]) / (nu * nv)

        cosab = _unidot((ax, ay), (bx, by))
        cosac = _unidot((ax, ay), (cx, cy))
        cosbc = _unidot((bx, by), (cx, cy))
        if cosac >= cosab and cosbc >= cosab:
            return False  # blocked
    return True


def station_clear_mask(
    xc_km: float,
    yc_km: float,
    xs_km: np.ndarray,
    ys_km: np.ndarray,
    bar: BarrierSet | None,
    layer_k: int = 0,
) -> np.ndarray:
    """Boolean mask of stations clear of barriers for grid cell (layer 0-based)."""
    xs = np.atleast_1d(np.asarray(xs_km, dtype=np.float64))
    ys = np.atleast_1d(np.asarray(ys_km, dtype=np.float64))
    if bar is None or bar.nbar <= 0:
        return np.ones(xs.shape, dtype=bool)
    # KBAR is 1-based top level; layers with k >= kbar are clear
    if layer_k + 1 > int(bar.kbar):
        return np.ones(xs.shape, dtype=bool)
    return np.array(
        [same_side(xc_km, yc_km, float(xs[s]), float(ys[s]), bar) for s in range(xs.size)],
        dtype=bool,
    )


@dataclass
class LakeBreezeBox:
    xg1: float
    xg2: float
    yg1: float
    yg2: float
    xbcst: float
    ybcst: float
    xecst: float
    yecst: float
    metbxid: list[int] = field(default_factory=list)


@dataclass
class LakeBreezeConfig:
    boxes: list[LakeBreezeBox]

    @classmethod
    def from_inp(cls, inp) -> "LakeBreezeConfig | None":
        if hasattr(inp, "get_bool"):
            enabled = bool(inp.get_bool("LLBREZE", False))
        else:
            raw = inp.get("LLBREZE", False) if hasattr(inp, "get") else False
            if isinstance(raw, str):
                enabled = raw.strip().upper() in ("T", "TRUE", ".TRUE.", "1", "YES")
            else:
                enabled = bool(raw)
        if not enabled:
            return None
        nbox = int(inp.get_int("NBOX", 0))
        if nbox <= 0:
            return None

        def _list(key, n):
            a = list(inp.get_list_float(key) or [])
            while len(a) < n:
                a.append(0.0)
            return a[:n]

        xg1 = _list("XG1", nbox)
        xg2 = _list("XG2", nbox)
        yg1 = _list("YG1", nbox)
        yg2 = _list("YG2", nbox)
        xbc = _list("XBCST", nbox)
        ybc = _list("YBCST", nbox)
        xec = _list("XECST", nbox)
        yec = _list("YECST", nbox)
        nlb = list(inp.get_list_int("NLB") or [0] * nbox)
        met = list(inp.get_list_int("METBXID") or [])
        boxes = []
        cursor = 0
        for i in range(nbox):
            n = int(nlb[i]) if i < len(nlb) else 0
            ids = met[cursor : cursor + n] if n > 0 else []
            cursor += n
            boxes.append(
                LakeBreezeBox(
                    xg1[i], xg2[i], yg1[i], yg2[i],
                    xbc[i], ybc[i], xec[i], yec[i],
                    metbxid=ids,
                )
            )
        return cls(boxes=boxes)


def _point_in_box(x: float, y: float, box: LakeBreezeBox) -> bool:
    return (min(box.xg1, box.xg2) <= x <= max(box.xg1, box.xg2)) and (
        min(box.yg1, box.yg2) <= y <= max(box.yg1, box.yg2)
    )


def _coast_normal(box: LakeBreezeBox) -> tuple[float, float]:
    """Unit vector roughly onshore (perpendicular to coastline segment)."""
    dx = box.xecst - box.xbcst
    dy = box.yecst - box.ybcst
    # Perpendicular
    nx, ny = -dy, dx
    n = np.hypot(nx, ny)
    if n < 1e-12:
        return 1.0, 0.0
    return float(nx / n), float(ny / n)


def apply_lake_breeze(
    U: np.ndarray,
    V: np.ndarray,
    *,
    xorig_km: float,
    yorig_km: float,
    dgrid_km: float,
    cfg: LakeBreezeConfig,
    stn_u: np.ndarray | None = None,
    stn_v: np.ndarray | None = None,
    stn_ids: Sequence[int] | None = None,
    blend: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Blend surface (layer 0) winds inside lake-breeze boxes toward onshore flow.

    When station UV are provided for METBXID members, use their mean; else use
    coastline-normal * |U,V| magnitude.
    """
    U = np.asarray(U, dtype=np.float64).copy()
    V = np.asarray(V, dtype=np.float64).copy()
    nz, ny, nx = U.shape
    for box in cfg.boxes:
        # Station-mean onshore vector
        ub, vb = 0.0, 0.0
        nuse = 0
        if stn_u is not None and stn_ids is not None and box.metbxid:
            id_list = list(stn_ids)
            for sid in box.metbxid:
                if sid in id_list:
                    k = id_list.index(sid)
                    ub += float(np.atleast_1d(stn_u)[k])
                    vb += float(np.atleast_1d(stn_v)[k])
                    nuse += 1
        nx_hat, ny_hat = _coast_normal(box)
        for j in range(ny):
            for i in range(nx):
                xc = xorig_km + (i + 0.5) * dgrid_km
                yc = yorig_km + (j + 0.5) * dgrid_km
                if not _point_in_box(xc, yc, box):
                    continue
                spd = float(np.hypot(U[0, j, i], V[0, j, i]))
                if nuse > 0:
                    tu, tv = ub / nuse, vb / nuse
                else:
                    # Orient normal toward box interior mid-point
                    mid_x = 0.5 * (box.xg1 + box.xg2)
                    mid_y = 0.5 * (box.yg1 + box.yg2)
                    # Coast midpoint
                    cx = 0.5 * (box.xbcst + box.xecst)
                    cy = 0.5 * (box.ybcst + box.yecst)
                    if (mid_x - cx) * nx_hat + (mid_y - cy) * ny_hat < 0:
                        nx_hat, ny_hat = -nx_hat, -ny_hat
                    tu, tv = nx_hat * max(spd, 1.0), ny_hat * max(spd, 1.0)
                U[0, j, i] = (1.0 - blend) * U[0, j, i] + blend * tu
                V[0, j, i] = (1.0 - blend) * V[0, j, i] + blend * tv
    return U, V
