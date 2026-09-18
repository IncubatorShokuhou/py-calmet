#!/bin/python3

from __future__ import annotations

import argparse
import struct
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import BinaryIO, Dict, List, Optional, Tuple

import numpy as np

from .fortran_bin import (
    FortranRecordError,
    FortranSequentialReader,
    decode_string,
    parse_payload,
    split_labeled_record,
)


@dataclass
class RunControl:
    """Parsed CALMET run-control record."""

    start_time: datetime
    end_time: datetime
    timezone: str
    irlg: int
    irtype: int
    nx: int
    ny: int
    nz: int
    dgrid: float
    xorigr: float
    yorigr: float
    iwfcod: int
    nssta: int
    nusta: int
    npsta: int
    nowsta: int
    nlu: int
    iwat1: int
    iwat2: int
    lcalgrd: bool
    pmap: str
    datum: str
    daten: str
    feast: float
    fnorth: float
    utmhem: str
    iutmzn: int
    rnlat0: float
    relon0: float
    xlat1: float
    xlat2: float


@dataclass
class CalmetDataset:
    """Reader for `CALMET.DAT` files."""

    path: Path
    dataset_name: str = ""
    dataset_version: str = ""
    data_model: str = ""
    comments: List[str] = field(default_factory=list)
    run_control: Optional[RunControl] = None
    static_fields: Dict[str, np.ndarray] = field(default_factory=dict)
    fields_2d: Dict[str, np.ndarray] = field(default_factory=dict)
    fields_3d: Dict[str, np.ndarray] = field(default_factory=dict)
    time_bounds: List[Tuple[datetime, datetime]] = field(default_factory=list)

    @property
    def nx(self) -> int:
        """Return the grid size in the x direction."""
        return 0 if self.run_control is None else self.run_control.nx

    @property
    def ny(self) -> int:
        """Return the grid size in the y direction."""
        return 0 if self.run_control is None else self.run_control.ny

    @property
    def nz(self) -> int:
        """Return the number of vertical levels."""
        return 0 if self.run_control is None else self.run_control.nz

    @property
    def nt(self) -> int:
        """Return the number of loaded timesteps."""
        return len(self.time_bounds)

    @property
    def dx(self) -> float:
        """Return the x cell size in map units."""
        return 0.0 if self.run_control is None else self.run_control.dgrid

    @property
    def dy(self) -> float:
        """Return the y cell size in map units."""
        return 0.0 if self.run_control is None else self.run_control.dgrid

    @property
    def x(self) -> np.ndarray:
        """Return x coordinates for cell centers."""
        if self.run_control is None:
            raise ValueError("Dataset is not initialized")
        # XORIGR is the SW corner of cell (1,1); cell centers are +0.5 cell.
        return self.run_control.xorigr + (np.arange(self.nx, dtype=np.float32) + 0.5) * self.run_control.dgrid

    @property
    def y(self) -> np.ndarray:
        """Return y coordinates for cell centers."""
        if self.run_control is None:
            raise ValueError("Dataset is not initialized")
        return self.run_control.yorigr + (np.arange(self.ny, dtype=np.float32) + 0.5) * self.run_control.dgrid

    @classmethod
    def read(cls, path: str | Path) -> "CalmetDataset":
        """Read a `CALMET.DAT` file and return a populated dataset."""
        dataset = cls(path=Path(path))

        with FortranSequentialReader(dataset.path) as reader:
            dataset._read_header(reader)
            dataset._read_static_fields(reader)
            dataset._read_time_steps(reader)

        return dataset

    def _read_header(self, reader: FortranSequentialReader) -> None:
        record = reader.read_record()
        self.dataset_name = decode_string(record[0:16])
        self.dataset_version = decode_string(record[16:32])
        self.data_model = decode_string(record[32:96])

        ncom_record = reader.read_record()
        ncom = struct.unpack(f"{reader.endian}i", ncom_record)[0]

        self.comments = [decode_string(reader.read_record()) for _ in range(ncom)]

        run_record = reader.read_record()
        self.run_control = _parse_run_control(run_record, reader.endian)

    def _read_static_fields(self, reader: FortranSequentialReader) -> None:
        if self.run_control is None:
            raise ValueError("Run control must be loaded before static fields")

        specs: List[Tuple[str, str, Tuple[int, ...]]] = [
            ("ZFACE", "real", (self.nz + 1,)),
        ]

        if self.run_control.nssta >= 1:
            specs.extend(
                [
                    ("XSSTA", "real", (self.run_control.nssta,)),
                    ("YSSTA", "real", (self.run_control.nssta,)),
                ]
            )
        if self.run_control.nusta >= 1:
            specs.extend(
                [
                    ("XUSTA", "real", (self.run_control.nusta,)),
                    ("YUSTA", "real", (self.run_control.nusta,)),
                ]
            )
        if self.run_control.npsta >= 1:
            specs.extend(
                [
                    ("XPSTA", "real", (self.run_control.npsta,)),
                    ("YPSTA", "real", (self.run_control.npsta,)),
                ]
            )

        specs.extend(
            [
                ("Z0", "real", (self.ny, self.nx)),
                ("ILANDU", "int", (self.ny, self.nx)),
                ("ELEV", "real", (self.ny, self.nx)),
                ("XLAI", "real", (self.ny, self.nx)),
            ]
        )

        if self.run_control.nssta >= 1:
            specs.append(("NEARS", "int", (self.ny, self.nx)))

        for expected_label, dtype_name, shape in specs:
            label, _, payload = split_labeled_record(reader.read_record(), reader.endian)
            if label != expected_label:
                raise FortranRecordError(
                    f"Expected static record {expected_label}, found {label}"
                )
            self.static_fields[label] = parse_payload(payload, dtype_name, shape, reader.endian)

    def _read_time_steps(self, reader: FortranSequentialReader) -> None:
        if self.run_control is None:
            raise ValueError("Run control must be loaded before time steps")

        step_specs = self._time_step_specs()
        layered_targets = {"U", "V", "W", "T"}
        field2d_store: Dict[str, List[np.ndarray]] = {}
        field3d_store: Dict[str, List[np.ndarray]] = {}

        while True:
            try:
                first_record = reader.read_record()
            except EOFError:
                break

            step_arrays_3d: Dict[str, List[np.ndarray]] = {}
            step_arrays_2d: Dict[str, np.ndarray] = {}
            step_start: Optional[datetime] = None
            step_end: Optional[datetime] = None

            all_records = [first_record]
            all_records.extend(reader.read_record() for _ in range(len(step_specs) - 1))

            for raw_record, (expected_label, dtype_name, shape, target_name) in zip(
                all_records, step_specs
            ):
                label, time_info, payload = split_labeled_record(raw_record, reader.endian)
                if label != expected_label:
                    raise FortranRecordError(
                        f"Expected time-step record {expected_label}, found {label}"
                    )

                begin_time = _parse_yyyyjjjhh(time_info[0], time_info[1])
                end_time = _parse_yyyyjjjhh(time_info[2], time_info[3])
                if step_start is None:
                    step_start = begin_time
                    step_end = end_time

                array = parse_payload(payload, dtype_name, shape, reader.endian)
                if target_name in layered_targets:
                    step_arrays_3d.setdefault(target_name, []).append(array)
                else:
                    step_arrays_2d[target_name] = array

            assert step_start is not None and step_end is not None
            self.time_bounds.append((step_start, step_end))

            for name, layers in step_arrays_3d.items():
                field3d_store.setdefault(name, []).append(np.stack(layers, axis=0))
            for name, array in step_arrays_2d.items():
                field2d_store.setdefault(name, []).append(array)

        self.fields_3d = {
            name: np.stack(step_values, axis=0).astype(np.float32, copy=False)
            for name, step_values in field3d_store.items()
        }
        self.fields_2d = {
            name: np.stack(step_values, axis=0) for name, step_values in field2d_store.items()
        }

    def _time_step_specs(self) -> List[Tuple[str, str, Tuple[int, ...], str]]:
        if self.run_control is None:
            raise ValueError("Run control must be loaded before reading time steps")

        specs: List[Tuple[str, str, Tuple[int, ...], str]] = []

        for level in range(1, self.nz + 1):
            level_tag = f"{level:3d}"
            specs.append((f"U-LEV{level_tag}", "real", (self.ny, self.nx), "U"))
            specs.append((f"V-LEV{level_tag}", "real", (self.ny, self.nx), "V"))
            if self.run_control.lcalgrd:
                specs.append((f"WFACE{level_tag}", "real", (self.ny, self.nx), "W"))

        if self.run_control.irtype != 0 and self.run_control.lcalgrd:
            for level in range(1, self.nz + 1):
                specs.append((f"T-LEV{level:3d}", "real", (self.ny, self.nx), "T"))

        if self.run_control.irtype != 0:
            specs.extend(
                [
                    ("IPGT", "int", (self.ny, self.nx), "IPGT"),
                    ("USTAR", "real", (self.ny, self.nx), "USTAR"),
                    ("ZI", "real", (self.ny, self.nx), "ZI"),
                    ("EL", "real", (self.ny, self.nx), "EL"),
                    ("WSTAR", "real", (self.ny, self.nx), "WSTAR"),
                ]
            )
            if self.run_control.npsta != 0:
                specs.append(("RMM", "real", (self.ny, self.nx), "RMM"))
            specs.extend(
                [
                    ("TEMPK", "real", (self.ny, self.nx), "TEMPK"),
                    ("RHO", "real", (self.ny, self.nx), "RHO"),
                    ("QSW", "real", (self.ny, self.nx), "QSW"),
                    ("IRH", "int", (self.ny, self.nx), "IRH"),
                ]
            )
            if self.run_control.npsta != 0:
                specs.append(("IPCODE", "int", (self.ny, self.nx), "IPCODE"))

        return specs

    def info(self) -> str:
        """Return a human-readable metadata summary."""
        lines = [
            f"File: {self.path}",
            f"Dataset: {self.dataset_name} v{self.dataset_version}",
            f"Model: {self.data_model}",
        ]

        if self.run_control is not None:
            rc = self.run_control
            lines.extend(
                [
                    f"Grid: {rc.nx} x {rc.ny} x {rc.nz}",
                    f"Spacing: {rc.dgrid} m",
                    f"Origin: ({rc.xorigr}, {rc.yorigr})",
                    f"Run type: {rc.irtype}",
                    f"CALGRID fields: {rc.lcalgrd}",
                    f"Timezone: {rc.timezone}",
                    f"Run start: {rc.start_time.isoformat(sep=' ')}",
                    f"Run end: {rc.end_time.isoformat(sep=' ')}",
                    f"Timesteps loaded: {self.nt}",
                    f"2D fields: {sorted(self.fields_2d)}",
                    f"3D fields: {sorted(self.fields_3d)}",
                ]
            )

        return "\n".join(lines)

    def get_time_bounds(self) -> List[Tuple[datetime, datetime]]:
        """Return begin/end datetimes for each timestep."""
        return list(self.time_bounds)

    def get_static_field(self, label: str) -> np.ndarray:
        """Return one static field by label, for example ``Z0`` or ``ELEV``."""
        return self.static_fields[label]

    def get_2d_field(self, label: str) -> np.ndarray:
        """Return one time-varying 2D field by label."""
        return self.fields_2d[label]

    def get_3d_field(self, label: str, level: Optional[int] = None) -> np.ndarray:
        """Return one time-varying 3D field by label.

        If `level` is provided, returns the single level with shape
        `[time, y, x]`. Otherwise returns the full field with shape
        `[time, level, y, x]`.
        """
        data = self.fields_3d[label]
        if level is None:
            return data
        if level < 1 or level > data.shape[1]:
            raise IndexError(f"Level must be in 1..{data.shape[1]}")
        return data[:, level - 1, :, :]

    def plot_wind_animation(
        self,
        level: int = 1,
        stride: int = 3,
        interval: int = 200,
        scale: Optional[float] = None,
        cmap: str = "viridis",
    ) -> tuple[object, object]:
        """Build a matplotlib animation of wind speed and vectors."""
        if "U" not in self.fields_3d or "V" not in self.fields_3d:
            raise ValueError("This file does not contain U/V wind fields")

        if level < 1 or level > self.nz:
            raise IndexError(f"Level must be in 1..{self.nz}")

        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation

        u = self.get_3d_field("U", level=level)
        v = self.get_3d_field("V", level=level)
        speed = np.hypot(u, v)

        x = self.x
        y = self.y
        x2d, y2d = np.meshgrid(x, y)
        qx = x2d[::stride, ::stride]
        qy = y2d[::stride, ::stride]

        fig, ax = plt.subplots()
        mesh = ax.pcolormesh(x2d, y2d, speed[0], shading="auto", cmap=cmap)
        quiver = ax.quiver(
            qx,
            qy,
            u[0, ::stride, ::stride],
            v[0, ::stride, ::stride],
            scale=scale,
        )
        title = ax.set_title(self._time_title(0, level))
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        fig.colorbar(mesh, ax=ax, label="Wind speed")

        def update(frame: int):
            mesh.set_array(speed[frame].ravel())
            quiver.set_UVC(
                u[frame, ::stride, ::stride],
                v[frame, ::stride, ::stride],
            )
            title.set_text(self._time_title(frame, level))
            return mesh, quiver, title

        anim = FuncAnimation(fig, update, frames=self.nt, interval=interval, blit=False)
        return fig, anim

    def _time_title(self, frame: int, level: int) -> str:
        begin_time, end_time = self.time_bounds[frame]
        return (
            f"CALMET winds level {level} | "
            f"{begin_time.isoformat(sep=' ')} -> {end_time.isoformat(sep=' ')}"
        )

def _parse_run_control(payload: bytes, endian: str) -> RunControl:
    offset = 0

    ints_10 = struct.unpack_from(f"{endian}10i", payload, offset)
    offset += 40
    axtz = decode_string(payload[offset : offset + 8])
    offset += 8
    irlg, irtype = struct.unpack_from(f"{endian}2i", payload, offset)
    offset += 8
    nx, ny, nz = struct.unpack_from(f"{endian}3i", payload, offset)
    offset += 12
    dgrid, xorigr, yorigr = struct.unpack_from(f"{endian}3f", payload, offset)
    offset += 12
    iwfcod, nssta, nusta, npsta, nowsta, nlu, iwat1, iwat2 = struct.unpack_from(
        f"{endian}8i", payload, offset
    )
    offset += 32
    lcalgrd_raw = struct.unpack_from(f"{endian}i", payload, offset)[0]
    offset += 4
    pmap = decode_string(payload[offset : offset + 8])
    offset += 8
    datum = decode_string(payload[offset : offset + 8])
    offset += 8
    daten = decode_string(payload[offset : offset + 12])
    offset += 12
    feast, fnorth = struct.unpack_from(f"{endian}2f", payload, offset)
    offset += 8
    utmhem = decode_string(payload[offset : offset + 4])
    offset += 4
    iutmzn = struct.unpack_from(f"{endian}i", payload, offset)[0]
    offset += 4
    rnlat0, relon0, xlat1, xlat2 = struct.unpack_from(f"{endian}4f", payload, offset)

    start_time = datetime(
        ints_10[0], ints_10[1], ints_10[2], ints_10[3]
    ) + timedelta(seconds=ints_10[4])
    end_time = datetime(ints_10[5], ints_10[6], ints_10[7], ints_10[8]) + timedelta(
        seconds=ints_10[9]
    )

    return RunControl(
        start_time=start_time,
        end_time=end_time,
        timezone=axtz,
        irlg=irlg,
        irtype=irtype,
        nx=nx,
        ny=ny,
        nz=nz,
        dgrid=dgrid,
        xorigr=xorigr,
        yorigr=yorigr,
        iwfcod=iwfcod,
        nssta=nssta,
        nusta=nusta,
        npsta=npsta,
        nowsta=nowsta,
        nlu=nlu,
        iwat1=iwat1,
        iwat2=iwat2,
        lcalgrd=bool(lcalgrd_raw),
        pmap=pmap,
        datum=datum,
        daten=daten,
        feast=feast,
        fnorth=fnorth,
        utmhem=utmhem,
        iutmzn=iutmzn,
        rnlat0=rnlat0,
        relon0=relon0,
        xlat1=xlat1,
        xlat2=xlat2,
    )

def _parse_yyyyjjjhh(date_code: int, seconds: int) -> datetime:
    """Decode CALMET NDATHR = YYYY*100000 + JJJ*100 + HH.

    The remainder after the 4-digit year is *five* digits (JJJHH), not six.
    Using ``% 1_000_000`` mis-parses years whose last digit is not 0
    (e.g. 2005-08-28 → Julian 5240 / 2019-05-07).
    """
    year = date_code // 100_000
    jday_hour = date_code % 100_000
    jday = jday_hour // 100
    hour = jday_hour % 100
    if jday < 1 or jday > 366 or hour > 24:
        raise ValueError(f"Invalid YYYYJJJHH date code: {date_code}")
    base = datetime(year, 1, 1) + timedelta(days=jday - 1, hours=hour)
    return base + timedelta(seconds=seconds)


def main() -> None:
    parser = argparse.ArgumentParser(description="Read CALMET.DAT files")
    parser.add_argument("path", help="Path to CALMET.DAT")
    parser.add_argument(
        "--animate",
        action="store_true",
        help="Display a matplotlib animation for the wind fields",
    )
    parser.add_argument(
        "--level",
        type=int,
        default=1,
        help="Vertical level used by --animate-winds (1-based)",
    )
    args = parser.parse_args()

    dataset = CalmetDataset.read(args.path)
    print(dataset.info())

    if args.animate:
        import matplotlib.pyplot as plt

        _, anim = dataset.plot_wind_animation(level=args.level)
        # Keep a live reference so matplotlib does not garbage-collect the animation.
        _ = anim
        plt.show()


if __name__ == "__main__":
    main()


def read_calmet_dat(path):
    """Read a CALMET.DAT file into a CalmetDataset."""
    return CalmetDataset.read(path)


# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

class FortranSequentialWriter:
    """Write unformatted Fortran sequential records (4-byte markers)."""

    def __init__(self, path: str | Path, endian: str = "<"):
        self.path = Path(path)
        self.endian = endian
        self.handle: Optional[BinaryIO] = None

    def __enter__(self) -> "FortranSequentialWriter":
        self.handle = self.path.open("wb")
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        if self.handle is not None:
            self.handle.close()
            self.handle = None

    def write_record(self, payload: bytes) -> None:
        assert self.handle is not None
        n = len(payload)
        marker = struct.pack(f"{self.endian}i", n)
        self.handle.write(marker + payload + marker)


def _pack_string(s: str, width: int) -> bytes:
    b = s.encode("ascii", errors="replace")[:width]
    return b + b" " * (width - len(b))


def _yyyyjjjhh(dt: datetime) -> tuple[int, int]:
    jday = int(dt.timetuple().tm_yday)
    code = dt.year * 100_000 + jday * 100 + dt.hour
    return code, dt.minute * 60 + dt.second


def _pack_run_control(rc: RunControl) -> bytes:
    endian = "<"
    parts = []
    st, et = rc.start_time, rc.end_time
    parts.append(
        struct.pack(
            f"{endian}10i",
            st.year, st.month, st.day, st.hour, st.minute * 60 + st.second,
            et.year, et.month, et.day, et.hour, et.minute * 60 + et.second,
        )
    )
    parts.append(_pack_string(rc.timezone, 8))
    parts.append(struct.pack(f"{endian}2i", rc.irlg, rc.irtype))
    parts.append(struct.pack(f"{endian}3i", rc.nx, rc.ny, rc.nz))
    parts.append(struct.pack(f"{endian}3f", rc.dgrid, rc.xorigr, rc.yorigr))
    parts.append(
        struct.pack(
            f"{endian}8i",
            rc.iwfcod, rc.nssta, rc.nusta, rc.npsta, rc.nowsta, rc.nlu, rc.iwat1, rc.iwat2,
        )
    )
    parts.append(struct.pack(f"{endian}i", 1 if rc.lcalgrd else 0))
    parts.append(_pack_string(rc.pmap, 8))
    parts.append(_pack_string(rc.datum, 8))
    parts.append(_pack_string(rc.daten, 12))
    parts.append(struct.pack(f"{endian}2f", rc.feast, rc.fnorth))
    parts.append(_pack_string(rc.utmhem, 4))
    parts.append(struct.pack(f"{endian}i", rc.iutmzn))
    parts.append(struct.pack(f"{endian}4f", rc.rnlat0, rc.relon0, rc.xlat1, rc.xlat2))
    return b"".join(parts)


def _pack_labeled(
    label: str,
    payload: bytes,
    time_info: tuple[int, int, int, int] = (0, 0, 0, 0),
    endian: str = "<",
) -> bytes:
    return _pack_string(label, 8) + struct.pack(f"{endian}4i", *time_info) + payload


def _pack_grid(arr: np.ndarray, dtype_name: str, endian: str = "<") -> bytes:
    """Pack [ny, nx] array in Fortran order (i fastest).

    Reader does: frombuffer → reshape((nx, ny), order='F').T → [ny, nx],
    so we must emit values with index i + j*nx = arr[j, i], i.e. C-ravel of [ny, nx].
    """
    a = np.asarray(arr)
    if a.ndim != 2:
        raise ValueError(f"Expected 2D grid, got {a.shape}")
    dt = np.dtype(f"{endian}f4" if dtype_name == "real" else f"{endian}i4")
    flat = np.asarray(a, dtype=dt).ravel(order="C")
    return flat.tobytes()


def write_calmet_dat(
    path: str | Path,
    *,
    result: "object",
    run_control: RunControl,
    z0: np.ndarray,
    landuse: np.ndarray,
    elev: np.ndarray,
    xlai: np.ndarray | None = None,
    xssta: np.ndarray | None = None,
    yssta: np.ndarray | None = None,
    xusta: np.ndarray | None = None,
    yusta: np.ndarray | None = None,
    nears: np.ndarray | None = None,
    comments: List[str] | None = None,
    rmm: np.ndarray | None = None,
    xpsta: np.ndarray | None = None,
    ypsta: np.ndarray | None = None,
) -> None:
    """Write a CALMET.DAT file readable by :class:`CalmetDataset`.

    ``result`` must expose U,V,W,T [nt,nz,ny,nx] and 2D fields
    IPGT,USTAR,ZI,EL,WSTAR,TEMPK,RHO,QSW,IRH [nt,ny,nx], plus zface.
    """
    path = Path(path)
    endian = "<"
    nx, ny, nz = run_control.nx, run_control.ny, run_control.nz
    nt = int(result.U.shape[0])
    if xlai is None:
        xlai = np.ones((ny, nx), dtype=np.float32)
    if comments is None:
        comments = [
            "Produced by py-calmet write_calmet_dat",
            f"mode={getattr(result, 'mode', 'unknown')}",
        ]

    with FortranSequentialWriter(path, endian=endian) as w:
        w.write_record(
            _pack_string("CALMET.DAT", 16)
            + _pack_string("2.1", 16)
            + _pack_string("No-Obs file structure with embedded control file", 64)
        )
        w.write_record(struct.pack(f"{endian}i", len(comments)))
        for c in comments:
            w.write_record(_pack_string(c, max(len(c), 1)))
        w.write_record(_pack_run_control(run_control))

        # static
        zface = np.asarray(result.zface, dtype=np.float32)
        w.write_record(
            _pack_labeled("ZFACE", np.asarray(zface, dtype=f"{endian}f4").tobytes())
        )
        if run_control.nssta >= 1:
            xs = np.asarray(xssta if xssta is not None else [0.0], dtype=np.float32)
            ys = np.asarray(yssta if yssta is not None else [0.0], dtype=np.float32)
            w.write_record(_pack_labeled("XSSTA", np.asarray(xs, dtype=f"{endian}f4").tobytes()))
            w.write_record(_pack_labeled("YSSTA", np.asarray(ys, dtype=f"{endian}f4").tobytes()))
        if run_control.nusta >= 1:
            xu = np.asarray(xusta if xusta is not None else [0.0], dtype=np.float32)
            yu = np.asarray(yusta if yusta is not None else [0.0], dtype=np.float32)
            w.write_record(_pack_labeled("XUSTA", np.asarray(xu, dtype=f"{endian}f4").tobytes()))
            w.write_record(_pack_labeled("YUSTA", np.asarray(yu, dtype=f"{endian}f4").tobytes()))
        if run_control.npsta >= 1:
            xp = np.asarray(xpsta if xpsta is not None else [0.0], dtype=np.float32)
            yp = np.asarray(ypsta if ypsta is not None else [0.0], dtype=np.float32)
            w.write_record(_pack_labeled("XPSTA", np.asarray(xp, dtype=f"{endian}f4").tobytes()))
            w.write_record(_pack_labeled("YPSTA", np.asarray(yp, dtype=f"{endian}f4").tobytes()))

        for label, arr, kind in (
            ("Z0", z0, "real"),
            ("ILANDU", landuse, "int"),
            ("ELEV", elev, "real"),
            ("XLAI", xlai, "real"),
        ):
            w.write_record(_pack_labeled(label, _pack_grid(arr, kind, endian)))
        if run_control.nssta >= 1:
            near = nears if nears is not None else np.ones((ny, nx), dtype=np.int32)
            w.write_record(_pack_labeled("NEARS", _pack_grid(near, "int", endian)))

        # timesteps
        dt = timedelta(seconds=3600)  # default; override via time_bounds if present
        start = run_control.start_time
        # Prefer evenly spaced from start using irlg
        step_sec = int((run_control.end_time - run_control.start_time).total_seconds() // max(nt, 1))
        if step_sec <= 0:
            step_sec = 3600

        for t in range(nt):
            t0 = start + timedelta(seconds=t * step_sec)
            t1 = t0 + timedelta(seconds=step_sec)
            c0, s0 = _yyyyjjjhh(t0)
            c1, s1 = _yyyyjjjhh(t1)
            tinfo = (c0, s0, c1, s1)

            for lev in range(1, nz + 1):
                w.write_record(
                    _pack_labeled(
                        f"U-LEV{lev:3d}",
                        _pack_grid(result.U[t, lev - 1], "real", endian),
                        tinfo,
                        endian,
                    )
                )
                w.write_record(
                    _pack_labeled(
                        f"V-LEV{lev:3d}",
                        _pack_grid(result.V[t, lev - 1], "real", endian),
                        tinfo,
                        endian,
                    )
                )
                if run_control.lcalgrd:
                    w.write_record(
                        _pack_labeled(
                            f"WFACE{lev:3d}",
                            _pack_grid(result.W[t, lev - 1], "real", endian),
                            tinfo,
                            endian,
                        )
                    )
            if run_control.irtype != 0 and run_control.lcalgrd:
                for lev in range(1, nz + 1):
                    w.write_record(
                        _pack_labeled(
                            f"T-LEV{lev:3d}",
                            _pack_grid(result.T[t, lev - 1], "real", endian),
                            tinfo,
                            endian,
                        )
                    )
            if run_control.irtype != 0:
                for label, arr, kind in (
                    ("IPGT", result.IPGT[t], "int"),
                    ("USTAR", result.USTAR[t], "real"),
                    ("ZI", result.ZI[t], "real"),
                    ("EL", result.EL[t], "real"),
                    ("WSTAR", result.WSTAR[t], "real"),
                ):
                    w.write_record(
                        _pack_labeled(label, _pack_grid(arr, kind, endian), tinfo, endian)
                    )
                if run_control.npsta != 0:
                    rr = rmm[t] if rmm is not None else np.zeros((ny, nx), dtype=np.float32)
                    w.write_record(
                        _pack_labeled("RMM", _pack_grid(rr, "real", endian), tinfo, endian)
                    )
                for label, arr, kind in (
                    ("TEMPK", result.TEMPK[t], "real"),
                    ("RHO", result.RHO[t], "real"),
                    ("QSW", result.QSW[t], "real"),
                    ("IRH", result.IRH[t], "int"),
                ):
                    w.write_record(
                        _pack_labeled(label, _pack_grid(arr, kind, endian), tinfo, endian)
                    )
                if run_control.npsta != 0:
                    w.write_record(
                        _pack_labeled(
                            "IPCODE",
                            _pack_grid(np.zeros((ny, nx), dtype=np.int32), "int", endian),
                            tinfo,
                            endian,
                        )
                    )


def write_calmet_netcdf(path: str | Path, result: "object", run_control: RunControl | None = None) -> None:
    """Write a CF-ish NetCDF of CalmetResult fields (optional dependency)."""
    try:
        from netCDF4 import Dataset
    except ImportError as exc:  # pragma: no cover
        raise ImportError("netCDF4 required for write_calmet_netcdf") from exc

    path = Path(path)
    nt, nz, ny, nx = result.U.shape
    with Dataset(path, "w") as ds:
        ds.createDimension("time", nt)
        ds.createDimension("level", nz)
        ds.createDimension("y", ny)
        ds.createDimension("x", nx)
        ds.createDimension("zface", nz + 1)
        ds.createVariable("zface", "f4", ("zface",))[:] = np.asarray(result.zface, dtype=np.float32)
        for name in ("U", "V", "W", "T"):
            var = ds.createVariable(name, "f4", ("time", "level", "y", "x"))
            var[:] = np.asarray(getattr(result, name), dtype=np.float32)
        for name in ("IPGT", "USTAR", "ZI", "EL", "WSTAR", "TEMPK", "RHO", "QSW", "IRH"):
            var = ds.createVariable(name, "f4", ("time", "y", "x"))
            var[:] = np.asarray(getattr(result, name), dtype=np.float32)
        ds.setncattr("title", "py-calmet output")
        ds.setncattr("mode", getattr(result, "mode", ""))
