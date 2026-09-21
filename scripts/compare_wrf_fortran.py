#!/usr/bin/env python3
"""Compare py-calmet vs Fortran on the real-WRF wrf_demo case (+ GEO.DAT).

Primary validation entrypoint. Exits non-zero if any wrf_demo mode fails
thresholds in ``tests/thresholds.py`` (WRF_DEMO_THRESH).

Reference preference
--------------------
1. Live Fortran ``calmet.x`` when available (``CALMET_X`` env, or
   ``vendor/calmet-fortran/src/calmet.x``, or per-mode symlink under
   ``cases/wrf_demo/{mode}/calmet.x``), re-run into a work dir and compare.
2. Else archived Fortran-from-WRF-case outputs under
   ``cases/wrf_demo/goldens/{mode}/CALMET.DAT`` — still the Katrina
   mountain-window case, **not** synthetic small_domain goldens.

Usage
-----
  python scripts/compare_wrf_fortran.py
  python scripts/compare_wrf_fortran.py --modes noobs,obs_model
  python scripts/compare_wrf_fortran.py --json /tmp/wrf_compare.json
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "tests"))

from py_calmet import run_calmet, read_calmet_dat  # noqa: E402
from py_calmet.core.met_utils import relative_rmse  # noqa: E402
from thresholds import WRF_DEMO_THRESH, WRF_DEMO_UV_FLOOR  # noqa: E402

CASE = REPO / "cases" / "wrf_demo"
GOLDENS = CASE / "goldens"
INPUTS = GOLDENS / "inputs"
SHARED = CASE / "shared"


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    if a.std() < 1e-9 or b.std() < 1e-9:
        return 0.0
    return float(np.corrcoef(a, b)[0, 1])


def find_calmet_x() -> Path | None:
    env = os.environ.get("CALMET_X")
    if env and Path(env).is_file():
        return Path(env)
    vendored = REPO / "vendor" / "calmet-fortran" / "src" / "calmet.x"
    if vendored.is_file():
        return vendored
    for mode in ("noobs", "obs_model", "obs"):
        p = CASE / mode / "calmet.x"
        if p.is_file() or p.is_symlink():
            try:
                if p.resolve().is_file():
                    return p.resolve()
            except FileNotFoundError:
                pass
    return None


def ensure_geo() -> None:
    geo = INPUTS / "geo.dat"
    if not geo.is_file():
        raise SystemExit(f"GEO.DAT missing on wrf_demo path: {geo}")
    text = geo.read_text()
    if "Katrina" not in text and "wrfout" not in text.lower():
        raise SystemExit(f"{geo} does not look like Katrina wrfout GEO")
    # mountain check
    elev_max = 0.0
    saw_ht = False
    for ln in text.splitlines():
        if "HTFAC" in ln:
            saw_ht = True
            continue
        if not saw_ht:
            continue
        parts = ln.split()
        if not parts:
            break
        try:
            elev_max = max(elev_max, max(float(x) for x in parts))
        except ValueError:
            break
    if elev_max < 1000.0:
        raise SystemExit(f"geo.dat max elev {elev_max} < 1000 m — expected mountains")


def run_live_fortran(mode: str, calmet_x: Path, work: Path) -> Path:
    """Re-run Fortran for one mode; return path to CALMET.DAT."""
    work.mkdir(parents=True, exist_ok=True)
    inp_src = GOLDENS / mode / "calmet.inp"
    if not inp_src.is_file():
        raise SystemExit(f"missing {inp_src} for live Fortran")
    shutil.copy2(inp_src, work / "calmet.inp")
    for name in ("geo.dat", "surf.dat", "up.dat", "3d.dat"):
        src = SHARED / name
        dst = work / name
        if dst.exists() or dst.is_symlink():
            dst.unlink()
        dst.symlink_to(src.resolve())
    cx = work / "calmet.x"
    if cx.exists() or cx.is_symlink():
        cx.unlink()
    cx.symlink_to(calmet_x.resolve())
    log = work / "run.log"
    with log.open("w") as fh:
        proc = subprocess.run(
            [str(cx), "calmet.inp"],
            cwd=work,
            stdout=fh,
            stderr=subprocess.STDOUT,
            check=False,
        )
    out = work / "CALMET.DAT"
    if not out.is_file():
        out = work / "calmet.dat"
    if proc.returncode != 0 or not out.is_file():
        raise SystemExit(
            f"live Fortran failed mode={mode} rc={proc.returncode}; see {log}"
        )
    return out


def compare_mode(mode: str, gold_dat: Path) -> dict:
    thr = WRF_DEMO_THRESH[mode]
    fl = WRF_DEMO_UV_FLOOR
    res = run_calmet(GOLDENS / mode, mode=mode, inputs_dir=INPUTS)
    gold = read_calmet_dat(gold_dat)
    Ug = gold.get_3d_field("U")
    Vg = gold.get_3d_field("V")
    Zi = gold.get_2d_field("ZI")
    stats = {
        "U": relative_rmse(res.U, Ug, floor=fl),
        "V": relative_rmse(res.V, Vg, floor=fl),
        "U_corr": _corr(res.U, Ug),
        "V_corr": _corr(res.V, Vg),
        "ZI": relative_rmse(res.ZI, Zi),
        "USTAR": relative_rmse(res.USTAR, gold.get_2d_field("USTAR")),
        "SPD": relative_rmse(np.hypot(res.U, res.V), np.hypot(Ug, Vg), floor=fl),
        "abs_U_rmse": float(np.sqrt(np.mean((res.U - Ug) ** 2))),
        "abs_V_rmse": float(np.sqrt(np.mean((res.V - Vg) ** 2))),
        "abs_ZI_rmse": float(np.sqrt(np.mean((res.ZI - Zi) ** 2))),
    }
    fails = []
    for key in ("U", "V", "ZI", "USTAR", "SPD"):
        if stats[key] > thr[key]:
            fails.append(f"{key}={stats[key]:.4f}>{thr[key]}")
    if stats["U_corr"] < thr["U_corr"]:
        fails.append(f"U_corr={stats['U_corr']:.3f}<{thr['U_corr']}")
    if stats["V_corr"] < thr["V_corr"]:
        fails.append(f"V_corr={stats['V_corr']:.3f}<{thr['V_corr']}")
    return {"mode": mode, "stats": stats, "thr": thr, "fails": fails, "ok": not fails}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--modes",
        default="noobs,obs_model,obs",
        help="comma-separated modes (default: all three)",
    )
    ap.add_argument("--json", type=Path, help="optional JSON report path")
    ap.add_argument(
        "--work",
        type=Path,
        default=CASE / "work" / "live_fortran",
        help="scratch dir for live Fortran re-runs",
    )
    args = ap.parse_args()
    modes = [m.strip() for m in args.modes.split(",") if m.strip()]

    ensure_geo()
    calmet_x = find_calmet_x()
    if calmet_x:
        ref_label = f"live Fortran ({calmet_x})"
    else:
        ref_label = (
            "archived Fortran-from-WRF-case goldens "
            "(cases/wrf_demo/goldens/*/CALMET.DAT; not synthetic domain)"
        )

    print("=== wrf_demo PRIMARY validation ===")
    print(f"GEO.DAT: {INPUTS / 'geo.dat'}")
    print(f"3D.DAT:  {INPUTS / '3d.dat'}")
    print(f"Fortran reference: {ref_label}")
    print(f"UV floor (rel RMSE): {WRF_DEMO_UV_FLOOR} m/s")
    print()

    report = {
        "reference": ref_label,
        "live_fortran": bool(calmet_x),
        "geo": str(INPUTS / "geo.dat"),
        "modes": [],
    }
    failed = False
    for mode in modes:
        if calmet_x:
            gold_dat = run_live_fortran(mode, calmet_x, args.work / mode)
        else:
            gold_dat = GOLDENS / mode / "CALMET.DAT"
            if not gold_dat.is_file():
                print(f"FAIL {mode}: missing golden {gold_dat}")
                failed = True
                continue
        row = compare_mode(mode, gold_dat)
        report["modes"].append(row)
        s = row["stats"]
        status = "PASS" if row["ok"] else "FAIL"
        if not row["ok"]:
            failed = True
        print(
            f"{status} {mode:10s}  "
            f"U={s['U']:.4f} (corr {s['U_corr']:.3f})  "
            f"V={s['V']:.4f} (corr {s['V_corr']:.3f})  "
            f"ZI={s['ZI']:.4f}  "
            f"|abs U/V/ZI RMSE {s['abs_U_rmse']:.3f}/{s['abs_V_rmse']:.3f}/{s['abs_ZI_rmse']:.1f}"
        )
        if row["fails"]:
            print(f"       exceeds: {', '.join(row['fails'])}")

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2) + "\n")
        print(f"\nwrote {args.json}")

    print()
    if failed:
        print("PRIMARY GATE FAILED — wrf_demo vs Fortran")
        return 1
    print("PRIMARY GATE PASSED — wrf_demo vs Fortran (real WRF + GEO.DAT)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
