"""Compare annual rotation branches with a cold-start continuous simulation.

Use a prepared transient-management namelist (including fallow records if desired).
All runs and outputs are isolated under --workdir; source inputs are read-only.
The test branches at January checkpoints, which precede annual management reads.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np
from netCDF4 import Dataset

from check_restart_exactness import compare


ROOT = Path(__file__).resolve().parents[1]
STAGER = ROOT / ".agents/skills/ecosim-restart-run/scripts/stage_restart_run.py"
spec = importlib.util.spec_from_file_location("stage_restart", STAGER)
stage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage)


def compare_history(reference, branch, offset):
    differences = []
    with Dataset(reference) as a, Dataset(branch) as b:
        a.set_auto_maskandscale(False)
        b.set_auto_maskandscale(False)
        if set(a.variables) != set(b.variables):
            differences.append({"variable_sets_differ": True})
        for name in a.variables.keys() & b.variables.keys():
            av, bv = a[name], b[name]
            x, y = av[:], bv[:]
            if "time" in av.dimensions:
                if av.dimensions[0] != "time":
                    raise ValueError(f"Unexpected time ordering: {name}")
                x = x[offset:offset + len(y)]
            if name == "time_bounds":
                y = y + offset  # Branch history is relative to its new start date.
            if x.shape != y.shape or x.dtype != y.dtype or av.dimensions != bv.dimensions:
                differences.append({"variable": name, "schema_mismatch": True})
                continue
            equal = np.array_equal(x, y, equal_nan=True) if x.dtype.kind in "fc" else np.array_equal(x, y)
            if not equal:
                differences.append({"variable": name, "different_elements": int(np.count_nonzero(x != y))})
            if x.dtype.kind in "fc":
                for label, values, var in (("continuous", x, av), ("branch", y, bv)):
                    valid = np.ones(values.shape, dtype=bool)
                    for attr in ("_FillValue", "missing_value"):
                        if attr in var.ncattrs():
                            fill = var.getncattr(attr)
                            valid &= ~np.isnan(values) if np.any(np.isnan(fill)) else values != fill
                    if np.any(~np.isfinite(values[valid])):
                        differences.append({"variable": name, "nonfinite": label})
    return differences


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--namelist", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--years", type=int, default=3)
    args = parser.parse_args()
    if args.years < 2:
        parser.error("at least two years are required")
    executable = args.executable.resolve(strict=True)
    namelist = args.namelist.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = namelist.read_text()
    # Existing example namelists use inline comments on solver-count arrays.
    source = re.sub(r"(?m)![^\n]*$", "", source)

    def absolute_input(match):
        candidate = namelist.parent / match[1].rstrip()
        return "'" + str(candidate.resolve()) + "'" if candidate.is_file() else match[0]

    source = re.sub(r"'([^'\n]+)'", absolute_input, source)
    start = datetime.strptime(stage.get_value(source, "start_date").strip("'\""), "%Y%m%d%H%M%S")
    if (start.month, start.day, start.hour, start.minute, start.second) != (1, 1, 0, 0, 0):
        parser.error("the continuous simulation must start at midnight on January 1")
    for key, value, group in (
        ("case_name", "'annual'", "ecosim"), ("finidat", "''", "ecosim"),
        ("continue_run", ".false.", "ecosim"), ("hist_nhtfrq", "-24", "ecosim"),
        ("hist_mfilt", str(args.years * 366 + 1), "ecosim"),
        ("stop_option", "'nyears'", "ecosim_time"), ("stop_n", str(args.years), "ecosim_time"),
        ("rest_opt", "'nyears'", "ecosim_time"), ("rest_frq", "1", "ecosim_time"),
    ):
        source = stage.set_value(source, key, value, group)
    parent = work / "continuous"
    parent.mkdir()
    (parent / "runfile.namelist").write_text(source)

    def run(folder):
        with (folder / "run.log").open("w") as log:
            result = subprocess.run([str(executable), "runfile.namelist"], cwd=folder,
                                    stdout=log, stderr=subprocess.STDOUT)
        body = (folder / "run.log").read_text(errors="replace")
        if result.returncode or "ENDRUN" in body or "Fortran runtime error" in body:
            raise RuntimeError(f"Model failed: {folder / 'run.log'}")
        print(f"Completed {folder.name}", flush=True)

    def branch(year, continuous):
        checkpoint = parent / f"annual.ecosim.r.{year}-01-01-000000.nc"
        pointer = parent / "rpointer.esim"
        while True:
            # The pointer is written after the checkpoint has been closed.
            if checkpoint.is_file() and pointer.is_file() and pointer.stat().st_mtime >= checkpoint.stat().st_mtime:
                if pointer.read_text().strip():
                    break
            if continuous.done():
                continuous.result()
                raise RuntimeError(f"Missing annual checkpoint: {checkpoint}")
            time.sleep(1)
        folder = work / f"branch-{year}"
        subprocess.run([sys.executable, str(STAGER), "--namelist", str(parent / "runfile.namelist"),
                        "--restart", str(checkpoint), "--out", str(folder),
                        "--case-name", f"annual.b{year}", "--stop-n", str(start.year + args.years - year)],
                       check=True, stdout=subprocess.DEVNULL)
        run(folder)
        return folder

    with ThreadPoolExecutor(max_workers=args.years) as pool:
        continuous = pool.submit(run, parent)
        futures = [pool.submit(branch, year, continuous)
                   for year in range(start.year + 1, start.year + args.years)]
        continuous.result()
        branches = [future.result() for future in futures]
    report = {"executable": str(executable), "namelist": str(namelist), "years": args.years, "branches": {}}
    reference_history, = parent.glob("*.ecosim.h0.*.nc")
    end = start.year + args.years
    for folder in branches:
        year = int(folder.name.split("-")[-1])
        offset = (datetime(year, 1, 1) - start).days
        branch_history, = folder.glob("*.ecosim.h0.*.nc")
        with Dataset(branch_history) as nc:
            expected_days = (datetime(end, 1, 1) - datetime(year, 1, 1)).days
            if len(nc.dimensions["time"]) != expected_days:
                raise RuntimeError(f"Incomplete branch history: {branch_history}")
        reference_checkpoint, = parent.glob(f"*.ecosim.r.{end}-01-01-*.nc")
        branch_checkpoint, = folder.glob(f"*.ecosim.r.{end}-01-01-*.nc")
        checkpoint_differences = compare(reference_checkpoint, branch_checkpoint)
        # Branch origin and history filenames intentionally differ. Every physical
        # state is compared; these are routing metadata, not model variables.
        checkpoint_differences = [d for d in checkpoint_differences
                                  if d.get("variable") not in
                                  ("start_date", "timemgr_rst_start_ymd", "locfnh", "locfnhr")]
        result = {"checkpoint": checkpoint_differences,
                  "history": compare_history(reference_history, branch_history, offset)}
        report["branches"][str(year)] = result
        print(f"{folder.name}: checkpoint differences={len(result['checkpoint'])}, "
              f"history differences={len(result['history'])}", flush=True)
    (work / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    return int(any(items for result in report["branches"].values() for items in result.values()))


if __name__ == "__main__":
    raise SystemExit(main())
