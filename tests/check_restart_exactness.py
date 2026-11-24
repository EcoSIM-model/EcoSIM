"""Compare a short uninterrupted EcoSIM run with a stopped/resumed run.

Requires numpy and netCDF4. Run from the repository root, for example:
  python3 tests/check_restart_exactness.py --executable build/<config>/drivers/ecosim/ecosim.f90.x \
      --workdir /tmp/ecosim-restart-check

The work directory must not exist. Inputs and existing outputs are read-only.
Use --split-hours 53 to exercise an intraday restart, or --seed CHECKPOINT to
start both branches from an existing mature checkpoint. Old checkpoint recovery
is then common to both branches; only the newly written checkpoint is tested
for exact continuation. No reference baselines are updated.
"""

import argparse
from datetime import datetime, timedelta
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from netCDF4 import Dataset


ROOT = Path(__file__).resolve().parents[1]


def setting(text, key, value):
    pattern = rf"(?m)^(\s*{key}\s*=)[^\n]*"
    updated, count = re.subn(pattern, lambda match: match[1] + value, text)
    if count != 1:
        raise ValueError(f"Expected one {key} assignment, found {count}")
    return updated


def compare(left, right):
    differences = []
    with Dataset(left) as a, Dataset(right) as b:
        a.set_auto_maskandscale(False)
        b.set_auto_maskandscale(False)
        if {name: len(dim) for name, dim in a.dimensions.items()} != {
                name: len(dim) for name, dim in b.dimensions.items()}:
            differences.append({"dimension_sizes_mismatch": True})
        if set(a.variables) != set(b.variables):
            differences.append({"variable_sets": [sorted(a.variables), sorted(b.variables)]})
        for name in a.variables.keys() & b.variables.keys():
            x, y = a[name][:], b[name][:]
            if x.shape != y.shape or x.dtype != y.dtype:
                differences.append({"variable": name, "schema_mismatch": True})
                continue
            equal = np.array_equal(x, y, equal_nan=True) if x.dtype.kind in "fc" else np.array_equal(x, y)
            if not equal:
                entry = {"variable": name, "different_elements": int(np.count_nonzero(x != y))}
                if x.dtype.kind in "fiu":
                    entry["max_abs_difference"] = float(np.max(np.abs(x.astype(float) - y.astype(float))))
                differences.append(entry)
            if x.dtype.kind in "fc":
                for label, values, variable in (("reference", x, a[name]), ("resumed", y, b[name])):
                    valid = np.ones(values.shape, dtype=bool)
                    for attr in ("_FillValue", "missing_value"):
                        if attr in variable.ncattrs():
                            fill = variable.getncattr(attr)
                            valid &= ~np.isnan(values) if np.any(np.isnan(fill)) else values != fill
                    if np.any(~np.isfinite(values[valid])):
                        differences.append({"variable": name, "nonfinite": label})
    return sorted(differences, key=lambda item: item.get("variable", ""))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--namelist", type=Path, default=ROOT / "examples/run_dir/blodgett/Blodget.dfir.namelist")
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--days", type=int, default=5)
    parser.add_argument("--split-hours", type=int, default=48)
    parser.add_argument("--seed", type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    namelist = args.namelist.resolve(strict=True)
    work = args.workdir.resolve()
    if not 0 < args.split_hours < args.days * 24:
        parser.error("split-hours must lie within the short run")
    source = namelist.read_text()
    # Preserve each case's input selection while moving execution to scratch.
    def absolute_input(match):
        value = match[1].rstrip()  # Fortran trims trailing blanks in filenames.
        candidate = namelist.parent / value
        return "'" + str(candidate.resolve()) + "'" if candidate.is_file() else match[0]
    source = re.sub(r"'([^'\n]+)'", absolute_input, source)
    if float(re.search(r"(?m)^\s*delta_time\s*=\s*([\d.]+)", source)[1]) != 3600:
        parser.error("this driver test requires hourly timesteps")
    origin = datetime.strptime(re.search(r"(?m)^\s*start_date\s*=\s*'(\d+)'", source)[1], "%Y%m%d%H%M%S")
    start = origin
    seed = args.seed.resolve(strict=True) if args.seed else None
    if seed:
        match = re.search(r"\.r\.(\d{4}-\d{2}-\d{2})-(\d{6})\.nc$", seed.name)
        start = datetime.strptime(match[1], "%Y-%m-%d") + timedelta(seconds=int(match[2]))
    # The standalone driver counts steps from January 1 of the start year.
    offset = int((start - datetime(origin.year, 1, 1)).total_seconds() // 3600)
    if offset < 0:
        parser.error("seed precedes the configured start_date")
    work.mkdir(parents=True, exist_ok=False)

    def copy_checkpoint(checkpoint, destination):
        shutil.copy2(checkpoint, destination / checkpoint.name)
        prefix, stamp = checkpoint.name.split(".ecosim.r.")
        for path in checkpoint.parent.glob(prefix + ".ecosim.rh*." + stamp):
            shutil.copy2(path, destination / path.name)
        for path in checkpoint.parent.glob(prefix + ".ecosim.h*.nc"):
            shutil.copy2(path, destination / path.name)
        (destination / "rpointer.esim").write_text(checkpoint.name + "\n")

    def run(name, hours, checkpoint=None):
        folder = work / name
        folder.mkdir()
        if checkpoint:
            copy_checkpoint(checkpoint, folder)
        text = source
        for key, value in {"stop_option": "'nsteps'", "stop_n": str(offset + hours),
                           "rest_opt": "'never'", "rest_frq": "1",
                           "continue_run": ".true." if checkpoint else ".false."}.items():
            text = setting(text, key, value)
        (folder / "case.namelist").write_text(text)
        with (folder / "run.log").open("w") as log:
            subprocess.run([str(executable), "case.namelist"], cwd=folder,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        steps = re.findall(r"wrote out restart data at nstep\s*=\s*(\d+)",
                           (folder / "run.log").read_text(errors="replace"))
        if not steps or int(steps[-1]) != offset + hours:
            raise RuntimeError(f"Incomplete run; inspect {folder / 'run.log'}")
        date = start + timedelta(hours=hours)
        stamp = date.strftime("%Y-%m-%d-") + f"{date.hour * 3600:06d}.nc"
        files = list(folder.glob("*.ecosim.r." + stamp))
        if len(files) != 1:
            raise RuntimeError(f"Missing final checkpoint; inspect {folder / 'run.log'}")
        # ENDRUN can return status zero after a NetCDF error. Validate the file.
        with Dataset(files[0]) as dataset:
            if not dataset.variables:
                raise RuntimeError(f"Empty checkpoint; inspect {folder / 'run.log'}")
        return files[0]

    reference = run("continuous", args.days * 24, seed)
    split = run("first_segment", args.split_hours, seed)
    resumed = run("resumed", args.days * 24, split)
    results = {"checkpoint": compare(reference, resumed)}
    # History includes partially accumulated samples across intraday restarts.
    for path in resumed.parent.glob("*.ecosim.h*.nc"):
        if not (reference.parent / path.name).exists():
            results[path.name] = [{"unexpected_file": str(path)}]
    for path in reference.parent.glob("*.ecosim.h*.nc"):
        other = resumed.parent / path.name
        results[path.name] = compare(path, other) if other.exists() else [{"missing_file": str(other)}]
    report = {"executable": str(executable), "namelist": str(namelist), "seed": str(seed) if seed else None,
              "days": args.days, "split_hours": args.split_hours, "comparisons": results}
    (work / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    count = sum(map(len, results.values()))
    print(f"{'PASS' if count == 0 else 'FAIL'}: {count} differences; {work / 'comparison.json'}")
    return int(count != 0)


if __name__ == "__main__":
    raise SystemExit(main())
