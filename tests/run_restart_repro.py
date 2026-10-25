#!/usr/bin/env python3
"""Compare fresh and split EcoSIM runs, including partial-day history recovery.

Requires numpy/netCDF4. Run from the repository root, for example:
  python tests/run_restart_repro.py --executable build/CONFIG/drivers/ecosim/ecosim.f90.x \
      --namelist examples/run_dir/blodgett/Blodget.ndlf.lig.namelist \
      --output build/runs/restart-regression
All model writes go to a new scratch mirror. No baselines are updated.
"""
from pathlib import Path
import argparse
import datetime as dt
import hashlib
import json
import math
import re
import shutil
import subprocess

import numpy as np
from netCDF4 import Dataset


def set_value(text, name, value):
    pattern = rf"(?m)^\s*{name}\s*=.*$"
    if not re.search(pattern, text):
        raise ValueError(f"Missing namelist assignment: {name}")
    return re.sub(pattern, f"{name}={value}", text)


def compare_numeric(reference, candidate):
    different = []
    count = 0
    with Dataset(reference) as a, Dataset(candidate) as b:
        if set(a.variables) != set(b.variables):
            raise AssertionError(f"Different variable sets: {candidate}")
        for name, var in a.variables.items():
            if np.issubdtype(var.dtype, np.number):
                count += 1
                ax = np.ma.asarray(var[:])
                ay = np.ma.asarray(b[name][:])
                finite = np.isfinite(ax.compressed()).all() and np.isfinite(ay.compressed()).all()
                same_schema = var.dtype == b[name].dtype and var.dimensions == b[name].dimensions
                same_mask = np.array_equal(np.ma.getmaskarray(ax), np.ma.getmaskarray(ay))
                same_data = np.array_equal(ax.filled(0), ay.filled(0), equal_nan=True)
                if not finite or not same_schema or not same_mask or not same_data:
                    different.append(name)
    return {"compared_fields": count, "different_fields": different}


def copy_history_prefix(source, destination, elapsed_days):
    """Exclude future records so copied baseline output cannot conceal a failure."""
    with Dataset(source) as src:
        keep = int(np.count_nonzero(np.asarray(src['time_bounds'][:])[:, 1]
                                    <= elapsed_days))
        if not keep:
            return
        with Dataset(destination, 'w', format=src.data_model) as dst:
            dst.setncatts({n: src.getncattr(n) for n in src.ncattrs()})
            for name, dimension in src.dimensions.items():
                dst.createDimension(name, None if dimension.isunlimited() else len(dimension))
            for name, var in src.variables.items():
                kwargs = {'fill_value': var.getncattr('_FillValue')} if '_FillValue' in var.ncattrs() else {}
                out = dst.createVariable(name, var.dtype, var.dimensions, **kwargs)
                out.setncatts({n: var.getncattr(n) for n in var.ncattrs() if n != '_FillValue'})
                index = tuple(slice(0, keep) if n == 'time' else slice(None)
                              for n in var.dimensions)
                out[:] = var[index]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--namelist', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--days', type=int, default=20)
    parser.add_argument('--split-hours', type=int, nargs='+', default=[240, 252])
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    executable = args.executable.resolve()
    original = args.namelist.resolve().read_text()
    out = args.output.resolve()
    if out.exists():
        raise ValueError(f"Use a new output directory: {out}")
    if any(h <= 0 or h >= args.days * 24 for h in args.split_hours):
        raise ValueError('Splits must lie inside the simulation window')
    code = '\n'.join(line.split('!')[0] for line in original.splitlines())
    if re.search(r"finidat\s*=\s*'[^ '\n]+", code, re.I):
        raise ValueError('This regression runner expects a fresh-start case')
    start = re.search(r"start_date\s*=\s*'(\d{14})'", code).group(1)
    start = dt.datetime.strptime(start, '%Y%m%d%H%M%S')
    delta = float(re.search(r'delta_time\s*=\s*([\d.]+)', code).group(1))
    if delta != 3600 or start.hour or start.minute or start.second:
        raise ValueError('This runner requires hourly steps beginning at midnight')
    text = set_value(original, 'continue_run', '.false.')
    text = set_value(text, 'stop_option', "'nsteps'")
    text = set_value(text, 'stop_n', str(args.days * 24))
    text = set_value(text, 'rest_opt', "'nsteps'")
    interval = math.gcd(*args.split_hours)
    text = set_value(text, 'rest_frq', str(interval))

    def stage(name):
        mirror = out / name / 'eco'
        case = mirror / 'examples/run_dir' / args.namelist.parent.name
        case.mkdir(parents=True)
        (mirror / 'input_data').symlink_to(root / 'input_data', target_is_directory=True)
        (mirror / 'examples/inputs').symlink_to(root / 'examples/inputs', target_is_directory=True)
        return case

    def run(case, namelist):
        (case / 'test.namelist').write_text(namelist)
        log_path = case.parents[3] / 'run.log'
        with log_path.open('w') as log:
            result = subprocess.run([str(executable), 'test.namelist'], cwd=case,
                                    stdout=log, stderr=subprocess.STDOUT)
        log = log_path.read_text(errors='replace')
        steps = re.findall(r'wrote out restart data at nstep\s*=\s*(\d+)', log)
        if result.returncode or 'ENDRUN:' in log or not steps or int(steps[-1]) != args.days * 24:
            raise AssertionError(f"Run failed or stopped early; see {log_path}")

    baseline = stage('continuous')
    run(baseline, text)
    def restart_stamp(date):
        seconds = date.hour * 3600 + date.minute * 60 + date.second
        return f'{date:%Y-%m-%d}-{seconds:06d}'

    final_date = restart_stamp(start + dt.timedelta(days=args.days))
    final = next(baseline.glob(f'*.r.{final_date}.nc'))
    results = {'executable': str(executable), 'executable_sha256': hashlib.sha256(executable.read_bytes()).hexdigest(),
               'namelist': str(args.namelist.resolve()), 'days': args.days, 'splits': []}
    for hour in args.split_hours:
        stamp = restart_stamp(start + dt.timedelta(hours=hour))
        checkpoint = next(baseline.glob(f'*.r.{stamp}.nc'))
        case = stage(f'restart_{hour:06d}h')
        shutil.copy2(checkpoint, case / checkpoint.name)
        for path in baseline.glob(f'*.rh*.{stamp}.nc'):
            shutil.copy2(path, case / path.name)
        for path in baseline.glob('*.h0.*.nc'):
            copy_history_prefix(path, case / path.name, hour / 24)
        (case / 'rpointer.esim').write_text('./' + checkpoint.name + '\n')
        run(case, set_value(text, 'continue_run', '.true.'))
        record = {'split_hour': hour, 'restart': compare_numeric(final, case / final.name),
                  'history': {p.name: compare_numeric(p, case / p.name)
                              for p in baseline.glob('*.h0.*.nc')}}
        record['pass'] = not record['restart']['different_fields'] and all(
            not v['different_fields'] for v in record['history'].values())
        results['splits'].append(record)
        (out / 'results.json').write_text(json.dumps(results, indent=2) + '\n')
        print(f"Split at {hour} h: {'PASS' if record['pass'] else 'FAIL'}", flush=True)
    if not all(r['pass'] for r in results['splits']):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
