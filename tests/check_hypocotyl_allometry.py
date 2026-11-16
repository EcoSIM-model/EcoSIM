"""Integration check of explicit HYPA loading, legacy equivalence and emergence response.

Use a soybean checkpoint before planting and a parent namelist with hourly output.
HYPA is required and uses the standard trait reader, like WDLF.
All simulations and modified parameter fixtures live in a NEW scratch directory.
The optional legacy executable tests exact numerical agreement with the old code.
Example (from repository root):
  python tests/check_hypocotyl_allometry.py --executable build/CONFIG/drivers/ecosim/ecosim.f90.x \
    --namelist PARENT/runfile.namelist --restart PARENT/CASE.ecosim.r.2022-05-01-000000.nc \
    --workdir /tmp/hypa-check --legacy-executable /tmp/ecosim_before.x
"""
import argparse
import json
import re
from pathlib import Path
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

from check_restart_exactness import compare

ROOT = Path(__file__).resolve().parents[1]
STAGER = ROOT / '.agents/skills/ecosim-restart-run/scripts/stage_restart_run.py'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--legacy-executable', type=Path)
    parser.add_argument('--namelist', type=Path, required=True)
    parser.add_argument('--restart', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--days', type=int, default=20)
    args = parser.parse_args()
    args.workdir.mkdir(parents=True, exist_ok=False)
    results = {}
    runs = {}
    fixtures = [('explicit100', 100.), ('lower25', 25.)]
    if args.legacy_executable:
        fixtures.insert(0, ('legacy', 100.))
    for label, value in fixtures:
        folder = args.workdir / label
        stage = [sys.executable, str(STAGER), '--namelist', str(args.namelist.resolve()),
                 '--restart', str(args.restart.resolve()), '--out', str(folder),
                 '--case-name', 'hypa.check', '--stop-n', str(args.days),
                 '--stop-option', 'ndays', '--copy-input', 'pft_file_in',
                 '--add-fields', 'EmergHeight_pft,POPN_pft']
        subprocess.run(stage, check=True, capture_output=True, text=True)
        manifest = json.loads((folder / 'restart_run.json').read_text())
        # The copied input is the only parameter fixture edited by this test.
        text = (folder / 'runfile.namelist').read_text()
        pft_path = Path(re.search(r"(?im)^\s*pft_file_in\s*=\s*'([^']+)'", text)[1])
        if not pft_path.is_absolute():
            pft_path = folder / pft_path
        with Dataset(pft_path, 'a') as ds:
            if 'HYPA' in ds.variables:
                var = ds['HYPA']
            else:
                var = ds.createVariable('HYPA', 'f8', ds['WDLF'].dimensions)
                var.long_name = 'Hypocotyl-leaf area allometry coefficient'
                var.units = '1'
            var[:] = value
        subprocess.run([sys.executable, str(ROOT / 'ecosim-agent/Tools/case_preflight.py'),
                        str(folder / 'runfile.namelist'), '--output', str(folder / 'preflight.json')],
                       check=True, capture_output=True, text=True)
        exe = args.legacy_executable if label == 'legacy' else args.executable
        with (folder / 'model.log').open('w') as log:
            proc = subprocess.run([str(exe.resolve()), 'runfile.namelist'], cwd=folder,
                                  stdout=log, stderr=subprocess.STDOUT)
        log = (folder / 'model.log').read_text(errors='replace')
        checkpoints = sorted(folder.glob('*.ecosim.r.*.nc'))
        steps = re.findall(r'wrote out restart data at nstep\s*=\s*(\d+)', log)
        complete = (proc.returncode == 0 and 'ENDRUN' not in log
                    and 'restFile_open Successfully' in log
                    and bool(steps) and int(steps[-1]) == manifest['expected_nstep']
                    and any(manifest['expected_end'][:10] in p.name for p in checkpoints))
        results[label] = {'completed': complete}
        if not complete:
            raise RuntimeError(f'{label}: incomplete simulation; see {folder}/model.log')
        desc = '\n'.join(p.read_text() for p in folder.glob('*.plant_trait.desc.*'))
        if label != 'legacy' and not re.search(r'HYPA.*allometry coefficient.*(?:100|25)', desc):
            raise AssertionError(f'{label}: resolved HYPA absent from trait description')
        balance = {}
        for hist in folder.glob('*.ecosim.h0.*.nc'):
            with Dataset(hist) as ds:
                for element in 'CNP':
                    name = f'Plant_BALANCE_{element}_pft'
                    values = np.ma.asarray(ds[name][:]).compressed()
                    if not np.isfinite(values).all():
                        raise AssertionError(f'{label}: nonfinite {name}')
                    balance[element] = max(balance.get(element, 0.), float(np.max(np.abs(values), initial=0.)))
        results[label]['maximum_absolute_plant_balance_g_per_m2'] = balance
        if any(value > 1e-7 for value in balance.values()):
            raise AssertionError(f'{label}: plant C/N/P balance exceeds 1e-7 g/m2')
        runs[label] = folder
    for label in ('legacy',):
        if label not in runs:
            continue
        diffs = {}
        for ref in list(runs['explicit100'].glob('*.ecosim.h0.*.nc')) + list(runs['explicit100'].glob('*.ecosim.r.*.nc')):
            diff = compare(ref, runs[label] / ref.name)
            if diff:
                diffs[ref.name] = diff
        results[label]['differences_from_explicit100'] = diffs
        if diffs:
            raise AssertionError(f'{label}: numerical differences from legacy-default path')
    def height(label):
        with Dataset(next(runs[label].glob('*.ecosim.h0.*.nc'))) as ds:
            result = np.asarray(ds['EmergHeight_pft'][:, 0])
            return np.where(np.asarray(ds['POPN_pft'][:, 0]) > 0, result, -np.inf)
    original, lower = height('explicit100'), height('lower25')
    changed = not np.array_equal(original, lower)
    results['lower25']['emergence_height_responds'] = changed
    # This requires an emergence window in the selected checkpoint/interval.
    first = lambda x: int(np.flatnonzero(x > -0.001)[0]) if np.any(x > -0.001) else None
    results['explicit100']['first_height_threshold_record'] = first(original)
    results['lower25']['first_height_threshold_record'] = first(lower)
    if first(original) is None or (first(lower) is not None and first(lower) <= first(original)):
        raise AssertionError('Lower HYPA did not delay the emergence-height threshold')
    if not changed:
        raise AssertionError('No emergence response: select a checkpoint/interval spanning planting')
    (args.workdir / 'summary.json').write_text(json.dumps(results, indent=2) + '\n')
    print(json.dumps(results, indent=2))


if __name__ == '__main__':
    main()
