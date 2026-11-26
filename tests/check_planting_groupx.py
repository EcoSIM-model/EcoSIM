"""Exercise planting GROUPX overrides through the complete model input path.

Use a transient crop example and a warm January checkpoint for the next year.
All edited inputs, runs and diagnostics are isolated under a new --workdir.
Checks enabled/disabled flags, legacy records, annual reset, derived maturity,
and rejection of malformed/invalid overrides. No source inputs are modified.
"""
import argparse
import importlib.util
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from netCDF4 import Dataset, chartostring

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    'stage', ROOT/'.agents/skills/ecosim-restart-run/scripts/stage_restart_run.py')
stage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--executable', type=Path, required=True)
    ap.add_argument('--namelist', type=Path, required=True)
    ap.add_argument('--annual-restart', type=Path, required=True,
                    help='January 1 checkpoint at the next year, before annual management refresh')
    ap.add_argument('--workdir', type=Path, required=True)
    args = ap.parse_args()
    exe = args.executable.resolve(strict=True)
    nml = args.namelist.resolve(strict=True)
    restart = args.annual_restart.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    text = '\n'.join(stage.strip_inline_comment(line) for line in nml.read_text().splitlines())

    def absolute_input(match):
        p = nml.parent/match[1].rstrip()
        return "'"+str(p.resolve())+"'" if p.is_file() else match[0]

    text = re.sub(r"'([^'\n]+)'", absolute_input, text)
    pft = Path(stage.get_value(text, 'pft_mgmt_in').strip("'\""))
    year0 = int(stage.get_value(text, 'start_date').strip("'\"")[:4])
    with Dataset(pft) as ds:
        assert int(ds['pft_dflag'][...]) == 1, 'Use annual transient planting information'
        assert {year0, year0+1} <= set(ds['year'][:]), 'Need two consecutive years'
    for key, val, group in [
        ('case_name', "'groupx'", 'ecosim'), ('finidat', "''", 'ecosim'),
        ('continue_run', '.false.', 'ecosim'), ('disp_planttrait', '.true.', 'ecosim'),
        ('start_date', f"'{year0}0101000000'", 'ecosim'),
        ('forc_periods', f'{year0},{year0+1},1,{year0},{year0+1},0,{year0},{year0+1},0', 'ecosim'),
        ('hist_mfilt', '2000', 'ecosim'), ('hist_nhtfrq', '-24', 'ecosim'),
        ('stop_option', "'ndays'", 'ecosim_time'), ('stop_n', '20', 'ecosim_time'),
        ('rest_opt', "'nyears'", 'ecosim_time'), ('rest_frq', '1', 'ecosim_time'),
    ]:
        text = stage.set_value(text, key, val, group)
    results = {}

    def run(label, suffix, expected_error=None, annual=False):
        folder = work/label
        folder.mkdir()
        copy = folder/'planting.nc'
        shutil.copy2(pft, copy)
        with Dataset(copy, 'r+') as ds:
            var = ds['pft_pltinfo']
            records = chartostring(var[:])
            years = ds['year'][:]
            for iy, y in enumerate(years):
                for it in range(records.shape[1]):
                    for ip in range(records.shape[2]):
                        line = str(records[iy, it, ip]).strip()
                        if not line:
                            continue
                        extra = suffix(int(y)) if callable(suffix) else suffix
                        fields = re.split(r'[\s,]+', line)[:3]
                        # Explicit field count is part of the format under test.
                        records[iy, it, ip] = ' '.join(fields)+extra
            width = var.shape[-1]
            padded = np.asarray([str(s).ljust(width) for s in records.flat],
                                dtype=f'S{width}').reshape(records.shape)
            var[:] = padded.view('S1').reshape(*records.shape, width)
        source = stage.set_value(text, 'pft_mgmt_in', "'"+str(copy)+"'", 'ecosim')
        if annual:
            # Restore the previous year's state, then apply new annual management.
            # Avoid a cold start into a winter soil profile without spinup.
            source = stage.set_value(source, 'finidat', "'"+str(restart)+"'", 'ecosim')
            source = stage.set_value(source, 'start_date', f"'{year0+1}0101000000'", 'ecosim')
            source = stage.set_value(source, 'forc_periods',
                f'{year0+1},{year0+1},1,{year0+1},{year0+1},0,{year0+1},{year0+1},0', 'ecosim')
            source = stage.set_value(source, 'stop_n', '2', 'ecosim_time')
        (folder/'runfile.namelist').write_text(source)
        with (folder/'run.log').open('w') as log:
            proc = subprocess.run([str(exe), 'runfile.namelist'], cwd=folder,
                                  stdout=log, stderr=subprocess.STDOUT, timeout=600)
        body = (folder/'run.log').read_text(errors='replace')
        if expected_error:
            assert expected_error in body, (label, proc.returncode, body[-3000:])
            results[label] = 'rejected as expected'
            return
        assert proc.returncode == 0 and 'ENDRUN' not in body, (label, body[-3000:])
        steps = 48 if annual else 480
        assert re.search(r'nstep\s*=?\s*'+str(steps)+r'\b', body), (label, body[-3000:])
        values = {}
        for f in sorted(folder.glob('groupx.plant_trait.desc.*')):
            year = int(f.suffix[1:])
            desc = f.read_text()
            values[year] = [float(v) for v in re.findall(
                r'\bGROUPX\s*\|[^\n]*:\s*([^\s]+)', desc)]
            checkpoint = folder/f'groupx.ecosim.r.{year+1}-01-01-000000.nc'
            if not checkpoint.exists():
                checkpoint = folder/f'groupx.ecosim.r.{year}-01-03-000000.nc'
            if checkpoint.exists():
                xtli = [float(v) for v in re.findall(r'\bXTLI\s*\|[^\n]*:\s*([^\s]+)', desc)]
                ibtyp = [int(v) for v in re.findall(r'(?m)\bIBTYP:[^\n]*\s([0-9]+)\s*$', desc)]
                assert len(xtli) == len(ibtyp) == len(values[year])
                max_nodes = int(re.search(r'(?i)\bMaxNodesPerBranch\s*=\s*(\d+)',
                    (ROOT/'f90src/Mesh/GridConsts.F90').read_text())[1])
                expected = (np.array(values[year])-xtli)/np.where(np.array(ibtyp) == 0, 1, max_nodes)
                with Dataset(checkpoint) as state:
                    # Restart storage reserves JP slots, including inactive fill.
                    np.testing.assert_allclose(state['MatureGroup_pft'][:len(expected)],
                                               expected, rtol=0, atol=1.e-12)
        assert values and all(values.values()), (label, 'missing GROUPX diagnostic')
        results[label] = values
        print(label, values, flush=True)
        return values

    baseline = run('legacy3', '')[year0]
    assert run('legacy4', ' 999')[year0] == baseline
    assert run('flag0', ' 0 99')[year0] == baseline
    assert run('flag1', ' 1 13.25')[year0] == [13.25]*len(baseline)
    annual = run('annual', lambda y: ' 1 13.25' if y == year0 else ' 1 14.5', annual=True)
    assert annual[year0+1] == [14.5]*len(annual[year0+1])
    # The crop may rotate, so compare fallback against a legacy run in that year.
    legacy = run('annual_legacy', '', annual=True)
    reset = run('annual_reset_flag0', ' 0 99', annual=True)
    assert reset[year0+1] == legacy[year0+1]
    for label, suffix, error in [
        ('bad_flag', ' 2 12', 'Planting GROUPX flag must be 0 or 1'),
        ('negative', ' 1 -1', 'Planting GROUPX must be finite and nonnegative'),
        ('nan', ' 1 NaN', 'Planting GROUPX must be finite and nonnegative'),
        ('infinite', ' 1 Infinity', 'Planting GROUPX must be finite and nonnegative'),
        ('missing', ' 1 twelve', 'Invalid planting record'),
        ('extra', ' 1 12 99', 'Plant records require explicit values'),
    ]:
        run(label, suffix, error)
    (work/'results.json').write_text(json.dumps(results, indent=2)+'\n')
    print('PASS: legacy forms, flags, yearly overrides/reset and invalid planting GROUPX')


if __name__ == '__main__':
    main()
