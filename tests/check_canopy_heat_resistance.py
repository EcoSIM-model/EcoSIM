"""Exercise canopy heat/H2O resistance mapping, history and legacy restart reads.

Run an active prescribed canopy continuously for 48 hours and as 24+24 hours.
The resumed variants cover current, missing heat, and legacy water cache names.
All inputs and outputs are isolated under a new --workdir; requires netCDF4/numpy.
"""
import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from netCDF4 import Dataset

from check_annual_restart import ROOT, stage
from check_restart_exactness import compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', required=True, type=Path)
    parser.add_argument('--workdir', required=True, type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = ROOT / 'examples/run_dir/SatePhenol/SatePhenol_test.namelist'
    text = re.sub(r'(?m)![^\n]*$', '', source.read_text())

    def absolute_input(match):
        candidate = source.parent / match[1].rstrip()
        return "'" + str(candidate.resolve()) + "'" if candidate.is_file() else match[0]

    text = re.sub(r"'([^'\n]+)'", absolute_input, text)
    phenology = work / 'phenology.nc'
    shutil.copyfile(ROOT / 'examples/inputs/blodgett/Blodget_sp_phenol_20261103.nc', phenology)
    with Dataset(phenology, 'r+') as nc:
        for field, value in [('LAI', 1.2), ('SAI', .3), ('HTOP', 17.)]:
            nc[field][:] = 0.
            nc[field][0, 0, 0, :] = value
    for key, value, group in [
        ('case_name', "'heat-resistance'", 'ecosim'),
        ('sp_phenol_in', "'" + str(phenology) + "'", 'ecosim'),
        ('pft_mgmt_in', "'NO'", 'ecosim'),
        ('forc_periods', '2013,2014,1', 'ecosim'),
        ('finidat', "''", 'ecosim'),
        ('hist_nhtfrq', '1', 'ecosim'), ('hist_mfilt', '8760', 'ecosim'),
        ('hist_fincl1', "'BLYR_RSC_H2O_pft','BLYR_RSC_HEAT_pft','CAN_H_pft'", 'ecosim'),
        ('stop_option', "'nsteps'", 'ecosim_time'),
        ('rest_opt', "'never'", 'ecosim_time')]:
        text = stage.set_value(text, key, value, group)

    def run(name, hours, seed=None, legacy=None):
        folder = work / name
        folder.mkdir()
        if seed:
            for pattern in ['*.ecosim.r.*.nc', '*.ecosim.rh*.nc', '*.ecosim.h*.nc']:
                for path in seed.parent.glob(pattern):
                    shutil.copy2(path, folder / path.name)
            checkpoint = folder / seed.name
            if legacy:
                with Dataset(checkpoint, 'r+') as nc:
                    nc.renameVariable('RawHeatCanopy2Atm_pft', 'unused_heat_cache')
                    if legacy == 'old_name':
                        nc.renameVariable('RawH2OCanopy2Atm_pft', 'RawCanopy2Atm_pft')
            (folder / 'rpointer.esim').write_text(seed.name + '\n')
        case = stage.set_value(text, 'stop_n', str(hours), 'ecosim_time')
        case = stage.set_value(case, 'continue_run', '.true.' if seed else '.false.', 'ecosim')
        (folder / 'case.namelist').write_text(case)
        with (folder / 'run.log').open('w') as log:
            subprocess.run([str(executable), 'case.namelist'], cwd=folder,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        log = (folder / 'run.log').read_text(errors='replace')
        steps = re.findall(r'wrote out restart data at nstep\s*=\s*(\d+)', log)
        assert 'ENDRUN' not in log and steps and int(steps[-1]) == hours, folder
        checkpoints = sorted(folder.glob('*.ecosim.r.*.nc'))
        checkpoint = checkpoints[-1]
        assert not compare(checkpoint, checkpoint), folder
        with Dataset(checkpoint) as nc:
            water = nc['RawH2OCanopy2Atm_pft'][:]
            heat = nc['RawHeatCanopy2Atm_pft'][:]
            np.testing.assert_array_equal(np.ma.getmaskarray(water), np.ma.getmaskarray(heat))
            np.testing.assert_allclose(heat.compressed(), water.compressed() / .92, rtol=1.e-14)
            assert np.any(heat.compressed() > 0.), folder
        history, = folder.glob('*.ecosim.h0.*.nc')
        assert not compare(history, history), folder
        with Dataset(history) as nc:
            heat = nc['BLYR_RSC_HEAT_pft'][:]
            water = nc['BLYR_RSC_H2O_pft'][:]
            assert nc['BLYR_RSC_HEAT_pft'].units == 's/m'
            np.testing.assert_array_equal(np.ma.getmaskarray(water), np.ma.getmaskarray(heat))
            np.testing.assert_allclose(heat.compressed(), water.compressed() / .92, rtol=2.e-7)
            assert np.any(np.abs(nc['CAN_H_pft'][:].compressed()) > 1.e-6), folder
        return checkpoint, history

    reference = run('continuous', 48)
    first = run('first_segment', 24)
    results = {}
    for name, legacy in [('current', None), ('missing_heat', 'missing_heat'), ('old_name', 'old_name')]:
        resumed = run(name, 48, first[0], legacy)
        differences = {'checkpoint': compare(reference[0], resumed[0]),
                       'history': compare(reference[1], resumed[1])}
        results[name] = differences
    (work / 'comparison.json').write_text(json.dumps(results, indent=2) + '\n')
    assert all(not differences for result in results.values() for differences in result.values()), results
    print('PASS: active canopy state/history ratios; exact 24+24 hour continuation for all three restart schemas')
    print(work / 'comparison.json')


if __name__ == '__main__':
    main()
