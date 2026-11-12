"""Check current reference-air exchange in the existing bare-soil example.

Runs in a new scratch directory; requires numpy/netCDF4. Default: 10 days.
Use --days 365 to include snow-free diurnal cycles. No source inputs are edited.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from netCDF4 import Dataset

from check_annual_restart import ROOT, stage
from check_restart_exactness import compare


def validate(folder, hours):
    """Check endpoints, finite state, inactive feedback and heat-flux direction."""
    log = (folder / 'run.log').read_text(errors='replace')
    steps = re.findall(r'wrote out restart data at nstep\s*=\s*(\d+)', log)
    assert steps and int(steps[-1]) == hours, folder
    history, = folder.glob('*.ecosim.h0.*.nc')
    assert not compare(history, history), history
    with Dataset(history) as dataset:
        def column(name):
            return np.asarray(dataset[name][:, 0], dtype=float)
        air = column('AIR_TEMP_col')
        endpoint = column('TKQ_col') - 273.15
        soil = np.asarray(dataset['TEMP_vr'][:, 0, 0], dtype=float)
        sensible = column('ECO_HeatS_col')  # positive into soil
        assert len(air) == hours
        np.testing.assert_allclose(endpoint, air, rtol=0., atol=2.e-5)
        np.testing.assert_array_equal(column('VPQ_col'), column('HUM_col'))
        bare = (column('Frac_Snow_Ground_col') < 1.e-6) & (column('Frac_Litr_Ground_col') < 1.e-6)
        previous_soil = np.roll(soil, 1)
        consecutive = bare & np.roll(bare, 1)
        consecutive[0] = False
        # Hourly H and end temperatures need not agree during a gradient
        # reversal within the hour. Screen only separated, same-side endpoints.
        wrong = consecutive & (
            ((np.maximum(soil, previous_soil) < air - 1.) & (sensible < -10.))
            | ((np.minimum(soil, previous_soil) > air + 1.) & (sensible > 10.)))
        assert not np.any(wrong), np.flatnonzero(wrong).tolist()
        report = {'hours': hours, 'snow_litter_free_hours': int(np.count_nonzero(bare)),
                  'opposite_gradient_hours': int(np.count_nonzero(wrong)),
                  'max_air_endpoint_difference_K': float(np.max(np.abs(endpoint-air)))}
    checkpoint = sorted(folder.glob('*.ecosim.r.*.nc'))[-1]
    assert not compare(checkpoint, checkpoint), checkpoint
    with Dataset(checkpoint) as dataset:
        assert np.all(dataset['LeafStalkAreaAll_col'][:] == 0.)
        for name in ('Air_Heat_Sens_store_col', 'Air_Heat_Latent_store_col',
                     'dSensHeatdTKQ_col', 'dLatHeatdVPQ_col'):
            assert np.all(dataset[name][:] == 0.), name
    (folder / 'surface-exchange-check.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', required=True, type=Path)
    parser.add_argument('--workdir', required=True, type=Path)
    parser.add_argument('--days', default=10, type=int)
    args = parser.parse_args()
    if args.days <= 0:
        parser.error('--days must be positive')
    executable = args.executable.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = ROOT / 'examples/run_dir/bare_soil/BareSoil.namelist'
    text = re.sub(r'(?m)![^\n]*$', '', source.read_text())
    def absolute_input(match):
        candidate = source.parent / match[1].rstrip()
        return repr(str(candidate.resolve())) if candidate.is_file() else match[0]
    text = re.sub(r"'([^'\n]+)'", absolute_input, text)
    for key, value, group in (
        ('case_name', "'bare-exchange'", 'ecosim'),
        ('hist_fincl1', "'TKQ_col','VPQ_col','Frac_Snow_Ground_col','Frac_Litr_Ground_col'", 'ecosim'),
        ('hist_nhtfrq', '1', 'ecosim'), ('hist_mfilt', str(args.days*24), 'ecosim'),
        ('continue_run', '.false.', 'ecosim'), ('finidat', "''", 'ecosim'),
        ('stop_option', "'nsteps'", 'ecosim_time'), ('stop_n', str(args.days*24), 'ecosim_time'),
        ('rest_opt', "'never'", 'ecosim_time'),
    ):
        text = stage.set_value(text, key, value, group)
    (work / 'case.namelist').write_text(text)
    with (work / 'run.log').open('w') as log:
        subprocess.run([str(executable), 'case.namelist'], cwd=work, stdout=log,
                       stderr=subprocess.STDOUT, check=True)
    print('PASS:', json.dumps(validate(work, args.days*24)))


if __name__ == '__main__':
    main()
