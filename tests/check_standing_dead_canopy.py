"""Exercise standing-dead canopy exchange in isolated full-model restart branches.

The supplied seed is copied before constructing controlled dead-only, mixed and
snow-buried initial states. These are thermal-path fixtures, not mortality
calibration runs. Original inputs, checkpoints and histories remain read-only.
Requires NumPy and netCDF4.
"""
import argparse
import json
from pathlib import Path
import shutil
import re
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

from check_annual_restart import ROOT


def set_slot(nc, name, value, slot=0):
    """Assign one named PFT on the variable's actual NetCDF axis order."""
    var = nc[name]
    index = [slice(None)] * var.ndim
    index[var.dimensions.index('pft')] = slot
    var[tuple(index)] = value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--seed', type=Path, required=True)
    parser.add_argument('--namelist', type=Path, required=True)
    parser.add_argument('--workdir', type=Path, required=True)
    args = parser.parse_args()
    exe = args.executable.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    parent_source = args.namelist.resolve()
    text = re.sub(r'(?m)![^\n]*$', '', parent_source.read_text())
    # Use one short forcing block for these controlled thermal experiments;
    # expanding the example's long spin-up cycles exceeds its namelist array.
    text = re.sub(r'(?m)^\s*forc_periods\s*=.*$', 'forc_periods=2013,2014,1', text)
    # Normalize the existing multiline history list before using the branch
    # stager, whose field-list replacement operates on a single assignment line.
    text = re.sub(r'(?ms)^\s*hist_fincl1\s*=.*?(?=^\s*\w+\s*=|^\s*/)',
                  "hist_fincl1='CAN_H_pft'\n", text)
    def absolute_input(match):
        candidate = parent_source.parent / match[1].rstrip()
        return repr(str(candidate.resolve())) if candidate.is_file() else match[0]
    text = re.sub(r"'([^'\n]+)'", absolute_input, text)
    parent = work/'parent.namelist'
    parent.write_text(text)
    report = {}
    for case in ('dead_only', 'mixed', 'buried'):
        seed_dir = work / (case + '-seed')
        seed_dir.mkdir()
        seed = seed_dir / args.seed.name
        shutil.copyfile(args.seed, seed)
        with Dataset(seed, 'r+') as nc:
            set_slot(nc, 'StandDeadCompKElms_pft', 0.)
            v = nc['StandDeadCompKElms_pft']
            for element, mass in enumerate((100., 1., .1)):
                idx = [0] * v.ndim
                idx[v.dimensions.index('elmnts')] = element
                v[tuple(idx)] = mass
            set_slot(nc, 'StandDeadStrutElms_pft', np.array([100., 1., .1]))
            set_slot(nc, 'PlantPopuDead_pft', .015)
            set_slot(nc, 'CanopyHeightDead_pft', .001 if case == 'buried' else 2.)
            if case != 'mixed':
                for name in ('PlantPopuLive_pft', 'IsPlantActive_pft', 'CanopyHeightLive_pft',
                             'CanopyLeafArea_pft', 'CanopyStemSurfArea_pft', 'LeafStalkAreaAct_pft',
                             'CanopyLeafAreaZ_pft', 'CanopyStemSurfAreaZ_pft',
                             'CanopyLeafArea_lnode', 'CanopyStalkSurfArea_lbrch',
                             'LeafAreaZsec_brch', 'StemAreaZsec_brch', 'CanopyBiomWater_pft'):
                    if name in nc.variables:
                        set_slot(nc, name, 0.)
                for name in ('CanopyLeafArea_col', 'StemArea_col', 'CanopyLeafAareZ_col', 'CanopyStemAareZ_col'):
                    if name in nc.variables:
                        nc[name][:] = 0.
            if case == 'buried':
                # A consistent one-layer dry-snow fixture: 0.08 m SWE at density
                # 0.1 gives 0.8 m bulk depth over the existing one-m2 column.
                for name in ('VLDrySnoWE', 'VLIceSnow', 'VLWatSnow', 'VLSnoDWIprev_snvr',
                             'SnowThickL_snvr', 'VLHeatCapSnow_snvr'):
                    nc[name][:] = 0.
                nc['VLDrySnoWE'][0, 0] = .08
                nc['SnoDensL'][:] = .1
                nc['SnowThickL_snvr'][0, 0] = .8
                nc['VLSnoDWIprev_snvr'][0, 0] = .8
                nc['VLHeatCapSnow_snvr'][0, 0] = 2.095*.08
                nc['TKSnow'][:] = 260.
                nc['TCSnow'][:] = 260.-273.15
                nc['SnowDepth'][:] = .8
                nc['VcumDrySnoWE_col'][:] = .08
                nc['VcumSnoDWI_col'][:] = .8
            assert nc['StandDeadStrutElms_pft'][0, 0] == 100.
        staged = work / (case + '-stage')
        cmd = [sys.executable, str(ROOT/'.agents/skills/ecosim-restart-run/scripts/stage_restart_run.py'),
               '--namelist', str(parent), '--restart', str(seed), '--out', str(staged),
               '--case-name', case, '--stop-n', '1', '--stop-option', 'nsteps',
               '--set', 'hist_nhtfrq=1', '--set', 'hist_mfilt=8760',
               '--add-fields', 'CAN_H_pft,CAN_LE_pft,CAN_RN_pft,CAN_G_pft,CAN_TEMPC_pft,QVegTransp_pft,STANDING_DEAD_C_pft',
               '--exe', str(exe)]
        with (seed_dir/'stage.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True)
        first, = staged.glob(case + '.ecosim.r.*.nc')
        target = work / case
        with (seed_dir/'comparison.log').open('w') as log:
            subprocess.run([sys.executable, str(ROOT/'tests/check_restart_exactness.py'),
                            '--executable', str(exe), '--namelist', str(staged/'runfile.namelist'),
                            '--seed', str(first), '--workdir', str(target), '--days', '2', '--split-hours', '24'],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        checkpoint = sorted((target/'continuous').glob(case + '.ecosim.r.*.nc'))[-1]
        with Dataset(checkpoint) as nc:
            if case != 'mixed':
                assert nc['PlantPopuLive_pft'][0] == 0.
                assert nc['Transpiration_pft'][0] == 0.
            assert nc['StandDeadStrutElms_pft'][0, 0] > 0.
            assert 100. < nc['TKC_pft'][0] < 400.
            assert 0. < nc['VHeatCapCanopy_pft'][0] < 10.
        history, = (target/'continuous').glob('*.ecosim.h0.*.nc')
        with Dataset(history) as nc:
            h = np.ma.asarray(nc['CAN_H_pft'][:, 0])
            assert not np.any(np.ma.getmaskarray(h)), 'dead canopy heat was hidden in history'
            if case == 'dead_only':
                assert np.max(np.abs(h)) > 1.e-6, 'dead material did not exchange sensible heat'
            if case == 'buried':
                assert np.all(h == 0.), 'buried material exchanged with air'
            if case != 'mixed':
                assert np.all(nc['QVegTransp_pft'][:, 0] == 0.)
            report[case] = {'maximum_absolute_sensible_W_m2': float(np.max(np.abs(h))), 'restart_differences': 0}
        print('PASS:', case, json.dumps(report[case]), flush=True)
    (work/'summary.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
