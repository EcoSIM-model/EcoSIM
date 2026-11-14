"""Check completed EcoSIM runs for a plant carbon-balance jump at harvest.

The history must contain daily samples before and after the specified harvest.
Example: python tests/check_harvest_balance.py RUN --harvest-date 2022-10-04 --output check.json
"""
import argparse
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from netCDF4 import Dataset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run', type=Path)
    parser.add_argument('--harvest-date', required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    text = (args.run/'runfile.namelist').read_text()
    year = re.search(r"start_date\s*=\s*'(\d{4})", text)[1]
    bounds, values = [], {k: [] for k in ['C','N','P']}
    for path in sorted(args.run.glob('*.ecosim.h0.*.nc')):
        with Dataset(path) as ds:
            bounds.append(np.asarray(ds['time_bounds'][:]))
            for element in values:
                v = np.ma.asarray(ds['Plant_BALANCE_'+element+'_pft'][:])[:,0]
                if np.ma.getmaskarray(v).any() or not np.isfinite(v).all():
                    raise ValueError('Invalid active plant-balance history')
                values[element].append(np.asarray(v))
    b = np.concatenate(bounds)
    assert len(b) and np.all(b[:,1]-b[:,0] == 1) and np.array_equal(b[1:,0],b[:-1,1])
    frame = pd.DataFrame({k: np.concatenate(v) for k,v in values.items()},
                         index=pd.Timestamp(year+'-01-01')+pd.to_timedelta(b[:,0],unit='D'))
    harvest = pd.Timestamp(args.harvest_date)
    jumps = frame.loc[harvest+pd.Timedelta(days=1)]-frame.loc[harvest-pd.Timedelta(days=1)]
    maximum = frame.abs().max()
    log = (args.run/'run.log').read_text(errors='replace')
    stage = json.loads((args.run/'restart_run.json').read_text())
    steps = re.findall(r'nstep\s*=\s*(\d+)',log)
    checks = dict(completed=bool(steps and int(steps[-1])==stage['expected_nstep'] and 'ENDRUN' not in log),
                  carbon_residual_below_1e_6_g_per_m2=bool(maximum.C<1e-6),
                  all_elements_below_1e_6_g_per_m2=bool((maximum<1e-6).all()),
                  carbon_harvest_jump_below_1e_8_g_per_m2=bool(abs(jumps.C)<1e-8))
    result = dict(checks=checks,harvest_date=args.harvest_date,
                  max_absolute_cumulative_residual_g_per_m2=maximum.to_dict(),
                  cumulative_residual_change_across_harvest_g_per_m2=jumps.to_dict(),
                  daily_samples=len(frame))
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
