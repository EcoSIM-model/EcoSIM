"""Check paired, live-crop XRLA runs for realized leaf-node progression.

Supply run directories from a common checkpoint with identical inputs except
XRLA. End both runs during active vegetative growth, before node-slot wrapping
or harvest. This checks realized biomass and LAI rather than just leaf counters.
Example: python tests/check_leaf_node_progression.py --normal RUN002 --slower RUN001 --output result.json
"""
import argparse
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def read_run(folder):
    """Read active PFT 0, main branch and paired daily histories."""
    checkpoints = sorted(folder.glob('*.ecosim.r.*.nc'),
                         key=lambda path: path.name.split('.ecosim.r.')[1])
    if not checkpoints:
        raise ValueError(f'No output checkpoint in {folder}')
    checkpoint = checkpoints[-1]
    with Dataset(checkpoint) as ds:
        def at_pft(name):
            var = ds[name]
            return np.ma.asarray(var[:]).take(0, axis=var.dimensions.index('pft'))
        branch = int(at_pft('MainBranchNum_pft')) - 1
        if branch < 0 or branch >= len(at_pft('KVSTG')):
            raise ValueError(f'No valid live main branch in {checkpoint}')
        frontier = int(at_pft('KVSTG')[branch])
        carbon = at_pft('LeafElmntNode_brch')[branch, :, 0]
        occupied = np.flatnonzero(np.ma.filled(carbon > 1e-6, False)).tolist()
        leaves = float(at_pft('NumOfLeaves_brch')[branch])
    values = {}
    bounds = []
    for path in sorted(folder.glob('*.ecosim.h0.*.nc')):
        with Dataset(path) as ds:
            bounds.append(np.asarray(ds['time_bounds'][:]))
            for name in ['ECO_LAI_col', 'Plant_BALANCE_C_pft',
                         'Plant_BALANCE_N_pft', 'Plant_BALANCE_P_pft']:
                v = np.ma.asarray(ds[name][:])[:, 0]
                if np.ma.getmaskarray(v).any() or not np.isfinite(v).all():
                    raise ValueError(f'Missing or nonfinite {name} in {path}')
                values.setdefault(name, []).append(np.asarray(v))
    if not bounds:
        raise ValueError(f'No daily history in {folder}')
    values = {name: np.concatenate(parts) for name, parts in values.items()}
    balance = {name: float(np.max(np.abs(v))) for name, v in values.items()
               if name.startswith('Plant_BALANCE')}
    return dict(checkpoint=str(checkpoint), frontier=frontier,
                occupied_leaf_nodes=occupied, main_branch_appeared_leaves=leaves,
                max_absolute_balance_residual=balance), values['ECO_LAI_col'], np.concatenate(bounds)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--normal', type=Path, required=True)
    parser.add_argument('--slower', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    normal, lai_n, bounds_n = read_run(args.normal)
    slower, lai_s, bounds_s = read_run(args.slower)
    checks = {
        'same_checkpoint_date': normal['checkpoint'].split('.ecosim.r.')[1] == slower['checkpoint'].split('.ecosim.r.')[1],
        'same_daily_intervals': np.array_equal(bounds_n, bounds_s),
        'normal_has_new_leaf_nodes': normal['frontier'] > 1 and any(k > 1 for k in normal['occupied_leaf_nodes']),
        'slower_has_new_leaf_nodes': slower['frontier'] > 1 and any(k > 1 for k in slower['occupied_leaf_nodes']),
        'slower_appearance_changes_actual_frontier': slower['frontier'] < normal['frontier'],
        'slower_appearance_changes_leaf_count': slower['main_branch_appeared_leaves'] < normal['main_branch_appeared_leaves'],
        'leaf_area_responds': lai_n.shape == lai_s.shape and not np.array_equal(lai_n, lai_s),
        'plant_balance_within_1e_7_g_per_m2': all(v < 1e-7 for run in [normal, slower]
                                               for v in run['max_absolute_balance_residual'].values()),
    }
    result = dict(normal=normal, slower=slower, checks=checks,
                  max_absolute_LAI_difference=float(np.max(np.abs(lai_n-lai_s))))
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
    return 0 if all(checks.values()) else 1


if __name__ == '__main__':
    raise SystemExit(main())
