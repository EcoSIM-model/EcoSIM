#!/usr/bin/env python3
"""Check emitted soil+litter+snow exchange against total ground evaporation.

Accept history files or directories. The tolerance scales with the stored
floating-point precision and component magnitudes, including cancellation.
This checks history accounting, not the physical accuracy of evaporation.
"""
import argparse
import json
from pathlib import Path
import re

import numpy as np
from netCDF4 import Dataset

FIELDS = ('EVAPSoil_col', 'EVAPLitr_col', 'EVAPSnow_col', 'EVAPGrnd_col')


def check_file(path):
    """Validate units, masks, finiteness and the component identity in one tape."""
    with Dataset(path) as dataset:
        missing = [name for name in FIELDS if name not in dataset.variables]
        if missing:
            raise ValueError(f'{path}: missing fields {missing}')
        variables = [dataset.variables[name] for name in FIELDS]
        if any(v.dimensions != variables[0].dimensions for v in variables):
            raise ValueError(f'{path}: component dimensions differ')
        if any(v.units != 'mm H2O/hr' for v in variables):
            raise ValueError(f'{path}: component units differ')
        arrays = [np.ma.asarray(v[:]) for v in variables]
        masks = [np.ma.getmaskarray(a) for a in arrays]
        if any(not np.array_equal(mask, masks[0]) for mask in masks):
            raise ValueError(f'{path}: component masks differ')
        values = [np.asarray(a.compressed(), dtype=np.float64) for a in arrays]
        if not values[0].size:
            raise ValueError(f'{path}: no valid exchange samples')
        if any(not np.isfinite(v).all() for v in values):
            raise ValueError(f'{path}: nonfinite exchange')
        error = np.abs(values[0]+values[1]+values[2]-values[3])
        eps = max(np.finfo(v.dtype).eps for v in variables)
        scale = sum(np.abs(v) for v in values)
        bound = 8*eps*scale + np.finfo(np.float64).tiny
        if np.any(error > bound):
            raise ValueError(f'{path}: partition mismatch; max error {error.max()} mm/hr')
        return {'file': str(path.resolve()), 'samples': len(error),
                'max_abs_error_mm_hr': float(error.max()),
                'max_precision_scaled_error': float(np.max(error/bound)),
                'nonzero_samples': {name: int(np.count_nonzero(v)) for name, v in zip(FIELDS, values)}}


def main():
    """Check selected tapes and optionally save a machine-readable summary."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('paths', nargs='+', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    files = set()
    for path in args.paths:
        candidates = path.rglob('*.nc') if path.is_dir() else [path]
        files.update(p for p in candidates if re.search(r'\.ecosim\.h\d+\.', p.name))
    if not files:
        parser.error('no history tapes selected')
    records = [check_file(path) for path in sorted(files)]
    result = {'status': 'passed', 'files_checked': len(records), 'files': records}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': 'passed', 'files_checked': len(records),
                      'max_abs_error_mm_hr': max(x['max_abs_error_mm_hr'] for x in records)}, indent=2))


if __name__ == '__main__':
    main()
