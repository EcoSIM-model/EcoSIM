"""Check that canopy layers represent every emerged leaf node's C/N/P and area.

Use an emerged, live canopy checkpoint before harvest. A missing layer projection
can silently destroy node biomass when harvest reconstructs the remaining pools.
Example: python tests/check_leaf_layer_closure.py CHECKPOINT --output result.json
"""
import argparse
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset


def check(path, pft):
    """Compare authoritative node pools with their summed spatial projections."""
    results = {}
    with Dataset(path) as ds:
        for node_name, layer_name in [('LeafElmntNode_brch', 'LeafLayerElms_node'),
                                      ('ARLF', 'CanopyLeafArea_lnode')]:
            node_var, layer_var = ds[node_name], ds[layer_name]
            node_dims = list(node_var.dimensions)
            layer_dims = list(layer_var.dimensions)
            node = np.ma.asarray(node_var[:]).take(pft, axis=node_dims.index('pft'))
            layer = np.ma.asarray(layer_var[:]).take(pft, axis=layer_dims.index('pft'))
            node_dims.remove('pft')
            layer_dims.remove('pft')
            layer = layer.sum(axis=layer_dims.index('levcan'))
            layer_dims.remove('levcan')
            layer = layer.transpose([layer_dims.index(d) for d in node_dims])
            if np.ma.getmaskarray(node).any() or np.ma.getmaskarray(layer).any():
                raise ValueError('Select an active PFT with fully populated node fields')
            if not np.isfinite(node).all() or not np.isfinite(layer).all():
                raise ValueError('Nonfinite leaf pool or layer projection')
            delta = np.asarray(node-layer)
            bad = np.abs(delta) > 1e-9 + 1e-12*np.abs(node)
            entries = []
            for index in np.argwhere(bad):
                ix = tuple(index)
                entries.append(dict(index=dict(zip(node_dims, map(int,index))),
                                    node=float(node[ix]), layer_sum=float(layer[ix]),
                                    difference=float(delta[ix])))
            results[node_name] = dict(max_absolute_difference=float(np.max(np.abs(delta))),
                                      total_difference=float(delta.sum()), mismatches=entries)
    return dict(checkpoint=str(path), pft=pft, passed=not any(x['mismatches'] for x in results.values()),
                comparisons=results)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoint', type=Path)
    parser.add_argument('--pft', type=int, default=0, help='zero-based active PFT slot')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    report = check(args.checkpoint, args.pft)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(('PASS' if report['passed'] else 'FAIL')+': '+str(args.output))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
