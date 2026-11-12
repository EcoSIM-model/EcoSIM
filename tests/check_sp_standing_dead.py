"""Check completed output from the SatePhenolStandingDead example.

Requires NumPy and netCDF4. Reads outputs without changing them. This checks
passive exchange, zero living/root activity, and the represented heat inventory;
restart equality is checked separately by the ecosim-restart-test runner.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

from check_restart_exactness import compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir', type=Path, required=True)
    parser.add_argument('--hours', type=int, default=480)
    args = parser.parse_args()
    folder = args.run_dir
    checkpoint = sorted(folder.glob('*.ecosim.r.*.nc'))[-1]
    assert not compare(checkpoint, checkpoint), 'nonfinite checkpoint values'
    with Dataset(checkpoint) as nc:
        assert nc['NP'][:].item() == nc['NP0'][:].item() == 1
        for field in ('PlantPopuLive_pft', 'CanopyLeafArea_pft', 'CanopyStemSurfArea_pft',
                      'CanopyLeafSheathC_pft', 'CanopySapwoodC_pft', 'CanopyBiomWater_pft',
                      'ShootElms_pft', 'PopuRootMycoC_pvr', 'RootAbsorbLenPerPlant_pvr',
                      'RootTotLenPerPlant_pvr', 'Transpiration_pft'):
            assert np.all(nc[field][:] == 0.), field
        assert nc['CanopyHeightDead_pft'][0] == 17.
        assert nc['CanopyHeight_col'][:].item() == 17.
        capacity = float(nc['VHeatCapCanopy_pft'][0])
        assert capacity > 0., 'missing native canopy heat capacity'
        # The only PFT is passive, so the column inventory must be C*T.
        energy = float(nc['CanopyHeatStor_col'][:].item())
        np.testing.assert_allclose(energy, capacity*float(nc['TKC_pft'][0]), rtol=1.e-12)
    histories = sorted(folder.glob('*.ecosim.h0.*.nc'))
    samples = 0
    values = {}
    for history in histories:
        assert not compare(history, history), 'nonfinite history values'
        with Dataset(history) as nc:
            samples += len(nc.dimensions['time'])
            for field in ('CAN_H_pft', 'CAN_LE_pft', 'CAN_RN_pft', 'CAN_G_pft',
                          'CAN_TEMPC_pft', 'QVegTransp_pft', 'STANDING_DEAD_C_pft'):
                active = np.ma.asarray(nc[field][:, 0])
                assert not np.any(np.ma.getmaskarray(active)), field
                values.setdefault(field, []).extend(active.tolist())
            roots = np.ma.asarray(nc['RootC_2nd_pvr'][:, :, 0])
            assert not np.any(np.ma.getmaskarray(roots)), 'missing active root history'
            assert np.all(roots == 0.), 'dead canopy acquired live fine roots'
    assert samples == args.hours, (samples, args.hours)
    for field in ('CAN_H_pft', 'CAN_LE_pft', 'CAN_RN_pft', 'CAN_G_pft'):
        assert np.max(np.abs(values[field])) > 1.e-6, field
    heat = np.asarray([values[field] for field in ('CAN_RN_pft', 'CAN_H_pft', 'CAN_LE_pft', 'CAN_G_pft')])
    tolerance = 8.*np.finfo(np.float32).eps*np.maximum(1., np.abs(heat).sum(axis=0))
    assert np.all(np.abs(heat.sum(axis=0)) <= tolerance), 'passive canopy energy does not close'
    assert np.all(np.asarray(values['QVegTransp_pft']) == 0.)
    np.testing.assert_allclose(values['STANDING_DEAD_C_pft'], 5000., rtol=0., atol=0.)
    result = {'hours': samples, 'living_population_and_roots': 'zero',
              'transpiration': 'zero', 'standing_dead_C_gC_m2': 5000.,
              'canopy_capacity_MJ_K': capacity,
              'sensible_W_m2_range': [min(values['CAN_H_pft']), max(values['CAN_H_pft'])],
              'latent_W_m2_range': [min(values['CAN_LE_pft']), max(values['CAN_LE_pft'])]}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
