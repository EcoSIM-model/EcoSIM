"""Guard against false passes in the numerical restart regression harness."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

import numpy as np
from netCDF4 import Dataset

spec = importlib.util.spec_from_file_location('restart_runner', Path(__file__).with_name('run_restart_repro.py'))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class RestartComparisonTests(unittest.TestCase):
    def setUp(self):
        self.scratch = tempfile.TemporaryDirectory()
        self.addCleanup(self.scratch.cleanup)
        self.root = Path(self.scratch.name)

    def make_file(self, path):
        with Dataset(path, 'w', format='NETCDF3_64BIT_OFFSET') as d:
            d.createDimension('time', None)
            d.createDimension('nbnd', 2)
            b = d.createVariable('time_bounds', 'f8', ('time', 'nbnd'))
            b[:] = [[0, 1], [1, 2], [2, 3]]
            v = d.createVariable('C', 'f8', ('time',), fill_value=1.e36)
            v.units = 'g C m-2'
            v[:] = [1., 2., 1.e36]

    def test_perturbation_is_detected_and_fill_is_allowed(self):
        a, b = self.root/'a.nc', self.root/'b.nc'
        self.make_file(a)
        self.make_file(b)
        self.assertEqual(runner.compare_numeric(a, b)['different_fields'], [])
        with Dataset(b, 'a') as d:
            d['C'][0] = 1.01
        self.assertEqual(runner.compare_numeric(a, b)['different_fields'], ['C'])

    def test_identical_real_nans_do_not_pass(self):
        a, b = self.root/'a.nc', self.root/'b.nc'
        self.make_file(a)
        self.make_file(b)
        for p in [a, b]:
            with Dataset(p, 'a') as d:
                d['C'][0] = np.nan
        self.assertEqual(runner.compare_numeric(a, b)['different_fields'], ['C'])

    def test_large_integer_difference_is_not_rounded_away(self):
        a, b = self.root/'a.nc', self.root/'b.nc'
        for path, value in [(a, 2**53), (b, 2**53 + 1)]:
            with Dataset(path, 'w') as d:
                v = d.createVariable('counter', 'i8')
                v.assignValue(value)
        self.assertEqual(runner.compare_numeric(a, b)['different_fields'], ['counter'])

    def test_intraday_copy_excludes_unfinished_and_future_records(self):
        a, b = self.root/'a.nc', self.root/'b.nc'
        self.make_file(a)
        runner.copy_history_prefix(a, b, 2.5)
        with Dataset(b) as d:
            self.assertEqual(len(d.dimensions['time']), 2)
            np.testing.assert_array_equal(d['C'][:], [1., 2.])
            self.assertEqual(d['C'].units, 'g C m-2')
            self.assertEqual(d['C']._FillValue, 1.e36)


if __name__ == '__main__':
    unittest.main()
