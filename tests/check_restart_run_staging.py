"""Regression checks for staging a continuation from a real checkpoint layout."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from datetime import datetime

import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / '.agents/skills/ecosim-restart-run/scripts/stage_restart_run.py'
spec = importlib.util.spec_from_file_location('restart_staging', SCRIPT)
staging = importlib.util.module_from_spec(spec)
spec.loader.exec_module(staging)


class RestartStagingTest(unittest.TestCase):
    def stage(self, mode, seconds=0, extra=()):
        """Stage a checkpoint whose history names exist only in the master file."""
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        parent = Path(tmp.name) / 'parent'
        parent.mkdir()
        namelist = parent / 'lake.namelist'
        namelist.write_text("""&ecosim
case_name='Lake'
start_date='18000101000000'
forc_periods=1980,1989,2
npxs=30 !water/heat substeps
npys=10 !gas substeps
sp_phenol_in = 'phenology.nc'
/
&ecosim_time
delta_time=3600.
/
""")
        (parent / 'phenology.nc').write_bytes(b'phenology')
        original = namelist.read_bytes()
        stamp = f'1800-07-30-{seconds:06d}.nc'
        restart = parent / ('Lake.ecosim.r.' + stamp)
        history_name = 'Lake.ecosim.h0.1800-01-01-00000.nc'
        with Dataset(restart, 'w') as nc:
            nc.createDimension('ntapes', 1)
            nc.createDimension('max_chars', 256)
            v = nc.createVariable('locfnh', 'S1', ('ntapes', 'max_chars'))
            v[0, :] = np.frombuffer(('./' + history_name).ljust(256).encode(), dtype='S1')
        rh = parent / ('Lake.ecosim.rh0.' + stamp)
        with Dataset(rh, 'w') as nc:
            nc.createDimension('scalar', 1)
            nc.createVariable('ntimes', 'i4', ('scalar',))[:] = 7
        (parent / history_name).write_bytes(b'preserved history content')
        out = Path(tmp.name) / 'staged'
        args = [sys.executable, str(SCRIPT), '--namelist', str(namelist),
                '--restart', str(restart), '--out', str(out), '--mode', mode,
                '--stop-n', '10', '--stop-option', 'ndays', *extra]
        result = subprocess.run(args, capture_output=True, text=True)
        self.assertEqual(namelist.read_bytes(), original)
        return parent, out, result

    def test_continue_copies_master_named_history_and_adds_duration(self):
        parent, out, result = self.stage('continue')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((out / 'Lake.ecosim.h0.1800-01-01-00000.nc').read_bytes(),
                         b'preserved history content')
        record = json.loads((out / 'restart_run.json').read_text())
        self.assertEqual(record['expected_nstep'], 5280)
        self.assertEqual(record['expected_end'], '1800-08-09T00:00:00')
        self.assertIn("stop_option = 'nsteps'", (out / 'runfile.namelist').read_text())
        self.assertTrue((parent / 'Lake.ecosim.r.1800-07-30-000000.nc').exists())

    def test_branch_retains_intraday_checkpoint_time(self):
        _, out, result = self.stage('branch', seconds=18000)
        self.assertEqual(result.returncode, 0, result.stderr)
        text = (out / 'runfile.namelist').read_text()
        self.assertIn("start_date = '18000730050000'", text)
        record = json.loads((out / 'restart_run.json').read_text())
        self.assertEqual(record['expected_nstep'], 5285)
        self.assertEqual(record['expected_end'], '1800-08-09T05:00:00')

    def test_phenology_input_is_resolved_and_copied_when_requested(self):
        parent, out, result = self.stage('branch')
        self.assertEqual(result.returncode, 0, result.stderr)
        text = (out / 'runfile.namelist').read_text()
        self.assertIn(f"sp_phenol_in = '{parent / 'phenology.nc'}'", text)

        _, copied_out, result = self.stage('branch', extra=('--copy-input', 'sp_phenol_in'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((copied_out / 'phenology.nc').read_bytes(), b'phenology')
        self.assertIn("sp_phenol_in = './phenology.nc'", (copied_out / 'runfile.namelist').read_text())

    def test_calendar_endpoints(self):
        self.assertEqual(staging.run_end_time(datetime(2000, 2, 29, 5), 1, 'nyears', 3600),
                         datetime(2001, 2, 28, 5))
        self.assertEqual(staging.run_end_time(datetime(2000, 1, 31, 5), 1, 'nmonths', 3600),
                         datetime(2000, 2, 29, 5))

    def test_comment_markers_inside_strings_are_preserved(self):
        self.assertEqual(staging.get_value("file='a''!b.nc' !comment", 'file'), "'a''!b.nc'")

    def test_short_run_is_rejected_even_with_success_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            fake = Path(tmp) / 'fake-model'
            fake.write_text('#!/bin/sh\necho "restFile_open Successfully"\n'
                            'echo "Successfully wrote out restart data at nstep = 1"\n')
            fake.chmod(0o755)
            _, _, result = self.stage('continue', extra=('--exe', str(fake)))
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('completed=False', result.stdout)


if __name__ == '__main__':
    unittest.main()
