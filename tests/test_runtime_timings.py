"""Exercise the real Fortran runtime profiler without model inputs."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
SOURCE = '''program test_timings
  use iso_fortran_env, only: int64
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use timings
  implicit none
  type(timer_stamp_type) :: outer, inner, synthetic
  integer :: i
  character(len=32) :: name
  character(len=512) :: destination
  call get_command_argument(1,destination)
  synthetic = timer_stamp_type(98_int64,10_int64,100_int64)
  if (abs(timer_elapsed(synthetic,2_int64)-0.5_r8)>1.e-12_r8) stop 1
  synthetic = timer_stamp_type(10_int64,10_int64,huge(0_int64))
  if (abs(timer_elapsed(synthetic,15_int64)-0.5_r8)>1.e-12_r8) stop 2
  synthetic = timer_stamp_type(huge(0_int64)-1_int64,10_int64,huge(0_int64))
  if (abs(timer_elapsed(synthetic,1_int64)-0.3_r8)>1.e-12_r8) stop 3
  call init_timer(trim(destination),.true.,'20261008-120000')
  call start_timer(inner)
  call end_timer('setup',inner)
  call flush_timer_detail()
  call start_timer(outer)
  call wait_ticks(0.02_r8)
  call start_timer(inner)
  call wait_ticks(0.01_r8)
  call end_timer('inner',inner)
  call wait_ticks(0.02_r8)
  call end_timer('outer',outer)
  call end_timer_loop()
  do i=1,2
    call start_timer(inner)
    call end_timer('inner',inner)
  enddo
  do i=1,220
    write(name,'(A,I0)') 'region_',i
    call start_timer(inner)
    call end_timer(trim(name),inner)
  enddo
  call start_timer(inner)
  call end_timer('a very long label, with a "quoted" name',inner)
  call end_timer_loop()
  call finalize_timer()
  call finalize_timer()
  if (len(timer_launch_label())/=15) stop 4
  call init_timer(trim(destination)//'/second')
  call start_timer(inner)
  call end_timer('new_run',inner)
  call end_timer_loop()
  call finalize_timer()
contains
  subroutine wait_ticks(seconds)
    real(r8), intent(in) :: seconds
    type(timer_stamp_type) :: stamp
    integer(int64) :: count
    call start_timer(stamp)
    do
      call system_clock(count)
      if (timer_elapsed(stamp,count)>=seconds) exit
    enddo
  end subroutine
end program
'''


class RuntimeTimingTest(unittest.TestCase):
    def test_profiler(self):
        compiler = os.environ.get('FC') or shutil.which('gfortran')
        if not compiler:
            self.skipTest('Set FC to a Fortran compiler or install gfortran')
        with tempfile.TemporaryDirectory(prefix="ecosim timer '") as directory:
            work = Path(directory)
            source = work / 'test.F90'
            source.write_text(SOURCE)
            exe = work / 'test'
            subprocess.run([compiler, '-cpp', '-fcheck=all', '-Wall', '-Wextra',
                            str(ROOT / 'f90src/Utils/data_kind_mod.F90'),
                            str(ROOT / 'f90src/Utils/timings.F90'), str(source),
                            '-o', str(exe)], cwd=work, check=True)
            subprocess.run([str(exe), str(work)], cwd=work, check=True)
            def read_summary(path):
                with path.open() as stream:
                    return {row['timer']: row for row in csv.DictReader(
                        line for line in stream if not line.startswith('#'))}
            rows = read_summary(work / 'timing/summary.csv')
            self.assertEqual(len(rows), 224)
            self.assertEqual(int(rows['inner']['calls']), 3)
            self.assertGreater(float(rows['outer']['total_s']), float(rows['inner']['total_s']) + .025)
            self.assertIn('a very long label, with a "quoted" name', rows)
            for row in rows.values():
                total, mean, minimum, maximum = map(float, (row['total_s'], row['mean_s'], row['min_s'], row['max_s']))
                self.assertAlmostEqual(total / int(row['calls']), mean)
                self.assertLessEqual(minimum, maximum)
                self.assertGreaterEqual(minimum, 0)
                self.assertGreaterEqual(float(row['run_percent']), 0)
            with (work / 'timing/time-20261008-120000.txt').open() as stream:
                detail = list(csv.DictReader(stream))
            self.assertEqual([row['step'] for row in detail if row['timer'] == 'setup'], ['-1'])
            inner = [row for row in detail if row['timer'] == 'inner']
            self.assertEqual([(row['step'], row['calls']) for row in inner], [('0', '1'), ('1', '2')])
            self.assertEqual(len([row for row in detail if row['timer'] == 'outer']), 1)
            self.assertEqual(set(read_summary(work / 'second/timing/summary.csv')), {'new_run'})
            self.assertEqual(list((work / 'second/timing').glob('time*.txt')), [])


if __name__ == '__main__':
    unittest.main()
