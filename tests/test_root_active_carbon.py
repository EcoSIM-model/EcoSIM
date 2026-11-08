"""Exercise compiled root-pool aggregation and the production FWTC/FWTS block.

Build EcoSIM first. ECOSIM_TEST_BUILD may select a configured build directory.
The compiler, ABI flags, and libraries are taken from that build, not PATH.
"""
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BUILD = Path(os.environ.get('ECOSIM_TEST_BUILD', ROOT / 'build/Darwin-arm64-double-Release'))


class RootActiveCarbonTest(unittest.TestCase):
    def test_current_pools_and_paired_weights(self):
        linkfile = BUILD / 'drivers/ecosim/CMakeFiles/ecosim.f90.x.dir/link.txt'
        self.assertTrue(linkfile.is_file(), 'Build EcoSIM before running this test')
        source = (ROOT / 'f90src/Plant_bgc/PlantNonstElmDynMod.F90').read_text()
        start = source.index('  !canopy sink weighting factor', source.index('  subroutine ShootRootElmTransfer'))
        end = source.index('  IF(RootSinkC(ipltroot)', start)
        weights = source[start:end]
        harness = '''program test_root_active
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use ElmIDMod, only: ielmc
  use PlantBiomassAPIData, only: plt_biom
  use PlantMorphologyAPIData, only: plt_morph
  use PlantSiteAPIData, only: plt_site
  use PlantBalMod, only: UpdateRootActiveC
  implicit none
  integer :: NZ, i
  real(r8) :: FWTC, FWTS
  real(r8) :: CanopyLeafSheathC_pft(1), ZERO4Groth_pft(1), RootElms_pft(3,1)
  real(r8) :: shoots(4), expected_c(4), expected_s(4)
  allocate(plt_morph%Myco_pft(1),plt_morph%NumStructuralRootAxes_pft(1))
  allocate(plt_biom%RootActiveC_pft(1))
  allocate(plt_biom%RootMycoNonstElms_rpvr(3,2,2,1))
  allocate(plt_biom%RootMyco2ndStrutElms_rpvr(3,2,2,2,1))
  allocate(plt_biom%Root1stActStructElms_rpvr(3,2,2,1))
  allocate(plt_biom%Root1stLigStructElms_rpvr(3,2,2,1))
  allocate(plt_biom%RootMediumStructElms_rpvr(3,2,2,1))
  plt_site%MaxNumRootLays=2
  plt_morph%Myco_pft=2
  plt_morph%NumStructuralRootAxes_pft=2
  associate(active => plt_biom%RootActiveC_pft, &
    reserve => plt_biom%RootMycoNonstElms_rpvr, &
    fine => plt_biom%RootMyco2ndStrutElms_rpvr, &
    primary => plt_biom%Root1stActStructElms_rpvr, &
    inactive => plt_biom%Root1stLigStructElms_rpvr, &
    medium => plt_biom%RootMediumStructElms_rpvr)
    reserve=0._r8;fine=0._r8;primary=0._r8;inactive=0._r8;medium=0._r8
    reserve(ielmc,:,:,:)=4._r8
    fine(ielmc,:,:,:,:)=1._r8
    primary(ielmc,:,:,:)=2._r8
    medium(ielmc,:,:,:)=3._r8
    active=-999._r8
    call UpdateRootActiveC(1)
    if(abs(active(1)-44._r8)>1.e-12_r8)stop 1
    inactive(ielmc,:,:,:)=1.e6_r8
    call UpdateRootActiveC(1)
    if(abs(active(1)-44._r8)>1.e-12_r8)stop 2
    primary(ielmc,1,1,1)=7._r8
    call UpdateRootActiveC(1)
    if(abs(active(1)-49._r8)>1.e-12_r8)stop 3
    primary(ielmc,1,1,1)=2._r8
  end associate
  call UpdateRootActiveC(1)
  ZERO4Groth_pft=1.e-12_r8
  RootElms_pft=1.e6_r8  ! Total/inactive C must not drive either weight.
  shoots=(/100._r8,10._r8,10._r8,0._r8/)
  expected_c=(/0.29348_r8,1._r8,0._r8,1._r8/)
  expected_s=(/1._r8,10._r8/29.348_r8,1._r8,1._r8/)
  NZ=1
  associate(RootActiveC_pft => plt_biom%RootActiveC_pft)
    do i=1,4
      CanopyLeafSheathC_pft(1)=shoots(i)
      if(i>=3)RootActiveC_pft=0._r8
''' + weights + '''
      if(abs(FWTC-expected_c(i))>1.e-12_r8)stop 4
      if(abs(FWTS-expected_s(i))>1.e-12_r8)stop 5
    enddo
  end associate
  print *, 'PASS: current active pools, inactive-C exclusion, both weights, and empty pools'
end program test_root_active
'''
        command = shlex.split(linkfile.read_text())
        if '--' in command:
            command = command[command.index('--') + 1:]
        cleaned = []
        i = 0
        while i < len(command):
            if command[i] == '-o':
                i += 2
                continue
            if not command[i].endswith('.o'):
                cleaned.append(command[i])
            i += 1
        includes = []
        for mod in BUILD.rglob('*.mod'):
            inc = '-I' + str(mod.parent.resolve())
            if inc not in includes:
                includes.append(inc)
        with tempfile.TemporaryDirectory(prefix='ecosim_rootactive_') as work:
            program = Path(work) / 'test.F90'
            exe = Path(work) / 'test.x'
            program.write_text(harness)
            proc = subprocess.run(cleaned[:1] + includes + [str(program), '-o', str(exe)] + cleaned[1:],
                                  cwd=BUILD / 'drivers/ecosim', capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            proc = subprocess.run([str(exe)], cwd=work, capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            self.assertIn('PASS:', proc.stdout)


if __name__ == '__main__':
    unittest.main()
