"""Exercise the production cohort integrator with a standalone Fortran driver."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class MediumRootCohorts(unittest.TestCase):
    def test_investment_geometry_and_conservation(self):
        compiler = shutil.which('agf13') or shutil.which('gfortran')
        if not compiler:
            self.skipTest('Fortran compiler unavailable')
        stub = '''module data_kind_mod
use iso_fortran_env, only: DAT_KIND_R8=>real64
end module
module DebugToolMod
contains
subroutine PrintInfo(message)
character(len=*),intent(in)::message
end subroutine
end module
'''
        driver = '''program test
use MediumRootCohortMod
use data_kind_mod,only:r8=>DAT_KIND_R8
implicit none
real(r8)::n,l,e(3),yn,yl,ye(3),g(3),nn(2),ll(2),rr(2),before(3)
integer::i
! A vast count target with tiny positive growth must recruit only paid length.
n=10;l=100;e=[100._r8,2._r8,0.1_r8];yn=0;yl=0;ye=0
g=[0.001_r8,0.00002_r8,0.000001_r8];before=e
call GrowMediumCohorts(g,10000._r8,1._r8,0._r8,0.01_r8,0.1_r8,0.001_r8,n,l,e,yn,yl,ye)
if(n-10._r8>0.10000001_r8.or.n<=10._r8)stop 1
if(abs((l-100._r8)-g(1))>1.e-12_r8)stop 2
if(any(abs(e-before-g)>1.e-12_r8))stop 3
if(any(ye<0._r8).or.any(ye>e))stop 4
! Recruitment cannot spend old length or shorten the established cohort.
call MediumCohortGeometry(n,l,e,yn,yl,ye,1._r8,nn,ll,rr)
if(ll(2)<100._r8.or.abs(ll(2)/nn(2)-10._r8)>1.e-4_r8)stop 5
! Zero growth has no recruitment even with an arbitrarily large target.
before=e;g=0
call GrowMediumCohorts(g,1.e8_r8,1._r8,0._r8,0.01_r8,0.1_r8,0.001_r8,n,l,e,yn,yl,ye)
if(any(e/=before))stop 6
! Repeated recruitment, extension and maturation must conserve aggregate CNP.
n=0;l=0;e=0;yn=0;yl=0;ye=0
g=[0.1_r8,0.002_r8,0.0001_r8]
do i=1,1000
before=e
call GrowMediumCohorts(g,100._r8,1._r8,0.2_r8,0.01_r8,0.1_r8,0.001_r8,n,l,e,yn,yl,ye)
if(any(abs(e-before-g)>1.e-12_r8))stop 7
if(yn<0._r8.or.yn>n.or.yl<0._r8.or.yl>l.or.any(ye<0._r8).or.any(ye>e))stop 8
if(n>100.000001_r8.or.l<yn*0.01_r8)stop 9
enddo
if(n-yn<=0._r8.or.yn<=0._r8)stop 10
call MediumCohortGeometry(n,l,e,yn,yl,ye,1._r8,nn,ll,rr)
if(abs(sum(nn)-n)>1.e-12_r8.or.abs(sum(ll)-l)>1.e-12_r8)stop 11
if(abs(sum(acos(-1._r8)*rr**2*ll)-e(1))>1.e-10_r8)stop 12
print *, 'PASS: recruitment cost, extension competition, CNP closure, cohort geometry, maturation'
end program
'''
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / 'stubs.f90').write_text(stub)
            (tmp / 'driver.f90').write_text(driver)
            command = [compiler, '-ffree-line-length-none', '-fcheck=all',
                       '-ffpe-trap=invalid,zero,overflow', str(tmp / 'stubs.f90'),
                       str(ROOT / 'f90src/Plant_bgc/MediumRootCohortMod.F90'),
                       str(tmp / 'driver.f90'), '-o', str(tmp / 'test')]
            subprocess.run(command, cwd=tmp, check=True, capture_output=True, text=True)
            result = subprocess.run([str(tmp / 'test')], cwd=tmp, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('PASS:', result.stdout)


if __name__ == '__main__':
    unittest.main()
