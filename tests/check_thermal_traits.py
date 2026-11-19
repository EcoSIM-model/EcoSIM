"""Check production thermal-trait loading and independent temperature responses.

Compile the production reader and three growth functions with a small NetCDF
fixture. Test absent/default traits, per-PFT overrides, invalid records, and
legacy response equivalence. All generated files stay in a temporary directory.
Example: python3 tests/check_thermal_traits.py --fc agf13 \
  --nf-config build/CONFIG/bin/nf-config
"""
import argparse
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

import numpy as np
from netCDF4 import Dataset, default_fillvals

ROOT = Path(__file__).resolve().parents[1]
TRAITS = dict(RG_ELOW=197500., RG_EHIGH=222500., CG_ELOW=197500.,
              CG_EHIGH=222500., LG_ELOW=197500., LG_EHIGH=218500.)


def procedure(source, name):
    match = re.search(r'^  (?:pure )?(?:subroutine|function) ' + name +
                      r'\b.*?^  end (?:subroutine|function) ' + name + r'\b[^\n]*',
                      source, re.M | re.S | re.I)
    assert match, name
    return match[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fc', default='gfortran')
    parser.add_argument('--nf-config', type=Path, required=True)
    args = parser.parse_args()
    info = (ROOT/'f90src/IOutils/PlantInfoMod.F90').read_text()
    math = (ROOT/'f90src/Plant_bgc/PlantMathFuncMod.F90').read_text()
    source = '''module fixture
use iso_fortran_env, only: r8=>real64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use netcdf
implicit none
real(r8), parameter :: RGASC=8.3143_r8
type file_desc
 integer :: fh
end type
type(file_desc) :: pft_nfid
contains
subroutine PrintInfo(message)
 character(len=*), intent(in) :: message
end subroutine
subroutine endrun(message)
 character(len=*), intent(in) :: message
 print *,message
 error stop 1
end subroutine
subroutine check_ret(status,message)
 integer, intent(in) :: status
 character(len=*), intent(in) :: message
 if(status/=nf90_noerr)call endrun(message//nf90_strerror(status))
end subroutine
'''
    source += procedure(info, 'ReadThermalTrait') + '\n'
    for part in ('root', 'canopy', 'leave'):
        source += procedure(math, f'calc_{part}_grow_tempf') + '\n'
    source += '''end module
program check
use fixture
implicit none
character(len=1024) :: filename
real(r8) :: values(2),single(1),energies(6,2),t,rt,st,legacy(3),actual(3)
integer :: i,j
call get_command_argument(1,filename)
call check_ret(nf90_open(trim(filename),nf90_nowrite,pft_nfid%fh),'open')
'''
    for i, (name, default) in enumerate(TRAITS.items(), 1):
        source += f"call ReadThermalTrait('{name}',values,1,{default}_r8)\n"
        source += f"call ReadThermalTrait('{name}',single,2,{default}_r8)\n"
        source += "if(single(1)/=values(2))error stop 'per-PFT/catalog mismatch'\n"
        source += f"energies({i},:)=values\nwrite(*,'(A,2ES25.16)')'{name}',values\n"
    source += '''call check_ret(nf90_close(pft_nfid%fh),'close')
! PFT 1 always has defaults; PFT 2 may have one independent override.
do j=1,1401
 t=253.15_r8+real(j-1,r8)*.05_r8
 rt=RGASC*t;st=710._r8*t
 legacy(1)=exp(25.229_r8-62500._r8/rt)/(1+exp((197500._r8-st)/rt)+exp((st-222500._r8)/rt))
 legacy(2)=legacy(1)
 legacy(3)=exp(24.269_r8-60000._r8/rt)/(1+exp((197500_r8-st)/rt)+exp((st-218500._r8)/rt))
 actual=[calc_root_grow_tempf(t,energies(1,1),energies(2,1)), &
         calc_canopy_grow_tempf(t,energies(3,1),energies(4,1)), &
         calc_leave_grow_tempf(t,energies(5,1),energies(6,1))]
 if(any(actual/=legacy))error stop 'legacy default response changed'
enddo
t=298.15_r8
write(*,'(A,3ES25.16)')'FACTORS', &
 calc_root_grow_tempf(t,energies(1,2),energies(2,2)), &
 calc_canopy_grow_tempf(t,energies(3,2),energies(4,2)), &
 calc_leave_grow_tempf(t,energies(5,2),energies(6,2))
end program
'''
    with tempfile.TemporaryDirectory(prefix='ecosim-thermal-') as temp:
        work = Path(temp)
        (work/'check.f90').write_text(source)
        config = str(args.nf_config.resolve())
        flags = lambda option: shlex.split(subprocess.check_output([config, option], text=True))
        subprocess.run([*shlex.split(args.fc), '-ffree-line-length-none', '-fbounds-check',
                        '-ffpe-trap=invalid,zero,overflow', '-O2', *flags('--fflags'),
                        'check.f90', '-o', 'check', *flags('--flibs')], cwd=work, check=True)

        def run(label, overrides=None, bad=None):
            path = work/(label+'.nc')
            with Dataset(path, 'w') as ds:
                ds.createDimension('pft', 2)
                ds.createVariable('WDLF', 'f8', ('pft',))[:] = 4.
                for key, value in (overrides or {}).items():
                    ds.createVariable(key, 'f8', ('pft',))[:] = [TRAITS[key], value]
                if bad:
                    kind, value = bad
                    dims = ('pft',)
                    if kind == 'axis':
                        ds.createDimension('wrong', 2)
                        dims = ('wrong',)
                    elif kind == 'rank':
                        ds.createDimension('extra', 1)
                        dims = ('pft', 'extra')
                    kw = {'fill_value': -777.} if kind == 'fill' else {}
                    var = ds.createVariable('LG_ELOW', 'f8', dims, **kw)
                    if kind == 'missing':
                        var.missing_value = value
                    if kind != 'unwritten':
                        var[:] = value
            proc = subprocess.run([str(work/'check'), str(path)], text=True, capture_output=True)
            if bad:
                assert proc.returncode != 0, label
                assert 'LG_ELOW' in proc.stdout, (label, proc.stdout, proc.stderr)
                return
            assert proc.returncode == 0, (label, proc.stdout, proc.stderr)
            rows = {line.split()[0]: np.array(line.split()[1:], float)
                    for line in proc.stdout.splitlines() if line.strip()}
            for key, default in TRAITS.items():
                np.testing.assert_array_equal(rows[key], [default, (overrides or {}).get(key, default)])
            return rows['FACTORS']

        baseline = run('absent')
        np.testing.assert_array_equal(run('explicit_defaults', TRAITS), baseline)
        for i, (key, value) in enumerate(TRAITS.items()):
            changed = run(key, {key: value+1000.})
            np.testing.assert_array_equal(np.delete(changed, i//2), np.delete(baseline, i//2))
            assert changed[i//2] != baseline[i//2], key
        invalid = [('nan', np.nan), ('inf', np.inf), ('zero', 0.), ('negative', -1.),
                   ('fill', -777.), ('missing', 999.), ('unwritten', 0.),
                   ('default_float_fill', default_fillvals['f4']),
                   ('axis', 197500.), ('rank', 197500.)]
        for label, value in invalid:
            run(label, bad=(label, value))
    print('PASS: 1401-temperature legacy equivalence, catalog/scalar loading, six independent overrides, 10 invalid inputs')


if __name__ == '__main__':
    main()
