"""Check production FL_POPT/FL_PSENS loading on catalog and single-PFT paths.

The new optional traits use hours and inverse hours respectively. Missing inputs
explicitly select resolved XDL and zero photoperiod sensitivity. No cultivar
sensitivity is inferred from XPPD. Run with --fc agf13 --nf-config BUILD/bin/nf-config.
"""
import argparse
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np
from netCDF4 import Dataset
from check_phenology_phases import procedure

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fc', default='gfortran')
    parser.add_argument('--nf-config', type=Path, required=True)
    args = parser.parse_args()
    reader = procedure((ROOT/'f90src/IOutils/PlantInfoMod.F90').read_text(), 'ReadFloralTrait')
    source = '''module fixture
use iso_fortran_env,only:r8=>real64
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use netcdf
implicit none
type file_desc
 integer :: fh
end type
type(file_desc) :: pft_nfid
contains
subroutine PrintInfo(message)
 character(len=*),intent(in)::message
end subroutine
subroutine endrun(message)
 character(len=*),intent(in)::message
 print *,message
 stop 1
end subroutine
subroutine check_ret(status,message)
 integer,intent(in)::status
 character(len=*),intent(in)::message
 if(status/=nf90_noerr)call endrun(message//nf90_strerror(status))
end subroutine
'''+reader+'''
end module
program check
use fixture
implicit none
character(len=1024)::filename
real(r8)::opt(2),sens(2),one(1)
call get_command_argument(1,filename)
call check_ret(nf90_open(trim(filename),nf90_nowrite,pft_nfid%fh),'open')
call ReadFloralTrait('FL_POPT',opt,1,-1._r8,-1._r8,24._r8)
call ReadFloralTrait('FL_PSENS',sens,1,0._r8,0._r8,1.e6_r8)
call ReadFloralTrait('FL_POPT',one,2,-1._r8,-1._r8,24._r8)
if(one(1)/=opt(2))error stop 'single-PFT photoperiod differs from catalog'
call ReadFloralTrait('FL_PSENS',one,2,0._r8,0._r8,1.e6_r8)
if(one(1)/=sens(2))error stop 'single-PFT sensitivity differs from catalog'
call check_ret(nf90_close(pft_nfid%fh),'close')
write(*,*) opt,sens
end program
'''
    config = str(args.nf_config.resolve())
    flags = lambda option: shlex.split(subprocess.check_output([config, option], text=True))
    with tempfile.TemporaryDirectory(prefix='ecosim-floral-traits-') as temp:
        work = Path(temp)
        (work/'check.f90').write_text(source)
        subprocess.run([*shlex.split(args.fc), '-ffree-line-length-none', '-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow', *flags('--fflags'),
                        'check.f90', '-o', 'check', *flags('--flibs')], cwd=work, check=True)

        def run(label, opt=None, sens=None, bad=None):
            path = work/(label+'.nc')
            with Dataset(path, 'w') as ds:
                ds.createDimension('pft', 2)
                ds.createVariable('WDLF', 'f8', ('pft',))[:] = 4.
                for name, values in [('FL_POPT', opt), ('FL_PSENS', sens)]:
                    if values is None:
                        continue
                    dims = ('pft',)
                    if bad == 'axis':
                        ds.createDimension('wrong', 2)
                        dims = ('wrong',)
                    if bad == 'rank':
                        ds.createDimension('extra', 1)
                        dims = ('pft', 'extra')
                    var = ds.createVariable(name, 'f8', dims)
                    if bad != 'unwritten':
                        var[:] = np.array(values).reshape(var.shape)
                    if bad == 'missing':
                        var.missing_value = values[1]
            proc = subprocess.run([str(work/'check'), str(path)], capture_output=True, text=True)
            if bad or label.startswith('invalid'):
                assert proc.returncode != 0, (label, proc.stdout)
            else:
                assert proc.returncode == 0, (label, proc.stdout, proc.stderr)
                values = np.fromstring(proc.stdout, sep=' ')
                np.testing.assert_array_equal(values, [*(opt if opt is not None else [-1., -1.]),
                                                      *(sens if sens is not None else [0., 0.])])

        run('defaults')
        run('overrides', [11.5, 13.], [.2, .4])
        run('sensitivity_only', sens=[0., .3])
        for label, values in [('negative', [0., -.1]), ('nan', [0., np.nan]),
                              ('infinity', [0., np.inf]), ('huge', [0., 1.e100])]:
            run('invalid_'+label, sens=values)
        run('invalid_daylength', opt=[12., 25.])
        for bad in ['axis', 'rank', 'unwritten', 'missing']:
            run(bad, sens=[0., .4], bad=bad)
    print('PASS: floral defaults, overrides, catalog/single-PFT equivalence and invalid NetCDF records')


if __name__ == '__main__':
    main()
