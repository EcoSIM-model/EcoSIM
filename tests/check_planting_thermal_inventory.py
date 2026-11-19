"""Compile production planting initialization against retained dead/wet fixtures.

Verify that planting preserves retained donor heat and books net water/energy
imports, including cold-start and mixed old/new material. No simulation input
files or model outputs are modified.
"""
import argparse
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def procedure(path, name):
    text = path.read_text()
    match = re.search(r'^  (?:pure )?(?:subroutine|function) '+name+
                      r'\b.*?^  end (?:subroutine|function) '+name+r'\b[^\n]*',
                      text, re.M | re.S | re.I)
    assert match, name
    return match[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--fc', required=True)
    args = ap.parse_args()
    init = procedure(ROOT/'f90src/Plant_bgc/InitPlantMod.F90', 'InitPlantHeatWater')
    fields = {}
    for instance, field in re.findall(r'(plt_\w+)%(\w+)', init):
        fields.setdefault(instance, set()).add(field)
    declarations = []
    for instance, names in fields.items():
        declarations.append(f'type {instance}_type')
        for name in sorted(names):
            dims = '' if name == 'ATCA' or name.endswith('_col') else '(1)'
            if name in ('ShootElms_pft', 'StandDeadStrutElms_pft'):
                dims = '(3,1)'
            declarations.append(f' real(r8) :: {name}{dims}=0._r8')
        declarations += ['end type', f'type({instance}_type) :: {instance}']
    constants = (ROOT/'f90src/Utils/EcoSimConst.F90').read_text()
    selected = []
    for name in ('cpw', 'cpo', 'cps', 'ppmc', 'OMCMassFrac', 'DensitySolidOM'):
        selected.append(re.search(r'^.*:: '+name+r'\s*=.*$', constants, re.M)[0])
    helpers = '\n'.join(procedure(ROOT/'f90src/Plant_bgc/PlantMathFuncMod.F90', name)
                        for name in ('CanopyStemThermalCarbon', 'get_FDM'))
    helpers += '\n'+procedure(ROOT/'f90src/Minimath/MiniFuncMod.F90', 'gOC_to_m3_OM')
    source = '''module fixture
use iso_fortran_env, only: r8=>real64,iulog=>error_unit
use ieee_arithmetic, only: ieee_is_finite
implicit none
integer, parameter :: ielmc=1
integer :: j
real(r8) :: before_energy,before_water,retained,added,temperature
real(r8), parameter :: tol=1.e-12_r8
character(len=1),parameter :: mod_filename='x'
type units_type
contains
 procedure,nopass :: Celcius2Kelvin
 procedure,nopass :: Kelvin2Celcius
end type
type(units_type) :: units
'''+ '\n'.join(selected+declarations)+'''
contains
pure function Celcius2Kelvin(x) result(y)
 real(r8),intent(in)::x
 real(r8)::y
 y=x+273.15_r8
end function
pure function Kelvin2Celcius(x) result(y)
 real(r8),intent(in)::x
 real(r8)::y
 y=x-273.15_r8
end function
pure function AZMAX1(x) result(y)
 real(r8),intent(in)::x
 real(r8)::y
 y=max(0._r8,x)
end function
subroutine PrintInfo(x)
 character(len=*),intent(in)::x
end subroutine
subroutine endrun(x,line)
 character(len=*),intent(in)::x
 integer,intent(in)::line
 error stop 'unexpected endrun'
end subroutine
'''+helpers+'\n'+init+'''
end module
program test
use fixture
implicit none
! Three fixtures: cold seed; dead/wet residue; residue plus newly added tissue.
do j=1,3
 plt_site%ATCA=7._r8
 plt_morph%StemSpecVolume_pft=1.e-5_r8
 plt_biom%StandDeadStrutElms_pft=0._r8
 plt_morph%StandDeadSurfArea_pft=0._r8
 plt_ew%WatHeldOnCanopy_pft=0._r8
 plt_ew%SnowOnCanopy_pft=0._r8
 plt_ew%CanopyBiomWater_pft=0._r8
 plt_biom%CanopyLeafSheathC_pft=10._r8
 plt_biom%CanopySapwoodC_pft=0._r8
 plt_ew%TKC_pft=260._r8
 retained=0._r8
 if(j>1)then
   plt_biom%StandDeadStrutElms_pft(1,1)=1000._r8
   plt_morph%StandDeadSurfArea_pft=3._r8
   plt_ew%WatHeldOnCanopy_pft=.002_r8
   plt_ew%SnowOnCanopy_pft=.004_r8
   plt_ew%CanopyBiomWater_pft=.003_r8
   retained=cpo*(750.e-6_r8/(.55_r8*1.30_r8))+cpw*.002_r8+cps*.004_r8
 endif
 if(j==2)plt_biom%CanopyLeafSheathC_pft=0._r8
 before_water=plt_ew%CanopyBiomWater_pft(1)
 before_energy=(retained+cpw*before_water)*260._r8
 plt_ew%ENGYX_pft=before_energy
 plt_ew%CanopyHeatLoss2Dist_col=0._r8
 plt_ew%QCanopyWatLoss2Dist_col=0._r8
 added=cpo*(plt_biom%CanopyLeafSheathC_pft(1)*1.e-6_r8/(.55_r8*1.30_r8)) &
   +cpw*ppmc*plt_biom%CanopyLeafSheathC_pft(1)/get_FDM(real(-1.0e-03,r8))
 temperature=280.15_r8
 if(retained>0._r8)temperature=(retained*260._r8+added*280.15_r8)/(retained+added)
 call InitPlantHeatWater(1)
 if(abs(plt_ew%TKC_pft(1)-temperature)>tol)error stop 'retained temperature'
 if(abs(plt_ew%ENGYX_pft(1)-(retained*260._r8+added*280.15_r8))>tol)error stop 'mixed heat'
 if(abs(plt_ew%ENGYX_pft(1)-before_energy+plt_ew%CanopyHeatLoss2Dist_col)>tol)error stop 'heat ledger'
 if(abs(plt_ew%CanopyBiomWater_pft(1)-before_water+plt_ew%QCanopyWatLoss2Dist_col)>tol) &
   error stop 'water ledger'
enddo
print *, 'planting thermal inventories passed'
end program
'''
    with tempfile.TemporaryDirectory(prefix='planting-heat-') as tmp:
        path = Path(tmp)/'test.F90'
        path.write_text(source)
        exe = Path(tmp)/'test.x'
        subprocess.run([args.fc, '-cpp', '-ffree-line-length-none', '-fcheck=all',
                        str(path), '-o', str(exe)], cwd=tmp, check=True)
        subprocess.run([str(exe)], cwd=tmp, check=True)


if __name__ == '__main__':
    main()
