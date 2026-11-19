"""Compile production crop induction routines and test their developmental behavior.

Run with --fc agf13. Checks fixed/increasing/decreasing photoperiods, both response
signs, exposure history, node competence, stress increments, branch isolation,
physical-clock invariance, frozen completed induction, and the legacy perennial
bypass. Integration/restart tests separately exercise state transfer and NetCDF.
"""
import argparse
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

from check_phenology_phases import procedure

ROOT = Path(__file__).resolve().parents[1]


def fixture():
    source = (ROOT/'f90src/Plant_bgc/PlantPhenolMod.F90').read_text()
    routines = '\n'.join(procedure(source, name) for name in
                         ['UpdateBranchNodeNumber', 'CheckBranchNodeState', 'InitiateBranchFloral'])
    apis = {'plt_pheno': 'PlantPhenologyAPIData', 'plt_site': 'PlantSiteAPIData',
            'plt_morph': 'PlantMorphologyAPIData', 'plt_distb': 'PlantDisturbanceAPIData',
            'plt_rbgc': 'PlantRootBGCAPIData', 'plt_ew': 'PlantEnergyWaterAPIData'}
    declarations = []
    for instance in sorted(set(re.findall(r'(plt_\w+)%', routines))):
        api = (ROOT/'f90src/APIData'/f'{apis[instance]}.F90').read_text()
        fields = []
        for field in sorted(set(re.findall(instance+r'%(\w+)', routines))):
            decl = re.search(r'^.*::\s*'+field+r'\b[^\n]*', api, re.M | re.I)[0]
            kind, initial = ('integer', '0') if re.search(r'^\s*integer', decl, re.I) else ('real(r8)', '0._r8')
            alloc = re.search(r'allocate\(this%'+field+r'\((.*?)\)\)', api, re.I)
            dims = '('+alloc[1]+')' if alloc else ''
            fields.append(f'  {kind} :: {field}{dims}={initial}')
        declarations += [f'type {instance}_type', *fields, 'end type',
                         f'type({instance}_type) :: {instance}']
    mini = (ROOT/'f90src/Minimath/MiniMathMod.F90').read_text()
    math = (ROOT/'f90src/Plant_bgc/PlantMathFuncMod.F90').read_text()
    helpers = '\n'.join(procedure(mini, name) for name in ['AZMAX1_s', 'AZMAX1_d'])
    helpers += '\n'+procedure(math, 'is_cold_deciduos')
    return '''module fixture
use iso_fortran_env, only: r8 => real64
use ElmIDMod
implicit none
integer, parameter :: JP1=1, MaxNumBranches=2, NumGrowthStages=10, NumCanopyLayers1=2
real(r8), parameter :: tiny_val=1.e-12_r8
real(r8), parameter :: GrothStageNorm4VegetaPheno=2._r8, GrothStageNorm4ReprodPheno=.667_r8
interface AZMAX1
 module procedure AZMAX1_s,AZMAX1_d
end interface
'''+ '\n'.join(declarations)+'''
contains
subroutine PrintInfo(message)
 character(len=*),intent(in) :: message
end subroutine
'''+helpers+'\n'+routines+'\n'+HELPERS+'\nend module\n'


HELPERS = '''subroutine setup(response,day,prev,sensitivity)
 integer,intent(in) :: response
 real(r8),intent(in) :: day,prev,sensitivity
 plt_pheno%iPlantCalendar_brch=0
 plt_pheno%iPlantPhenolPattern_pft=iplt_annual
 plt_pheno%iPlantPhenolType_pft=iphenotyp_evgreen
 plt_pheno%iPlantPhotoperiodType_pft=response
 plt_pheno%iPlantDevelopPattern_pft=ideterminate
 plt_pheno%MatureGroup_pft=2._r8
 plt_pheno%MatureGroup_brch=2._r8
 plt_pheno%FloralDevNodes_brch=0._r8
 plt_pheno%RefNodeInitRate_pft=.25_r8
 plt_pheno%RateRefLeafAppearance_pft=.2_r8
 plt_pheno%FloralPhotoOpt_pft=12._r8
 plt_pheno%FloralPhotoSens_pft=sensitivity
 plt_pheno%Hours4Leafout_brch=0._r8
 plt_pheno%HourReq4LeafOut_brch=999999._r8
 plt_pheno%CriticPhotoPeriod_pft=12._r8
 plt_pheno%PhotoPeriodSens_pft=.25_r8
 plt_pheno%LeafNumberAtFloralInit_brch=0._r8
 plt_morph%ShootNodeNum_brch=4._r8
 plt_morph%ShootNodeNumAtInitFloral_brch=4._r8
 plt_morph%NumOfLeaves_brch=0._r8
 plt_morph%CanopyHeightLive_pft=0._r8
 plt_site%DayLenthCurrent=day
 plt_site%DayLenthPrev=prev
 plt_site%iYearCurrent=2001
 plt_distb%iDayPlanting_pft=1
 plt_distb%iYearPlanting_pft=2001
 plt_ew%SnowDepth=0._r8
end subroutine
integer function flowering_hour(response,day,prev,sensitivity)result(hour)
 integer,intent(in) :: response
 real(r8),intent(in) :: day,prev,sensitivity
 call setup(response,day,prev,sensitivity)
 do hour=1,100
  call UpdateBranchNodeNumber(100,hour,1,1,1._r8,1._r8,1._r8)
  call InitiateBranchFloral(100,hour,1,1,day<prev)
  if(plt_pheno%iPlantCalendar_brch(ipltcal_InitFloral,1,1)/=0)return
 enddo
 error stop 'no floral initiation under constant finite photoperiod'
end function
'''

DRIVER = '''program check_crop_floral
use fixture
implicit none
integer :: h,i
real(r8) :: held,ln2
ln2=log(2._r8)
! Nine quarter-node hours exceed a two-node requirement at any daylength direction.
if(flowering_hour(iphotop_short,11._r8,10._r8,ln2)/=9)error stop 'increasing inductive days blocked'
if(flowering_hour(iphotop_short,11._r8,11._r8,ln2)/=9)error stop 'constant inductive days blocked'
if(flowering_hour(iphotop_short,11._r8,12._r8,ln2)/=9)error stop 'shortening days bypass competence'
if(flowering_hour(iphotop_short,13._r8,13._r8,ln2)/=17)error stop 'quantitative short-day delay lost'
if(flowering_hour(iphotop_long,11._r8,11._r8,ln2)/=17)error stop 'long-day response sign wrong'
if(flowering_hour(iphotop_long,13._r8,13._r8,ln2)/=9)error stop 'long-day inductive rate wrong'
if(flowering_hour(iphotop_neutral,20._r8,21._r8,100._r8)/=9)error stop 'neutral timing changed'
! Completed induction remains frozen; physical clocks continue unchanged.
held=plt_pheno%FloralDevNodes_brch(1,1)
call UpdateBranchNodeNumber(100,10,1,1,1._r8,1._r8,1._r8)
if(plt_pheno%FloralDevNodes_brch(1,1)/=held)error stop 'completed induction advanced'
if(plt_morph%ShootNodeNum_brch(1,1)/=6.5_r8)error stop 'photoperiod changed physical nodes'
if(abs(plt_morph%NumOfLeaves_brch(1,1)-2._r8)>1.e-12_r8)error stop 'physical leaf clock changed'
if(plt_pheno%FloralDevNodes_brch(2,1)/=0._r8)error stop 'branch histories coupled'
! Exposure history: four half-rate hours followed by seven full-rate hours.
call setup(iphotop_short,13._r8,13._r8,ln2)
do i=1,11
 if(i==5)plt_site%DayLenthCurrent=11._r8
 call UpdateBranchNodeNumber(100,i,1,1,1._r8,1._r8,1._r8)
 call InitiateBranchFloral(100,i,1,1,.true.)
 if(i<11.and.plt_pheno%iPlantCalendar_brch(ipltcal_InitFloral,1,1)/=0)error stop 'history accumulated too fast'
enddo
if(plt_pheno%iPlantCalendar_brch(ipltcal_InitFloral,1,1)/=100)error stop 'induction history not retained'
! Even a large induction state cannot bypass physical competence.
call setup(iphotop_short,11._r8,12._r8,ln2)
plt_pheno%FloralDevNodes_brch=100._r8
call InitiateBranchFloral(100,1,1,1,.true.)
if(plt_pheno%iPlantCalendar_brch(ipltcal_InitFloral,1,1)/=0)error stop 'competence bypassed'
! The accumulated increment includes the same temperature/water/oxygen factors.
call setup(iphotop_short,13._r8,13._r8,ln2)
call UpdateBranchNodeNumber(100,1,1,1,.5_r8,.8_r8,.5_r8)
if(abs(plt_pheno%FloralDevNodes_brch(1,1)-.025_r8)>1.e-12_r8)error stop 'stress increment wrong'
! Extreme sensitivity remains finite under the model's trapped-underflow build.
call setup(iphotop_short,24._r8,24._r8,1.e6_r8)
call UpdateBranchNodeNumber(100,1,1,1,.001_r8,1._r8,1._r8)
if(plt_pheno%FloralDevNodes_brch(1,1)<=0._r8)error stop 'numerical response floor lost'
! Noncrop perennial shortening-day behavior is retained.
call setup(iphotop_long,11._r8,12._r8,ln2)
plt_pheno%iPlantPhenolPattern_pft=iplt_perennial
plt_pheno%iPlantPhenolType_pft=iphenotyp_coldecid
call InitiateBranchFloral(100,1,1,1,.true.)
if(plt_pheno%iPlantCalendar_brch(ipltcal_InitFloral,1,1)/=100)error stop 'legacy perennial route changed'
print *, 'PASS: absolute photoperiod, induction memory, competence, clocks, stress, branch isolation and legacy perennial'
end program
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fc', default='gfortran')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='ecosim-crop-floral-') as temp:
        work = Path(temp)
        (work/'check.F90').write_text(fixture()+DRIVER)
        subprocess.run([*shlex.split(args.fc), '-ffree-line-length-none', '-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow,underflow',
                        str(ROOT/'f90src/Modelconfig/ElmIDMod.F90'), str(work/'check.F90'),
                        '-o', str(work/'check')], cwd=work, check=True)
        subprocess.run([str(work/'check')], cwd=work, check=True)


if __name__ == '__main__':
    main()
