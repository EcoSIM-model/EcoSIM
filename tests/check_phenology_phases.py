"""Exercise preparation and advancement boundaries using production procedures.

Compile the actual phase routines and temperature/rate equations against small
cohort fixtures. Verify planting/eligibility before exchange, no early clock
advance, current-temperature progression, and inactive/dormant behavior.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def procedure(text, name):
    pattern = (r'^  (?:(?:pure|integer)\s+)?(?:subroutine|function) ' + name +
               r'\b.*?^  end (?:subroutine|function) ' + name + r'\b[^\n]*')
    match = re.search(pattern, text, re.M | re.S | re.I)
    assert match, name
    return match[0]


def fixture(phenology):
    routines = '\n'.join(procedure(phenology, name) for name in
                         ['PreparePlantPhenology', 'PhenologyUpdate',
                          'CalcPhenolEnvfactor', 'UpdateBranchNodeNumber'])
    apis = {'plt_pheno': 'PlantPhenologyAPIData', 'plt_site': 'PlantSiteAPIData',
            'plt_morph': 'PlantMorphologyAPIData', 'plt_ew': 'PlantEnergyWaterAPIData',
            'plt_rbgc': 'PlantRootBGCAPIData'}
    declarations = []
    for instance in sorted(set(re.findall(r'(plt_\w+)%', routines))):
        api = (ROOT/'f90src/APIData'/f'{apis[instance]}.F90').read_text()
        fields = []
        for field in sorted(set(re.findall(instance+r'%(\w+)', routines))):
            decl = re.search(r'^.*::\s*'+field+r'\b[^\n]*', api, re.M | re.I)[0]
            if re.search(r'^\s*character', decl, re.I):
                kind, initial = 'character(len=16)', "''"
            elif re.search(r'^\s*integer', decl, re.I):
                kind, initial = 'integer', '0'
            else:
                kind, initial = 'real(r8)', '0._r8'
            alloc = re.search(r'allocate\(this%'+field+r'\((.*?)\)\)', api, re.I)
            dims = '('+alloc[1]+')' if alloc else ''
            fields.append(f'  {kind} :: {field}{dims}={initial}')
        declarations += [f'type {instance}_type', *fields, 'end type',
                         f'type({instance}_type) :: {instance}']
    mini = (ROOT/'f90src/Minimath/MiniMathMod.F90').read_text()
    math = (ROOT/'f90src/Plant_bgc/PlantMathFuncMod.F90').read_text()
    helpers = '\n'.join(procedure(mini, name) for name in ['AZMAX1_s', 'AZMAX1_d'])
    helpers += '\n'+procedure(math, 'calc_leave_grow_tempf')
    return '''module fixture
use iso_fortran_env, only: r8 => real64
use ElmIDMod
implicit none
integer, parameter :: JP1=3, MaxNumBranches=2, NumGrowthStages=10
real(r8), parameter :: RGASC=8.3143_r8
real(r8), parameter :: tiny_val=1.e-12_r8
real(r8), parameter :: GrothStageNorm4VegetaPheno=2._r8, GrothStageNorm4ReprodPheno=.667_r8
type yearIJ_type
 integer :: I=130,J=1
end type
integer :: flag_calls=0,stage_calls=0,emergence_calls=0,branch_calls=0,advance_calls=0
logical :: new_cohort=.true.,allow_emergence=.true.
interface AZMAX1
 module procedure AZMAX1_s,AZMAX1_d
end interface
'''+'\n'.join(declarations)+'''
contains
subroutine PrintInfo(message)
 character(len=*),intent(in) :: message
end subroutine
subroutine set_plant_flags(I,J,NZ)
 integer,intent(in) :: I,J,NZ
 flag_calls=flag_calls+1
 if(NZ==1.and.new_cohort)then
  plt_site%PlantPopuLive_pft(NZ)=12._r8
  plt_pheno%IsPlantActive_pft(NZ)=itrue
  new_cohort=.false.
 endif
end subroutine
subroutine FindMainBranchNumber(NZ)
 integer,intent(in) :: NZ
 plt_morph%MainBranchNum_pft(NZ)=1
end subroutine
subroutine StagePlantPhenology(I,J,NZ)
 integer,intent(in) :: I,J,NZ
 stage_calls=stage_calls+1
end subroutine
subroutine TestPlantEmergence(I,J,NZ)
 integer,intent(in) :: I,J,NZ
 emergence_calls=emergence_calls+1
 if(allow_emergence)plt_pheno%iPlantCalendar_brch(ipltcal_Emerge,1,NZ)=I
end subroutine
subroutine root_shoot_branching(I,J,NZ)
 integer,intent(in) :: I,J,NZ
 branch_calls=branch_calls+1
end subroutine
subroutine Emerged_plant_Phenology(I,J,NZ)
 integer,intent(in) :: I,J,NZ
 real(r8) :: t,w,o
 advance_calls=advance_calls+1
 call CalcPhenolEnvfactor(I,J,NZ,t,w,o)
 call UpdateBranchNodeNumber(I,J,1,NZ,t,w,o)
end subroutine
'''+helpers+'\n'+routines+'\nend module\n'


DRIVER = '''program check_phenology_phases
use fixture
implicit none
type(yearIJ_type) :: hour
real(r8) :: expected,before
plt_site%NP=3
plt_site%DATAP=['soyb41          ','maiz41          ','NO              ']
plt_site%PlantPopuLive_pft=[100._r8,0._r8,0._r8]
plt_site%PlantPopu_col=777._r8
plt_pheno%IsPlantActive_pft=ifalse
plt_pheno%dReproNodeNumNormByMatG_brch=99._r8
plt_pheno%iPlantPhenolPattern_pft=iplt_annual
plt_pheno%RefNodeInitRate_pft=.013_r8
plt_pheno%RateRefLeafAppearance_pft=.010_r8
plt_pheno%MatureGroup_pft=8._r8
plt_pheno%TKGroth_pft=283.15_r8
plt_pheno%TempOffset_pft=1.3335_r8
plt_pheno%PlantO2Stress_pft=.25_r8
plt_ew%PSICanopy_pft=-1._r8
call PreparePlantPhenology(hour)
if(flag_calls/=2.or.stage_calls/=1.or.emergence_calls/=1)error stop 'cohort preparation skipped'
if(plt_site%PlantPopu_col/=12._r8.or.plt_site%NumActivePlants/=1)error stop 'wrong initialized population'
if(plt_pheno%iPlantCalendar_brch(ipltcal_Emerge,1,1)/=hour%I)error stop 'emergence unavailable to hydraulics'
if(any(plt_morph%NumOfLeaves_brch/=0._r8).or.any(plt_morph%ShootNodeNum_brch/=0._r8))error stop 'preparation advanced development'
if(any(plt_pheno%dReproNodeNumNormByMatG_brch/=0._r8))error stop 'stale reproductive increment'
if(branch_calls/=0.or.advance_calls/=0)error stop 'early branching or clock update'
! Represent the intervening hydraulic/chemical solve with changed conditions.
plt_pheno%TKGroth_pft=298.15_r8
plt_ew%PSICanopy_pft=-2._r8
expected=calc_leave_grow_tempf(299.4835_r8)*exp(-.05_r8)*.5_r8
call PhenologyUpdate(hour)
if(branch_calls/=1.or.advance_calls/=1)error stop 'wrong advancement count'
if(abs(plt_morph%NumOfLeaves_brch(1,1)-.010_r8*expected)>1.e-12_r8)error stop 'leaf clock used old environment'
if(abs(plt_morph%ShootNodeNum_brch(1,1)-.013_r8*expected)>1.e-12_r8)error stop 'node clock used old environment'
if(any(plt_morph%NumOfLeaves_brch(:,2:)/=0._r8))error stop 'inactive cohort progressed'
if(flag_calls/=2.or.stage_calls/=1.or.emergence_calls/=1)error stop 'progression repeated initialization'
if(plt_site%PlantPopu_col/=12._r8)error stop 'progression changed population total'
! Next-hour preparation neither advances again nor accumulates population twice.
before=plt_morph%NumOfLeaves_brch(1,1)
hour%J=2
call PreparePlantPhenology(hour)
if(plt_site%PlantPopu_col/=12._r8.or.plt_site%NumActivePlants/=1)error stop 'hourly totals accumulated'
if(plt_morph%NumOfLeaves_brch(1,1)/=before.or.advance_calls/=1)error stop 'preparation advanced twice'
! Unemerged managed crops wait; dormant-seed path retains its explicit exception.
allow_emergence=.false.
plt_pheno%iPlantCalendar_brch(ipltcal_Emerge,1,1)=0
plt_pheno%doInitPlant_pft=ifalse
call PreparePlantPhenology(hour)
call PhenologyUpdate(hour)
if(advance_calls/=1.or.plt_morph%NumOfLeaves_brch(1,1)/=before)error stop 'unemerged crop clock advanced'
plt_pheno%doInitPlant_pft(1)=itrue
call PhenologyUpdate(hour)
if(advance_calls/=2)error stop 'dormant seed exception lost'
print *, 'PASS: preparation, initialized totals, phase isolation, current environment and dormant/inactive cohorts'
end program
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fc', default=os.environ.get('FC', 'gfortran'))
    parser.add_argument('--phenology-source', type=Path,
                        default=ROOT/'f90src/Plant_bgc/PlantPhenolMod.F90')
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix='ecosim_phenology_phases_') as work:
        folder = Path(work)
        source = folder/'check.F90'
        source.write_text(fixture(args.phenology_source.read_text())+DRIVER)
        subprocess.run(shlex.split(args.fc)+['-ffree-line-length-none','-fcheck=all',
                       '-ffpe-trap=invalid,zero,overflow',
                       str(ROOT/'f90src/Modelconfig/ElmIDMod.F90'),str(source),
                       '-o',str(folder/'check')],cwd=folder,check=True)
        subprocess.run([str(folder/'check')],cwd=folder,check=True)


if __name__ == '__main__':
    main()
