"""Check maintenance-first shoot respiration using complete production procedures.

Exercise both sides of emergence with nutrient shortage, carbon shortage, water,
oxygen, grain-fill shutdown and substrate caps; verify pool and respiration ledgers.
"""
import argparse
import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def procedure(text, name):
    match = re.search(r'^\s*(?:pure\s+)?(?:subroutine|function) ' + name +
                      r'\b.*?^\s*end (?:subroutine|function) ' + name + r'\b[^\n]*',
                      text, re.M | re.S | re.I)
    if match is None:
        raise ValueError(name)
    return match[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fc', default=os.environ.get('FC', 'gfortran'))
    args = parser.parse_args()
    source = (ROOT / 'f90src/Plant_bgc/PlantBranchMod.F90').read_text()
    routines = '\n'.join(procedure(source, n) for n in
                         ['ComputRAutoAfEmergence', 'ComputRAutoB4Emergence'])
    math = procedure((ROOT / 'f90src/Minimath/MiniMathMod.F90').read_text(), 'AZMAX1_s')
    bryophyte = procedure((ROOT / 'f90src/Plant_bgc/PlantMathFuncMod.F90').read_text(), 'is_root_bryophyte')
    pars = (ROOT / 'f90src/Modelpars/PlantBGCPars.F90').read_text()
    assignments = '\n'.join(re.search(r'^\s*' + n + r'\s*=.*$', pars, re.M)[0]
                            for n in ['VMXC', 'CNKI', 'CPKI', 'RmSpecPlant'])
    fixture = r'''
module branch_fixture
  ! Description:
  ! Provide isolated branch API state for the complete production procedures.
  use iso_fortran_env, only: r8=>real64
  use ElmIDMod
  implicit none
  integer,parameter :: JZ1=1
  real(r8),parameter :: tiny_val=1.e-50_r8
  real(r8) :: VMXC,CNKI,CPKI,RmSpecPlant
  type biomass_state
    real(r8) :: CanopyNonstElms_brch(3,1,1),LeafPetoNonstElmConc_brch(3,1,1)
  end type
  type site_state
    real(r8) :: ZERO=1.e-12_r8
  end type
  type phenology_state
    integer :: iPlantRootProfile_pft(1)=1,iPlantPhenolType_pft(1)=0
    real(r8) :: fTCanopyGroth_pft(1)=1._r8,fTgrowRootP_vr(1,1)=1._r8
  end type
  type photo_state
    real(r8) :: GrainFillDowreg_brch(1,1)=1._r8
  end type
  type morphology_state
    integer :: NGTopRootLayer_pft(1)=1
  end type
  type root_state
    real(r8) :: RootRespPotent_pvr(2,1,1)=0,RootCO2EmisPot_pvr(2,1,1)=0,RootCO2Autor_pvr(2,1,1)=0
    real(r8) :: RAutoRootO2Limter_rpvr(2,1,1)=1._r8
  end type
  type rates_state
    real(r8) :: CO2NetFix_pft(1)=0,CanopyGrosRCO2_pft(1)=0,CanopyResp_brch(1,1)=0
    real(r8) :: Eco_AutoR_CumYr_col=0,ECO_ER_col=0,GrossCO2Fix_pft(1)=0
    real(r8) :: RCanMaintDef_CO2_pft(1)=0,CanopyRespC_CumYr_pft(1)=0
    real(r8) :: CanopyNLimFactor_brch(1,1)=0,CanopyPLimFactor_brch(1,1)=0
  end type
  type(biomass_state) :: plt_biom
  type(site_state) :: plt_site
  type(phenology_state) :: plt_pheno
  type(photo_state) :: plt_photo
  type(morphology_state) :: plt_morph
  type(root_state) :: plt_rbgc
  type(rates_state) :: plt_bgcr
  interface AZMAX1
    module procedure AZMAX1_s
  end interface
contains
  subroutine PrintInfo(message)
    ! Description:
    ! Suppress model logging in the isolated test.
    character(len=*),intent(in) :: message ! Diagnostic message.
  end subroutine PrintInfo
''' + math + '\n' + bryophyte + '\n' + routines + r'''
end module branch_fixture
program test
  ! Description:
  ! Check maintenance, growth, deficit and C/N/P ledger invariants across emergence.
  use branch_fixture
  implicit none
  integer :: mode,cases=0
  real(r8) :: growth,deficit,cnpg,used(3),water,turg,structN,cns,cps
  real(r8) :: yield=.5_r8,tfn(1),autoresp,expected
''' + assignments + r'''
  do mode=1,2
    call reset()
    call run()
    call near(deficit,0._r8,'N/P shortage does not restrict maintenance funding')
    call near(growth*yield,.0007_r8,'N/P acts on post-maintenance growth')
    call near(cnpg,.05_r8,'nutrient multiplier remains available for leaf stoichiometry')

    call reset()
    plt_biom%CanopyNonstElms_brch(2:3,1,1)=0._r8
    plt_biom%LeafPetoNonstElmConc_brch(:,1,1)=plt_biom%CanopyNonstElms_brch(:,1,1)
    call run()
    call near(deficit,0._r8,'zero N/P still permits maintenance')
    call near(growth,0._r8,'zero N/P prevents new tissue')
    call near(used(ielmc),.001_r8,'only maintenance C is consumed without nutrients')

    call reset()
    plt_photo%GrainFillDowreg_brch=0._r8
    call run()
    call near(deficit,0._r8,'grain-fill shutdown preserves maintenance funding')
    call near(growth,0._r8,'grain-fill shutdown blocks growth')

    call reset()
    plt_biom%CanopyNonstElms_brch(ielmc,1,1)=.01_r8
    call run()
    call near(deficit,.00085_r8,'carbon shortage remains a maintenance deficit')
    call near(growth,0._r8,'carbon deficit blocks growth')

    call reset()
    plt_biom%CanopyNonstElms_brch(:,1,1)=0._r8
    plt_biom%LeafPetoNonstElmConc_brch(:,1,1)=0._r8
    call run()
    call near(deficit,.001_r8,'empty reserves preserve the carbon deficit')
    call near(used(ielmc),0._r8,'empty reserves cannot supply C')

    call reset()
    water=.5_r8
    call run()
    call near(growth*yield,.000325_r8,'water restriction remains on mobilization')

    call reset()
    plt_pheno%fTCanopyGroth_pft=.5_r8;plt_pheno%fTgrowRootP_vr=.5_r8
    call run()
    call near(growth*yield,.000325_r8,'temperature restriction remains on mobilization')

    call reset()
    turg=.5_r8
    call run()
    call near(growth*yield,.00035_r8,'turgor restriction remains on growth')

    call reset()
    plt_biom%LeafPetoNonstElmConc_brch(:,1,1)=0._r8
    call run()
    call near(growth*yield,.014_r8,'unrestricted nutrient fallback recovers original allowance')

    call reset()
    cns=1000._r8;cps=100._r8
    expected=plt_biom%CanopyNonstElms_brch(ielmp,1,1)*yield/(cps+.002_r8+.001_r8*.05_r8)
    call run()
    call near(growth*yield,expected,'stoichiometric N/P substrate cap remains active')

    call reset()
    plt_photo%GrainFillDowreg_brch=0._r8
    plt_biom%CanopyNonstElms_brch(ielmc,1,1)=.01_r8
    plt_pheno%fTCanopyGroth_pft=1000._r8;plt_pheno%fTgrowRootP_vr=1000._r8
    structN=20._r8
    call run()
    call near(used(ielmc),.01_r8,'mobilization is bounded by available C')
    call near(deficit,.19_r8,'mobilization cap preserves remaining maintenance demand')
  enddo

  mode=2
  call reset()
  plt_rbgc%RAutoRootO2Limter_rpvr=.02_r8
  call run()
  call near(deficit,.0007_r8,'pre-emergence oxygen shortage remains a deficit')
  call near(growth,0._r8,'oxygen deficit blocks growth')

  call reset()
  plt_rbgc%RAutoRootO2Limter_rpvr=.5_r8
  call run()
  call near(growth*yield,.000325_r8,'oxygen acts before maintenance subtraction')
  print *, 'PASS:',cases,'complete branch respiration scenarios, with C/N/P accounting'
contains
  subroutine reset()
    ! Description:
    ! Restore a nutrient-limited branch with sufficient carbohydrate reserves.
    plt_bgcr=rates_state();plt_rbgc=root_state()
    plt_biom%CanopyNonstElms_brch(:,1,1)=[1._r8,.02_r8,.01_r8/19._r8]
    plt_biom%LeafPetoNonstElmConc_brch=plt_biom%CanopyNonstElms_brch
    plt_pheno%fTCanopyGroth_pft=1._r8;plt_pheno%fTgrowRootP_vr=1._r8
    plt_photo%GrainFillDowreg_brch=1._r8;plt_rbgc%RAutoRootO2Limter_rpvr=1._r8
    water=1._r8;turg=1._r8;structN=.1_r8;cns=.1_r8;cps=.01_r8;tfn=1._r8
  end subroutine reset
  subroutine run()
    ! Description:
    ! Call a production routine and check elemental and respiration ledgers.
    if(mode==1)then
      call ComputRAutoAfEmergence(1,1,1,1,yield,.02_r8,.002_r8,cns,cps,.01_r8,.001_r8, &
        0._r8,0._r8,1._r8,water,turg,structN,used,cnpg,deficit,growth)
      autoresp=-plt_bgcr%CanopyResp_brch(1,1)
    else
      call ComputRAutoB4Emergence(1,1,1,1,tfn,yield,.02_r8,.002_r8,cns,cps,.01_r8,.001_r8, &
        structN,water,turg,used,cnpg,deficit,growth)
      autoresp=-plt_rbgc%RootCO2Autor_pvr(ipltroot,1,1)
    endif
    call near(used(ielmc),autoresp+growth*(1._r8-yield),'C withdrawal = respiration + structural growth')
    if(any(used< -1.e-12_r8) .or. any(used>plt_biom%CanopyNonstElms_brch(:,1,1)+1.e-12_r8)) &
      error stop 'Consumption exceeds available C/N/P'
    cases=cases+1
  end subroutine run
  subroutine near(actual,target,label)
    ! Description:
    ! Compare independent expectations with production results.
    real(r8),intent(in) :: actual,target ! Quantities on the same carbon or ratio basis.
    character(len=*),intent(in) :: label ! Scientific property under test.
    if(abs(actual-target)>1.e-12_r8)then
      print *, 'FAIL:',label,actual,target
      error stop 1
    endif
  end subroutine near
end program test
'''
    with tempfile.TemporaryDirectory(prefix='ecosim_branch_maint_') as tmp:
        path = Path(tmp)
        f = path / 'test.F90'
        f.write_text(fixture)
        cmd = [args.fc, '-cpp', '-O0', '-g', '-fdefault-real-8', '-fdefault-double-8',
               '-fcheck=all', '-ffpe-trap=invalid,zero,overflow', '-ffree-line-length-none',
               str(ROOT / 'f90src/Modelconfig/ElmIDMod.F90'), str(f), '-o', str(path / 'test.x')]
        subprocess.run(cmd, check=True, cwd=path)
        subprocess.run([str(path / 'test.x')], check=True, cwd=path)


if __name__ == '__main__':
    main()
