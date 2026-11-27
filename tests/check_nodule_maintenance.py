"""Exercise complete canopy/root nodule procedures and their C/N/P ledgers.

Check maintenance priority under structural nutrient deficiency, fixation energy,
carbon/water/oxygen limits, substrate saturation, host exchange and initial infection.
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
    text = (ROOT / 'f90src/Plant_bgc/NoduleBGCMod.F90').read_text()
    routines = '\n'.join(procedure(text, n) for n in ['CanopyNoduleBiochemistry', 'RootNodulBiochemistry'])
    math = (ROOT / 'f90src/Minimath/MiniMathMod.F90').read_text()
    funcs = '\n'.join(procedure(math, n) for n in ['AZMAX1_s', 'AZERO', 'safe_adb'])
    plantmath = (ROOT / 'f90src/Plant_bgc/PlantMathFuncMod.F90').read_text()
    funcs += '\n' + '\n'.join(procedure(plantmath, n) for n in ['is_root_N2fix', 'is_canopy_N2fix'])
    pars = (ROOT / 'f90src/Modelpars/PlantBGCPars.F90').read_text()
    names = ['VMXO','SPNDLK','RmSpecPlant','EN2F','SPNDL','RCCZN','RCCYN','RCCXN','RCCQN',
             'CNKI','CPKI','CZKM','CPKM','NodulBiomCatInfection','CCNGB','CCNGR','FXRN']
    assignments = '\n'.join(re.search(r'^\s*' + n + r'\s*=.*$', pars, re.M)[0] for n in names)
    declarations = 'real(r8) :: ' + ','.join(n for n in names if n != 'FXRN') + ',FXRN(6)\n'
    fixture = r'''
module nodule_fixture
  ! Description:
  ! Isolate API state while compiling complete production nodule procedures.
  use iso_fortran_env, only: r8=>real64
  use ElmIDMod
  implicit none
  integer,parameter :: JZ1=1,jsken=1
  real(r8),parameter :: tiny_val=1.e-50_r8
''' + declarations + r'''
  type config_state
    integer :: k_fine_comp=1,ifoliar=1,iroot=2
  end type
  type site_state
    integer :: NU=1,MaxNumRootLays=1
    real(r8) :: AREA3(1)=1._r8,ZERO=1.e-12_r8
  end type
  type morphology_state
    integer :: iPlantNfixType_pft(1)=1,NMaxRootBotLayer_pft(1)=1
  end type
  type phenology_state
    real(r8) :: fTCanopyGroth_pft(1)=1._r8,fTgrowRootP_vr(1,1)=1._r8,fRootGrowPSISense_pvr(2,1,1)=1._r8
  end type
  type chemistry_state
    real(r8) :: PlantElmAllocMat4Litr(3,2,1,1)=1._r8
  end type
  type allometry_state
    real(r8) :: NoduGrowthYield_pft(1)=.5_r8,rNCNodule_pft(1)=.1_r8,rPCNoduler_pft(1)=.01_r8
  end type
  type biomass_state
    real(r8) :: CanopyLeafSheathC_brch(1,1)=10._r8,ZERO4Groth_pft(1)=1.e-12_r8,ZERO4LeafVar_pft(1)=1.e-12_r8
    real(r8) :: CanopyNonstElms_brch(3,1,1)=0,CanopyNodulNonstElms_brch(3,1,1)=0,CanopyNodulStrutElms_brch(3,1,1)=0
    real(r8) :: RootNodulNonstElms_rpvr(3,1,1)=0,RootNodulStrutElms_rpvr(3,1,1)=0
    real(r8) :: RootMycoNonstElms_rpvr(3,2,1,1)=0,PopuRootMycoC_pvr(2,1,1)=10._r8
  end type
  type rates_state
    real(r8) :: CanopyGrosRCO2_pft(1)=0,ECO_ER_col=0,CanopyRespC_CumYr_pft(1)=0,Eco_AutoR_CumYr_col=0,CO2NetFix_pft(1)=0
    real(r8) :: LitrfallElms_pvr(3,1,1,0:1,1)=0,NodulInfectElms_pft(3,1)=0,RNFixCO2_pft(1)=0
    real(r8) :: RCO2Nodule_pvr(1,1)=0,RootN2Fix_pvr(1,1)=0
  end type
  type root_state
    real(r8) :: CanopyN2Fix_pft(1)=0,RootN2Fix_pft(1)=0,RAutoRootO2Limter_rpvr(2,1,1)=1._r8
    real(r8) :: RootRespPotent_pvr(2,1,1)=0,RootCO2EmisPot_pvr(2,1,1)=0,RootCO2Autor_pvr(2,1,1)=0
  end type
  type(config_state) :: pltpar
  type(site_state) :: plt_site
  type(morphology_state) :: plt_morph
  type(phenology_state) :: plt_pheno
  type(chemistry_state) :: plt_soilchem
  type(allometry_state) :: plt_allom
  type(biomass_state) :: plt_biom
  type(rates_state) :: plt_bgcr
  type(root_state) :: plt_rbgc
  interface AZMAX1
    module procedure AZMAX1_s
  end interface
contains
  subroutine PrintInfo(message)
    ! Description:
    ! Suppress model logging in the isolated test.
    character(len=*),intent(in) :: message ! Diagnostic message.
  end subroutine PrintInfo
''' + funcs + '\n' + routines + r'''
end module nodule_fixture
program test
  ! Description:
  ! Check both nodule locations, fixed-N energy, host transfer and elemental closure.
  use nodule_fixture
  implicit none
  integer :: mode,cases=0
  real(r8) :: before(3),after(3),litter(3),fixN,resp,water,tfn(1),b0,expected,energy,need,factor,concC
''' + assignments + r'''
  do mode=1,2
    call reset()
    call run()
    call near(fixN,0._r8,'N-saturated nodules do not fix N2')

    call reset()
    call structural_P(0._r8)
    call run()
    call near(structural_C(),1._r8-SPNDL,'P restriction does not create maintenance senescence')
    call near(fixN,0._r8,'zero structural P still blocks growth/fixation activity')

    call reset()
    call structural_P(1.e-8_r8)
    call run()
    if(structural_C()<1._r8-SPNDL-1.e-12_r8)error stop 'Low structural P restricted maintenance'

    call reset()
    call structural_N(.025_r8)
    concC=1._r8
    factor=.5_r8
    energy=(VMXO*concC/(concC+SPNDLK)-RmSpecPlant*.025_r8)*factor
    need=(.1_r8-.025_r8)/EN2F
    expected=EN2F*energy*need/(energy+need)
    call run()
    call near(fixN,expected,'fixation retains its demand and carbon energy controls')

    call reset()
    call reserve_C(0._r8)
    call run()
    if(structural_C()>=1._r8-SPNDL)error stop 'Real carbon shortage must retain maintenance senescence'
    call near(fixN,0._r8,'no carbon supply means no fixation')

    call reset()
    call reserve_C(1.e-6_r8)
    call run()
    if(structural_C()>=1._r8-SPNDL)error stop 'Carbon-limited substrate gate must remain'

    call reset()
    water=0._r8;plt_pheno%fRootGrowPSISense_pvr=0._r8
    call run()
    if(structural_C()>=1._r8)error stop 'Water-limited supply must retain maintenance deficit'

    call reset()
    plt_pheno%fTCanopyGroth_pft=0._r8;plt_pheno%fTgrowRootP_vr=0._r8
    call run()
    if(structural_C()>=1._r8)error stop 'Cold-limited mobilization must remain'

    call reset()
    plt_pheno%fTCanopyGroth_pft=100._r8;plt_pheno%fTgrowRootP_vr=100._r8
    call run()

    call reset()
    plt_morph%iPlantNfixType_pft=iN2fixtyp_none
    call run()
    call near(resp,0._r8,'nonfixers do not run nodule metabolism')
    call near(fixN,0._r8,'nonfixers do not acquire fixed N')

    call reset()
    plt_biom%CanopyNodulStrutElms_brch=0._r8;plt_biom%RootNodulStrutElms_rpvr=0._r8
    call run()
    call near(plt_bgcr%NodulInfectElms_pft(ielmc,1),NodulBiomCatInfection,'initial infection source remains accounted')
  enddo
  mode=2
  call reset()
  plt_rbgc%RAutoRootO2Limter_rpvr=.001_r8
  call run()
  if(structural_C()>=1._r8-SPNDL)error stop 'Root nodule oxygen shortage must remain a deficit'

  call reset()
  call structural_N(.025_r8)
  plt_rbgc%RAutoRootO2Limter_rpvr=.5_r8
  concC=1._r8;factor=.5_r8
  energy=(.5_r8*VMXO*concC/(concC+SPNDLK)-RmSpecPlant*.025_r8)*factor
  need=(.1_r8-.025_r8)/EN2F
  expected=EN2F*energy*need/(energy+need)
  call run()
  call near(fixN,expected,'root oxygen restricts supply before maintenance and N fixation')
  print *, 'PASS:',cases,'complete nodule scenarios with C/N/P closure'
contains
  subroutine reset()
    ! Description:
    ! Restore an active nodule and its host before each independent test.
    plt_biom=biomass_state();plt_bgcr=rates_state();plt_rbgc=root_state();plt_pheno=phenology_state()
    plt_morph%iPlantNfixType_pft=merge(in2fixtyp_canopy_fast,in2fixtyp_root_fast,mode==1)
    plt_biom%CanopyNonstElms_brch(:,1,1)=[1._r8,.05_r8,.005_r8]
    plt_biom%RootMycoNonstElms_rpvr(:,ipltroot,1,1)=[1._r8,.05_r8,.005_r8]
    plt_biom%CanopyNodulNonstElms_brch(:,1,1)=[1._r8,.02_r8,.002_r8]
    plt_biom%RootNodulNonstElms_rpvr(:,1,1)=[1._r8,.02_r8,.002_r8]
    plt_biom%CanopyNodulStrutElms_brch(:,1,1)=[1._r8,.1_r8,.01_r8]
    plt_biom%RootNodulStrutElms_rpvr(:,1,1)=[1._r8,.1_r8,.01_r8]
    water=1._r8;tfn=1._r8
  end subroutine reset
  subroutine structural_P(value)
    ! Description:
    ! Set nodule structural P without changing its carbon reserve.
    real(r8),intent(in) :: value ! Structural P, g P.
    plt_biom%CanopyNodulStrutElms_brch(ielmp,1,1)=value
    plt_biom%RootNodulStrutElms_rpvr(ielmp,1,1)=value
  end subroutine structural_P
  subroutine structural_N(value)
    ! Description:
    ! Set structural N to test the independent demand for nitrogen fixation.
    real(r8),intent(in) :: value ! Structural N, g N.
    plt_biom%CanopyNodulStrutElms_brch(ielmn,1,1)=value
    plt_biom%RootNodulStrutElms_rpvr(ielmn,1,1)=value
  end subroutine structural_N
  subroutine reserve_C(value)
    ! Description:
    ! Change nodule carbohydrate availability.
    real(r8),intent(in) :: value ! Nonstructural C, g C.
    plt_biom%CanopyNodulNonstElms_brch(ielmc,1,1)=value
    plt_biom%RootNodulNonstElms_rpvr(ielmc,1,1)=value
  end subroutine reserve_C
  real(r8) function structural_C()
    ! Description:
    ! Return structural C for the active nodule location.
    if(mode==1)then
      structural_C=plt_biom%CanopyNodulStrutElms_brch(ielmc,1,1)
    else
      structural_C=plt_biom%RootNodulStrutElms_rpvr(ielmc,1,1)
    endif
  end function structural_C
  function total_elements()result(total)
    ! Description:
    ! Include both partners so host-nodule transfers cancel in elemental budgets.
    real(r8) :: total(3) ! Host reserve plus nodule structural and reserve C/N/P, g.
    if(mode==1)then
      total=plt_biom%CanopyNonstElms_brch(:,1,1)+plt_biom%CanopyNodulStrutElms_brch(:,1,1) &
        +plt_biom%CanopyNodulNonstElms_brch(:,1,1)
    else
      total=plt_biom%RootMycoNonstElms_rpvr(:,ipltroot,1,1)+plt_biom%RootNodulStrutElms_rpvr(:,1,1) &
        +plt_biom%RootNodulNonstElms_rpvr(:,1,1)
    endif
  end function total_elements
  subroutine run()
    ! Description:
    ! Execute a production procedure and require C/N/P budgets to close.
    before=total_elements()
    if(mode==1)then
      call CanopyNoduleBiochemistry(1,1,1,1,1._r8,water)
      fixN=plt_rbgc%CanopyN2Fix_pft(1)
    else
      call RootNodulBiochemistry(1,1,1,tfn)
      fixN=plt_rbgc%RootN2Fix_pft(1)
    endif
    after=total_elements();resp=plt_bgcr%RNFixCO2_pft(1)
    litter=sum(sum(sum(sum(plt_bgcr%LitrfallElms_pvr,dim=5),dim=4),dim=3),dim=2)
    call near(after(ielmc)+litter(ielmc)+resp,before(ielmc)+plt_bgcr%NodulInfectElms_pft(ielmc,1),'C = pools + litter + respiration')
    call near(after(ielmn)+litter(ielmn),before(ielmn)+fixN+plt_bgcr%NodulInfectElms_pft(ielmn,1),'N gains equal fixation and infection')
    call near(after(ielmp)+litter(ielmp),before(ielmp)+plt_bgcr%NodulInfectElms_pft(ielmp,1),'P conserved through growth and turnover')
    if(minval(after)< -1.e-12_r8 .or. resp<0._r8 .or. fixN<0._r8)error stop 'Invalid elemental budget'
    if(mode==1)then
      if(minval(plt_biom%CanopyNodulStrutElms_brch)< -1.e-12_r8 .or. &
        minval(plt_biom%CanopyNodulNonstElms_brch)< -1.e-12_r8 .or. &
        minval(plt_biom%CanopyNonstElms_brch)< -1.e-12_r8)error stop 'Negative canopy pools'
    else
      if(minval(plt_biom%RootNodulStrutElms_rpvr)< -1.e-12_r8 .or. &
        minval(plt_biom%RootNodulNonstElms_rpvr)< -1.e-12_r8 .or. &
        minval(plt_biom%RootMycoNonstElms_rpvr)< -1.e-12_r8)error stop 'Negative root pools'
      if(plt_rbgc%RootRespPotent_pvr(ipltroot,1,1)+1.e-12_r8<resp)error stop 'Actual respiration exceeds oxygen-unlimited potential'
    endif
    cases=cases+1
  end subroutine run
  subroutine near(actual,target,label)
    ! Description:
    ! Compare production output with an independent budget or expected behavior.
    real(r8),intent(in) :: actual,target ! Quantities with matching elemental units.
    character(len=*),intent(in) :: label ! Tested physical property.
    if(abs(actual-target)>1.e-12_r8)then
      print *, 'FAIL:',mode,cases,label,actual,target
      error stop 1
    endif
  end subroutine near
end program test
'''
    with tempfile.TemporaryDirectory(prefix='ecosim_nodule_maint_') as tmp:
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
