"""Exercise production fine-root carbon funding under nutrient and carbon stress.

Compile the actual Grow2ndRootAxes mobilization/growth block with isolated state.
Check maintenance priority, post-maintenance growth limitation, substrate caps,
oxygen/carbon shortages and saturation for both roots and mycorrhizae.
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fc', default=os.environ.get('FC', 'gfortran'))
    args = parser.parse_args()
    text = (ROOT / 'f90src/Plant_bgc/RootMod.F90').read_text()
    start = text.index('  ! Fine-root and mycorrhizal maintenance use potential C mobilization')
    stop = text.index('!     TOTAL NON-STRUCTURAL C,N,P USED IN SECONDARY ROOT GROWTH', start)
    block = text[start:stop]
    fixture = r'''
program fine_root_maintenance_test
  ! Description:
  ! Exercise production carbon mobilization and growth statements with fixed inputs.
  use iso_fortran_env, only: r8 => real64
  implicit none
  integer, parameter :: ipltroot=1,ielmc=1,ielmn=2,ielmp=3
  integer :: N=1,L=1,NZ=1,root_type
  real(r8) :: VMXC=.015_r8, VMXCPot,FineRootGrowthFactor
  real(r8) :: RootMycoNonstElms_rpvr(3,2,1,1),fTGrowRootP_vr(1,1)
  real(r8) :: ROOTNLim_rpvr(2,1,1),ROOTPLim_rpvr(2,1,1),GrainFillDowreg_brch(1,1)
  integer :: MainBranchNum_pft(1)=1
  real(r8) :: fRootGrowPSISense,FracRoot2ndCSinkL,RAutoRootO2Limter_rpvr(2,1,1)
  real(r8) :: RNonstCO2_OUltd,RNonstCO2_Oltd,RPotCO2Groth2nd_OUltd,RPotCO2Groth2nd_Oltd
  real(r8) :: Rmaint2nd_CO2,RGrowCO2_OUltd,RGrowCO2_Oltd,RespElongWatSens
  real(r8) :: RootMaintDef_CO2_pvr(2,1,1),ZPOOLB,PPOOLB,FNP,CNRTS_pft(1),CPRTS_pft(1),DMRespEff

  do root_type=1,2
  call reset()
  call calculate()
  call near(RNonstCO2_Oltd,.015_r8,'potential supply bypasses N/P growth restriction')
  call near(RPotCO2Groth2nd_Oltd,.014_r8,'maintenance paid before growth')
  call near(RGrowCO2_Oltd,.0007_r8,'N/P restricts only post-maintenance remainder')
  call near(RootMaintDef_CO2_pvr(N,1,1),0._r8,'nutrient shortage alone does not create maintenance deficit')

  call reset()
  ROOTNLim_rpvr=0._r8;ROOTPLim_rpvr=0._r8;RootMycoNonstElms_rpvr(2:3,:,:,:)=0._r8
  call calculate()
  call near(RNonstCO2_Oltd,.015_r8,'zero nutrients retain maintenance access to C')
  call near(RGrowCO2_Oltd,0._r8,'zero nutrients prevent growth')

  call reset()
  RootMycoNonstElms_rpvr(1,:,:,:)=.01_r8
  call calculate()
  call near(RNonstCO2_Oltd,.00015_r8,'limited reserve C bounds mobilization')
  call near(RootMaintDef_CO2_pvr(N,1,1),-.00085_r8,'real carbon shortage remains visible')

  call reset()
  RAutoRootO2Limter_rpvr=.02_r8
  call calculate()
  call near(RNonstCO2_Oltd,.0003_r8,'oxygen still restricts supply')
  call near(RGrowCO2_Oltd,0._r8,'oxygen deficit prevents growth')
  call near(RootMaintDef_CO2_pvr(N,1,1),-.0007_r8,'oxygen shortage remains visible')

  call reset()
  GrainFillDowreg_brch=0._r8
  call calculate()
  call near(RNonstCO2_Oltd,.015_r8,'grain-fill shutdown preserves maintenance supply')
  call near(RGrowCO2_Oltd,0._r8,'grain-fill shutdown prevents growth')

  call reset()
  ROOTNLim_rpvr=1._r8;ROOTPLim_rpvr=1._r8
  call calculate()
  call near(RGrowCO2_Oltd,.014_r8,'unlimited nutrients recover original growth allowance')

  call reset()
  CNRTS_pft=1000._r8;CPRTS_pft=1000._r8
  call calculate()
  call near(RGrowCO2_Oltd,RootMycoNonstElms_rpvr(3,N,1,1)/1000._r8*.5_r8, &
    'independent stoichiometric substrate cap remains active')

  call reset()
  fTGrowRootP_vr=1000._r8
  call calculate()
  call near(RNonstCO2_Oltd,1._r8,'mobilization cannot exceed available reserve C')

  enddo

  root_type=2
  call reset()
  RootMycoNonstElms_rpvr(:,1,:,:)=0._r8
  ROOTNLim_rpvr(1,:,:)=0._r8;ROOTPLim_rpvr(1,:,:)=0._r8
  call calculate()
  call near(RNonstCO2_Oltd,.015_r8,'fungi use their own carbon reserve')
  call near(RGrowCO2_Oltd,.0007_r8,'fungi use their own nutrient factors')
  print *, 'PASS: 17 production root/mycorrhizal maintenance/growth scenarios'
contains
  subroutine reset()
    ! Description:
    ! Initialize a nutrient-limited root with sufficient carbohydrate reserves.
    N=root_type;RootMycoNonstElms_rpvr(1,:,:,:)=1._r8
    RootMycoNonstElms_rpvr(2,:,:,:)=.02_r8
    RootMycoNonstElms_rpvr(3,:,:,:)=.01_r8/19._r8
    fTGrowRootP_vr=1._r8;fRootGrowPSISense=1._r8;FracRoot2ndCSinkL=1._r8
    ROOTNLim_rpvr=1._r8/6._r8;ROOTPLim_rpvr=.05_r8
    GrainFillDowreg_brch=1._r8;RAutoRootO2Limter_rpvr=1._r8
    Rmaint2nd_CO2=.001_r8;RespElongWatSens=1._r8;RootMaintDef_CO2_pvr=0._r8
    CNRTS_pft=.1_r8;CPRTS_pft=.01_r8;DMRespEff=.5_r8
  end subroutine reset
  subroutine near(actual,expected,label)
    ! Description:
    ! Compare a production result against an independently specified expectation.
    real(r8),intent(in) :: actual,expected ! Carbon amounts for the one-hour step, g C.
    character(len=*),intent(in) :: label ! Scientific behavior being checked.
    if(abs(actual-expected)>1.e-12_r8)then
      print *, 'FAIL: ',label,actual,expected
      error stop 1
    endif
  end subroutine near
  real(r8) function AZMAX1(value)
    ! Description:
    ! Supply the production nonnegative scalar operation for the isolated fixture.
    real(r8),intent(in) :: value ! Scalar to bound below by zero.
    AZMAX1=max(0._r8,value)
  end function AZMAX1
  subroutine calculate()
    ! Description:
    ! Run the unmodified statements extracted from production Grow2ndRootAxes.
''' + block + '''
  end subroutine calculate
end program fine_root_maintenance_test
'''
    with tempfile.TemporaryDirectory(prefix='ecosim_fine_root_') as tmp:
        path = Path(tmp)
        source = path / 'test.F90'
        source.write_text(fixture)
        exe = path / 'test.x'
        subprocess.run([args.fc, '-O0', '-g', '-fdefault-real-8', '-fdefault-double-8',
                        '-fcheck=all', '-ffpe-trap=invalid,zero,overflow',
                        '-ffree-line-length-none', str(source), '-o', str(exe)], check=True, cwd=path)
        subprocess.run([str(exe)], check=True, cwd=path)


if __name__ == '__main__':
    main()
