"""Exercise production boundary outlets with a shared, nearly exhausted donor.

Extract the private production procedures to avoid adding a test-only public API.
Requires a GNU-compatible FC. No model inputs or reference results are modified.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'f90src/HydroTherm/SoilPhys/WatsubMod.F90'


def procedure(source, name):
    match = re.search(rf'^  subroutine {name}\b.*?^  end subroutine {name}\b',
                      source.read_text(), re.M | re.S | re.I)
    if not match:
        raise ValueError(f'Missing {name} in {source}')
    return match[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fc', default=os.environ.get('FC', 'gfortran'))
    parser.add_argument('--source', type=Path, default=SOURCE)
    args = parser.parse_args()
    # The current budget helper also allows testing the old outlet procedures.
    code = procedure(SOURCE, 'RemainingLiquidAfterFlows') + '\n'
    code += '\n'.join(procedure(args.source, n) for n in
                      ('DischargeOverWaterTable', 'VertBoundaryDrainM'))
    fixture = '''
module fixture
  use iso_fortran_env, only: r8 => real64
  implicit none
  integer, parameter :: iWestEastDirection=1,iNorthSouthDirection=2,iVerticalDirection=3
  integer :: FlowDirIndicator_col(2,2)=1,NL_col(2,2)=2,NUM_col(2,2)=1
  real(r8) :: VLWatMicP1_vr(0:3,2,2),VLWatMacP1_vr(0:3,2,2)
  real(r8) :: WaterFlow2Micpt_3D(3,0:3,2,2),WaterFlow2Macpt_3D(3,0:3,2,2)
  real(r8) :: WaterFlow2MicptX_3D(3,0:3,2,2),HeatFlow2Soili_3D(3,0:3,2,2)
  real(r8) :: VLSoilPoreMicP_vr(0:3,2,2)=1._r8,ZEROS2(2,2)=1.e-16_r8
  real(r8) :: SLOPE_col(3,2,2)=1._r8,DLYR_3D(3,0:3,2,2)=1._r8
  real(r8) :: WaterTable2GroundSlope_col(2,2)=1._r8
  real(r8) :: SoilDepthMidLay_vr(0:3,2,2)=0._r8,WaterTableDepZi_col(2,2)=0._r8
  real(r8) :: PSISoilMatricPtmp_vr(0:3,2,2)=0._r8,PSI_at_Sautration_vr(0:3,2,2)=0._r8
  real(r8) :: SoilWatAirDry_vr(0:3,2,2)=0._r8,POROS_vr(0:3,2,2)=1._r8
  real(r8) :: VLSoilMicP_vr(0:3,2,2)=1._r8,HydroCond_3D(3,1,0:3,2,2)=1._r8
  real(r8) :: HydCondSoil_3D(3,0:3,2,2),AREA_3D(3,0:3,2,2)=1._r8
  real(r8) :: VLMacroPores1_vr(0:3,2,2)=1._r8,HydroCondMacP1_vr(0:3,2,2)=1._r8
  real(r8) :: TKSoil1_vr(0:3,2,2)=280._r8
  real(r8) :: QDrainM_col(2,2),QDrain_col(2,2),HeatDrain_col(2,2)
  real(r8) :: QDrain_cum_col(2,2),QDrainloss_vr(0:3,2,2)
  real(r8), parameter :: cpw=4.186_r8,dts_HeatWatTP=1._r8/30._r8
  real(r8), parameter :: dts_wat=0.33333_r8,mGravAccelerat=0.00981_r8
  logical :: use_ecosys_drainage=.false.
  interface AMAX1
    module procedure max2
  end interface
  interface AMIN1
    module procedure min2
  end interface
contains
  subroutine PrintInfo(text)
    character(*),intent(in) :: text
  end subroutine
  real(r8) function min2(a,b)
    real(r8),intent(in) :: a,b
    min2=min(a,b)
  end function
  real(r8) function max2(a,b)
    real(r8),intent(in) :: a,b
    max2=max(a,b)
  end function
  real(r8) function AZMAX1(a)
    real(r8),intent(in) :: a
    AZMAX1=max(0._r8,a)
  end function
  real(r8) function AZMIN1(a)
    real(r8),intent(in) :: a
    AZMIN1=min(0._r8,a)
  end function
  real(r8) function AZERO(a)
    real(r8),intent(in) :: a
    AZERO=a
  end function
  logical function isclose(a,b)
    real(r8),intent(in) :: a,b
    isclose=abs(a-b)<1.e-20_r8
  end function
  real(r8) function safe_adb(a,b)
    real(r8),intent(in) :: a,b
    safe_adb=a/b
  end function
  integer function getMoistK(a,b)
    real(r8),intent(in) :: a,b
    getMoistK=1
  end function
  integer function SoilFaceSlot(n,src,dst)
    integer,intent(in) :: n,src,dst
    SoilFaceSlot=src
    if(n==3)SoilFaceSlot=dst
  end function
  subroutine SoilFaceIncomingRange(n,l,ny,nx,first,last)
    integer,intent(in) :: n,l,ny,nx
    integer,intent(out) :: first,last
    ! Single-column fixture: incoming external/vertical faces retain layer L.
    first=l;last=l
  end subroutine
  subroutine reset()
    VLWatMicP1_vr=0._r8; VLWatMacP1_vr=0._r8
    WaterFlow2Micpt_3D=0._r8;WaterFlow2Macpt_3D=0._r8
    WaterFlow2MicptX_3D=0._r8;HeatFlow2Soili_3D=0._r8
    QDrainM_col=0._r8;QDrain_col=0._r8;HeatDrain_col=0._r8
    QDrain_cum_col=0._r8;QDrainloss_vr=0._r8
  end subroutine
''' + code + '''
end module
program check_outlets
  use fixture
  implicit none
  integer :: mode,n,side,tx,ty,tile
  real(r8) :: xn,mic,mac,before,expected
  do mode=1,5
    call reset()
    use_ecosys_drainage=mode==4
    VLWatMicP1_vr(1,1,1)=1.e-6_r8
    VLWatMacP1_vr(1,1,1)=2.5e-8_r8
    ! Include the already accepted vertical inflow and outflow in the budget.
    WaterFlow2Macpt_3D(3,1,1,1)=1.67e-8_r8
    WaterFlow2Macpt_3D(3,2,1,1)=0.2e-8_r8
    before=1.e-6_r8+2.5e-8_r8+1.67e-8_r8-0.2e-8_r8
    if(mode==5)then
      call reset()
      before=0._r8
    endif
    do n=1,2
      do side=1,2
        xn=real(2*side-3,r8)
        tx=1;ty=1
        if(side==1.and.n==1)tx=2
        if(side==1.and.n==2)ty=2
        if(mode==3)then
          call VertBoundaryDrainM(60,13,n,1,1,1,tx,ty,1,1._r8,xn,1._r8,1._r8)
        else
          do tile=1,2
            expected=1._r8
            if(mode==2)expected=1.e-10_r8
            call DischargeOverWaterTable(60,13,n,1,1,1,tx,ty,1,.true.,.true., &
              1._r8,1._r8,0._r8,xn,expected,0._r8)
          enddo
        endif
      enddo
    enddo
    call RemainingLiquidAfterFlows(1,1,1,mic,mac)
    if(abs(mic+mac+QDrain_col(1,1)-before)>1.e-18_r8) &
      error stop 'outlets overdraw the shared water inventory'
    if(abs(HeatDrain_col(1,1)-cpw*280._r8*QDrain_col(1,1))>1.e-15_r8) &
      error stop 'heat export does not match accepted liquid'
    if(abs(QDrainM_col(1,1)-QDrain_col(1,1))>1.e-18_r8.or. &
       abs(QDrain_cum_col(1,1)-QDrain_col(1,1))>1.e-18_r8.or. &
       abs(QDrainloss_vr(1,1,1)-QDrain_col(1,1))>1.e-18_r8)error stop 'drainage ledgers disagree'
    if(mode==2)then
      expected=16._r8*mGravAccelerat*1.e-10_r8/2._r8*dts_HeatWatTP
      if(abs(QDrain_col(1,1)-expected)>1.e-20_r8)error stop 'unlimited flux changed'
    endif
  enddo
  print *, 'Boundary drainage shared donor and heat budgets passed'
end program
'''
    with tempfile.TemporaryDirectory(prefix='ecosim-boundary-budget-') as work:
        path = Path(work)
        (path / 'test.f90').write_text(fixture)
        subprocess.run(shlex.split(args.fc) + ['-ffree-line-length-none', '-fcheck=all',
                       str(path / 'test.f90'), '-o', str(path / 'test.x')], cwd=path, check=True)
        subprocess.run([str(path / 'test.x')], cwd=path, check=True)


if __name__ == '__main__':
    main()
