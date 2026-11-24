module UptakePars
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSimConst
  use ElmIDMod
  implicit none
  public
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
!
!     MaxIterNum=max number of cycles in convergence soln for water uptake
!     DIFFX,DIFFY=acceptance criteria in convergence soln
!     FMN=min PFT:total population ratio
!     RawMinH2OCanopy2Atm_pft,RACX=min,max canopy boundary layer resistance [h m-1]
!     ResistLeafSurf=surface resistance to evaporation [h m-1]
!     EMODW=wood modulus of elasticity (MPa)
!     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition on root,myco N,P uptake(g g-1)
!     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation (h-1)
!

  integer :: MaxIterNum
  real(r8) :: DIFFX
  real(r8) :: DIFFY
  real(r8) :: FMN
  real(r8) :: RawMinH2OCanopy2Atm_pft !this value varies with vegetation type, forest (0.0014) < long grass (0.0083) < short grass (0.0278) < moss (0.0417), [h/m]
  real(r8) :: RACX
  real(r8) :: ResistLeafSurf
  real(r8) :: EMMC         !canopy emissivity of longwave radiation
  real(r8) :: EMODW
  real(r8) :: ZCKI
  real(r8) :: PCKI
  real(r8) :: ZPKI
  real(r8) :: PZKI
  real(r8) :: FEXUDE(NumPlantChemElms)

  contains
  subroutine InitUptakePars
  implicit none

  MaxIterNum              = 200
  DIFFX                   = 1.0E-09_r8
  DIFFY                   = 0.5E-02_r8
  FMN                     = 1.e-6_r8 !ppmc
  RawMinH2OCanopy2Atm_pft    = 0.00139_r8  ! [h/m]
  RACX                    = 0.0278_r8
  ResistLeafSurf          = 0.0139_r8

  EMMC          = 0.97_r8
  EMODW         = 50.0_r8
  ZCKI          = 0.5E-01_r8
  PCKI          = 0.5E-02_r8
  ZPKI          = ZCKI/PCKI
  PZKI          = PCKI/ZCKI
  FEXUDE(ielmc) = 0.5E-03_r8
  FEXUDE(ielmn) = 1.0E-02_r8
  FEXUDE(ielmp) = 1.0E-02_r8
  
  end subroutine InitUptakePars
end module UptakePars
