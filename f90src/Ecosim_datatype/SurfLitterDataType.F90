module SurfLitterDataType
  ! Own native surface-litter state, transfer coefficients and flux accumulators.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod, only: idg_beg,idg_NH3
  use EcoSIMConfig, only : jcplx=> jcplxc, jcplx1=>jcplxcm1
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8) ,target,allocatable ::   BulkDensLitR(:)                   !surface litter bulk density,	[Mg m-3]
  real(r8) ,target,allocatable ::   GasXferCoeffLitr2Atm_col(:,:,:,:) !litter-area/time-scaled gas coefficient (gas,M,NY,NX) [m3]
  integer  ,target,allocatable ::   iLitrType_col(:,:,:)              !surface litter type:1 = plant, 2 = manure,[-]
  real(r8) ,target,allocatable ::   XTillCorp_col(:,:)                !factor for surface litter incorporation and soil mixing
  real(r8) ,target,allocatable ::   WatFLoLitr2SoilM_col(:,:,:)        !water transfer between soil surface and surface litter, [g d-2 t-1]
  real(r8) ,target,allocatable ::   WatFLoLitr2MacPM_col(:,:,:)        !accepted litter-to-macropore water volume per soil substep, [m3 d-2]
  real(r8) ,target,allocatable ::   WatFlowSno2LitRM_col(:,:,:)       !meltwater flux into surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   FracSurfByLitR_col(:,:)           !fraction of soil surface covered by surface litter, [-]
  real(r8) ,allocatable ::   FracSurfByLitRM_col(:,:,:)               !litter fraction used by gas/liquid interface coefficients (M,NY,NX) [-]
  real(r8) ,target,allocatable ::   HeatFLoByWat2LitR_col(:,:)        !net heat transfer to surface litter, [MJ d-2 t-1]
  real(r8) ,target,allocatable ::   VLitR_col(:,:)                    !surface litter volume, [m3 d-2]
  real(r8) ,target,allocatable ::   VHeatCapLitRMin_col(:,:)          !threshold surface litter heat capacity, [MJ d-2 K-1]
  real(r8) ,target,allocatable ::   VWatLitRHoldCapcity_col(:,:)      !surface litter water holding capacity, [m3 d-2]
  real(r8) ,target,allocatable ::   WatFLo2LitR_col(:,:)              !net water transfer to surface litter, [MJ d-2 t-1]
  real(r8) ,target,allocatable ::   TLitrIceFlxThaw_col(:,:)          !water from ice thaw in surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   TLitrIceHeatFlxFrez_col(:,:)      !latent heat released from water freeze in surface litter, [MJ d-2 h-1]
  real(r8) ,target,allocatable ::   Rain2LitRSurf_col(:,:)            !precipitation flux into surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   Irrig2LitRSurf_col(:,:)               !irrigation flux into surface litter, [m3 d-2 h-1]
  real(r8) ,target,allocatable ::   POROS0_col(:,:)                   !litter porosity, [m3 pore m-3 litr]
  real(r8) ,target,allocatable ::   SurfLitterC_col(:,:,:)                        !surface litter OM in each complex,	[g d-2]
  real(r8) ,target,allocatable ::   RC0ff_col(:,:)                        !surface litter OM in the autotrophic complex,[g d-2]
  real(r8) ,target,allocatable ::   LitWatMassBeg_col(:,:)            !total inital water mass in litter layer, [m3 H2O d-2]
  real(r8) ,target,allocatable ::   LitWatMassEnd_col(:,:)            !total inital water mass in litter layer, [m3 H2O d-2]
  real(r8) ,target,allocatable ::   Rain2LitR_col(:,:)                !total precipiation reaches the litter layer, [m3 H3O d-2 h-1]
  private :: InitAllocate
  contains
!----------------------------------------------------------------------

  subroutine InitSurfLitter(NumOfLitrCmplxs)
  ! Initialize litter state and gas-transfer work arrays for the configured grid.

  implicit none
  integer, intent(in) :: NumOfLitrCmplxs

  call InitAllocate(NumOfLitrCmplxs)
  end subroutine InitSurfLitter
!----------------------------------------------------------------------

  subroutine InitAllocate(NumOfLitrCmplxs)
  ! Allocate and zero litter state, fluxes and species/substep gas coefficients.

  implicit none
  integer, intent(in) :: NumOfLitrCmplxs

  allocate(Rain2LitR_col(JY,JX)); Rain2LitR_col=0._r8
  allocate(LitWatMassBeg_col(JY,JX)); LitWatMassBeg_col = 0._r8
  allocate(LitWatMassEnd_col(JY,JX)); LitWatMassEnd_col = 0._r8
  allocate(BulkDensLitR(1:NumOfLitrCmplxs));    BulkDensLitR =0._r8
  allocate(GasXferCoeffLitr2Atm_col(idg_beg:idg_NH3,60,JY,JX)); GasXferCoeffLitr2Atm_col=0._r8
  allocate(iLitrType_col(2,JY,JX));      iLitrType_col=0
  allocate(XTillCorp_col(JY,JX));        XTillCorp_col=0._r8
  allocate(WatFLoLitr2SoilM_col(60,JY,JX));     WatFLoLitr2SoilM_col=0._r8
  allocate(WatFLoLitr2MacPM_col(60,JY,JX));     WatFLoLitr2MacPM_col=0._r8
  allocate(WatFlowSno2LitRM_col(60,JY,JX));     WatFlowSno2LitRM_col=0._r8
  allocate(FracSurfByLitR_col(JY,JX));         FracSurfByLitR_col=0._r8
  allocate(FracSurfByLitRM_col(60,JY,JX));     FracSurfByLitRM_col=0._r8
  allocate(HeatFLoByWat2LitR_col(JY,JX));        HeatFLoByWat2LitR_col=0._r8
  allocate(VLitR_col(JY,JX));         VLitR_col=0._r8
  allocate(VHeatCapLitRMin_col(JY,JX));       VHeatCapLitRMin_col=0._r8
  allocate(VWatLitRHoldCapcity_col(JY,JX));       VWatLitRHoldCapcity_col=0._r8
  allocate(WatFLo2LitR_col(JY,JX));         WatFLo2LitR_col=0._r8
  allocate(TLitrIceFlxThaw_col(JY,JX));        TLitrIceFlxThaw_col=0._r8
  allocate(TLitrIceHeatFlxFrez_col(JY,JX));       TLitrIceHeatFlxFrez_col=0._r8
  allocate(Rain2LitRSurf_col(JY,JX));        Rain2LitRSurf_col=0._r8
  allocate(Irrig2LitRSurf_col(JY,JX));        Irrig2LitRSurf_col=0._r8
  allocate(POROS0_col(JY,JX));       POROS0_col=0._r8
  allocate(SurfLitterC_col(1:NumOfLitrCmplxs,JY,JX));      SurfLitterC_col=0._r8
  allocate(RC0ff_col(JY,JX)); RC0ff_col=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSurfLitter
  ! Release native litter state, fluxes and gas-transfer work arrays.
  use abortutils, only : destroy
  implicit none

  call destroy(Rain2LitR_col)
  call destroy(LitWatMassBeg_col)
  call destroy(LitWatMassEnd_col)
  call destroy(BulkDensLitR)
  call destroy(GasXferCoeffLitr2Atm_col)
  call destroy(iLitrType_col)
  call destroy(XTillCorp_col)
  call destroy(WatFLoLitr2SoilM_col)
  call destroy(WatFLoLitr2MacPM_col)
  call destroy(WatFlowSno2LitRM_col)
  call destroy(FracSurfByLitR_col)
  call destroy(FracSurfByLitRM_col)
  call destroy(HeatFLoByWat2LitR_col)
  call destroy(VLitR_col)
  call destroy(VHeatCapLitRMin_col)
  call destroy(VWatLitRHoldCapcity_col)
  call destroy(WatFLo2LitR_col)
  call destroy(TLitrIceFlxThaw_col)
  call destroy(TLitrIceHeatFlxFrez_col)
  call destroy(Rain2LitRSurf_col)
  call destroy(Irrig2LitRSurf_col)
  call destroy(POROS0_col)
  call destroy(SurfLitterC_col)
  call destroy(RC0ff_col)

  end subroutine DestructSurfLitter

!----------------------------------------------------------------------
  pure subroutine LitterInfiltrationFractions(WaterMic,WaterMac,WaterRunoff,WaterDonor,WaterTol,MaxFraction,FracMic,FracMac)
  ! Share a litter tracer donor between infiltration routes. Reserve the
  ! existing runoff proposal's capped fraction; do not count incoming runon.
  real(r8), intent(in) :: WaterMic,WaterMac,WaterRunoff,WaterDonor,WaterTol,MaxFraction
  real(r8), intent(out) :: FracMic,FracMac
  real(r8) :: TotalWater,Fraction,RunoffFraction,Limit

  FracMic=0._r8; FracMac=0._r8
  TotalWater=MAX(0._r8,WaterMic)+MAX(0._r8,WaterMac)
  IF(WaterDonor.LE.WaterTol .OR. TotalWater.LE.0._r8)RETURN
  Limit=MIN(1._r8,MAX(0._r8,MaxFraction))
  RunoffFraction=MIN(Limit,MAX(0._r8,WaterRunoff)/WaterDonor)
  Fraction=MIN(Limit,TotalWater/WaterDonor,1._r8-RunoffFraction)
  FracMic=Fraction*MAX(0._r8,WaterMic)/TotalWater
  FracMac=Fraction*MAX(0._r8,WaterMac)/TotalWater
  end subroutine LitterInfiltrationFractions

end module SurfLitterDataType
