module HydrologyDiagMod  
  ! Description:
  ! Diagnose water/ice contents, connected water-table depth and active-layer depth
  ! from model state without changing water inventories.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: isclose, AZMAX1, AZMIN1, safe_adb
  use EcosimConst,   only: mGravAccelerat
  use DebugToolMod  
  use SoilWaterDataType
  use GridDataType  
  use EcoSIMCtrlDataType
  use FlagDataType
  use SoilPropertyDataType
  use SoilPhysDataType
  implicit none
  private

  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: UpdateSoilWaterContents
  public :: ComputeWaterTableEquilContent
  public :: DiagWaterTableDepZi
  public :: DiagActiveLayerDepth

  contains
!------------------------------------------------------------------------------------------

  subroutine ComputeWaterTableEquilContent(L,NY,NX,Potential,MaxWaterContent,WaterContent)
  ! Description:
  ! Preserve the equilibrium-content approximation used by water-table diagnosis.
  ! It extrapolates the log-linear wet relation without SRP or a dry-branch switch;
  ! replacing it with a full retention inverse would change existing diagnostics.
  implicit none
  integer, intent(in) :: L       ! soil layer above the diagnosed saturated zone [-]
  integer, intent(in) :: NY,NX   ! row and column indices [-]
  real(r8), intent(in) :: Potential       ! equilibrium matric potential [MPa]
  real(r8), intent(in) :: MaxWaterContent ! liquid capacity after ice/air allowance [m3 m-3 soil]
  real(r8), intent(out) :: WaterContent   ! bounded equilibrium liquid content [m3 m-3 soil]
  character(len=*), parameter :: subname='ComputeWaterTableEquilContent'

  call PrintInfo('beg '//subname)
  IF(Potential.GE.-EXP(LOGPSIAtSat(NY,NX)))THEN
    WaterContent=MaxWaterContent
  ELSE
    WaterContent=AMIN1(MaxWaterContent,EXP((LOGPSIAtSat(NY,NX)-LOG(-Potential)) &
      *PSD_vr(L,NY,NX)/LOGPSIMXD_col(NY,NX)+LOGPOROS_vr(L,NY,NX)))
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine ComputeWaterTableEquilContent
!------------------------------------------------------------------------------------------

  subroutine UpdateSoilWaterContents(NY,NX)
  ! Update diagnostic water, ice and air contents without changing inventories.
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L
  character(len=*), parameter :: subname='UpdateSoilWaterContents'

  call PrintInfo('beg '//subname)
  ! Refresh fractions from current inventories before hydraulic pressures are computed.
  DO L=NU_col(NY,NX),NL_col(NY,NX)
    IF(VLSoilMicP_vr(L,NY,NX).LE.ZEROS(NY,NX))THEN
      THETW_vr(L,NY,NX)    = POROS_vr(L,NY,NX)
      THETI_vr(L,NY,NX)    = 0._r8
      ThetaAir_vr(L,NY,NX) = 0._r8
    ELSE
      THETW_vr(L,NY,NX)    = AZMAX1(AMIN1(POROS_vr(L,NY,NX),VLWatMicP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)))
      THETI_vr(L,NY,NX)    = AZMAX1(AMIN1(POROS_vr(L,NY,NX),VLiceMicP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)))
      ThetaAir_vr(L,NY,NX) = AZMAX1(VLsoiAirP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX))
    ENDIF

  ENDDO
  call PrintInfo('end '//subname)
  end subroutine UpdateSoilWaterContents
!------------------------------------------------------------------------------------------

  subroutine DiagWaterTableDepZi(I,J,NY,NX)
  ! Description:
  !Diagnose water table depth in the soil column
  implicit none
  integer, intent(in) :: I,J,NY,NX

  real(r8) :: THETAirP_vr(JZ)  !air-filled porosity in layer
  real(r8) :: PSIEquil         !hydrostatic equilibrium matric potential [MPa]
  real(r8) :: THETW1
  real(r8) :: THETWM
  real(r8) :: THETPX
  real(r8), parameter :: THETPW =1.e-3_r8   !maximum air-filled porosity for saturation (m3 m-3)
  integer :: LL,L
  logical :: FoundWaterTable
  character(len=*), parameter :: subname='DiagWaterTableDepZi'

!     begin_execution
  call PrintInfo('beg '//subname)

  ! Call after UpdateSoilWaterContents and GetSoilHydraulicVars so fractions
  ! and matric pressures describe the same current soil state.
  DO L=NU_col(NY,NX),NL_col(NY,NX)
    THETAirP_vr(L)=AZMAX1(POROS_vr(L,NY,NX)-THETW_vr(L,NY,NX)-THETI_vr(L,NY,NX))
  ENDDO

  IF(IDWaterTable_col(NY,NX).NE.0)THEN
    ! If no connected saturated zone is found, estimate a below-profile table
    ! from bottom-layer pressure, assuming psi=0 at the table and hydrostatic
    ! equilibrium. This is a diagnostic estimate, not a water transfer.
    L=NL_col(NY,NX)
    WaterTableDepZi_col(NY,NX)=AMAX1(CumDepz2LayBottom_vr(L,NY,NX), &
      SoilDepthMidLay_vr(L,NY,NX)-PSISoilMatricP_vr(L,NY,NX)/mGravAccelerat)

    DO L=NU_col(NY,NX),NL_col(NY,NX)
      IF(THETAirP_vr(L).GT.THETPW)cycle
      FoundWaterTable=.true.
      IF(SoilDepthMidLay_vr(L,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN
        ! Reject a perched candidate separated from the external table (or
        ! profile bottom) by an unsaturated layer, including a dry bottom.
        DO LL=L+1,NL_col(NY,NX)
          IF(THETAirP_vr(LL).GT.THETPW)THEN
            FoundWaterTable=.false.
            exit
          ENDIF
          IF(SoilDepthMidLay_vr(LL,NY,NX).GE.ExtWaterTable_col(NY,NX))exit
        ENDDO
      ENDIF
      IF(.not.FoundWaterTable)cycle

      IF(L.GT.NU_col(NY,NX))THEN
        ! Interpolate in the layer immediately above the accepted saturated
        ! layer. Use the same absolute air tolerance as the detection test;
        ! ice occupies pores but is not part of the liquid-water threshold.
        LL=L-1
        PSIEquil=PSISoilMatricP_vr(L,NY,NX)-mGravAccelerat* &
          (SoilDepthMidLay_vr(L,NY,NX)-SoilDepthMidLay_vr(LL,NY,NX))
        THETWM=AZMAX1(POROS_vr(LL,NY,NX)-THETI_vr(LL,NY,NX)-THETPW)
        call ComputeWaterTableEquilContent(LL,NY,NX,PSIEquil,THETWM,THETW1)
        THETPX=0._r8
        IF(THETWM.GT.THETW1)THEN
          THETPX=AMIN1(1._r8,AZMAX1((THETWM-THETW_vr(LL,NY,NX))/(THETWM-THETW1)))
        ENDIF
        WaterTableDepZi_col(NY,NX)=CumDepz2LayBottom_vr(LL,NY,NX) &
          -DLYR_3D(3,LL,NY,NX)*(1._r8-THETPX)
      ELSE
        WaterTableDepZi_col(NY,NX)=CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)
      ENDIF
      exit
    ENDDO
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine DiagWaterTableDepZi
!------------------------------------------------------------------------------------------

  subroutine DiagActiveLayerDepth(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  real(r8) :: VLiceTot,VLPoreTot,VLiceTotL,VOLWTL,VLPoreTotL
  integer :: LL,L
  logical :: FoundActiveLayer
  logical :: goto5701
  character(len=*), parameter :: subname='DiagActiveLayerDepth'

!     begin_execution
  call PrintInfo('beg '//subname)
  ! No frozen barrier unless the current profile establishes one below.
  ! Reset each hour so complete thaw cannot retain an earlier active-layer depth.
  ActiveLayDepZ_col(NY,NX)=9999._r8
  FoundActiveLayer = .false.
  DO L=NUI_col(NY,NX),NLI_col(NY,NX)
!
!     VOLI,VOLIH=ice volume in micropores,macropores
!     VOLW,VOLWH=water volume in micropores,macropores
!     VOLA,VOLAH=total volume in micropores,macropores
!     ActiveLayDepZ_col=active layer depth
!     CDPTH,DLYR=depth to bottom,thickness of soil layer
!
    IF(FoundActiveLayer)exit
    !current layer
    VLiceTot=VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)
    VLPoreTot=VLMicP_vr(L,NY,NX)+VLMacroPores_vr(L,NY,NX)

    IF(VLPoreTot.GT.ZEROS2(NY,NX).AND.VLiceTot.GT.0.01*VLPoreTot)THEN
      !significant ice and pore
      D5700: DO LL=MIN(L+1,NL_col(NY,NX)),NL_col(NY,NX)
        VLiceTotL  = VLiceMicP_vr(LL,NY,NX)+VLiceMacP_vr(LL,NY,NX)
        VOLWTL     = VLWatMicP_vr(LL,NY,NX)+VLWatMacP_vr(LL,NY,NX)
        VLPoreTotL = VLMicP_vr(LL,NY,NX)+VLMacroPores_vr(LL,NY,NX)
        !defined as significant increase in ice content
        goto5701=(VLPoreTotL.GT.ZEROS2(NY,NX).AND.VLiceTotL.LT.ZERO2*VLPoreTotL)
        if(goto5701)exit
      ENDDO D5700

      if(.not. goto5701)then
        IF(VLPoreTot.GT.ZEROS2(NY,NX))THEN
          ActiveLayDepZ_col(NY,NX)=CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)*AMIN1(1.0_r8,VLiceTot/VLPoreTot)
        ELSE
          ActiveLayDepZ_col(NY,NX)=CumDepz2LayBottom_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)
        ENDIF
        FoundActiveLayer=.true.
      else
        !not active layer
        ActiveLayDepZ_col(NY,NX)=9999.0_r8
      endif
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine DiagActiveLayerDepth
!------------------------------------------------------------------------------------------

end module HydrologyDiagMod
