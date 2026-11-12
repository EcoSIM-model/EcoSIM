module SnowBalanceMod
  ! Preserve snow mass and energy through hourly updates, removal and relayering.
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use data_const_mod,   only: spval => DAT_CONST_SPVAL
  use abortutils,       only: endrun
  use SnowTransportMod, only: SoluteTransportThruSnow
  use EcoSIMCtrlMod,    only: lverb,  snowRedist_model, fixWaterLevel
  use minimathmod,      only: AZMAX1, isclose, AZMIN1, AZMAX1d, AZERO
  use MiniFuncMod,      only: gOC_to_m3_OM
  use DebugToolMod
  use SoilPropertyDataType
  use SurfLitterDataType
  use SOMDataType
  USE SnowDataType
  use GridConsts
  use GridDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use SoilHeatDataType
  use ClimForcDataType
  use EcosimConst
  USE SoilPhysDataType
  USE EcoSimSumDataType
  USE FlagDataType
  use SnowPhysData
  use AqueChemDatatype
  use SoilBGCDataType
  use ChemTranspDataType

  use UnitMod, only : units
implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: SnowMassUpdate     !called in redist
  public :: SnowpackLayering   !called in redist after SnowMassUpdate
  public :: ZeroSnowArrays
  public :: DebugSnowPrint
  contains

  subroutine DebugSnowPrint(I,J,NY,NX,info)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  character(len=*), intent(in) :: info
  integer :: L
    if(I==152 .and. J==1 .and. NX==7)then
    L=1
    write(113,*)trim(info)//'snoadd',L,CumSno2SnowL_snvr(L,NY,NX),CumWat2SnowL_snvr(L,NY,NX), &
      CumIce2SnowL_snvr(L,NY,NX),CumHeat2SnowL_snvr(L,NY,NX)
    endif

  end subroutine DebugSnowPrint
!------------------------------------------------------------------------------------------

  subroutine SnowMassUpdate(I,J,NY,NX,QWatinfl2Mic,QHeatInfl2Soil)
  !
  !Description: update snow mass
  !
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8),optional,intent(out) :: QWatinfl2Mic
  real(r8),optional,intent(out) :: QHeatInfl2Soil
  character(len=*), parameter :: subname='SnowMassUpdate'
  real(r8) :: QWatinfl2Mic_loc,QHeatInfl2Soil_loc
  real(r8) :: VOLSWI,ENGYW,TKSnow_pre
  real(r8) :: ENGY,TKSX
  integer :: L,L2
  logical :: test_exist
      !
      !     CALCULATE SNOWPACK TEMPERATURE FROM ITS CHANGE
      !     IN HEAT STORAGE
      !
  call PrintInfo('beg '//subname)

  QSnowH2Oloss_col(NY,NX) = QSnowH2Oloss_col(NY,NX)-TDrysnoByRedist_col(NY,NX)  &
    -TIceSnoByRedist_col(NY,NX)*DENSICE-TWatSnoByRedist_col(NY,NX)
  QSnowHeatLoss_col(NY,NX)= QSnowHeatLoss_col(NY,NX)-THeatSnoByRedist_col(NY,NX)

  !update snow layer from top to bottom
  VOLSWI     = 0.0_r8
  TKSnow_pre = TKSnow_snvr(1,NY,NX)
!  if(TKSnow_snvr(1,NY,NX)<230._r8)then
!    write(*,*)I,J,TKSnow_pre,TairK_col(NY,NX)
!    call endrun(trim(mod_filename)//' at line',__LINE__)
!  endif

  test_exist=present(QWatinfl2Mic).OR.present(QHeatInfl2Soil)

  D9780: DO L=1,JS
    L2 = MIN(JS,L+1)
    IF(abs(SnoXfer2SnoLay_snvr(L,NY,NX))>0._r8)THEN
      if(L .LT. JS)then
        CumSno2SnowL_snvr(L,NY,NX)  = CumSno2SnowL_snvr(L,NY,NX)+SnoXfer2SnoLay_snvr(L,NY,NX)-SnoXfer2SnoLay_snvr(L2,NY,NX)
        CumWat2SnowL_snvr(L,NY,NX)  = CumWat2SnowL_snvr(L,NY,NX)+WatXfer2SnoLay_snvr(L,NY,NX)-WatXfer2SnoLay_snvr(L2,NY,NX)
        CumIce2SnowL_snvr(L,NY,NX)  = CumIce2SnowL_snvr(L,NY,NX)+IceXfer2SnoLay_snvr(L,NY,NX)-IceXfer2SnoLay_snvr(L2,NY,NX)
        CumHeat2SnowL_snvr(L,NY,NX) = CumHeat2SnowL_snvr(L,NY,NX)+HeatXfer2SnoLay_snvr(L,NY,NX)-HeatXfer2SnoLay_snvr(L2,NY,NX)
      ELSE
        CumSno2SnowL_snvr(L,NY,NX)  = CumSno2SnowL_snvr(L,NY,NX)+SnoXfer2SnoLay_snvr(L,NY,NX)
        CumWat2SnowL_snvr(L,NY,NX)  = CumWat2SnowL_snvr(L,NY,NX)+WatXfer2SnoLay_snvr(L,NY,NX)
        CumIce2SnowL_snvr(L,NY,NX)  = CumIce2SnowL_snvr(L,NY,NX)+IceXfer2SnoLay_snvr(L,NY,NX)
        CumHeat2SnowL_snvr(L,NY,NX) = CumHeat2SnowL_snvr(L,NY,NX)+HeatXfer2SnoLay_snvr(L,NY,NX)
      ENDIF
    ENDIF

    call UpdateSnowLayerL(I,J,L,NY,NX,VOLSWI)

!    call SoluteTransportThruSnow(I,J,L,NY,NX)

  ENDDO D9780



  VcumDrySnoWE_col(NY,NX) = sum(VLDrySnoWE_snvr(1:JS,NY,NX))
  VcumWatSnow_col(NY,NX)  = sum(VLWatSnow_snvr(1:JS,NY,NX))
  VcumIceSnow_col(NY,NX)  = sum(VLIceSnow_snvr(1:JS,NY,NX))
  VcumSnoDWI_col(NY,NX)   = sum(VLSnoDWIprev_snvr(1:JS,NY,NX))
  SnowDepth_col(NY,NX)    = sum(SnowThickL_snvr(1:JS,NY,NX))
  VcumSnowWE_col(NY,NX)   = VcumDrySnoWE_col(NY,NX)+VcumIceSnow_col(NY,NX)*DENSICE+VcumWatSnow_col(NY,NX)

!
! IF SNOWPACK DISAPPEARS

  !if(NY.EQ.1)then
  !  write(*,*) "After summation: "
  !  write(*,*) "   VcumDrySnoWE_col(NY,NX): ", VcumDrySnoWE_col(NY,NX)
  !  write(*,*) "   VcumWatSnow_col(NY,NX):  ", VcumWatSnow_col(NY,NX)
  !  write(*,*) "   VcumIceSnow_col(NY,NX):  ", VcumIceSnow_col(NY,NX)
  !  write(*,*) "   VcumSnoDWI_col(NY,NX):   ", VcumSnoDWI_col(NY,NX)
  !  write(*,*) "   SnowDepth_col(NY,NX):    ", SnowDepth_col(NY,NX)
  !  write(*,*) "   VcumSnowWE_col(NY,NX):   ", VcumSnowWE_col(NY,NX)
  !endif
!
! IF SNOWPACK DISAPPEARS

! intermediate disappearance, jyt commented out the following, Mar 14, 2026
!  IF(SoilBulkDensity_vr(NUM_col(NY,NX),NY,NX).LE.ZERO .or. SoilOrgM_vr(ielmc,0,NY,NX)<=1.e-2_r8)THEN
!    VLWatMicP_vr(NUM_col(NY,NX),NY,NX) = VLWatMicP_vr(NUM_col(NY,NX),NY,NX)+QSnoWatXfer2Soil_col(NY,NX)
!    VLiceMicP_vr(NUM_col(NY,NX),NY,NX) = VLiceMicP_vr(NUM_col(NY,NX),NY,NX)+QSnoIceXfer2Soil_col(NY,NX)

!    TKSX = TKS_vr(NUM_col(NY,NX),NY,NX)
!    ENGY = VHeatCapacity_vr(NUM_col(NY,NX),NY,NX)*TKSX

!    VHeatCapacity_vr(NUM_col(NY,NX),NY,NX) = VHeatCapacitySoilM_vr(NUM_col(NY,NX),NY,NX) &
!      +cpw*(VLWatMicP_vr(NUM_col(NY,NX),NY,NX)+VLWatMacP_vr(NUM_col(NY,NX),NY,NX)) &
!      +cpi*(VLiceMicP_vr(NUM_col(NY,NX),NY,NX)+VLiceMacP_vr(NUM_col(NY,NX),NY,NX))

!    IF(VHeatCapacity_vr(NUM_col(NY,NX),NY,NX).GT.ZEROS(NY,NX) .and. abs(QSnoHeatXfer2Soil_col(NY,NX))>ZEROS(NY,NX))THEN
!      TKS_vr(NUM_col(NY,NX),NY,NX) = (ENGY+QSnoHeatXfer2Soil_col(NY,NX))/VHeatCapacity_vr(NUM_col(NY,NX),NY,NX)
!    ELSEIF(VHeatCapacity_vr(NUM_col(NY,NX),NY,NX).LE.ZEROS(NY,NX))then
!      TKS_vr(NUM_col(NY,NX),NY,NX)=TairK_col(NY,NX)
!    ENDIF

!  endif

  call SnowpackDisapper(I,J,NY,NX,test_exist,QWatinfl2Mic_loc,QHeatInfl2Soil_loc)

  if(present(QWatinfl2Mic))QWatinfl2Mic=QWatinfl2Mic_loc
  if(present(QHeatInfl2Soil))QHeatInfl2Soil=QHeatInfl2Soil_loc
  TCSnow_snvr(1,NY,NX)=units%Kelvin2Celcius(TKSnow_snvr(1,NY,NX))

  call PrintInfo('end '//subname)
  end subroutine SnowMassUpdate

!------------------------------------------------------------------------------------------

  subroutine SnowpackDisapper(I,J,NY,NX,test_exist,QWatinfl2Mic,QHeatInfl2Soil)
  ! Transfer a removed top snow store's constituents and sensible energy to
  ! the receiving surface, retaining frozen water as ice in the native model.
  implicit none
  integer, intent(in) :: NY,NX,I,J
  logical, intent(in) :: test_exist
  real(r8),intent(out):: QWatinfl2Mic,QHeatInfl2Soil
  character(len=*), parameter :: subname='SnowpackDisapper'
  integer :: L
  real(r8) :: tksx
  real(r8) :: FLWS,FLWW,FLWI
  real(r8) :: HeatFlo2Surface,ENGYS,ENGY1,ENGY
!     begin_execution
!

  call PrintInfo('beg '//subname)
  QWatinfl2Mic=0._r8
  QHeatInfl2Soil=0._r8
  IF(VLHeatCapSnow_snvr(1,NY,NX).GT.0.0_r8 .AND. VLHeatCapSnow_snvr(1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX) &
    .AND. TairK_col(NY,NX).GT.TFICE)THEN

    !air temperature above freezing, surface snow layer heat insignificant, so it is merged
    !to the surface layer, and all varaibles are reset
    ENGYS           = TKSnow_snvr(1,NY,NX)*VLHeatCapSnow_snvr(1,NY,NX)
    ENGY1           = TKS_vr(NUM_col(NY,NX),NY,NX)*VHeatCapacity_vr(NUM_col(NY,NX),NY,NX)
    FLWS            = VLDrySnoWE_snvr(1,NY,NX)
    FLWW            = VLWatSnow_snvr(1,NY,NX)
    FLWI            = VLIceSnow_snvr(1,NY,NX)
    HeatFlo2Surface = (cpw*FLWW+cps*FLWS+cpi*FLWI)*TKSnow_snvr(1,NY,NX)

    !reset snow layer variables
    VLDrySnoWE_snvr(1,NY,NX)    = 0.0_r8
    VLWatSnow_snvr(1,NY,NX)     = 0.0_r8
    VLIceSnow_snvr(1,NY,NX)     = 0.0_r8
    VLHeatCapSnow_snvr(1,NY,NX) = 0.0_r8

    ! Remove only the transferred layer; occupied lower layers retain energy.
    VLSnoDWIprev_snvr(1,NY,NX)=0._r8
    SnowThickL_snvr(1,NY,NX)=0._r8
    SnoDens_snvr(1,NY,NX)=NewSnowDens_col(NY,NX)
    VcumDrySnoWE_col(NY,NX)=SUM(VLDrySnoWE_snvr(:,NY,NX))
    VcumWatSnow_col(NY,NX)=SUM(VLWatSnow_snvr(:,NY,NX))
    VcumIceSnow_col(NY,NX)=SUM(VLIceSnow_snvr(:,NY,NX))
    VcumSnoDWI_col(NY,NX)=SUM(VLSnoDWIprev_snvr(:,NY,NX))
    SnowDepth_col(NY,NX)=SUM(SnowThickL_snvr(:,NY,NX))
    VcumSnowWE_col(NY,NX)=VcumDrySnoWE_col(NY,NX)+VcumWatSnow_col(NY,NX)+DENSICE*VcumIceSnow_col(NY,NX)

    IF(SoilBulkDensity_vr(NUM_col(NY,NX),NY,NX).GT.ZERO .and. SoilOrgM_vr(ielmc,0,NY,NX)>1.e-2_r8)THEN
      ENGY                      = TKS_vr(0,NY,NX)*VHeatCapacity_vr(0,NY,NX)
      VLWatMicP_vr(0,NY,NX)     = VLWatMicP_vr(0,NY,NX)+FLWW
      VLiceMicP_vr(0,NY,NX)     = VLiceMicP_vr(0,NY,NX)+FLWI+FLWS/DENSICE
      VHeatCapacity_vr(0,NY,NX) = cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,0,NY,NX)) &
        +cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)
      IF(abs(HeatFlo2Surface)>ZEROS(NY,NX))THEN
        TKS_vr(0,NY,NX)           = (ENGY+HeatFlo2Surface)/VHeatCapacity_vr(0,NY,NX)
      ENDIF

      WatFLo2LitR_col(NY,NX)    = WatFLo2LitR_col(NY,NX) + FLWW+FLWI*DENSICE+FLWS
    else
      if(test_exist)then
        !This is for coupling with external models other than EcoSIM itself
        !The disappeared snow becomes water/heat flux to soil

        QHeatInfl2Soil = HeatFlo2Surface
        !dump all ice to micropore, assuming no macropore of the host model
        QWatinfl2Mic   = FLWW+FLWI*DENSICE+FLWS
      else
        write(*,*)I*1000+J,'SnowpackDisapper2'
        if(.not.fixWaterLevel)then
          QWatinfl2Mic   = 0._r8
          QHeatInfl2Soil = 0._r8
          !update top soil layer variables
          VLWatMicP_vr(NUM_col(NY,NX),NY,NX) = VLWatMicP_vr(NUM_col(NY,NX),NY,NX)+FLWW
          VLiceMicP_vr(NUM_col(NY,NX),NY,NX) = VLiceMicP_vr(NUM_col(NY,NX),NY,NX)+FLWI+FLWS/DENSICE

          ENGY  = VHeatCapacity_vr(NUM_col(NY,NX),NY,NX)*TKS_vr(NUM_col(NY,NX),NY,NX)

          VHeatCapacity_vr(NUM_col(NY,NX),NY,NX) = VHeatCapSolidSoil_vr(NUM_col(NY,NX),NY,NX) &
            +cpw*(VLWatMicP_vr(NUM_col(NY,NX),NY,NX)+VLWatMacP_vr(NUM_col(NY,NX),NY,NX)) &
            +cpi*(VLiceMicP_vr(NUM_col(NY,NX),NY,NX)+VLiceMacP_vr(NUM_col(NY,NX),NY,NX))

          IF(VHeatCapacity_vr(NUM_col(NY,NX),NY,NX).GT.ZEROS(NY,NX) .and. abs(HeatFlo2Surface)>ZEROS(NY,NX))THEN
            TKS_vr(NUM_col(NY,NX),NY,NX) = (ENGY+HeatFlo2Surface)/VHeatCapacity_vr(NUM_col(NY,NX),NY,NX)
          ELSEIF(VHeatCapacity_vr(NUM_col(NY,NX),NY,NX).LE.ZEROS(NY,NX))THEN
            TKS_vr(NUM_col(NY,NX),NY,NX)=TairK_col(NY,NX)
          ENDIF
        endif
        Qinflx2Soil_col(NY,NX)  = Qinflx2Soil_col(NY,NX)+FLWW+FLWI*DENSICE+FLWS
        QSnowH2Oloss_col(NY,NX) = QSnowH2Oloss_col(NY,NX)+FLWW+FLWI*DENSICE+FLWS
        QSnowHeatLoss_col(NY,NX)= QSnowHeatLoss_col(NY,NX)+HeatFlo2Surface
      endif
    endif
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine SnowpackDisapper
!------------------------------------------------------------------------------------------
  subroutine DealHighTempSnow(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8) :: SnoIce,dHeat,dIce,ENGY,dHPhaseChange
  real(r8) :: dVice,fi,fs,cpold,TKX

  SnoIce=VLDrySnoWE_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
  if(SnoIce<ZEROS(NY,NX))return
!    write(*,*)'not sufficient ice to counter high temperature drysnow ice',I+J/24.,VLDrySnoWE_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX)
!  endif
  TKX=TKSnow_snvr(L,NY,NX)
  !the starting enthalpy
  ENGY  = VLHeatCapSnow_snvr(L,NY,NX)*TKSnow_snvr(L,NY,NX)
  dHeat = VLHeatCapSnow_snvr(L,NY,NX)*(TKSnow_snvr(L,NY,NX)-TFICE)
  fs    = AZMAX1(VLDrySnoWE_snvr(L,NY,NX))/SnoIce
  fi    = (1._r8-fs)/DENSICE

  dIce=dHeat/(LtHeatIceMelt+(cpw-cpi*fi-cps*fs)*TFICE)

  if(dIce<SnoIce)then
  !some ice melt
  !  print*,VLDrySnoWE_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),fs,fi
    dHPhaseChange               = -dIce*LtHeatIceMelt  !thaw, cooling
    dVice                       = dIce/SnoIce
    VLDrySnoWE_snvr(L,NY,NX)    = VLDrySnoWE_snvr(L,NY,NX)-fs*dIce
    VLIceSnow_snvr(L,NY,NX)     = VLIceSnow_snvr(L,NY,NX)-fi*dIce
    VLWatSnow_snvr(L,NY,NX)     = VLWatSnow_snvr(L,NY,NX)+dIce
    cpold                       = VLHeatCapSnow_snvr(L,NY,NX)
    VLHeatCapSnow_snvr(L,NY,NX) = cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)
    TKSnow_snvr(L,NY,NX)=(ENGY+dHPhaseChange)/VLHeatCapSnow_snvr(L,NY,NX)
  else
    !all ice melt
    dHPhaseChange               = -SnoIce*LtHeatIceMelt  !thaw, cooling
    VLDrySnoWE_snvr(L,NY,NX)    = 0._r8
    VLIceSnow_snvr(L,NY,NX)     = 0._r8
    VLWatSnow_snvr(L,NY,NX)     = VLWatSnow_snvr(L,NY,NX)+SnoIce
    VLHeatCapSnow_snvr(L,NY,NX) = cpw*VLWatSnow_snvr(L,NY,NX)
    TKSnow_snvr(L,NY,NX)        = (ENGY+dHPhaseChange)/VLHeatCapSnow_snvr(L,NY,NX)
  endif
  if(TKSnow_snvr(L,NY,NX)<200._r8)then
    write(*,*)'High temp',L,TKSnow_snvr(L,NY,NX),TKX
    call endrun(trim(mod_filename)//' at line',__LINE__)
  endif
  end subroutine DealHighTempSnow
!------------------------------------------------------------------------------------------

  subroutine DealNegativeSnowMass(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8) :: dHPhaseChange,VLDrySnoWEtmp
  real(r8) :: SnoWEtot,ENGYW,dVLDrySnoWEtmp
  real(r8) :: tKNew,vlheatnew,dIce,dHeat,fs,fi
  real(r8) :: dwat,cphwat,vcphsnw


  cphwat  = cpw*VLWatSnow_snvr(L,NY,NX)
  vcphsnw = cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)
  if(abs(cphwat/vcphsnw)<1.e-3_r8)return

  !total snow mass
  SnoWEtot=VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
  if(TKSnow_snvr(L,NY,NX)<TFICE)then
!     write(*,*)'negative water',L,TKSnow_snvr(L,NY,NX)
!    call endrun('Temeprature too low to resolve negative water '//trim(mod_filename)//' at line',__LINE__)
  endif
  if(SnoWEtot<0._r8)then
    call endrun('Negative snow mass '//trim(mod_filename)//' at line',__LINE__)
  endif
  !the starting enthalpy
  ENGYW=VLHeatCapSnow_snvr(L,NY,NX)*TKSnow_snvr(L,NY,NX)

  !thaw all ice + snow, absorb heat/cooling (<0)
!  if(I>=138.and.I<=139)write(149,*)'neg',I+J/24.,VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX),VLDrySnoWE_snvr(L,NY,NX),&
!    TKSnow_snvr(L,NY,NX),(CumHeat2SnowL_snvr(L,NY,NX)+XPhaseChangeHeatL_snvr(L,NY,NX)+ENGYW)/(cpi*TFICE)
  dHPhaseChange=-LtHeatIceMelt*(VLDrySnoWE_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE)

  vlheatnew=cpw*SnoWEtot
  !compute potential temperature
  tkNew=(ENGYW+CumHeat2SnowL_snvr(L,NY,NX)+XPhaseChangeHeatL_snvr(L,NY,NX)+dHPhaseChange)/vlheatnew

  if(tkNew>TFICE)then
    !all ice & snow are melt
    XPhaseChangeHeatL_snvr(L,NY,NX) = XPhaseChangeHeatL_snvr(L,NY,NX)+dHPhaseChange
    VLWatSnow_snvr(L,NY,NX)         = SnoWEtot
    VLDrySnoWE_snvr(L,NY,NX)        = 0._r8
    VLIceSnow_snvr(L,NY,NX)         = 0._r8
  else
    !some snow and ice remains, refreeze/release heat
    dHeat=-(tkNew-TFICE)*vlheatnew    !>0, enthalpy excess for freeze thaw
    VLDrySnoWEtmp=AZMAX1(VLDrySnoWE_snvr(L,NY,NX))+AZMAX1(VLIceSnow_snvr(L,NY,NX))*DENSICE
    fs=AZMAX1(VLDrySnoWE_snvr(L,NY,NX))/VLDrySnoWEtmp
    fi=(1._r8-fs)/DENSICE
    dIce=dHeat/(LtHeatIceMelt+(cpw-cpi*fi-cps*fs)*TFICE)          !>0

    if(dIce <= SnoWEtot)then
      if(VLDrySnoWEtmp>0._r8)then
        dVLDrySnoWEtmp           = dIce/VLDrySnoWEtmp
        VLDrySnoWE_snvr(L,NY,NX) = dVLDrySnoWEtmp*AZMAX1(VLDrySnoWE_snvr(L,NY,NX))
        VLIceSnow_snvr(L,NY,NX)  = dVLDrySnoWEtmp*AZMAX1(VLIceSnow_snvr(L,NY,NX))
        VLWatSnow_snvr(L,NY,NX)  = AZMAX1(SnoWEtot-VLDrySnoWE_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX)*DENSICE)
        !freeze releases heat
        dHPhaseChange=(VLDrySnoWE_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX)*DENSICE-VLDrySnoWEtmp)*LtHeatIceMelt
      else
        !meaning water is positive, freezing occured, heat release
        VLIceSnow_snvr(L,NY,NX)  = dIce/DENSICE
        VLDrySnoWE_snvr(L,NY,NX) = 0._r8
        VLWatSnow_snvr(L,NY,NX)  = AZMAX1(SnoWEtot-VLDrySnoWE_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX)*DENSICE)
        dHPhaseChange            = dIce*LtHeatIceMelt
      endif

      XPhaseChangeHeatL_snvr(L,NY,NX)=XPhaseChangeHeatL_snvr(L,NY,NX)+dHPhaseChange
    else
      !negative water is made by ice, no phase change, but there is conversion
      dVLDrySnoWEtmp           = SnoWEtot/VLDrySnoWEtmp
      VLDrySnoWE_snvr(L,NY,NX) = dVLDrySnoWEtmp*AZMAX1(VLDrySnoWE_snvr(L,NY,NX))
      VLIceSnow_snvr(L,NY,NX)  = dVLDrySnoWEtmp*AZMAX1(VLIceSnow_snvr(L,NY,NX))
      VLWatSnow_snvr(L,NY,NX)  = 0._r8
    endif
  endif

  end subroutine DealNegativeSnowMass

!------------------------------------------------------------------------------------------

  subroutine UpdateSnowLayerL(I,J,L,NY,NX,VOLSWI)
  ! Update snow geometry and diagnostics while retaining each positive store's energy.
  implicit none
  integer, intent(in) :: L,NY,NX,I,J
  real(r8), intent(inout) :: VOLSWI
  character(len=*), parameter :: subname='UpdateSnowLayerL'
  real(r8) :: TKWX,VHCPWZ(JZ,JY,JX)
  real(r8) :: ENGYW
  real(r8) :: VOLSF,TCASF
  real(r8) :: DENSF,DDENS2,DDENS1
  real(r8) :: CVISC                        !snow compaction factor
  real(r8) :: DENSX
  real(r8) :: frcnew,SnowIceMass,VLwatNet
  real(r8) :: vwat,vice,vdry
  logical :: watchk,icechk
  real(r8), parameter :: tinyw=1.e-13_r8
! begin_execution

  call PrintInfo('beg '//subname)
  ! Phase change was accepted by the substep solver even if the layer has
  ! since disappeared or transferred its remaining constituents to the ground.
  THeatSnowThaw_col(NY,NX)=THeatSnowThaw_col(NY,NX)+XPhaseChangeHeatL_snvr(L,NY,NX)
  !the line below is a hack, and likely a better snow layering scheme is needed.
  if(L.eq.1.and.isclose(TCSnow_snvr(1,NY,NX),spval))then
    TCSnow_snvr(1,NY,NX)=units%Kelvin2Celcius(TairK_col(NY,NX))
  endif
!
! ADD CHANGES IN SNOW, WATER AND ICE
!
! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
! CumSno2SnowLay,CumWat2SnowLay,CumIce2SnowLay=net fluxes of snow,water,ice in snowpack
! XSnowThawMassL,XIceThawMassL=freeze-thaw flux from watsub.f
!
!
  vwat = VLWatSnow_snvr(L,NY,NX);vice = VLIceSnow_snvr(L,NY,NX);vdry = VLDrySnoWE_snvr(L,NY,NX)

!
! ACCUMULATE SNOW MASS FOR CALCULATING COMPRESSION
!
! VOLWSI=accumulated water equivalent volume in snowpack
! XFLWS=snow transfer from watsub.f
! VOLSF=snowfall volume
! DENSS=snow density in layer
!

  IF(L.EQ.1)THEN
    VOLSWI=VOLSWI+(VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE)

    if(VOLSWI<0._r8)then
      write(*,*)'VOLSWI=',VOLSWI,VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX)*DENSICE
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
!
!   RESET SNOW SURFACE DENSITY FOR SNOWFALL
!
    IF(SnoXfer2SnoLay_snvr(L,NY,NX).GT.0.0_r8)THEN
      !DENSF: fresh snow density
      DENSX = SnoDens_snvr(L,NY,NX)
      TCASF = AMAX1(-15.0_r8,AMIN1(2.0_r8,TCA_col(NY,NX)))
      DENSF = 0.05_r8+1.7E-03_r8*(TCASF+15.0_r8)**1.5_r8
      VOLSF = AMIN1(SnoXfer2SnoLay_snvr(L,NY,NX),VLDrySnoWE_snvr(L,NY,NX))/DENSF + &
        AZMAX1(VLDrySnoWE_snvr(L,NY,NX)-SnoXfer2SnoLay_snvr(L,NY,NX))/SnoDens_snvr(L,NY,NX)
      if(VOLSF>0._r8)SnoDens_snvr(L,NY,NX)=VLDrySnoWE_snvr(L,NY,NX)/VOLSF
      !write(*,*)'xVOLSSL=',VLDrySnoWE_snvr(L,NY,NX),SnoXfer2SnoLay_snvr(L,NY,NX),SnoDens_snvr(L,NY,NX),VOLSF
    ENDIF
  ELSE
    VOLSWI=VOLSWI+0.5_r8*(VLDrySnoWE_snvr(L-1,NY,NX)+VLWatSnow_snvr(L-1,NY,NX) &
      +VLIceSnow_snvr(L-1,NY,NX)*DENSICE+VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX) &
      +VLIceSnow_snvr(L,NY,NX)*DENSICE)

    if(VOLSWI<0._r8)then
      write(*,*)'iVOLSWI=',VOLSWI,VLDrySnoWE_snvr(L-1,NY,NX)+VLWatSnow_snvr(L-1,NY,NX) &
        +VLIceSnow_snvr(L-1,NY,NX)*DENSICE,VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX) &
        +VLIceSnow_snvr(L,NY,NX)*DENSICE
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDIF
  !
  ! SNOWPACK COMPRESSION
  !
  ! DDENS1 = Temperature effect on snow density
  ! DDENS2 = Compression effect on snow density
  ! DENSS  = Snow density in layer
  ! VOLSSL,VOLWSL,VOLISL=snow water equivalent,water,ice volume in snowpack layer
  ! VLSnoDWIprev_snvr=snowpack layer volume
  ! DLYRS=snowpack layer depth
  ! cumSnowDepz_col=cumulative depth to bottom of snowpack layer
  ! VHCPW=snowpack layer heat capacity
  ! TKW,TCSnow=snowpack layer temperature K,oC
  ! CumHeat2SnowLay=convective heat fluxes of snow,water,ice in snowpack
  ! XPhaseChangeHeatL=latent heat flux from freeze-thaw from watsub.f
  ! HEATIN_lnd=cumulative net surface heat transfer
  ! VOLSS,VOLWS,VOLIS=total snow water equivalent, water, ice content of snowpack
  ! VOLS,SnowDepth=total snowpack volume, depth
  !
  if(TKSnow_snvr(L,NY,NX)<0.)call endrun('too low snow temp '//trim(mod_filename)//' at line',__LINE__)

  if(SnoDens_snvr(L,NY,NX)>0._r8)then
    IF(SnoDens_snvr(L,NY,NX).LT.0.25_r8)THEN
      DDENS1=SnoDens_snvr(L,NY,NX)*1.0E-05_r8*EXP(0.04_r8*TCSnow_snvr(L,NY,NX))
    ELSE
      DDENS1=0.0_r8
    ENDIF

    CVISC  = 0.25_r8*EXP(-0.08_r8*TCSnow_snvr(L,NY,NX)+23.0_r8*SnoDens_snvr(L,NY,NX))
    DDENS2 = SnoDens_snvr(L,NY,NX)*VOLSWI/(AREA_3D(3,NU_col(NY,NX),NY,NX)*CVISC)

    SnoDens_snvr(L,NY,NX)=SnoDens_snvr(L,NY,NX)+DDENS1+DDENS2

    if(SnoDens_snvr(L,NY,NX)<0._r8)then
      write(*,*)'DDENS1=',SnoDens_snvr(L,NY,NX),DDENS1,DDENS2,L
      write(*,*)SnoXfer2SnoLay_snvr(L,NY,NX),VLDrySnoWE_snvr(L,NY,NX)
      call endrun("negative snow dens")
    endif
  endif

  !there is snow in layer L
  VLSnoDWIprev_snvr(L,NY,NX)=VLDrySnoWE_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)

  IF(VLSnoDWIprev_snvr(L,NY,NX) .GT. 0._r8)THEN
    SnowThickL_snvr(L,NY,NX)    = AZMAX1(VLSnoDWIprev_snvr(L,NY,NX))/AREA_3D(3,NU_col(NY,NX),NY,NX)
    cumSnowDepz_col(L,NY,NX)    = cumSnowDepz_col(L-1,NY,NX)+SnowThickL_snvr(L,NY,NX)
    VHCPWZ(L,NY,NX)             = VLHeatCapSnow_snvr(L,NY,NX)
    TKWX                        = TKSnow_snvr(L,NY,NX)
    ENGYW                       = VLHeatCapSnow_snvr(L,NY,NX)*TKWX
    VLHeatCapSnow_snvr(L,NY,NX) = cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)

    ! Retain the layer's stored sensible energy, including for thin snow.
    ! Phase-change energy was accepted by the substep solve; no temperature
    ! reset or compensating external heat input is introduced here.
    IF(VLHeatCapSnow_snvr(L,NY,NX).GT.0._r8)THEN
      TKSnow_snvr(L,NY,NX)=ENGYW/VLHeatCapSnow_snvr(L,NY,NX)
    ENDIF
    !there is no snow in layer L, shrink, and add those mass to layer below?
  ELSE

    VLDrySnoWE_snvr(L,NY,NX)   = 0.0_r8
    VLWatSnow_snvr(L,NY,NX)    = 0.0_r8
    VLIceSnow_snvr(L,NY,NX)    = 0.0_r8
    VLSnoDWIprev_snvr(L,NY,NX) = 0.0_r8
    SnowThickL_snvr(L,NY,NX)   = 0.0_r8
    cumSnowDepz_col(L,NY,NX)   = cumSnowDepz_col(L-1,NY,NX)

    !top snow layer
    IF(L.EQ.1)THEN
      TKSnow_snvr(L,NY,NX) = TairK_col(NY,NX)
    ELSE
      TKSnow_snvr(L,NY,NX) = TKSnow_snvr(L-1,NY,NX)
    ENDIF
  ENDIF
  TCSnow_snvr(L,NY,NX)=units%Kelvin2Celcius(TKSnow_snvr(L,NY,NX))
  call PrintInfo('end '//subname)
  end subroutine UpdateSnowLayerL

!------------------------------------------------------------------------------------------

  subroutine SnowpackLayering(I,J,NY,NX)
  !
  !Description:
  !Relayering (by division or combination) snow after mass update
  !
  implicit none
  integer, intent(in) :: I,J,NY,NX

  character(len=*), parameter :: subname='SnowpackLayering'
  real(r8) :: FX,FY
  integer :: L,L1,L0
  real(r8) :: ENGY0X,ENGY0,ENGY1X,ENGY1
  real(r8) :: DDLYXS,DDLYRS
  real(r8) :: DDLYXX,VOLSLX
  real(r8) :: DrySnowWE,DrySnowVol

  integer :: IFLGLS,idn,idg,idsalt
  integer, parameter :: inochange=0
  integer, parameter :: iexpand=1
  integer, parameter :: ishrink=2
  integer :: nsnol_beg_col

!     begin_execution
! from surface to bottom, and modify the bottom layer
! there is snow

  call PrintInfo('beg '//subname)
  nsnol_beg_col=nsnol_col(NY,NX)

!  write(111,*)'nsnol_col=',nsnol_col(NY,NX)
  IF(VLHeatCapSnow_snvr(1,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
    D325: DO L=1,JS-1

      IF(VLSnoDWIprev_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
        !compute excessive thickness for layer L: DDLYXS
        DDLYXS=(VLSnoDWIMax_snvr(L,NY,NX)-VLDrySnoWE_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX) &
          -VLWatSnow_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX))/AREA_3D(3,L,NY,NX)

        !snow is expanding from layer L into L+1, or layer L is shrinking, but L+1 has snow
        IF(DDLYXS.LT.-ZERO .OR. SnowThickL_snvr(L+1,NY,NX).GT.ZERO)THEN
          !current volume is greater than allowed, or next layer exists
          !case 1: DDLYXS< 0, layer L extends into layer L+1, DDLYRS<0: amount of expand layer L into layer L+1
          !case 2: 0<DDLYXS<SnowThickL_snvr(L+1,NY,NX), layer L still has space, can take DDLYRS of layer L+1
          !case 3: DDLYXS>SnowThickL_snvr(L+1,NY,NX), layer L still has space, and can even hold DDLYRS of layer L+1
          !
          DDLYRS = AMIN1(DDLYXS,SnowThickL_snvr(L+1,NY,NX))
          IFLGLS = iexpand         !expand
        ELSE
          !layer L has space and layer L+1 disappears, so layer L+1 combines into layer L
          !volume less than allowed, and no next layer
          !DDLYXS: is the depth change of layer L
          DDLYXS=(VLSnoDWIprev_snvr(L,NY,NX)-VLDrySnoWE_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX) &
            -VLWatSnow_snvr(L,NY,NX)-VLIceSnow_snvr(L,NY,NX))/AREA_3D(3,L,NY,NX)
          DDLYRS = DDLYXS
          IFLGLS = ishrink         !shrink
        ENDIF
      !layer L previously has nothing
      ELSE
        !current layer is empty, do nothing
        DDLYRS = 0.0_r8      !no change
        IFLGLS = inochange
      ENDIF
      !
      !     RESET SNOW LAYER DEPTHS
      !
      cumSnowDepz_col(L,NY,NX) = cumSnowDepz_col(L,NY,NX)+DDLYRS
      SnowThickL_snvr(L,NY,NX) = cumSnowDepz_col(L,NY,NX)-cumSnowDepz_col(L-1,NY,NX)
!
      !     TRANSFER STATE VARIABLES BETWEEN LAYERS
      !
      IF(ABS(DDLYRS).GT.ZERO)THEN
        IF(DDLYRS.GT.0.0_r8)THEN
          !incoporating L+1 into L
          L1 = L     !dest
          L0 = L+1   !source
          IF(DDLYRS.LT.DDLYXS)THEN
            !full L0 into L1
            FX=1.0_r8
          ELSE
            !partial L0 into L1
            FX=AMIN1(1.0_r8,DDLYRS*AREA_3D(3,L0,NY,NX)/VLSnoDWIprev_snvr(L0,NY,NX))
          ENDIF
        ELSE
          !expanding L into L+1
          L1 = L+1  !dest
          L0 = L    !source
          IF(VLSnoDWIprev_snvr(L0,NY,NX).LT.VLSnoDWIMax_snvr(L0,NY,NX))THEN
            FX=0.0_r8
          ELSE
            !FX fraction to be donated from L0 to L1
            FX=AMIN1(1.0_r8,-DDLYRS*AREA_3D(3,L0,NY,NX)/VLSnoDWIprev_snvr(L0,NY,NX))
          ENDIF
        ENDIF
!   donor L0, target L1
        IF(FX.GT.0.0_r8)THEN
          FY=1.0_r8-FX
!
!     TARGET SNOW LAYER
!         volume/mass
          ! Relayering preserves dry-snow bulk volume; compression is separate.
          DrySnowWE = VLDrySnoWE_snvr(L1,NY,NX)+FX*VLDrySnoWE_snvr(L0,NY,NX)
          DrySnowVol = VLDrySnoWE_snvr(L1,NY,NX)/SnoDens_snvr(L1,NY,NX) &
            +FX*VLDrySnoWE_snvr(L0,NY,NX)/SnoDens_snvr(L0,NY,NX)
          if(DrySnowVol>0._r8)then
            SnoDens_snvr(L1,NY,NX)=DrySnowWE/DrySnowVol
          endif
          VLDrySnoWE_snvr(L1,NY,NX)   = VLDrySnoWE_snvr(L1,NY,NX)+FX*VLDrySnoWE_snvr(L0,NY,NX)
          VLWatSnow_snvr(L1,NY,NX)    = VLWatSnow_snvr(L1,NY,NX)+FX*VLWatSnow_snvr(L0,NY,NX)
          VLIceSnow_snvr(L1,NY,NX)    = VLIceSnow_snvr(L1,NY,NX)+FX*VLIceSnow_snvr(L0,NY,NX)
          VLSnoDWIprev_snvr(L1,NY,NX) = VLDrySnoWE_snvr(L1,NY,NX)/SnoDens_snvr(L1,NY,NX)+VLWatSnow_snvr(L1,NY,NX)+VLIceSnow_snvr(L1,NY,NX)

!         energy
          ENGY1X                       = VLHeatCapSnow_snvr(L1,NY,NX)*TKSnow_snvr(L1,NY,NX)
          ENGY0X                       = VLHeatCapSnow_snvr(L0,NY,NX)*TKSnow_snvr(L0,NY,NX)
          ENGY1                        = ENGY1X+FX*ENGY0X
          VLHeatCapSnow_snvr(L1,NY,NX) = cps*VLDrySnoWE_snvr(L1,NY,NX)+cpw*VLWatSnow_snvr(L1,NY,NX)+cpi*VLIceSnow_snvr(L1,NY,NX)

          IF(VLHeatCapSnow_snvr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
            TKSnow_snvr(L1,NY,NX)=ENGY1/VLHeatCapSnow_snvr(L1,NY,NX)
          ELSE
            TKSnow_snvr(L1,NY,NX)=TKSnow_snvr(L0,NY,NX)
          ENDIF
          if(TKSnow_snvr(L1,NY,NX)<200._r8)then
            write(*,*)'wwwerid temp',L1,TKSnow_snvr(L1,NY,NX),TKSnow_snvr(L0,NY,NX)
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          TCSnow_snvr(L1,NY,NX)=units%Kelvin2Celcius(TKSnow_snvr(L1,NY,NX))
          !------------------------------------------------------------
!          chemicals
          !gas
          DO idg=idg_beg,idg_NH3
            trcg_solsml_snvr(idg,L1,NY,NX)=AZERO(trcg_solsml_snvr(idg,L1,NY,NX)+FX*trcg_solsml_snvr(idg,L0,NY,NX))
          ENDDO

          !nutrients
          DO idn=ids_nut_beg,ids_nuts_end
            trcn_solsml_snvr(idn,L1,NY,NX)=AZERO(trcn_solsml_snvr(idn,L1,NY,NX)+FX*trcn_solsml_snvr(idn,L0,NY,NX))
          ENDDO
          !salt
          IF(salt_model)THEN
            DO idsalt=idsalt_beg,idsalt_end
              trcSalt_ml_snvr(idsalt,L1,NY,NX)=AZERO(trcSalt_ml_snvr(idsalt,L1,NY,NX)+FX*trcSalt_ml_snvr(idsalt,L0,NY,NX))
            ENDDO
          ENDIF
          !------------------------------------------------------------
!
!     SOURCE SNOW LAYER
!         volume/mass
          VLDrySnoWE_snvr(L0,NY,NX)   = FY*VLDrySnoWE_snvr(L0,NY,NX)
          VLWatSnow_snvr(L0,NY,NX)    = FY*VLWatSnow_snvr(L0,NY,NX)
          VLIceSnow_snvr(L0,NY,NX)    = FY*VLIceSnow_snvr(L0,NY,NX)
          VLSnoDWIprev_snvr(L0,NY,NX) = VLDrySnoWE_snvr(L0,NY,NX)/SnoDens_snvr(L0,NY,NX)+VLWatSnow_snvr(L0,NY,NX)+VLIceSnow_snvr(L0,NY,NX)
!     energy
          ENGY0=FY*ENGY0X
          VLHeatCapSnow_snvr(L0,NY,NX)=cps*VLDrySnoWE_snvr(L0,NY,NX)+cpw*VLWatSnow_snvr(L0,NY,NX)+cpi*VLIceSnow_snvr(L0,NY,NX)
          IF(VLHeatCapSnow_snvr(L0,NY,NX).GT.ZEROS(NY,NX))THEN
            TKSnow_snvr(L0,NY,NX)=ENGY0/VLHeatCapSnow_snvr(L0,NY,NX)
          ELSE
            TKSnow_snvr(L0,NY,NX)=TKSnow_snvr(L1,NY,NX)
          ENDIF
          if(TKSnow_snvr(L0,NY,NX)<200._r8)then
            write(*,*)'eeeewird',L0,TKSnow_snvr(L0,NY,NX),TKSnow_snvr(L1,NY,NX)
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          TCSnow_snvr(L0,NY,NX)=units%Kelvin2Celcius(TKSnow_snvr(L0,NY,NX))

          !------------------------------------------------------------
          !   chemicals
          DO idg=idg_beg,idg_NH3
            trcg_solsml_snvr(idg,L0,NY,NX)=FY*trcg_solsml_snvr(idg,L0,NY,NX)
          ENDDO
          DO idn=ids_nut_beg,ids_nuts_end
            trcn_solsml_snvr(idn,L0,NY,NX)=FY*trcn_solsml_snvr(idn,L0,NY,NX)
          ENDDO

          IF(salt_model)THEN
            DO idsalt=idsalt_beg,idsalt_end
              trcSalt_ml_snvr(idsalt,L0,NY,NX)=FY*trcSalt_ml_snvr(idsalt,L0,NY,NX)
            ENDDO
          ENDIF
          !------------------------------------------------------------
        ENDIF
      ENDIF
    ENDDO D325
  ENDIF

  !update volumetric snow heat capacity
  nsnol_col(NY,NX)=0
  DO L=1,JS
    if(TKSnow_snvr(L,NY,NX)/=spval)then
      VLHeatCapSnow_snvr(L,NY,NX)=cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)
      if(VLHeatCapSnow_snvr(L,NY,NX)>ZEROS(NY,NX))then
        nsnol_col(NY,NX)=nsnol_col(NY,NX)+1
      endif
    else
      VLDrySnoWE_snvr(L,NY,NX)    = 0._r8
      VLWatSnow_snvr(L,NY,NX)     = 0._r8
      VLIceSnow_snvr(L,NY,NX)     = 0._r8
      VLHeatCapSnow_snvr(L,NY,NX) = 0._r8
      SnoDens_snvr(L,NY,NX)       = NewSnowDens_col(NY,NX)
    endif
    if(VLHeatCapSnow_snvr(L,NY,NX)<=0._r8)then
      if(L>1)then
        TKSnow_snvr(L,NY,NX)=spval
      else
        TKSnow_snvr(L,NY,NX)=TairK_col(NY,NX)
      endif
    endif

  ENDDO

  ! Keep occupied layers contiguous from the surface for tracer transport.
  ! A positive store remains occupied even below the heat-solver threshold.
  ! Residual chemicals in an empty destination remain with the incoming snow.
  DO L=2,JS
    if(VLHeatCapSnow_snvr(L,NY,NX)>0._r8)then
      L1=1
      do while(L1<L)
        if(VLHeatCapSnow_snvr(L1,NY,NX)<=0._r8)exit
        L1=L1+1
      enddo
      if(L1<L)then
        VLDrySnoWE_snvr(L1,NY,NX)    = VLDrySnoWE_snvr(L,NY,NX)
        VLWatSnow_snvr(L1,NY,NX)     = VLWatSnow_snvr(L,NY,NX)
        VLIceSnow_snvr(L1,NY,NX)     = VLIceSnow_snvr(L,NY,NX)
        VLHeatCapSnow_snvr(L1,NY,NX) = VLHeatCapSnow_snvr(L,NY,NX)
        TKSnow_snvr(L1,NY,NX)        = TKSnow_snvr(L,NY,NX)
        SnoDens_snvr(L1,NY,NX)       = SnoDens_snvr(L,NY,NX)

        VLDrySnoWE_snvr(L,NY,NX)    = 0._r8
        VLWatSnow_snvr(L,NY,NX)     = 0._r8
        VLIceSnow_snvr(L,NY,NX)     = 0._r8
        VLHeatCapSnow_snvr(L,NY,NX) = 0._r8
        TKSnow_snvr(L,NY,NX)        = spval
        SnoDens_snvr(L,NY,NX)       = NewSnowDens_col(NY,NX)

        DO idg=idg_beg,idg_NH3
          trcg_solsml_snvr(idg,L1,NY,NX) = trcg_solsml_snvr(idg,L1,NY,NX)+trcg_solsml_snvr(idg,L,NY,NX)
          trcg_solsml_snvr(idg,L,NY,NX)   = 0._r8
        ENDDO
        DO idn=ids_nut_beg,ids_nuts_end
          trcn_solsml_snvr(idn,L1,NY,NX) = trcn_solsml_snvr(idn,L1,NY,NX)+trcn_solsml_snvr(idn,L,NY,NX)
          trcn_solsml_snvr(idn,L,NY,NX)   = 0._r8
        ENDDO

        IF(salt_model)THEN
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_ml_snvr(idsalt,L1,NY,NX) = trcSalt_ml_snvr(idsalt,L1,NY,NX)+trcSalt_ml_snvr(idsalt,L,NY,NX)
            trcSalt_ml_snvr(idsalt,L,NY,NX)   = 0._r8
          ENDDO
        ENDIF
      endif
    endif
  ENDDO

  ! Rebuild geometry and diagnostics from the final layer order.
  nsnol_col(NY,NX)=0
  DO L=1,JS
    VLSnoDWIprev_snvr(L,NY,NX) = VLDrySnoWE_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX) &
      +VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)
    SnowThickL_snvr(L,NY,NX) = VLSnoDWIprev_snvr(L,NY,NX)/AREA_3D(3,L,NY,NX)
    cumSnowDepz_col(L,NY,NX) = cumSnowDepz_col(L-1,NY,NX)+SnowThickL_snvr(L,NY,NX)
    if(TKSnow_snvr(L,NY,NX)/=spval)then
      TCSnow_snvr(L,NY,NX)=units%Kelvin2Celcius(TKSnow_snvr(L,NY,NX))
    else
      TCSnow_snvr(L,NY,NX)=spval
    endif
    if(VLHeatCapSnow_snvr(L,NY,NX)>0._r8)then
      nsnol_col(NY,NX)=nsnol_col(NY,NX)+1
    endif
  ENDDO
  VcumSnoDWI_col(NY,NX) = sum(VLSnoDWIprev_snvr(1:JS,NY,NX))
  SnowDepth_col(NY,NX)  = sum(SnowThickL_snvr(1:JS,NY,NX))
!  if(I==19 .and. I>=10)write(115,*)I*1000+J,'nsnol_beg_col',nsnol_beg_col,nsnol_col(NY,NX),trcs_solml_vr(idg_O2,0,NY,NX),&
!    trcg_solsml_snvr(idg_O2,1,NY,NX)
  if(nsnol_col(NY,NX)==0)then
    ! Include residual tracers in empty slots beyond the former active count.
    ! Add all tracers to the litter layer.
    DO L=1,JS
      DO idg=idg_beg,idg_NH3
        if(trcg_solsml_snvr(idg,L,NY,NX)>0._r8)then
          trcs_solml_vr(idg,0,NY,NX)       = trcs_solml_vr(idg,0,NY,NX)+trcg_solsml_snvr(idg,L,NY,NX)
          trcg_snowMassloss_col(idg,NY,NX) = trcg_snowMassloss_col(idg,NY,NX)+trcg_solsml_snvr(idg,L,NY,NX)
          trcg_solsml_snvr(idg,L,NY,NX)    = 0._r8
        endif
      ENDDO

      do idn=ids_nut_beg,ids_nuts_end
        if(trcn_solsml_snvr(idn,L,NY,NX)>0._r8)then
          trcs_solml_vr(idn,0,NY,NX)       = trcs_solml_vr(idn,0,NY,NX)+trcn_solsml_snvr(idn,L,NY,NX)
          trcn_snowMassloss_col(idn,NY,NX) = trcn_snowMassloss_col(idn,NY,NX)+trcn_solsml_snvr(idn,L,NY,NX)
          trcn_solsml_snvr(idn,L,NY,NX)    = 0._r8
        endif
      ENDDO

      IF(salt_model)THEN
        DO idsalt=idsalt_beg,idsalt_end
          if(trcSalt_ml_snvr(idsalt,L,NY,NX)>0._r8) then
            trcSalt_solml_vr(idsalt,0,NY,NX)       = trcSalt_solml_vr(idsalt,0,NY,NX)+trcSalt_ml_snvr(idsalt,L,NY,NX)
            trcSalt_snowMassloss_col(idsalt,NY,NX) = trcSalt_snowMassloss_col(idsalt,NY,NX)+trcSalt_ml_snvr(idsalt,L,NY,NX)
            trcSalt_ml_snvr(idsalt,L,NY,NX)        = 0._r8
          endif
        ENDDO
      ENDIF
    ENDDO
  endif

!  if(I==19 .and. I>=10)write(115,*)I*1000+J,'nsnol_afe_col',nsnol_beg_col,nsnol_col(NY,NX),trcs_solml_vr(idg_O2,0,NY,NX),&
!    trcg_solsml_snvr(idg_O2,1,NY,NX)
!  write(111,*)'nsnol_col=',nsnol_col(NY,NX)
  call PrintInfo('end '//subname)
  end subroutine SnowpackLayering

!------------------------------------------------------------------------------------------
  subroutine ZeroSnowArrays(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  trcg_LossXSnowRedist_col(idg_beg:idg_NH3,NY,NX)          = 0.0_r8
  trcn_LossXSnowRedist_col(ids_nut_beg:ids_nuts_end,NY,NX) = 0.0_r8

  DO  L=1,JS
    trcg_AquaAdv_NetFlx_snvr(idg_beg:idg_NH3,L,NY,NX)          = 0.0_r8
    trcn_AquaAdv_NetFlx_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX) = 0.0_r8
  ENDDO

  IF(salt_model)THEN
!     INITIALIZE NET SOLUTE AND GAS FLUXES FROM SNOWPACK DRIFT
!
    trcSalt_LossXSnowRedist_col(idsalt_beg:idsalt_end,NY,NX)=0.0_r8
    DO  L=1,JS
      trcSalt_AquaAdv_NetFlx_snvr(idsalt_beg:idsalt_end,L,NY,NX)=0.0_r8
    ENDDO
  endif
  end subroutine ZeroSnowArrays

end module SnowBalanceMod
