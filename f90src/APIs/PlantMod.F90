module PlantMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8,yearIJ_type
  use grosubsMod,        only: GrowPlants
  use PlantPhenolMod,    only: PhenologyUpdate
  use UptakesMod,        only: RootUptakes
  use PlantDisturbMod,   only: PrepLandscapeGrazing
  use PlantMgmtDataType, only: NP_col, NP0_col
  use MiniMathMod,       only: fixEXConsumpFlux
  use EcoSIMCtrlMod,     only: lverb, ldo_sp_mode
  use PlantDebugMod,     only: PrintRootTracer
  use LitterFallMod,     only: ReSeedPlants
  use DebugToolMod,      only: PrintInfo
  use InitPlantMod,      only: SeedDeposition  
  use PlantAPI4Uptake
  use TracerIDMod
  use GridDataType
  use PlantDataRateType
  use SoilBGCDataType
  use EcoSimSumDataType
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantAPI
  use ClimForcDataType
  use PlantTraitDataType
  use PlantCanAPI
  use ExtractsMod
  use PlantBalMod
implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PlantModel
  public :: PlantCanopyRadsModel
  public :: UpdateDailyGPP

  contains

  subroutine PlantModel(yearIJ,NHW,NHE,NVN,NVS)
  !
  !run the plant biogeochemistry model
  !
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8) :: t1,tvegE(NumPlantChemElms)
  real(r8) :: SeedCDeposition_1d(JP1), PPmax_1d(JP1)

  integer :: NY,NX,NZ
  character(len=*), parameter :: subname='PlantModel'
333   FORMAT(A8)

  call PrintInfo('beg '//subname)
  if(.not.ldo_sp_mode)then
    call PrepLandscapeGrazing(yearIJ%I,yearIJ%J,NHW,NHE,NVN,NVS)
    plt_site%PlantElemntStoreLandscape(:)=PlantElemntStoreLandscape(:)
  endif

  DO NX=NHW,NHE
    DO NY=NVN,NVS

      if(ldo_sp_mode)then
        !do prescribed phenolgoy
        call PlantUptakeAPISend(yearIJ%I,yearIJ%J,NY,NX)        

        CALL ROOTUPTAKES(yearIJ)
        
        call extracts(yearIJ)
        
        call PlantUPtakeAPIRecv(yearIJ%I,yearIJ%J,NY,NX)
      else
      
        call  PlantAPISend(yearIJ,NY,NX)

        call EnterPlantBalance(yearIJ,NP0_col(NY,NX))

        if(yearIJ%J.eq.1) then          
          SeedCDeposition_1d = SeedCDeposition_pft(:,NY,NX)
          PPmax_1d = PPmax_pft(:,NY,NX)
          call SeedDeposition(yearIJ%I,SeedCDeposition_1d,PPmax_1d,NP_col(NY,NX))          
        endif

        !Phenological update, determine living/active branches      
        CALL PhenologyUpdate(yearIJ)

        !Predict uptake fluxes of nutrients and O2        
        CALL ROOTUPTAKES(yearIJ)

        !Do growth of active branches and roots
        CALL GROWPLANTS(yearIJ)

        !aggregate varaibles
        CALL EXTRACTs(yearIJ)
        
        Call ReSeedPlants(yearIJ,NP_col(NY,NX))
        
        call ExitPlantBalance(yearIJ,NP0_col(NY,NX))

        call PlantAPIRecv(yearIJ%I,yearIJ%J,NY,NX)
        
      endif

    ENDDO
  ENDDO
  
  PlantElemntStoreLandscape(:)=plt_site%PlantElemntStoreLandscape(:)
  call PrintInfo('end '//subname)

  end subroutine PlantModel
!------------------------------------------------------------------------------------------
  subroutine UpdateDailyGPP(J,NHW,NHE,NVN,NVS,PlantGPPAvailable)
  ! Called once after each hourly plant update, independently of warming dates.
  implicit none
  integer, intent(in) :: J,NHW,NHE,NVN,NVS
  logical, intent(in) :: PlantGPPAvailable !false when no current prognostic GPP was computed
  character(len=*), parameter :: subname='UpdateDailyGPP'
  integer :: NY,NX

  call PrintInfo('beg '//subname)
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      ! Keep the daily GPP resets here rather than moving them to DAY.
      ! DAY precedes checkpoint loading: intraday restarts must retain the
      ! restored partial sum/count, while hour 1 starts a fresh accumulation.
      if(J==1)then
        GPP_DaySum_col(NY,NX)=0._r8
        GPP_DayHours_col(NY,NX)=0
      endif
      if(PlantGPPAvailable)then
        GPP_DaySum_col(NY,NX)=GPP_DaySum_col(NY,NX) &
          +sum(max(0._r8,GrossCO2Fix_pft(1:NP_col(NY,NX),NY,NX)))
      endif
      GPP_DayHours_col(NY,NX)=GPP_DayHours_col(NY,NX)+1
      if(J==24)then
        ! Publish only after the day's last warming/plant calculations. Keep
        ! yesterday's total fixed during the following day's weather and soil calls.
        GPP_PrevDay_col(NY,NX)=0._r8
        if(GPP_DayHours_col(NY,NX)==24)GPP_PrevDay_col(NY,NX)=GPP_DaySum_col(NY,NX)
        ! Clear before checkpoint writing so midnight restarts load zero
        ! partial state, while retaining the completed total in GPP_PrevDay_col.
        GPP_DaySum_col(NY,NX)=0._r8
        GPP_DayHours_col(NY,NX)=0
      endif
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine UpdateDailyGPP
!------------------------------------------------------------------------------------------
!  function test_active_plant(NY,NX)result(activeroot)
  !
  !Description:
  !Check the existence of active plants
!  implicit none
!  integer, intent(in) :: NY,NX
!  logical :: activeroot
!  integer :: NZ

!  activeroot=.false.
!  DO NZ=1,NP_col(NY,NX)
!    activeroot=IsPlantActive_pft(NZ,NY,NX).EQ.iTrue .and. PlantPopuLive_pft(NZ,NY,NX)>ZEROS(NY,NX)
!    if(activeroot)return
!  ENDDO
!  end function test_active_plant

!------------------------------------------------------------------------------------------

  subroutine PlantCanopyRadsModel(I,J,NY,NX,DepthSurfWatIce)
  use SurfaceRadiationMod
  use CanopyDataType, only : CanopyLeafArea_lnode
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(in) :: DepthSurfWatIce

  call PlantAPICanMSend(NY,NX)

  call CanopyConditionModel(I,J,DepthSurfWatIce)

  call PlantAPICanMRecv(NY,NX)

  end subroutine PlantCanopyRadsModel

end module PlantMod
