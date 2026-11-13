module WatsubMod
!!
! Description:
! Do water and enerby balance calculation.
! The module diagnoses the mass and energy fluxes associated
! with soil/snow water (vapor, liquid and ice) and energy, and update
! them in redistmod.F90

  use SoilFaceMappingMod, only: SoilFaceSlot,SoilFaceIncomingRange, &
    RebuildSoilFaceMap,DestroySoilFaceMap
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use data_const_mod,    only: GravAcceleration=>DAT_CONST_G
  use abortutils,        only: endrun,                      print_info
  use ElmIDMod,          only: iWestEastDirection,          iNorthSouthDirection, iVerticalDirection
  use SurfPhysData,      only: InitSurfPhysData,            DestructSurfPhysData
  use EcosysWarmingMod,  only: check_warming_dates,         is_warming_layerL,    apply_soil_cable_warming
  use PlantDataRateType, only: TWaterPlantRoot2SoilPrev_vr, THeatPlantRoot2SoilPrev_vr
  use MiniFuncMod,       only: gOC_to_m3_OM  
  use DebugToolMod
  use EcoSIMCtrlMod  
  use minimathmod      
  use RootDataType
  use EcosimConst
  use MiniFuncMod
  use EcoSIMSolverPar
  use SOMDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilWaterDataType
  use SoilHeatDatatype
  use EcoSIMCtrlDataType
  use LandSurfDataType
  use ClimForcDataType
  use FertilizerDataType
  use SnowDataType
  use PlantTraitDataType
  use SurfLitterDataType
  use SurfSoilDataType
  use SurfSoilDataType
  use CanopyDataType
  use ChemTranspDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use IrrigationDataType
  use GridDataType
  use WatsubDataMod
  use PhysPars    
  use SnowPhysMod
  use HydroThermData
  use SurfPhysMod
  use SoilPhysParaMod
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__
  real(r8), parameter :: tiny_wat=1.e-12_r8   !threshold water flux

  integer :: curday,curhour
  logical :: do_warming

  public :: watsub,InitWatsub
  public :: DestroyWatsub
  public :: ExitSoilHydroUpdate
  contains
!------------------------------------------------------------------------------------------
  subroutine InitWatsub

  implicit none

  call InitSurfPhysData  

  call InitWatSubData

  end subroutine InitWatsub
!------------------------------------------------------------------------------------------
  subroutine DestroyWatsub
  implicit none
  call DestroySoilFaceMap
  call DestructSurfPhysData
  end subroutine DestroyWatsub
!------------------------------------------------------------------------------------------

  SUBROUTINE watsub(I,J,NHW,NHE,NVN,NVS)
! Description:
!     THIS SUBROUTINE CACULATES ENERGY BALANCES OF SNOW, RESIDUE
!     AND SOIL SURFACES, FREEZING, THAWING, AND HEAT AND WATER
!     TRANSFER THROUGH SOIL PROFILES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='watsub'
  integer :: K0,K1
  integer :: KL,L,L2,LL,M,MM,M1,M2,M3,M4,M5,M6,NX,NY
  integer :: N,N1,N2,N3,N4,N5,N6,NN,N4B,N5B,NUX
  real(r8):: ResistVapThroughLitr(JY,JX),HeatFluxAir2Soi(JY,JX)
  REAL(R8):: RainEkReducedKsat(JY,JX)

  REAL(R8) :: TopLayWatVol_col(JY,JX)
  real(r8) :: Qinfl2MicP_col(JY,JX)
  real(r8) :: HeatInfl2Soil(JY,JX)
  real(r8) :: Qinfl2MacP_col(JY,JX)
  real(r8) :: dHGrnddTKsoi1(JY,JX)    !step-integrated radiation+sensible heat derivative [MJ/K]
  real(r8) :: twatmass0(JY,JX),twatmasstM(JY,JX)
  integer  :: NUX0(JY,JX)
!  real(r8) :: dtime
! begin_execution

  call PrintInfo('beg '//subname)
  if(fixWaterLevel)HeatAdv_scal=0._r8

  do_warming=check_warming_dates(iYearCurrent,I,J)

  if(do_warming)then
    call apply_soil_cable_warming(I,J,NHW,NHE,NVN,NVS)  
  endif

  ! Shared surface staging retires snow before copying native working states.
  call StageSurfacePhysModel(I,J,NHW,NHE,NVN,NVS,ResistVapThroughLitr)
  !
  call BeginMassCheck(NHW,NHE,NVN,NVS,twatmass0)
  twatmasstM=twatmass0
  ! Hour1 reset infiltration. Staging's snow transfer is already in the copied
  ! soil, but remains part of the cumulative hourly infiltration ledger.
  twatmass0(NVN:NVS,NHW:NHE)=twatmass0(NVN:NVS,NHW:NHE)-Qinflx2Soil_col(NVN:NVS,NHW:NHE)
  !
  call InitSoilHydrauics(NHW,NHE,NVN,NVS)

!  dtime=0._r8
  ! iterate for NPH times
  D3320: DO M=1,NPH
    !begin mass check for ieration M

    ! dtime=dtime+dts_HeatWatTP
    call PrepHydroThermIterM(M,NHW,NHE,NVN,NVS,TopLayWatVol_col,NUX0)

    call RunSurfacePhysModelM(I,J,M,NHE,NHW,NVS,NVN,ResistVapThroughLitr,RainEkReducedKsat,&
      TopLayWatVol_col,HeatFluxAir2Soi,dHGrnddTKsoi1,Qinfl2MicP_col,HeatInfl2Soil,Qinfl2MacP_col)
    
    call CopySoilWatVolIterateM(I,J,M,NHW,NHE,NVN,NVS,TopLayWatVol_col)

    if(lverb)write(*,*)'subsurface 3D flow'    
    call Subsurface3DInternalFlowM(I,J,M,NHW,NHE,NVN,NVS,RainEkReducedKsat,HeatFluxAir2Soi)

    call XBoundaryFlowM(I,J,M,NHW,NHE,NVN,NVS,RainEkReducedKsat)

    call Summarize3DFlowM(I,J,M,NHW,NHE,NVN,NVS,Qinfl2MicP_col,HeatInfl2Soil,Qinfl2MacP_col)

    if(snowRedist_model)call AccumulateSnowRedisFluxM(I,J,M,NHW,NHE,NVN,NVS)

    call UpdateSoilMoistTempM(I,J,M,NHW,NHE,NVN,NVS,dHGrnddTKsoi1)

    call AggregateSurfRunoffFluxM(I,J,M,NHW,NHE,NVN,NVS)

    call UpdateSurfaceAtM(I,J,M,NHW,NHE,NVN,NVS)

    IF(M.NE.NPH)THEN
      call UpdateStateFluxAtM(I,J,M,NHW,NHE,NVN,NVS)
    ELSE
      call UpdateFluxAtExit(I,J,NHW,NHE,NVN,NVS)
    ENDIF

    call ExitMassBalanceM(I,J,M,NHW,NHE,NVN,NVS,twatmass0,NUX0,twatmasstM)

  ENDDO D3320

  call ExitSoilHydroUpdate(NHW,NHE,NVN,NVS)

  call PrintInfo('end '//subname)
  END subroutine watsub
!------------------------------------------------------------------------------------------

  subroutine ExitSoilHydroUpdate(NHW,NHE,NVN,NVS)
  ! Description:
  ! Diagnose the frozen boundary and commit accepted water, ice, temperature
  ! and heat capacity before later processes update the material inventories.
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NX,NY,L
  logical :: found_frozen
  real(r8) :: VLSoilPoreMicPX
  character(len=*), parameter :: subname='ExitSoilHydroUpdate'

  call PrintInfo('beg '//subname)
  !diagnose surface ice
  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      ! Without a frozen barrier, ebullition must scan upward from the current
      ! bottom layer. Reset here to avoid zero at startup or a stale thawed barrier.
      NLF_col(NY,NX)=NL_col(NY,NX)
      found_frozen=.false.
      DO L=NUM_col(NY,NX),NL_col(NY,NX)
        VLSoilPoreMicPX       = AREA_3D(3,L,NY,NX)*DLYR_3D(3,L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
        ! Macropore volume is already pore space; porosity applies only to the
        ! micropore-domain soil volume.
        VOLTX_vr(L,NY,NX) = POROS_vr(L,NY,NX)*VLSoilPoreMicPX+VLMacroPores_vr(L,NY,NX)
        IF(VOLTX_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
          ThetaICEZ_vr(L,NY,NX) = AZMAX1(AMIN1(1._r8, &
            (VLiceMicP1_vr(L,NY,NX)+AMIN1(VLMacroPores_vr(L,NY,NX), &
            VLiceMacP1_vr(L,NY,NX)))/VOLTX_vr(L,NY,NX)))
        ELSE
          ThetaICEZ_vr(L,NY,NX)=0._r8
        ENDIF
        if(ThetaICEZ_vr(L,NY,NX)>0.9_r8 .and. .not. found_frozen)then
          found_frozen=.true.
          NLF_col(NY,NX)=L
        endif
        !micropore water change due to freeze-thaw:=final-initial
        DVLiceMicP_vr(L,NY,NX)       = VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX)-VLiceMicP_vr(L,NY,NX)-VLiceMacP_vr(L,NY,NX)

        ! Commit accepted capacity with temperature before biology/disturbance
        ! can change the organic pools used by the later diagnostic refresh.
        VHeatCapacity_vr(L,NY,NX) = VHeatCapacity1_vr(L,NY,NX)
        TKS_vr(L,NY,NX)        = TKSoil1_vr(L,NY,NX)
        VLWatMicP_vr(L,NY,NX)  = VLWatMicP1_vr(L,NY,NX)
        VLWatMicPX_vr(L,NY,NX) = VLWatMicPX1_vr(L,NY,NX)
        VLiceMicP_vr(L,NY,NX)  = VLiceMicP1_vr(L,NY,NX)
        VLWatMacP_vr(L,NY,NX)  = VLWatMacP1_vr(L,NY,NX)
        VLiceMacP_vr(L,NY,NX)  = VLiceMacP1_vr(L,NY,NX)
      ENDDO

      DO L=NU_col(NY,NX),NUM_col(NY,NX)-1
        DVLiceMicP_vr(L,NY,NX) = VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)
      ENDDO
      !PAR radiation to Soil surface
      RadPAR2Soil_col(NY,NX) = RadPARGrnd_col(NY,NX)*FracSurfSnoFree_col(NY,NX)*FracSurfBareSoil_col(NY,NX)
      !PAR radiation to litter surface
      RadPAR2LitR_col(NY,NX) = RadPARGrnd_col(NY,NX)*FracSurfSnoFree_col(NY,NX)*FracSurfByLitR_col(NY,NX)

    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine ExitSoilHydroUpdate

!------------------------------------------------------------------------------------------  

  subroutine ExitMassBalanceM(I,J,M,NHW,NHE,NVN,NVS,twatmass0,NUX0,twatmasstM)

  implicit none
  integer, intent(in)    :: I,J,M
  integer, intent(in)    :: NHW,NHE,NVN,NVS
  real(r8),intent(in)    :: twatmass0(JY,JX)
  integer, intent(in)    :: NUX0(JY,JX)
  real(r8),intent(inout) :: twatmasstM(JY,JX)  

  real(r8)  :: twatmass1(JY,JX)
  character(len=*), parameter :: subname='ExitMassBalanceM'
  real(r8) :: dwat,dwat0,dwat1
  integer :: NY,NX,L

  call PrintInfo('beg '//subname)
  DO NX=NHW,NHE
    DO  NY=NVN,NVS

     twatmass1(NY,NX)=0._r8
     dwat0=0._r8
     DO L=NUX0(NY,NX),NUM_col(NY,NX)-1
       dwat0=dwat0 + VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE
     ENDDO

      D131: DO L=NUM_col(NY,NX),NL_col(NY,NX)
        twatmass1(NY,NX)=twatmass1(NY,NX)+VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE
      ENDDO D131

      if(fixWaterLevel)then
        dwat=twatmass0(NY,NX)-twatmass1(NY,NX)
      else
        dwat=twatmass0(NY,NX)-twatmass1(NY,NX)+Qinflx2Soil_col(NY,NX)+QLaterFlow2Cell_col(NY,NX)-QDischarg2WTBL_col(NY,NX) &
          + TWaterPlantRoot2SoilX_col(NY,NX)-QDrain_col(NY,NX)
      endif

      dwat1=twatmasstM(NY,NX)-twatmass1(NY,NX)+Qinflx2SoilM_col(NY,NX)+QLaterFlow2CellM_col(NY,NX)-QDischarM_col(NY,NX) &
          + TWaterPlantRoot2SoilXM_col(NY,NX)-QDrainM_col(NY,NX)

      if(abs(dwat1)>1.e-4_r8)then
        write(994,*)('-',L=1,100)
        write(994,*)I*1000+J/24.,NY,NX,NUM_col(NY,NX),NU_col(NY,NX)
        write(994,*)'dwat.    =',dwat1,twatmasstM(NY,NX)-twatmass1(NY,NX)
        write(994,*)'infl.    =',Qinflx2SoilM_col(NY,NX)
        write(994,*)'laterflow=',QLaterFlow2CellM_col(NY,NX)
        write(994,*)'discharg2=',QDischarM_col(NY,NX)
        write(994,*)'qdrain.  =',QDrainM_col(NY,NX)
        write(994,*)'root2soi =',TWaterPlantRoot2SoilXM_col(NY,NX)              

        write(996,*)'dwat.    =',dwat1,twatmasstM(NY,NX)-twatmass1(NY,NX)
        write(996,*)'infl.    =',Qinflx2SoilM_col(NY,NX)
        write(996,*)'laterflow=',QLaterFlow2CellM_col(NY,NX)
        write(996,*)'disx2wtbl=',QDischarM_col(NY,NX)
        write(996,*)'qdrain.  =',QDrainM_col(NY,NX)
        write(996,*)'root2soi =',TWaterPlantRoot2SoilXM_col(NY,NX)              

        write(995,*)('-',L=1,100)
        write(995,*)I*1000+J/24.,NY,NX,NUM_col(NY,NX),NU_col(NY,NX)
        DO L=NUM_col(NY,NX),NL_col(NY,NX)
          write(995,*)L,VLWatMicP1_vr(L,NY,NX),VLWatMacP1_vr(L,NY,NX),VLiceMicP1_vr(L,NY,NX),VLiceMacP1_vr(L,NY,NX)
        ENDDO
        call endrun('soil H2O error test failure in '//trim(mod_filename)//' at line',__LINE__)        
      endif
      if(abs(dwat)>1.e-4_r8)then
        write(993,*)I*1000+J/24.,M,NY,NX, NUM_col(NY,NX),NU_col(NY,NX)
        write(993,*)'watbegend=',twatmass0(NY,NX),twatmass1(NY,NX)
        write(993,*)'dwat.    =',dwat,twatmass0(NY,NX)-twatmass1(NY,NX)
        write(993,*)'infl.    =',Qinflx2Soil_col(NY,NX)
        write(993,*)'laterflow=',QLaterFlow2Cell_col(NY,NX)
        write(993,*)'disx2wtbl=',QDischarg2WTBL_col(NY,NX)
        write(993,*)'qdrain.  =',QDrain_col(NY,NX)
        write(993,*)'root2soi =',TWaterPlantRoot2SoilX_col(NY,NX)      
 
        if(abs(dwat)>1.e-4_r8)call endrun('soil H2O error test failure in '//trim(mod_filename)//' at line',__LINE__)
      endif
      twatmasstM(NY,NX) = twatmass1(NY,NX)
!      twatmass0(NY,NX)=twatmass1(NY,NX)
    ENDDO
  ENDDO  
  call PrintInfo('end '//subname)
  end subroutine ExitMassBalanceM
!------------------------------------------------------------------------------------------  

  subroutine PrepHydroThermIterM(M,NHW,NHE,NVN,NVS,TopLayWatVol_col,NUX0)
  !
  !make a copy of toplayer water content
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: TopLayWatVol_col(JY,JX)
  integer, intent(out) :: NUX0(JY,JX)
  character(len=*), parameter :: subname='PrepHydroThermIterM'
  integer :: NY,NX,L,N

  call PrintInfo('beg '//subname)
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      NUX0(NY,NX)                 = NUM_col(NY,NX)
      QDischarM_col(NY,NX)        = 0._r8
      QDrainM_col(NY,NX)          = 0._r8
      Qinflx2SoilM_col(NY,NX)     = 0._r8
      QLaterFlow2CellM_col(NY,NX) = 0._r8
      TopLayWatVol_col(NY,NX)= VLWatMicP1_vr(NUM_col(NY,NX),NY,NX)     
    ENDDO
  ENDDO

  TMLiceThawedMicP_vr           = 0.0_r8
  TMLiceThawedMacP_vr           = 0.0_r8
  TLPhaseChangeHeat2Soi1_vr     = 0.0_r8
  WatNetFlow2MicP_3DM_vr        = 0.0_r8
  WatNetFlow2MicptX_3DM_vr      = 0.0_r8
  WatNetFlow2Macpt_3DM_vr       = 0.0_r8
  THeatFlow2Soil_3DM_vr         = 0.0_r8
  WaterFlow2Micpt_3D            = 0.0_r8
  WaterFlow2MicptX_3D           = 0.0_r8
  WaterFlow2Macpt_3D            = 0.0_r8
  HeatFlow2Soili_3D             = 0.0_r8
  WaterFlow2MicPM_3D(M,:,:,:,:) = 0.0_r8
  WaterFlow2MacPM_3D(M,:,:,:,:) = 0.0_r8


  call PrintInfo('end '//subname)
  end subroutine PrepHydroThermIterM

!------------------------------------------------------------------------------------------  

  subroutine CopySoilWatVolIterateM(I,J,M,NHW,NHE,NVN,NVS,TopLayWatVol_col)

  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(in) :: TopLayWatVol_col
  integer :: L,NY,NX
  
  !VOLW2 will be updated in soil surface model

  DO NX=NHW,NHE
    DO  NY=NVN,NVS  
      DO L=NUM_col(NY,NX)+1,NL_col(NY,NX)        
        VLWatMicP2_vr(L,NY,NX)=VLWatMicP1_vr(L,NY,NX)
      ENDDO  
      VLWatMicP2_vr(NUM_col(NY,NX),NY,NX)=TopLayWatVol_col(NY,NX)

    ENDDO
  ENDDO
  end subroutine CopySoilWatVolIterateM

!------------------------------------------------------------------------------------------  
  subroutine BeginMassCheck(NHW,NHE,NVN,NVS,twatmass0)
  implicit none
  integer,  intent(in) :: NHW,NHE,NVN,NVS
  real(r8), intent(out) :: twatmass0(JY,JX)  
  character(len=*), parameter :: subname='BeginMassCheck'
  integer :: NY,NX, L

  call PrintInfo('beg '//subname)
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      !initialize for mass tracer
      twatmass0(NY,NX)                 = 0._r8
      TWaterPlantRoot2SoilX_col(NY,NX) = 0._r8
      DO L=NUM_col(NY,NX),NL_col(NY,NX)
        twatmass0(NY,NX) = twatmass0(NY,NX)+VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE
      ENDDO
    ENDDO    
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine BeginMassCheck
!------------------------------------------------------------------------------------------  


!------------------------------------------------------------------------------------------

  subroutine InitSoilHydrauics(NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: N,N1,N2,N3,N4,N5,N6
  integer :: NY,NX,L,LL
!     begin_execution
!
!     INITIALIZE SOIL HYDRAULIC PARAMETERS IN LOCAL ARRAYS
!     FOR LATER USE IN WATER TRANSFER ALGORITHMS
!
!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  call RebuildSoilFaceMap(NHW,NHE,NVN,NVS)
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      D35: DO L=NUM_col(NY,NX),NL_col(NY,NX)
        dWaterPlantRoot2SoilPrev_vr(L,NY,NX)=TWaterPlantRoot2SoilPrev_vr(L,NY,NX)*dts_HeatWatTP
        dHeatPlantRoot2SoilPrev_vr(L,NY,NX) =THeatPlantRoot2SoilPrev_vr(L,NY,NX)*dts_HeatWatTP
        THeatPlantRoot2SoilPrev_vr(L,NY,NX)=0._r8
        if(VLSoilPoreMicP_vr(L,NY,NX).LE.ZEROS2(NY,NX))cycle D35
        D40: DO N=FlowDirIndicator_col(NY,NX),3
          N1=NX;N2=NY;N3=L
! in the EW direction
          IF(N.EQ.iWestEastDirection)THEN
            IF(NX.EQ.NHE)THEN
              cycle
            ELSE
              N4=NX+1
              N5=NY
              N6=L
            ENDIF
! in the NS direction
          ELSEIF(N.EQ.iNorthSouthDirection)THEN
            IF(NY.EQ.NVS)THEN
              cycle
            ELSE
              N4=NX
              N5=NY+1
              N6=L
            ENDIF
! in the vertical direction
          ELSEIF(N.EQ.iVerticalDirection)THEN
            IF(L.EQ.NL_col(NY,NX))THEN
              cycle
            ELSE
              N4=NX
              N5=NY
              N6=L+1
            ENDIF
          ENDIF
          ! Match the destination layer used by Subsurface3DInternalFlowM.
          do LL=N6,NL_col(N5,N4)
            if(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))then
              N6=LL
              exit
            endif
          enddo
          if(LL>NL_col(N5,N4))cycle D40
          !
          !     MACROPORE CONDUCTIVITY FROM 'HOUR1' AND GRAVITATIONAL
          !     GRADIENT USED TO CALCULATE MACROPORE FLOW FOR USE BELOW
          !
          !     HydroCondMacP1=macropore hydraulic conductivity
          !     AVCNHL_3D=macropore hydraulic conductance
          !     DLYR=layer depth
          !
          IF(HydroCondMacP1_vr(N3,N2,N1).GT.ZERO .AND. HydroCondMacP1_vr(N6,N5,N4).GT.ZERO)THEN
            !when both src and dest have macro pores
            AVCNHL_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)=2.0_r8*HydroCondMacP1_vr(N3,N2,N1)*HydroCondMacP1_vr(N6,N5,N4) &
              /(HydroCondMacP1_vr(N3,N2,N1)*DLYR_3D(N,N6,N5,N4)+HydroCondMacP1_vr(N6,N5,N4)*DLYR_3D(N,N3,N2,N1))
          ELSE
            AVCNHL_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)=0.0_r8
          ENDIF
        ENDDO D40
      ENDDO D35
    ENDDO
  ENDDO
  end subroutine InitSoilHydrauics

!------------------------------------------------------------------------------------------

  subroutine Subsurface3DInternalFlowM(I,J,M,NHW,NHE,NVN,NVS,KSatRedusByRainKinetEnergy,HeatFluxAir2Soi)
  !
  !Description
  !3D flow between internal grids 
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in)  :: M,NHW,NHE,NVN,NVS
  real(r8), dimension(:,:),intent(in) :: KSatRedusByRainKinetEnergy(:,:)
  real(r8), dimension(:,:),intent(in) :: HeatFluxAir2Soi(:,:)

  character(len=*), parameter :: subname='Subsurface3DInternalFlowM'
  integer :: N,N1,N2,N3,N4,N5,N6,L,LL,K1,KL,NY,NX
  real(r8) :: WTHET1,FCDX,FCLX,FCX
  real(r8) :: PSISV1,TKY,PSDX
  real(r8) :: WPLX,WPX,WTHET2

  real(r8) :: ConvectVapFlux,HeatByConvectVapFlux,PSISVL
  real(r8) :: PredDarcyFlowMax
  real(r8) :: TKLX
  real(r8) :: WatDarcyFlowMicP,HeatByDarcyFlowMicP,WaterMacpFlow
  real(r8) :: HeatByFlowMacP,HeatByWatFlowMicP,HFLWS,THETA1,THETAL  
  logical  :: LInvalidMacP     !to disable/enable macropore
  logical  :: checkLayerSrc,checkLayerTgt

  !     begin_execution

  call PrintInfo('beg '//subname)
  !
  !     WATER AND ENERGY TRANSFER THROUGH SOIL PROFILE
  !
  !     N3,N2,N1=L,NY,NX of source grid cell
  !     N6,N5,N4=L,NY,NX of destination grid cell
  !

  call InitSoil3DModelIterateM(I,J,M,NHW,NHE,NVN,NVS)
  
  !If waterlevel is fixed, like shallow lake, without considering the change
  !of water depth due to hydrological fluxes 

  DO NX=NHW,NHE    !sweep from west to east
    DO NY=NVN,NVS  !sweep from north to south

      D4400: DO L=1,NL_col(NY,NX)  !sweep from top to bottom
        N1=NX;N2=NY;N3=L
        !
        !     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
        !
        D4320: DO N=FlowDirIndicator_col(N2,N1),3
          IF(N.EQ.iWestEastDirection)THEN            
            !west to east
            IF(NX.EQ.NHE)THEN  !skip eastern boundary              
              cycle
            ELSE
              N4=NX+1;N5=NY;N6=L  
            ENDIF
          ELSEIF(N.EQ.iNorthSouthDirection)THEN
            !north to south
            IF(NY.EQ.NVS)THEN  !skip southern boundary
              cycle
            ELSE
              N4=NX;N5=NY+1;N6=L
            ENDIF
          ELSEIF(N.EQ.iVerticalDirection)THEN
            !bottom layer
            IF(L.EQ.NL_col(NY,NX))THEN !skip the bottom boundary
              cycle
            ELSE
              !destination is layer below
              N4=NX;N5=NY;N6=L+1
            ENDIF
          ENDIF
          !
          !'SKIP NON-EXISTENT DESTINATION SOIL LAYERS'
          ! identified by soil volume
          D1100: DO LL=N6,NL_col(N5,N4)
            IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
              N6=LL
              exit
            ENDIF
          ENDDO D1100
          ! No internal face exists without an active destination layer.
          if(LL>NL_col(N5,N4))cycle D4320

          ! source grid at surface, make a copy of the destination layer number
          IF(N3.EQ.NU_col(N2,N1))N6X(N2,N1)=N6
          !
          !     POROSITIES 'THETP*', WATER CONTENTS 'THETA*', AND POTENTIALS
          !     'PSIS*' FOR EACH GRID CELL
          !
          !     THETA1,THETAL=micropore water concn in source,destination cells
          !     THETY=hygroscopic water concentration
          !     POROS=soil porosity
          !     VLSoilPoreMicP_vrI=soil volume excluding rock, macropore
          !
          IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
            !make sure both source and destination grids are legitimate
            !and match the source top layer to the destination top layer 
            checkLayerSrc=N3.GE.NUM_col(N2,N1) .AND. N3.LE.NL_col(N2,N1)
            checkLayerTgt=N6.GE.NUM_col(N5,N4) .AND. N6.LE.NL_col(N5,N4)

            IF(checkLayerSrc .AND. checkLayerTgt)THEN

              call CalcSoilWatPotential(NY,NX,N1,N2,N3,PSISoilMatricPtmp_vr(N3,N2,N1),THETA1)

              call CalcSoilWatPotential(NY,NX,N4,N5,N6,PSISoilMatricPtmp_vr(N6,N5,N4),THETAL)
              !
              !     ACCOUNT FOR WETTING FRONTS WHEN CALCULATING WATER CONTENTS,
              !     MATRIC WATER POTENTIALS AND HYDRAULIC CONDUCTIVITIES USED
              !     IN WATER FLUX CALCULATIONS
              !
              !     HydCondSrc,CNDL=hydraulic conductivities in source,destination cells
              !     KSatRedusByRainKinetEnergy=reduction in soil surface Ksat from rainfall energy impact
              !     PSISM1=soil matric potential
              !     VLWatMicPX1=VLWatMicP1 accounting for wetting front
              !
              !     (CURRENT WATER POTENTIAL > AIR ENTRY WATER POTENTIAL)
              
              call MicropXGridDarcyFlow(N,N1,N2,N3,N4,N5,N6,THETA1,THETAL,          &
                KSatRedusByRainKinetEnergy(NY,NX),PredDarcyFlowMax,WatDarcyFlowMicP,&
                HeatByDarcyFlowMicP,PSISV1,PSISVL)          
                
              call MacropXgridFLow(N,N1,N2,N3,N4,N5,N6,WaterMacpFlow,HeatByFlowMacP,LMacPoreOff)

              call WaterVaporXgridFlow(M,N,N1,N2,N3,N4,N5,N6,PSISV1,PSISVL,ConvectVapFlux,HeatByConvectVapFlux)
          
              !     WaterFlow2MicptX_3D=total unsaturated water+vapor flux to destination
              !     HeatByWatFlowMicP=total convective heat flux from water+vapor flux          
              !flow into cell (N6,N5,N4)
              WaterFlow2Micpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)  = WatDarcyFlowMicP+ConvectVapFlux
              WaterFlow2MicptX_3D(N,SoilFaceSlot(N,N3,N6),N5,N4) = PredDarcyFlowMax+ConvectVapFlux
              WaterFlow2Macpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)  = WaterMacpFlow
              HeatFlow2Soili_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)   = HeatByDarcyFlowMicP+HeatByConvectVapFlux+HeatByFlowMacP
              WaterFlow2Micptl_3D(N,SoilFaceSlot(N,N3,N6),N5,N4) = WaterFlow2Micpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
              WaterFlow2Macptl_3D(N,SoilFaceSlot(N,N3,N6),N5,N4) = WaterFlow2Macpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)

              !   compute heat flux by conduction
              call SolveXgridHeatConduction(I,J,M,N,N1,N2,N3,N4,N5,N6,ConvectVapFlux,HeatByConvectVapFlux,HeatFluxAir2Soi(NY,NX))

              !
              !     TOTAL WATER, VAPOR AND HEAT FLUXES
              !
              !     FLW,FLWX,FLWH=total water flux through micropores,macropores
              !     HFLW=total heat flux
              !     WaterFlow2MicPM_3D=water flux used for solute flux calculations in TranspNoSalt.f
              !
              if(.not.fixWaterLevel)then
                WaterFlowSoiMicP_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)  = WaterFlowSoiMicP_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)+WaterFlow2Micpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                WaterFlowSoiMicPX_3D(N,SoilFaceSlot(N,N3,N6),N5,N4) = WaterFlowSoiMicPX_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)+WaterFlow2MicptX_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                WaterFlowSoiMacP_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)  = WaterFlowSoiMacP_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)+WaterFlow2Macpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
              endif
              HeatFlow2Soil_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)     = HeatFlow2Soil_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)+HeatFlow2Soili_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)

              WaterFlow2MicPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4) = WaterFlow2Micpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
              WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4) = WaterFlow2Macpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)

              IF(N.EQ.iVerticalDirection)THEN
                FILMM_vr(M,N6,N5,N4)=FilmThickness(PSISoilMatricPtmp_vr(N6,N5,N4))
              ENDIF
            ENDIF
          ENDIF
        ENDDO D4320
      ENDDO D4400
    ENDDO
  ENDDO  

  call PrintInfo('end '//subname)
  end subroutine Subsurface3DInternalFlowM
!------------------------------------------------------------------------------------------

  subroutine XBoundaryFlowM(I,J,M,NHW,NHE,NVN,NVS,RainEkReducedKsat)
  !
  !Description
  ! boundary flow involes exchange with external water table, and through tile drainage
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8),intent(in) :: RainEkReducedKsat(JY,JX)
  character(len=*), parameter :: subname='XBoundaryFlowM'

  integer :: NY,NX
  integer :: L,LL
  integer :: N,NN,N1,N2,N3
  integer :: N4,N5,N4B,N5B,N6           !internal grids, grid indices for snow redistribution
  integer :: M1,M2,M3,M4,M5,M6          !for boundary grids
  integer :: K1,KL
  logical :: DoMacPDischarg2Tile
  logical :: DoMicPDischarg2ExtWTBL
  logical :: DoMacPDischarg2ExtWTBL
  logical :: DoMicPDischarg2Tile
  real(r8) :: FLWT
  real(r8) :: RechargSurf
  real(r8) :: RechargDist2WTBL
  real(r8) :: Recharg2WTBLScal   !multiplier to water exchange with external water table.
  real(r8) :: Surf2MacPDepz,CNDL
  real(r8) :: THETAX  
  real(r8) :: AirfMicP,VOLPX2,AirfMacP
  real(r8) :: XN,THETA1,HeatFlx

  character(len=20) :: str_dir
  logical :: donot_drain,lZeroRunoff,checkActiveFlow
  character(len=2), parameter :: direcs(2)=(/'ew','ns'/)
  ! begin_execution
  call PrintInfo('beg '//subname)
  ! AirfMicP,AirfMacP=air-filled porosity in micropores,macropores
  D9595: DO  NX=NHW,NHE                !from west to east
    D9590: DO  NY=NVN,NVS              !from north to south
      D9585: DO L=NUM_col(NY,NX),NL_col(NY,NX) !from top to bottom

        AirfMicP      = VLMicP1_vr(L,NY,NX)-VLWatMicP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)
        VOLPX2        = AirfMicP
        AirfMacP      = VLMacroPores1_vr(L,NY,NX)-VLWatMacP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)
        Surf2MacPDepz = get_Surf2MacPDepz(L,NY,NX)

        !-----------------------------------
        !
        !     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
        !     ENTERED IN 'READS'
        !
        !     N3,N2,N1 for source grid cell
        !     M6,M5,M4 for destination grid cell
        !       
        N1=NX;N2=NY;N3=L
        !
        !     LOCATE EXTERNAL BOUNDARIES
        !
        D9580: DO N=FlowDirIndicator_col(NY,NX),3

          D9575: DO NN=1,2   !go from iFront=1, to iBehind=2
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !            !Locate the boundary 
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
            str_dir='none'       
            IF(N.EQ.iWestEastDirection)THEN
              ! along the W-E direction, NX              
              N4  = NX+1; N5  = NY  !eastern
              N4B = NX-1; N5B = NY  !western
              N6  = L
              IF(NN.EQ.iFront)THEN            !target grid is on the front of (N2,N1)    
                IF(NX.EQ.NHE)THEN             !eastern boundary, outflow  |---|->             
                  M1 = NX;  M2  = NY; M3 = L  !source
                  M4 = NX+1;M5 = NY;  M6 = L  !target
                  XN = -1.0_r8                !going out of eastern boundary
                  RechargSurf       = RechargEastSurf_col(M2,M1)
                  RechargDist2WTBL  = RechrgDistEastSubSurf_col(M2,M1)      !distance to eastern external water table [m]
                  Recharg2WTBLScal  = RechargRateEastWTBL_col(M2,M1)  
                  str_dir='east'
                  !do nothing if not on the boundary  
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN                
                IF(NX.EQ.NHW)THEN            !western boundary, inflow   -|-> |               
                  M4 = NX; M5 = NY; M6 = L   !target on the boundary
                  XN = 1.0_r8                !coming in from western boundary
                  RechargSurf       = RechargWestSurf_col(M5,M4)
                  RechargDist2WTBL  = RechrgDistWestSubSurf_col(M5,M4)      !distance to western external water table [m]
                  Recharg2WTBLScal  = RechargRateWestWTBL_col(M5,M4)
                  str_dir='west'
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              ! along the N-S direction          , NY    
              N4  = NX; N5  = NY+1 !south
              N4B = NX; N5B = NY-1 !north
              N6  = L
              IF(NN.EQ.iFront)THEN             !-----                
                IF(NY.EQ.NVS)THEN              !southern boundary, outflow     \|/ !-----
                  M1 = NX; M2 = NY; M3   = L   !source
                  M4 = NX; M5 = NY+1; M6 = L   !target
                  XN                = -1.0_r8    !going out of south boundary
                  RechargSurf       = RechargSouthSurf_col(M2,M1)
                  RechargDist2WTBL  = RechrgDistSouthSubSurf_col(M2,M1)   !distance to southern external water table [m]
                  Recharg2WTBLScal  = RechargRateSouthWTBL_col(M2,M1)
                  str_dir='south'
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN      !current grid (NY,NX) is behind 
                IF(NY.EQ.NVN)THEN            !northern boundary, inflow     ----
                  M4 = NX;M5  = NY; M6   = L !target grid on the boundary
                  XN          = 1.0_r8       !coming in from north boundary
                  RechargSurf      = RechargNorthSurf_col(M5,M4)
                  RechargDist2WTBL = RechrgDistNorthSubSurf_col(M5,M4)    !distance to northern external water table [m]
                  Recharg2WTBLScal = RechargRateNorthWTBL_col(M5,M4)
                  str_dir='north'
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN              
              N4 = NX;N5 = NY;N6 = L+1
              IF(NN.EQ.iFront)THEN
                !bottom
                IF(L.EQ.NL_col(NY,NX))THEN
                  M1 = NX; M2 = NY; M3  = L   !source
                  M4 = NX; M5 = NY; M6 = L+1  !target
                  XN              = -1.0_r8    !going out lower boundary
                  RechargDist2WTBL  = 1._r8
                  Recharg2WTBLScal = RechargBottom_col(M2,M1)
                  str_dir='vert'
                ELSE
                  cycle
                ENDIF
              !surface infiltration is done in another case
              ELSEIF(NN.EQ.iBehind)THEN
                cycle
              ENDIF
            ENDIF
            !-----------------------------------
            ! configure water table info, should the exchange to water table be limited only to boundary grids?
            call Config4XWaterTableDrain(I,J,N,L,NY,NX,XN,Surf2MacPDepz,DoMicPDischarg2ExtWTBL,DoMacPDischarg2ExtWTBL)

            call Config4TileDrainage(I,J,N,L,NY,NX,XN,Surf2MacPDepz,DoMicPDischarg2Tile,DoMacPDischarg2Tile)

            !
            !     BOUNDARY SURFACE RUNOFF DEPENDING ON ASPECT, SLOPE
            !     VELOCITY, HYDRAULIC RADIUS AND SURFACE WATER STORAGE
            !
            !     CDPTH,CumLitRDepzInit_col=current,initial surface elevation
            !     BKDS=bulk density
            !     XGridRunoffFlag_2DH,RCHQ*=runoff boundary flags

            !     top soil layer and surface soil layer is active, litter layer is lower than its initial thickness 
            !     or the grid is a soil
            !    surface lateral flow 
            IF(L.EQ.NUM_col(N2,N1) .AND. N.NE.iVerticalDirection                                & ! do lateral flow at surface
              .AND. (CumDepz2LayBottom_vr(NU_col(N2,N1)-1,N2,N1).LE.CumLitRDepzInit_col(N2,N1)  & ! in the soil
              .OR. SoilBulkDensity_vr(NUI_col(N2,N1),N2,N1).GT.ZERO))THEN                         ! it is soil, not a pond

              !  NO runoff
              lZeroRunoff=.not.XGridRunoffFlag_2DH(NN,N,N2,N1)   &              !Runoff flag is off
                .OR. isclose(RechargSurf,0._r8)                                 !Allowed runoff rate is zero

              !do runoff
              IF(.not.LZeroRunoff)then
                call XBoundSurfaceRunoffM(I,J,M,N,NN,N1,N2,M4,M5,RechargSurf,XN)
              ENDIF
            ENDIF
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! BOUNDARY SUBSURFACE WATER AND HEAT TRANSFER DEPENDING
            ! ON LEVEL OF WATER TABLE
            !
            IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(NY,NX))THEN       !current grid has legitimate soil pore
              checkActiveFlow=FlowDirIndicator_col(N2,N1).NE.iVerticalDirection .OR. N.EQ.iVerticalDirection

              IF(checkActiveFlow)THEN                  
                !IF NO WATER TABLE THEN FLUX IS THROUGH BOTTOM BOUNDARY                
                IF((IDWaterTable_col(N2,N1).EQ.0 .OR. N.EQ.iVerticalDirection) .AND. .not.isclose(Recharg2WTBLScal,0.0_r8))THEN    
                  !handle vertical flow
                  call VertBoundaryDrainM(I,J,N,N1,N2,N3,M4,M5,M6,RainEkReducedKsat(NY,NX),XN,RechargDist2WTBL,Recharg2WTBLScal)
                  !
                  !lateral flow
                ELSE

                  if((DoMicPDischarg2ExtWTBL .or. DoMacPDischarg2ExtWTBL) .and. (.not.isclose(Recharg2WTBLScal,0._r8)))then
                    !discharge to external water table
                    CALL DischargeOverWaterTable(I,J,N,N1,N2,N3,M4,M5,M6,DoMicPDischarg2ExtWTBL,DoMacPDischarg2ExtWTBL,&
                      RechargDist2WTBL,Recharg2WTBLScal,Surf2MacPDepz,XN,ExtWaterTable_col(N2,N1),FracLayVolBelowExtWTBL_vr(N3,N2,N1))
                  endif

                  if((DoMicPDischarg2Tile .or. DoMacPDischarg2Tile) .and. (.not.isclose(Recharg2WTBLScal,0._r8)))then
                    !tile drainage
                    CALL DischargeOverWaterTable(I,J,N,N1,N2,N3,M4,M5,M6,DoMicPDischarg2Tile,DoMacPDischarg2Tile,&
                      RechargDist2WTBL,Recharg2WTBLScal,Surf2MacPDepz,XN,TileWaterTable_col(N2,N1),FracLayVolBelowTileWTBL_vr(N3,N2,N1))                      
                  endif

                  call RechargeFromExtWaterTBL(I,J,M,N,N1,N2,N3,M4,M5,M6,Surf2MacPDepz,RechargDist2WTBL,&
                    Recharg2WTBLScal,XN,AirfMicP,VOLPX2,AirfMacP,str_dir)!,donot_drain)

                ENDIF
                !
                !     SUBSURFACE HEAT SOURCE/SINK
                !     HFLWL=heat flux across lower boundary
                !     TK1=lower boundary soil temperature
                !     TKSD=deep source/sink temperature from geothermal flux
                !     TCNDG=thermal conductivity below lower boundary
                !     SoilHeatSrcDepth_col,CDPTH=depth of thermal sink/source, lower boundary
                !     KoppenClimZone=Koppen climate zone
                IF(N.EQ.iVerticalDirection .AND. KoppenClimZone_col(N2,N1).NE.-2)THEN
                  !heat flux going out (<0)
                  HeatFlx=(TKSoil1_vr(N3,N2,N1)-TKSD_col(N2,N1))*TCNDG/(SoilHeatSrcDepth_col(N2,N1)-CumDepz2LayBottom_vr(N3,N2,N1)) &
                    *AREA_3D(N,N3,N2,N1)*dts_HeatWatTP
                  HeatFlow2Soili_3D(N,M6,M5,M4) = HeatFlow2Soili_3D(N,M6,M5,M4)+heatFlx                  
                  HeatSource_col(N2,N1)         = HeatSource_col(N2,N1)-HeatFlx
                ENDIF

                if(.not.fixWaterLevel)then
                  WaterFlowSoiMicP_3D(N,M6,M5,M4)  = WaterFlowSoiMicP_3D(N,M6,M5,M4)+WaterFlow2Micpt_3D(N,M6,M5,M4)
                  WaterFlowSoiMicPX_3D(N,M6,M5,M4) = WaterFlowSoiMicPX_3D(N,M6,M5,M4)+WaterFlow2MicptX_3D(N,M6,M5,M4)
                  WaterFlowSoiMacP_3D(N,M6,M5,M4)  = WaterFlowSoiMacP_3D(N,M6,M5,M4)+WaterFlow2Macpt_3D(N,M6,M5,M4)
                  WaterFlow2MicPM_3D(M,N,M6,M5,M4) = WaterFlow2Micpt_3D(N,M6,M5,M4)
                  WaterFlow2MacPM_3D(M,N,M6,M5,M4) = WaterFlow2Macpt_3D(N,M6,M5,M4)
                endif
                HeatFlow2Soil_3D(N,M6,M5,M4)     = HeatFlow2Soil_3D(N,M6,M5,M4)+HeatFlow2Soili_3D(N,M6,M5,M4)                
              ENDIF
            ENDIF
          ENDDO D9575  !loop NN
          !
          !     NET WATER AND HEAT FLUXES IN RUNOFF AND SNOW DRIFT
          !
          !     TQR1,THQR1=net runoff,convective heat from runoff
          !     TQS1,TQW1,TQI1,THQS1=net snow,water,ice, heat from snowpack runoff
          !     WatFlx2LitRByRunoff,HeatFlx2LitRByRunoff=runoff, convective heat from runoff
          !     DrySnoFlxBySnowRedistribut,WatFlxBySnowRedistribut,IceFlxBySnowRedistribut=snow,water,ice transfer
          !     HeatFlxBySnowRedistribut=convective heat transfer from snow,water,ice transfer

          IF(L.EQ.NUM_col(N2,N1) .AND. N.NE.iVerticalDirection)THEN
            if(snowRedist_model)call SumSnowDriftByRunoffM(M,N,N1,N2,N4,N5,N4B,N5B)
          ENDIF

        ENDDO D9580  !loop boundaries through 3 directions

        !
        !     DISSIPATE WETTING FRONT
        !
        !     VLWatMicP1=soil micropore water content
        !     VLWatMicPX1=soil micropore water content behind wetting front
        !     FLWVL=water flux from wetted to drier soil
        !
      ENDDO D9585  !loop vertically through layers
    ENDDO D9590
  ENDDO D9595
  call PrintInfo('end '//subname)
  end subroutine XBoundaryFlowM
!------------------------------------------------------------------------------------------
  subroutine Summarize3DFlowM(I,J,M,NHW,NHE,NVN,NVS,Qinfl2MicP_col,HeatInfl2Soil,Qinfl2MacP_col)
  ! Description:
  ! Assemble accepted soil face transfers and freeze/thaw, tallying internal
  ! lateral water and heat separately from external drainage and vertical flow.
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: Qinfl2MicP_col(JY,JX)
  real(r8),intent(in) :: HeatInfl2Soil(JY,JX)
  real(r8),intent(in) :: Qinfl2MacP_col(JY,JX)  
  character(len=*), parameter :: subname='Summarize3DFlowM'
  integer :: NY,NX,L,N,LL
  integer :: first_face,last_face
  integer :: N3,N2,N1
  integer :: N6,N5,N4
  real(r8):: FINHX
  real(r8) :: VLWatMicP1X
  real(r8) :: VOLP1X,VLWatMacP1X,VOLPH1X  
  logical :: CHECKGRID
  character(len=2), parameter :: direcs(2)=(/'ew','ns'/)

  ! begin_execution
  call PrintInfo('beg '//subname)
  
  D95951: DO NX=NHW,NHE                !from west to east
    D95901: DO NY=NVN,NVS              !from north to south
      HeatFlx2Grnd_col(NY,NX) = HeatFlx2Grnd_col(NY,NX)+HeatInfl2Soil(NY,NX)
      Qinflx2Soil_col(NY,NX)  = Qinflx2Soil_col(NY,NX)+Qinfl2MicP_col(NY,NX)+Qinfl2MacP_col(NY,NX)
      Qinflx2SoilM_col(NY,NX) = Qinflx2SoilM_col(NY,NX)+Qinfl2MicP_col(NY,NX)+Qinfl2MacP_col(NY,NX)

      D95851: DO L=NUM_col(NY,NX),NL_col(NY,NX) !from top to bottom

        N1=NX;N2=NY;N3=L
        !
        !     LOCATE EXTERNAL BOUNDARIES
        !
        D95801: DO N=FlowDirIndicator_col(NY,NX),3
          IF(N.EQ.iWestEastDirection)THEN
            ! along the W-E direction, NX              
            N6=L; N4=NX+1;N5  = NY  !eastern
            
          ELSEIF(N.EQ.iNorthSouthDirection)THEN
            ! along the N-S direction , NY    
            N6=L;N4=NX; N5=NY+1 !south
            
          ELSEIF(N.EQ.iVerticalDirection)THEN              
            N4 = NX;N5 = NY;N6 = L+1
          ENDIF  
          !
          !     NET WATER AND HEAT FLUXES THROUGH SOIL AND SNOWPACK
          !
          !     TFLWL,THFLWL=net water micropore,macropore flux
          !     THFLWL=net convective+conductive heat flux
          !     FLWL =micropore water,heat flux
          !     WaterFlow2Macpt_3D=macropore water,heat flux
          !     HFLWL=soil heat flux
          !
          IF(FlowDirIndicator_col(N2,N1).NE.iVerticalDirection .OR. N.EQ.iVerticalDirection)THEN
            D12003: DO LL=N6,NL_col(N5,N4)
              IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
                N6=LL
                exit
              ENDIF
            ENDDO D12003

  call SoilFaceIncomingRange(N,N3,N2,N1,first_face,last_face)
            !compute the net flux 
            IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
              if(.not.fixWaterLevel)then
                WatNetFlow2MicP_3DM_vr(N3,N2,N1)   = WatNetFlow2MicP_3DM_vr(N3,N2,N1)+sum(WaterFlow2Micpt_3D(N,first_face:last_face,N2,N1))-WaterFlow2Micpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                WatNetFlow2MicptX_3DM_vr(N3,N2,N1) = WatNetFlow2MicptX_3DM_vr(N3,N2,N1)+sum(WaterFlow2MicptX_3D(N,first_face:last_face,N2,N1))-WaterFlow2MicptX_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                WatNetFlow2Macpt_3DM_vr(N3,N2,N1)  = WatNetFlow2Macpt_3DM_vr(N3,N2,N1)+sum(WaterFlow2Macpt_3D(N,first_face:last_face,N2,N1))-WaterFlow2Macpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)

              endif
              THeatFlow2Soil_3DM_vr(N3,N2,N1)    = THeatFlow2Soil_3DM_vr(N3,N2,N1)+sum(HeatFlow2Soili_3D(N,first_face:last_face,N2,N1))-HeatFlow2Soili_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)

              if(N.NE.iVerticalDirection)then     

                CHECKGRID=(N1.GT.NHW .AND. N1.LT.NHE) .AND. (N2.GT.NVN .AND. N2.LT.NVS) .OR.                      & !innter grid
                  (N.EQ.iWestEastDirection .AND. N1.GT.NHW .AND. N1.LT.NHE .AND. (N2.EQ.NVN .OR. N2.EQ.NVS)) .OR. & !north/south boundary
                  (N.EQ.iNorthSouthDirection .AND. N2.GT.NVN .AND. N2.LT.NVS .AND. (N1.EQ.NHW .OR. N1.EQ.NHE))      !west/east boundary

                IF((N.EQ.iWestEastDirection .AND. NX.EQ.NHW .AND. NHW.LT.NHE) .OR. & !western boundary
                  (N.EQ.iNorthSouthDirection .AND. NY.EQ.NVN .AND. NVN.LT.NVS)) THEN !northern boundary 
                  ! External boundary heat is already in the drainage ledger.
                  ! Only exchanges with another active column enter this tally.
                  HeatLateral2Cell_col(N2,N1)=HeatLateral2Cell_col(N2,N1) &
                    -HeatFlow2Soili_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                  QLaterFlow2Cell_col(N2,N1)  = QLaterFlow2Cell_col(N2,N1)-WaterFlow2Micpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)-WaterFlow2Macpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                  QLaterFlow2CellM_col(N2,N1) = QLaterFlow2CellM_col(N2,N1)-WaterFlow2Micpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)-WaterFlow2Macpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                  
                ELSEIF((N.EQ.iWestEastDirection .AND. NX.EQ.NHE .AND. NHW.LT.NHE) .OR. & !eastern boundary
                  (N.EQ.iNorthSouthDirection .AND. NY.EQ.NVS .AND. NVN.LT.NVS)) THEN !southern boundary 
                  HeatLateral2Cell_col(N2,N1)=HeatLateral2Cell_col(N2,N1) &
                    +sum(HeatFlow2Soili_3D(N,first_face:last_face,N2,N1))
                  QLaterFlow2Cell_col(N2,N1)  = QLaterFlow2Cell_col(N2,N1)+sum(WaterFlow2Micpt_3D(N,first_face:last_face,N2,N1))+sum(WaterFlow2Macpt_3D(N,first_face:last_face,N2,N1))
                  QLaterFlow2CellM_col(N2,N1) = QLaterFlow2CellM_col(N2,N1)+sum(WaterFlow2Micpt_3D(N,first_face:last_face,N2,N1))+sum(WaterFlow2Macpt_3D(N,first_face:last_face,N2,N1))
                ELSEIF(CHECKGRID)THEN
                  HeatLateral2Cell_col(N2,N1)=HeatLateral2Cell_col(N2,N1) &
                    +sum(HeatFlow2Soili_3D(N,first_face:last_face,N2,N1)) &
                    -HeatFlow2Soili_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                  QLaterFlow2Cell_col(N2,N1)  = QLaterFlow2Cell_col(N2,N1)+sum(WaterFlow2Micpt_3D(N,first_face:last_face,N2,N1))-WaterFlow2Micpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4) &
                    +sum(WaterFlow2Macpt_3D(N,first_face:last_face,N2,N1))-WaterFlow2Macpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                  QLaterFlow2CellM_col(N2,N1) = QLaterFlow2CellM_col(N2,N1)+sum(WaterFlow2Micpt_3D(N,first_face:last_face,N2,N1))-WaterFlow2Micpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4) &
                    +sum(WaterFlow2Macpt_3D(N,first_face:last_face,N2,N1))-WaterFlow2Macpt_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
                  
                ENDIF
              endif
            ENDIF
          ENDIF
        ENDDO D95801
        !
        !     INFILTRATION OF WATER FROM MACROPORES INTO MICROPORES
        !
        !     VLWatMacP1=macropore volume
        !     FINHX,FINHL=macro-micropore transfer unltd,ltd by water,air volume
        !     FWatExMacP2MicPM_vr=macro-micropore transfer for use in TranspNoSalt.f
        !     HydroCond_3D=hydraulic conductivity
        !     PSISE,PSISoilAirEntry_vr=air entry,matric water potentials
        !     PHOL,MacPRadius=path length between,radius of macropores from hour1.f
        !     dts_HeatWatTP=time step
        !     VLWatMicP1X,VOLP1X=current micropore water,air volume
        !     VLWatMacP1X,VOLPH1X=current macropore water,air volume
        !     PSI_at_Sautration_vr < 0,    
        IF(VLWatMacP1_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
          !south-north direction, is it true?
          FINHX = TwoPiCON*HydroCond_3D(2,1,N3,N2,N1)*AREA_3D(3,N3,N2,N1) &
            *(PSI_at_Sautration_vr(N3,N2,N1)-PSISoilMatricPtmp_vr(N3,N2,N1)) &
            /LOG(PathLenMacPore_vr(N3,N2,N1)/MacPoreRadius_vr(N3,N2,N1))*dts_HeatWatTP
          VLWatMicP1X = VLWatMicP1_vr(N3,N2,N1)+WatNetFlow2MicP_3DM_vr(N3,N2,N1)+FWatIrrigate2MicP1_vr(N3,N2,N1)
          VOLP1X      = AZMAX1(VLMicP1_vr(N3,N2,N1)-VLWatMicP1X-VLiceMicP1_vr(N3,N2,N1))
          VLWatMacP1X = VLWatMacP1_vr(N3,N2,N1)+WatNetFlow2Macpt_3DM_vr(N3,N2,N1)
          VOLPH1X     = AZMAX1(VLMacroPores1_vr(N3,N2,N1)-VLWatMacP1X-VLiceMacP1_vr(N3,N2,N1))

          IF(FINHX.GT.0.0_r8)THEN
            FWatExMacP2MicPiM_vr(N3,N2,N1)=AZMAX1(AMIN1(FINHX,VLWatMacP1X,VOLP1X))
          ELSE
            FWatExMacP2MicPiM_vr(N3,N2,N1)=AZMIN1(AMAX1(FINHX,-VOLPH1X,-VLWatMicP1X))
          ENDIF
          FWatExMacP2MicPM_vr(M,N3,N2,N1) = FWatExMacP2MicPiM_vr(N3,N2,N1)
          FWatExMacP2MicP_vr(N3,N2,N1)    = FWatExMacP2MicP_vr(N3,N2,N1)+FWatExMacP2MicPiM_vr(N3,N2,N1)
        ELSE
          FWatExMacP2MicPiM_vr(N3,N2,N1)   = 0.0_r8
          FWatExMacP2MicPM_vr(M,N3,N2,N1) = 0.0_r8
        ENDIF

        call FreezeThawIterateM(I,J,M,NY,NX,L,N1,N2,N3)
      ENDDO D95851
    ENDDO D95901
  ENDDO D95951
  call PrintInfo('end '//subname)
  end subroutine Summarize3DFlowM
!------------------------------------------------------------------------------------------

  subroutine VertBoundaryDrainM(I,J,N,N1,N2,N3,M4,M5,M6,RainEkReducedKsat,XN,RechargDist2WTBL,Recharg2WTBLScal)    
  !
  !Drainage across the lateral or bottom boundaries.      
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: N
  integer,  intent(in) :: N1,N2,N3  !source
  integer,  intent(in) :: M4,M5,M6  !dest, when XN=1, N2=M5,N1=M4
  real(r8), intent(in) :: RainEkReducedKsat 
  real(r8), intent(in) :: XN        !direction
  real(r8), intent(in) :: RechargDist2WTBL    !lateral distance to external water table, [m]
  real(r8), intent(in) :: Recharg2WTBLScal

  character(len=*), parameter :: subname='VertBoundaryDrainM'
  real(r8) :: THETA1,HydcondSrc,RechargRate
  real(r8) :: watflx,heatflx,HydGrad,VOLWZ
  integer :: K1 !saturation index, saturate=1
!

  call PrintInfo('beg '//subname)
! IDWaterTable=water table flag

  RechargRate=Recharg2WTBLScal*dts_HeatWatTP*RechargDist2WTBL
  !involve no water table or in vertical direction
  !water content
  THETA1 = AMAX1(SoilWatAirDry_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),safe_adb(VLWatMicP1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))
  K1     = getMoistK(THETA1,POROS_vr(N3,N2,N1))

  !compute hydraulic conductivity for micropore
  HydcondSrc=HydroCond_3D(N,K1,N3,N2,N1)
  IF(N3.EQ.NUM_col(N2,N1))THEN  !on the surface for WE, NS boundaries.
    HydcondSrc=HydcondSrc*RainEkReducedKsat    
  ENDIF
  HydCondSoil_3D(N,N3,N2,N1) = HydcondSrc

  !x-section area scaled hydraulic gradient, HydGrad > 0, when XN=-1, ES; HydGrad < 0, when XN=1, WN
  HydGrad=XN*mGravAccelerat*(-ABS(SLOPE_col(N,N2,N1)))*AREA_3D(3,N3,N2,N1)

  VOLWZ=AZMAX1(VLWatMicP1_vr(N3,N2,N1)*dts_wat) 

  WaterFlow2Micpt_3D(N,M6,M5,M4)=AZERO(AMIN1(VOLWZ, AMAX1(-VOLWZ,RechargRate*HydGrad*HydcondSrc)))

  VOLWZ=AZMAX1(VLWatMacP1_vr(N3,N2,N1)*dts_wat) 
  WaterFlow2MicptX_3D(N,M6,M5,M4) = WaterFlow2Micpt_3D(N,M6,M5,M4)  
  WaterFlow2Macpt_3D(N,M6,M5,M4)  = AZERO(AMIN1(VOLWZ,AMAX1(-VOLWZ,RechargRate*HydGrad*HydroCondMacP1_vr(N3,N2,N1))))
    
  watflx  = WaterFlow2Micpt_3D(N,M6,M5,M4)+WaterFlow2Macpt_3D(N,M6,M5,M4)

  !from (N3,N2,N1) to (M6,M5,M4)
  if(watflx>0._r8)then  
    heatflx = cpw*TKSoil1_vr(N3,N2,N1)*watflx
    !from (M6,M5,M4) to (N3,N2,N1)
  else
    heatflx = cpw*TKSoil1_vr(M6,M5,M4)*watflx
  endif

  HeatFlow2Soili_3D(N,M6,M5,M4) = heatflx

  QDrain_col(N2,N1)     = QDrain_col(N2,N1) - XN*watflx
  QDrain_cum_col(N2,N1) = QDrain_cum_col(N2,N1) - XN*watflx
  QDrainloss_vr(N3,N2,N1)=QDrainloss_vr(N3,N2,N1) - XN*watflx

  QDrainM_col(N2,N1)   = QDrainM_col(N2,N1) - XN*watflx
  HeatDrain_col(N2,N1) = HeatDrain_col(N2,N1)-XN*heatflx

  call PrintInfo('end '//subname)
  end subroutine VertBoundaryDrainM

!------------------------------------------------------------------------------------------
  subroutine UpdateSoilMoistTempM(I,J,M,NHW,NHE,NVN,NVS,dHGrnddTKsoi1)
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  ! Step-integrated net surface-heat derivative [MJ/K], at the old topsoil
  ! temperature. Its accepted correction is recorded after the implicit solve.
  real(r8), intent(in):: dHGrnddTKsoi1(JY,JX)
  character(len=*), parameter :: subname='UpdateSoilMoistTempM'
  integer :: NY,NX,L
  real(r8) :: tk1pres,tk1l,NetSurfaceHeatCorrection
  real(r8) :: ENGY1,VLTSoiPore,VHXX
  real(r8) :: TKXX,VLWMicPre,VLHeatCapacityPre,dHeat
  integer :: it   !current model step in the simulation year
  real(r8) :: dwat,dwat1,rootC,dflow_mic,dflow_mac,pscal,VLWatMicP1t_vr
  real(r8),parameter :: tau=100._r8
  it=(I-1)*24+J

! begin_execution
  call PrintInfo('beg '//subname)
  D11: DO NX=NHW,NHE
    D12: DO NY=NVN,NVS
      TWaterPlantRoot2SoilXM_col(NY,NX) = 0._r8
      dwat=0._r8;dwat1=0._r8;dflow_mic=0._r8;dflow_mac=0._r8
      D13: DO L = NUM_col(NY,NX), NL_col(NY,NX)

        IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN

          VLWMicPre              = VLWatMicP1_vr(L,NY,NX)
          dwat1=dwat1+VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE

          VLWatMicP1_vr(L,NY,NX) = VLWatMicP1_vr(L,NY,NX)+WatNetFlow2MicP_3DM_vr(L,NY,NX)+FWatExMacP2MicPiM_vr(L,NY,NX) &
            +TMLiceThawedMicP_vr(L,NY,NX)+FWatIrrigate2MicP1_vr(L,NY,NX)

          VLWatMicPX1_vr(L,NY,NX) = VLWatMicPX1_vr(L,NY,NX)+WatNetFlow2MicptX_3DM_vr(L,NY,NX)+FWatExMacP2MicPiM_vr(L,NY,NX) &
             +TMLiceThawedMicP_vr(L,NY,NX)+FWatIrrigate2MicP1_vr(L,NY,NX)

          VLiceMicP1_vr(L,NY,NX)  = VLiceMicP1_vr(L,NY,NX)-TMLiceThawedMicP_vr(L,NY,NX)/DENSICE
          
          VLWatMacP1_vr(L,NY,NX) = VLWatMacP1_vr(L,NY,NX)+WatNetFlow2Macpt_3DM_vr(L,NY,NX)-FWatExMacP2MicPiM_vr(L,NY,NX) &
            +TMLiceThawedMacP_vr(L,NY,NX)
          VLiceMacP1_vr(L,NY,NX) = VLiceMacP1_vr(L,NY,NX)-TMLiceThawedMacP_vr(L,NY,NX)/DENSICE

          dflow_mic=dflow_mic+WatNetFlow2MicP_3DM_vr(L,NY,NX)
          dflow_mac=dflow_mac+WatNetFlow2Macpt_3DM_vr(L,NY,NX)
          IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
            ! air-filled space
            VLairMicP_vr(L,NY,NX)  = VLMicP1_vr(L,NY,NX)-VLWatMicP1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)
            VLairMicP1_vr(L,NY,NX) = AZMAX1(VLairMicP_vr(L,NY,NX))
            VLairMacP_vr(L,NY,NX)  = VLMacroPores1_vr(L,NY,NX)-VLWatMacP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)
            VLairMacP1_vr(L,NY,NX) = AZMAX1(VLairMacP_vr(L,NY,NX))
            VLMacroPores1_vr(L,NY,NX)    = AZMAX1(VLMacroPores_vr(L,NY,NX)-FVOLAH*CCLAY_vr(L,NY,NX) &
              *(safe_adb(VLWatMicP1_vr(L,NY,NX),VLSoilMicP_vr(L,NY,NX))-WiltPoint_vr(L,NY,NX))*VGeomLayer_vr(L,NY,NX))
          ELSE
            VLairMicP_vr(L,NY,NX)  = 0.0_r8
            VLairMicP1_vr(L,NY,NX) = 0.0_r8
            VLairMacP_vr(L,NY,NX)  = 0.0_r8
            VLairMacP1_vr(L,NY,NX) = 0.0_r8
            VLMicP1_vr(L,NY,NX)    = VLWatMicP1_vr(L,NY,NX)+VLiceMicP1_vr(L,NY,NX)
            VLMacroPores1_vr(L,NY,NX)    = 0.0_r8
          ENDIF

          dwat=dwat+VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE
          !total soil volume with pores
          VLTSoiPore                        = VLSoilMicP_vr(L,NY,NX)+VLMacroPores1_vr(L,NY,NX)
          FracSoiPAsWat_vr(L,NY,NX)         = AZMAX1t((VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX))/VLTSoiPore)
          FracSoiPAsIce_vr(L,NY,NX)         = AZMAX1t((VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))/VLTSoiPore)
          FracAirFilledSoilPore_vr(L,NY,NX) = AZMAX1t((VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX))/VLTSoiPore)
          
          IF(VLMicP1_vr(L,NY,NX)+VLMacroPores1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            FracSoilAsAirt(L,NY,NX)=AZMAX1((VLairMicP1_vr(L,NY,NX)+ &
              VLairMacP1_vr(L,NY,NX))/(VLMicP1_vr(L,NY,NX)+VLMacroPores1_vr(L,NY,NX)))
          ELSE
            FracSoilAsAirt(L,NY,NX)=0.0_r8
          ENDIF
          IF(VLMacroPores1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            SoilFracAsMacP1_vr(L,NY,NX) = SoilFracAsMacP_vr(L,NY,NX)*VLMacroPores1_vr(L,NY,NX)/VLMacroPores_vr(L,NY,NX)
            HydroCondMacP1_vr(L,NY,NX)  = HydroCondMacP_vr(L,NY,NX)*(VLMacroPores1_vr(L,NY,NX)/VLMacroPores_vr(L,NY,NX))**2
          ELSE
            SoilFracAsMacP1_vr(L,NY,NX) = 0.0_r8
            HydroCondMacP1_vr(L,NY,NX)  = 0.0_r8
          ENDIF
          SoilFracAsMicP_vr(L,NY,NX) = 1.0_r8-SoilFracAsMacP1_vr(L,NY,NX)
          TKXX                       = TKSoil1_vr(L,NY,NX)
          VHXX                       = VHeatCapacity1_vr(L,NY,NX)
          ENGY1                      = VHXX*TKXX
          
          VLHeatCapacityA_vr(L,NY,NX) = VHeatCapSolidSoil_vr(L,NY,NX)+cpw*VLWatMicP1_vr(L,NY,NX)+cpi*VLiceMicP1_vr(L,NY,NX)
          if(plantOM4Heat .and. VLHeatCapacityA_vr(L,NY,NX)>0._r8)then
            rootC=sum(RootMycoMassElm_vr(ielmc,:,L,NY,NX))
            VLHeatCapacityA_vr(L,NY,NX)=VLHeatCapacityA_vr(L,NY,NX)+cpo*gOC_to_m3_OM(rootC)
          endif
          VLHeatCapacityB_vr(L,NY,NX) = cpw*VLWatMacP1_vr(L,NY,NX)+cpi*VLiceMacP1_vr(L,NY,NX)
          VHeatCapacity1_vr(L,NY,NX)  = VLHeatCapacityA_vr(L,NY,NX)+VLHeatCapacityB_vr(L,NY,NX)

         IF(VHeatCapacity1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN

            VLWatMicP1t_vr  = VLWatMicP1_vr(L,NY,NX)+dWaterPlantRoot2SoilPrev_vr(L,NY,NX)

            if(VLWatMicP1t_vr<0._r8 .and. VLWatMicP1_vr(L,NY,NX)>0._r8)then
              pscal=-VLWatMicP1_vr(L,NY,NX)/dWaterPlantRoot2SoilPrev_vr(L,NY,NX)
              TWaterPlantRoot2SoilPrev_vr(L,NY,NX) = TWaterPlantRoot2SoilPrev_vr(L,NY,NX)-(1._r8-pscal)*dWaterPlantRoot2SoilPrev_vr(L,NY,NX)
              VLWatMicP1_vr(L,NY,NX)               = VLWatMicP1_vr(L,NY,NX)+dWaterPlantRoot2SoilPrev_vr(L,NY,NX)*pscal
            else
              pscal=1._r8
              VLWatMicP1_vr(L,NY,NX)=VLWatMicP1t_vr
            endif

            VLWatMicPX1_vr(L,NY,NX) = VLWatMicPX1_vr(L,NY,NX)+dWaterPlantRoot2SoilPrev_vr(L,NY,NX)*pscal
            VLWatMicPX1_vr(L,NY,NX) = AMIN1(VLWatMicP1_vr(L,NY,NX),VLWatMicPX1_vr(L,NY,NX))

            TWaterPlantRoot2SoilXM_col(NY,NX) = TWaterPlantRoot2SoilXM_col(NY,NX)+dWaterPlantRoot2SoilPrev_vr(L,NY,NX)*pscal
            TWaterPlantRoot2SoilX_col(NY,NX) = TWaterPlantRoot2SoilX_col(NY,NX)+dWaterPlantRoot2SoilPrev_vr(L,NY,NX)*pscal

            if(dWaterPlantRoot2SoilPrev_vr(L,NY,NX)>0._r8)then
              ! The previous-hour transfer is applied here. Tally its accepted
              ! donor enthalpy in both directions for the delayed-transfer ledger.
              THeatPlantRoot2SoilPrev_vr(L,NY,NX)=THeatPlantRoot2SoilPrev_vr(L,NY,NX) &
                +dHeatPlantRoot2SoilPrev_vr(L,NY,NX)
              VHeatCapacity1_vr(L,NY,NX) = VHeatCapacity1_vr(L,NY,NX)+cpw*dWaterPlantRoot2SoilPrev_vr(L,NY,NX)
              if(L.EQ.NUM_col(NY,NX))then
                TKSoil1_vr(L,NY,NX)  = (ENGY1+THeatFlow2Soil_3DM_vr(L,NY,NX)+HeatIrrigation1_vr(L,NY,NX)&
                  +TLPhaseChangeHeat2Soi1_vr(L,NY,NX)+dHeatPlantRoot2SoilPrev_vr(L,NY,NX)-dHGrnddTKsoi1(NY,NX)*TKSoil1_vr(L,NY,NX)) &
                  /(VHeatCapacity1_vr(L,NY,NX)-dHGrnddTKsoi1(NY,NX))                 
              else
                TKSoil1_vr(L,NY,NX)        = (ENGY1+THeatFlow2Soil_3DM_vr(L,NY,NX)+HeatIrrigation1_vr(L,NY,NX)&
                  +TLPhaseChangeHeat2Soi1_vr(L,NY,NX)+dHeatPlantRoot2SoilPrev_vr(L,NY,NX))/VHeatCapacity1_vr(L,NY,NX)   
              endif
            else
              if(L.EQ.NUM_col(NY,NX))then
                TKSoil1_vr(L,NY,NX)  = (ENGY1+THeatFlow2Soil_3DM_vr(L,NY,NX)+HeatIrrigation1_vr(L,NY,NX)&
                  +TLPhaseChangeHeat2Soi1_vr(L,NY,NX)-dHGrnddTKsoi1(NY,NX)*TKSoil1_vr(L,NY,NX)) &
                  /(VHeatCapacity1_vr(L,NY,NX)-dHGrnddTKsoi1(NY,NX))                              
              else
                TKSoil1_vr(L,NY,NX)  = (ENGY1+THeatFlow2Soil_3DM_vr(L,NY,NX)+HeatIrrigation1_vr(L,NY,NX)&
                  +TLPhaseChangeHeat2Soi1_vr(L,NY,NX))/VHeatCapacity1_vr(L,NY,NX)                              
              endif  
              VHeatCapacity1_vr(L,NY,NX)          = VHeatCapacity1_vr(L,NY,NX)+cpw*dWaterPlantRoot2SoilPrev_vr(L,NY,NX)*pscal
              THeatPlantRoot2SoilPrev_vr(L,NY,NX) = THeatPlantRoot2SoilPrev_vr(L,NY,NX)+ &
                cpw*TKSoil1_vr(L,NY,NX)*dWaterPlantRoot2SoilPrev_vr(L,NY,NX)*pscal
            endif

            ! The stabilizing derivative represents real accepted boundary heat.
            ! Record it before any separately accounted prescribed warming.
            if(L.eq.NUM_col(NY,NX))then
              call CorrectSoilSurfaceHeatM(NY,NX,TKXX,TKSoil1_vr(L,NY,NX), &
                dHGrnddTKsoi1(NY,NX),NetSurfaceHeatCorrection)
              THeatFlow2Soil_3DM_vr(L,NY,NX)=THeatFlow2Soil_3DM_vr(L,NY,NX)+NetSurfaceHeatCorrection
            endif

            if(TKSoil1_vr(L,NY,NX)>400._r8.or.TKSoil1_vr(L,NY,NX)<100._r8)then
              write(*,*)'======'
              write(*,*)'VLWatMicP1_vr(L,NY,NX)=',VLWMicPre,VLWatMicP1_vr(L,NY,NX),L,VLairMicP_vr(L,NY,NX),&
                VLiceMicP1_vr(L,NY,NX)
              write(*,*)WatNetFlow2MicP_3DM_vr(L,NY,NX),FWatExMacP2MicPiM_vr(L,NY,NX), &
                TMLiceThawedMicP_vr(L,NY,NX),FWatIrrigate2MicP1_vr(L,NY,NX)
              write(*,*)'M, L=',M,L,NY,NX,NUM_col(NY,NX),VGeomLayer_vr(L,NY,NX),ZEROS2(NY,NX),VLWatMicPX1_vr(L,NY,NX)
              write(*,*)'SoilBulkDensity_vr(L,NY,NX)=',SoilBulkDensity_vr(L,NY,NX),TKS_vr(L,NY,NX),ZEROS(NY,NX)
              write(*,*)'VHeatCapacity1_vr(L,NY,NX),TKSoil1_vr(L,NY,NX),TKXX',L,VHeatCapacity1_vr(L,NY,NX),TKSoil1_vr(L,NY,NX),TKXX
              write(*,*)VLMicP1_vr(L,NY,NX),VLMacroPores1_vr(L,NY,NX)
              if(TKSoil1_vr(L,NY,NX)>1.e3_r8.or.TKSoil1_vr(L,NY,NX)<0._r8)call endrun(trim(mod_filename)//' at line',__LINE__)
            endif

            ! do soil warming
            if(do_warming)then
              if(is_warming_layerL(L,NY,NX))then     
                tk1l                   = TKSoil1_vr(L,NY,NX)
                dHeat                  = (TKS_ref_vr(it,L,NY,NX)-tk1l)*VHeatCapacity1_vr(L,NY,NX)
                HeatSource_vr(L,NY,NX) = HeatSource_vr(L,NY,NX)+dHeat
                TKSoil1_vr(L,NY,NX)    = tk1l+dHeat/(VHeatCapacity1_vr(L,NY,NX))
              endif
            endif
            
            if(TKSoil1_vr(L,NY,NX)>400.0_r8)then
              write(*,*)'TKSoil1_vr(L,NY,NX)',I+J/24.,M,L,TKSoil1_vr(L,NY,NX),tk1l,VHXX,VHeatCapacity1_vr(L,NY,NX),TKS_vr(L,NY,NX)
              write(*,*)ENGY1/VHeatCapacity1_vr(L,NY,NX),THeatFlow2Soil_3DM_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX),&
               TLPhaseChangeHeat2Soi1_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX),HeatIrrigation1_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX)
              write(*,*)VLWatMicPM_vr(M,L,NY,NX),VLWatMicPM_vr(M+1,L,NY,NX),THeatFlow2Soil_3DM_vr(L,NY,NX)/VHeatCapacity1_vr(L,NY,NX)
              call endrun('crazy temperature in '//trim(mod_filename)//' at line',__LINE__)
            endif

          ELSEIF(L.EQ.1)THEN
            TKSoil1_vr(L,NY,NX)=TairK_col(NY,NX)
          ELSE
            TKSoil1_vr(L,NY,NX)=TKSoil1_vr(L-1,NY,NX)
          ENDIF
        ENDIF
      ENDDO D13

    ENDDO D12
  ENDDO D11
  call PrintInfo('end '//subname)
  end subroutine UpdateSoilMoistTempM

!------------------------------------------------------------------------------------------

  subroutine UpdateStateFluxAtM(I,J,M,NHW,NHE,NVN,NVS)
  !
  !Description
  ! Early exit of the watsub solver
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  integer :: NY,NX
  integer :: L,NUX,LL,Ls
  real(r8):: VHCPX,HFLX,dCPX,dwat

! begin_execution
  D9795: DO NX=NHW,NHE
    D9790: DO NY=NVN,NVS
      dwat=0._r8
      !
      ! SOIL LAYER WATER, ICE AND TEMPERATURE
      !
      D9785: DO L=NUM_col(NY,NX),NL_col(NY,NX)
        IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN

          !record intermediate variables for bgc calculation
          VLWatMicPM_vr(M+1,L,NY,NX) = VLWatMicP1_vr(L,NY,NX)
          VLWatMacPM_vr(M+1,L,NY,NX) = VLWatMacP1_vr(L,NY,NX)
          VLsoiAirPM_vr(M+1,L,NY,NX) = VLairMicP1_vr(L,NY,NX)+VLairMacP1_vr(L,NY,NX)+THETPI*(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))

          !change in soil air volume
          ReductVLsoiAirPM_vr(M,L,NY,NX)     = VLsoiAirPM_vr(M,L,NY,NX)-VLsoiAirPM_vr(M+1,L,NY,NX)
          FracAirFilledSoilPoreM_vr(M+1,L,NY,NX) = FracAirFilledSoilPore_vr(L,NY,NX)
 !
        ELSE
          !layer L disappears
          VLWatMicPM_vr(M+1,L,NY,NX)         = 0.0_r8
          VLWatMacPM_vr(M+1,L,NY,NX)         = 0.0_r8
          VLsoiAirPM_vr(M+1,L,NY,NX)         = 0.0_r8
          ReductVLsoiAirPM_vr(M,L,NY,NX)     = VLsoiAirPM_vr(M,L,NY,NX)
          FracAirFilledSoilPoreM_vr(M+1,L,NY,NX) = 0.0_r8
        ENDIF
        dwat=dwat+VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)+(VLiceMicP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))*DENSICE
      ENDDO D9785

      !
      !       RESET SURFACE LAYER NUMBER AND TRANSFER ALL WATER TO SOIL SURFACE LAYER
      !       IF LAKE SURFACE LAYER IS LOST TO EVAPORATION
      !
      !       NUM=new surface layer number after CO2CompenPoint_nodeete lake evaporation
      !       LakeSurfFlowMicP,LakeSurfFlowMacP,LakeSurfHeatFlux=lake surface water flux, heat flux if lake surface disappears
      ! NUM layer is pure water, and its heat capacity and thereby water is too little
      IF(SoilBulkDensity_vr(NUM_col(NY,NX),NY,NX).LE.ZERO .AND. VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX).LE.VHCPNX_col(NY,NX))THEN
        !the soil/water profile moves down        
        NUX=NUM_col(NY,NX)                 
        D9970: DO  LL=NUX+1,NL_col(NY,NX)
          !the first meaningful layer
          IF(VHeatCapacity1_vr(LL,NY,NX).GT.0._r8)THEN
            NUM_col(NY,NX) = LL
            !heat capacity of a too thin layer
            dCPX           = cpw*(VLWatMicP1_vr(NUX,NY,NX)+VLWatMacP1_vr(NUX,NY,NX)) &
              +cpi*(VLiceMicP1_vr(NUX,NY,NX)+VLiceMacP1_vr(NUX,NY,NX))
            !enthalpy
            HFLX           = dCPX*TKSoil1_vr(NUX,NY,NX)

            !NUX layer is gone
            WaterFlowSoiMicP_3D(3,NUM_col(NY,NX),NY,NX)  = WaterFlowSoiMicP_3D(3,NUM_col(NY,NX),NY,NX) +VLWatMicP1_vr(NUX,NY,NX)
            WaterFlowSoiMicPX_3D(3,NUM_col(NY,NX),NY,NX) = WaterFlowSoiMicPX_3D(3,NUM_col(NY,NX),NY,NX)+VLWatMicP1_vr(NUX,NY,NX)
            WaterFlowSoiMacP_3D(3,NUM_col(NY,NX),NY,NX)  = WaterFlowSoiMacP_3D(3,NUM_col(NY,NX),NY,NX) +VLWatMacP1_vr(NUX,NY,NX)
            HeatFlow2Soil_3D(3,NUM_col(NY,NX),NY,NX)     = HeatFlow2Soil_3D(3,NUM_col(NY,NX),NY,NX) + HFLX

            WatFlow2TopSoiMicP_col(NY,NX)                = WaterFlowSoiMicP_3D(3,NUM_col(NY,NX),NY,NX)
            WatFlow2TopSoiMicPX_col(NY,NX)               = WaterFlowSoiMicPX_3D(3,NUM_col(NY,NX),NY,NX)
            WatFlow2TopSoiMacP_col(NY,NX)                = WaterFlowSoiMacP_3D(3,NUM_col(NY,NX),NY,NX)
            HeatFlow2TopSoi_col(NY,NX)                   = HeatFlow2Soil_3D(3,NUM_col(NY,NX),NY,NX)

            VHCPX                                   = VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX)
            QIceInflx_vr(NUM_col(NY,NX),NY,NX)      = QIceInflx_vr(NUM_col(NY,NX),NY,NX) &
              +VLiceMicP1_vr(NUX,NY,NX)+VLiceMacP1_vr(NUX,NY,NX)
            VLWatMicP1_vr(NUM_col(NY,NX),NY,NX)     = VLWatMicP1_vr(NUM_col(NY,NX),NY,NX)+ VLWatMicP1_vr(NUX,NY,NX)
            VLWatMacP1_vr(NUM_col(NY,NX),NY,NX)     = VLWatMacP1_vr(NUM_col(NY,NX),NY,NX)+ VLWatMacP1_vr(NUX,NY,NX)
            VLiceMicP1_vr(NUM_col(NY,NX),NY,NX)     = VLiceMicP1_vr(NUM_col(NY,NX),NY,NX) &
              +VLiceMicP1_vr(NUX,NY,NX)+VLiceMacP1_vr(NUX,NY,NX)
            VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX) = VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX)+dCPX
            !update the temperature by merging the two layers
            TKSoil1_vr(NUM_col(NY,NX),NY,NX)        = (VHCPX*TKSoil1_vr(NUM_col(NY,NX),NY,NX)+HFLX)/VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX)
            exit
          ENDIF
        ENDDO D9970
        !set the layer content to zero
        VLWatMicP1_vr(NUX,NY,NX) = 0._r8
        VLWatMacP1_vr(NUX,NY,NX) = 0._r8
        VLiceMicP1_vr(NUX,NY,NX) = 0._r8
        VLiceMacP1_vr(NUX,NY,NX) = 0._r8
      ENDIF
    ENDDO D9790
  ENDDO D9795

  end subroutine UpdateStateFluxAtM

!------------------------------------------------------------------------------------------

  subroutine UpdateFluxAtExit(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX
! begin_execution
  if(fixWaterLevel)return

  D9695: DO NX=NHW,NHE
    D9690: DO NY=NVN,NVS
      !there is not change in surface layer      
      IF(NUM_col(NY,NX).EQ.NU_col(NY,NX))THEN
        !vertical flux, incoming
        !top layer is the same
        LakeSurfFlowMicP_col(NY,NX)  = WaterFlowSoiMicP_3D(3,N6X(NY,NX),NY,NX)
        LakeSurfFlowMicPX_col(NY,NX) = WaterFlowSoiMicPX_3D(3,N6X(NY,NX),NY,NX)
        LakeSurfFlowMacP_col(NY,NX)  = WaterFlowSoiMacP_3D(3,N6X(NY,NX),NY,NX)
        LakeSurfHeatFlux_col(NY,NX)  = HeatFlow2Soil_3D(3,N6X(NY,NX),NY,NX)

        !the top soil/water layer has changed, e.g. for ponding water
      ELSE
        LakeSurfFlowMicP_col(NY,NX)  = WatFlow2TopSoiMicP_col(NY,NX)
        LakeSurfFlowMicPX_col(NY,NX) = WatFlow2TopSoiMicPX_col(NY,NX)
        LakeSurfFlowMacP_col(NY,NX)  = WatFlow2TopSoiMacP_col(NY,NX)
        LakeSurfHeatFlux_col(NY,NX)  = HeatFlow2TopSoi_col(NY,NX)
      ENDIF
!      call writeSurfDiagnosis(I,J,NY,NX)
    ENDDO D9690
  ENDDO D9695
  
  end subroutine UpdateFluxAtExit

!--------------------------------------------------------------------------
  pure function get_Surf2MacPDepz(L,NY,NX)result(Surf2MacPDepz)
  !
  !Description:
  !get the cumulative depth of macropore water 
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8) :: Surf2MacPDepz     !depth from soil surface to macropore in layer L

  !has macropore water
  IF(VLMacroPores1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    Surf2MacPDepz= CumDepz2LayBottom_vr(L,NY,NX)-(VLWatMacP1_vr(L,NY,NX)+VLiceMacP1_vr(L,NY,NX))/VLMacroPores1_vr(L,NY,NX)*DLYR_3D(3,L,NY,NX)
  ELSE
    !has no macropore water
    Surf2MacPDepz=CumDepz2LayBottom_vr(L,NY,NX)
  ENDIF

  end function get_Surf2MacPDepz

!--------------------------------------------------------------------------
  subroutine Config4TileDrainage(I,J,N,L,NY,NX,XN,Surf2MacPDepz,DoMicPDischarg2Tile,DoMacPDischarg2Tile)

  implicit none
  integer , intent(in) :: I,J,N,L,NY,NX
  real(r8), intent(in) :: XN
  real(r8),intent(in) :: Surf2MacPDepz                 !depth of macropore water
  logical, intent(out) :: DoMicPDischarg2Tile
  logical, intent(out) :: DoMacPDischarg2Tile
  character(len=*), parameter :: subname='Config4TileDrainage'
  integer :: LL
  real(r8) :: TileWaterTableX,PSISWD,PSISWDC
  real(r8) :: TileEquilPotential ![MPa]

  call PrintInfo('beg '//subname)
  !
  !     IDWaterTable=water table flag
  !     DPTH,TileWaterTable_col=depth to layer midpoint, artificial water table
  !     PSISM1,PSISE=soil,air entry matric potential
  !     TileWaterTableX=equilibrium water potential with artificial water table
  !     DoMicPDischarg2Tile=micropore discharge flag to artificial water table
  ! 
  DoMicPDischarg2Tile=.false.
  DoMacPDischarg2Tile=.false.
  if(N.eq.iVerticalDirection)return

  IF(IDWaterTable_col(NY,NX).GE.3 .AND. &
    SoilDepthMidLay_vr(L,NY,NX).LT.TileWaterTable_col(NY,NX))THEN !above tile

    PSISWDC = XN*mGravAccelerat*SLOPE_col(N,NY,NX)*(1.0_r8-WaterTable2GroundSlope_col(NY,NX))
    PSISWD  = PSISWDC*DLYR_3D(N,L,NY,NX)*0.5_r8
    !calc equilibrium tile water potential using hydrostatic equilbrium
    TileEquilPotential=mGravAccelerat*(SoilDepthMidLay_vr(L,NY,NX)-TileWaterTable_col(NY,NX))-PSISWD

    !equilibrium tile water table below soil middle layer depth
    IF(PSISM1_vr(L,NY,NX).GT.TileEquilPotential)THEN
      DoMicPDischarg2Tile=.true.

      IF(L.LT.NL_col(NY,NX))THEN
        D9568: DO  LL=L+1,NL_col(NY,NX)
          PSISWD = PSISWDC*DLYR_3D(N,LL,NY,NX)*0.5_r8
          TileWaterTableX=TileWaterTable_col(NY,NX)+(PSI_at_Sautration_vr(LL,NY,NX)+PSISWD)/mGravAccelerat
          IF(SoilDepthMidLay_vr(LL,NY,NX).LT.TileWaterTableX)THEN !above tile
            IF((PSISM1_vr(LL,NY,NX).LE.mGravAccelerat*(SoilDepthMidLay_vr(LL,NY,NX)-TileWaterTableX) &
             .AND. LL.NE.NL_col(NY,NX)) & !not bottom layer
             .OR. SoilDepthMidLay_vr(LL,NY,NX).GT.ActiveLayDepZ_col(NY,NX))THEN !layer below active layer
              DoMicPDischarg2Tile=.false.
            ENDIF
          ENDIF
        ENDDO D9568
      ENDIF
    ENDIF
  ENDIF
  !
  !     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO TILE DRAIN
  !
  !     VLMacP1,VLWatMacP1,VLiceMacP1=macropore volume,water,ice content
  !     CDPTH=depth to layer bottom
  !     DLYR=layer thickness
  !     DoMacPDischarg2Tile=macropore discharge flag to artificial water table
  !
  if(LMacPoreOff)return
  DoMacPDischarg2Tile=IDWaterTable_col(NY,NX).GE.3 .AND.              & !title drainage is on
    Surf2MacPDepz.LT.TileWaterTable_col(NY,NX)                        & !layer above tile
    .AND. VLWatMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX)                       !has water in macropore   

  IF(L.LT.NL_col(NY,NX) .and. DoMacPDischarg2Tile)THEN
    D9569: DO  LL=L+1,NL_col(NY,NX)
      IF(SoilDepthMidLay_vr(LL,NY,NX).LT.TileWaterTable_col(NY,NX))THEN
        ! the layer is above tile drain water table
          IF(VLMacroPores1_vr(LL,NY,NX).LE.ZEROS(NY,NX))THEN
          ! no macropore discharge to tile
            DoMacPDischarg2Tile=.false.
        ENDIF
      ENDIF
    ENDDO D9569
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine Config4TileDrainage
!------------------------------------------------------------------

  subroutine Config4XWaterTableDrain(I,J,N,L,NY,NX,XN,Surf2MacPDepz,DoMicPDischarg2ExtWTBL,DoMacPDischarg2ExtWTBL)
  !
  !configure drainage to external water table
  implicit none
  integer , intent(in) :: I,J,N,L,NY,NX
  real(r8), intent(in) :: XN
  real(r8), intent(in):: Surf2MacPDepz                         !depth to layer macropore water
  logical , intent(out):: DoMicPDischarg2ExtWTBL       !true to do water discharge from micropore to external water table 
  logical , intent(out):: DoMacPDischarg2ExtWTBL       !true to do water discharge from macropore to external water table

  character(len=*), parameter :: subname='Config4XWaterTableDrain' 
  real(r8) :: ExtWaterTableEquil
  real(r8) :: PSISWDC,PSISWD !pressure drop due to elevation change (following the slope) [MPa] 
  real(r8) :: WaterTableXEquilPotential ![MPa]
  integer :: LL

  DoMicPDischarg2ExtWTBL = .false.
  DoMacPDischarg2ExtWTBL = .false.
  if(N.eq.iVerticalDirection)return  
  ! ECOSIM_ZM: legacy ecosys watsub.f (5283--5310) activates micropore
  ! discharge only above its conductivity-derived wetness threshold, and
  ! requires lower layers to meet that threshold too. Compute the equivalent
  ! threshold from EcoSIM's OWN hydraulic table; do not change the retention
  ! curve or its air-entry parameter used by other processes.

  if (use_ecosys_drainage) then
    if (IDWaterTable_col(NY,NX).ne.0 .and. SoilDepthMidLay_vr(L,NY,NX).lt.ExtWaterTable_col(NY,NX)) then
      if (PSISoilMatricPtmp_vr(L,NY,NX).gt.ZmOutletWetnessPotential(L,NY,NX)) then
        DoMicPDischarg2ExtWTBL = .true.
        if (L.lt.NL_col(NY,NX)) then
          do LL=L+1,NL_col(NY,NX)
            if (SoilDepthMidLay_vr(LL,NY,NX).lt.ExtWaterTable_col(NY,NX)) then
              if (PSISoilMatricPtmp_vr(LL,NY,NX).le.ZmOutletWetnessPotential(LL,NY,NX) .or. &
                  SoilDepthMidLay_vr(LL,NY,NX).gt.ActiveLayDepZ_col(NY,NX)) then
                DoMicPDischarg2ExtWTBL = .false.
                exit
              endif
            endif
          enddo
        endif
      endif
    endif
    ! Macropore activation retains the upstream/legacy water-table test.
    DoMacPDischarg2ExtWTBL=IDWaterTable_col(NY,NX).ne.0 .and. Surf2MacPDepz.lt.ExtWaterTable_col(NY,NX) &
      .and. VLWatMacP1_vr(L,NY,NX).gt.ZEROS2(NY,NX)
    return
  endif

  !     IDENTIFY CONDITIONS FOR MICROPRE DISCHARGE TO WATER TABLE
  !
  !     IDWaterTable=water table flag
  !     DPTH,ExtWaterTable=depth to layer midpoint,natural water table
  !     PSISM1_vr(<0),PSI_at_Sautration_vr(<0)=matric,air entry water potential [MPa], from water height to MPa, 1.e-6*9.8*1.e3*H
  !     ExtWaterTableEquil=equilibrium water potential with natural water table
  !     ActiveLayDepZ_col=active layer depth
  !     DoMicPDischarg2ExtWTBL=micropore discharge flag to natural water table
  !    total water potential psi+rho*grav*h
  ! note that depth is counting downward, so that smaller value designates shallower depth.
  ! layer 
  IF(IDWaterTable_col(NY,NX).NE.0 .AND. SoilDepthMidLay_vr(L,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN !layer above external water table        

    PSISWDC = XN*mGravAccelerat*SLOPE_col(N,NY,NX)*(1.0_r8-WaterTable2GroundSlope_col(NY,NX))
    PSISWD  = PSISWDC*DLYR_3D(N,L,NY,NX)*0.5_r8

    !assuming zero suction pressure at water table z_w
    !psi_ei-g*z_i=-g*z_w, psi_ei=g(z_i-z_w), or, z_i=z_w+psi_{ei}/g

    WaterTableXEquilPotential= mGravAccelerat*(SoilDepthMidLay_vr(L,NY,NX)-ExtWaterTable_col(NY,NX))-PSISWD

    if(PSISM1_vr(L,NY,NX).GT.WaterTableXEquilPotential)then
      DoMicPDischarg2ExtWTBL=.true.
      D9565: DO LL=L+1,NL_col(NY,NX)
        PSISWD = PSISWDC*DLYR_3D(N,LL,NY,NX)*0.5_r8      
        !water level difference        
        ExtWaterTableEquil=ExtWaterTable_col(NY,NX)+(PSI_at_Sautration_vr(LL,NY,NX)+PSISWD)/mGravAccelerat  

        !soil layer mid-depth above water table
        IF(SoilDepthMidLay_vr(LL,NY,NX).LT.ExtWaterTableEquil)THEN
          !if a lower layer LL does meet with the drainage requirement, turn off the discharge and ensure the
          !water move internally instead
          IF((PSISM1_vr(LL,NY,NX).LE.mGravAccelerat*(SoilDepthMidLay_vr(LL,NY,NX)-ExtWaterTableEquil) &   !pressure head is too low
            .AND. LL.NE.NL_col(NY,NX)) .OR. SoilDepthMidLay_vr(LL,NY,NX).GT.ActiveLayDepZ_col(NY,NX))THEN !layer below active layer
            DoMicPDischarg2ExtWTBL=.false.
          ENDIF
        ENDIF
      ENDDO D9565
    ENDIF
  ENDIF

  if(LMacPoreOff)return
  !
  !     IDENTIFY CONDITIONS FOR MACROPORE DISCHARGE TO WATER TABLE
  !
  !     VLMacP1,VLWatMacP1,VLiceMacP1=macropore volume,water,ice content
  !     CDPTH=depth to layer bottom
  !     DLYR=layer thickness
  !     DoMacPDischarg2ExtWTBL=macropore discharge flag to natural water table
  !  
  DoMacPDischarg2ExtWTBL=IDWaterTable_col(NY,NX).NE.0 .AND. Surf2MacPDepz.LT.ExtWaterTable_col(NY,NX) &
    .AND. VLWatMacP1_vr(L,NY,NX).GT.ZEROS2(NY,NX)
  
  IF(L.LT.NL_col(NY,NX) .and. DoMacPDischarg2ExtWTBL)THEN
    DO  LL=L+1,NL_col(NY,NX)
      IF(SoilDepthMidLay_vr(LL,NY,NX).LT.ExtWaterTable_col(NY,NX))THEN
        ! the layer is above tile drain water table
          IF(VLMacroPores1_vr(LL,NY,NX).LE.ZEROS(NY,NX))THEN
          ! no macropore discharge to tile
            DoMacPDischarg2ExtWTBL=.false.
        ENDIF
      ENDIF
    ENDDO 
  ENDIF

  end subroutine Config4XWaterTableDrain
!------------------------------------------------------------------------------------------
  function ZmOutletWetnessPotential(L,NY,NX) result(psi_threshold)
    ! ECOSIM_ZM: reproduce legacy hour1.f's FSCNV=0.1 crossing of the
    ! vertical conductivity lookup. K=1 is saturated, increasing K is drier.
    ! Binary search is equivalent to scanning this monotonic lookup and avoids
    ! rebuilding the hydraulic table inside every solver substep.
    integer, intent(in) :: L,NY,NX
    integer :: lo,hi,mid,nk
    real(r8) :: psi_threshold,cutoff,theta_threshold

    nk=size(HydroCond_3D,2)
    cutoff=0.1_r8*SatHydroCondVert_vr(L,NY,NX)
    psi_threshold=PSI_at_Sautration_vr(L,NY,NX)
    if (cutoff.le.0._r8 .or. HydroCond_3D(3,nk,L,NY,NX).ge.cutoff) return
    lo=2; hi=nk
    do while (lo.lt.hi)
      mid=(lo+hi)/2
      if (HydroCond_3D(3,mid,L,NY,NX).lt.cutoff) then
        hi=mid
      else
        lo=mid+1
      endif
    enddo
    theta_threshold=POROS_vr(L,NY,NX)*(1._r8-real(lo-2,r8)/real(nk,r8))
    call ComputePsiMCM(L,NY,NX,theta_threshold,psi_threshold)
  end function ZmOutletWetnessPotential

!------------------------------------------------------------------------------------------

  subroutine InitSoil3DModelIterateM(I,J,M,NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  character(len=*),parameter :: subname='InitSoil3DModelIterateM'
  integer :: NY,NX,L
  real(r8) :: THETWT,TScal4Aquadifsvity,THETWA
  real(r8) :: VLSoiPorAvail   !avaiable pore space for air + water, excluding ice
  real(r8) :: VLWatSoi,scalar
  real(r8) :: THETWH,Z3S

  call PrintInfo('beg '//subname)
  WaterFlow2Micptl_3D=0._r8  
  WaterFlow2Macptl_3D=0._r8        

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      D9885: DO L=NUM_col(NY,NX),NL_col(NY,NX)
    !
    !   GAS EXCHANGE COEFFICIENTS SOIL LAYERS
    !
    !   VOLA1,VLiceMicP1,VLWatMicP1=total,ice-,water-filled microporosity
    !   VLMacP1,VLiceMacP1,VLWatMacP1=total,ice-,water-filled macroporosity
    !   VLsoiAirPM=air-filled porosity
    !   TScal4Aquadifsvity=temperature effect on gas diffusivity
    !   DiffusivitySolutEff=rate constant for air-water gas exchange
    !   Z1S,Z2SW,Z2SD,Z3SX=parameters for soil air-water gas transfers
    !   XNPD=time step for gas transfer calculations
    !   TORT,TortMacPM_vr=tortuosity for aqueous diffn in micropores,macropres
    !
        VLWatSoi      = VLWatMicP1_vr(L,NY,NX)+VLWatMacP1_vr(L,NY,NX)
        VLSoiPorAvail = VLMicP1_vr(L,NY,NX)+VLMacroPores1_vr(L,NY,NX)-VLiceMicP1_vr(L,NY,NX)-VLiceMacP1_vr(L,NY,NX)

        !unsaturated soil
        IF(VLSoiPorAvail.GT.ZEROS2(NY,NX) .AND. VLsoiAirPM_vr(M,L,NY,NX).GT.ZEROS2(NY,NX))THEN
          THETWA                             = AZMAX1(AMIN1(1.0_r8,VLWatSoi/VLSoiPorAvail))
          TScal4Aquadifsvity                 = TEFAQUDIF(TKSoil1_vr(L,NY,NX))
          Z3S                                = FieldCapacity_vr(L,NY,NX)/POROS_vr(L,NY,NX)
          scalar                             = TScal4Aquadifsvity*XNPD
          DiffusivitySolutEffM_vr(M,L,NY,NX) = fDiffusivitySolutEff(scalar,THETWA,Z3S)
        ELSE
          DiffusivitySolutEffM_vr(M,L,NY,NX)=0.0_r8
        ENDIF

        !Pure water layer
        IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
          THETWT                  = safe_adb(VLWatMicPM_vr(M,L,NY,NX),VLSoilMicP_vr(L,NY,NX))
          TortMicPM_vr(M,L,NY,NX) = TortMicporew(THETWT)*(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
        !soil layer  
        ELSE
          TortMicPM_vr(M,L,NY,NX)=0.7_r8
        ENDIF
        !
        IF(VLMacroPores1_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          THETWH                  = VLWatMacPM_vr(M,L,NY,NX)/VLMacroPores1_vr(L,NY,NX)
          TortMacPM_vr(M,L,NY,NX) = TortMacporew(THETWH)*SoilFracAsMacP_vr(L,NY,NX)
        ELSE
          TortMacPM_vr(M,L,NY,NX)=0.0_r8
        ENDIF
      ENDDO D9885
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine InitSoil3DModelIterateM

!------------------------------------------------------------------------------------------
  subroutine MicropXGridDarcyFlow(N,N1,N2,N3,N4,N5,N6,THETA1,THETAL,KSatRedusByRainKinetEnergy,&
    PredDarcyFlowMax,WatDarcyFlowMicP,HeatByDarcyFlowMicP,PSISV1,PSISVL)     
  !
  !Description:
  !micropore flow across grids     
  implicit none
  integer, intent(in)  :: N                !flow direction
  integer, intent(in)  :: N1,N2,N3         !source grid
  integer, intent(in)  :: N4,N5,N6         !destination grid
  real(r8), intent(in) :: THETA1,THETAL    !micropore water concn in source,destination cells
  real(r8), intent(in) :: KSatRedusByRainKinetEnergy
  real(r8), intent(out):: PredDarcyFlowMax
  real(r8), intent(out):: WatDarcyFlowMicP
  real(r8), intent(out):: HeatByDarcyFlowMicP
  real(r8), intent(out):: PSISV1,PSISVL

  character(len=*), parameter :: subname='MicropXGridDarcyFlow'
  real(r8) :: AVE_CONDUCTANCE,HydCondSrc,HydCondDest
  real(r8) :: PSISTL,PSIST1,THETW1
  real(r8) :: PtWatDarcyFlux,PredDarcyFlow,PSISTMP,THETWL
  integer :: K1,KL

  call PrintInfo('beg '//subname)
   PredDarcyFlowMax    = 0._r8
   WatDarcyFlowMicP    = 0._r8
   HeatByDarcyFlowMicP = 0._r8
   PSISV1              = 0._r8
   PSISVL              = 0._r8

  IF(DLYR_3D(N,N6,N5,N4).LE.0._r8 .OR. DLYR_3D(N,N3,N2,N1).LE.0._r8 .OR. fixWaterLevel)return

  !     THETW1,THETWL=water concentrations in source,destination cells
  !air entry pressure is the pressure when large-pores started to be air-filled, when soil matric pressure is less
  !than air entry pressure, it is unsaturated
  IF(PSISoilMatricPtmp_vr(N3,N2,N1).GT.PSISoilAirEntry_vr(N3,N2,N1) .AND. &
    PSISoilMatricPtmp_vr(N6,N5,N4).GT.PSISoilAirEntry_vr(N6,N5,N4))THEN
    !Both grids are saturated
    THETW1 = THETA1
    THETWL = THETAL
    K1     = getMoistK(THETW1,POROS_vr(N3,N2,N1))
    KL     = getMoistK(THETWL,POROS_vr(N6,N5,N4))
    PSISM1_vr(N3,N2,N1) = PSISoilMatricPtmp_vr(N3,N2,N1)
    PSISM1_vr(N6,N5,N4) = PSISoilMatricPtmp_vr(N6,N5,N4)
    !
    !     GREEN-AMPT FLOW IF ONE LAYER IS SATURATED
    !     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POENTIAL)
    !
    !     GREEN-AMPT FLOW IF SOURCE CELL SATURATED
    !ThetaSat_vr=micropore soil water content
  ELSEIF(PSISoilMatricPtmp_vr(N3,N2,N1).GT.PSISoilAirEntry_vr(N3,N2,N1))THEN
    !Grid (N3,N2,N1) is saturated, while grid (N6,N5,N4) is not
    THETW1 = THETA1
    THETWL = AMAX1(SoilWatAirDry_vr(N6,N5,N4),AMIN1(POROS_vr(N6,N5,N4),safe_adb(VLWatMicPX1_vr(N6,N5,N4),VLSoilMicP_vr(N6,N5,N4))))
    K1     = getMoistK(THETW1,POROS_vr(N3,N2,N1))
    KL     = getMoistK(AMIN1(ThetaSat_vr(N6,N5,N4),THETWL),POROS_vr(N6,N5,N4))

    PSISM1_vr(N3,N2,N1) = PSISoilMatricPtmp_vr(N3,N2,N1)
    
    IF(VLSoilMicPMass_vr(N6,N5,N4).GT.ZEROS(N2,N1))THEN

      call ComputePsiMCM(N6,N5,N4,THETWL,PSISM1_vr(N6,N5,N4))

    ELSE
      !pond
      THETWL              = POROS_vr(N6,N5,N4)
      PSISM1_vr(N6,N5,N4) = PSI_at_Sautration_vr(N6,N5,N4)
    ENDIF
    !
    !     GREEN-AMPT FLOW IF ADJACENT CELL SATURATED
    
  ELSEIF(PSISoilMatricPtmp_vr(N6,N5,N4).GT.PSISoilAirEntry_vr(N6,N5,N4))THEN
    !Grid (N3,N2,N1) is not saturated, while grid (N6,N5,N4) is saturated
    THETW1 = AMAX1(SoilWatAirDry_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),safe_adb(VLWatMicPX1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))
    THETWL = THETAL
    K1     = getMoistK(AMIN1(ThetaSat_vr(N3,N2,N1),THETW1),POROS_vr(N3,N2,N1))
    KL     = getMoistK(THETWL,POROS_vr(N6,N5,N4))
    PSISM1_vr(N6,N5,N4) = PSISoilMatricPtmp_vr(N6,N5,N4)

    IF(VLSoilMicPMass_vr(N3,N2,N1).GT.ZEROS(N2,N1))THEN
      call ComputePsiMCM(N3,N2,N1,THETW1,PSISM1_vr(N3,N2,N1)) 
    ELSE
      THETW1              = POROS_vr(N3,N2,N1)
      PSISM1_vr(N3,N2,N1) = PSI_at_Sautration_vr(N3,N2,N1)
    ENDIF
    !
    !     RICHARDS FLOW IF NEITHER CELL IS SATURATED
    !     (CURRENT WATER POTENTIAL < AIR ENTRY WATER POTENTIAL)
    !
  ELSE
    !neither is saturated
    THETW1              = THETA1
    THETWL              = THETAL
    K1                  = getMoistK(THETW1,POROS_vr(N3,N2,N1)) 
    KL                  = getMoistK(THETWL,POROS_vr(N6,N5,N4))
    PSISM1_vr(N3,N2,N1) = PSISoilMatricPtmp_vr(N3,N2,N1)
    PSISM1_vr(N6,N5,N4) = PSISoilMatricPtmp_vr(N6,N5,N4)
  ENDIF
  !
  !  'HYDRAULIC CONUCTIVITY'
  !
  !     HydCondSrc,HydCondDest=hydraulic conductivity of source,destination layer
  !     HydroCond_3D=lateral(1,2),vertical(3) micropore hydraulic conductivity
  !
  HydCondSrc=HydroCond_3D(N,K1,N3,N2,N1)  
  IF(N3.EQ.NUM_col(N2,N1))THEN
    !surface soil
    HydCondSrc=HydCondSrc*KSatRedusByRainKinetEnergy
  ENDIF
  !
  HydCondDest                = HydroCond_3D(N,KL,N6,N5,N4)
  HydCondSoil_3D(N,N3,N2,N1) = HydCondSrc
  HydCondSoil_3D(N,N6,N5,N4) = HydCondDest
  !
  !     TOTAL SOIL WATER POTENTIAL = MATRIC, GRAVIMETRIC + OSMOTIC
  !
  !     PSISM1,PSISH,PSISO=soil matric,gravitational,osmotic potentials
  !src
  PSIST1 = PSISM1_vr(N3,N2,N1)+PSIGrav_vr(N3,N2,N1)+PSISoilOsmotic_vr(N3,N2,N1)
  PSISV1 = PSISM1_vr(N3,N2,N1)+PSISoilOsmotic_vr(N3,N2,N1)

  !dest
  PSISTL = PSISM1_vr(N6,N5,N4)+PSIGrav_vr(N6,N5,N4)+PSISoilOsmotic_vr(N6,N5,N4)
  PSISVL = PSISM1_vr(N6,N5,N4)+PSISoilOsmotic_vr(N6,N5,N4)
  !
  !     'HYDRAULIC CONDUCTIVITY FROM CURRENT WATER CONTENT'
  !     AND LOOKUP ARRAY GENERATED IN 'HOUR1'
  !
  !     HydCondSrc,HydCondDest=hydraulic conductivities in source,destination cells
  !     KSatRedusByRainKinetEnergy=reduction in soil surface Ksat from rainfall energy impact
  !     AVE_CONDUCTANCE=source-destination hydraulic conductance
  !     DLYR=layer thickness
  !
  IF(HydCondSrc.GT.ZERO .AND. HydCondDest.GT.ZERO)THEN
    AVE_CONDUCTANCE=2.0_r8*HydCondSrc*HydCondDest/(HydCondSrc*DLYR_3D(N,N6,N5,N4)+HydCondDest*DLYR_3D(N,N3,N2,N1))
  ELSE
    AVE_CONDUCTANCE=0.0_r8
  ENDIF
  !
  !     WATER FLUX FROM WATER POTENTIALS, HYDRAULIC CONDUCTIVITY
  !     CONSTRAINED BY WATER POTENTIAL GRADIENT, COUPLED WITH
  !     CONVECTIVE HEAT FLUX FROM WATER FLUX
  !
  !     dts_HeatWatTP=time step of flux calculations
  !     VOLW2,VOLP1=water,air contents of source,destination micropores
  !     HeatByWatFlowMicP=convective heat flux from micropore water flux
  !     VLairMicP=excess water+ice relative to porosity
  !
  PtWatDarcyFlux=AVE_CONDUCTANCE*(PSIST1-PSISTL)*AREA_3D(N,N3,N2,N1)*dts_HeatWatTP

  IF(PtWatDarcyFlux.GE.0.0_r8)THEN
    !water flow from (N3,N2,N1) to (N6,N5,N4)
    !
    IF(THETW1.GT.ThetaSat_vr(N3,N2,N1))THEN
      !water more than air-entry saturation
      PredDarcyFlow=PtWatDarcyFlux+AMIN1((THETW1-ThetaSat_vr(N3,N2,N1))*VLSoilMicP_vr(N3,N2,N1),&
        AZMAX1((ThetaSat_vr(N6,N5,N4)-THETWL)*VLSoilMicP_vr(N6,N5,N4)))*dts_wat
    ELSE
      PredDarcyFlow=PtWatDarcyFlux
    ENDIF

    WatDarcyFlowMicP = AZMAX1(AMIN1(PredDarcyFlow,VLWatMicP2_vr(N3,N2,N1)*dts_wat,VLairMicP1_vr(N6,N5,N4)*dts_wat))
    PredDarcyFlowMax = AZMAX1(AMIN1(PtWatDarcyFlux,VLWatMicP2_vr(N3,N2,N1)*dts_wat,VLairMicP1_vr(N6,N5,N4)*dts_wat))

  ELSE
    !water flow from (N6,N5,N4) to (N3,N2,N1)
    !
    IF(THETWL.GT.ThetaSat_vr(N6,N5,N4))THEN
      PredDarcyFlow=PtWatDarcyFlux+AMAX1((ThetaSat_vr(N6,N5,N4)-THETWL)*VLSoilMicP_vr(N6,N5,N4), &
        AZMIN1((THETW1-ThetaSat_vr(N3,N2,N1))*VLSoilMicP_vr(N3,N2,N1)))*dts_wat
    ELSE
      PredDarcyFlow=PtWatDarcyFlux
    ENDIF

    WatDarcyFlowMicP=AZMIN1(AMAX1(PredDarcyFlow,-VLWatMicP2_vr(N6,N5,N4)*dts_wat,-VLairMicP1_vr(N3,N2,N1)*dts_wat))
    PredDarcyFlowMax=AZMIN1(AMAX1(PtWatDarcyFlux,-VLWatMicP2_vr(N6,N5,N4)*dts_wat,-VLairMicP1_vr(N3,N2,N1)*dts_wat))

  ENDIF

  IF(N.EQ.iVerticalDirection .AND. VLairMicP_vr(N6,N5,N4).LT.0.0_r8)THEN
    !flow in the vertical direction, destination grid is saturated
    WatDarcyFlowMicP = WatDarcyFlowMicP+AZMIN1(AMAX1(-VLWatMicP2_vr(N6,N5,N4)*dts_wat,VLairMicP_vr(N6,N5,N4)))
    PredDarcyFlowMax = PredDarcyFlowMax+AZMIN1(AMAX1(-VLWatMicP2_vr(N6,N5,N4)*dts_wat,VLairMicP_vr(N6,N5,N4)))
  ENDIF

  !water flow from (N3,N2,N1) to (N6,N5,N4)
  IF(WatDarcyFlowMicP.GT.0.0_r8)THEN
    HeatByDarcyFlowMicP=cpw*TKSoil1_vr(N3,N2,N1)*WatDarcyFlowMicP*HeatAdv_scal
  ELSE
    !water flow from (N6,N5,N4) to (N3,N2,N1)
    HeatByDarcyFlowMicP=cpw*TKSoil1_vr(N6,N5,N4)*WatDarcyFlowMicP*HeatAdv_scal
  ENDIF

  VLWatMicP2_vr(N3,N2,N1) = VLWatMicP2_vr(N3,N2,N1)-WatDarcyFlowMicP
  VLWatMicP2_vr(N6,N5,N4) = VLWatMicP2_vr(N6,N5,N4)+WatDarcyFlowMicP
  call PrintInfo('end '//subname)
  end subroutine MicropXGridDarcyFlow
!------------------------------------------------------------------------------------------
  subroutine MacropXgridFLow(N,N1,N2,N3,N4,N5,N6,WaterMacpFlow,HeatByFlowMacP,LInvalidMacP)
  !
  ! Description:
  ! Compute internal macropore water and donor-temperature heat transfers.
  ! Limit upward overflow to liquid water available in the lower layer.
  implicit none
  integer, intent(in) :: N          !flow direction
  integer, intent(in) :: N1,N2,N3   !source grid
  integer, intent(in) :: N4,N5,N6   !destination grid
  real(r8), intent(out) :: HeatByFlowMacP   ! >0, flux into dst, 
  real(r8), intent(out) :: WaterMacpFlow
  logical, intent(inout) :: LInvalidMacP   !==1, dest grid invalid
  real(r8) :: FLWHX,PSISH1,PSISHL
  character(len=*), parameter :: subname='MacropXgridFLow'
  !     PSISH1,PSISHL=macropore total water potl in source,destination
  !     DLYR=layer thickness
  !     VLWatMacP1,VOLPH1=macropore water,air content

  HeatByFlowMacP = 0._r8
  WaterMacpFlow  = 0._r8

  if(fixWaterLevel)return

  IF(VLMacroPores1_vr(N3,N2,N1).GT.ZEROS2(N2,N1) .AND. VLMacroPores1_vr(N6,N5,N4).GT.ZEROS2(N5,N4) .AND. (.not.LInvalidMacP))THEN
    !macropore flow is on
    PSISH1=PSIGrav_vr(N3,N2,N1)+mGravAccelerat*DLYR_3D(3,N3,N2,N1)*(AMIN1(1.0_r8,&
      AZMAX1(VLWatMacP1_vr(N3,N2,N1)/VLMacroPores1_vr(N3,N2,N1)))-0.5_r8)
    PSISHL=PSIGrav_vr(N6,N5,N4)+mGravAccelerat*DLYR_3D(3,N6,N5,N4)*(AMIN1(1.0_r8,&
      AZMAX1(VLWatMacP1_vr(N6,N5,N4)/VLMacroPores1_vr(N6,N5,N4)))-0.5_r8)
    !
    !     MACROPORE FLOW IF GRAVITATIONAL GRADIENT IS POSITIVE
    !     AND MACROPORE POROSITY EXISTS IN ADJACENT CELL
    !
    !     FLWHX,WaterFlow2Macpt_3D=macropore water flux unltd,ltd by source water
    !     dts_HeatWatTP=time step of flux calculations
    !     VOLW2,VOLP1=water,air contents of source,destination micropores
    !     HeatByFlowMacP=convective heat flux from micropore water flux
    ! from (N3,N2,N1) to (N6,N5,N4)
    FLWHX=AVCNHL_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)*(PSISH1-PSISHL)*AREA_3D(N,N3,N2,N1)*dts_HeatWatTP
    IF(N.NE.iVerticalDirection)THEN
      !horizontal direction
      IF(PSISH1.GT.PSISHL)THEN
        !flow from source to dest grid
        WaterMacpFlow=AZMAX1(AMIN1(AMIN1(VLWatMacP1_vr(N3,N2,N1),VLairMacP1_vr(N6,N5,N4))*dts_wat,FLWHX))
      ELSEIF(PSISH1.LT.PSISHL)THEN
        !flow from dest grid to source
        WaterMacpFlow=AZMIN1(AMAX1(AMAX1(-VLWatMacP1_vr(N6,N5,N4),-VLairMacP1_vr(N3,N2,N1))*dts_wat,FLWHX))
      ELSE
        WaterMacpFlow=0.0_r8
      ENDIF
    ELSE
      !vertical direction
      WaterMacpFlow=AZMAX1(AMIN1(AMIN1(VLWatMacP1_vr(N3,N2,N1)*dts_wat &
        +WaterFlow2Macpt_3D(N,N3,N2,N1),VLairMacP1_vr(N6,N5,N4)*dts_wat),FLWHX))
    ENDIF

    IF(N.EQ.iVerticalDirection)THEN
      ! Negative air space includes ice expansion; it is not a liquid donor.
      ! Bound the accepted upward transfer before computing its heat and solute
      ! transport, as for the other Darcy fluxes. Retain excess ice in place.
      WaterMacpFlow=AMAX1(-AZMAX1(VLWatMacP1_vr(N6,N5,N4)*dts_wat), &
        WaterMacpFlow+AZMIN1(VLairMacP_vr(N6,N5,N4)))
    ENDIF

    if(abs(WaterMacpFlow).LT.tiny_wat)then
      WaterMacpFlow=0._r8
    endif
  ELSE
    WaterMacpFlow   = 0.0_r8        
  ENDIF

  IF(WaterMacpFlow.GT.0.0_r8)THEN
    HeatByFlowMacP=cpw*TKSoil1_vr(N3,N2,N1)*WaterMacpFlow
  ELSE
    HeatByFlowMacP=cpw*TKSoil1_vr(N6,N5,N4)*WaterMacpFlow
  ENDIF

  end subroutine MacropXgridFLow
!------------------------------------------------------------------------------------------

  subroutine WaterVaporXgridFlow(M,N,N1,N2,N3,N4,N5,N6,PSISV1,PSISVL,ConvectVapFlux,HeatByConvectVapFlux)
  !
  !Description:
  !water vapor flow across internal grids
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), intent(in) :: PSISV1,PSISVL
  real(r8), intent(out) :: ConvectVapFlux     !vapor flux from (N3,N2,N1) to (N6,N5,N4)
  real(r8), intent(out) :: HeatByConvectVapFlux

  character(len=*), parameter :: subname='WaterVaporXgridFlow'
  real(r8) :: TK11,TK12,VP1,VPL,VPY,CNV1,CNVL
  REAL(R8) :: ATCNVL,PotentialVaporFlux,MaxVaporFlux
  !     VAPOR PRESSURE AND DIFFUSIVITY IN EACH GRID CELL
  !
  
  !     THETPM,AirFillPore_Min=current, minimum air-filled porosity
  !     TK11,TK12=interim soil temperature in source,destination
  !     VP1,VPL=vapor concentration in source,destination
  !     PSISV1,PSISVL=matric+osmotic water potl in source,destination
  !     CNV1,CNVL=vapor conductivities of source, destination
  !     POROS,POROQ=porosity, tortuosity
  !     WVapDifusvitySoil_vr=vapor diffusivity
  !     ATCNVL=source,destination vapor conductance
  !     DLYR=soil layer depth
  !     PotentialVaporFlux,MaxVaporFlux=vapor flux unlimited,limited by vapor
  !     VPY=equilibrium vapor concentration
  !     dts_wat=time step for flux calculations
  !     ConvectVapFlux,ConvectiveHeatFlux=vapor flux and its convective heat flux
!
  call PrintInfo('beg '//subname)

  ConvectVapFlux       = 0._r8
  HeatByConvectVapFlux = 0._r8
  if(fixWaterLevel)return
  IF(FracAirFilledSoilPoreM_vr(M,N3,N2,N1).GT.AirFillPore_Min .AND. FracAirFilledSoilPoreM_vr(M,N6,N5,N4).GT.AirFillPore_Min)THEN
    TK11   = TKSoil1_vr(N3,N2,N1)
    TK12   = TKSoil1_vr(N6,N5,N4)
    VP1    = vapsat(TK11)*EXP(18.0_r8*PSISV1/(RGASC*TK11))
    VPL    = vapsat(TK12)*EXP(18.0_r8*PSISVL/(RGASC*TK12))
    CNV1   = WVapDifusvitySoil_vr(N3,N2,N1)*FracAirFilledSoilPoreM_vr(M,N3,N2,N1)**2*POROQ/POROS_vr(N3,N2,N1)
    CNVL   = WVapDifusvitySoil_vr(N6,N5,N4)*FracAirFilledSoilPoreM_vr(M,N6,N5,N4)**2*POROQ/POROS_vr(N6,N5,N4)
    ATCNVL = 2.0_r8*CNV1*CNVL/(CNV1*DLYR_3D(N,N6,N5,N4)+CNVL*DLYR_3D(N,N3,N2,N1))
    !
    !     VAPOR FLUX FROM VAPOR PRESSURE AND DIFFUSIVITY,
    !     AND CONVECTIVE HEAT FLUX FROM VAPOR FLUX
    !
    PotentialVaporFlux = ATCNVL*(VP1-VPL)*AREA_3D(N,N3,N2,N1)*dts_HeatWatTP
    VPY                = (VP1*VLsoiAirPM_vr(M,N3,N2,N1)+VPL*VLsoiAirPM_vr(M,N6,N5,N4))/(VLsoiAirPM_vr(M,N3,N2,N1)+VLsoiAirPM_vr(M,N6,N5,N4))
    MaxVaporFlux       = (VP1-VPY)*VLsoiAirPM_vr(M,N3,N2,N1)*dts_wat

    !from grid (N3,N2,N1) to (N6,N5,N4)
    IF(PotentialVaporFlux.GE.0.0_r8)THEN
      ConvectVapFlux       = AZMAX1d(AMIN1(PotentialVaporFlux,MaxVaporFlux,(VLWatMacP1_vr(N3,N2,N1)+VLWatMicP1_vr(N3,N2,N1))*dts_wat),tiny_wat)
      HeatByConvectVapFlux = (cpw*TKSoil1_vr(N3,N2,N1)+EvapLHTC)*ConvectVapFlux
    !from grid (N6,N5,N4) to (N3,N2,N1)
    ELSE
      ConvectVapFlux       = AZMIN1d(AMAX1(PotentialVaporFlux,MaxVaporFlux,-(VLWatMacP1_vr(N6,N5,N4)+VLWatMicP1_vr(N6,N5,N4))*dts_wat),tiny_wat)
      HeatByConvectVapFlux = (cpw*TKSoil1_vr(N6,N5,N4)+EvapLHTC)*ConvectVapFlux
    ENDIF
  ELSE
    ConvectVapFlux       = 0.0_r8
    HeatByConvectVapFlux = 0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine WaterVaporXgridFlow
!------------------------------------------------------------------------------------------

  subroutine SolveXgridHeatConduction(I,J,M,N,N1,N2,N3,N4,N5,N6,ConvectVapFlux,HeatByConvectVapFlux,HeatFluxAir2Soi)
  !
  !Description
  !Solve for heat flux due to conduction along temperature gradient
  implicit none
  integer , intent(in) :: I,J,M
  integer , intent(in) :: N
  integer , intent(in) :: N1,N2,N3  !source
  integer , intent(in) :: N4,N5,N6  !dest
  real(r8), intent(in) :: ConvectVapFlux          !
  real(r8), intent(in) :: HeatByConvectVapFlux    !convective heat flux out of source grid into dest grid
  real(r8), intent(in) :: HeatFluxAir2Soi        !ground heat into surface layer
  real(r8) :: TK1X,TKLX,TKY,HFLWC,HFLWX,HeatCondSoi
  real(r8) :: ATCNDL,DTKX,dCPv
  real(r8) :: ThermCondDst,ThermCondSrc
  !

  !     THERMAL CONDUCTIVITY IN EACH GRID CELL
  !
  !     DTH*,RYL*,DNU*,TRB*,XNUS*=turbulence effects on thermal conductivity
  !     THETWX,FracSoiPAsAir=water,air concentration
  !     TCNDW*,TCNDA*=thermal conductivity of water,air
  !     ThermCondSrc,TCNDL=soil thermal conductivity in source,destination
  !     WTHET*=multiplier for air concn in thermal conductivity
  !     ATCNDL=source-destination thermal conductance
  !
  dCPv=cpw*ConvectVapFlux
  DTKX=ABS(TKSoil1_vr(N3,N2,N1)-TKSoil1_vr(N6,N5,N4))

  call CalcSoilThermConductivity(N1,N2,N3,DTKX,ThermCondSrc)

  call CalcSoilThermConductivity(N4,N5,N6,DTKX,ThermCondDst)

  ATCNDL=safe_adb(2.0_r8*ThermCondSrc*ThermCondDst,(ThermCondSrc*DLYR_3D(N,N6,N5,N4)+ThermCondDst*DLYR_3D(N,N3,N2,N1)))

  !     HEAT FLOW FROM THERMAL CONDUCTIVITY AND TEMPERATURE GRADIENT
  !
  !     VHeatCapacity1_vr,VHCPW=volumetric heat capacity of soil,snowpack
  !     TK1X,TKLX=interim temperatures of source,destination
  !     ConvectiveHeatFlux,HeatFluxAir2Soi=convective heat from soil vapor flux
  !     HeatFluxAir2Soi=storage heat flux from snowpack
  !     TKY=equilibrium source-destination temperature
  !     HFLWC,HFLWX=source-destination heat flux unltd,ltd by heat
  !     ATCNDL=source-destination thermal conductance
  !     HeatCondSoi=source-destination conductive heat flux
  !     HFLWL=total conductive+convective source-destination heat flux
  !
  
  IF(VHeatCapacity1_vr(N3,N2,N1).GT.VHCPNX_col(N2,N1))THEN
    !surface layer and does not have significant snowpack
    IF(N3.EQ.NUM_col(N2,N1) .AND. VLHeatCapSnow_snvr(1,N2,N1).LE.VLHeatCapSnowMin_col(N2,N1))THEN
      
      TK1X=TKSoil1_vr(N3,N2,N1)-(HeatByConvectVapFlux-HeatFluxAir2Soi)/(VHeatCapacity1_vr(N3,N2,N1)-dCPv)
      
      if(abs(TK1X)>1.e5_r8)then
        write(*,*)'TKSoil1_vr(N3,N2,N1)-HeatByConvectVapFlux/VHeatCapacity1_vr(N3,N2,N1)',&
          TKSoil1_vr(N3,N2,N1),HeatByConvectVapFlux,HeatFluxAir2Soi,VHeatCapacity1_vr(N3,N2,N1)
        write(*,*)'N1,n2,n3',N1,N2,N3
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif

    ELSE

      TK1X=TKSoil1_vr(N3,N2,N1)-HeatByConvectVapFlux/(VHeatCapacity1_vr(N3,N2,N1)-dCPv)

      if(abs(TK1X)>1.e5_r8)then
        write(*,*)'TKSoil1_vr(N3,N2,N1)-HeatByConvectVapFlux/VHeatCapacity1_vr(N3,N2,N1)',&
          TKSoil1_vr(N3,N2,N1),HeatByConvectVapFlux,VHeatCapacity1_vr(N3,N2,N1)
        write(*,*)'N1,n2,n3',N1,N2,N3
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif
    ENDIF
  ELSE
    TK1X=TKSoil1_vr(N3,N2,N1)
  ENDIF

  IF(VHeatCapacity1_vr(N6,N5,N4).GT.ZEROS(N2,N1))THEN
    TKLX=TKSoil1_vr(N6,N5,N4)+HeatByConvectVapFlux/(VHeatCapacity1_vr(N6,N5,N4)+dCPv)
  ELSE
    TKLX=TKSoil1_vr(N6,N5,N4)
  ENDIF
  
  if(VHeatCapacity1_vr(N3,N2,N1)+VHeatCapacity1_vr(N6,N5,N4)>0._r8)then
    !temperature at the interface between (N3,N2,N1) and (N6,N5,N4)
    TKY=((VHeatCapacity1_vr(N3,N2,N1)-dCPv)*TK1X+(VHeatCapacity1_vr(N6,N5,N4)+dCPv)*TKLX) &
      /(VHeatCapacity1_vr(N3,N2,N1)+VHeatCapacity1_vr(N6,N5,N4))
  ELSE
    TKY=(TK1X+TKLX)/2._r8
  endif 

  HFLWX = (TK1X-TKY)*VHeatCapacity1_vr(N3,N2,N1)*dts_wat
  HFLWC = ATCNDL*(TK1X-TKLX)*AREA_3D(N,N3,N2,N1)*dts_HeatWatTP
  IF(HFLWC.GE.0.0_r8)THEN
    !from (N3,N2,N1) into (N6,N5,N4)
    HeatCondSoi=AZMAX1(AMIN1(HFLWX,HFLWC))
  ELSE
    !from (N6,N5,N4) into (N3,N2,N1)
    HeatCondSoi=AZMIN1(AMAX1(HFLWX,HFLWC))
  ENDIF

  HeatFlow2Soili_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)=HeatFlow2Soili_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)+HeatCondSoi

  end subroutine SolveXgridHeatConduction
!------------------------------------------------------------------------------------------
  subroutine DischargeOverWaterTable(I,J,N,N1,N2,N3,M4,M5,M6,DoMicPDischarg2ExtWTBL,DoMacPDischarg2ExtWTBL,&
    RechargDist2WTBL,Recharg2WTBLScal,Surf2MacPDepz,XN,TargetWaterTBL,FracLayVolBelowWTBL)
  !
  !Description:
  !  
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: N
  integer, intent(in) :: N1,N2,N3  !source grid
  integer, intent(in) :: M4,M5,M6  !dest grid
  logical, intent(in) :: DoMicPDischarg2ExtWTBL
  logical, intent(in) :: DoMacPDischarg2ExtWTBL

  real(r8), intent(in):: RechargDist2WTBL              !distance to external water table [m]
  real(r8),intent(in) :: Recharg2WTBLScal
  real(r8), intent(in):: XN                            !flow direction, <0 out of source grid
  real(r8), intent(in):: Surf2MacPDepz                         !depth to layer macropore water 
  real(r8), intent(in):: TargetWaterTBL                !vertical location of the external water table for tile
  real(r8), intent(in) :: FracLayVolBelowWTBL          !Fraction of soil layer below water table 

  character(len=*), parameter :: subname='DischargeOverWaterTable'
  real(r8) :: FLWT,MacPoreDischarg,FLWTH
  real(r8) :: PSISWD            !pressure head drop from the grid to the external 
                                !water table due to elevation change along the slope [MPa]
  real(r8) :: PSISWT,PSISWTH
  real(r8) :: watflx, heatflx
  real(r8) :: outlet_theta,outlet_conductivity
  integer :: outlet_k
  
  character(len=2), parameter :: direcs(2)=(/'ew','ns'/)

  call PrintInfo('beg '//subname)

  !     MICROPORE DISCHARGE ABOVE WATER TABLE
  !
  !     DoMicPDischarg2ExtWTBL=micropore discharge flag to natural water table
  !     PSISWD=water potential drop from water table slope
  !     XN,RCHG*=direction indicator,boundary flag
  !     SLOPE=tan(lateral slope), projected onto the grid direction
  !     DLYR=layer width
  !     WaterTable2GroundSlope_col=slope of water table relative to surface slope
  !     PSISWT=water potential driving micropore discharge
  !     PSISoilMatricPtmp_vr,PSISO=matric,osmotic water potential
  !     DPTH,ExtWaterTable=depth to layer midpoint,natural water table
  !     WaterTableDepZi_col=internal water table depth
  !     FLWL=micropore discharge to natural water table
  !     HFLWL=convective heat from discharge to natural water table
  !     HydroCond_3D=saturated hydraulic conductivity
  !     FracLayVolBelowExtWTBL_vr=fraction of layer below natural water table
  !
  IF(DoMicPDischarg2ExtWTBL .AND. (.not.isclose(Recharg2WTBLScal,0._r8)))THEN        
    !pressure drop due to the change of disturbance following the slope
    PSISWD = XN*0.5_r8*mGravAccelerat*SLOPE_col(N,N2,N1)*DLYR_3D(N,N3,N2,N1)*(1.0_r8-WaterTable2GroundSlope_col(N2,N1))
    
    if (use_ecosys_drainage) then
      ! ECOSIM_ZM: legacy watsub.f 5819--5834 includes the saturation
      ! potential and selects K(theta), not Ksat. Preserve MPa, m and the
      ! upstream timestep, area and flux-accounting conventions.
      PSISWT=AZMIN1(-PSISoilMatricPtmp_vr(N3,N2,N1)+PSI_at_Sautration_vr(N3,N2,N1) &
        +mGravAccelerat*(AMIN1(SoilDepthMidLay_vr(N3,N2,N1),WaterTableDepZi_col(N2,N1))-TargetWaterTBL))
        
      outlet_theta=AMAX1(SoilWatAirDry_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1), &
        safe_adb(VLWatMicP1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))
      outlet_k            = getMoistK(outlet_theta,POROS_vr(N3,N2,N1))
      outlet_conductivity = HydroCond_3D(N,outlet_k,N3,N2,N1)
    else
      !psi-g(z_i-zw), pressure gradient, assuming zero suction pressure at the external water table
      PSISWT = AZMIN1(-PSISoilMatricPtmp_vr(N3,N2,N1) + &
        mGravAccelerat*(AMIN1(SoilDepthMidLay_vr(N3,N2,N1),WaterTableDepZi_col(N2,N1))-TargetWaterTBL))   
      !set conductivity at saturation value  
      outlet_conductivity=HydroCond_3D(N,1,N3,N2,N1)
    endif

    PSISWT=PSISWT-PSISWD
    if(PSISWT<0._r8)then
      ! Upstream Ksat for tiles; ECOSIM_ZM K(theta) for natural outlets only.
      FLWT=PSISWT*outlet_conductivity*AREA_3D(N,N3,N2,N1)*(1.0_r8-FracLayVolBelowWTBL)/(RechargDist2WTBL+1.0_r8) &
        *Recharg2WTBLScal*dts_HeatWatTP
    else
      FLWT = 0._r8
    endif

    watflx                          = XN*FLWT
    heatflx                         = cpw*TKSoil1_vr(N3,N2,N1)*watflx
    WaterFlow2Micpt_3D(N,M6,M5,M4)  = WaterFlow2Micpt_3D(N,M6,M5,M4)+watflx
    WaterFlow2MicptX_3D(N,M6,M5,M4) = WaterFlow2MicptX_3D(N,M6,M5,M4)+watflx
    HeatFlow2Soili_3D(N,M6,M5,M4)   = HeatFlow2Soili_3D(N,M6,M5,M4)+heatflx

    QDrainM_col(N2,N1) = QDrainM_col(N2,N1)-XN*watflx
    QDrain_col(N2,N1)      = QDrain_col(N2,N1) -XN*watflx
    HeatDrain_col(N2,N1)   = HeatDrain_col(N2,N1)-XN*heatflx
    QDrain_cum_col(N2,N1)  = QDrain_cum_col(N2,N1) -XN*watflx
    QDrainloss_vr(N3,N2,N1) = QDrainloss_vr(N3,N2,N1)-XN*watflx

  ENDIF
  !
  !     MACROPORE DISCHARGE ABOVE WATER TABLE
  !
  !     DoMacPDischarg2ExtWTBL=macropore discharge flag to natural water table
  !     PSISWD=water potential from water table slope
  !     XN,RCHG*=direction indicator,boundary flag
  !     SLOPE=tan(lateral slope), projected onto the grid direction
  !     DLYR=layer width
  !     WaterTable2GroundSlope_col=water table slope
  !     PSISWTH=water potential driving macropore discharge
  !     PSISO=osmotic water potential
  !     ExtWaterTable=natural water table
  !     WaterTableDepZi_col=depth to internal water table
  !     FLWTH,MacPoreDischarg=macropore discharge unltd,ltd by macropore water
  !     HydroCondMacP1=macropore hydraulic conductivity
  !     WaterFlow2Macpt_3D=macropore discharge to natural water table
  !     HFLWL=convective heat from discharge to natural water table
  !     HydroCond_3D=saturated hydraulic conductivity
  !     FracLayVolBelowExtWTBL_vr=fraction of layer below natural water table
  !
  IF(DoMacPDischarg2ExtWTBL .AND. (.not.isclose(Recharg2WTBLScal,0._r8)) .AND. VLMacroPores1_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    
    PSISWD = XN*0.5_r8*mGravAccelerat*SLOPE_col(N,N2,N1)*DLYR_3D(N,N3,N2,N1)*(1.0_r8-WaterTable2GroundSlope_col(N2,N1))
    
    PSISWTH = mGravAccelerat*(AMIN1(Surf2MacPDepz,WaterTableDepZi_col(N2,N1))-TargetWaterTBL)

    PSISWTH = PSISWTH-PSISWD
    if(PSISWTH<0._r8)then
      FLWTH=PSISWTH*HydroCondMacP1_vr(N3,N2,N1)*AREA_3D(N,N3,N2,N1) &
        *(1.0_r8-FracLayVolBelowWTBL)/(RechargDist2WTBL+1.0_r8)*Recharg2WTBLScal*dts_HeatWatTP
    else
      FLWTH=0._r8
    endif

    MacPoreDischarg=AMAX1(FLWTH,AZMIN1(-(VLWatMacP1_vr(N3,N2,N1)*dts_wat &
      +WaterFlow2Macpt_3D(3,N3,N2,N1)-WaterFlow2Macpt_3D(3,N3+1,N2,N1))))

    watflx  = XN*MacPoreDischarg
    heatflx = cpw*TKSoil1_vr(N3,N2,N1)*watflx

    WaterFlow2Macpt_3D(N,M6,M5,M4) = WaterFlow2Macpt_3D(N,M6,M5,M4)+watflx
    HeatFlow2Soili_3D(N,M6,M5,M4)  = HeatFlow2Soili_3D(N,M6,M5,M4)+heatflx

    QDrainM_col(N2,N1) = QDrainM_col(N2,N1)-XN*watflx    
    QDrain_col(N2,N1)      = QDrain_col(N2,N1) -XN*watflx
    HeatDrain_col(N2,N1)   = HeatDrain_col(N2,N1)-XN*heatflx
    QDrain_cum_col(N2,N1)  = QDrain_cum_col(N2,N1) -XN*watflx
    QDrainloss_vr(N3,N2,N1) = QDrainloss_vr(N3,N2,N1) -XN*watflx
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine DischargeOverWaterTable
!------------------------------------------------------------------------------------------

  SUBROUTINE RechargeFromExtWaterTBL(I,J,M,N,N1,N2,N3,M4,M5,M6,Surf2MacPDepz,RechargDist2WTBL,&
    Recharg2WTBLScal,XN,AirfMicP,VOLPX2,AirfMacP,str_dir)!,donot_drain)
  !
  !subsurface recharge to soil micropore and macropores from external water table
  !it considers the existence of frozen layers
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: N
  integer, intent(in) :: N1,N2,N3   !source grid
  integer, intent(in) :: M4,M5,M6   !dest grid
  real(r8),intent(in) :: Recharg2WTBLScal
  real(r8),intent(in) :: Surf2MacPDepz         !vertical location of the macropore layer [m]
  real(r8),intent(in) :: RechargDist2WTBL      !distance between in soil water table and external water table [m]
  real(r8),intent(in) :: XN                    !flow direction, > 0 (in), < 0 (out)
  real(r8),intent(inout):: AirfMicP            !air-filled space in grid (N3,N2,N1)
  real(r8),intent(inout) :: VOLPX2,AirfMacP
  character(len=*),intent(in) :: str_dir  
!  logical, intent(in) :: donot_drain

  character(len=*), parameter :: subname='RechargeFromExtWaterTBL'
  real(r8) :: FLWU,FLWUL,FLWUX,FLWUH,FLWUHL  
  real(r8) :: PSISUT,PSISUTH,PSISWD
  real(r8) :: heatflx,watflx,watflxu

  call PrintInfo('beg '//subname)
  !     MICROPORE RECHARGE BELOW WATER TABLE

  IF(SoilDepthMidLay_vr(N3,N2,N1).GE.ExtWaterTable_col(N2,N1)                      & !deeper than water table
    .AND. ActiveLayDepZ_col(N2,N1).GT.ExtWaterTable_col(N2,N1)                     & !active layer deepd enough
    .AND. SoilDepthMidLay_vr(N3,N2,N1).LT.ActiveLayDepZ_col(N2,N1)                 & !inside active layer
    .AND. (AirfMicP.GT.ZEROS2(N2,N1) .OR. SoilBulkDensity_vr(N3,N2,N1).LE.ZERO)    & !able to accept water
    .AND. VLairMicP_vr(N3,N2,N1).GT.0.0_r8                                         & !grid is unsaturated
    .AND. (.not.isclose(Recharg2WTBLScal,0._r8)))THEN                                !water exchange with watertable enabled
    
    PSISWD = XN*0.5_r8*mGravAccelerat*SLOPE_col(N,N2,N1)*DLYR_3D(N,N3,N2,N1)*(1.0_r8-WaterTable2GroundSlope_col(N2,N1))
    PSISUT = AZMAX1(-PSISoilMatricPtmp_vr(N3,N2,N1) + &
      mGravAccelerat*(SoilDepthMidLay_vr(N3,N2,N1)-ExtWaterTable_col(N2,N1)))

    !outflow/drainage  (>0.)
    PSISUT=PSISUT+PSISWD
    IF(PSISUT.GT.0.0_r8)then
      FLWU = PSISUT*HydroCond_3D(N,1,N3,N2,N1)*AREA_3D(N,N3,N2,N1)*FracLayVolBelowExtWTBL_vr(N3,N2,N1) &
        /(RechargDist2WTBL+1.0_r8)*Recharg2WTBLScal*dts_HeatWatTP
    else
      FLWU = 0._r8
    endif

    !within a time step, the incoming flow cannot exceed avaiable air-filled pores   
    IF(SoilBulkDensity_vr(N3,N2,N1).GT.ZERO)THEN
      FLWUL = AMIN1(FLWU,AirfMicP)
      FLWUX = AMIN1(FLWU,VOLPX2)
    ELSE
      FLWUL = FLWU
      FLWUX = FLWU
    ENDIF
    watflx  = XN*FLWUL
    watflxu = XN*FLWUX
    heatflx = cpw*TKSoil1_vr(N3,N2,N1)*watflx
        
    WaterFlow2Micpt_3D(N,M6,M5,M4)  = WaterFlow2Micpt_3D(N,M6,M5,M4)+watflx
    WaterFlow2MicptX_3D(N,M6,M5,M4) = WaterFlow2MicptX_3D(N,M6,M5,M4)+watflxu
    HeatFlow2Soili_3D(N,M6,M5,M4)   = HeatFlow2Soili_3D(N,M6,M5,M4)+heatflx
    AirfMicP                        = AirfMicP-XN*watflx
    VOLPX2                          = VOLPX2-XN*watflxu

    QDischarg2WTBL_col(N2,N1) = QDischarg2WTBL_col(N2,N1)-XN*watflx
    QDischarM_col(N2,N1)      = QDischarM_col(N2,N1)-XN*watflx
    HeatDischar_col(N2,N1)    = HeatDischar_col(N2,N1)-XN*heatflx

  ENDIF

  !     MACROPORE RECHARGE BELOW WATER TABLE
  !
  IF(Surf2MacPDepz.GT.ExtWaterTable_col(N2,N1)                     & !deeper than water table
    .AND. ActiveLayDepZ_col(N2,N1).GT.ExtWaterTable_col(N2,N1)     & !active layer below water table
    .AND. SoilDepthMidLay_vr(N3,N2,N1).LT.ActiveLayDepZ_col(N2,N1) & !midlayer depth above active water layer
    .AND. AirfMacP.GT.ZEROS2(N2,N1)                                & !macropore has air-filled fraction
    .AND. (.not.isclose(Recharg2WTBLScal,0.0_r8)))THEN               !recharge is on

    PSISWD  = XN*0.5_r8*mGravAccelerat*SLOPE_col(N,N2,N1)*DLYR_3D(N,N3,N2,N1)*(1.0_r8-WaterTable2GroundSlope_col(N2,N1))
    PSISUTH = mGravAccelerat*(Surf2MacPDepz-ExtWaterTable_col(N2,N1))

    !outflow/drainage
    PSISUTH = PSISUTH+PSISWD
    IF(PSISUTH.GT.0.0_r8)then
      FLWUH  = PSISUTH*HydroCondMacP1_vr(N3,N2,N1)*AREA_3D(N,N3,N2,N1)*FracLayVolBelowExtWTBL_vr(N3,N2,N1) &
        /(RechargDist2WTBL+1.0_r8)*Recharg2WTBLScal*dts_HeatWatTP
    else
      FLWUH  = 0._r8
    endif

    FLWUHL                         = AMIN1(FLWUH,AirfMacP)
    watflx                         = XN*FLWUHL
    heatflx                        = cpw*TKSoil1_vr(N3,N2,N1)*watflx
    WaterFlow2Macpt_3D(N,M6,M5,M4) = WaterFlow2Macpt_3D(N,M6,M5,M4)+watflx
    HeatFlow2Soili_3D(N,M6,M5,M4)  = HeatFlow2Soili_3D(N,M6,M5,M4)+heatflx
    AirfMacP                       = AirfMacP-XN*watflx

    QDischarg2WTBL_col(N2,N1) = QDischarg2WTBL_col(N2,N1)-XN*watflx
    QDischarM_col(N2,N1)      = QDischarM_col(N2,N1)-XN*watflx
    HeatDischar_col(N2,N1)    = HeatDischar_col(N2,N1)-XN*heatflx
  ENDIF
  call PrintInfo('end '//subname)    
  end SUBROUTINE RechargeFromExtWaterTBL
!------------------------------------------------------------------------------------------

  subroutine FreezeThawIterateM(I,J,M,NY,NX,L,N1,N2,N3)
  !
  !Description
  !Applying freezethaw in cell (N3,N2,N1)
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NY,NX,L,N1,N2,N3

  character(len=*), parameter :: subname='FreezeThawIterateM'
  real(r8) :: VLHeatCapacityBX
  real(r8) :: PSISMX
  real(r8) :: TFREEZ           !freezing point 
  real(r8) :: VLWatMicP1X      !micropore liquid water before freeze/thaw
  real(r8) :: VLWatMacP1X      !macropore liquid water before freeze/thaw
  real(r8) :: VLHeatCapacityAX !
  real(r8) :: TK1App           !apparent temperature before freeze/thaw
  real(r8) :: MacPIceHeatFlxFrezPt,MacPIceHeatFlxFrez
  real(r8) :: ENGY1
  real(r8) :: MicPIceHeatFlxFrez   !actual latent heat from micropore freeze-thaw
  real(r8) :: MicPIceHeatFlxFrezPt !potential heat flux supporting micropore freeze-thaw
  real(r8) :: VLHeatCapacityX  !heat capacity before freeze/thaw
  real(r8) :: dcpo,rootC             !heat capacity due to living root [MJ/K/d-2]

  !
  !     FREEZE-THAW IN SOIL LAYER MICROPORE FROM NET CHANGE IN SOIL
  !     LAYER HEAT STORAGE
  !
  !     TFREEZ=micropore freezing temperature
  !     PSISoilAirEntry_vr,PSISO=micropore matric,osmotic potential
  !     VLWatMicP1*,VLiceMicP1=micropore water,ice volume
  !     VLWatMacP1*,VLiceMacP1=macropore water,ice volume
  !     VLHeatCapacityX,VLHeatCapacityAX,VLHeatCapacityBX=total soil,micropore,macropore heat capacity
  !     VHCM=soil solid volumetric heat capacity
  !     TK1*=soil temperature
  !     THFLWL=total soil conductive, convective heat flux
  !     HeatIrrigation1_vr=subsurface convective heat flux
  !     MicPIceHeatFlxFrezPt,MicPIceHeatFlxFrez=latent heat from micro freeze-thaw unltd,ltd by water,ice
  !     FIceThawedMicP_vr=soil water flux from micropore freeze-thaw
  !     MacPIceHeatFlxFrezPt,SoilWatFrezHeatRelease_vr=latent heat from macro freeze-thaw unltd,ltd by water,ice
  !     FIceThawedMacP_vr=soil water flux from macropore freeze-thaw
  !
  call PrintInfo('beg '//subname)
  if(plantOM4Heat)then
    rootC=sum(RootMycoMassElm_vr(ielmc,:,N3,N2,N1))
    dcpo = cpo*gOC_to_m3_OM(rootC)
  else
    dcpo=0._r8
  endif
  PSISMX          = PSISoilMatricPtmp_vr(N3,N2,N1)+PSISoilOsmotic_vr(N3,N2,N1)
  
  TFREEZ=get_Tfrez(PSISMX)  
  !obtain the potential water
  VLWatMicP1X     = VLWatMicP1_vr(N3,N2,N1)+WatNetFlow2MicP_3DM_vr(N3,N2,N1)+FWatExMacP2MicPiM_vr(N3,N2,N1)+FWatIrrigate2MicP1_vr(N3,N2,N1)
  VLWatMacP1X     = VLWatMacP1_vr(N3,N2,N1)+WatNetFlow2Macpt_3DM_vr(N3,N2,N1)-FWatExMacP2MicPiM_vr(N3,N2,N1)
  ENGY1           = VHeatCapacity1_vr(N3,N2,N1)*TKSoil1_vr(N3,N2,N1)
  VLHeatCapacityX = VHeatCapSolidSoil_vr(N3,N2,N1)+cpw*(VLWatMicP1X+VLWatMacP1X)+cpi*(VLiceMicP1_vr(N3,N2,N1)+VLiceMacP1_vr(N3,N2,N1))

  IF(VLHeatCapacityX.GT.ZEROS(NY,NX))THEN
    VLHeatCapacityX=VLHeatCapacityX+dcpo
    TK1App=(ENGY1+THeatFlow2Soil_3DM_vr(N3,N2,N1)+HeatIrrigation1_vr(N3,N2,N1))/VLHeatCapacityX
!
!     FREEZE-THAW IN SOIL LAYER MACROPORE FROM NET CHANGE IN SOIL
!     LAYER HEAT STORAGE
!   TFice=frozen temperature
!   There is freeze-thaw in macropore, it is assumed there is no freezing temperature depression in 
!   macropores, therefore, water will first freeze in macropores.
    IF((TK1App.LT.TFREEZ .AND. VLWatMacP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)) &
      .OR.(TK1App.GT.TFREEZ .AND. VLiceMacP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)))THEN
      
      VLHeatCapacityBX     = cpw*VLWatMacP1X+cpi*VLiceMacP1_vr(L,NY,NX)
      !where is the following equation come from?
      MacPIceHeatFlxFrezPt = VLHeatCapacityBX*(TFREEZ-TK1App)/((1.0_r8+6.2913E-03_r8*(TFREEZ-TFICE)))*dts_wat 
!          /(1.0_r8-0._r8*0.10_r8*PSISMX)

      !Ice thawed, absorb heat
      IF(MacPIceHeatFlxFrezPt.LT.0.0_r8)THEN        
        MacPIceHeatFlxFrez=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMacP1_vr(N3,N2,N1)*dts_wat,MacPIceHeatFlxFrezPt)        
      !Water frozen, release heat
      ELSE        
        MacPIceHeatFlxFrez=AMIN1(LtHeatIceMelt*VLWatMacP1X*dts_wat,MacPIceHeatFlxFrezPt)
      ENDIF
      FIceThawedMacP_vr(N3,N2,N1)=-MacPIceHeatFlxFrez/LtHeatIceMelt
    !No freeze-thaw in macropore 
    ELSE
      MacPIceHeatFlxFrez     = 0.0_r8
      FIceThawedMacP_vr(N3,N2,N1) = 0.0_r8
    ENDIF

    !freeze-thaw in micropore
    IF((TK1App.LT.TFREEZ .AND. VLWatMicP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)) &
      .OR.(TK1App.GT.TFREEZ .AND. VLiceMicP1_vr(N3,N2,N1).GT.ZERO*VGeomLayer_vr(N3,N2,N1)))THEN
      
      VLHeatCapacityAX     = VHeatCapSolidSoil_vr(N3,N2,N1)+cpw*VLWatMicP1X+cpi*VLiceMicP1_vr(N3,N2,N1)
      if(VLHeatCapacityAX>0._r8)VLHeatCapacityAX=VLHeatCapacityAX+dcpo
      MicPIceHeatFlxFrezPt = VLHeatCapacityAX*(TFREEZ-TK1App)/((1.0_r8+6.2913E-03_r8*(TFREEZ-TFice)))*dts_wat

      !fusion energy absorb (<0) in thaw
      IF(MicPIceHeatFlxFrezPt.LT.0.0_r8)THEN
        MicPIceHeatFlxFrez=AMAX1(-LtHeatIceMelt*DENSICE*VLiceMicP1_vr(N3,N2,N1)*dts_wat,MicPIceHeatFlxFrezPt)
      !fusion energy release (>0) in freeze
      ELSE
        MicPIceHeatFlxFrez=AMIN1(LtHeatIceMelt*VLWatMicP1X*dts_wat,MicPIceHeatFlxFrezPt)
      ENDIF
      FIceThawedMicP_vr(N3,N2,N1)=-MicPIceHeatFlxFrez/LtHeatIceMelt
    ELSE
      MicPIceHeatFlxFrez     = 0.0_r8
      FIceThawedMicP_vr(N3,N2,N1) = 0.0_r8
    ENDIF

  ELSE
    MicPIceHeatFlxFrez          = 0.0_r8
    FIceThawedMicP_vr(N3,N2,N1) = 0.0_r8
    MacPIceHeatFlxFrez          = 0.0_r8
    FIceThawedMacP_vr(N3,N2,N1) = 0.0_r8
  ENDIF
  !total heat released due to phase change 
  SoilWatFrezHeatRelease_vr(N3,N2,N1)=MicPIceHeatFlxFrez+MacPIceHeatFlxFrez !+ &
    !(cpw-cpi*DENSICE)*(FIceThawedMicP_vr(N3,N2,N1)+FIceThawedMacP_vr(N3,N2,N1))*TFREEZ
  
  !
  !     TOTAL AND HOURLY ACCUMULATED FREEZE-THAW FLUXES
  !
  !     THAW,TLIceThawMacP=hourly accumulated freeze-thaw flux in micropores,macropores
  !     HTHAW=hourly accumulated freeze-thaw latent heat flux
  !     TWFLXL,TMLiceThawMacP=total accumulated freeze-thaw in micropores,macropores
  !     TLPhaseChangeHeat2Soi1_vr=total latent heat flux from melting

  ! Summarize fluxes for the current iteration
  TLIceThawMicP_vr(N3,N2,N1)          = TLIceThawMicP_vr(N3,N2,N1)+FIceThawedMicP_vr(N3,N2,N1)
  TLIceThawMacP_vr(N3,N2,N1)          = TLIceThawMacP_vr(N3,N2,N1)+FIceThawedMacP_vr(N3,N2,N1)
  TLPhaseChangeHeat2Soi_vr(N3,N2,N1)  = TLPhaseChangeHeat2Soi_vr(N3,N2,N1)+SoilWatFrezHeatRelease_vr(N3,N2,N1)
  TMLiceThawedMicP_vr(N3,N2,N1)       = TMLiceThawedMicP_vr(N3,N2,N1)+FIceThawedMicP_vr(N3,N2,N1)
  TMLiceThawedMacP_vr(N3,N2,N1)       = TMLiceThawedMacP_vr(N3,N2,N1)+FIceThawedMacP_vr(N3,N2,N1)
  TLPhaseChangeHeat2Soi1_vr(N3,N2,N1) = TLPhaseChangeHeat2Soi1_vr(N3,N2,N1)+SoilWatFrezHeatRelease_vr(N3,N2,N1)

  call PrintInfo('end '//subname)
  end subroutine FreezeThawIterateM

end module WatsubMod
