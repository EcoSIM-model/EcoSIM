module SnowPhysMod
! Description:
! Evolve snow water, ice, temperature and surface exchange, and accumulate
! accepted mass and energy transfers for the surface and hourly balances.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL
  use abortutils,    only: endrun, iulog
  use EcoSIMCtrlMod, only: fixWaterLevel, mod_snow_albedo
  use DebugToolMod,  only: PrintInfo
  use SnowDataType
  use SurfLitterDataType
  use SOMDataType
  use GridDataType
  use MiniMathMod
  use ClimForcDataType
  use EcoSimConst
  use SnowPhysData
  use EcoSIMSolverPar
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use HydroThermData
  use SurfaceAirExchangeMod, only: SurfaceExchangeResistance
  use SurfLitterDataType
  use SurfSoilDataType
  use SoilHeatDataType
  use SoilPropertyDataType
  use SoilPhysDataType
  use SoilPhysParaMod, only: ComputePsiMCM
  use ChemTranspDataType
  use LandSurfDataType
  use PhysPars
  use UnitMod, only : units
  use ChemTracerParsMod
  use MiniFuncMod
  use EcoSIMCtrlMod
implicit none
  private
  character(len=*), parameter :: mod_filename=&
  __FILE__
  integer :: MMit
  public :: InitSnowLayers
  public :: SnowRedistributionM
  public :: CalcSnowBNDResistance
  public :: ZeroSnowFluxM
  public :: PrepIterSnowLayerM
  public :: InitSnowAccumsM
  public :: StageSnowModel
  public :: LimitSnowPhaseChange
  public :: SolveSnowpackM
  public :: SumSnowDriftByRunoffM
  public :: UpdateSnowAtM
  public :: UpdateSnowPack1M
  public :: AccumulateSnowRedisFluxM
contains
!------------------------------------------------------------------------------------------
  subroutine InitSnowLayers(NY,NX,XSW)

  implicit none
  integer, intent(in) :: NY,NX
  real(r8), optional, intent(out) :: XSW  !total water equivalent snow
  character(len=*), parameter :: subname='InitSnowLayers'
  real(r8) :: DLYRSI
  real(r8) :: VOLSWI
  real(r8) :: GridArea
  !reference snow thickness [m]
  real(r8), parameter :: cumSnowDepzRef_col(JS)=(/0.05_r8,0.15_r8,0.30_r8,0.60_r8,1.00_r8/)
  ! The maximum initial snowpack depth is 1.0 m.
  integer :: L

  ! begin_execution
  call PrintInfo('beg '//subname)
  !
  ! cumSnowDepz_col=depth to bottom
  ! DENSICE=ice density, 0.92 Mg m-3
  ! NewSnowDens=snow density (Mg m-3)
  ! VOLSS,VcumWatSnow_col,VOLIS,VOLS=snow,water,ice,total snowpack volume(m3)

! cumSnowDepzRef_col=depth to bottom of snowpack layers (m), i.e. from current layere to surface
! DLYRS=snowpack layer thickness (m)
! VOLSSL,VcumWatSnow_colL,VOLISL,VLSnoDWIprev_snvr=snow,water,ice,total layer volume(m3), water equivalent snow
! DENSS=layer density (Mg m-3)
! TKW,TCSnow=later temperature K,oC
! VHCPW=layer volumetric heat capacity (MJ m-3 K-1)
! SnowDepth=total snow height in the column

  ! Truncate before deriving any stores or summaries so initialization agrees.
  if(SnowDepth_col(NY,NX)>cumSnowDepzRef_col(JS))then
    write(iulog,*) 'WARNING: InitSnowLayers NY,NX=',NY,NX, &
      ': initial snow depth (m)=',SnowDepth_col(NY,NX), &
      ' truncated to (m)=',cumSnowDepzRef_col(JS)
    SnowDepth_col(NY,NX)=cumSnowDepzRef_col(JS)
  endif

  GridArea=DH_col(NY,NX)*DV_col(NY,NX)
  if(present(XSW))XSW=0._r8
  cumSnowDepz_col(0,NY,NX) = 0.0_r8
  NewSnowDens_col(NY,NX)   = 0.10_r8   ![Mg m-3]
  VcumDrySnoWE_col(NY,NX)  = SnowDepth_col(NY,NX)*NewSnowDens_col(NY,NX)*GridArea
  VcumWatSnow_col(NY,NX)   = 0.0_r8
  VcumIceSnow_col(NY,NX)   = 0.0_r8
  VcumSnoDWI_col(NY,NX)    = VcumDrySnoWE_col(NY,NX)/NewSnowDens_col(NY,NX)+VcumWatSnow_col(NY,NX)+VcumIceSnow_col(NY,NX)
!  VOLSWI=0.0_r8

  !build the snow profile, topdown
  nsnol_col(NY,NX)=0
  D9580: DO L=1,JS
    IF(L.EQ.1)THEN
      !top snow layer
      DLYRSI                   = cumSnowDepzRef_col(L)
      SnowThickL_snvr(L,NY,NX) = AMIN1(DLYRSI,SnowDepth_col(NY,NX))
    ELSE
      DLYRSI                   = cumSnowDepzRef_col(L)-cumSnowDepzRef_col(L-1)
      SnowThickL_snvr(L,NY,NX) = AMIN1(DLYRSI,AZMAX1(SnowDepth_col(NY,NX)-cumSnowDepzRef_col(L-1)))
    ENDIF
    VLSnoDWIMax_snvr(L,NY,NX)=DLYRSI*GridArea      !maximum snow volume alloed in layer L

    if(SnowThickL_snvr(L,NY,NX) .GT. ZEROS(NY,NX))then
      !there is meaningful snow
      nsnol_col(NY,NX)=nsnol_col(NY,NX)+1
      VLDrySnoWE_snvr(L,NY,NX) = SnowThickL_snvr(L,NY,NX)*NewSnowDens_col(NY,NX)*GridArea
      VLWatSnow_snvr(L,NY,NX)  = 0.0_r8
      VLIceSnow_snvr(L,NY,NX)  = 0.0_r8

  !    IF(L.EQ.1)THEN
  !      VOLSWI=VOLSWI+0.5_r8*(VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE)
  !    ELSE
  !      VOLSWI=VOLSWI+0.5_r8*(VLDrySnoWE_snvr(L-1,NY,NX)+VLWatSnow_snvr(L-1,NY,NX) &
  !        +VLIceSnow_snvr(L-1,NY,NX)*DENSICE+VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX) &
  !        +VLIceSnow_snvr(L,NY,NX)*DENSICE)
  !    ENDIF

      SnoDens_snvr(L,NY,NX)      = NewSnowDens_col(NY,NX)
      VLSnoDWIprev_snvr(L,NY,NX) = VLDrySnoWE_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)
      cumSnowDepz_col(L,NY,NX)    = cumSnowDepz_col(L-1,NY,NX)+SnowThickL_snvr(L,NY,NX)
      TKSnow_snvr(L,NY,NX)        = AMIN1(Tref,TairKClimMean_col(NY,NX))
      TCSnow_snvr(L,NY,NX)        = TKSnow_snvr(L,NY,NX)-Tref
      VLHeatCapSnow_snvr(L,NY,NX) = cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)
    ELSE
      VLDrySnoWE_snvr(L,NY,NX)    = 0._r8
      VLWatSnow_snvr(L,NY,NX)     = 0.0_r8
      VLIceSnow_snvr(L,NY,NX)     = 0.0_r8
      SnoDens_snvr(L,NY,NX)       = NewSnowDens_col(NY,NX)
      VLSnoDWIprev_snvr(L,NY,NX)  = 0._r8
      cumSnowDepz_col(L,NY,NX)    = cumSnowDepz_col(L-1,NY,NX)+SnowThickL_snvr(L,NY,NX)
      TKSnow_snvr(L,NY,NX)        = TKS_vr(NU_col(NY,NX),NY,NX)
      TCSnow_snvr(L,NY,NX)        = TKSnow_snvr(L,NY,NX)-Tref
      VLHeatCapSnow_snvr(L,NY,NX) = 0._r8
    ENDIF
    if(present(XSW))XSW=XSW+VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
    
  ENDDO D9580
  !
  ! Retain the historical threshold for layer merging/removal. Thermal
  ! exchange uses the numerical zero so thin positive stores remain active.
  !
  VLHeatCapSnowMin_col(NY,NX)=VLHeatCapSnoMin*AREA_3D(3,NU_col(NY,NX),NY,NX)
  call PrintInfo('end '//subname)
  end subroutine InitSnowLayers

!------------------------------------------------------------------------------------------

  subroutine InitSnowAccumsM(I,J,M,NY,NX)
  !
  !called each iterations
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NY,NX
  character(len=*), parameter :: subname='InitSnowAccumsM'
  integer :: L

  !  TFLWS,TFLWW,TFLWI=net fluxes of snow,water,ice in snowpack
  !  CumHeat2SnowLay=convective heat fluxes of snow,water,ice in snowpack
  !
  call PrintInfo('beg '//subname)
  D9875: DO L=1,JS
    CumSno2SnowLM_snvr(L,NY,NX)      = 0.0_r8
    CumWat2SnowLM_snvr(L,NY,NX)      = 0.0_r8
    CumIce2SnowLM_snvr(L,NY,NX)      = 0.0_r8
    CumHeat2SnowLM_snvr(L,NY,NX)     = 0.0_r8
    XSnowThawMassLM_snvr(L,NY,NX)    = 0.0_r8
    XIceThawMassLM_snvr(L,NY,NX)     = 0.0_r8
    XPhaseChangeHeatLM_snvr(L,NY,NX) = 0.0_r8

    SnoX2SnoLay_snvr(L,NY,NX)      = 0.0_r8
    WatX2SnoLay_snvr(L,NY,NX)      = 0.0_r8
    IceX2SnoLay_snvr(L,NY,NX)      = 0.0_r8
    HeatX2SnoLay_snvr(L,NY,NX)     = 0.0_r8
    WatFlowInSnowM_snvr(M,L,NY,NX) = 0.0_r8
  ENDDO D9875
  !
  ! SURFACE FLUX ACCUMULATORS
  !
  ! TWFLXL,TMLiceThawMacP=total freeze-thaw in micropores,macropores
  ! TLPhaseChangeHeat2Soi=total latent heat from freeze-thaw
  ! TFLWL,TFLWHL=net water flux in micropores,macropores
  ! THFLWL=net heat flux
  !
  !
  ! ENERGY EXCHANGE VARIABLES AT SNOW SURFACE IF PRESENT
  !
  ! Radnet2Snow,LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatSensAir2Snow,HeatNetFlx2Snow=netradn,latent,convective,sensible
  ! and storage heat fluxes
  ! CumWatFlx2SoiMacP=water from snowpack to soil micropores,macropores
  ! cumNetHeatFlow2Soil=conv heat from snowpack to soil micropores,macropores
  !
  ! EVAPS,EVAPW=evaporation from soil,snowpack surfaces
  ! WatFlowSno2LitRM_col,WatFlowSno2MicPM_col,WatFlowSno2MacPM_col=water into litter,soil micropores,micropores for use in TranspNoSalt.f
  !
  VapXAir2Sno_col(NY,NX)=0._r8;EVAPS_col(NY,NX)=0.0_r8;EVAPW_col(NY,NX)=0.0_r8
  WatFlowSno2LitRM_col(M,NY,NX)=0.0_r8
  WatFlowSno2MicPM_col(M,NY,NX)=0.0_r8
  WatFlowSno2MacPM_col(M,NY,NX)=0.0_r8

  call PrintInfo('end '//subname)
  end subroutine InitSnowAccumsM


!------------------------------------------------------------------------------------------

  subroutine SnowPackIterationMM(dt_SnoHeat,I,J,M,NY,NX,TotWatXFlx2SoiMicP,TotHeatFlow2Soi,WatFlowSno2MacP,&
    TotSnoWatFlow2Litr,TotSnoHeatFlow2Litr,CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,&
    CumWatXFlx2SoiMicP,CumSnowWatFLow2LitR,CumNetHeatFlow2LitR,cumNetHeatFlow2Soil)
  implicit none
  real(r8), intent(in) :: dt_SnoHeat             !time step size for snow iteration, [h]
  integer,  intent(in) :: M,NY,NX,I,J
  real(r8), intent(out) :: TotWatXFlx2SoiMicP,TotHeatFlow2Soi,WatFlowSno2MacP
  real(r8), intent(out) :: TotSnoWatFlow2Litr,TotSnoHeatFlow2Litr
  real(r8), intent(inout) :: CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP
  real(r8), intent(inout) :: CumSnowWatFLow2LitR
  real(r8), intent(inout) :: CumNetHeatFlow2LitR
  real(r8), intent(inout) :: cumNetHeatFlow2Soil
  character(len=*), parameter :: subname='SnowPackIterationMM'

  real(r8) :: TCND1W,TCNDR
  real(r8) :: ATCNDW,VapCondSnoWeited,TCNDS
  real(r8) :: H2OVapFlx,H2OVapFlxMax

  real(r8) :: VapConvFlxInSnow,HeatbyVapConvInSnow
  real(r8) :: HeatCnduct,HeatCnductMax
  real(r8) :: PSISV1,VLairSno1,VLairSno2,VPY
  real(r8) :: TKY,FracAsAirSno1,AirFilledFracSnow
  real(r8) :: FCDX,LOGFCX,FCX,VapFlxSno2Soi1,HeatConvFlxSno2Soi1
  real(r8) :: VapCond1,VapCond2,CNVR,DENSW1,DENSW2
  real(r8) :: WatVapFloInSnow,HeatByWatVapFloInSnow,WatFlowSno2Soil,WatFlowSno2MicP,PtWatFlowSno2Soi
  real(r8) :: WatFloInSno,HeatFlxByWatFloInSno
  real(r8) :: WatFlowSno2LitR        !water flow from bottom snow layer to litter
  real(r8) :: HeatFlowSno2LitrByWat  !heat flow from bottom snow layer to litter due to water flow
  real(r8) :: WatFloInSnoMax
  real(r8) :: LiquidBudget(JS)
  real(r8), parameter :: mscal=0.99999_r8
  real(r8) :: HeatFlowSno2SoiByWat   !heat flow from bottom snow layer to soil due to water flow
  !real(r8) :: HAltFlowSno2SoiByWat
  real(r8) :: HeatCndFlxSno2Soi,HeatCndFlxInSno,PSDX,TCND2W,THETRR
  real(r8) :: TK0X,TK1X,VapSnoSrc,VapSnoDest,WPX
  integer :: L,L2
  logical :: lchkBottomL

  ! begin_execution
  ! PHYSICAL AND HYDRAULIC PROPERTIES OF SNOWPACK INCLUDING
  ! AIR AND WATER-FILLED POROSITY, WATER POTENTIAL OF UNDERLYING
  ! SOIL SURFACE USED IN FLUX CALCULATIONS
  !
  ! VHCPW,VLHeatCapSnowMin_col=current, minimum snowpack heat capacities
  ! VLDrySnoWE0M,VOLI0M,VOLW0M,VLSnoDWI1=snow,ice,water,total snowpack volume
  ! DENSS,DENSICE,NewSnowDens=snow,ice,minimum snow density
  ! AREA=area of grid cell
  ! VLairSno1=snowpack air volume
  ! FracAsAirSno1=snowpack air concentration
  ! VapCond1=snowpack vapor conductivity
  ! VapSnoSrc=snowpack vapor concentration
  ! TKSnow1=snowpack temperature
  ! WGSGW=vapor diffusivity
  ! DENSW1=snowpack density
  ! VHCPWMM=previous snowpack heat capacity
  call PrintInfo('beg '//subname)

  lchkBottomL=.false.
  ! Share each donor's liquid among all outgoing faces before accumulating
  ! water, heat or tracer fluxes. Reserve atmospheric evaporation first;
  ! positive incoming transfers are not borrowed for this iteration's exports.
  LiquidBudget    = AZMAX1(VLWatSnow0M_snvr(1:JS,NY,NX))
  LiquidBudget(1) = AZMAX1(LiquidBudget(1)+MIN(0._r8,WatX2SnoLay_snvr(1,NY,NX)))
  LiquidBudget    = LiquidBudget*mscal
  !loop from surface (L=1) to bottom (L=JS, maximum)
  D9880: DO L=1,JS

    IF(VLHeatCapSnowM1_snvr(L,NY,NX).GT.ZEROS(NY,NX))THEN
      ! active snow layer
      VLSnoDWI1_snvr(L,NY,NX)   = VLDrySnoWE0M_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX)+VLWatSnow0M_snvr(L,NY,NX)+VLIceSnow0M_snvr(L,NY,NX)
      SnowThickL0_snvr(L,NY,NX) = VLSnoDWI1_snvr(L,NY,NX)/AREA_3D(3,NUM_col(NY,NX),NY,NX)
      VLairSno1                 = AZMAX1(VLSnoDWI1_snvr(L,NY,NX)-VLDrySnoWE0M_snvr(L,NY,NX)-VLIceSnow0M_snvr(L,NY,NX)-VLWatSnow0M_snvr(L,NY,NX))
      FracAsAirSno1             = AMAX1(THETPI,VLairSno1/VLSnoDWI1_snvr(L,NY,NX))
      VapCond1                  = FracAsAirSno1**2*H2OVapDifsc_snvr(L,NY,NX)
      VapSnoSrc                 = vapsat(TKSnow1_snvr(L,NY,NX))  !*dssign(VLWatSnow0M_snvr(L,NY,NX))

      IF(VLSnoDWI1_snvr(L,NY,NX).GT.0._r8)THEN
        !maximum snow density is 0.6 g/cm3
        DENSW1=AMIN1(0.6_r8,(VLDrySnoWE0M_snvr(L,NY,NX)+VLWatSnow0M_snvr(L,NY,NX)+VLIceSnow0M_snvr(L,NY,NX)*DENSICE)/VLSnoDWI1_snvr(L,NY,NX))
      ELSE
        DENSW1=NewSnowDens_col(NY,NX)
      ENDIF
      !
      ! SNOW THERMAL CONDUCTIVITY FROM J GLACIOL 43:26-41
      !
      ! TCND1W=snow thermal conductivity
      ! WatFloInSnoMax=porosity-unconstrained snow water flux
      !
      TCND1W=0.0036_r8*10._r8**(2.650_r8*DENSW1-1.652_r8)
      !
      ! DISCHARGE OF MELTWATER AND ITS HEAT FROM SNOWPACK LAYER
      ! TO LOWER SNOWPACK LAYER
      ! the value 0.05 below is arbitrary,
      WatFloInSnoMax=AZMAX1(AZMAX1(VLWatSnow0M_snvr(L,NY,NX))-0.05_r8*AZMAX1(VLDrySnoWE0M_snvr(L,NY,NX)))*dts_sno
      WatFloInSnoMax=MIN(WatFloInSnoMax,LiquidBudget(L))
      !
      ! WATER AND HEAT FLUXES IN SNOWPACK
      !
      ! SnowThickL0_snvr=snow layer thickness
      ! WatFloInSno=porosity-constrained snow water flux
      ! HeatFlxByWatFloInSno=convective heat flux from water flux
      ! id of next/target snow layer
      ! maximum JS snow layers
      ! next layer is L2
      L2=MIN(JS,L+1)

      !not bottom layer, and it is heat significant
      IF(L.LT.JS .AND. VLHeatCapSnowM1_snvr(L2,NY,NX).GT.ZEROS(NY,NX))THEN
        !if L==JS-1, L2==JS, so top layer is treated here.

        VLSnoDWI1_snvr(L2,NY,NX)   = VLDrySnoWE0M_snvr(L2,NY,NX)/SnoDens_snvr(L2,NY,NX)+VLWatSnow0M_snvr(L2,NY,NX)+VLIceSnow0M_snvr(L2,NY,NX)
        SnowThickL0_snvr(L2,NY,NX) = VLSnoDWI1_snvr(L2,NY,NX)/AREA_3D(3,NUM_col(NY,NX),NY,NX)
        VLairSno2                  = AZMAX1(VLSnoDWI1_snvr(L2,NY,NX)-VLDrySnoWE0M_snvr(L2,NY,NX)-VLIceSnow0M_snvr(L2,NY,NX)-VLWatSnow0M_snvr(L2,NY,NX))
        AirFilledFracSnow          = AMAX1(THETPI,VLairSno2/VLSnoDWI1_snvr(L2,NY,NX))

        !water flows only in air-filled pores
        WatFloInSno          = AMIN1(VLairSno2,WatFloInSnoMax)
        LiquidBudget(L)      = AZMAX1(LiquidBudget(L)-WatFloInSno)
        HeatFlxByWatFloInSno = cpw*TKSnow1_snvr(L,NY,NX)*WatFloInSno
        !
        ! VAPOR FLUX IN SNOWPACK
        !
        ! VLairSno1,VLairSno2=air-filled volumes of source, destination layers
        ! L2=destination layer
        ! VapCond1,VapCond2=vapor conductivities of source, destination layers
        ! VapSnoSrc,VapSnoDest=vapor concentrations of source, destination layers
        ! TKSnow1=soil temperature
        ! VapCondSnoWeited=snow vapor conductance
        ! SnowThickL0_snvr=snow layer thickness
        ! H2OVapFlx,H2OVapFlxMax=vapor-unconstrained,vapor-constrained vapor flux
        ! VapConvFlxInSnow,HeatbyVapConvInSnow=vapor flux and its convective heat flux
        !
        IF(VLairSno1.GT.ZEROS2(NY,NX) .AND. VLairSno2.GT.ZEROS2(NY,NX))THEN
          !both layer has air-filled pores
          VapCond2         = AirFilledFracSnow**2*H2OVapDifsc_snvr(L2,NY,NX)
          VapSnoDest       = vapsat(TKSnow1_snvr(L2,NY,NX))  !*dssign(VLWatSnow0M_snvr(L2,NY,NX))
          VapCondSnoWeited = 2.0_r8*VapCond1*VapCond2/(VapCond1*SnowThickL0_snvr(L2,NY,NX)+VapCond2*SnowThickL0_snvr(L,NY,NX))
          H2OVapFlx        = VapCondSnoWeited*(VapSnoSrc-VapSnoDest)*AREA_3D(3,NUM_col(NY,NX),NY,NX)*FracSurfSnowCovered_col(NY,NX)*dt_SnoHeat
          VPY              = (VapSnoSrc*VLairSno1+VapSnoDest*VLairSno2)/(VLairSno1+VLairSno2)
          !
          !H2OVapFlxMax>0, moves out from layer L
          H2OVapFlxMax=(VapSnoSrc-VPY)*VLairSno1*dts_sno
          !out of layer L
          IF(H2OVapFlx.GE.0.0_r8)THEN
            VapConvFlxInSnow    = AZMAX1(AMIN1(H2OVapFlx,H2OVapFlxMax,VLWatSnow0M_snvr(L,NY,NX)*dts_wat,LiquidBudget(L)))
            HeatbyVapConvInSnow = (cpw*TKSnow1_snvr(L,NY,NX)+EvapLHTC)*VapConvFlxInSnow
          !into layer L
          ELSE
            VapConvFlxInSnow    = AZMIN1(AMAX1(H2OVapFlx,H2OVapFlxMax,-VLWatSnow0M_snvr(L2,NY,NX)*dts_wat,-LiquidBudget(L2)))
            HeatbyVapConvInSnow = (cpw*TKSnow1_snvr(L2,NY,NX)+EvapLHTC)*VapConvFlxInSnow
          ENDIF
        ELSE
          VapConvFlxInSnow    = 0.0_r8
          HeatbyVapConvInSnow = 0.0_r8
        ENDIF
        LiquidBudget(L)=MAX(0._r8,LiquidBudget(L)-MAX(0._r8,VapConvFlxInSnow))
        LiquidBudget(L2)=MAX(0._r8,LiquidBudget(L2)+MIN(0._r8,VapConvFlxInSnow))
        !
        ! HEAT FLUX IN SNOWPACK
        !
        ! DENSW2,TCNDW2=density,thermal conductivity in destination layer
        ! ATCNDW=thermal conductance
        ! SnowThickL0_snvr=layer thickness
        ! TKY=equilibrium temperature
        ! HeatCnductMax,HeatCnduct=heat-constrained,heat-unconstrained heat fluxes
        ! VHCPWMM,TKSnow1=volumetric heat capacity,temperature
        ! dts_wat=time step for flux calculations
        ! FSNW=snow cover fraction
        ! dt_SnoHeat=time step for snowpack flux calculations
        ! HeatCndFlxInSno=snowpack heat flux
        ! FLW0S,FLQ0I,FLQ0W=snow,ice,water fluxes through snowpack
        ! HFLW0W=convective heat flux snow,water,ice fluxes
        !
        IF(VLSnoDWI1_snvr(L2,NY,NX).GT.ZEROS2(NY,NX))THEN
          DENSW2=AMIN1(0.6_r8,(VLDrySnoWE0M_snvr(L2,NY,NX)+VLWatSnow0M_snvr(L2,NY,NX) &
            +VLIceSnow0M_snvr(L2,NY,NX)*DENSICE)/VLSnoDWI1_snvr(L2,NY,NX))
        ELSE
          !initialize new snow layer
          DENSW2=NewSnowDens_col(NY,NX)
        ENDIF
        TCND2W = 0.0036_r8*10._r8**(2.650_r8*DENSW2-1.652_r8)
        ATCNDW = 2.0_r8*TCND1W*TCND2W/(TCND1W*SnowThickL0_snvr(L2,NY,NX)+TCND2W*SnowThickL0_snvr(L,NY,NX))
        TKY    = (TKSnow1_snvr(L,NY,NX)*VLHeatCapSnowM1_snvr(L,NY,NX)+TKSnow1_snvr(L2,NY,NX) &
          *VLHeatCapSnowM1_snvr(L2,NY,NX))/(VLHeatCapSnowM1_snvr(L,NY,NX)+VLHeatCapSnowM1_snvr(L2,NY,NX))

        HeatCnductMax = (TKSnow1_snvr(L,NY,NX)-TKY)*VLHeatCapSnowM1_snvr(L,NY,NX)*dts_sno
        HeatCnduct    = ATCNDW*(TKSnow1_snvr(L,NY,NX)-TKSnow1_snvr(L2,NY,NX))*AREA_3D(3,NUM_col(NY,NX),NY,NX) &
          *FracSurfSnowCovered_col(NY,NX)*dt_SnoHeat

        IF(HeatCnduct.GE.0.0_r8)THEN
          HeatCndFlxInSno=AZMAX1(AMIN1(HeatCnductMax,HeatCnduct))
        ELSE
          HeatCndFlxInSno=AZMIN1(AMAX1(HeatCnductMax,HeatCnduct))
        ENDIF

        WatVapFloInSnow                 = WatFloInSno+VapConvFlxInSnow
        HeatByWatVapFloInSnow           = HeatFlxByWatFloInSno+HeatbyVapConvInSnow+HeatCndFlxInSno
        SnoX2SnoLay_snvr(L2,NY,NX)      = 0.0_r8
        WatX2SnoLay_snvr(L2,NY,NX)      = WatVapFloInSnow
        IceX2SnoLay_snvr(L2,NY,NX)      = 0.0_r8
        HeatX2SnoLay_snvr(L2,NY,NX)     = HeatByWatVapFloInSnow
        WatFlowInSnowM_snvr(M,L2,NY,NX) = WatFlowInSnowM_snvr(M,L2,NY,NX)+WatFloInSno

        if(abs(WatVapFloInSnow)>1.e10)call endrun(trim(mod_filename)//' at line',__LINE__)
        !
        ! DISCHARGE OF MELTWATER AND ITS HEAT FROM LOWEST SNOWPACK LAYER
        ! TO RESIDUE, SURFACE SOIL MICROPORES AND MACROPORES
      ELSE
        !dealing with residual snow
        !L==JS, or the layer L2 has insignificant snow, so layer L is directly interacting with litter and soil surface
        IF(.not.lchkBottomL)THEN
          !
          !interaction with respect to litter layer and topsoil
          !potential flow to soil
          !first  flow to micro- and macropores
          if(FracSurfByLitR_col(NY,NX)>ZEROL)then
            PtWatFlowSno2Soi      = WatFloInSnoMax*FracSurfBareSoil_col(NY,NX)
            WatFlowSno2MicP       = AMIN1(VLairMicP1_vr(NUM_col(NY,NX),NY,NX)*dts_wat,PtWatFlowSno2Soi*SoilFracAsMicP_vr(NUM_col(NY,NX),NY,NX))
            WatFlowSno2MacP       = AMIN1(VLairMacP1_vr(NUM_col(NY,NX),NY,NX)*dts_wat,PtWatFlowSno2Soi*SoilFracAsMacP1_vr(NUM_col(NY,NX),NY,NX))
            WatFlowSno2Soil       = WatFlowSno2MicP+WatFlowSno2MacP
            HeatFlowSno2SoiByWat  = cpw*TKSnow1_snvr(L,NY,NX)*WatFlowSno2Soil*HeatAdv_scal
            WatFlowSno2LitR       = WatFloInSnoMax-WatFlowSno2Soil
            HeatFlowSno2LitrByWat = cpw*TKSnow1_snvr(L,NY,NX)*WatFlowSno2LitR
          else
            PtWatFlowSno2Soi      = WatFloInSnoMax
            WatFlowSno2MicP       = PtWatFlowSno2Soi*SoilFracAsMicP_vr(NUM_col(NY,NX),NY,NX)
            WatFlowSno2MacP       = PtWatFlowSno2Soi*SoilFracAsMacP1_vr(NUM_col(NY,NX),NY,NX)
            WatFlowSno2Soil       = WatFlowSno2MicP+WatFlowSno2MacP
            HeatFlowSno2SoiByWat  = cpw*TKSnow1_snvr(L,NY,NX)*WatFlowSno2Soil*HeatAdv_scal
            !HeatFlowSno2SoiByWat  = cpw*(TKSnow1_snvr(L,NY,NX)-TKSoil1_vr(NUM_col(NY,NX),NY,NX))*WatFlowSno2Soil*HeatAdv_scal
            WatFlowSno2LitR       = 0._r8
            HeatFlowSno2LitrByWat = 0._r8
          endif
          !
          LiquidBudget(L)=MAX(0._r8,LiquidBudget(L)-WatFlowSno2Soil-WatFlowSno2LitR)
          call SnowTopSoilExch(dt_SnoHeat,M,L,NY,NX,VapCond1,VapSnoSrc,VLairSno1,TCND1W,&
            VapFlxSno2Soi1,HeatConvFlxSno2Soi1,HeatCndFlxSno2Soi,VapCond2,PSISV1,TCNDS,LiquidBudget(L))
          !  
          TotWatXFlx2SoiMicP = WatFlowSno2MicP+VapFlxSno2Soi1
          TotHeatFlow2Soi     = HeatFlowSno2SoiByWat+HeatConvFlxSno2Soi1+HeatCndFlxSno2Soi
          TotSnoWatFlow2Litr  = WatFlowSno2LitR
          TotSnoHeatFlow2Litr = HeatFlowSno2LitrByWat
          !
          !
          ! HEAT FLUX AMONG SNOWPACK, SURFACE RESIDUE AND SURFACE SOIL
          !
          !
          !surface litter layer is active
          IF(VHeatCapacity1_vr(0,NY,NX).GT.VHeatCapLitRMin_col(NY,NX) .and. FracSurfByLitR_col(NY,NX).GT.ZEROL)THEN
            !should the litter heat/water states be updated here?
            !
            call SnowSurLitterExch(I,J,dt_SnoHeat,M,L,NY,NX,VapCond1,VapCond2,TCND1W,VLairSno1,PSISV1,TCNDS,&
              TotWatXFlx2SoiMicP,TotHeatFlow2Soi,TotSnoWatFlow2Litr,TotSnoHeatFlow2Litr,LiquidBudget(L))

          ENDIF
          !
          ! GATHER WATER, VAPOR AND HEAT FLUXES INTO FLUX ARRAYS
          ! FOR LATER UPDATES TO STATE VARIABLES
          !
          CumWatFlx2SoiMicP  = CumWatFlx2SoiMicP+TotWatXFlx2SoiMicP
          CumWatXFlx2SoiMicP  = CumWatXFlx2SoiMicP+WatFlowSno2MicP
          CumWatFlx2SoiMacP   = CumWatFlx2SoiMacP+WatFlowSno2MacP
          !TotHeatFlow2Soi     = HeatFlowSno2SoiByWat+HeatConvFlxSno2Soi1+HeatCndFlxSno2Soi+cumHeatConvFlxLitr2Soi1+cumHeatCndFlxLitr2Soi
          
          !Heat flux diagnostics
          HeatFlowSno2SoiByWat_col(NY,NX) = HeatFlowSno2SoiByWat_col(NY,NX)+HeatFlowSno2SoiByWat
          !HAltFlowSno2SoiByWat_col(NY,NX) = HAltFlowSno2SoiByWat_col(NY,NX)+HAltFlowSno2SoiByWat
          HeatConvFlxSno2Soi_col(NY,NX) = HeatConvFlxSno2Soi_col(NY,NX)+HeatConvFlxSno2Soi1
          HeatCndFlxSno2Soi_col(NY,NX) = HeatCndFlxSno2Soi_col(NY,NX)+HeatCndFlxSno2Soi
          !HeatConvFlxLitr2Soi_col(NY,NX) = HeatConvFlxLitr2Soi_col(NY,NX)+cumHeatConvFlxLitr2Soi1
          !HeatCndFlxLitr2Soi_col(NY,NX) = HeatConvFlxLitr2Soi_col(NY,NX)+cumHeatCndFlxLitr2Soi
          
          !cumNetHeatFlow2Soil  = cumNetHeatFlow2Soil+TotHeatFlow2Soi
          !TotSnoWatFlow2Litr  = WatFlowSno2LitR+CumVapFlxSno2Litr-CumVapFlxLitr2Soi
          !TotSnoHeatFlow2Litr = HeatFlowSno2LitrByWat+CumHeatConvFlxSno2Litr+CumHeatCndFlxSno2Litr &
          !  -cumHeatConvFlxLitr2Soi1-cumHeatCndFlxLitr2Soi

          cumNetHeatFlow2Soil  = cumNetHeatFlow2Soil+TotHeatFlow2Soi
          CumSnowWatFLow2LitR = CumSnowWatFLow2LitR+TotSnoWatFlow2Litr
          CumNetHeatFlow2LitR = CumNetHeatFlow2LitR+TotSnoHeatFlow2Litr
          QSnowH2Oloss_col(NY,NX)   = QSnowH2Oloss_col(NY,NX)+TotSnoWatFlow2Litr+TotWatXFlx2SoiMicP+WatFlowSno2MacP

          WatFlowSno2LitRM_col(M,NY,NX) = WatFlowSno2LitRM_col(M,NY,NX)+WatFlowSno2LitR
          WatFlowSno2MicPM_col(M,NY,NX) = WatFlowSno2MicPM_col(M,NY,NX)+WatFlowSno2MicP
          WatFlowSno2MacPM_col(M,NY,NX) = WatFlowSno2MacPM_col(M,NY,NX)+WatFlowSno2MacP
          !do diagnoses
          WatConvSno2MicP_snvr(L,NY,NX)  = WatConvSno2MicP_snvr(L,NY,NX)+TotWatXFlx2SoiMicP
          WatConvSno2MacP_snvr(L,NY,NX)  = WatConvSno2MacP_snvr(L,NY,NX)+WatFlowSno2MacP
          HeatConvSno2Soi_snvr(L,NY,NX)  = HeatConvSno2Soi_snvr(L,NY,NX)+TotHeatFlow2Soi
          WatConvSno2LitR_snvr(L,NY,NX)  = WatConvSno2LitR_snvr(L,NY,NX)+TotSnoWatFlow2Litr
          HeatConvSno2LitR_snvr(L,NY,NX) = HeatConvSno2LitR_snvr(L,NY,NX)+TotSnoHeatFlow2Litr
          lchkBottomL                    = .true.
        ENDIF
      ENDIF
    ENDIF
  ENDDO D9880
  call PrintInfo('end '//subname)
  end subroutine SnowPackIterationMM

!------------------------------------------------------------------------------------------

  subroutine LimitSnowPhaseChange(RequestedHeat,DrySnowWE,IceVolume,LiquidVolume,PhaseHeat,SnowMeltWE,IceMeltWE)
  ! Description:
  ! Limit melting/freezing by positive donor inventories, including tiny films.
  ! Derive latent heat from the accepted mass so no heat is booked without its
  ! corresponding phase transfer. Negative ice melt denotes freezing.
  implicit none
  real(r8), intent(in) :: RequestedHeat !requested latent heat, positive for freezing [MJ/cell]
  real(r8), intent(in) :: DrySnowWE     !available dry snow [m3 water equivalent/cell]
  real(r8), intent(in) :: IceVolume     !available ice [m3 ice/cell]
  real(r8), intent(in) :: LiquidVolume  !available liquid [m3 water/cell]
  real(r8), intent(out) :: PhaseHeat   !accepted latent heat [MJ/cell]
  real(r8), intent(out) :: SnowMeltWE  !dry snow converted to liquid [m3 water/cell]
  real(r8), intent(out) :: IceMeltWE   !ice converted to liquid, negative for freezing [m3 water/cell]
  character(len=*), parameter :: subname='LimitSnowPhaseChange'
  real(r8) :: SolidWE,MeltWE,DryWE,FrozenWE

  call PrintInfo('beg '//subname)
  SnowMeltWE=0._r8
  IceMeltWE=0._r8
  if(RequestedHeat<0._r8)then
    DryWE=max(0._r8,DrySnowWE)
    FrozenWE=DENSICE*max(0._r8,IceVolume)
    SolidWE=DryWE+FrozenWE
    if(SolidWE>0._r8)then
      MeltWE=min(SolidWE,-RequestedHeat/LtHeatIceMelt)
      if(MeltWE==SolidWE)then
        SnowMeltWE=DryWE
        IceMeltWE=FrozenWE
      else
        SnowMeltWE=MeltWE*(DryWE/SolidWE)
        IceMeltWE=min(FrozenWE,max(0._r8,MeltWE-SnowMeltWE))
      endif
    endif
  else
    IceMeltWE=-min(max(0._r8,LiquidVolume),RequestedHeat/LtHeatIceMelt)
  endif
  PhaseHeat=-LtHeatIceMelt*(SnowMeltWE+IceMeltWE)
  call PrintInfo('end '//subname)
  end subroutine LimitSnowPhaseChange

!------------------------------------------------------------------------------------------
  subroutine SolveSnowpackM(I,J,M,NY,NX,LatentHeatAir2Sno,Radnet2Snow,HeatSensEvapAir2Snow,HeatSensAir2Snow,&
    HeatNetFlx2Snow,CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP,CumSnowWatFLow2LitR,&
    CumNetHeatFlow2LitR,cumNetHeatFlow2Soil)
  ! Description:
  ! Subcycle snow mass, heat and phase change with linearized longwave and
  ! sensible exchange; accumulate the same accepted fluxes for outer updates.
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M     !soil heat-flow iteration id
  integer, intent(in) :: NY,NX
  real(r8), intent(inout) :: HeatSensEvapAir2Snow
  real(r8), intent(out) :: CumSnowWatFLow2LitR     !cumulative snow water flow to litter
  real(r8), intent(out) :: HeatNetFlx2Snow
  real(r8), intent(out) :: LatentHeatAir2Sno
  real(r8), intent(out) :: Radnet2Snow
  real(r8), intent(out) :: HeatSensAir2Snow
  real(r8), intent(out) :: CumNetHeatFlow2LitR
  real(r8), intent(out) :: cumNetHeatFlow2Soil
  real(r8), intent(out) :: CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,CumWatXFlx2SoiMicP
  character(len=*), parameter :: subname='SolveSnowpackM'
  integer :: MM,L,L2

  real(r8) :: TotWatXFlx2SoiMicP,TotHeatFlow2Soi
  real(r8) :: ENGY0
  real(r8) :: TKApp,TKXR,TK1X
  real(r8) :: WatFlowSno2MacP,TotSnoWatFlow2Litr,TotSnoHeatFlow2Litr
  real(r8) :: NetIce2LayL,NetSno2LayL,NetWat2LayL,HeatByFrezThaw,TFLX1
  real(r8) :: NetHeat2LayL,VLHeatCapSnowMX
  real(r8) :: VOLI0X,VOLS0X,VOLW0X,IceThawMass,SnowThawMass
  real(r8) :: tNetWat2LayL,vwat,vdry,vice
  real(r8) :: dt_SnoHeat   !time step size for snow model iteration
  real(r8) :: cphwat,vhcp0
  real(r8) :: LongwaveDeriv,LongwaveCorrection
  real(r8) :: SnowSensibleDeriv,SensibleDeriv,BoundaryHeatDeriv !step-integrated coefficients [MJ/K]
  real(r8) :: SensibleCorrection                             !additional outgoing sensible heat [MJ/cell]
  real(r8) :: SnoFallWE,Rainfall,IceFallVol,HeatSnofall2Snow
  real(r8) :: TopLaySnoVol      !sum of constituent-volume weights for empirical albedo [m3]
  real(r8), parameter :: Kirchhoff_cf=1._r8+TFICE*(cpw-cps)/LtHeatIceMelt
  !     begin_execution
  !     HEAT AND VAPOR FLUXES BETWEEN SNOWPACK AND ATMOSPHERE
  !
  !     VHCPWM=volumetric heat capacity of snowpack
  !     NPS=number of cycles for solving snowpack heat and water fluxes
  !     SnowAlbedo=snowpack albedo
  !     VLDrySnoWE0M,VOLI0M,VOLW0M=snow,ice,water volumes
  !     RFLX0=net radiation input
  !     RadSWonSno=shortwave radiation at snowpack surface
  !     LWRad2Snow_col=longwave radn incident at snowpack surface
  !     LWRadSno1=longwave radn emitted by snowpack surface
  !     TKSnow1=snowpack surface temperature
  !     RadNet2Sno2=net radiation

  call PrintInfo('beg '//subname)
  CumWatFlx2SoiMacP   = 0._r8
  cumNetHeatFlow2Soil = 0._r8
  CumWatFlx2SoiMicP   = 0._r8
  CumNetHeatFlow2LitR = 0.0_r8
  CumWatXFlx2SoiMicP  = 0._r8
  HeatSensAir2Snow    = 0.0_r8
  Radnet2Snow         = 0.0_r8
  LatentHeatAir2Sno   = 0.0_r8
  HeatNetFlx2Snow     = 0.0_r8
  CumSnowWatFLow2LitR = 0._r8
  tNetWat2LayL        = 0._r8    !local diagnostics
  MMit                = 0
  dt_snoHeat          = dts_HeatWatTP*XNPS
  SnoFallWE           = SnowFallt_col(NY,NX)*XNPS
  Rainfall            = Rain2Snowt_col(NY,NX)*XNPS
  IceFallVol          = Ice2Snowt_col(NY,NX)*XNPS
  HeatSnofall2Snow    = PrecHeat2Snowt_col(NY,NX)*XNPS


  D3000: DO MM = 1, NPS

    SnowSensibleDeriv=0._r8
    TopLaySnoVol=(VLDrySnoWE0M_snvr(1,NY,NX)+SnoFallWE)+(VLIceSnow0M_snvr(1,NY,NX)+IceFallVol) &
      +VLWatSnow0M_snvr(1,NY,NX)+Rainfall

    if(TopLaySnoVol.GT.0._r8)then

      call SnowAtmosExchangeMM(I,J,M,NY,NX,SnoFallWE,Rainfall,IceFallVol,TopLaySnoVol,HeatSnofall2Snow,LatentHeatAir2Sno,&
        HeatSensEvapAir2Snow,HeatNetFlx2Snow,Radnet2Snow,HeatSensAir2Snow,SnowSensibleDeriv)

      call SnowPackIterationMM(dt_snoHeat,I,J,M,NY,NX,TotWatXFlx2SoiMicP,TotHeatFlow2Soi,WatFlowSno2MacP,&
        TotSnoWatFlow2Litr,TotSnoHeatFlow2Litr,CumWatFlx2SoiMacP,CumWatFlx2SoiMicP,&
        CumWatXFlx2SoiMicP,CumSnowWatFLow2LitR,CumNetHeatFlow2LitR,cumNetHeatFlow2Soil)
    endif
!
!     ACCUMULATE SNOWPACK FLUXES TO LONGER TIME STEP FOR
!     LITTER, SOIL FLUX CALCULATIONS
!
!     LWRadBySurf_col=total longwave emission
!     XFLWS,WatXfer2SnoLay,IceXfer2SnoLay=hourly accumulated snow,water,ice transfer
!     HeatXfer2SnoLay=hourly convective heat flux from snow,water,ice transfer
!     NetSno2LayL,NetWat2LayL,NetIce2LayL=net snow,water,ice transfer
!     NetHeat2LayL=convective heat flux from net snow,water,ice transfer
!     TFLWS,TFLWW,TFLWI=accumulated net snow,water,ice transfer
!     CumHeat2SnowLay=convective heat flux from accumd snow,water,ice transfer
!
   D9860: DO L=1,JS

      VOLS0X = VLDrySnoWE0M_snvr(L,NY,NX)
      VOLW0X = VLWatSnow0M_snvr(L,NY,NX)
      VOLI0X = VLIceSnow0M_snvr(L,NY,NX)

      !next layer
      L2 = MIN(JS,L+1)
!
!     IF WITHIN SNOWPACK
!
      IF(L.LT.JS .AND. VLHeatCapSnowM1_snvr(L2,NY,NX).GT.ZEROS(NY,NX))THEN
        NetSno2LayL  = SnoX2SnoLay_snvr(L,NY,NX)-SnoX2SnoLay_snvr(L2,NY,NX)
        NetWat2LayL  = WatX2SnoLay_snvr(L,NY,NX)-WatX2SnoLay_snvr(L2,NY,NX)
        NetIce2LayL  = IceX2SnoLay_snvr(L,NY,NX)-IceX2SnoLay_snvr(L2,NY,NX)
        NetHeat2LayL = HeatX2SnoLay_snvr(L,NY,NX)-HeatX2SnoLay_snvr(L2,NY,NX)
        !
        ! Donor limits were applied to the shared faces before accumulation.
        !layer L2 has insignificant snow, or L==JS, and layer L has significant snow
      ELSEIF(VLHeatCapSnowM1_snvr(L,NY,NX).GT.ZEROS(NY,NX))THEN

        NetSno2LayL  = SnoX2SnoLay_snvr(L,NY,NX)
        NetWat2LayL  = WatX2SnoLay_snvr(L,NY,NX)-TotSnoWatFlow2Litr-TotWatXFlx2SoiMicP-WatFlowSno2MacP
        NetIce2LayL  = IceX2SnoLay_snvr(L,NY,NX)
        NetHeat2LayL = HeatX2SnoLay_snvr(L,NY,NX)-TotSnoHeatFlow2Litr-TotHeatFlow2Soi

      ELSE
        NetSno2LayL  = 0.0_r8
        NetWat2LayL  = 0.0_r8
        NetIce2LayL  = 0.0_r8
        NetHeat2LayL = 0.0_r8
      ENDIF

      ! Phase change uses the material remaining after all accepted transfers.
      VOLS0X = VOLS0X+NetSno2LayL
      VOLW0X = VOLW0X+NetWat2LayL
      VOLI0X = VOLI0X+NetIce2LayL
      SnoXfer2SnoLay_snvr(L,NY,NX)  = SnoXfer2SnoLay_snvr(L,NY,NX)+SnoX2SnoLay_snvr(L,NY,NX)
      WatXfer2SnoLay_snvr(L,NY,NX)  = WatXfer2SnoLay_snvr(L,NY,NX)+WatX2SnoLay_snvr(L,NY,NX)
      IceXfer2SnoLay_snvr(L,NY,NX)  = IceXfer2SnoLay_snvr(L,NY,NX)+IceX2SnoLay_snvr(L,NY,NX)
      HeatXfer2SnoLay_snvr(L,NY,NX) = HeatXfer2SnoLay_snvr(L,NY,NX)+HeatX2SnoLay_snvr(L,NY,NX)
!
!     FREEZE-THAW IN SNOWPACK FROM NET CHANGE IN SNOWPACK
!     HEAT STORAGE
!
!     VLDrySnoWE0M,VOLW0M,VOLI0M=snow,water,ice volume
!     VHCPWMM,VLHeatCapSnowMX,VLHeatCapSnowMin_col=previous,current,minimum heat capacity
!     NetHeat2LayL=net conductive+convective heat flux
!     TFLX1=unconstrained latent heat flux from freeze-thaw
!     HeatByFrezThaw=source-limited latent heat flux from freeze-thaw
!     SnowThawMass,IceThawMass=freeze-thaw changes in water,ice
!     SnowThawMassL,IceThawMassL=accumulated freeze-thaw
!     PhaseChangeHeatL=accumulated latent heat flux from freeze-thaw
!     XSnowThawMassL,XIceThawMassL=hourly accumulated freeze-thaw
!     XPhaseChangeHeatL=hourly accumulated latent heat flux from freeze-thaw
!
      ENGY0           = VLHeatCapSnowM1_snvr(L,NY,NX)*TKSnow1_snvr(L,NY,NX)
      VLHeatCapSnowMX = cps*VOLS0X+cpw*VOLW0X+cpi*VOLI0X

      ! Freeze resistance and local air during this substep. Both boundary
      ! heat losses respond to the new top-snow temperature [MJ/K].
      SensibleDeriv=0._r8
      if(L==1)SensibleDeriv=SnowSensibleDeriv
      LongwaveDeriv=0._r8
      if(L==1 .and. VLHeatCapSnowMX>ZEROS(NY,NX))then
        LongwaveDeriv=4._r8*LWEmscefSnow_col(NY,NX)*TKSnow1_snvr(L,NY,NX)**3/real(NPS,kind=r8)
      endif
      BoundaryHeatDeriv=LongwaveDeriv+SensibleDeriv

      IF(VLHeatCapSnowMX.GT.ZEROS(NY,NX))THEN
        !apparent temperature before freeze-thaw
        TKApp=(ENGY0+NetHeat2LayL+BoundaryHeatDeriv*TKSnow1_snvr(L,NY,NX))/(VLHeatCapSnowMX+BoundaryHeatDeriv)

        IF((TKApp.LT.TFice .AND. VOLW0X.GT.ZERO*VcumSnoDWI_col(NY,NX)) &
          .OR.(TKApp.GT.TFice .AND. VOLI0X+VOLS0X.GT.ZERO*VcumSnoDWI_col(NY,NX)))THEN
          !freeze-thaw condition met
          TFLX1=(VLHeatCapSnowMX+BoundaryHeatDeriv)*(TFice-TKApp)/Kirchhoff_cf*dts_wat

          call LimitSnowPhaseChange(TFLX1,VOLS0X,VOLI0X,VOLW0X,HeatByFrezThaw,SnowThawMass,IceThawMass)
        ELSE
          TFLX1          = 0.0_r8
          HeatByFrezThaw = 0.0_r8
          SnowThawMass   = 0.0_r8
          IceThawMass    = 0.0_r8
        ENDIF
      ELSE
        HeatByFrezThaw = 0.0_r8
        SnowThawMass   = 0.0_r8
        IceThawMass    = 0.0_r8
      ENDIF
      !
      !     INTERNAL SNOWPACK SNOW, WATER, ICE, TEMPERATURE
      !
      !     VLDrySnoWE0M,VOLW0M,VOLI0M=snow water eqv,water,ice volume
      !     NetSno2LayL,NetWat2LayL,NetIce2LayL=net snow,water,ice transfer
      !     NetHeat2LayL=conductive+convective heat from snow,water,ice transfer
      !     SnowThawMass,IceThawMass=freeze-thaw changes in water,ice
      !     HeatByFrezThaw=source-limited latent heat flux from freeze-thaw
      !     DENSICE=ice density
      !     TKSnow1,TairK=snowpack,air temperature
      !     VHCPWMM,VLHeatCapSnowMin_col=snowpack, minimum heat capacity
      !
      if(abs(NetWat2LayL)>1.e10)call endrun(trim(mod_filename)//' at line',__LINE__)

      vdry                       = VLDrySnoWE0M_snvr(L,NY,NX)
      vwat                       = VLWatSnow0M_snvr(L,NY,NX)
      vice                       = VLIceSnow0M_snvr(L,NY,NX)
      ! Retain positive residuals instead of deleting them at a size cutoff.
      VLDrySnoWE0M_snvr(L,NY,NX) = MAX(0._r8,VOLS0X-SnowThawMass)
      VLWatSnow0M_snvr(L,NY,NX)  = MAX(0._r8,VOLW0X+SnowThawMass+IceThawMass)
      VLIceSnow0M_snvr(L,NY,NX)  = MAX(0._r8,VOLI0X-IceThawMass/DENSICE)

      cphwat                        = cpw*VLWatSnow0M_snvr(L,NY,NX)
      tNetWat2LayL                  = tNetWat2LayL+NetWat2LayL
      ENGY0                         = VLHeatCapSnowM1_snvr(L,NY,NX)*TKSnow1_snvr(L,NY,NX)
      vhcp0                         = VLHeatCapSnowM1_snvr(L,NY,NX)
      VLHeatCapSnowM1_snvr(L,NY,NX) = cps*VLDrySnoWE0M_snvr(L,NY,NX)+cpw*VLWatSnow0M_snvr(L,NY,NX)+cpi*VLIceSnow0M_snvr(L,NY,NX)
      TK1X                          = TKSnow1_snvr(L,NY,NX)

      ! A depleted layer can still finish its accepted boundary exchange.
      ! Keeping the boundary derivative closes energy even when C_new is zero.
      IF(VLHeatCapSnowM1_snvr(L,NY,NX)+BoundaryHeatDeriv.GT.0._r8)THEN
        TKSnow1_snvr(L,NY,NX)=(ENGY0+NetHeat2LayL+HeatByFrezThaw+BoundaryHeatDeriv*TK1X) &
          /(VLHeatCapSnowM1_snvr(L,NY,NX)+BoundaryHeatDeriv)
      ELSEIF(L.EQ.1)THEN
        TKSnow1_snvr(L,NY,NX)=TairK_col(NY,NX)
      ELSE
        TKSnow1_snvr(L,NY,NX)=TKSnow1_snvr(L-1,NY,NX)
      ENDIF
      ! Record the same implicit radiation exchange used by the temperature solve.
      ! Positive correction is additional outgoing longwave; net heat loses it.
      if(L==1)then
        LongwaveCorrection=LongwaveDeriv*(TKSnow1_snvr(L,NY,NX)-TK1X)
        NetHeat2LayL=NetHeat2LayL-LongwaveCorrection
        HeatX2SnoLay_snvr(1,NY,NX)=HeatX2SnoLay_snvr(1,NY,NX)-LongwaveCorrection
        HeatXfer2SnoLay_snvr(1,NY,NX)=HeatXfer2SnoLay_snvr(1,NY,NX)-LongwaveCorrection
        LWRadBySurf_col(NY,NX)=LWRadBySurf_col(NY,NX)+LongwaveCorrection
        Radnet2Snow=Radnet2Snow-LongwaveCorrection
        HeatNetFlx2Snow=HeatNetFlx2Snow-LongwaveCorrection
        ! Replace K_H*(Tq-Told) by K_H*(Tq-Tnew) in every accepted heat
        ! ledger. The outer reconstruction consumes these corrected fluxes.
        SensibleCorrection=SensibleDeriv*(TKSnow1_snvr(L,NY,NX)-TK1X)
        NetHeat2LayL=NetHeat2LayL-SensibleCorrection
        HeatX2SnoLay_snvr(1,NY,NX)=HeatX2SnoLay_snvr(1,NY,NX)-SensibleCorrection
        HeatXfer2SnoLay_snvr(1,NY,NX)=HeatXfer2SnoLay_snvr(1,NY,NX)-SensibleCorrection
        HeatSensAir2Snow=HeatSensAir2Snow-SensibleCorrection
        HeatNetFlx2Snow=HeatNetFlx2Snow-SensibleCorrection
      endif
      CumSno2SnowLM_snvr(L,NY,NX)      = CumSno2SnowLM_snvr(L,NY,NX)+NetSno2LayL
      CumWat2SnowLM_snvr(L,NY,NX)      = CumWat2SnowLM_snvr(L,NY,NX)+NetWat2LayL
      CumIce2SnowLM_snvr(L,NY,NX)      = CumIce2SnowLM_snvr(L,NY,NX)+NetIce2LayL
      XSnowThawMassLM_snvr(L,NY,NX)    = XSnowThawMassLM_snvr(L,NY,NX)+SnowThawMass
      XIceThawMassLM_snvr(L,NY,NX)     = XIceThawMassLM_snvr(L,NY,NX)+IceThawMass
      CumHeat2SnowLM_snvr(L,NY,NX)     = CumHeat2SnowLM_snvr(L,NY,NX)+NetHeat2LayL
      XPhaseChangeHeatLM_snvr(L,NY,NX) = XPhaseChangeHeatLM_snvr(L,NY,NX)+HeatByFrezThaw
      tEnGYM_snvr(L,NY,NX)=ENGY0+NetHeat2LayL+HeatByFrezThaw

      if(VLHeatCapSnowM1_snvr(L,NY,NX).GT.ZEROS(NY,NX) .and. &
        (TK1X/=spval .and. abs(TK1X-TKSnow1_snvr(L,NY,NX))>20._r8 .or. TKSnow1_snvr(L,NY,NX)<200._r8))then
        write(iulog,*)'dh, TK1X, TKSnow1',I+J/24.,M,L,TK1X,TKSnow1_snvr(L,NY,NX)
        write(iulog,*)'VH, VHMX, tair,TKSnow0',VLHeatCapSnowM1_snvr(L,NY,NX), VLHeatCapSnowMX,TairK_col(NY,NX),TKSnow_snvr(L,NY,NX)
        write(iulog,*)'TKS_1, TKS_0',TKS_vr(NUM_col(NY,NX),NY,NX),TKS_vr(0,NY,NX),ZEROS(NY,NX)
        write(iulog,*)'energy(ENGY0, NetHeat2LayL, HeatByFT ', ENGY0, NetHeat2LayL, HeatByFrezThaw
        write(iulog,*)'mass0 (dry, wat, ice): ',vdry,vwat,vice
        write(iulog,*)'masst (dry, wat, ice): ',VLDrySnoWE0M_snvr(L,NY,NX),VLWatSnow0M_snvr(L,NY,NX),VLIceSnow0M_snvr(L,NY,NX)
        write(iulog,*)'NYNX',NY,NX
        call endrun('crazy snow temperature '//trim(mod_filename),__LINE__)
      endif

!      if(M.EQ.1.AND.MM.EQ.1.AND.L.EQ.1)then
!         write(99,*) TKSnow1_snvr(1,NY,NX), " " , ENGY0, " " , NetHeat2LayL, " " , HeatByFrezThaw
!      endif

      if(cphwat < 0._r8)then
        if(abs(cphwat/VLHeatCapSnowM1_snvr(L,NY,NX))>1.e-3_r8)then
          write(*,*)'too great negative water content',VLWatSnow0M_snvr(L,NY,NX),&
            VLDrySnoWE0M_snvr(L,NY,NX),VLIceSnow0M_snvr(L,NY,NX)
          write(*,*)vwat,vdry,vice,cphwat/VLHeatCapSnowM1_snvr(L,NY,NX)
          call endrun('negative snow water',__LINE__)
        endif
      endif

    ENDDO D9860
  ENDDO D3000

  !if(M.EQ.1)then
  !  write(99,*) TKSnow1_snvr(1,NY,NX), " " , ENGY0, " " , NetHeat2LayL, " " , HeatByFrezThaw
  !endif

  call PrintInfo('end '//subname)
 end subroutine SolveSnowpackM
!------------------------------------------------------------------------------------------
  subroutine SnowAtmosExchangeMM(I,J,M,NY,NX,SnoFallWE,Rainfall,IceFallVol,TopLaySnoVol,HeatSnofall2Snow,&
    LatentHeatAir2Sno,HeatSensEvapAir2Snow,HeatNetFlx2Snow,Radnet2Snow,HeatSensAir2Snow,SnowSensibleDeriv)
  ! Description:
  ! Evaluate snow-atmosphere exchange at the old state and return the frozen
  ! sensible coefficient for the semi-implicit snow-temperature correction.
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M    !soil heat-flow iteration id
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: HeatSnofall2Snow,IceFallVol,SnoFallWE,Rainfall
  real(r8),intent(in) :: TopLaySnoVol
  real(r8), intent(inout) :: LatentHeatAir2Sno
  real(r8), intent(inout) :: HeatSensEvapAir2Snow    !cumulated heat by vapor advection from air to snow [MJ]
  real(r8), intent(inout) :: HeatNetFlx2Snow,Radnet2Snow
  real(r8), intent(inout) :: HeatSensAir2Snow
  real(r8), intent(out) :: SnowSensibleDeriv !area/time-scaled sensible coefficient for this inner step [MJ/K]
  character(len=*), parameter :: subname='SnowAtmosExchangeMM'
  real(r8) :: SnowAlbedo    !snow surface albedo for reflecting short wave radiation
  real(r8) :: RFLX0         !net radiation (from sun + canopy) on snow surface [MJ]
  real(r8) :: RI            !Richardson number
  real(r8) :: LWRadSno1     !snow emitted long wave radiation, [MJ]
  real(r8) :: RadNet2Sno2   !net radiation on snow surface [MJ]
  real(r8) :: RAGX
  real(r8):: Raa
  real(r8) :: CdSnoEvap        !time scaled conductance for snowpack latent heat flux, [m2 h]/[h/m] =[m3],
  real(r8) :: CdSnoHSens       !area/time-scaled snow sensible coefficient [MJ/K]
  real(r8) :: VPSno0           !vapor pressure, kg H2O/m3
  real(r8) :: EVAPW2           !snow loss as evaporation [kg H2O]
  real(r8) :: EVAPX2
  real(r8) :: HeatSensAir2Sno2      !sensible heat flux [MJ]
  real(r8) :: LatentHeatAir2Sno2    !latent heat flux [MJ]
  real(r8) :: EvapSublimation2                           !snow loss as sublimation
  real(r8) :: MaxVapXAir2Sno                             ![kg H2O]
  real(r8) :: HeatNetFlx2Sno1,HeatNetFlx2Sno2
  real(r8) :: HeatAdvAir2SnoByEvap2                      !convective heat flux
  real(r8) :: NetHeatAir2Snow,SnofallRain
  real(r8) :: SnofallDry,Snofallice
  real(r8) :: RadSWbySnow                                !shortwave radiation absorbed by snow [MJ]


  call PrintInfo('beg '//subname)

  if(mod_snow_albedo)then
    SnowAlbedo = SnowAlbedo_col(NY,NX)
  else
    !note TopLaySnoVol below is neither a water equivalent volume nor geometric volume, it is only for empirical
    !calculation of surface albedo. The albedo calculation can be done better
    SnowAlbedo=(0.85_r8*(VLDrySnoWE0M_snvr(1,NY,NX)+SnoFallWE)+0.30_r8*(VLIceSnow0M_snvr(1,NY,NX)+IceFallVol) &
    +0.06_r8*(VLWatSnow0M_snvr(1,NY,NX)+Rainfall))/TopLaySnoVol
  endif

  RadSWbySnow          = (1.0_r8-SnowAlbedo)*RadSW2Sno_col(NY,NX)
  RFLX0                = RadSWbySnow+LWRad2Snow_col(NY,NX)                            !incoming radiation, short + longwave
  LWRadSno1            = LWEmscefSnow_col(NY,NX)*TKSnow1_snvr(1,NY,NX)**4/real(NPS,kind=r8)         !emitting longwave radiation,
  RadNet2Sno2          = RFLX0-LWRadSno1                                              !net radiation
  Eco_RadSW_col(NY,NX) = Eco_RadSW_col(NY,NX) + RadSWbySnow

  !
  !     AERODYNAMIC RESISTANCE ABOVE SNOWPACK INCLUDING
  !     RESISTANCE IMPOSED BY PLANT CANOPY
  !
  !     RI=Richardsons number
  !     RIB=isothermal RI
  !     TKQ=canopy air temperature
  !     RAGX,RA=snowpack blr, h/m
  !     RAG,RAGW=isothermal blrs at ground,snowpack surfaces
  !
  RI   = RichardsonNumber(RIB_col(NY,NX),TKQ_col(NY,NX),TKSnow1_snvr(1,NY,NX))
  RAGX = AMAX1(RAM,0.8_r8*RawH2OSnow2CanAir_col(NY,NX), &
    AMIN1(1.2_r8*RawH2OSnow2CanAir_col(NY,NX),RawIsoTSoil2CanAir_col(NY,NX)/(1.0_r8-10.0_r8*RI)))
  RawH2OSnow2CanAir_col(NY,NX) = RAGX
  RAa = SurfaceExchangeResistance(NY,NX,TKSnow1_snvr(1,NY,NX),RAGX)
  !
  ! PARAMETERS FOR CALCULATING LATENT AND SENSIBLE HEAT FLUXES
  !
  !     CdSnoEvap,CdSnoHSens=conductance for snowpack latent,sensible heat fluxes
  !     PAREW,PARSW=conductances for latent,sensible heat fluxes
  !     RZ=surface resistance
  !     VPSno0,VPQ_col=vapor pressure at snowpack surface, canopy air, ton H2O/m3
  !     MaxVapXAir2Sno,EVAPW2,EvapSublimation2=evaporation total, water,snow
  !     XNPS=1/NPS
  !     LatentHeatAir2Sno2=latent heat flux
  !     VAP,VAPS=latent heat of evaporation,sublimation
  !     HeatAdvAir2SnoByEvap2=convective heat of evaporation flux
  !
  CdSnoEvap      = AScaledCdWOverSnow_col(NY,NX)/(RAa+RZ)
  CdSnoHSens     = AScaledCdHOverSnow_col(NY,NX)/RAa
  SnowSensibleDeriv=CdSnoHSens
  VPSno0         = vapsat(TKSnow1_snvr(1,NY,NX))  !*dssign(VLWatSnow0M_snvr(1,NY,NX))
  MaxVapXAir2Sno = CdSnoEvap*(VPQ_col(NY,NX)-VPSno0)

  !first: evaporation from snow held water
  EVAPW2 = AMAX1(MaxVapXAir2Sno,-AZMAX1(VLWatSnow0M_snvr(1,NY,NX)*dts_sno))
  !second: sublimation from dry snow
  EVAPX2 = AZMIN1(MaxVapXAir2Sno-EVAPW2)
  !then the loss is from dry snow
  EvapSublimation2   = AMAX1(EVAPX2,-AZMAX1(VLDrySnoWE0M_snvr(1,NY,NX)*dts_sno))
  LatentHeatAir2Sno2 = EVAPW2*EvapLHTC+EvapSublimation2*SublmHTC
  ! Differentiate the sequential liquid/snow donor limits at fixed stores.
  ! Match the latent-energy ledger, including the sublimation latent heat.
  if(MaxVapXAir2Sno>-AZMAX1(VLWatSnow0M_snvr(1,NY,NX)*dts_sno))then
    SurfLatHeatVapDeriv_col(NY,NX)=SurfLatHeatVapDeriv_col(NY,NX)+EvapLHTC*CdSnoEvap
  elseif(MaxVapXAir2Sno-EVAPW2>-AZMAX1(VLDrySnoWE0M_snvr(1,NY,NX)*dts_sno))then
    SurfLatHeatVapDeriv_col(NY,NX)=SurfLatHeatVapDeriv_col(NY,NX)+SublmHTC*CdSnoEvap
  endif

  IF(MaxVapXAir2Sno.LT.0.0_r8)THEN
    !snow is losing water/heat
    HeatAdvAir2SnoByEvap2=(EVAPW2*cpw+EvapSublimation2*cps)*TKSnow1_snvr(1,NY,NX)
  ELSE
    !snow is gaining water/heat, condensation/deposition
    HeatAdvAir2SnoByEvap2=(EVAPW2*cpw+EvapSublimation2*cps)*TKQ_col(NY,NX)
  ENDIF
!
!     SOLVE FOR SNOWPACK SURFACE TEMPERATURE AT WHICH ENERGY
!     BALANCE OCCURS, SOLVE AND ACCUMULATE LATENT, SENSIBLE
!     STORAGE HEAT FLUXES AND EVAPORATION
!
!     HeatSensAir2Sno2,LatentHeatAir2Sno2,RadNet2Sno2=sensible,latent heat fluxes, net radiation
!     HeatAdvAir2SnoByEvap2=convective heat flux from LatentHeatAir2Sno2
!     HeatNetFlx2Sno1=storage heat flux
!     SnoFall,Rainfall,IceFall=snow,water,ice input to snowpack
!     HeatSnofall2Snow=convective heat from snow,water,ice input to snowpack
!

  HeatSensAir2Sno2 = CdSnoHSens*(TKQ_col(NY,NX)-TKSnow1_snvr(1,NY,NX))
  ! The hourly air predictor holds accepted surface temperatures fixed.
  ! Its partial derivative remains K_H after the implicit flux correction.
  SurfSensHeatConduct_col(NY,NX)=SurfSensHeatConduct_col(NY,NX)+CdSnoHSens
  !occasionally, RadNet2Sno2 and HeatSensAir2Sno2 go to infinity
  HeatNetFlx2Sno1      = RadNet2Sno2+LatentHeatAir2Sno2+HeatSensAir2Sno2

  HeatNetFlx2Sno2      = HeatNetFlx2Sno1+HeatAdvAir2SnoByEvap2
  Radnet2Snow          = Radnet2Snow+RadNet2Sno2
  LatentHeatAir2Sno    = LatentHeatAir2Sno+LatentHeatAir2Sno2
  HeatSensEvapAir2Snow = HeatSensEvapAir2Snow+HeatAdvAir2SnoByEvap2
  HeatSensAir2Snow     = HeatSensAir2Snow+HeatSensAir2Sno2
  HeatNetFlx2Snow      = HeatNetFlx2Snow+HeatNetFlx2Sno2

  EVAPS_col(NY,NX)       = EVAPS_col(NY,NX)+EvapSublimation2
  EVAPW_col(NY,NX)       = EVAPW_col(NY,NX)+EVAPW2
  VapXAir2Sno_col(NY,NX) = VapXAir2Sno_col(NY,NX)+EvapSublimation2+EVAPW2
  SnofallDry             = SnoFallWE+EvapSublimation2
  SnofallRain            = Rainfall+EVAPW2
  Snofallice             = IceFallVol
  NetHeatAir2Snow        = HeatSnofall2Snow+HeatNetFlx2Sno2

  SnoX2SnoLay_snvr(1,NY,NX)      = SnofallDry
  WatX2SnoLay_snvr(1,NY,NX)      = SnofallRain
  IceX2SnoLay_snvr(1,NY,NX)      = Snofallice
  HeatX2SnoLay_snvr(1,NY,NX)     = NetHeatAir2Snow
  WatFlowInSnowM_snvr(M,1,NY,NX) = WatFlowInSnowM_snvr(M,1,NY,NX)+SnoFallWE+Rainfall+IceFallVol

  Prec2Snow_col(NY,NX)     = Prec2Snow_col(NY,NX)+SnoFallWE+Rainfall+IceFallVol
  PrecHeat2Snow_col(NY,NX) = PrecHeat2Snow_col(NY,NX)+HeatSnofall2Snow
  ! Record the actual nested-substep input, also when snow later disappears.
  PrecHeat_col(NY,NX) = PrecHeat_col(NY,NX)+HeatSnofall2Snow
  QSnowH2Oloss_col(NY,NX)  = QSnowH2Oloss_col(NY,NX)-EvapSublimation2-EVAPW2
  RainPrec2Sno_col(NY,NX)  = RainPrec2Sno_col(NY,NX)+Rainfall

  if(WatFlowInSnowM_snvr(M,1,NY,NX)>0._r8 .and. isclose(TCSnow_snvr(1,NY,NX),spval))then
    TCSnow_snvr(1,NY,NX)=units%Kelvin2Celcius(TairK_col(NY,NX))
  endif

  LWRadBySurf_col(NY,NX)=LWRadBySurf_col(NY,NX)+LWRadSno1
  call PrintInfo('end '//subname)
  end subroutine SnowAtmosExchangeMM
!------------------------------------------------------------------------------------------
  subroutine PrepIterSnowLayerM(I,J,M,NY,NX)
  !
  !called at the begining of each soil heat-flow iteration
  implicit none
  integer, intent(in) :: M   !soil heat-flow iteration id
  integer, intent(in) :: NY,NX,I,J
  character(len=*), parameter :: subname='PrepIterSnowLayerM'
  integer :: L

  ! FLUX VARIABLES IN SNOWPACK
  !
  ! PhaseChangeHeatL=latent heat from freeze-thaw
  ! SnowThawMassL,IceThawMassL=freeze-thaw between snow,ice and water
  ! FLW0S,FLW0W,FLW0I=snow,water,ice fluxes
  ! HFLW0W=convective heat flux from snow,water,ice fluxes
  ! WatFlowInSnowM=snowpack water flux
  ! VLDrySnoWE0M,VOLW0M,VOLI0M=snow,water,ice contents
  ! VHCPWMM=volumetric heat capacity
  ! TKSnow1=snow temperature
  !
  call PrintInfo('beg '//subname)
  D9765: DO L=1,JS
    VLDrySnoWE0M_snvr(L,NY,NX)    = VLDrySnoWE0_snvr(L,NY,NX)
    VLWatSnow0M_snvr(L,NY,NX)     = VLWatSnow0_snvr(L,NY,NX)
    VLIceSnow0M_snvr(L,NY,NX)     = VLIceSnow0_snvr(L,NY,NX)
    VLHeatCapSnowM1_snvr(L,NY,NX) = cps*VLDrySnoWE0_snvr(L,NY,NX)+cpw*VLWatSnow0_snvr(L,NY,NX)+cpi*VLIceSnow0_snvr(L,NY,NX)
    TKSnow1_snvr(L,NY,NX)         = TKSnow0_snvr(L,NY,NX)
  ENDDO D9765
  call PrintInfo('end '//subname)
  end subroutine PrepIterSnowLayerM
!------------------------------------------------------------------------------------------

  subroutine UpdateSnowPack1M(I,J,M,NY,NX,RadHeatCorrection)
  ! Reconstruct the outer-step snow state from accepted mass and energy,
  ! retaining thin stores and returning the accepted longwave correction.
  implicit NONE
  integer, intent(in) :: I,J,M,NY,NX
  real(r8), intent(out) :: RadHeatCorrection !additional net radiation into snow [MJ/cell/step]

  integer :: L
  real(r8) :: ENGY0,vdry,vwat,vice
  real(r8) :: dHPhaseChange,VLDrySnoWEtmp
  real(r8) :: TKX,tENGY,dENGY,LongwaveStabilizer
  logical :: ActiveSnow_test
!
!     VOLS0,VOLW0,VOLI0=snow,water,ice volumes in snowpack
!     TFLWS,TFLWW=net snow,water flux
!     SnowThawMassL,IceThawMassL=snow-water,ice-water freeze-thaw flux
!     DENSICE=ice density
!     VHCPWM=snowpack volumetric heat capacity
!     TK0=snowpack temperature
!     CumHeat2SnowLay=total snowpack conductive+convective heat flux
!     PhaseChangeHeatL=snowpack latent heat flux from freeze-thaw
!

  RadHeatCorrection=0._r8

  D9780: DO L=1,JS
    vdry=VLDrySnoWE0_snvr(L,NY,NX); vwat=VLWatSnow0_snvr(L,NY,NX); vice=VLIceSnow0_snvr(L,NY,NX)
    ENGY0                     = VLSnowHeatCapM_snvr(M,L,NY,NX)*TKSnow0_snvr(L,NY,NX)
    TKX                       = TKSnow0_snvr(L,NY,NX)
    VLDrySnoWE0_snvr(L,NY,NX) = MAX(0._r8,VLDrySnoWE0_snvr(L,NY,NX)+CumSno2SnowLM_snvr(L,NY,NX)-XSnowThawMassLM_snvr(L,NY,NX))
    VLIceSnow0_snvr(L,NY,NX)  = MAX(0._r8,VLIceSnow0_snvr(L,NY,NX)+CumIce2SnowLM_snvr(L,NY,NX)-XIceThawMassLM_snvr(L,NY,NX)/DENSICE)
    VLWatSnow0_snvr(L,NY,NX)  = MAX(0._r8,VLWatSnow0_snvr(L,NY,NX)+CumWat2SnowLM_snvr(L,NY,NX) &
      +XSnowThawMassLM_snvr(L,NY,NX)+XIceThawMassLM_snvr(L,NY,NX))

    VLSnowHeatCapM_snvr(M+1,L,NY,NX) = cps*VLDrySnoWE0_snvr(L,NY,NX)+cpw*VLWatSnow0_snvr(L,NY,NX)+cpi*VLIceSnow0_snvr(L,NY,NX)

    LongwaveStabilizer=0._r8
    if(L==1)LongwaveStabilizer=4._r8*LWEmscefSnow_col(NY,NX)*TKX**3
    ActiveSnow_test = VLSnowHeatCapM_snvr(M+1,L,NY,NX)+LongwaveStabilizer.GT.0._r8
    IF(ActiveSnow_test)THEN
      ! Retain the historical outer longwave stabilizer as additional damping
      ! of the full-step reconstruction; the inner snow solve is also implicit.
      TKSnow0_snvr(L,NY,NX) = (ENGY0+CumHeat2SnowLM_snvr(L,NY,NX)+XPhaseChangeHeatLM_snvr(L,NY,NX) &
        +LongwaveStabilizer*TKX)/(VLSnowHeatCapM_snvr(M+1,L,NY,NX)+LongwaveStabilizer)
      if(L==1)then
        ! Net incoming correction is negative for warming [MJ/cell]. Return it
        ! before the caller accumulates Radnet2Snow; outgoing LW is stored here.
        RadHeatCorrection=-LongwaveStabilizer*(TKSnow0_snvr(L,NY,NX)-TKX)
        CumHeat2SnowLM_snvr(L,NY,NX)=CumHeat2SnowLM_snvr(L,NY,NX)+RadHeatCorrection
        HeatXfer2SnoLay_snvr(L,NY,NX)=HeatXfer2SnoLay_snvr(L,NY,NX)+RadHeatCorrection
        tEnGYM_snvr(L,NY,NX)=tEnGYM_snvr(L,NY,NX)+RadHeatCorrection
        LWRadBySurf_col(NY,NX)=LWRadBySurf_col(NY,NX)-RadHeatCorrection
      endif
      
      ! Only an empty store without boundary exchange needs a placeholder.
    ELSEIF(L.EQ.1)THEN
      TKSnow0_snvr(L,NY,NX)=TairK_col(NY,NX)
    ELSE
     TKSnow0_snvr(L,NY,NX)=TKSnow0_snvr(L-1,NY,NX)
    ENDIF

    if(TKSnow0_snvr(L,NY,NX)<200._r8)then
      write(iulog,*)'weird temp',I+J/24.,M,L,TKX,TKSnow0_snvr(L,NY,NX),TairK_col(NY,NX),ActiveSnow_test,VLSnowHeatCapM_snvr(M+1,L,NY,NX)
      write(iulog,*)L,VLDrySnoWE0_snvr(L,NY,NX),VLWatSnow0_snvr(L,NY,NX),VLIceSnow0_snvr(L,NY,NX)
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
    dENGY=CumHeat2SnowLM_snvr(L,NY,NX)+XPhaseChangeHeatLM_snvr(L,NY,NX)

  ENDDO D9780

  if(VLWatSnow0_snvr(1,NY,NX)<0._r8)then
    L=1
    write(*,*)'!make up the negative water from dry snow',VLWatSnow0_snvr(1,NY,NX),NY,NX
    call endrun(trim(mod_filename)//' at line',__LINE__)
    VLDrySnoWEtmp=VLDrySnoWE0_snvr(L,NY,NX)+VLWatSnow0_snvr(L,NY,NX)
    !this is equivalent to thaw, which releases latent heat
    if(VLDrySnoWEtmp>0._r8)then
      VLDrySnoWE0_snvr(L,NY,NX)        = VLDrySnoWEtmp
      ENGY0                            = VLSnowHeatCapM_snvr(M,L,NY,NX)*TKSnow0_snvr(L,NY,NX)
      VLSnowHeatCapM_snvr(M+1,L,NY,NX) = cps*VLDrySnoWE0_snvr(L,NY,NX)+cpw*VLWatSnow0_snvr(L,NY,NX)+cpi*VLIceSnow0_snvr(L,NY,NX)
      dHPhaseChange                    = -LtHeatIceMelt*VLWatSnow0_snvr(L,NY,NX)
      TKSnow0_snvr(L,NY,NX)            = (ENGY0+dHPhaseChange)/VLSnowHeatCapM_snvr(M+1,L,NY,NX)
      XPhaseChangeHeatLM_snvr(L,NY,NX) = XPhaseChangeHeatLM_snvr(L,NY,NX)+dHPhaseChange
      VLWatSnow0_snvr(L,NY,NX)         = 0._r8
    else
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif

  endif

  !update diagnostics
  DO L=1,JS
    CumSno2SnowL_snvr(L,NY,NX)      = CumSno2SnowL_snvr(L,NY,NX)+CumSno2SnowLM_snvr(L,NY,NX)
    CumWat2SnowL_snvr(L,NY,NX)      = CumWat2SnowL_snvr(L,NY,NX)+CumWat2SnowLM_snvr(L,NY,NX)
    CumIce2SnowL_snvr(L,NY,NX)      = CumIce2SnowL_snvr(L,NY,NX)+CumIce2SnowLM_snvr(L,NY,NX)
    CumHeat2SnowL_snvr(L,NY,NX)     = CumHeat2SnowL_snvr(L,NY,NX)+CumHeat2SnowLM_snvr(L,NY,NX)
    XSnowThawMassL_snvr(L,NY,NX)    = XSnowThawMassL_snvr(L,NY,NX)+XSnowThawMassLM_snvr(L,NY,NX)
    XIceThawMassL_snvr(L,NY,NX)     = XIceThawMassL_snvr(L,NY,NX)+XIceThawMassLM_snvr(L,NY,NX)
    XPhaseChangeHeatL_snvr(L,NY,NX) = XPhaseChangeHeatL_snvr(L,NY,NX)+XPhaseChangeHeatLM_snvr(L,NY,NX)
  ENDDO

  end subroutine UpdateSnowPack1M

!------------------------------------------------------------------------------------------
  subroutine AccumulateSnowRedisFluxM(I,J,M,NHW,NHE,NVN,NVS)
  ! Description:
  ! Accumulate accepted lateral snow constituents and sensible heat for native
  ! snow updates and the column heat balance, with positive heat into the cell.
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      CumSno2SnowL_snvr(1,NY,NX)  = CumSno2SnowL_snvr(1,NY,NX)+cumDrySnoByRedistM_col(NY,NX)
      CumWat2SnowL_snvr(1,NY,NX)  = CumWat2SnowL_snvr(1,NY,NX)+cumWatSnoByRedistM_col(NY,NX)
      CumIce2SnowL_snvr(1,NY,NX)  = CumIce2SnowL_snvr(1,NY,NX)+cumIceSnoByRedistM_col(NY,NX)
      CumHeat2SnowL_snvr(1,NY,NX) = CumHeat2SnowL_snvr(1,NY,NX)+cumHeatSnoByRedistM_col(NY,NX)
      ! This is internal transport, separate from atmospheric heat exchange.
      HeatLateral2Cell_col(NY,NX)=HeatLateral2Cell_col(NY,NX)+cumHeatSnoByRedistM_col(NY,NX)
    ENDDO
  ENDDO
  end subroutine AccumulateSnowRedisFluxM
!------------------------------------------------------------------------------------------
  subroutine UpdateSnowAtM(I,J,M,NY,NX)
  ! Description:
  ! Update snow after redistribution and transfer unresolved native films to
  ! the receiving surface with their phases and energy intact.
  implicit none
  integer, intent(in) :: M,NY,NX,I,J

  integer :: L
  real(r8) :: tk1pres
  real(r8) :: ENGY0,ENGY1
  real(r8) :: FLWI,FLWW,FLWS
  real(r8) :: HFLWS
  real(r8) :: dWatMac,dwatMic
!     SNOWPACK WATER, ICE, SNOW AND TEMPERATURE

  !      if(curday>=176)then
  !        write(*,*)'line',__LINE__,'tk1',TKSoil1_vr(8,ny,nx),TKSoil1_vr(9,ny,nx),M
  !      endif
  !
  !

  if (M.NE.NPH)then

    VLSnowHeatCapM_snvr(M+1,1,NY,NX) = cps*VLDrySnoWE0_snvr(1,NY,NX)+cpw*VLWatSnow0_snvr(1,NY,NX)+cpi*VLIceSnow0_snvr(1,NY,NX)

    if(any((/cumDrySnoByRedistM_col(NY,NX)/=0._r8,cumWatSnoByRedistM_col(NY,NX)/=0._r8, &
      cumIceSnoByRedistM_col(NY,NX)/=0._r8/)))then
      ENGY0   = VLSnowHeatCapM_snvr(M+1,1,NY,NX)*TKSnow0_snvr(1,NY,NX)
      VLDrySnoWE0_snvr(1,NY,NX) = VLDrySnoWE0_snvr(1,NY,NX)+cumDrySnoByRedistM_col(NY,NX)
      VLWatSnow0_snvr(1,NY,NX)  = VLWatSnow0_snvr(1,NY,NX)+cumWatSnoByRedistM_col(NY,NX)
      VLIceSnow0_snvr(1,NY,NX)  = VLIceSnow0_snvr(1,NY,NX)+cumIceSnoByRedistM_col(NY,NX)

      VLSnowHeatCapM_snvr(M+1,1,NY,NX) = cps*VLDrySnoWE0_snvr(1,NY,NX)+cpw*VLWatSnow0_snvr(1,NY,NX)+cpi*VLIceSnow0_snvr(1,NY,NX)

      IF(VLSnowHeatCapM_snvr(M+1,1,NY,NX).GT.0._r8)THEN
        TKSnow0_snvr(1,NY,NX)=(ENGY0+cumHeatSnoByRedistM_col(NY,NX))/VLSnowHeatCapM_snvr(M+1,1,NY,NX)
      ELSE
        TKSnow0_snvr(1,NY,NX)=TairK_col(NY,NX)
      ENDIF
    endif
  ELSE
    VLHeatCapSnow_snvr(1,NY,NX) = cps*VLDrySnoWE0_snvr(1,NY,NX)+cpw*VLWatSnow0_snvr(1,NY,NX)+cpi*VLIceSnow0_snvr(1,NY,NX)
    if(any((/cumDrySnoByRedistM_col(NY,NX)/=0._r8,cumWatSnoByRedistM_col(NY,NX)/=0._r8,&
      cumIceSnoByRedistM_col(NY,NX)/=0._r8/)))then
      ENGY0   = VLHeatCapSnow_snvr(1,NY,NX)*TKSnow0_snvr(1,NY,NX)
      VLDrySnoWE0_snvr(1,NY,NX) = VLDrySnoWE0_snvr(1,NY,NX)+cumDrySnoByRedistM_col(NY,NX)
      VLWatSnow0_snvr(1,NY,NX)  = VLWatSnow0_snvr(1,NY,NX)+cumWatSnoByRedistM_col(NY,NX)
      VLIceSnow0_snvr(1,NY,NX)  = VLIceSnow0_snvr(1,NY,NX)+cumIceSnoByRedistM_col(NY,NX)

      VLHeatCapSnow_snvr(1,NY,NX) = cps*VLDrySnoWE0_snvr(1,NY,NX)+cpw*VLWatSnow0_snvr(1,NY,NX)+cpi*VLIceSnow0_snvr(1,NY,NX)

      IF(VLHeatCapSnow_snvr(1,NY,NX).GT.0._r8)THEN
        TKSnow0_snvr(1,NY,NX)=(ENGY0+cumHeatSnoByRedistM_col(NY,NX))/VLHeatCapSnow_snvr(1,NY,NX)
      ELSE
        TKSnow0_snvr(1,NY,NX)=TairK_col(NY,NX)
      ENDIF
    endif
  ENDIF

  if(VLWatSnow0_snvr(1,NY,NX)<0._r8)then
    write(*,*)'UpdateSnowAtM NY,NX ',NY,NX,M,VLWatSnow0_snvr(1,NY,NX),cumWatSnoByRedistM_col(NY,NX)
    call endrun(trim(mod_filename)//' at line',__LINE__)
  endif
    !
    !     IF SNOWPACK DISAPPEARS ALL MATERIAL,HEAT TRANSFERRED TO SOIL SURFACE
    !
    !     VHeatCapacity1_vr,VLHeatCapacityA,VLHeatCapacityP=total soil,soil+micropore,macropore heat capacity
    !     TK1=soil surface temperature, why not to litter layer
    !
  ! Test the current constituents, including redistribution, on every substep.
  VLHeatCapSnow_snvr(1,NY,NX) = cps*VLDrySnoWE0_snvr(1,NY,NX) &
    +cpw*VLWatSnow0_snvr(1,NY,NX)+cpi*VLIceSnow0_snvr(1,NY,NX)
  IF(VLHeatCapSnow_snvr(1,NY,NX).LE.VLHeatCapSnowMin_col(NY,NX) .AND. &
    (TairK_col(NY,NX).GT.TFice.or.(.not.ats_cpl_mode.and..not.fixWaterLevel.and. &
      sum(VLDrySnoWE0_snvr(:,NY,NX))+sum(VLWatSnow0_snvr(:,NY,NX))+DENSICE*sum(VLIceSnow0_snvr(:,NY,NX)) &
      <=MinResolvedSnowSWE*AREA_3D(3,NUM_col(NY,NX),NY,NX))))THEN
    FLWS                      = VLDrySnoWE0_snvr(1,NY,NX)
    FLWW                      = VLWatSnow0_snvr(1,NY,NX)
    FLWI                      = VLIceSnow0_snvr(1,NY,NX)
    HFLWS                     = (cpw*FLWW+cps*FLWS+cpi*FLWI)*TKSnow0_snvr(1,NY,NX)
    VLDrySnoWE0_snvr(1,NY,NX) = 0._r8
    VLWatSnow0_snvr(1,NY,NX)  = 0._r8
    VLIceSnow0_snvr(1,NY,NX)  = 0._r8

    !add ice and water to litter layer
    tk1pres        = TKSoil1_vr(0,NY,NX)
    ENGY1          = TKSoil1_vr(0,NY,NX)*VHeatCapacity1_vr(0,NY,NX)

    !it is not a pond
    IF(SoilBulkDensity_vr(NUM_col(NY,NX),NY,NX).GT.ZERO .and. SoilOrgM_vr(ielmc,0,NY,NX)>1.e-2_r8)THEN
      VLWatMicP1_vr(0,NY,NX)     = VLWatMicP1_vr(0,NY,NX)+FLWW
      VLiceMicP1_vr(0,NY,NX)     = VLiceMicP1_vr(0,NY,NX)+FLWI+FLWS/DENSICE
      VHeatCapacity1_vr(0,NY,NX) = cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,0,NY,NX))+cpw*VLWatMicP1_vr(0,NY,NX)+cpi*VLiceMicP1_vr(0,NY,NX)  !update heat capacity
      TKSoil1_vr(0,NY,NX)        = (ENGY1+HFLWS)/VHeatCapacity1_vr(0,NY,NX)
      QSnowH2Oloss_col(NY,NX)    = QSnowH2Oloss_col(NY,NX)+FLWW+FLWI*DENSICE+FLWS
      QSnowHeatLoss_col(NY,NX)   = QSnowHeatLoss_col(NY,NX)+HFLWS

      WatFLo2LitR_col(NY,NX)     = WatFLo2LitR_col(NY,NX) + FLWW+FLWI*DENSICE+FLWS
    !pond
    else
      if(.not. fixWaterLevel)then
        !in ATS mode the soil state is owned by ATS, so only the fluxes below are kept
        if(.not. ats_cpl_mode)then
          VLWatMicP1_vr(NUM_col(NY,NX),NY,NX)      = VLWatMicP1_vr(NUM_col(NY,NX),NY,NX)+FLWW
          VLiceMicP1_vr(NUM_col(NY,NX),NY,NX)      = VLiceMicP1_vr(NUM_col(NY,NX),NY,NX)+FLWI+FLWS/DENSICE
          ENGY1                                = VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX)*TKSoil1_vr(NUM_col(NY,NX),NY,NX)
          VLHeatCapacityA_vr(NUM_col(NY,NX),NY,NX) = VHeatCapSolidSoil_vr(NUM_col(NY,NX),NY,NX) &
            +cpw*VLWatMicP1_vr(NUM_col(NY,NX),NY,NX)+cpi*VLiceMicP1_vr(NUM_col(NY,NX),NY,NX)
          VLHeatCapacityB_vr(NUM_col(NY,NX),NY,NX)=cpw*VLWatMacP1_vr(NUM_col(NY,NX),NY,NX) &
            +cpi*VLiceMacP1_vr(NUM_col(NY,NX),NY,NX)
          VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX) = VLHeatCapacityA_vr(NUM_col(NY,NX),NY,NX)+VLHeatCapacityB_vr(NUM_col(NY,NX),NY,NX)

          ! topsoil layer is there
          IF(VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
            tk1pres                      = TKSoil1_vr(NUM_col(NY,NX),NY,NX)
            TKSoil1_vr(NUM_col(NY,NX),NY,NX) = (ENGY1+HFLWS)/VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX)
          ELSE
            TKSoil1_vr(NUM_col(NY,NX),NY,NX) = TairK_col(NY,NX)
          ENDIF
        else
          !ATS receives the residual as liquid water through Qinflx2Soil_col and its
          !energy through HeatFlx2Grnd_col. The export subtracts cpw*TFice per unit of
          !Qinflx2Soil_col (ATS energy reference), so add that back here and send the
          !residual's energy relative to TFice, including the latent heat needed to
          !melt its ice and dry snow, which ATS cannot see.
          HeatFlx2Grnd_col(NY,NX) = HeatFlx2Grnd_col(NY,NX) &
            +cpw*(TKSnow0_snvr(1,NY,NX)-TFice)*FLWW &
            +cpi*(TKSnow0_snvr(1,NY,NX)-TFice)*FLWI &
            +cps*(TKSnow0_snvr(1,NY,NX)-TFice)*FLWS &
            -LtHeatIceMelt*(FLWI*DENSICE+FLWS) &
            +cpw*TFice*(FLWW+FLWI*DENSICE+FLWS)
        endif

        QSnoWatXfer2Soil_col(NY,NX) = QSnoWatXfer2Soil_col(NY,NX)+FLWW
        QSnoIceXfer2Soil_col(NY,NX) = QSnoIceXfer2Soil_col(NY,NX)+FLWI+FLWS/DENSICE
        QSnoHeatXfer2Soil_col(NY,NX)= QSnoHeatXfer2Soil_col(NY,NX)+HFLWS

        Qinflx2Soil_col(NY,NX)   = Qinflx2Soil_col(NY,NX)+FLWW+FLWI*DENSICE+FLWS
        QSnowH2Oloss_col(NY,NX)  = QSnowH2Oloss_col(NY,NX)+FLWW+FLWI*DENSICE+FLWS
        QSnowHeatLoss_col(NY,NX) = QSnowHeatLoss_col(NY,NX)+HFLWS
        Qinflx2SoilM_col(NY,NX)  = Qinflx2SoilM_col(NY,NX)+FLWW+FLWI*DENSICE+FLWS
      endif
    endif
  ENDIF

  ! Keep the next substep capacity consistent if the snow was cleared above.
  if(M.NE.NPH)then
    VLSnowHeatCapM_snvr(M+1,1,NY,NX) = cps*VLDrySnoWE0_snvr(1,NY,NX) &
      +cpw*VLWatSnow0_snvr(1,NY,NX)+cpi*VLIceSnow0_snvr(1,NY,NX)
  endif

  !update the snow state variables after iteration M
  D60: DO L=1,JS
    VLDrySnoWE_snvr(L,NY,NX)    = VLDrySnoWE0_snvr(L,NY,NX)
    VLIceSnow_snvr(L,NY,NX)     = VLIceSnow0_snvr(L,NY,NX)
    VLWatSnow_snvr(L,NY,NX)     = VLWatSnow0_snvr(L,NY,NX)
    VLSnoDWIprev_snvr(L,NY,NX)  = VLDrySnoWE_snvr(L,NY,NX)/SnoDens_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)
    VLHeatCapSnow_snvr(L,NY,NX) = cps*VLDrySnoWE_snvr(L,NY,NX)+cpw*VLWatSnow_snvr(L,NY,NX)+cpi*VLIceSnow_snvr(L,NY,NX)
    SnowThickL_snvr(L,NY,NX)    = VLSnoDWIprev_snvr(L,NY,NX)/AREA_3D(3,NUM_col(NY,NX),NY,NX)
    TKSnow_snvr(L,NY,NX)        = TKSnow0_snvr(L,NY,NX)
    TCSnow_snvr(L,NY,NX)        = units%Kelvin2Celcius(TKSnow_snvr(L,NY,NX))
  ENDDO D60

  end subroutine UpdateSnowAtM

!------------------------------------------------------------------------------------------
  subroutine SnowRedistributionM(M,NY,NX,NHE,NHW,NVS,NVN)
  ! Route top-layer constituents together, sharing one donor budget across
  ! all neighbors. Incoming snow is available only on the next substep.
  implicit none
  integer, intent(in) :: M,NY,NX,NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='SnowRedistributionM'
  integer :: N,NN,N4,N5,FaceSide,Step
  real(r8) :: Fraction(2,2),SnowDepthGradH,QSX,DonorScale,FaceDistance

  call PrintInfo('beg '//subname)
  Fraction=0._r8
  DO N=1,2
    DO NN=1,2
      Step=1
      IF(NN.EQ.iBehind)Step=-1
      N4=NX; N5=NY
      IF(N.EQ.iWestEastDirection)THEN
        N4=NX+Step
      ELSE
        N5=NY+Step
      ENDIF
      IF(N4.LT.NHW .OR. N4.GT.NHE .OR. N5.LT.NVN .OR. N5.GT.NVS)CYCLE
      ! Horizontal center distance is independent of soil-layer numbering.
      IF(N.EQ.iWestEastDirection)THEN
        FaceDistance=0.5_r8*(DH_col(NY,NX)+DH_col(N5,N4))
      ELSE
        FaceDistance=0.5_r8*(DV_col(NY,NX)+DV_col(N5,N4))
      ENDIF
      ! Retain the snow-depth gradient formulation, using the actual neighbor
      ! in both directions. Donor states have completed local physics.
      SnowDepthGradH=(Altitude_col(NY,NX)+SnowDepth_col(NY,NX) &
        -Altitude_col(N5,N4)-SnowDepth_col(N5,N4))/FaceDistance
      IF(SnowDepthGradH.GT.0._r8 .AND. SnowDepth_col(NY,NX).GT.MinSnowDepth)THEN
        Fraction(N,NN)=MIN(1._r8,SnowDepthGradH/MAX(1._r8,FaceDistance)*dts_HeatWatTP)
      ENDIF
    ENDDO
  ENDDO

  DonorScale=1._r8/MAX(1._r8,SUM(Fraction))
  DO N=1,2
    DO NN=1,2
      Step=1
      FaceSide=iBehind
      IF(NN.EQ.iBehind)THEN
        Step=-1
        FaceSide=iFront
      ENDIF
      N4=NX; N5=NY
      IF(N.EQ.iWestEastDirection)THEN
        N4=NX+Step
      ELSE
        N5=NY+Step
      ENDIF
      IF(N4.LT.NHW .OR. N4.GT.NHE .OR. N5.LT.NVN .OR. N5.GT.NVS)CYCLE
      QSX=Fraction(N,NN)*DonorScale
      DrySnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4)=QSX*AZMAX1(VLDrySnoWE0_snvr(1,NY,NX))
      WatSnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4)=QSX*AZMAX1(VLWatSnow0_snvr(1,NY,NX))
      IceSnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4)=QSX*AZMAX1(VLIceSnow0_snvr(1,NY,NX))
      HeatSnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4)=TKSnow0_snvr(1,NY,NX) &
        *(cps*DrySnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4)+cpw*WatSnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4) &
        +cpi*IceSnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4))
      IFLBSM_2DH(M,N,FaceSide,N5,N4)=-1
      IF(QSX.GT.0._r8)THEN
        IFLBSM_2DH(M,N,FaceSide,N5,N4)=0
        IF(FaceSide.EQ.iFront)IFLBSM_2DH(M,N,FaceSide,N5,N4)=1
      ENDIF
      ! Accumulate only accepted fluxes; the state and tracer consumers use
      ! these same substep faces, including the common donor reduction.
      DrySnoByRedistrib_2DH(N,FaceSide,N5,N4)=DrySnoByRedistrib_2DH(N,FaceSide,N5,N4)+DrySnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4)
      WatSnoByRedist_2DH(N,FaceSide,N5,N4)=WatSnoByRedist_2DH(N,FaceSide,N5,N4)+WatSnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4)
      IceSnoBySnowRedist_2DH(N,FaceSide,N5,N4)=IceSnoBySnowRedist_2DH(N,FaceSide,N5,N4)+IceSnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4)
      HeatSnoByRedist_2DH(N,FaceSide,N5,N4)=HeatSnoByRedist_2DH(N,FaceSide,N5,N4)+HeatSnoFlxByRedistM_2DH(M,N,FaceSide,N5,N4)
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SnowRedistributionM

!------------------------------------------------------------------------------------------
  subroutine SumSnowDriftByRunoffM(M,N,N1,N2,N4,N5,N4B,N5B)
  use SoilWaterDataType, only : IFLBM_2DH
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2
  integer, intent(in) :: N4,N5    !forward dest grid
  integer, intent(in) :: N4B,N5B  !backward dest grid
  integer :: NN

  DO NN=1,2
    cumDrySnoByRedistM_col(N2,N1)  = cumDrySnoByRedistM_col(N2,N1)+DrySnoFlxByRedistM_2DH(M,N,NN,N2,N1)
    cumWatSnoByRedistM_col(N2,N1)  = cumWatSnoByRedistM_col(N2,N1)+WatSnoFlxByRedistM_2DH(M,N,NN,N2,N1)
    cumIceSnoByRedistM_col(N2,N1)  = cumIceSnoByRedistM_col(N2,N1)+IceSnoFlxByRedistM_2DH(M,N,NN,N2,N1)
    cumHeatSnoByRedistM_col(N2,N1) = cumHeatSnoByRedistM_col(N2,N1)+HeatSnoFlxByRedistM_2DH(M,N,NN,N2,N1)

    TDrysnoByRedist_col(N2,N1)  = TDrysnoByRedist_col(N2,N1)+DrySnoFlxByRedistM_2DH(M,N,NN,N2,N1)
    TWatSnoByRedist_col(N2,N1)  = TWatSnoByRedist_col(N2,N1)+WatSnoFlxByRedistM_2DH(M,N,NN,N2,N1)
    TIceSnoByRedist_col(N2,N1)  = TIceSnoByRedist_col(N2,N1)+IceSnoFlxByRedistM_2DH(M,N,NN,N2,N1)
    THeatSnoByRedist_col(N2,N1) = THeatSnoByRedist_col(N2,N1)+HeatSnoFlxByRedistM_2DH(M,N,NN,N2,N1)

  ENDDO

  if(IFLBSM_2DH(M,N,iBehind,N5,N4).EQ.0)then
    NN                             = iBehind
    cumDrySnoByRedistM_col(N2,N1)  = cumDrySnoByRedistM_col(N2,N1)-DrySnoFlxByRedistM_2DH(M,N,NN,N5,N4)
    cumWatSnoByRedistM_col(N2,N1)  = cumWatSnoByRedistM_col(N2,N1)-WatSnoFlxByRedistM_2DH(M,N,NN,N5,N4)
    cumIceSnoByRedistM_col(N2,N1)  = cumIceSnoByRedistM_col(N2,N1)-IceSnoFlxByRedistM_2DH(M,N,NN,N5,N4)
    cumHeatSnoByRedistM_col(N2,N1) = cumHeatSnoByRedistM_col(N2,N1)-HeatSnoFlxByRedistM_2DH(M,N,NN,N5,N4)

    TDrysnoByRedist_col(N2,N1)  = TDrysnoByRedist_col(N2,N1)-DrySnoFlxByRedistM_2DH(M,N,NN,N5,N4)
    TWatSnoByRedist_col(N2,N1)  = TWatSnoByRedist_col(N2,N1)-WatSnoFlxByRedistM_2DH(M,N,NN,N5,N4)
    TIceSnoByRedist_col(N2,N1)  = TIceSnoByRedist_col(N2,N1)-IceSnoFlxByRedistM_2DH(M,N,NN,N5,N4)
    THeatSnoByRedist_col(N2,N1) = THeatSnoByRedist_col(N2,N1)-HeatSnoFlxByRedistM_2DH(M,N,NN,N5,N4)
  ENDIF

  if(N4B.GT.0 .and. N5B.GT.0)then
    NN                             = iFront
    cumDrySnoByRedistM_col(N2,N1)  = cumDrySnoByRedistM_col(N2,N1)-DrySnoFlxByRedistM_2DH(M,N,NN,N5B,N4B)
    cumWatSnoByRedistM_col(N2,N1)  = cumWatSnoByRedistM_col(N2,N1)-WatSnoFlxByRedistM_2DH(M,N,NN,N5B,N4B)
    cumIceSnoByRedistM_col(N2,N1)  = cumIceSnoByRedistM_col(N2,N1)-IceSnoFlxByRedistM_2DH(M,N,NN,N5B,N4B)
    cumHeatSnoByRedistM_col(N2,N1) = cumHeatSnoByRedistM_col(N2,N1)-HeatSnoFlxByRedistM_2DH(M,N,NN,N5B,N4B)

    TDrysnoByRedist_col(N2,N1)  = TDrysnoByRedist_col(N2,N1)-DrySnoFlxByRedistM_2DH(M,N,NN,N5B,N4B)
    TWatSnoByRedist_col(N2,N1)  = TWatSnoByRedist_col(N2,N1)-WatSnoFlxByRedistM_2DH(M,N,NN,N5B,N4B)
    TIceSnoByRedist_col(N2,N1)  = TIceSnoByRedist_col(N2,N1)-IceSnoFlxByRedistM_2DH(M,N,NN,N5B,N4B)
    THeatSnoByRedist_col(N2,N1) = THeatSnoByRedist_col(N2,N1)-HeatSnoFlxByRedistM_2DH(M,N,NN,N5B,N4B)
  ENDIF

  END subroutine SumSnowDriftByRunoffM
!------------------------------------------------------------------------------------------
  subroutine ZeroSnowFluxM(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX

  cumDrySnoByRedistM_col(NY,NX)  = 0.0_r8
  cumWatSnoByRedistM_col(NY,NX)  = 0.0_r8
  cumIceSnoByRedistM_col(NY,NX)  = 0.0_r8
  cumHeatSnoByRedistM_col(NY,NX) = 0.0_r8

  end subroutine ZeroSnowFluxM
!------------------------------------------------------------------------------------------
  FUNCTION CalcSnowBNDResistance(NY,NX)result(RAS)
  ! Snow-pore gas resistance using column-mean depth; the gas-path assembler
  ! converts it to snow-covered-patch depth before mixing surface pathways.
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: RAS
  real(r8) :: RASL,RASX
  real(r8) :: FracAsAirSno
  real(r8) :: TFACW
  integer :: L

!     RAS,RASL=blrs of snowpack,snowpack layer
!     VOLS,VLSnoDWI1=volume of snowpack,snowpack layer
!     DLYRS=snowpack later depth
!     WGSGW=vapor diffusivity in snowpack
!     FracAsAirSno=snowpack air-filled porosity
!     THETPI=air content of ice
!     VOLS0,VOLI0,VOLW0,VLSnoDWI1=snow,ice,water,total volumes of snowpack

  RAS=0.0_r8
  IF(VcumSnoDWI_col(NY,NX).GT.ZEROS2(NY,NX))THEN
    D9775: DO L=1,JS
      IF(VLSnoDWI1_snvr(L,NY,NX).GT.0._r8)THEN
        if(ats_cpl_mode)then
          !In the coupler I need these defined which are typically done in hour1
          TFACW=TEFGASDIF(TKSnow_snvr(L,NY,NX))
          H2OVapDifsc_snvr(L,NY,NX)=WGSG*TFACW
        endif
        RASX         = SnowThickL_snvr(L,NY,NX)/H2OVapDifsc_snvr(L,NY,NX)
        FracAsAirSno = AMAX1(THETPI,1.0_r8-(VLDrySnoWE0_snvr(L,NY,NX)+VLIceSnow0_snvr(L,NY,NX) &
            +VLWatSnow0_snvr(L,NY,NX))/VLSnoDWI1_snvr(L,NY,NX))
        RASL = RASX/AMAX1(ZERO,FracAsAirSno)**2
        RAS  = RAS+RASL
      ENDIF
    ENDDO D9775
  ENDIF
  END FUNCTION CalcSnowBNDResistance

!------------------------------------------------------------------------------------------

  subroutine SnowSurfLitRIteration(I,J,dt_SnoHeat,L,M,NY,NX,AvgThermCondctSnoLitR,AvgThermCondctSoilLitR,&
    AvgVaporCondctSoilLitR,AvgVaporCondctSnowLitR,PSISV1,VLairSno1, &
    TotWatXFlx2SoiMicP,TotHeatFlow2Soi,TotSnoWatFlow2Litr,TotSnoHeatFlow2Litr,SnowLiquidBudget)
  implicit none
  integer , intent(in) :: I,J
  real(r8), intent(in) :: dt_SnoHeat            !time step size for snow iteration
  integer, intent(in) :: M   !iteration id
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: AvgThermCondctSnoLitR,AvgThermCondctSoilLitR,AvgVaporCondctSnowLitR
  real(r8), intent(in) :: AvgVaporCondctSoilLitR,PSISV1,VLairSno1
  real(r8), intent(inout) :: TotWatXFlx2SoiMicP
  real(r8), intent(inout) :: TotHeatFlow2Soi
  real(r8), intent(inout) :: TotSnoWatFlow2Litr
  real(r8), intent(inout) :: TotSnoHeatFlow2Litr
  real(r8), intent(inout) :: SnowLiquidBudget

  character(len=*), parameter :: subname='SnowSurfLitRIteration'
  real(r8)  :: CumVapFlxSno2Litr        !cumulative water flux from snow to litter
  real(r8)  :: CumVapFlxLitr2Soi        !cumulative water flux from litter to soil
  real(r8)  :: cumHeatConvFlxLitr2Soi1  !convective heat flux (associated with water flux) from litter to soil
  real(r8)  :: CumHeatConvFlxSno2Litr   !convective heat flux (associated with water flux) from snow to litter
  real(r8)  :: CumHeatCndFlxSno2Litr    !conductive heat flux from snow to litter
  real(r8)  :: cumHeatCndFlxLitr2Soi    !conductive heat flux from litter to soil
  integer :: NN
  real(r8)  :: TK0X         !snowpack temperature, [K]
  real(r8)  :: TKXR         !litter temperature,   [K]
  real(r8)  :: TK1X         !surface soil temperature, [K]
  real(r8) :: VP1,VPR,VapSnow0,VPY
  real(r8) :: H2OVapFlx,VapConvFlxLitr2Soi
  real(r8) :: HeatConvFlxLitr2Soi
  real(r8) :: VapFlxSno2Litr,HeatConvFlxSno2Litr,HeatCndFlxSno2Litr
  real(r8) :: H2OVapFlxMax,HeatCnduct,HeatCndFlxLitr2Soi
  real(r8) :: HeatCnductMax,TKY
  real(r8) :: VHeatCapacityLitR
  real(r8) :: VHeatCapacitySoil
  real(r8) :: VLWatLitR
  real(r8) :: VLWatSoil
  real(r8) :: ENGYS,ENGYR
  real(r8) :: dts_litrvapht          !time step for litter vapor flux calculation embedded in each snow iteration
  real(r8) :: dLWRaddTLitR
  real(r8) :: dLWRaddTSoil

! begin_execution
  call PrintInfo('beg '//subname)

  TK0X = TKSnow1_snvr(L,NY,NX)
  TKXR = TKSoil1_vr(0,NY,NX)
  TK1X = TKSoil1_vr(NUM_col(NY,NX),NY,NX)

  VHeatCapacityLitR = VHeatCapacity1_vr(0,NY,NX)
  VLWatLitR         = VLWatMicP1_vr(0,NY,NX)
  VHeatCapacitySoil = VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX)
  VLWatSoil         = VLWatMicP1_vr(NUM_col(NY,NX),NY,NX)
  dts_litrvapht     = dt_SnoHeat/(real(NPR,kind=r8))

  CumVapFlxSno2Litr       = 0.0_r8
  CumVapFlxLitr2Soi       = 0.0_r8
  cumHeatConvFlxLitr2Soi1 = 0.0_r8    
  CumHeatConvFlxSno2Litr  = 0.0_r8
  CumHeatCndFlxSno2Litr   = 0.0_r8
  cumHeatCndFlxLitr2Soi   = 0.0_r8

  D4000: DO NN = 1, NPR
    dLWRaddTLitR = -4._r8*LWEmscefLitR_col(NY,NX)*TKXR**3/(real(NPR,kind=r8))
    dLWRaddTSoil = -4._r8*LWEmscefSoil_col(NY,NX)*TK1X**3/(real(NPR,kind=r8))
    !
    MMit=MMit+1
    ! VAPOR FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
    !
    ! VapSnow0,VPR,VPY=snowpack,litter, equilibrium vapor concentration, [ton H2O]
    ! TK0X,TKXR=snowpack,litter temperature
    ! PSISM1=litter matric water potential
    ! H2OVapFlx,H2OVapFlxMax=vapor-unconstrained,vapor-constrained vapor flux
    ! AREA=area of grid cell
    ! FSNW,CVRD=snow,litter cover fraction
    ! dts_litrvapht=time step for flux calculation
    ! FLVRSX=snow-litter vapor flux
    ! HFLVRSX=convective heat flux from snow-litter vapor flux
    !
    !residue vapor pressure
    VPR=vapsat(TKXR)*EXP(18.0_r8*PSISM1_vr(0,NY,NX)/(RGASC*TKXR))
    
    if(abs(VPR)>1.e20_r8)then
      write(*,*)'TKXR=',TKXR,TKSoil1_vr(0,NY,NX),TKSoil1_vr(NUM_col(NY,NX),NY,NX),NN
      write(*,*)'PSISM1_vr(0,NY,NX)=',PSISM1_vr(0,NY,NX)
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif

    !AvgVaporCondctSnowLitR=snowpack-litter conductance
    !there is litter layer with sufficient moisture
    IF(VLairSno1.GT.ZEROS2(NY,NX) .AND. FracAirFilledSoilPoreM_vr(M,0,NY,NX).GT.AirFillPore_Min)THEN
      VapSnow0=vapsat(TK0X)  !*dssign(VLWatSnow0M_snvr(L,NY,NX)) !snow vapor pressure, saturated
      !snow -> residue vapor flux (>0)
      H2OVapFlx=AvgVaporCondctSnowLitR*(VapSnow0-VPR)*AREA_3D(3,NUM_col(NY,NX),NY,NX)*FracSurfSnowCovered_col(NY,NX) &
        *FracSurfByLitR_col(NY,NX)*dts_litrvapht
      !equilibrium vapor pressure
      VPY          = (VapSnow0*VLairSno1+VPR*VLsoiAirPM_vr(M,0,NY,NX))/(VLairSno1+VLsoiAirPM_vr(M,0,NY,NX))
      H2OVapFlxMax = (VapSnow0-VPY)*VLairSno1*dt_watvap  !from snow to litter

      IF(H2OVapFlx.GE.0.0_r8)THEN
        !vapor flux from snow to litter
        VapFlxSno2Litr      = AZMAX1(AMIN1(H2OVapFlx,H2OVapFlxMax,VLWatSnow0M_snvr(L,NY,NX)*dts_litrvapht,SnowLiquidBudget))      !water flux
        HeatConvFlxSno2Litr = (cpw*TK0X+EvapLHTC)*VapFlxSno2Litr         !enthalpy flux associated with water flux
      ELSE
        !vapor flux from litter to snow
        VapFlxSno2Litr      = AZMIN1(AMAX1(H2OVapFlx,H2OVapFlxMax,-AZMAX1(VLWatLitR)*dts_litrvapht))
        HeatConvFlxSno2Litr = (cpw*TKXR+EvapLHTC)*VapFlxSno2Litr
      ENDIF
    ELSE
      VapFlxSno2Litr      = 0.0_r8
      HeatConvFlxSno2Litr = 0.0_r8
    ENDIF
    SnowLiquidBudget=MAX(0._r8,SnowLiquidBudget-MAX(0._r8,VapFlxSno2Litr))
    !
    ! HEAT FLUX BETWEEN SNOWPACK AND SURFACE RESIDUE
    !
    ! TKY=snow-litter equilibrium temperature, assuming simple mass summation
    ! HeatCnduct,HeatCnductMax=snow-litter heat flux unltd,ltd by heat
    ! HeatCndFlxSno2Litr=snow-litter heat flux
    ! VLSnowHeatCapMM= volumetric heat capacity in snow layer
    TKY           = (TK0X*VLHeatCapSnowM1_snvr(L,NY,NX)+TKXR*VHeatCapacityLitR)/(VLHeatCapSnowM1_snvr(L,NY,NX)+VHeatCapacityLitR)
    HeatCnductMax = (TK0X-TKY)*VLHeatCapSnowM1_snvr(L,NY,NX)*dt_watvap
    HeatCnduct    = AvgThermCondctSnoLitR*(TK0X-TKXR)*AREA_3D(3,NUM_col(NY,NX),NY,NX)*FracSurfSnowCovered_col(NY,NX)*FracSurfByLitR_col(NY,NX)*dts_litrvapht
    IF(HeatCnduct.GE.0.0_r8)THEN
      !sensible heat flux from snow to litter
      HeatCndFlxSno2Litr=AZMAX1(AMIN1(HeatCnductMax,HeatCnduct))
    ELSE
      !sensible heat flux from litter to snow
      HeatCndFlxSno2Litr=AZMIN1(AMAX1(HeatCnductMax,HeatCnduct))
    ENDIF
    !
    !     VAPOR FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
    !
    !     THETPM,VsoiPM=air-filled porosity,volume
    !     VP1,VPY=soil,litter-soil equilibrium vapor concentration
    !     TK1X=soil temperature
    !     PSISV1=soil matric+osmotic water potentials
    !     H2OVapFlx,H2OVapFlxMax=vapor-unconstrained,vapor-constrained vapor flux
    !     VapConvFlxLitr2Soi=litter-soil vapor flux
    !     HeatConvFlxLitr2Soi=convective heat of litter-soil vapor flux
    !     TKXR,TK1X=interim calculation of litter,soil temperatures
    !
    !both litter layer and topsoil are not-saturated
    IF(VLsoiAirPM_vr(M,0,NY,NX).GT.ZEROS(NY,NX) .AND. VLsoiAirPM_vr(M,NUM_col(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
      VP1       = vapsat(TK1X)*EXP(18.0_r8*PSISV1/(RGASC*TK1X))
      H2OVapFlx = AvgVaporCondctSoilLitR*(VPR-VP1)*AREA_3D(3,NUM_col(NY,NX),NY,NX)*FracSurfSnowCovered_col(NY,NX)*FracSurfByLitR_col(NY,NX)*dts_litrvapht

      VPY=(VPR*VLsoiAirPM_vr(M,0,NY,NX)+VP1*VLsoiAirPM_vr(M,NUM_col(NY,NX),NY,NX)) &
        /(VLsoiAirPM_vr(M,0,NY,NX)+VLsoiAirPM_vr(M,NUM_col(NY,NX),NY,NX))
      H2OVapFlxMax=(VPR-VPY)*VLsoiAirPM_vr(M,0,NY,NX)*dt_watvap

      IF(H2OVapFlx.GE.0.0_r8)THEN
        !vapor flux from litter to soil
        VapConvFlxLitr2Soi  = AZMAX1(AMIN1(H2OVapFlx,H2OVapFlxMax,VLWatLitR*dts_litrvapht, &
            MAX(0._r8,VLWatLitR+MIN(0._r8,VapFlxSno2Litr))))
        HeatConvFlxLitr2Soi = (cpw*TKXR*HeatAdv_scal+EvapLHTC)*VapConvFlxLitr2Soi  !enthalpy flux
      ELSE
        !vapor flux from soil to litter
        VapConvFlxLitr2Soi  = AZMIN1(AMAX1(H2OVapFlx,H2OVapFlxMax,-AZMAX1(VLWatSoil)*dts_litrvapht))
        HeatConvFlxLitr2Soi = (cpw*TK1X*HeatAdv_scal+EvapLHTC)*VapConvFlxLitr2Soi
      ENDIF
    ELSE
      VapConvFlxLitr2Soi  = 0.0_r8
      HeatConvFlxLitr2Soi = 0.0_r8
    ENDIF
    !update litter layer temperature
    ENGYR             = TKXR*VHeatCapacityLitR
    VLWatLitR         = VLWatLitR+VapFlxSno2Litr-VapConvFlxLitr2Soi
    VHeatCapacityLitR = VHeatCapacityLitR+cpw*(VapFlxSno2Litr-VapConvFlxLitr2Soi)
    TKXR              = (ENGYR+HeatConvFlxSno2Litr-HeatConvFlxLitr2Soi-dLWRaddTLitR*TKXR)/(VHeatCapacityLitR-dLWRaddTLitR)

    !update top soil layer temperature
    ENGYS             = TK1X*VHeatCapacitySoil
    VHeatCapacitySoil = VHeatCapacitySoil+cpw*VapConvFlxLitr2Soi
    VLWatSoil         = VLWatSoil + VapConvFlxLitr2Soi
    TK1X              = (ENGYS+HeatConvFlxLitr2Soi-dLWRaddTSoil*TK1X)/(VHeatCapacitySoil-dLWRaddTSoil)
    !
    !     sensible HEAT FLUX BETWEEN SURFACE RESIDUE AND SOIL SURFACE
    !
    !     TKY=litter-soil equilibrium temperature
    !     HeatCnduct,HeatCnductMax=litter-soil heat flux unltd,ltg by heat
    !     HeatCndFlxLitr2Soi=litter-soil heat flux
    !
    TKY           = (TKXR*VHeatCapacityLitR+TK1X*VHeatCapacitySoil)/(VHeatCapacityLitR+VHeatCapacitySoil)
    HeatCnductMax = (TKXR-TKY)*VHeatCapacityLitR*dt_watvap
    HeatCnduct    = AvgThermCondctSoilLitR*(TKXR-TK1X)*AREA_3D(3,NUM_col(NY,NX),NY,NX)*FracSurfSnowCovered_col(NY,NX)*FracSurfByLitR_col(NY,NX)*dts_litrvapht

    IF(HeatCnduct.GE.0.0_r8)THEN
      !sensible heat flux from litter to soil
      HeatCndFlxLitr2Soi=AZMAX1(AMIN1(HeatCnductMax,HeatCnduct))
    ELSE
      !sensible heat flux from soil to litter
      HeatCndFlxLitr2Soi=AZMIN1(AMAX1(HeatCnductMax,HeatCnduct))
    ENDIF
!
!     ACCUMULATE SNOW-LITTER, LITTER-SOIL HEAT FLUXES
!     WITHIN LONGER TIME STEP FOR SNOWPACK FLUX CALCULATIONS
!
    CumVapFlxSno2Litr      = CumVapFlxSno2Litr+VapFlxSno2Litr
    CumHeatConvFlxSno2Litr = CumHeatConvFlxSno2Litr+HeatConvFlxSno2Litr
    CumHeatCndFlxSno2Litr  = CumHeatCndFlxSno2Litr+HeatCndFlxSno2Litr

    CumVapFlxLitr2Soi       = CumVapFlxLitr2Soi+VapConvFlxLitr2Soi
    cumHeatConvFlxLitr2Soi1 = cumHeatConvFlxLitr2Soi1+HeatConvFlxLitr2Soi
    cumHeatCndFlxLitr2Soi   = cumHeatCndFlxLitr2Soi+HeatCndFlxLitr2Soi

    TK0X = TK0X-HeatConvFlxSno2Litr/VLHeatCapSnowM1_snvr(L,NY,NX)
    TKXR = (TKXR*VHeatCapacityLitR-dLWRaddTLitR*TKXR-HeatCndFlxLitr2Soi)/(VHeatCapacityLitR-dLWRaddTLitR)
    TK1X = (TK1X*VHeatCapacitySoil-dLWRaddTSoil*TK1X+HeatCndFlxLitr2Soi)/(VHeatCapacitySoil-dLWRaddTSoil)

    if(TKXR<0._r8)then
      write(*,*)'unphysical litter temeprature'
      write(*,*)'TKXR',TKXR,HeatConvFlxSno2Litr,HeatCndFlxLitr2Soi,VHeatCapacityLitR
      call endrun(trim(mod_filename)//'at line',__LINE__)
    endif

  ENDDO D4000

  TotWatXFlx2SoiMicP = TotWatXFlx2SoiMicP+CumVapFlxLitr2Soi
  TotHeatFlow2Soi     = TotHeatFlow2Soi+cumHeatConvFlxLitr2Soi1+cumHeatCndFlxLitr2Soi
  TotSnoWatFlow2Litr  = TotSnoWatFlow2Litr+CumVapFlxSno2Litr-CumVapFlxLitr2Soi       
  TotSnoHeatFlow2Litr = TotSnoHeatFlow2Litr+CumHeatConvFlxSno2Litr+CumHeatCndFlxSno2Litr &
    -cumHeatConvFlxLitr2Soi1-cumHeatCndFlxLitr2Soi

  call PrintInfo('end '//subname)
  end subroutine SnowSurfLitRIteration

!------------------------------------------------------------------------------------------
  subroutine StageSnowModel(I,J,NY,NX)
  implicit none
  integer, intent(in) :: NY,NX,I,J
  character(len=*), parameter :: subname='StageSnowModel'

  integer :: L


  ! VOLS0,VOLSSL=snowpack snow content (water equivalent)
  ! VOLI0,VOLISSL=snowpack ice content
  ! VOLW0,VOLWSL=snowpack water content
  ! VLSnoDWI1,VLSnoDWIprev_snvr=snowpack volume
  ! SnowThickL0_snvr,DLYRS=snowpack depth
  ! VLSnowHeatCapM,VHCPW=snowpack heat capacity
  ! TK0,TKW=snowpack temperature
  call PrintInfo('beg '//subname)
  TDrysnoByRedist_col(NY,NX)  = 0._r8
  TWatSnoByRedist_col(NY,NX)  = 0._r8
  TIceSnoByRedist_col(NY,NX)  = 0._r8
  THeatSnoByRedist_col(NY,NX) = 0._r8

  D60: DO L=1,JS
    VLDrySnoWE0_snvr(L,NY,NX)      = VLDrySnoWE_snvr(L,NY,NX)
    VLIceSnow0_snvr(L,NY,NX)       = VLIceSnow_snvr(L,NY,NX)
    VLWatSnow0_snvr(L,NY,NX)       = VLWatSnow_snvr(L,NY,NX)
    VLSnoDWI1_snvr(L,NY,NX)        = VLSnoDWIprev_snvr(L,NY,NX)
    SnowThickL0_snvr(L,NY,NX)      = SnowThickL_snvr(L,NY,NX)
    VLSnowHeatCapM_snvr(1,L,NY,NX) = VLHeatCapSnow_snvr(L,NY,NX)
    TKSnow0_snvr(L,NY,NX)          = TKSnow_snvr(L,NY,NX)

    CumSno2SnowL_snvr(L,NY,NX)   = 0.0_r8
    CumWat2SnowL_snvr(L,NY,NX)   = 0.0_r8
    CumIce2SnowL_snvr(L,NY,NX)   = 0.0_r8
    CumHeat2SnowL_snvr(L,NY,NX)  = 0.0_r8
    XSnowThawMassL_snvr(L,NY,NX) = 0._r8
    XIceThawMassL_snvr(L,NY,NX)  = 0._r8
  ENDDO D60
  call PrintInfo('end '//subname)
  end subroutine StageSnowModel
!------------------------------------------------------------------------------------------

  subroutine UpdateSoilWaterPotential(NY,NX)
  ! Description:
  ! Refresh the working surface pressure for snow exchange. Ordinary soil shares
  ! the common retention evaluator; the existing pond/ice branches remain here.
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: THETW1,VPY,FCX,WPX
  real(r8) :: LOGFCX,LOGWPX,PSDX,FCDX
  character(len=*), parameter :: subname='UpdateSoilWaterPotential'

  call PrintInfo('beg '//subname)
  THETW1=AMAX1(SoilWatAirDry_vr(NUM_col(NY,NX),NY,NX),AMIN1(POROS_vr(NUM_col(NY,NX),NY,NX) &
    ,safe_adb(VLWatMicP1_vr(NUM_col(NY,NX),NY,NX),VLSoilMicP_vr(NUM_col(NY,NX),NY,NX))))
  IF(VLSoilMicPMass_vr(NUM_col(NY,NX),NY,NX).GT.ZEROS(NY,NX))THEN
    call ComputePsiMCM(NUM_col(NY,NX),NY,NX,THETW1,PSISM1_vr(NUM_col(NY,NX),NY,NX))
  ELSEIF(VLSoilPoreMicP_vr(NUM_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    FCX    = FCI*FracSoiPAsIce_vr(NUM_col(NY,NX),NY,NX)
    WPX    = WPI*FracSoiPAsIce_vr(NUM_col(NY,NX),NY,NX)
    LOGFCX = LOG(FCX)
    LOGWPX = LOG(WPX)
    PSDX   = LOGPOROS_vr(NUM_col(NY,NX),NY,NX)-LOGFCX
    FCDX   = LOGFCX-LOGWPX
    IF(FracSoiPAsWat_vr(NUM_col(NY,NX),NY,NX).LT.FCX)THEN
      PSISM1_vr(NUM_col(NY,NX),NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD_col(NY,NX) &
        +((LOGFCX-LOG(FracSoiPAsWat_vr(NUM_col(NY,NX),NY,NX)))/FCDX*LOGPSIMND_col(NY,NX))))
    ELSEIF(FracSoiPAsWat_vr(NUM_col(NY,NX),NY,NX).LT.POROS_vr(NUM_col(NY,NX),NY,NX)-DTHETW)THEN
      PSISM1_vr(NUM_col(NY,NX),NY,NX)=-EXP(LOGPSIAtSat(NY,NX) &
        +(((LOGPOROS_vr(NUM_col(NY,NX),NY,NX)-LOG(FracSoiPAsWat_vr(NUM_col(NY,NX),NY,NX)))/PSDX)*LOGPSIMXD_col(NY,NX)))
    ELSE
      THETW1                      = POROS_vr(NUM_col(NY,NX),NY,NX)
      PSISM1_vr(NUM_col(NY,NX),NY,NX) = PSI_at_Sautration_vr(NUM_col(NY,NX),NY,NX)
    ENDIF
  ELSE
    THETW1                      = POROS_vr(NUM_col(NY,NX),NY,NX)
    PSISM1_vr(NUM_col(NY,NX),NY,NX) = PSI_at_Sautration_vr(NUM_col(NY,NX),NY,NX)
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine UpdateSoilWaterPotential
!------------------------------------------------------------------------------------------
  subroutine SnowSurLitterExch(I,J,dt_SnoHeat,M,L,NY,NX,VapCond1,VapCond2,TCND1W,VLairSno1,PSISV1,TCNDS,&
    TotWatXFlx2SoiMicP,TotHeatFlow2Soi,TotSnoWatFlow2Litr,TotSnoHeatFlow2Litr,SnowLiquidBudget)
  implicit none
  integer  , intent(in) :: I,J
  real(r8) , intent(in) :: dt_SnoHeat
  integer  , intent(in) :: M,L,NY,NX
  real(r8) , intent(in) :: VapCond1,VapCond2
  real(r8) , intent(in) :: TCND1W             !heat conductivity in snow
  real(r8) , intent(in) :: VLairSno1,PSISV1
  real(r8) , intent(in) :: TCNDS              !heat conductivity in soil
  real(r8), intent(inout) :: TotWatXFlx2SoiMicP
  real(r8), intent(inout) :: TotHeatFlow2Soi
  real(r8), intent(inout) :: TotSnoWatFlow2Litr
  real(r8), intent(inout) :: TotSnoHeatFlow2Litr
  real(r8), intent(inout) :: SnowLiquidBudget

  character(len=*), parameter :: subname='SnowSurLitterExch'

  real(r8) :: TK0X,TKXR,TK1X
  real(r8) :: CNVR
  real(r8) :: AvgVaporCondctSnowLitR
  real(r8) :: AvgVaporCondctSoilLitR
  real(r8) :: THETRR                   !volume fraction as litter
  real(r8) :: TCNDR                    !heat conductivity in litter
  real(r8) :: AvgThermCondctSnoLitR    !weighted heat conductivity between snow and litter
  real(r8) :: AvgThermCondctSoilLitR   !weighted heat conductivity between soil and litter
  real(r8) :: ThetaWLitR               !volumetric water content in litter.

  ! ThetaWLitR,THETW1=litter, soil water concentration
  ! VWatLitRHoldCapcity=litter water retention capacity
  call PrintInfo('beg '//subname)
  IF(VLitR_col(NY,NX).LE.ZEROS(NY,NX))return

  IF(VLitR_col(NY,NX).GT.ZEROS(NY,NX).AND.VLWatMicP1_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    ThetaWLitR=AMIN1(VWatLitRHoldCapcity_col(NY,NX),VLWatMicP1_vr(0,NY,NX))/VLitR_col(NY,NX)
    IF(ThetaWLitR.LT.FieldCapacity_vr(0,NY,NX))THEN
      PSISM1_vr(0,NY,NX)=AMAX1(PSIHY,-EXP(LOGPSIFLD_col(NY,NX)+((LOGFldCapacity_vr(0,NY,NX)-LOG(ThetaWLitR))/FCD_vr(0,NY,NX) &
        *LOGPSIMND_col(NY,NX))))
    ELSEIF(ThetaWLitR.LT.POROS0_col(NY,NX))THEN
      PSISM1_vr(0,NY,NX)=-EXP(LOGPSIAtSat(NY,NX)+(((LOGPOROS_vr(0,NY,NX)-LOG(ThetaWLitR))/PSD_vr(0,NY,NX))**SRP_vr(0,NY,NX) &
        *LOGPSIMXD_col(NY,NX)))
    ELSE
      ThetaWLitR=POROS0_col(NY,NX)
      PSISM1_vr(0,NY,NX)=PSI_at_Sautration_vr(0,NY,NX)
    ENDIF
  ELSE
    ThetaWLitR         = POROS0_col(NY,NX)
    PSISM1_vr(0,NY,NX) = PSI_at_Sautration_vr(0,NY,NX)
  ENDIF

  CNVR = safe_adb(VaporDiffusivityLitR_col(NY,NX)*FracAirFilledSoilPoreM_vr(M,0,NY,NX)*POROQ*FracAirFilledSoilPoreM_vr(M,0,NY,NX),POROS_vr(0,NY,NX))

  IF(FracSurfByLitR_col(NY,NX).GT.ZEROL)THEN
    IF(VapCond1.GT.ZERO .AND. CNVR.GT.ZERO)THEN
      AvgVaporCondctSnowLitR=2.0_r8*CNVR*VapCond1/(VapCond1*DLYRR_COL(NY,NX)+CNVR*SnowThickL0_snvr(L,NY,NX))
    ELSE
      AvgVaporCondctSnowLitR=2.0_r8*VapCond1/(DLYRR_COL(NY,NX)+SnowThickL0_snvr(L,NY,NX))
    ENDIF

    IF(CNVR.GT.ZERO.AND.VapCond2.GT.ZERO)THEN
      AvgVaporCondctSoilLitR=2.0_r8*CNVR*VapCond2/(CNVR*DLYR_3D(3,NUM_col(NY,NX),NY,NX)+VapCond2*DLYRR_COL(NY,NX))
    ELSE
      AvgVaporCondctSoilLitR=2.0_r8*VapCond2/(DLYR_3D(3,NUM_col(NY,NX),NY,NX)+DLYRR_COL(NY,NX))
    ENDIF

    THETRR = AZMAX1(1.0_r8-FracAirFilledSoilPore_vr(0,NY,NX)-FracSoiPAsWat_vr(0,NY,NX)-FracSoiPAsIce_vr(0,NY,NX))
    TCNDR  = (0.779_r8*THETRR*9.050E-04_r8+0.622_r8*FracSoiPAsWat_vr(0,NY,NX) &
      *2.067E-03_r8+0.380_r8*FracSoiPAsIce_vr(0,NY,NX)*7.844E-03_r8+FracAirFilledSoilPore_vr(0,NY,NX) &
      *9.050E-05_r8)/(0.779_r8*THETRR+0.622_r8*FracSoiPAsWat_vr(0,NY,NX) &
      +0.380_r8*FracSoiPAsIce_vr(0,NY,NX)+FracAirFilledSoilPore_vr(0,NY,NX))

    IF(TCND1W.GT.ZERO.AND.TCNDR.GT.ZERO)THEN
      AvgThermCondctSnoLitR=2.0_r8*TCND1W*TCNDR/(TCND1W*DLYRR_COL(NY,NX)+TCNDR*SnowThickL0_snvr(L,NY,NX))
    ELSE
      AvgThermCondctSnoLitR=0.0_r8
    ENDIF

    IF(TCNDR.GT.ZERO.AND.TCNDS.GT.ZERO)THEN
      AvgThermCondctSoilLitR=2.0_r8*TCNDR*TCNDS/(TCNDR*DLYR_3D(3,NUM_col(NY,NX),NY,NX)+TCNDS*DLYRR_COL(NY,NX))
    ELSE
      AvgThermCondctSoilLitR=0.0_r8
    ENDIF
  ELSE
    AvgVaporCondctSnowLitR = 0.0_r8
    AvgVaporCondctSoilLitR = 0.0_r8
    AvgThermCondctSnoLitR  = 0.0_r8
    AvgThermCondctSoilLitR = 0.0_r8
  ENDIF

  ! do snow-litter-soil heat-water exchange
  ! SHORTER TIME STEP FOR SURFACE RESIDUE FLUX CALCULATIONS
  call SnowSurfLitRIteration(I,J,dt_SnoHeat,L,M,NY,NX,AvgThermCondctSnoLitR,AvgThermCondctSoilLitR,&
    AvgVaporCondctSoilLitR,AvgVaporCondctSnowLitR,PSISV1,VLairSno1,&
    TotWatXFlx2SoiMicP,TotHeatFlow2Soi,TotSnoWatFlow2Litr,TotSnoHeatFlow2Litr,SnowLiquidBudget)
   
  call PrintInfo('end '//subname)
  end subroutine SnowSurLitterExch

!------------------------------------------------------------------------------------------

  subroutine SnowTopSoilExch(dt_SnoHeat,M,L,NY,NX,VapCond1,VapSnoSrc,VLairSno1,TCND1W,&
    VapFlxSno2Soi1,HeatConvFlxSno2Soi1,HeatCndFlxSno2Soi,VapCond2,PSISV1,TCNDS,SnowLiquidBudget)
  implicit none
  real(r8), intent(in) :: dt_SnoHeat   !time step size for snow iteration, [h]
  integer , intent(in) :: M            !soil heat-flow iteration id
  integer , intent(in) :: L,NY,NX
  real(r8), intent(inout) :: SnowLiquidBudget
  real(r8), intent(in) :: VapCond1,VapSnoSrc
  real(r8), intent(in) :: VLairSno1,TCND1W
  real(r8), intent(out):: VapFlxSno2Soi1,HeatConvFlxSno2Soi1
  real(r8), intent(out) :: HeatCndFlxSno2Soi,VapCond2,PSISV1,TCNDS
  real(r8) :: AvgThermCondctSoilLitR
  real(r8) :: WTHET2,VPY,VapSoiDest,AvgVaporCondctSoilLitR,H2OVapFlx,H2OVapFlxMax
  real(r8) :: HeatCnduct,HeatCnductMax,TKWX1,TKY
  character(len=*), parameter :: subname='SnowTopSoilExch'

  call PrintInfo('beg '//subname)
  !in ATS mode keep the ATS matric potential (see SoilSRFEnerbyBalanceM)
  if(.not.ats_cpl_mode)call UpdateSoilWaterPotential(NY,NX)
  !micropore pressure, excluding gravitational pressure
  PSISV1=PSISM1_vr(NUM_col(NY,NX),NY,NX)+PSISoilOsmotic_vr(NUM_col(NY,NX),NY,NX)
!
  ! VAPOR FLUX BETWEEN SNOWPACK AND SOIL SURFACE
  !
  ! VLairSno1,THETPM=air volume,concentration
  ! VapCond1,VapCond2=vapor conductances of source, destination layers
  ! VapSnoSrc,VapSoiDest=vapor concentrations of source, destination layers
  ! POROS,POROQ=porosity, tortuosity
  ! WVapDifusvitySoil_vr=vapor diffusivity
  ! TKSnow1,TK1=snow,soil surface temperature
  ! PSISV1=soil matric+osmotic potential
  ! AvgVaporCondctSoilLitR=snow-soil vapor conductance
  ! DLYR=soil surface layer depth
  ! H2OVapFlx,H2OVapFlxMax=vapor flux unlimited,limited by vapor
  ! VPY=equilibrium vapor concentration
  ! dts_wat=time step for flux calculations
  ! VapFlxSno2Soi1,HeatConvFlxSno2Soi1=vapor flux and its convective heat flux
  !
  IF(VLairSno1.GT.ZEROS2(NY,NX).AND.FracAirFilledSoilPoreM_vr(M,NUM_col(NY,NX),NY,NX).GT.AirFillPore_Min)THEN
    VapSoiDest = vapsat(TKSoil1_vr(NUM_col(NY,NX),NY,NX))*EXP(18.0_r8*PSISV1/(RGASC*TKSoil1_vr(NUM_col(NY,NX),NY,NX)))
    VapCond2   = WVapDifusvitySoil_vr(NUM_col(NY,NX),NY,NX)*FracAirFilledSoilPoreM_vr(M,NUM_col(NY,NX),NY,NX)*POROQ &
      *FracAirFilledSoilPoreM_vr(M,NUM_col(NY,NX),NY,NX)/POROS_vr(NUM_col(NY,NX),NY,NX)

    AvgVaporCondctSoilLitR = 2.0_r8*VapCond1*VapCond2/(VapCond1*DLYR_3D(3,NUM_col(NY,NX),NY,NX)+VapCond2*SnowThickL0_snvr(L,NY,NX))
    H2OVapFlx              = AvgVaporCondctSoilLitR*(VapSnoSrc-VapSoiDest)*AREA_3D(3,NUM_col(NY,NX),NY,NX)*FracSurfSnowCovered_col(NY,NX)&
      *FracSurfBareSoil_col(NY,NX)*dt_SnoHeat
    VPY          = (VapSnoSrc*VLairSno1+VapSoiDest*VLsoiAirPM_vr(M,NUM_col(NY,NX),NY,NX))/(VLairSno1+VLsoiAirPM_vr(M,NUM_col(NY,NX),NY,NX))
    H2OVapFlxMax = (VapSnoSrc-VPY)*VLairSno1*dts_sno

    if(abs(H2OVapFlx)>1.e10)then
      print*,'hev',vapsat(TKSoil1_vr(NUM_col(NY,NX),NY,NX)),PSISV1,RGASC,TKSoil1_vr(NUM_col(NY,NX),NY,NX)
      print*,'h2ovapflx',AvgVaporCondctSoilLitR,VapSnoSrc,VapSoiDest,FracSurfSnowCovered_col(NY,NX)*FracSurfBareSoil_col(NY,NX)
    endif
    IF(H2OVapFlx.GE.0.0_r8)THEN
      !water/heat flux goes into soil
      VapFlxSno2Soi1      = AZMAX1(AMIN1(H2OVapFlx,H2OVapFlxMax,SnowLiquidBudget))
      HeatConvFlxSno2Soi1 = (cpw*TKSnow1_snvr(L,NY,NX)*HeatAdv_scal+EvapLHTC)*VapFlxSno2Soi1
    ELSE
      !water flux out of soil
      VapFlxSno2Soi1      = AZMIN1(AMAX1(H2OVapFlx,H2OVapFlxMax, &
        -AZMAX1(VLWatMicP1_vr(NUM_col(NY,NX),NY,NX))*dts_sno))
      HeatConvFlxSno2Soi1 = (cpw*TKSoil1_vr(NUM_col(NY,NX),NY,NX)*HeatAdv_scal+EvapLHTC)*VapFlxSno2Soi1
    ENDIF
  ELSE
    VapCond2            = 0.0_r8
    VapFlxSno2Soi1      = 0.0_r8
    HeatConvFlxSno2Soi1 = 0.0_r8
  ENDIF

  ! This exchange debits liquid water; solid sublimation remains in the
  ! atmospheric snow exchange, where its mass and latent heat are explicit.
  SnowLiquidBudget=MAX(0._r8,SnowLiquidBudget-MAX(0._r8,VapFlxSno2Soi1))
  !
  ! HEAT FLUX BETWEEN SNOWPACK AND SURFACE SOIL
  !
  ! WTHET2=multiplier for air concentration in thermal conductivity
  ! TCND1W,TCNDS=thermal conductivity of snowpack, soil surface
  ! STC,DTC=mineral component of thermal conductivity
  ! FracSoiPAsWat,FracSoiPAsIce,FracSoiPAsAir=soil surface water,ice,air concentrations
  ! BARE=soil surface fraction
  ! AvgThermCondctSoilLitR=snowpack-soil thermal conductance
  ! TKWX1=interim snowpack temperature
  ! TKY=equilibrium temperature
  ! HeatCnductMax,HeatCnduct=heat-constrained,heat-unconstrained heat fluxes
  ! dt_SnoHeat=time step for snowpack flux calculations
  ! HeatCndFlxSno2Soi=snowpack-soil heat flux
  ! FracSoilAsAirt=air-filled pores (including macropores and micropores)
  WTHET2 = 1.467_r8-0.467_r8*FracSoilAsAirt(NUM_col(NY,NX),NY,NX)
  TCNDS  = (NumerSolidThermCond_vr(NUM_col(NY,NX),NY,NX)+FracSoiPAsWat_vr(NUM_col(NY,NX),NY,NX) &
    *2.067E-03_r8+0.611_r8*FracSoiPAsIce_vr(NUM_col(NY,NX),NY,NX)*7.844E-03_r8 &
    +WTHET2*FracAirFilledSoilPore_vr(NUM_col(NY,NX),NY,NX)*9.050E-05_r8) &
    /(DenomSolidThermCond_vr(NUM_col(NY,NX),NY,NX)+FracSoiPAsWat_vr(NUM_col(NY,NX),NY,NX) &
    +0.611_r8*FracSoiPAsIce_vr(NUM_col(NY,NX),NY,NX) &
    +WTHET2*FracAirFilledSoilPore_vr(NUM_col(NY,NX),NY,NX))

  !the mean thermal conductivity
  IF(FracSurfBareSoil_col(NY,NX).GT.ZERO)THEN
    AvgThermCondctSoilLitR=2.0_r8*TCND1W*TCNDS/(TCND1W*DLYR_3D(3,NUM_col(NY,NX),NY,NX)+TCNDS*SnowThickL0_snvr(L,NY,NX))
  ELSE
    AvgThermCondctSoilLitR=0.0_r8
  ENDIF

  TKWX1 = TKSoil1_vr(NUM_col(NY,NX),NY,NX)+HeatConvFlxSno2Soi1/VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX)
  TKY   = (TKSnow1_snvr(L,NY,NX)*VLHeatCapSnowM1_snvr(L,NY,NX)+TKWX1*VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX)) &
    /(VLHeatCapSnowM1_snvr(L,NY,NX)+VHeatCapacity1_vr(NUM_col(NY,NX),NY,NX))
  HeatCnductMax = (TKSnow1_snvr(L,NY,NX)-TKY)*VLHeatCapSnowM1_snvr(L,NY,NX)*dts_sno
  HeatCnduct    = AvgThermCondctSoilLitR*(TKSnow1_snvr(L,NY,NX)-TKWX1)*AREA_3D(3,NUM_col(NY,NX),NY,NX) &
    *FracSurfSnowCovered_col(NY,NX)*FracSurfBareSoil_col(NY,NX)*dt_SnoHeat
  IF(HeatCnduct.GE.0.0_r8)THEN
    HeatCndFlxSno2Soi=AZMAX1(AMIN1(HeatCnductMax,HeatCnduct))
  ELSE
    HeatCndFlxSno2Soi=AZMIN1(AMAX1(HeatCnductMax,HeatCnduct))
  ENDIF
  call PrintInfo('beg '//subname)
  end subroutine SnowTopSoilExch
!------------------------------------------------------------------------------------------

end module SnowPhysMod
