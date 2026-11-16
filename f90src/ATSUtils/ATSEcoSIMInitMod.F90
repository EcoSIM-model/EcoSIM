module ATSEcoSIMInitMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilWaterDataType
  use SharedDataMod
  use GridDataType
  use SOMDataType
  USE SoilPhysDataType
  use LandSurfDataType
  use EcoSIMCtrlMod
  use HydroThermData, only : PSISM1_vr, TKSoil1_vr, VHeatCapacity1_vr, &
      SoilFracAsMicP_vr, VLWatMicP1_vr, VLiceMicP1_vr, FracSoiPAsWat_vr, &
      FracSoiPAsIce_vr, FracAirFilledSoilPore_vr!need the only as some vars
  use CanopyDataType, only: RadSWGrnd_col
  use ClimForcDataType, only : LWRadSky_col, TairK_col, &
      VPA_col, WindSpeedAtm_col, RainH, VPK_col
  use SoilPropertyDataType
  use SurfLitterDataType
  use SurfSoilDataType, only: VLWatHeldCapSurf_col
  use EcoSIMConfig
  use MiniMathMod
  use SnowDataType
  use EcoSiMParDataMod, only : micpar
implicit none
  character(len=*), private, parameter :: mod_filename=&
  __FILE__

  public :: THETRX
  real(r8), pointer :: THETRX(:)

  public :: Init_EcoSIM_Soil
  contains

  subroutine Init_EcoSIM_Soil(NYS)
  use EcoSimConst
  use GridMod           , only : SetMeshATS
  use InitAllocMod
  !use InitEcoSIM
  use UptakesMod
  !use Hour1Mod,         only: InitHour1
  use PlantBGCPars
  use PlantInfoMod
  use StartsMod, only : startsim, set_ecosim_solver, InitControlParms
  implicit none
  integer :: NY,NX,L,NHW,NHE,NVN,NVS
  integer, intent(in) :: NYS
  real(r8) :: YSIN(NumOfSkyAzimuthSects),YCOS(NumOfSkyAzimuthSects),SkyAzimuthAngle(NumOfSkyAzimuthSects)
  real(r8) :: DORGC(JZ),DVLiceMicP(JZ)
  real(r8) :: TXCO2(JY,JX),DORGE(JY,JX)
  real(r8) :: VOLISO,VOLPT,VOLTT
  real(r8) :: TFLWT
  real(r8) :: VPS(JY,JX)

  NHW=1;NHE=1;NVN=1;NVS=NYS
  !Setting some variables
  !That ecosim needs to recognize that it is running in coupled mode
  !with ATS and to turn off features unsupported in the coupler
  ATS_cpl_mode     = .true.
  column_mode      = .true.
  plant_model      = pheno_bool
  ldo_sp_mode      = pheno_bool
  microbial_model  = .false.
  soichem_model    = .false.
  snowRedist_model = .false.
  disp_planttrait  = .false.
  disp_modelconfig = .false.
  mod_snow_albedo  = a_bool

  !Set pft file manually for now
  pft_file_in = ecosim_pft_file_path
  !pft_file_in = '/home/ctg/Work/code/ats_dev_dir/amanzi_tpls-build-master-Debug/ecosim/ecosim-0.1.0-source/input_data/ecosim_pftpar_20260303.nc'
  
  
  !Calling some setup functions
  call SetMeshATS(NHW,NVN,NHE,NVS)
  call set_ecosim_solver(30, 10, 20, 20)

  !call InitModules()
  call InitAlloc()
  if(ldo_sp_mode)call ReadPlantTraitTable()
  call InitUptake
  allocate(THETRX(1:micpar%NumOfLitrCmplxs))
  THETRX=real((/4.0E-06,8.0E-06,8.0E-06/),r8)

  !setting a few variables
  FlowDirIndicator_col = 3 !Basically a switch, setting to 3 removes lateral flow
  MaxNumRootLays_col   = 1 !Is the number of layers down the roots go
  NX               = 1

  VOLISO = 0.0_r8
  TFLWT  = 0.0_r8
  VOLPT  = 0.0_r8
  VOLTT  = 0.0_r8

  !This fixes some issues with the initialization of the mesh
  !After this module startsim is run, which runs InitHGrid and
  !InitSoilLayerDepths. These draw their inital conditions from
  !DH_col, DV_col, and CumDepz2LayBottom_vr
  do NY = 1,NYS
    DH_col(NY,NX) = sqrt(column_area(NY))
    DV_col(NY,NX) = sqrt(column_area(NY))
    NU_col(NY,NX)               = a_NU(NY)
    NL_col(NY,NX)               = a_NL(NY)

    do L=NU_col(NY,NX),NL_col(NY,NX)
      CumDepz2LayBottom_vr(L,NY,NX)=a_CumDepz2LayBottom_vr(L,NY)
    enddo
  end do

  do NY=1,NYS
    TXCO2(NY,NX)            = 0.0_r8
    DORGE(NY,NX)            = 0.0_r8
    QRunSurf_col(NY,NX)     = 0.0_r8
    ASP_col(NY,NX)          = a_ASP(NY)
    TairK_col(NY,NX)=tairc(NY) !it's already in K??
    !convert VPA from ATS units (Pa) to EcoSIM (MPa)
    VPK_col(NY,NX)          = vpair(NY)/1.0e3 !vapor pressure in kPa
    VPA_col(NY,NX)              = VPK_col(NY,NX)*2.173E-03_r8/TairK_col(NY,NX)

    WindSpeedAtm_col(NY,NX)  = uwind(NY)*3600.0_r8
    SnowDepth_col(NY,NX) = surf_snow_depth(NY)
    
    !Need to check if litter area is set or not
    if(VGeomLayer_vr(0,NY,NX).EQ.0.0)then
      VGeomLayer_vr(0,NY,NX) = 0.1
    endif
    POROS0_col(NY,NX) = 0.5

    DO L=NU_col(NY,NX),NL_col(NY,NX)
      TKSoil1_vr(L,NY,NX) = a_TEMP(L,NY)
      POROS_vr(L,NY,NX)=a_PORO(L,NY)
      SoiBulkDensityt0_vr(L,NY,NX)=a_BKDSI(L,NY)/1.0e3_r8
      SoilBulkDensity_vr(L,NY,NX)=a_BKDSI(L,NY)/1.0e3_r8
      SoilFracAsMicP_vr(L,NY,NX) = 1.0
          
      !CSoilOrgM_vr(ielmc,L,NY,NX)=a_CORGC(L,NY)
      !CSoilOrgM_vr(ielmn,L,NY,NX)=a_CORGN(L,NY)
      !CSoilOrgM_vr(ielmp,L,NY,NX)=a_CORGP(L,NY)
      
      CSoilOrgM_vr(ielmc,L,NY,NX)=0.0_r8
      CSoilOrgM_vr(ielmn,L,NY,NX)=0.0_r8
      CSoilOrgM_vr(ielmp,L,NY,NX)=0.0_r8

    ENDDO
  ENDDO

  PSIAtFldCapacity_col = pressure_at_field_capacity
  PSIAtWiltPoint_col = pressure_at_wilting_point
  
  ! Coupled startup also needs native litter densities and numerical constants.
  ! Without this, the litter-free volume calculation divides zero by zero.
  call InitControlParms
  call startsim(NHW,NHE,NVN,NVS)
  do NY=1,NYS
    ! Native depression-storage geometry (Hour1Mod:ResetSurfResidualProperty).
    ! Its 0.02 m hydrological roughness is distinct from aerodynamic roughness.
    VLWatHeldCapSurf_col(NY,1)=max(0.001_r8,0.112_r8*0.020_r8 &
      +3.10_r8*0.020_r8**2-0.012_r8*0.020_r8*SLOPE_col(0,NY,1))*column_area(NY)
  enddo

  end subroutine Init_EcoSIM_Soil


end module ATSEcoSIMInitMod
