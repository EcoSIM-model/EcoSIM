module RestartMod
!
! Description:
! code to write restart/check point files
! code to read/write restart files
  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use data_const_mod    , only : spval => DAT_CONST_SPVAL, ispval => DAT_CONST_ISPVAL  
  use EcoSIMConfig      , only : jcplx=> jcplxc, NumMicbAFunGrupsPerCmplx, NumMicbHFunGrupsPerCmplx,nlbiomcp=>NumLiveMicrbCompts
  use EcoSIMConfig      , only : ndbiomcp=>NumDeadMicrbCompts,jsken=>jskenc, is_cold_run
  use EcoSIMConfig      , only : inst_suffix,ref_date,start_date, ctitle, finidat
  use EcoSIMConfig      , only : case_name,hostname,version,source,username
  use EcoSiMParDataMod  , only : micpar,pltpar
  use TracerIDMod       , only : trc_confs
  use EcoSIMCtrlMod     , only : etimer,do_budgets,plant_model,ldo_sp_mode,use_ecosys_drainage,DefaultSoilHydraulicScheme
  use SoilHydraulicSchemeMod, only: GetSoilHydraulicSchemeName,GetSoilHydraulicSchemeRevision
  use ecosim_log_mod   , only : errMsg => shr_log_errMsg
  use abortutils        , only : endrun,destroy,iulog
  use UnitMod           , only : units
  use DebugToolMod      , only : PrintInfo
  use HistFileMod       , only : hist_restart_ncd  
  use NumericalAuxMod
  use EcoSIMCtrlDataType  
  use restUtilMod  
  use ncdio_pio
  use MicrobialDataType
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use GridConsts
  use CanopyRadDataType
  use FlagDataType
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use ClimForcDataType
  use LandSurfDataType
  use PlantTraitDataType
  use PlantDataRateType
  use SnowDataType
  use SurfLitterDataType
  use SurfSoilDataType
  use CanopyDataType
  use SoilBGCDataType
  use PlantMgmtDataType
  use EcoSimSumDataType
  use RootDataType
  use EcosimBGCFluxType
  use EcoSIMHistMod
  use SoilPropertyDataType
  use IrrigationDataType
  use AqueChemDatatype
  use SedimentDataType
  use GridDataType
  use fileUtil
implicit none
  private
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  character(len=8),public :: rest_opt='year'
  integer         ,public :: rest_frq=-999999999
  integer,  parameter :: uninit_int = -999999999
  character(len=16), parameter :: nameg  = 'gridcell'     ! name of gridcells
  character(len=16), parameter :: namet  = 'topounit'     ! name of topographic units
  character(len=16), parameter :: namec  = 'column'       ! name of columns
  character(len=16), parameter :: namep  = 'pft'          ! name of patches

  public :: restFile
  public :: restartnc_soil_hydraulic_scheme
  public :: get_restart_date
  contains

!------------------------------------------------------------------------------------------
  subroutine restFile(flag)
  implicit none
  character(len=*), intent(in) :: flag
  character(len=256) :: fnamer       ! name of netcdf restart file
  character(len=256) :: filer                   ! restart file name
  character(len=256) :: path       ! name of netcdf restart file  
  character(len=18) :: rdate
  integer :: yr,mon,day,tod

  if (flag=='read')then
    call restFile_getfile(fnamer,path)
    call restFile_read( bounds, fnamer)
    ! At midnight the saved daily storage sum belongs to yesterday.
    ! Intraday restarts must retain it; daily setup resumes after this read.
    call etimer%get_curr_date(yr,mon,day,tod)
    if(plant_model.and.tod==0)then
      SeasonalNonstCDayAve_pft(:,bounds%NVN:bounds%NVS,bounds%NHW:bounds%NHE)=0._r8
    endif
  else if(flag=='write')then  
    call etimer%get_curr_date(yr,mon,day,tod)
    write(rdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i6.6)') yr,mon,day,tod
    filer = restFile_filename(rdate=rdate)
    call restFile_write(bounds, filer,rdate=rdate)
  endif
  end subroutine restFile

!------------------------------------------------------------------------------------------
  SUBROUTINE restartnc_plant(ncid,flag)  
!
!     THIS SUBROUTINE WRITES OUT ALL PLANT VARIABLES TO
!     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
!     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
!     TO INITIALIZE LATER MODEL RUNS
!
  use PlantTraitTableMod, only : RootVesselRadius_tab, AlphaVesselAxialResist_tab
  use PlantBGCPars, only : Rax_ref
  use EcoSimConst, only : BndlResistH2OToHeatRatio
  implicit none
  character(len=*) , intent(in) :: flag
  type(file_desc_t), intent(inout) :: ncid ! netcdf id

  integer :: NHW,NHE,NVN,NVS
  integer :: K,L,NX,NY,NZ,N,NB,NR,M,loc,ig
  character(len=64) :: root_cache_name
  integer :: root_conc_varid
  type(var_desc_t) :: root_conc_vardesc
  logical :: has_root_conc, has_water_resistance, has_heat_resistance
  real(r8), target :: datrc_1d(bounds%begc:bounds%endc)
  integer, target  :: datic_1d(bounds%begc:bounds%endc)
  real(r8), target :: datrp_1d(bounds%begp:bounds%endp)
  real(r8), target,allocatable :: datrp_2d(:,:)
  integer , target :: datip_1d(bounds%begp:bounds%endp)
  integer , target :: datip_2d(bounds%begp:bounds%endp,1:MaxNumRootAxes+1)
  real(r8), target, allocatable :: datrp_3d(:,:,:)
  real(r8), target, allocatable :: datrp_4d(:,:,:,:)  
  real(r8), target, allocatable :: datrp_5d(:,:,:,:,:)  
  integer , target, allocatable :: datip_3d(:,:,:)
  integer, pointer :: dat1pr(:),dat2pr(:,:),dat3pr(:,:,:),dat4pr(:,:,:,:)
  real(r8),pointer :: datpr1(:),datpr2(:,:),datpr3(:,:,:)
  real(r8),pointer :: datpr4(:,:,:,:),datpr5(:,:,:,:,:)
  integer :: ncols, npfts
  integer :: ic,ip,sz3,sz4,sz5,sz2
  integer :: routing_varid
  type(var_desc_t) :: routing_vardesc
  logical :: has_category_routing
  integer :: daily_storage_varid
  type(var_desc_t) :: daily_storage_vardesc
  logical :: has_daily_storage
  integer :: senesc_varid
  type(var_desc_t) :: senesc_vardesc
  logical :: has_senesc_snapshot
  integer :: np0_varid
  type(var_desc_t) :: np0_vardesc
  logical :: has_np0
  type(var_desc_t) :: single_grain_vardesc
  logical :: has_single_grain
  integer :: dead_height_varid
  type(var_desc_t) :: dead_height_vardesc
  logical :: has_dead_height
  logical :: retained_inactive_canopy
  integer :: canopy_state_schema
  logical :: retained_inactive_plant
  integer :: plant_state_schema
  integer :: restart_num_plants(JY,JX), restart_active_pft(JP,JY,JX)
  character(len=64) :: single_grain_varname

! execution begins here
  ! Legacy files masked these values in dormant/dead PFTs. New files retain
  ! the thermal/water state and the response used for standing-dead turnover.
  retained_inactive_canopy=.true.
  canopy_state_schema=1
  if(flag=='read')then
    call ncd_inqvid(ncid,'retained_inactive_canopy_state',dead_height_varid,dead_height_vardesc, &
      readvar=retained_inactive_canopy)
  endif
  if(flag/='read' .or. retained_inactive_canopy)then
    call restartvar(ncid,flag,varname='retained_inactive_canopy_state', &
      long_name='inactive canopy water, temperature and turnover response retained',units='none', &
      interpinic_flag='skip',data=canopy_state_schema)
  endif
  ! Dead and fallow slots still carry reserves, residue, dates and properties
  ! used at the next planting. Retain those values, including zero population.
  ! Older files masked this state; keep their filtered read for compatibility.
  retained_inactive_plant=.true.
  plant_state_schema=1
  if(flag=='read')then
    call ncd_inqvid(ncid,'retained_inactive_plant_state',dead_height_varid,dead_height_vardesc, &
      readvar=retained_inactive_plant)
  endif
  if(flag/='read' .or. retained_inactive_plant)then
    call restartvar(ncid,flag,varname='retained_inactive_plant_state', &
      long_name='inactive and fallow plant slot state retained',units='none', &
      interpinic_flag='skip',data=plant_state_schema)
  else
    write(iulog,*) 'Restart masks inactive plant state; exact recovery across replanting is not guaranteed'
  endif
  NHW = bounds%NHW;NVN = bounds%NVN
  NHE = bounds%NHE;NVS = bounds%NVS
  ncols = bounds%ncols
  npfts = bounds%npfts
  sz2=max(MaxNumBranches+1,JZ+1,NumCanopyLayers)
  sz3=max(jsken,NumPlantChemElms,JZ+1,MaxNodesPerBranch+1,NMaxRootSegs)
  sz4=max(MaxNodesPerBranch+1,NumCanopyLayers)
  sz5=max(NumLeafInclinationClasses,NumPlantChemElms)
  allocate(datrp_2d(bounds%begp:bounds%endp,sz2))
  allocate(datrp_3d(bounds%begp:bounds%endp,sz3,sz2))
  allocate(datrp_4d(bounds%begp:bounds%endp,sz4,sz3,sz2))  
  allocate(datrp_5d(bounds%begp:bounds%endp,sz5,sz4,sz3,sz2))  
  allocate(datip_3d(bounds%begp:bounds%endp,sz3,MaxNumBranches+1))

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='NP', dim1name='column',&
       long_name='number of plant species', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,NP_col,datic_1d)
  else 
    !print*,'NP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NP_col,datic_1d)
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='NP', dim1name='column',&
       long_name='number of plant species', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  !PFT slots including fallow ones (NP0_col >= NP_col); the per-PFT copies below cover all of them
  if(flag=='read')then
    call ncd_inqvid(ncid,'NP0',np0_varid,np0_vardesc,readvar=has_np0)
    if(has_np0)then
      dat1pr => datic_1d
      call restartvar(ncid, flag, varname='NP0', dim1name='column',&
         long_name='number of plant species slots, including fallow ones', units='none', interpinic_flag='skip', &
         data=dat1pr, missing_value=ispval, fill_value=ispval)
      call cpcol(flag,NHW,NHE,NVN,NVS,NP0_col,datic_1d)
    else
      !older files have no fallow slots
      NP0_col=NP_col
    endif
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NP0_col,datic_1d)
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='NP0', dim1name='column',&
       long_name='number of plant species slots, including fallow ones', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='LeafStalkAreaAll_col', dim1name='column',&
       long_name='leaf stalk area', units='m2 d-2', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,LeafStalkAreaAll_col,datrc_1d)
  else     
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LeafStalkAreaAll_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='LeafStalkAreaAll_col', dim1name='column',&
       long_name='leaf stalk area', units='m2 d-2', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='RoughnessLength', dim1name='column',&
       long_name='Roughness height', units='m', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,RoughnessLength_col,datrc_1d)
  else     
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RoughnessLength_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='RoughnessLength', dim1name='column',&
       long_name='Roughness height', units='m', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
  endif
  
  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='RootCO2Autor_col', dim1name='column',&
       long_name='Root autotrophic respiraiton', units='gC d-2 h-1', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,RootCO2Autor_col,datrc_1d)
  else     
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RootCO2Autor_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='RootCO2Autor_col', dim1name='column',&
       long_name='Root autotrophic respiraiton', units='gC d-2 h-1', interpinic_flag='skip', &
       data=datpr1, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    dat1pr=>datic_1d
    call restartvar(ncid, flag, varname='NumActivePlants', dim1name='column',&
       long_name='number of active PFT', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,NumActivePlants_col,datic_1d)

  else 
    !print*,'NumActivePlants'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NumActivePlants_col,datic_1d)
    dat1pr=>datic_1d
    call restartvar(ncid, flag, varname='NumActivePlants', dim1name='column',&
       long_name='number of active PFT', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='DATAPI', dim1name='pft',&
       long_name='ID for pft', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,DATAPI,datip_1d)
  else 
    !print*,'DATAPI'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,DATAPI,datip_1d)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='DATAPI', dim1name='pft',&
       long_name='ID for pft', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='IsPlantActive_pft', dim1name='pft',&
       long_name='flag for living pft', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,IsPlantActive_pft,datip_1d)

  else 
    !print*,'IsPlantActive_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,IsPlantActive_pft,datip_1d)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='IsPlantActive_pft', dim1name='pft',&
       long_name='flag for living pft', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iPlantStateLive_pft', dim1name='pft',&
       long_name='flag for species death, 0:live, 1:dead', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iPlantStateLive_pft,datip_1d)
  else 
    !print*,'iPlantStateLive_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iPlantStateLive_pft,datip_1d)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iPlantStateLive_pft', dim1name='pft',&
       long_name='flag for species death, 0:live, 1:dead', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif

  if(flag=='read')then
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='NumCogrowthNode_pft', dim1name='pft',&
       long_name='Number of growing plant nodes', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NumCogrowthNode_pft,datip_1d)
  else 
    !print*,'iPlantStateLive_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NumCogrowthNode_pft,datip_1d)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='NumCogrowthNode_pft', dim1name='pft',&
       long_name='Number of growing plant nodes', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:jsken)
    call restartvar(ncid, flag, varname='StandDeadCompKElms_pft', dim1name='pft',dim2name='nkinecmp',&
       dim3name='elmnts',long_name='standing dead element fraction', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StandDeadCompKElms_pft,datrp_3d)

  else
    !print*,'StandDeadCompKElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StandDeadCompKElms_pft,datrp_3d)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:jsken)    
    call restartvar(ncid, flag, varname='StandDeadCompKElms_pft', dim1name='pft',dim2name='elmnts',&
       dim3name='nkinecmp',long_name='standing dead element fraction', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)

  endif

  ! The copy filters describe saved slots, rather than living plants, in the
  ! complete schema. Do not change the model's actual plant activity flags.
  if(retained_inactive_plant)then
    restart_num_plants=NP0_col
    restart_active_pft=iTrue
  else
    restart_num_plants=NumActivePlants_col
    restart_active_pft=IsPlantActive_pft
  endif

  if(flag=='read')then  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='FireReSet_pft', dim1name='pft',&
       long_name='flag for pft reset due to fire', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FireReSet_pft,datip_1d,NumActivePlants=restart_num_plants)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FireReSet_pft,datip_1d,NumActivePlants=restart_num_plants)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='FireReSet_pft', dim1name='pft',&
       long_name='flag for pft reset due to fire', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif  

  if(flag=='read')then  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iYearPlanting_pft', dim1name='pft',&
       long_name='year of planting', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iYearPlanting_pft,datip_1d,NumActivePlants=restart_num_plants)
  else
    !print*,'iYearPlanting_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iYearPlanting_pft,datip_1d,NumActivePlants=restart_num_plants)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iYearPlanting_pft', dim1name='pft',&
       long_name='year of planting', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif  

  if(flag=='read')then  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iDayPlanting_pft', dim1name='pft',&
       long_name='day of planting', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iDayPlanting_pft,datip_1d,NumActivePlants=restart_num_plants)
  else
    !print*,'iDayPlanting_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iDayPlanting_pft,datip_1d,NumActivePlants=restart_num_plants)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iDayPlanting_pft', dim1name='pft',&
       long_name='day of planting', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif

  if(flag=='read')then  
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='iYearPlantHarvest_pft', dim1name='pft',&
       long_name='year of harvest', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iYearPlantHarvest_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)

  else
    !print*,'iYearPlantHarvest_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iYearPlantHarvest_pft,datip_1d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iYearPlantHarvest_pft', dim1name='pft',&
       long_name='year of harvest', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)    

  endif  

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iDayPlantHarvest_pft', dim1name='pft',&
     long_name='day of harvest', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iDayPlantHarvest_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else 
    !print*,'iDayPlantHarvest_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iDayPlantHarvest_pft,datip_1d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='iDayPlantHarvest_pft', dim1name='pft',&
     long_name='day of harvest', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif 

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NGTopRootLayer_pft', dim1name='pft',&
     long_name='top soil layer that has root', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NGTopRootLayer_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'NGTopRootLayer_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NGTopRootLayer_pft,datip_1d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NGTopRootLayer_pft', dim1name='pft',&
     long_name='top soil layer that has root', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        
  endif  

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='isPlantShootAlive_pft', dim1name='pft',&
     long_name='flag to detect canopy death', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)      
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,isPlantShootAlive_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)

  else
    !print*,'isPlantShootAlive_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,isPlantShootAlive_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='isPlantShootAlive_pft', dim1name='pft',&
     long_name='flag to detect canopy death', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)      
  endif

  if(flag=='read')then    
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='isPlantRootAlive_pft', dim1name='pft',&
     long_name='flag to detect root system death', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,isPlantRootAlive_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'isPlantRootAlive_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,isPlantRootAlive_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='isPlantRootAlive_pft', dim1name='pft',&
     long_name='flag to detect root system death', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        

  endif

  if(flag=='read')then   
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NumOfBranches_pft', dim1name='pft',&
     long_name='branch number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NumOfBranches_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else 
    !print*,'NumOfBranches_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NumOfBranches_pft,datip_1d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NumOfBranches_pft', dim1name='pft',&
     long_name='branch number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)            
  endif 

  if(flag=='read')then     
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='BranchNumber_pft', dim1name='pft',&
     long_name='branch number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,BranchNumber_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else 
    !print*,'BranchNumber_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,BranchNumber_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='BranchNumber_pft', dim1name='pft',&
     long_name='branch number', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)            
  endif  

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MainBranchNum_pft', dim1name='pft',&
     long_name='number of main branch', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MainBranchNum_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else 
    !print*,'MainBranchNum_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MainBranchNum_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MainBranchNum_pft', dim1name='pft',&
     long_name='number of main branch', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)     

  endif     

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='doInitPlant_pft', dim1name='pft',&
     long_name='PFT initialization:0=no,1=yes', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval, flag_values=(/0,1/))     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,doInitPlant_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'doInitPlant_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,doInitPlant_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='doInitPlant_pft', dim1name='pft',&
     long_name='PFT initialization:0=no,1=yes', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval, flag_values=(/0,1/))         
  endif  

  if(flag=='read')then  
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MaxSoilLays4Root_pft', dim1name='pft',&
     long_name='maximum root layer', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval, flag_values=(/0,1/))     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MaxSoilLays4Root_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MaxSoilLays4Root_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='MaxSoilLays4Root_pft', dim1name='pft',&
     long_name='maximum root layer', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval, flag_values=(/0,1/))         
  endif  

  if(flag=='read')then
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NumStructuralRootAxes_pft', dim1name='pft',&
     long_name='number of structural root axes', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NumStructuralRootAxes_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'NRT'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NumStructuralRootAxes_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NumStructuralRootAxes_pft', dim1name='pft',&
     long_name='number of structural root axes', units='none', interpinic_flag='skip', &
     data=dat1pr, missing_value=ispval, fill_value=ispval)    
  endif

  if(flag=='read')then
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NMaxRootBotLayer_pft', dim1name='pft',&
     long_name='Bottom soil layer number that has root', units='none', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, fill_value=ispval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NMaxRootBotLayer_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else 
    !print*,'NMaxRootBotLayer_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NMaxRootBotLayer_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='NMaxRootBotLayer_pft', dim1name='pft',&
     long_name='Bottom soil layer number that has root', units='none', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, fill_value=ispval)     

  endif 

  if(flag=='read')then
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='Myco_pft', dim1name='pft',&
     long_name='mycorrhizal type:0-no,1-yes', units='none', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
     fill_value=ispval,flag_values=(/0,1/))       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Myco_pft,datip_1d,NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'MY'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Myco_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d  
    call restartvar(ncid, flag, varname='Myco_pft', dim1name='pft',&
     long_name='mycorrhizal type:0-no,1-yes', units='none', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
     fill_value=ispval,flag_values=(/0,1/))       
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)    
    call restartvar(ncid, flag, varname='NRoot1stTipLay_raxes', dim1name='pft',dim2name='rootaxs',&
     long_name='maximum soil layer number for root axes', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, &
     fill_value=ispval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NRoot1stTipLay_raxes,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'NRoot1stTipLay_raxes'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NRoot1stTipLay_raxes,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)  
    call restartvar(ncid, flag, varname='NRoot1stTipLay_raxes', dim1name='pft',dim2name='rootaxs',&
     long_name='maximum soil layer number for root axes', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)     
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)  
    call restartvar(ncid, flag, varname='PlantExudElm_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total net root element uptake (+ve) - exudation (-ve)', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantExudElm_CumYr_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PlantExudElm_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantExudElm_CumYr_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)  
    call restartvar(ncid, flag, varname='PlantExudElm_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total net root element uptake (+ve) - exudation (-ve)', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)    
    call restartvar(ncid, flag, varname='LitrfallElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element LitrFall', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LitrfallElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LitrfallElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LitrfallElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)    
    call restartvar(ncid, flag, varname='LitrfallElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element LitrFall', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   

  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantN2Fix_CumYr_pft', dim1name='pft',&
     long_name='total plant N2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantN2Fix_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PlantN2Fix_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantN2Fix_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantN2Fix_CumYr_pft', dim1name='pft',&
     long_name='total plant N2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootUptk_N_CumYr_pft', dim1name='pft',&
     long_name='total plant N uptake', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootUptk_N_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootUptk_N_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootUptk_N_CumYr_pft', dim1name='pft',&
     long_name='total plant N uptake', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootUptk_P_CumYr_pft', dim1name='pft',&
     long_name='total plant P uptake', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootUptk_P_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootUptk_P_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootUptk_P_CumYr_pft', dim1name='pft',&
     long_name='total plant N uptake', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossResp_pft', dim1name='pft',&
     long_name='total plant respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrossResp_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'GrossResp_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrossResp_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossResp_pft', dim1name='pft',&
     long_name='total plant respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CNRTS_pft', dim1name='pft',&
     long_name='root N:C ratio x root growth yield', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CNRTS_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CNRTS_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CNRTS_pft', dim1name='pft',&
     long_name='root N:C ratio x root growth yield', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CPRTS_pft', dim1name='pft',&
     long_name='root P:C ratio x root growth yield', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CPRTS_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CPRTS_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CPRTS_pft', dim1name='pft',&
     long_name='root P:C ratio x root growth yield', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4Uptk_pft', dim1name='pft',&
     long_name='threshold for plant nutrient uptake', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ZERO4Uptk_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ZERO4Uptk_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4Uptk_pft', dim1name='pft',&
     long_name='threshold for plant nutrient uptake', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantO2Stress_pft', dim1name='pft',&
     long_name='plant O2 stress', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantO2Stress_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)

  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantO2Stress_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantO2Stress_pft', dim1name='pft',&
     long_name='plant O2 stress', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4LeafVar_pft', dim1name='pft',&
     long_name='threshold for plant leaf dynamics', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ZERO4LeafVar_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ZERO4LeafVar_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4LeafVar_pft', dim1name='pft',&
     long_name='threshold for plant leaf dynamics', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4Groth_pft', dim1name='pft',&
     long_name='threshold for plant growth', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ZERO4Groth_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ZERO4Groth_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ZERO4Groth_pft', dim1name='pft',&
     long_name='threshold for plant growth', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif


  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracPARads2Canopy_pft', dim1name='pft',&
     long_name='fraction of incoming PAR absorbed by canopy', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FracPARads2Canopy_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FracPARads2Canopy_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracPARads2Canopy_pft', dim1name='pft',&
     long_name='fraction of incoming PAR absorbed by canopy', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracPARads2LiveCanopy_pft', dim1name='pft',&
     long_name='fraction of incoming PAR absorbed by live canopy', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FracPARads2LiveCanopy_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FracPARads2LiveCanopy_pft,datrp_1d, &
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracPARads2LiveCanopy_pft', dim1name='pft',&
     long_name='fraction of incoming PAR absorbed by live canopy', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedTempSens_pft', dim1name='pft',&
     long_name='seed temperature sensitivity', units='1/oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedTempSens_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedTempSens_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedTempSens_pft', dim1name='pft',&
     long_name='seed temperature sensitivity', units='1/oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)         
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ETCanopy_CumYr_pft', dim1name='pft',&
     long_name='total transpiration (<0 into atmosphere)', units='m d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ETCanopy_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'ETCanopy_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ETCanopy_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ETCanopy_CumYr_pft', dim1name='pft',&
     long_name='total transpiration  (<0 into atmosphere)', units='m d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossCO2Fix_pft', dim1name='pft',&
     long_name='total gross CO2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrossCO2Fix_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'GrossCO2Fix_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrossCO2Fix_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='GrossCO2Fix_pft', dim1name='pft',&
     long_name='total gross CO2 fixation', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TdegCCanopy_pft', dim1name='pft',&
     long_name='canopy temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TdegCCanopy_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'TdegCCanopy_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TdegCCanopy_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TdegCCanopy_pft', dim1name='pft',&
     long_name='canopy temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='H2OCuticleResist_pft', dim1name='pft',&
     long_name='Cutitle resistance for H2O', units='h m-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,H2OCuticleResist_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,H2OCuticleResist_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='H2OCuticleResist_pft', dim1name='pft',&
     long_name='Cutitle resistance for H2O', units='h m-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootProteinCMax_pft', dim1name='pft',&
     long_name='Maximum protein concentration in root', units='', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootProteinCMax_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootProteinCMax_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootProteinCMax_pft', dim1name='pft',&
     long_name='Maximum protein concentration in root', units='', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKC_pft', dim1name='pft',&
     long_name='canopy temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    if(retained_inactive_canopy)then
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TKC_pft,datrp_1d)
    else
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TKC_pft,datrp_1d,NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    endif
  else
    !print*,'TKC'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TKC_pft,datrp_1d)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKC_pft', dim1name='pft',&
     long_name='canopy temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TCGroth_pft', dim1name='pft',&
     long_name='canopy growth temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TCGroth_pft,datrp_1d,NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'TCG'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TCGroth_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TCGroth_pft', dim1name='pft',&
     long_name='canopy growth temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKGroth_pft', dim1name='pft',&
     long_name='canopy growth temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TKGroth_pft,datrp_1d,NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'TKG'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TKGroth_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TKGroth_pft', dim1name='pft',&
     long_name='canopy growth temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           

  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='fTCanopyGroth_pft', dim1name='pft',&
     long_name='canopy temperature growth function', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)           
    if(retained_inactive_canopy)then
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,fTCanopyGroth_pft,datrp_1d)
    else
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,fTCanopyGroth_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    endif
  else
    !print*,'fTCanopyGroth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,fTCanopyGroth_pft,datrp_1d)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='fTCanopyGroth_pft', dim1name='pft',&
     long_name='canopy temperature growth function', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)               
  endif  

  if(flag=='read')then  
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopySapwoodC_pft', dim1name='pft',&
     long_name='canopy active stalk C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)               
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopySapwoodC_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopySapwoodC_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopySapwoodC_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopySapwoodC_pft', dim1name='pft',&
     long_name='canopy active stalk C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)               

  endif

  if(flag=='read')then    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyBiomWater_pft', dim1name='pft',&
     long_name='plant canopy water content in dry matter', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    if(retained_inactive_canopy)then
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyBiomWater_pft,datrp_1d)
    else
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyBiomWater_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    endif
  else
    !print*,'CanopyBiomWater_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyBiomWater_pft,datrp_1d)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyBiomWater_pft', dim1name='pft',&
     long_name='canopy water content in dry matter', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then    
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopy_pft', dim1name='pft',&
     long_name='canopy total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
    if(retained_inactive_canopy)then
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanopy_pft,datrp_1d)
    else
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanopy_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    endif
  else
    !print*,'PSICanP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanopy_pft,datrp_1d)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopy_pft', dim1name='pft',&
     long_name='canopy total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopyOsmo_pft', dim1name='pft',&
     long_name='plant canopy osmotic water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
    if(retained_inactive_canopy)then
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanopyOsmo_pft,datrp_1d)
    else
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanopyOsmo_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    endif
  else
    !print*,'PSICanopyOsmo_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanopyOsmo_pft,datrp_1d)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopyOsmo_pft', dim1name='pft',&
     long_name='plant canopy osmotic water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopyTurg_pft', dim1name='pft',&
     long_name='plant canopy turgor water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
    if(retained_inactive_canopy)then
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanopyTurg_pft,datrp_1d)
    else
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanopyTurg_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    endif
  else
    !print*,'PSICanopyTurg_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanopyTurg_pft,datrp_1d)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanopyTurg_pft', dim1name='pft',&
     long_name='plant canopy turgor water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2CuticleResist_pft', dim1name='pft',&
     long_name='CO2 cuticle resistance', units='h m-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CO2CuticleResist_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSICanopyTurg_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CO2CuticleResist_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2CuticleResist_pft', dim1name='pft',&
     long_name='CO2 cuticle resistance', units='h m-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootBiomCPerPlant_pft', dim1name='pft',&
     long_name='root C per plant', units='g g-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootBiomCPerPlant_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)

  else
    !print*,'RootBiomCPerPlant_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootBiomCPerPlant_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='RootBiomCPerPlant_pft', dim1name='pft',&
     long_name='root C per plant', units='g g-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyLeafArea_pft', dim1name='pft',&
     long_name='plant leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafArea_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyLeafArea_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafArea_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyLeafArea_pft', dim1name='pft',&
     long_name='plant leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyMassC_pft', dim1name='pft',&
     long_name='plant canopy biomass, leaf+PetolSheth+stalk mass', units='gC d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyMassC_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyMassC_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyMassC_pft', dim1name='pft',&
     long_name='plant canopy biomass, leaf+PetolSheth+stalk mass', units='gC d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyStemSurfArea_pft', dim1name='pft',&
     long_name='plant stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyStemSurfArea_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyStemSurfArea_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyStemSurfArea_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyStemSurfArea_pft', dim1name='pft',&
     long_name='plant stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif  

  ! Preserve the partial daily mean used by the annual false-break test.
  if(flag=='read')then
    call ncd_inqvid(ncid,'SeasonalNonstCDayAve_pft',daily_storage_varid,daily_storage_vardesc,&
      readvar=has_daily_storage)
    datrp_1d=0._r8
    if(has_daily_storage)then
      datpr1 => datrp_1d
      call restartvar(ncid, flag, varname='SeasonalNonstCDayAve_pft', dim1name='pft',&
        long_name='Accumulated hourly seasonal storage C divided by 24 for the current day', units='gC d-2', &
        interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    else if(plant_model)then
      ! Legacy checkpoints cannot reconstruct the hours preceding an intraday restart.
      write(iulog,*) 'Restart lacks SeasonalNonstCDayAve_pft; using zero. ',&
        'An intraday restart may change the first false-break daily test.'
    endif
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeasonalNonstCDayAve_pft,datrp_1d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeasonalNonstCDayAve_pft,datrp_1d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeasonalNonstCDayAve_pft', dim1name='pft',&
      long_name='Accumulated hourly seasonal storage C divided by 24 for the current day', units='gC d-2', &
      interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='Days4FalseBreak_pft', dim1name='pft',&
     long_name='Accumulated days to singifying false breakout of plants', units='# d-2', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, fill_value=ispval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Days4FalseBreak_pft,datip_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Days4FalseBreak_pft,datip_1d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    dat1pr => datip_1d
    call restartvar(ncid, flag, varname='Days4FalseBreak_pft', dim1name='pft',&
     long_name='Accumulated days to singifying false breakout of plants', units='# d-2', &
     interpinic_flag='skip', data=dat1pr, missing_value=ispval, fill_value=ispval)        
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedCMass_pft', dim1name='pft',&
     long_name='grain C at planting', units='gC seed-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedCMass_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedCMass_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedCMass_pft', dim1name='pft',&
     long_name='grain C at planting', units='gC seed-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPatSeeding_pft', dim1name='pft',&
     long_name='plant seeding population', units='# d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PPatSeeding_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PPatSeeding_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPatSeeding_pft', dim1name='pft',&
     long_name='plant seeding population', units='# d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
  endif
        
  if(flag=='read')then
    datpr1 => datrp_1d   
    call restartvar(ncid, flag, varname='CanopyLeafAreaMAX_pft', dim1name='pft',&
     long_name='Maximum canopy leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafAreaMAX_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafAreaMAX_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyLeafAreaMAX_pft', dim1name='pft',&
     long_name='Maximum canopy leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantPopuLive_pft', dim1name='pft',&
     long_name='plant population', units='# d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantPopuLive_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantPopuLive_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantPopuLive_pft', dim1name='pft',&
     long_name='plant population', units='# d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantPopuDead_pft', dim1name='pft',&
     long_name='Standing dead plant population', units='# d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantPopuDead_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantPopuDead_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PlantPopuDead_pft', dim1name='pft',&
     long_name='Standing dead plant population', units='# d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracGroth2Node_pft', dim1name='pft',&
     long_name='parameter for allocation of growth to nodes', units='', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FracGroth2Node_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FracGroth2Node_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='FracGroth2Node_pft', dim1name='pft',&
     long_name='parameter for allocation of growth to nodes', units='', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmnt_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EcoHavstElmnt_CumYr_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'EcoHavstElmnt_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EcoHavstElmnt_CumYr_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmnt_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='ShootElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant shoot elements', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'EcoHavstElmnt_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='ShootElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant shoot elements', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='Soil2RootMycoExudE_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root exudation', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Soil2RootMycoExudE_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Soil2RootMycoExudE_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='Soil2RootMycoExudE_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root exudation', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmntCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
     call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EcoHavstElmntCum_pft,datrp_2d,NumActivePlants=restart_num_plants,&
       IsPlantActive_pft=restart_active_pft)
  else 
    !print*,'EcoHavstElmntCum_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EcoHavstElmntCum_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmntCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total plant element harvest', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyRespC_CumYr_pft', dim1name='pft',&
     long_name='total autotrophic respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyRespC_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyRespC_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyRespC_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyRespC_CumYr_pft', dim1name='pft',&
     long_name='total autotrophic respiration', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='NetCumElmntFlx2Plant_pft', dim1name='pft',dim2name='elmnts',&
     long_name='effect of canopy element status on seed set', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NetCumElmntFlx2Plant_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'NetCumElmntFlx2Plant_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NetCumElmntFlx2Plant_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='NetCumElmntFlx2Plant_pft', dim1name='pft',dim2name='elmnts',&
     long_name='effect of canopy element status on seed set', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  

  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3Emis_CumYr_pft', dim1name='pft',&
     long_name='total canopy NH3 flux', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NH3Emis_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'NH3Emis_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NH3Emis_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3Emis_CumYr_pft', dim1name='pft',&
     long_name='total canopy NH3 flux', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SurfLitrfalStrutElms_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SurfLitrfalStrutElms_CumYr_pft,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SurfLitrfalStrutElms_CumYr_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SurfLitrfalStrutElms_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LitrfalStrutElms_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LitrfalStrutElms_CumYr_pft,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LitrfalStrutElms_CumYr_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LitrfalStrutElms_CumYr_pft', dim1name='pft',dim2name='elmnts',&
     long_name='total surface LitrFall element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPI_pft',dim1name='pft',&
     long_name='initial plant population over the grid', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PPI_pft,datrp_1d,NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PPI'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PPI_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPI_pft',dim1name='pft',&
     long_name='initial plant population over the grid', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       

  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPX_pft', dim1name='pft',&
     long_name='plant population over the grid', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PPX_pft,datrp_1d,NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PPX'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PPX_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PPX_pft', dim1name='pft',&
     long_name='plant population over the grid', units='# m-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StandDeadStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='standing dead element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StandDeadStrutElms_pft,datrp_2d)
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        DO NZ=1,NP0_col(NY,NX)
          if(any(StandDeadStrutElms_pft(:,NZ,NY,NX)==spval))then
            StandDeadStrutElms_pft(:,NZ,NY,NX)=sum(StandDeadCompKElms_pft(:,:,NZ,NY,NX),dim=2)
            write(iulog,*) 'Restart omitted inactive standing dead pool; reconstructed from components'
          endif
        ENDDO
      ENDDO
    ENDDO
  else
    !print*,'StandDeadStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StandDeadStrutElms_pft,datrp_2d)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StandDeadStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='standing dead element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyHeightLive_pft', dim1name='pft',&
     long_name='canopy height', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyHeightLive_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyHeight'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyHeightLive_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyHeightLive_pft', dim1name='pft',&
     long_name='canopy height', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif

  ! Dead canopy height carries the maximum of earlier live/dead heights.
  ! It cannot be reconstructed from current live height after mortality.
  if(flag=='read')then
    call ncd_inqvid(ncid,'CanopyHeightDead_pft',dead_height_varid,dead_height_vardesc, &
      readvar=has_dead_height)
    if(has_dead_height)then
      datpr1 => datrp_1d
      call restartvar(ncid, flag, varname='CanopyHeightDead_pft', dim1name='pft', &
        long_name='standing dead canopy height', units='m', interpinic_flag='skip', &
        data=datpr1, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyHeightDead_pft,datrp_1d)
    else
      write(iulog,*) 'Restart lacks CanopyHeightDead_pft; exact continuation is not guaranteed'
    endif
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyHeightDead_pft,datrp_1d)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyHeightDead_pft', dim1name='pft', &
      long_name='standing dead canopy height', units='m', interpinic_flag='skip', &
      data=datpr1, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SnowOnCanopy_pft', dim1name='pft',&
     long_name='plant canopy held snow content', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)            
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SnowOnCanopy_pft,datrp_1d)
    where(SnowOnCanopy_pft==spval)SnowOnCanopy_pft=0._r8
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SnowOnCanopy_pft,datrp_1d)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SnowOnCanopy_pft', dim1name='pft',&
     long_name='plant canopy held snow content', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='WatHeldOnCanopy_pft', dim1name='pft',&
     long_name='plant canopy held water content', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)            
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,WatHeldOnCanopy_pft,datrp_1d)
    where(WatHeldOnCanopy_pft==spval)WatHeldOnCanopy_pft=0._r8
  else
    !print*,'WatHeldOnCanopy_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,WatHeldOnCanopy_pft,datrp_1d)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='WatHeldOnCanopy_pft', dim1name='pft',&
     long_name='plant canopy held water content', units='m3 d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
  endif  
  
  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ClumpFactor_pft', dim1name='pft',&
     long_name='clumping factor for self-shading in canopy layer', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ClumpFactor_pft,datrp_1d,NumActivePlants=restart_num_plants,&
    IsPlantActive_pft=restart_active_pft)
  else
    !print*,'ClumpFactor'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ClumpFactor_pft,datrp_1d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ClumpFactor_pft', dim1name='pft',&
     long_name='clumping factor for self-shading in canopy layer', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant CO2 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CO2ByFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CO2ByFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CO2ByFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CO2ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant CO2 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
  endif

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CH4ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant CH4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CH4ByFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else  !define and write
    !print*,'CH4ByFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CH4ByFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CH4ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant CH4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='cNPP_pft', dim1name='pft',&
     long_name='plant cumulative NPP', units='gC d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,cumNPP_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else  !define and write
    !print*,'CH4ByFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,cumNPP_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='cNPP_pft', dim1name='pft',&
     long_name='plant cumulative NPP', units='gC d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='cGPP_pft', dim1name='pft',&
     long_name='plant cumulative GPP', units='gC d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrossCO2Fix_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else  !define and write
    !print*,'CH4ByFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrossCO2Fix_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='cGPP_pft', dim1name='pft',&
     long_name='plant cumulative GPP', units='gC d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)              
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ENGYX_pft', dim1name='pft',&
     long_name='canopy heat storage from previous time step', units='MJ d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ENGYX_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ENGYX_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ENGYX_pft', dim1name='pft',&
     long_name='canopy heat storage from previous time step', units='MJ d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)            

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='O2ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant O2 uptake from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)        
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,O2ByFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'O2ByFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,O2ByFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='O2ByFire_CumYr_pft', dim1name='pft',&
     long_name='plant O2 uptake from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)            

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3byFire_CumYr_pft', dim1name='pft',&
     long_name='plant NH3 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NH3byFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'NH3byFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NH3byFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='NH3byFire_CumYr_pft', dim1name='pft',&
     long_name='plant NH3 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='N2ObyFire_CumYr_pft', dim1name='pft',&
     long_name='plant N2O emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,N2ObyFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'N2ObyFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,N2ObyFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='N2ObyFire_CumYr_pft', dim1name='pft',&
     long_name='plant N2O emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PO4byFire_CumYr_pft', dim1name='pft',&
     long_name='plant PO4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PO4byFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PO4byFire_CumYr_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PO4byFire_CumYr_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PO4byFire_CumYr_pft', dim1name='pft',&
     long_name='plant PO4 emission from fire', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  
  
  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='dReproNodeNumNormByMatG_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant repro node normalized by maturity group', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,dReproNodeNumNormByMatG_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'GROUP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,dReproNodeNumNormByMatG_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='dReproNodeNumNormByMatG_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant repro node normalized by maturity group', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='MatureGroup_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant maturity group', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MatureGroup_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'GROUP'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MatureGroup_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='MatureGroup_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant maturity group', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootNodeNum_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootNodeNum_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'ShootNodeNum_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootNodeNum_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootNodeNum_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootNodeNumAtInitFloral_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at floral initiation', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootNodeNumAtInitFloral_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'ShootNodeNumAtInitFloral_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootNodeNumAtInitFloral_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootNodeNumAtInitFloral_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at floral initiation', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootNodeNumAtAnthesis_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at anthesis', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootNodeNumAtAnthesis_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'ShootNodeNumAtAnthesis_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootNodeNumAtAnthesis_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootNodeNumAtAnthesis_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at anthesis', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NumOfLeaves_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NumOfLeaves_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'NumOfLeaves_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NumOfLeaves_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NumOfLeaves_brch', dim1name='pft',dim2name='nbranches',&
     long_name='node number at anthesis', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafNumberAtFloralInit_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number at floral initiation', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafNumberAtFloralInit_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'VSTGX'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafNumberAtFloralInit_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafNumberAtFloralInit_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number at floral initiation', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NodeNumNormByMatgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during vegetative growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NodeNumNormByMatgrp_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'NodeNumNormByMatgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NodeNumNormByMatgrp_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='NodeNumNormByMatgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during vegetative growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ReprodNodeNumNormByMatrgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during reproductive growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ReprodNodeNumNormByMatrgrp_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'ReprodNodeNumNormByMatrgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ReprodNodeNumNormByMatrgrp_brch,datrp_2d,&
      NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ReprodNodeNumNormByMatrgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during reproductive growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='TotalNodeNumNormByMatgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during vegetative growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TotalNodeNumNormByMatgrp_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'TotalNodeNumNormByMatgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TotalNodeNumNormByMatgrp_brch,datrp_2d,&
      NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='TotalNodeNumNormByMatgrp_brch', dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during vegetative growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='TotReproNodeNumNormByMatrgrp_brch', &
     dim1name='pft',dim2name='nbranches',&
     long_name='normalized node number during reproductive growth stages', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TotReproNodeNumNormByMatrgrp_brch,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'TotReproNodeNumNormByMatrgrp_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TotReproNodeNumNormByMatrgrp_brch,&
      datrp_2d,NumActivePlants=restart_num_plants, &
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='TotReproNodeNumNormByMatrgrp_brch', &
      dim1name='pft',dim2name='nbranches',&
      long_name='normalized node number during reproductive growth stages', units='none', &
      interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4Leafout_brch', dim1name='pft',dim2name='nbranches',&
     long_name='heat requirement for spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4Leafout_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'Hours4Leafout_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4Leafout_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4Leafout_brch', dim1name='pft',dim2name='nbranches',&
     long_name='heat requirement for spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LeafOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='cold requirement for autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4LeafOff_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'Hours4LeafOff_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4LeafOff_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LeafOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='cold requirement for autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LenthenPhotoPeriod_brch', dim1name='pft',dim2name='nbranches',&
     long_name='initial heat requirement for spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4LenthenPhotoPeriod_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'Hours4LenthenPhotoPeriod_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4LenthenPhotoPeriod_brch,datrp_2d,&
      NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LenthenPhotoPeriod_brch', dim1name='pft',dim2name='nbranches',&
     long_name='initial heat requirement for spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4ShortenPhotoPeriod_brch', dim1name='pft',dim2name='nbranches',&
     long_name='initial cold requirement for autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4ShortenPhotoPeriod_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'Hours4ShortenPhotoPeriod_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4ShortenPhotoPeriod_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4ShortenPhotoPeriod_brch', dim1name='pft',dim2name='nbranches',&
     long_name='initial cold requirement for autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours2LeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='counter for mobilizing nonstructural C during spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours2LeafOut_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'Hours2LeafOut_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours2LeafOut_brch,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours2LeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='counter for mobilizing nonstructural C during spring leafout/dehardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HourFailGrainFill_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect physiological maturity from  grain fill', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HourFailGrainFill_brch,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'HourFailGrainFill_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HourFailGrainFill_brch,&
      datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HourFailGrainFill_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect physiological maturity from  grain fill', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HoursDoingRemob_brch', dim1name='pft',dim2name='nbranches',&
     long_name='counter for mobilizing nonstructural C during autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HoursDoingRemob_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'FLGZ'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HoursDoingRemob_brch,datrp_2d,&
      NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HoursDoingRemob_brch', dim1name='pft',dim2name='nbranches',&
     long_name='counter for mobilizing nonstructural C during autumn leafoff/hardening', units='h', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  


  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doInitLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,doInitLeafOut_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'doInitLeafOut_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,doInitLeafOut_brch,datip_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doInitLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)

  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='EnablePlantLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EnablePlantLeafOut_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'EnablePlantLeafOut_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EnablePlantLeafOut_brch,datip_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='EnablePlantLeafOut_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeaveOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,doPlantLeaveOff_brch,datip_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'doPlantLeaveOff_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,doPlantLeaveOff_brch,datip_2d,&
      NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='doPlantLeaveOff_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Prep4Literfall_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Prep4Literfall_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'Prep4Literfall_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Prep4Literfall_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Prep4Literfall_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LiterfalAftMature_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='h', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4LiterfalAftMature_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'Hours4LiterfalAftMature_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Hours4LiterfalAftMature_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='Hours4LiterfalAftMature_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch phenology flag', units='h', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='isPlantBranchAlive_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect branch death', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,isPlantBranchAlive_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'isPlantBranchAlive_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,isPlantBranchAlive_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='isPlantBranchAlive_brch', dim1name='pft',dim2name='nbranches',&
     long_name='flag to detect branch death', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='BranchNumerID_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch number', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,BranchNumerID_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'BranchNumerID_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,BranchNumerID_brch,datip_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='BranchNumerID_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch number', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLeafNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,KLeafNumber_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'KLEAF'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,KLeafNumber_brch,datip_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLeafNumber_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf number', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KVSTG', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,KHiestGroLeafNode_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'KVSTG'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,KHiestGroLeafNode_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KVSTG', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLowestGroLeafNode_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,KLowestGroLeafNode_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'KLowestGroLeafNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,KLowestGroLeafNode_brch,datip_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat2pr => datip_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='KLowestGroLeafNode_brch', dim1name='pft',dim2name='nbranches',&
     long_name='leaf growth stage counter', units='none', &
     interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    dat3pr => datip_3d(1:npfts,1:NumGrowthStages,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='iPlantCalendar_brch', dim1name='pft',dim2name='ngrstages', &
     dim3name='nbranches',long_name='plant growth stage', units='none', &
     interpinic_flag='skip', data=dat3pr, missing_value=ispval, fill_value=ispval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iPlantCalendar_brch,datip_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'iPlantCalendar'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,iPlantCalendar_brch,datip_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    dat3pr => datip_3d(1:npfts,1:NumGrowthStages,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='iPlantCalendar_brch', dim1name='pft',dim2name='ngrstages', &
     dim3name='nbranches',long_name='plant growth stage', units='none', &
     interpinic_flag='skip', data=dat3pr, missing_value=ispval, fill_value=ispval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='CanopyNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNonstElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyNonstElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNonstElms_pft,&
      datrp_2d,NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='CanopyNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif 

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='CanopyNodulNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy nodule nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNodulNonstElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyNodulNonstElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNodulNonstElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='CanopyNodulNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy nodule nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='VHeatCapCanopy_pft', dim1name='pft',&
     long_name='Volumetric canopy heat capacity', units='MJ d-2 K-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,VHeatCapCanopy_pft,datrp_1d,NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'VHeatCapCanopy_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,VHeatCapCanopy_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='VHeatCapCanopy_pft', dim1name='pft',&
     long_name='Volumetric canopy heat capacity', units='MJ d-2 K-1', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='DeltaTKC_pft', dim1name='pft',&
     long_name='change in canopy temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,DeltaTKC_pft,datrp_1d,NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
  else
    !print*,'DeltaTKC_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,DeltaTKC_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='DeltaTKC_pft', dim1name='pft',&
     long_name='change in canopy temperature', units='K', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyLeafSheathC_pft', dim1name='pft',&
     long_name='canopy leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafSheathC_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyLeafSheathC_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafSheathC_pft,datrp_1d, &
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='CanopyLeafSheathC_pft', dim1name='pft',&
     long_name='canopy leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='LeafStalkAreaAct_pft', dim1name='pft',&
     long_name='canopy leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafStalkAreaAct_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafStalkAreaAct_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafStalkAreaAct_pft,datrp_1d, &
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='LeafStalkAreaAct_pft', dim1name='pft',&
     long_name='canopy leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='MatureGroup_pft', dim1name='pft',&
     long_name='Plant maturity group', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MatureGroup_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafStalkAreaAct_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MatureGroup_pft,datrp_1d, &
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='MatureGroup_pft', dim1name='pft',&
     long_name='Plant maturity group', units='none', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanPDailyMin', dim1name='pft',&
     long_name='Daily minimum canopy water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanPDailyMin_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafStalkAreaAct_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSICanPDailyMin_pft,datrp_1d, &
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='PSICanPDailyMin', dim1name='pft',&
     long_name='Daily minimum canopy water potential', units='MPa', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='AvgCanopyBiomC2Graze_pft', dim1name='pft',&
     long_name='Mean canopy C for grazing', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,AvgCanopyBiomC2Graze_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafStalkAreaAct_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,AvgCanopyBiomC2Graze_pft,datrp_1d, &
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='AvgCanopyBiomC2Graze_pft', dim1name='pft',&
     long_name='Mean canopy C for grazing', units='g d-2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'RootElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='PlantElmBalCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='Cumulative plant element conservation error', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    if(retained_inactive_canopy)then
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantElmBalCum_pft,datrp_2d)
    else
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantElmBalCum_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    endif
  else
    !print*,'RootElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PlantElmBalCum_pft,datrp_2d)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='PlantElmBalCum_pft', dim1name='pft',dim2name='elmnts',&
     long_name='Cumulative plant element conservation error', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LeafStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafStrutElms_pft,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LeafStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='PetolShethStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolShethStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PetolShethStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolShethStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='PetolShethStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StalkStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy stalk element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'StalkStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StalkStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy stalk element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StalkRsrvElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy reserve element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkRsrvElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'StalkRsrvElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkRsrvElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='StalkRsrvElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy reserve element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='HuskStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy husk element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HuskStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'HuskStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HuskStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='HuskStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy husk element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EarStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy ear element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EarStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'EarStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EarStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EarStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy ear element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
    call restartvar(ncid, flag, varname='RootCRRadius0_rpvr', dim1name='pft',dim2name='levsoi',dim3name='rootaxs',&
     long_name='Initial radius for root that may undergo secondary growth', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootCRRadius0_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)

  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootCRRadius0_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
    call restartvar(ncid, flag, varname='RootCRRadius0_rpvr', dim1name='pft',dim2name='levsoi',dim3name='rootaxs',&
     long_name='Initial radius for root that may undergo secondary growth', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
    call restartvar(ncid, flag, varname='RootMediumXNum_rpvr', dim1name='pft',dim2name='levsoi',dim3name='rootaxs',&
     long_name='Medium root axes number in soil layer', units='# d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumXNum_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumXNum_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
    call restartvar(ncid, flag, varname='RootMediumXNum_rpvr', dim1name='pft',dim2name='levsoi',dim3name='rootaxs',&
     long_name='Medium root axes number in soil layer', units='# d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
    call restartvar(ncid, flag, varname='Root1stRadius_rpvr', dim1name='pft',dim2name='levsoi',dim3name='rootaxs',&
     long_name='Current radius for root that may undergo secondary growth', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stRadius_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)

  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stRadius_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
    call restartvar(ncid, flag, varname='Root1stRadius_rpvr', dim1name='pft',dim2name='levsoi',dim3name='rootaxs',&
     long_name='Current radius for root that may undergo secondary growth', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='GrainStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrainStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'GrainStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrainStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='GrainStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='canopy grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='ShootNoduleElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='Shoot total nodule mass', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootNoduleElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'ShootNoduleElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootNoduleElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='ShootNoduleElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='Shoot total nodule mass', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootNoduleElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='root total nodule mass', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNoduleElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'RootNoduleElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNoduleElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootNoduleElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='root total nodule mass', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SeasonalNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant stored nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeasonalNonstElms_pft,datrp_2d,NumActivePlants=restart_num_plants)
  else
    !print*,'NonStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeasonalNonstElms_pft,datrp_2d,NumActivePlants=restart_num_plants)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='SeasonalNonstElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant stored nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='HypocotHeight_pft', dim1name='pft',&
     long_name='Hypocotyledon height', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HypocotHeight_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'HypocotHeight_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HypocotHeight_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='HypocotHeight_pft', dim1name='pft',&
     long_name='Hypocotyledon height', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     

  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedDepth_pft', dim1name='pft',&
     long_name='seeding depth', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedDepth_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedDepth_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedDepth_pft', dim1name='pft',&
     long_name='seeding depth', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedMeanLen_pft', dim1name='pft',&
     long_name='mean seed length', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedMeanLen_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedMeanLen_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedMeanLen_pft', dim1name='pft',&
     long_name='mean seed length', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedVolumeMean_pft', dim1name='pft',&
     long_name='mean seed volume', units='m3', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedVolumeMean_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedVolumeMean_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedVolumeMean_pft', dim1name='pft',&
     long_name='mean seed volume', units='m3', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedAreaMean_pft', dim1name='pft',&
     long_name='mean seed area', units='m2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedAreaMean_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SeedAreaMean_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='SeedAreaMean_pft', dim1name='pft',&
     long_name='mean seed area', units='m2', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TC4LeafOut_pft', dim1name='pft',&
     long_name='leaf out temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TC4LeafOut_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TC4LeafOut_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TC4LeafOut_pft', dim1name='pft',&
     long_name='leaf out temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TC4LeafOff_pft', dim1name='pft',&
     long_name='leaf off temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TC4LeafOff_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SeedDepth_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,TC4LeafOff_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='TC4LeafOff_pft', dim1name='pft',&
     long_name='leaf off temperature', units='oC', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='HoursTooLowPsiCan_pft', dim1name='pft',&
     long_name='canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY)', units='h', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HoursTooLowPsiCan_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'HoursTooLowPsiCan_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HoursTooLowPsiCan_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='HoursTooLowPsiCan_pft', dim1name='pft',&
     long_name='canopy plant water stress indicator, number of hours PSICanopy_pft(< PSILY)', units='h', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ChillHours_pft', dim1name='pft',&
     long_name='chilling effect on CO2 fixation', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ChillHours_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CHILL'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ChillHours_pft,datrp_1d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr1 => datrp_1d
    call restartvar(ncid, flag, varname='ChillHours_pft', dim1name='pft',&
     long_name='chilling effect on CO2 fixation', units='m', &
     interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root structural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'RootStrutElms_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootStrutElms_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='RootStrutElms_pft', dim1name='pft',dim2name='elmnts',&
     long_name='plant root structural element', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SapwoodBiomassC_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch active stalk C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SapwoodBiomassC_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SapwoodBiomassC_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SapwoodBiomassC_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SapwoodBiomassC_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch active stalk C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyLeafSheathC_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant branch leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafSheathC_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyLeafSheathC_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafSheathC_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyLeafSheathC_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant branch leaf + sheath C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNonstElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNonstElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyNonstElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNonstElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNonstElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nonstructural element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNodulNonstElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nodule nonstructural C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNodulNonstElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyNodulNonstElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNodulNonstElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNodulNonstElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nodule nonstructural C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch shoot element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'ShootElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,ShootElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ShootElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch shoot element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNodulStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nodule element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNodulStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyNodulStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyNodulStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyNodulStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch nodule element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetolShethStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolShethStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PetolShethStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolShethStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetolShethStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch stalk element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'StalkStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch stalk element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkRsrvElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch reserve element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkRsrvElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'StalkRsrvElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkRsrvElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkRsrvElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch reserve element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HuskStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch husk element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)      
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HuskStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'HuskStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,HuskStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='HuskStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch husk element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='EarStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch ear element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EarStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'EarStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,EarStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='EarStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch ear element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='GrainStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrainStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'GrainStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GrainStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='GrainStrutElms_brch', dim1name='pft',dim2name='elmnts',&
      dim3name='nbranches',long_name='branch grain element', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PotentialSeedSites_brch', dim1name='pft',&
      dim2name='nbranches',long_name='branch seed sites', units=' d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)    
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PotentialSeedSites_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PotentialSeedSites_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PotentialSeedSites_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PotentialSeedSites_brch', dim1name='pft',&
      dim2name='nbranches',long_name='branch seed sites', units=' d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SetNumberSeeds_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch grain number', units='# d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SetNumberSeeds_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SetNumberSeeds_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SetNumberSeeds_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SetNumberSeeds_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch grain number', units='# d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    ! Accept the old restart field name without changing its per-grain values.
    single_grain_varname = 'SingleGrainMeanBiomC_brch'
    call check_var(ncid,trim(single_grain_varname),single_grain_vardesc,has_single_grain,print_err=.false.)
    if(.not.has_single_grain) single_grain_varname = 'GrainSeedBiomCMean_brch'
    call restartvar(ncid, flag, varname=trim(single_grain_varname), dim1name='pft',dim2name='nbranches',&
     long_name='Potential carbon mass per grain', units='gC seed-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SingleGrainMeanBiomC_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SingleGrainMeanBiomC_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SingleGrainMeanBiomC_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SingleGrainMeanBiomC_brch', dim1name='pft',dim2name='nbranches',&
     long_name='Potential carbon mass per grain', units='gC seed-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)

  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaLive_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafAreaLive_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafAreaLive_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafAreaLive_brch,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaLive_brch', dim1name='pft',dim2name='nbranches',&
     long_name='plant branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

!  if(flag=='read')then
!    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
!    call restartvar(ncid, flag, varname='LeafChemElmRemob_brch', dim1name='pft',dim2name='elmnts', &
!     dim3name='nbranches',long_name='branch leaf structural element', units='g d-2', &
!     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
!    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafChemElmRemob_brch,datrp_3d,NumActivePlants=restart_num_plants,&
!      IsPlantActive_pft=restart_active_pft)
!  else
    !print*,'LeafChemElmRemob_brch'
!    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafChemElmRemob_brch,datrp_3d,NumActivePlants=restart_num_plants,&
!      IsPlantActive_pft=restart_active_pft)
!    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
!    call restartvar(ncid, flag, varname='LeafChemElmRemob_brch', dim1name='pft',dim2name='elmnts', &
!     dim3name='nbranches',long_name='branch leaf structural element', units='g d-2', &
!     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
!  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaDying_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafAreaDying_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafAreaDying_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafAreaDying_brch,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaDying_brch', dim1name='pft',dim2name='nbranches',&
     long_name='branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanPBranchHeight', dim1name='pft',dim2name='nbranches',&
     long_name='branch height', units='m', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanPBranchHeight,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanPBranchHeight'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanPBranchHeight,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanPBranchHeight', dim1name='pft',dim2name='nbranches',&
     long_name='branch height', units='m', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  ! Senescence snapshots are needed between leaf appearances, including across restarts.
  if(plant_model)then
    if(flag=='read')then
      call ncd_inqvid(ncid,'LeafSenescInitialElms_brch',senesc_varid,senesc_vardesc,readvar=has_senesc_snapshot)
      datrp_3d=-1._r8
      if(has_senesc_snapshot)then
        datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
        call restartvar(ncid, flag, varname='LeafSenescInitialElms_brch', dim1name='pft',dim2name='elmnts',&
          dim3name='nbranches',long_name='Initial senescing leaf C/N/P mass; -1 means unset', units='g d-2', &
          interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      else
        ! Old files have no initial mass snapshot. Capture current tissue at the next senescence call.
        write(iulog,*) 'Restart lacks LeafSenescInitialElms_brch; senescence will restart from current tissue.'
      endif
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafSenescInitialElms_brch,datrp_3d,&
        NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafSenescInitialElms_brch,datrp_3d,&
        NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
      call restartvar(ncid, flag, varname='LeafSenescInitialElms_brch', dim1name='pft',dim2name='elmnts',&
        dim3name='nbranches',long_name='Initial senescing leaf C/N/P mass; -1 means unset', units='g d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif
    if(flag=='read')then
      call ncd_inqvid(ncid,'PetolSenescInitialElms_brch',senesc_varid,senesc_vardesc,readvar=has_senesc_snapshot)
      datrp_3d=-1._r8
      if(has_senesc_snapshot)then
        datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
        call restartvar(ncid, flag, varname='PetolSenescInitialElms_brch', dim1name='pft',dim2name='elmnts',&
          dim3name='nbranches',long_name='Initial senescing sheath/petiole C/N/P mass; -1 means unset', units='g d-2', &
          interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      else
        ! Old files have no initial mass snapshot. Capture current tissue at the next senescence call.
        write(iulog,*) 'Restart lacks PetolSenescInitialElms_brch; senescence will restart from current tissue.'
      endif
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolSenescInitialElms_brch,datrp_3d,&
        NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolSenescInitialElms_brch,datrp_3d,&
        NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
      call restartvar(ncid, flag, varname='PetolSenescInitialElms_brch', dim1name='pft',dim2name='elmnts',&
        dim3name='nbranches',long_name='Initial senescing sheath/petiole C/N/P mass; -1 means unset', units='g d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmntRemobFlx_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='Cached remobilizable leaf element mass', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafElmntRemobFlx_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafElmntRemobFlx_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafElmntRemobFlx_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmntRemobFlx_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='Cached remobilizable leaf element mass', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetolShethChemElmRemobFlx_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='Cached remobilizable sheath element mass', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolShethChemElmRemobFlx_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PetolShethChemElmRemobFlx_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolShethChemElmRemobFlx_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetolShethChemElmRemobFlx_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='Cached remobilizable sheath element mass', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SenecStalkStrutElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='branch stalk structural C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SenecStalkStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'SenecStalkStrutElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SenecStalkStrutElms_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumPlantChemElms,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='SenecStalkStrutElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nbranches',long_name='branch stalk structural C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyLeafAreaZ_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafAreaZ_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyLeafAreaZ_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafAreaZ_pft,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyLeafAreaZ_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif


  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyLeafCLyr_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafCLyr_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyLeafAreaZ_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafCLyr_pft,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyLeafCLyr_pft', dim1name='pft',dim2name='levcan',&
     long_name='pft canopy layer leaf C', units='g d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  !Like the other root-category state, routing exists only with the plant model.
  !Skip define/read/write when disabled: rootyps can be zero (NetCDF unlimited).
  if(plant_model)then
    if(flag=='read')then
      call ncd_inqvid(ncid,'RootFineFrac2Med_rpvr',routing_varid,routing_vardesc,readvar=has_category_routing)
      if(has_category_routing)then
        datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
        call restartvar(ncid, flag, varname='RootFineFrac2Med_rpvr', dim1name='pft',dim2name='rootyps',&
          dim3name='levsoi',long_name='Fraction of fine-root axes attached to medium roots by root category', units='-', &
          interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
        call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootFineFrac2Med_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
          IsPlantActive_pft=restart_active_pft)
      else
        ! Older restarts store one shared fraction. Retain it for both categories
        ! until the next root update recomputes the count-weighted fractions.
        datpr2 => datrp_2d(1:npfts,1:JZ)
        call restartvar(ncid, flag, varname='RootFineFrac2Med_pvr', dim1name='pft',dim2name='levsoi',&
          long_name='Fraction of fine roots attached to medium roots', units='-', &
          interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
        DO N=1,pltpar%jroots
          call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootFineFrac2Med_rpvr(N,:,:,:,:),datrp_2d,&
            NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
        ENDDO
      endif
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootFineFrac2Med_rpvr,datrp_3d,&
        NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootFineFrac2Med_rpvr', dim1name='pft',dim2name='rootyps',&
        dim3name='levsoi',long_name='Fraction of fine-root axes attached to medium roots by root category', units='-', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:JZ)
    call restartvar(ncid, flag, varname='RootMediumLength_pvr', dim1name='pft',dim2name='levsoi',&
     long_name='Medium size root mean length', units='m', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumLength_pvr,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumLength_pvr,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:JZ)
    call restartvar(ncid, flag, varname='RootMediumLength_pvr', dim1name='pft',dim2name='levsoi',&
     long_name='Medium size root mean length', units='m', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:JZ)
    call restartvar(ncid, flag, varname='CRootLumenArea_pvr', dim1name='pft',dim2name='levsoi',&
     long_name='Coarse root lumen area for xylem transport', units='-', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CRootLumenArea_pvr,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CRootLumenArea_pvr,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:JZ)
    call restartvar(ncid, flag, varname='CRootLumenArea_pvr', dim1name='pft',dim2name='levsoi',&
     long_name='Coarse root lumen area for xylem transport', units='-', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif


  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:JZ)
    call restartvar(ncid, flag, varname='MRootLumenArea_pvr', dim1name='pft',dim2name='levsoi',&
     long_name='Medium root lumen area for xylem transport', units='-', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MRootLumenArea_pvr,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,MRootLumenArea_pvr,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:JZ)
    call restartvar(ncid, flag, varname='MRootLumenArea_pvr', dim1name='pft',dim2name='levsoi',&
     long_name='Medium root lumen area for xylem transport', units='-', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:JZ)
    call restartvar(ncid, flag, varname='GroSrcRootStress_pvr', dim1name='pft',dim2name='levsoi',&
     long_name='Root growth stress due to nutrient and water', units='-', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GroSrcRootStress_pvr,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,GroSrcRootStress_pvr,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:JZ)
    call restartvar(ncid, flag, varname='GroSrcRootStress_pvr', dim1name='pft',dim2name='levsoi',&
     long_name='Root growth stress due to nutrient and water', units='-', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:NumCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyStemSurfAreaZ_pft', dim1name='pft',dim2name='levcan',&
     long_name='plant canopy layer stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyStemSurfAreaZ_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyStemSurfAreaZ_pft'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyStemSurfAreaZ_pft,datrp_2d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:NumCanopyLayers)
    call restartvar(ncid, flag, varname='CanopyStemSurfAreaZ_pft', dim1name='pft',dim2name='levcan',&
     long_name='plant canopy layer stem area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CPOOL4_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf nonstructural C4 content in C4 photosynthesi', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CPOOL4_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CPOOL4_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CPOOL4_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CPOOL4_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf nonstructural C4 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CPOOL3_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf nonstructural C3 content in C3 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CPOOL3_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CPOOL3_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CPOOL3_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CPOOL3_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='leaf nonstructural C3 content in C3 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CMassCO2BundleSheath_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='bundle sheath nonstructural C3 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CMassCO2BundleSheath_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CMassCO2BundleSheath_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CMassCO2BundleSheath_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CMassCO2BundleSheath_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='bundle sheath nonstructural C3 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)     
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CMassHCO3BundleSheath_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='bundle sheath nonstructural C3 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)     
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CMassHCO3BundleSheath_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CMassHCO3BundleSheath_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CMassHCO3BundleSheath_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CMassHCO3BundleSheath_node', dim1name='pft',dim2name='nodes',&
     dim3name='nbranches',long_name='bundle sheath nonstructural C3 content in C4 photosynthesis', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ARLF', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafArea_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'ARLF'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafArea_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='ARLF', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafElmntNode_brch,datrp_4d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafElmntNode_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafElmntNode_brch,datrp_4d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='leaf element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval) 
  endif  

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafProteinC_node', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='layer leaf protein C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafProteinC_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafProteinC_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafProteinC_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafProteinC_node', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='layer leaf protein C', units='g d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif

  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleLength_node', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='plant sheath height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetoleLength_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PetoleLength_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetoleLength_node,datrp_3d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleLength_node', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='plant sheath height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
  endif 

  if(flag=='read')then
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetolShethElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolShethElmntNode_brch,datrp_4d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'WGSHE'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetolShethElmntNode_brch,datrp_4d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetolShethElmntNode_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='sheath element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)      
  endif  
 
  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleProteinC_node', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='node sheath protein C', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetoleProteinC_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'WSSHE'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PetoleProteinC_node,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='PetoleProteinC_node', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='node sheath protein C', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)       
  endif  

  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkNodeHeight_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkNodeHeight_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'StalkNodeHeight_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkNodeHeight_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkNodeHeight_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='living internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkNodeVertLength_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='senescing internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkNodeVertLength_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'StalkNodeVertLength_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StalkNodeVertLength_brch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StalkNodeVertLength_brch', dim1name='pft',dim2name='nodes1',&
     dim3name='nbranches',long_name='senescing internode height', units='m', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then 
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StructInternodeElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='internode element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StructInternodeElms_brch,datrp_4d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'StructInternodeElms_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StructInternodeElms_brch,datrp_4d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StructInternodeElms_brch', dim1name='pft',dim2name='elmnts',&
     dim3name='nodes1',dim4name='nbranches',long_name='internode element', units='g d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then 
    datpr4 => datrp_4d(1:npfts,1:NumCanopyLayers,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyLeafArea_lnode', dim1name='pft',dim2name='levcan',&
     dim3name='nodes1',dim4name='nbranches',long_name='plant layer node branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafArea_lnode,datrp_4d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)

  else
    !print*,'CanopyLeafArea_lnode'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyLeafArea_lnode,datrp_4d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr4 => datrp_4d(1:npfts,1:NumCanopyLayers,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyLeafArea_lnode', dim1name='pft',dim2name='levcan',&
     dim3name='nodes1',dim4name='nbranches',long_name='plant layer node branch leaf area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr5 => datrp_5d(1:npfts,1:NumPlantChemElms,1:NumCanopyLayers,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafLayerElms_node', dim1name='pft',dim2name='elmnts',&
     dim3name='levcan',dim4name='nodes1',dim5name='nbranches',long_name='layer leaf element',&
     units='g d-2',interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafLayerElms_node,datrp_5d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafLayerElms_node'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafLayerElms_node,datrp_5d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr5 => datrp_5d(1:npfts,1:NumPlantChemElms,1:NumCanopyLayers,1:MaxNodesPerBranch+1,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafLayerElms_node', dim1name='pft',dim2name='elmnts',&
      dim3name='levcan',dim4name='nodes1',dim5name='nbranches',long_name='layer leaf element',&
      units='g d-2',interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr5 => datrp_5d(1:npfts,1:NumLeafInclinationClasses,1:NumCanopyLayers,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaZsec_brch', dim1name='pft',dim2name='cansecz',&
      dim3name='levcan',dim4name='nodes',dim5name='nbranches',long_name='leaf surface area',&
      units='m2 d-2',interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafAreaZsec_brch,datrp_5d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'LeafAreaZsec_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,LeafAreaZsec_brch,datrp_5d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr5 => datrp_5d(1:npfts,1:NumLeafInclinationClasses,1:NumCanopyLayers,1:MaxNodesPerBranch,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='LeafAreaZsec_brch', dim1name='pft',dim2name='cansecz',&
      dim3name='levcan',dim4name='nodes',dim5name='nbranches',long_name='leaf surface area',&
      units='m2 d-2',interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)
  endif  

  if(flag=='read')then 
    datpr3 => datrp_3d(1:npfts,1:NumCanopyLayers,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyStalkSurfArea_lbrch', dim1name='pft',dim2name='levcan',&
     dim3name='nbranches',long_name='plant canopy layer branch stem layer area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyStalkSurfArea_lbrch,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'CanopyStalkSurfArea_lbrch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CanopyStalkSurfArea_lbrch,datrp_3d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:NumCanopyLayers,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='CanopyStalkSurfArea_lbrch', dim1name='pft',dim2name='levcan',&
     dim3name='nbranches',long_name='plant canopy layer branch stem layer area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then 
    datpr4 => datrp_4d(1:npfts,1:NumLeafInclinationClasses,1:NumCanopyLayers,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StemAreaZsec_brch', dim1name='pft',dim2name='cansecz',&
      dim3name='levcan',dim4name='nbranches',long_name='stem surface area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)       
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StemAreaZsec_brch,datrp_4d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'StemAreaZsec_brch'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,StemAreaZsec_brch,datrp_4d,&
      NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
    datpr4 => datrp_4d(1:npfts,1:NumLeafInclinationClasses,1:NumCanopyLayers,1:MaxNumBranches)
    call restartvar(ncid, flag, varname='StemAreaZsec_brch', dim1name='pft',dim2name='cansecz',&
      dim3name='levcan',dim4name='nbranches',long_name='stem surface area', units='m2 d-2', &
     interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)       
  endif  

  if(plant_model)then
  if(flag=='read')then
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='PSIRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSIRoot_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSIRoot_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
    call restartvar(ncid, flag, varname='PSIRoot_pvr', dim1name='pft',dim2name='rootyps',&
     dim3name='levsoi',long_name='root total water potential', units='MPa', &
     interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
  endif  
  endif
  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='FineRootVolPerMassC_pft', dim1name='pft',dim2name='rootyps',&
     long_name='root volume per C mass', units='m3 g-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FineRootVolPerMassC_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,FineRootVolPerMassC_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='FineRootVolPerMassC_pft', dim1name='pft',dim2name='rootyps',&
     long_name='root volume per C mass', units='m3 g-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootPoreTortu4Gas_pft', dim1name='pft',dim2name='rootyps',&
     long_name='root pore tortuosity for gas transp', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootPoreTortu4Gas_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootPoreTortu4Gas_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootPoreTortu4Gas_pft', dim1name='pft',dim2name='rootyps',&
     long_name='root pore tortuosity for gas transp', units='none', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndXSecArea_pft', dim1name='pft',dim2name='rootyps',&
     long_name='fine root cross section area', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndXSecArea_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndXSecArea_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndXSecArea_pft', dim1name='pft',dim2name='rootyps',&
     long_name='fine root cross section area', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root1stXSecArea_pft', dim1name='pft',dim2name='rootyps',&
     long_name='primary root cross section area', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stXSecArea_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stXSecArea_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root1stXSecArea_pft', dim1name='pft',dim2name='rootyps',&
     long_name='primary root cross section area', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  
  
  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root1stMaxRadius1_pft', dim1name='pft',dim2name='rootyps',&
     long_name='primary root radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stMaxRadius1_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stMaxRadius1_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root1stMaxRadius1_pft', dim1name='pft',dim2name='rootyps',&
     long_name='primary root radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndMaxRadius1_pft', dim1name='pft',dim2name='rootyps',&
     long_name='fine root radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndMaxRadius1_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndMaxRadius1_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndMaxRadius1_pft', dim1name='pft',dim2name='rootyps',&
     long_name='fine root radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootRaidus_rpft', dim1name='pft',dim2name='rootyps',&
     long_name='root internal radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootRaidus_rpft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootRaidus_rpft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='RootRaidus_rpft', dim1name='pft',dim2name='rootyps',&
     long_name='root internal radius', units=' ', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndSpecLen_pft', dim1name='pft',dim2name='rootyps',&
     long_name='specific root length secondary axes', units='m g-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndSpecLen_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
  else
    !print*,'PSIRoot'
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndSpecLen_pft,datrp_2d,NumActivePlants=restart_num_plants,&
      IsPlantActive_pft=restart_active_pft)
    datpr2 => datrp_2d(1:npfts,1:pltpar%jroots)
    call restartvar(ncid, flag, varname='Root2ndSpecLen_pft', dim1name='pft',dim2name='rootyps',&
     long_name='specific root length secondary axes', units='m g-1', &
     interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(plant_model)then
    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PSIRootOSMO_vr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root osmotic water potential', units='MPa', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSIRootOSMO_vr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'PSIRootOSMO_vr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSIRootOSMO_vr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PSIRootOSMO_vr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root osmotic water potential', units='MPa', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  
  
    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PSIRootTurg_vr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root turgor water potential', units='MPa', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSIRootTurg_vr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'PSIRootTurg_vr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,PSIRootTurg_vr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PSIRootTurg_vr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root turgor water potential', units='MPa', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr2 => datrp_2d(1:npfts,1:JZ)
      call restartvar(ncid, flag, varname='RootMediumXNum_pvr', dim1name='pft',&
      dim2name='levsoi',long_name='soil layer root/myco number medium axes', units='# d-2', &
      interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumXNum_pvr,datrp_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else

      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumXNum_pvr,datrp_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr2 => datrp_2d(1:npfts,1:JZ)
      call restartvar(ncid, flag, varname='RootMediumXNum_pvr', dim1name='pft',&
      dim2name='levsoi',long_name='soil layer root/myco number medium axes', units='# d-2', &
      interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr2 => datrp_2d(1:npfts,1:JZ)
      call restartvar(ncid, flag, varname='Root1stXNumL_pvr', dim1name='pft',&
      dim2name='levsoi',long_name='soil layer root/myco number primary axes', units='# d-2', &
      interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stXNumL_pvr,datrp_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'Root1stXNumL_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stXNumL_pvr,datrp_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr2 => datrp_2d(1:npfts,1:JZ)
      call restartvar(ncid, flag, varname='Root1stXNumL_pvr', dim1name='pft',&
      dim2name='levsoi',long_name='soil layer root/myco number primary axes', units='# d-2', &
      interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndXNumL_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer 2nd order root number axes', units='# d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndXNumL_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'Root2ndXNumL_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndXNumL_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndXNumL_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer 2nd order root number axes', units='# d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)        
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootTotLenPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer length per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootTotLenPerPlant_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootTotLenPerPlant_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootTotLenPerPlant_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootTotLenPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer length per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootAbsorbLenPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer absorptive root length per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootAbsorbLenPerPlant_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootAbsorbLenPerPlant_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootAbsorbLenPerPlant_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootAbsorbLenPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer absorptive root length per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootLenDensPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer root length density', units='m m-3', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootLenDensPerPlant_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootLenDensPerPlant_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootLenDensPerPlant_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootLenDensPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer root length density', units='m m-3', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootPoreVol_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer volume air', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootPoreVol_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootPoreVol_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootPoreVol_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootPoreVol_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer volume air', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootVH2O_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer volume water', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootVH2O_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RTVLW'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootVH2O_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootVH2O_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer volume water', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMedTransptArea_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer medium size root cross section area for gas/water transport', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMedTransptArea_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMedTransptArea_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMedTransptArea_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer medium size cross section area for gas/water transport', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  


    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root1stTransptArea_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer primary root cross section area for gas/water transport', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stTransptArea_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stTransptArea_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root1stTransptArea_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer primary root cross section area for gas/water transport', units='m2 d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root1stRadius_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer primary root radius', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stRadius_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'Root1stRadius_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stRadius_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root1stRadius_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer primary root radius', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndRadius_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer secondary root radius', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndRadius_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'Root2ndRadius_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndRadius_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndRadius_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='soil layer secondary root radius', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootSAreaPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer area per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootSAreaPerPlant_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootSAreaPerPlant_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootSAreaPerPlant_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootSAreaPerPlant_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer area per plant', units='m p-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 

    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndEffLen4uptk_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='Mean 2nd root length', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndEffLen4uptk_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'Root2ndEffLen4uptk_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndEffLen4uptk_rpvr,datrp_3d,&
        NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='Root2ndEffLen4uptk_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='Mean 2nd root length', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootRespPotent_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root respiration unconstrained by O2', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootRespPotent_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootRespPotent_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootRespPotent_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootRespPotent_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root respiration unconstrained by O2', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootCO2Autor_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root respiration constrained by O2', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootCO2Autor_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootCO2Autor_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootCO2Autor_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootCO2Autor_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root respiration constrained by O2', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootCO2EmisPot_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root CO2 efflux unconstrained by root nonstructural C', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootCO2EmisPot_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootCO2EmisPot_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootCO2EmisPot_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootCO2EmisPot_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root CO2 efflux unconstrained by root nonstructural C', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then 
      datpr4 => datrp_4d(1:npfts,1:trc_confs%NGasTracers,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='trcg_rootml_pvr', dim1name='pft',dim2name='gastrcs',&
        dim3name='rootyps',dim4name='levsoi',long_name='root gaseous tracer contenta', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)       
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,trcg_rootml_pvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'trcg_rootml_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,trcg_rootml_pvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts,1:trc_confs%NGasTracers,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='trcg_rootml_pvr', dim1name='pft',dim2name='gastrcs',&
        dim3name='rootyps',dim4name='levsoi',long_name='root gaseous tracer contenta', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)     
    endif  

    if(flag=='read')then 
      datpr4 => datrp_4d(1:npfts,1:trc_confs%NGasTracers,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='trcs_rootml_pvr', dim1name='pft',dim2name='gastrcs',&
        dim3name='rootyps',dim4name='levsoi',long_name='root gaseous tracer contenta', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)     
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,trcs_rootml_pvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'trcs_rootml_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,trcs_rootml_pvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts,1:trc_confs%NGasTracers,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='trcs_rootml_pvr', dim1name='pft',dim2name='gastrcs',&
        dim3name='rootyps',dim4name='levsoi',long_name='root gaseous tracer contenta', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)     
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMycoActiveBiomC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer structural C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMycoActiveBiomC_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootMycoActiveBiomC_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMycoActiveBiomC_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMycoActiveBiomC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer structural C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PopuRootMycoC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root/mycorh layer C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col, PopuRootMycoC_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,' PopuRootMycoC_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col, PopuRootMycoC_pvr,datrp_3d,&
        NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='PopuRootMycoC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root/mycorh layer C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootProteinC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer protein C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootProteinC_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootProteinC_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootProteinC_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootProteinC_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root layer protein C', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif 

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootO2Dmnd4Resp_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root O2 demand from respiration', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootO2Dmnd4Resp_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootO2Dmnd4Resp_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootO2Dmnd4Resp_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootO2Dmnd4Resp_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root O2 demand from respiration', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNH4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NH4 non-band unconstrained by NH4', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNH4DmndSoil_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootNH4DmndSoil_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNH4DmndSoil_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNH4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NH4 non-band unconstrained by NH4', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNO3DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NH4 band unconstrained by NH4', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNO3DmndSoil_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootNO3DmndSoil_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNO3DmndSoil_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNO3DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NH4 band unconstrained by NH4', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH2PO4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of H2PO4 non-band', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootH2PO4DmndSoil_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootH2PO4DmndSoil_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootH2PO4DmndSoil_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH2PO4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of H2PO4 non-band', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH1PO4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='HPO4 demand in non-band by each root population', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootH1PO4DmndSoil_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootH1PO4DmndSoil_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootH1PO4DmndSoil_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH1PO4DmndSoil_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='HPO4 demand in non-band by each root population', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNH4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NO3 band unconstrained by NO3', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNH4DmndBand_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootNH4DmndBand_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNH4DmndBand_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNH4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NO3 band unconstrained by NO3', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNO3DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NO3 non-band unconstrained by NO3', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNO3DmndBand_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootNO3DmndBand_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNO3DmndBand_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNO3DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of NO3 non-band unconstrained by NO3', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH2PO4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of H2PO4 band', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootH2PO4DmndBand_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootH2PO4DmndBand_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootH2PO4DmndBand_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH2PO4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root uptake of H2PO4 band', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH1PO4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='HPO4 demand in band by each root population', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootH1PO4DmndBand_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootH1PO4DmndBand_pvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootH1PO4DmndBand_pvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootH1PO4DmndBand_pvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='HPO4 demand in band by each root population', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RAutoRootO2Limter_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='O2 constraint/stress to root respiration', units='none', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RAutoRootO2Limter_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'WFR'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RAutoRootO2Limter_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RAutoRootO2Limter_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='O2 constraint/stress to root respiration', units='none', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMycoNonstElms_rpvr', dim1name='pft',dim2name='elmnts', &
      dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col, RootMycoNonstElms_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,' RootMycoNonstElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col, RootMycoNonstElms_rpvr,datrp_4d,&
        NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootMycoNonstElms_rpvr', dim1name='pft',dim2name='elmnts', &
      dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootProteinConc_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root protein concentration', units='none', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootProteinConc_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'PSIRoot'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootProteinConc_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootProteinConc_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',long_name='root protein concentration', units='none', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)   
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
      root_cache_name='RootNonstructElmConc_pvr'
      call ncd_inqvid(ncid,trim(root_cache_name),root_conc_varid,root_conc_vardesc,readvar=has_root_conc)
      ! Accept checkpoints written before the component-array rename.
      if(.not.has_root_conc)root_cache_name='RootNonstructElmConc_rpvr'
      call restartvar(ncid, flag, varname=trim(root_cache_name), dim1name='pft',dim2name='elmnts', &
      dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element concentration', units='none', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col, RootNonstructElmConc_pvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col, RootNonstructElmConc_pvr,datrp_4d,&
        NumActivePlants=restart_num_plants,IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts,1:NumPlantChemElms, 1:pltpar%jroots,1:JZ)
      call restartvar(ncid, flag, varname='RootNonstructElmConc_pvr', dim1name='pft',dim2name='elmnts', &
      dim3name='rootyps',dim4name='levsoi',long_name='root layer nonstructural element concentration', units='none', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr2 => datrp_2d(1:npfts, 1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root1stDepz_raxes', dim1name='pft',&
      dim2name='rootaxs',long_name='primary root layer depth', units='m', &
      interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stDepz_raxes,datrp_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'Root1stDepz_raxes'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stDepz_raxes,datrp_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr2 => datrp_2d(1:npfts,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root1stDepz_raxes', dim1name='pft',&
      dim2name='rootaxs',long_name='primary root layer depth', units='m', &
      interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)     
      call restartvar(ncid, flag, varname='NActiveRootSegs_raxes', dim1name='pft',dim2name='rootaxs',&
        long_name='Number of active root segments for root axes', units='none', &
        interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)       
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NActiveRootSegs_raxes,datip_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,NActiveRootSegs_raxes,datip_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)  
      call restartvar(ncid, flag, varname='NActiveRootSegs_raxes', dim1name='pft',dim2name='rootaxs',&
       long_name='Number of active root segments for root axes', units='none', &
       interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)     
    endif

    if(flag=='read')then
      dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)     
      call restartvar(ncid, flag, varname='IndRootSegBase_raxes', dim1name='pft',dim2name='rootaxs',&
        long_name='Index for base segment under tracking for root axes', units='none', &
        interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)       
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,IndRootSegBase_raxes,datip_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,IndRootSegBase_raxes,datip_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)  
      call restartvar(ncid, flag, varname='IndRootSegBase_raxes', dim1name='pft',dim2name='rootaxs',&
       long_name='Index for base segment under tracking for root axes', units='none', &
       interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)     
    endif

    if(flag=='read')then
      dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)     
      call restartvar(ncid, flag, varname='IndRootSegTip_raxes', dim1name='pft',dim2name='rootaxs',&
        long_name='Index for tip segment under tracking for root axes', units='none', &
        interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)       
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,IndRootSegTip_raxes,datip_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,IndRootSegTip_raxes,datip_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      dat2pr => datip_2d(1:npfts,1:MaxNumRootAxes)  
      call restartvar(ncid, flag, varname='IndRootSegTip_raxes', dim1name='pft',dim2name='rootaxs',&
       long_name='Index for tip segment under tracking for root axes', units='none', &
       interpinic_flag='skip', data=dat2pr, missing_value=ispval, fill_value=ispval)     
    endif

    if(flag=='read')then
      datpr2 => datrp_2d(1:npfts,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootSegBaseDepth_raxes', dim1name='pft',dim2name='rootaxs',&
        long_name='base depth of different root axes', units='m', &
        interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootSegBaseDepth_raxes,datrp_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootSegBaseDepth_raxes,datrp_2d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr2 => datrp_2d(1:npfts,1:MaxNumRootAxes)  
      call restartvar(ncid, flag, varname='RootSegBaseDepth_raxes', dim1name='pft',dim2name='rootaxs',&
        long_name='base depth of different root axes', units='m', &
        interpinic_flag='skip', data=datpr2, missing_value=spval, fill_value=spval)     
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:NMaxRootSegs,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootSeglengths_raxes', dim1name='pft',dim2name='rootsegs',&
      dim3name='rootaxs',long_name='root length in each segment for root axes', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootSeglengths_raxes,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else      
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootSeglengths_raxes,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts, 1:NMaxRootSegs,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootSeglengths_raxes', dim1name='pft',dim2name='rootsegs',&
      dim3name='rootaxs',long_name='root length in each segment for root axes', units='m', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='SapFlowVLinear_rpvr', dim1name='pft',dim2name='levsoi',&
      dim3name='rootaxs',long_name='Linear sap flow for primary root axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SapFlowVLinear_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else      
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,SapFlowVLinear_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='SapFlowVLinear_rpvr', dim1name='pft',dim2name='levsoi',&
      dim3name='rootaxs',long_name='Linear sap flow for primary root axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:NMaxRootSegs,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootSegAges_raxes', dim1name='pft',dim2name='rootsegs',&
      dim3name='rootaxs',long_name='age of different active root segments for root axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootSegAges_raxes,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else      
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootSegAges_raxes,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts, 1:NMaxRootSegs,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootSegAges_raxes', dim1name='pft',dim2name='rootsegs',&
      dim3name='rootaxs',long_name='age of different active root segments for root axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms, 1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco1stElm_raxs', dim1name='pft',dim2name='elmnts',&
      dim3name='rootaxs',long_name='elmnts in primary root axes', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco1stElm_raxs,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootMyco1stElm_raxs'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco1stElm_raxs,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms, 1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco1stElm_raxs', dim1name='pft',dim2name='elmnts',&
      dim3name='rootaxs',long_name='elmnts in primary root axes', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco2ndSinkC_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='rootaxs',long_name='fine root carbon sink', units='gC d-2 h-1', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco2ndSinkC_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RTN2'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco2ndSinkC_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco2ndSinkC_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='rootaxs',long_name='fine root carbon sink', units='gC d-2 h-1', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='CytokininMRConc_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='cytokinin concentration in medium size roots', units='gC m-3 H2O', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CytokininMRConc_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CytokininMRConc_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='CytokininMRConc_rpvr', dim1name='pft',&
        dim2name='levsoi',dim3name='rootaxs',long_name='cytokinin concentration in medium size roots', units='gC m-3 H2O', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Cytokinin1stConc_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='cytokinin concentration in primary roots', units='gC m-3 H2O', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Cytokinin1stConc_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Cytokinin1stConc_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Cytokinin1stConc_rpvr', dim1name='pft',&
        dim2name='levsoi',dim3name='rootaxs',long_name='cytokinin concentration in primary roots', units='gC m-3 H2O', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Cytokinin2ndConc_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='rootaxs',long_name='cytokinin concentration in fine roots/mycorrhizae', units='gC m-3 H2O', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Cytokinin2ndConc_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Cytokinin2ndConc_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Cytokinin2ndConc_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='rootaxs',long_name='cytokinin concentration in fine roots/mycorrhizae', units='gC m-3 H2O', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root2ndXNum_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer secondary axes number', units='d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndXNum_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RTN2'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndXNum_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root2ndXNum_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer secondary axes number', units='d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='CRootLumenArea_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='coarse root axes lumen area', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CRootLumenArea_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,CRootLumenArea_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='CRootLumenArea_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='coarse root axes lumen area', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    endif

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMediumLength_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='root length of intermediate size root axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumLength_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumLength_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMediumLength_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='root length of intermediate size root axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    endif

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMediumRadius_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='root radius of intermediate size root axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumRadius_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumRadius_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMediumRadius_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='root radius of intermediate size root axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    endif

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootAge_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='root age of primary axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootAge_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootAge_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootAge_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='root age of primary axes', units='h', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    endif

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco1stSinkC_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='carbon sink by primary axes', units='m d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco1stSinkC_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco1stSinkC_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco1stSinkC_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='carbon sink by primary axes', units='m d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    endif

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root1stLenPP_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='primary root axis length in soil layer', units='m d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stLenPP_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'Root1stLenPP_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stLenPP_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root1stLenPP_rpvr', dim1name='pft',&
      dim2name='levsoi',dim3name='rootaxs',long_name='primary root axis length in soil layer', units='m d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)    
    endif

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco1stStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer primary axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco1stStrutElms_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootMyco1stStrutElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco1stStrutElms_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco1stStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer primary axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMediumStructElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer medium axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumStructElms_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMediumStructElms_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMediumStructElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer medium axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root1stActStructElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer primary axes active zone element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stActStructElms_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stActStructElms_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root1stActStructElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer primary axes active zone element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms, 1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root1stLigStructElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer primary axes lignified zone element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stLigStructElms_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root1stLigStructElms_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts, 1:NumPlantChemElms,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root1stLigStructElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer primary axes lignified zone element', units='g d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)
    endif  

    if(flag=='read')then
      datpr5 => datrp_5d(1:npfts, 1:NumPlantChemElms, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco2ndStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='rootyps',dim4name='levsoi',dim5name='rootaxs',long_name='root layer secondary axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco2ndStrutElms_rpvr,datrp_5d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootMyco2ndStrutElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootMyco2ndStrutElms_rpvr,datrp_5d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr5 => datrp_5d(1:npfts, 1:NumPlantChemElms, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='RootMyco2ndStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='rootyps',dim4name='levsoi',dim5name='rootaxs',long_name='root layer secondary axes element', units='g d-2', &
      interpinic_flag='skip', data=datpr5, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root2ndLen_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer length secondary axes', units='m d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndLen_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'Root2ndLen_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,Root2ndLen_rpvr,datrp_4d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr4 => datrp_4d(1:npfts, 1:pltpar%jroots,1:JZ,1:MaxNumRootAxes)
      call restartvar(ncid, flag, varname='Root2ndLen_rpvr', dim1name='pft',dim2name='rootyps',&
      dim3name='levsoi',dim4name='rootaxs',long_name='root layer length secondary axes', units='m d-2', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
      call restartvar(ncid, flag, varname='RootNodulNonstElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNodulNonstElms_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootNodulNonstElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNodulNonstElms_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
      call restartvar(ncid, flag, varname='RootNodulNonstElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',long_name='root layer nonstructural element', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  

    if(flag=='read')then
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
      call restartvar(ncid, flag, varname='RootNodulStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',long_name='root layer nodule element', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval)  
      call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNodulStrutElms_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
    else
      !print*,'RootNodulStrutElms_rpvr'
      if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,RootNodulStrutElms_rpvr,datrp_3d,NumActivePlants=restart_num_plants,&
        IsPlantActive_pft=restart_active_pft)
      datpr3 => datrp_3d(1:npfts, 1:NumPlantChemElms,1:JZ)
      call restartvar(ncid, flag, varname='RootNodulStrutElms_rpvr', dim1name='pft',dim2name='elmnts',&
      dim3name='levsoi',long_name='root layer nodule element', units='g d-2', &
      interpinic_flag='skip', data=datpr3, missing_value=spval, fill_value=spval) 
    endif  
  endif
  call destroy(datrp_2d)
  call destroy(datrp_3d)
  call destroy(datrp_4d)  
  call destroy(datrp_5d)  
  call destroy(datip_3d)

  ! Preserve geometry established at planting and cached branch activity.
  call restart_plant_profile(ncid,flag,'Root1stMaxRadius_pft',Root1stMaxRadius_pft,'m')
  call restart_plant_profile(ncid,flag,'Root2ndMaxRadius_pft',Root2ndMaxRadius_pft,'m')
  call restart_plant_profile(ncid,flag,'RubiscoActivity_brch',RubiscoActivity_brch,'none','nbranches')
  call restart_plant_profile(ncid,flag,'Root1stSpecLen_pft',Root1stSpecLen_pft,'m gC-1')
  ! These caches can outlive the current crop and are needed before planting
  ! and by the uptake API. STARTQ cannot reconstruct them in a fallow slot.
  call restart_pft_memory(ncid,flag,'RootSingleVesselArea_pft',RootSingleVesselArea_pft,'m2')
  call restart_pft_memory(ncid,flag,'RootSingleVesselRstaxial_pft',RootSingleVesselRstaxial_pft,'MPa h m-4')
  call restart_pft_memory(ncid,flag,'CRootActVolPerMassC_pft',CRootActVolPerMassC_pft,'m3 gC-1')
  call restart_pft_memory(ncid,flag,'CRootLigVolPerMassC_pft',CRootLigVolPerMassC_pft,'m3 gC-1')
  call restart_pft_memory(ncid,flag,'TempOffset_pft',TempOffset_pft,'degC')
  call restart_pft_memory(ncid,flag,'rPlantThermoAdaptZone_pft',rPlantThermoAdaptZone_pft,'none')
  call restart_pft_memory(ncid,flag,'HighTempLimitSeed_pft',HighTempLimitSeed_pft,'degC')
  ! Prior-step canopy, phenology and uptake caches remain observable at replanting.
  call restart_pft_memory(ncid,flag,'StalkHeight_pft',StalkHeight_pft,'m')
  call restart_pft_memory(ncid,flag,'StalkAveRadius_pft',StalkAveRadius_pft,'m')
  call restart_pft_memory(ncid,flag,'PlantinDepz_pft',PlantinDepz_pft,'m')
  call restart_pft_memory(ncid,flag,'RawIsoTCanopy2Atm_pft',RawIsoTCanopy2Atm_pft,'h m-1')
  call restart_pft_memory(ncid,flag,'CanopyHeight4WatUptake_pft',CanopyHeight4WatUptake_pft,'m')
  call restart_pft_memory(ncid,flag,'CanopySeedNumX_pft',CanopySeedNumX_pft,'d-2')
  call restart_pft_memory(ncid,flag,'CanPhenoMoistStress_pft',CanPhenoMoistStress_pft,'-')
  call restart_pft_memory(ncid,flag,'CanPhenoTempStress_pft',CanPhenoTempStress_pft,'-')
  call restart_pft_memory(ncid,flag,'PSICanPDailyMin_pft',PSICanPDailyMin_pft,'MPa')
  call restart_pft_memory(ncid,flag,'Root1stTipSinkWeight_pft',Root1stTipSinkWeight_pft,'d-2')
  call restart_pft_memory(ncid,flag,'NH3Dep2Can_pft',NH3Dep2Can_pft,'g d-2 h-1')
  call restart_pft_memory(ncid,flag,'PTSHTR_pft',PTSHTR_pft,'h-1')
  call restart_pft_memory(ncid,flag,'RootNH4Uptake_pft',RootNH4Uptake_pft,'g d-2 h-1')
  call restart_pft_memory(ncid,flag,'RootNO3Uptake_pft',RootNO3Uptake_pft,'g d-2 h-1')
  call restart_pft_memory(ncid,flag,'RootH2PO4Uptake_pft',RootH2PO4Uptake_pft,'g d-2 h-1')
  call restart_pft_memory(ncid,flag,'RootHPO4Uptake_pft',RootHPO4Uptake_pft,'g d-2 h-1')
  call restart_pft_memory(ncid,flag,'TurgEff4CanopyResp_pft',TurgEff4CanopyResp_pft,'-')
  call restart_pft_memory(ncid,flag,'CanPStomaResistH2O_pft',CanPStomaResistH2O_pft,'h m-1')
  call restart_pft_memory(ncid,flag,'CanopyMinStomaResistH2O_pft',CanopyMinStomaResistH2O_pft,'h m-1')
  call restart_pft_memory(ncid,flag,'O2I_pft',O2I_pft,'umol mol-1')
  call restart_pft_memory(ncid,flag,'LeafIntracellularCO2_pft',LeafIntracellularCO2_pft,'umol mol-1')
  call restart_pft_memory(ncid,flag,'AirConc_pft',AirConc_pft,'mol m-3')
  call restart_pft_memory(ncid,flag,'DiffCO2Atmos2Intracel_pft',DiffCO2Atmos2Intracel_pft,'umol m-3')
  call restart_pft_memory(ncid,flag,'CanopyGasCO2_pft',CanopyGasCO2_pft,'umol mol-1')
  call restart_pft_memory(ncid,flag,'aquCO2Intraleaf_pft',aquCO2Intraleaf_pft,'uM')
  call restart_pft_memory(ncid,flag,'O2L_pft',O2L_pft,'uM')
  call restart_pft_memory(ncid,flag,'CO2Solubility_pft',CO2Solubility_pft,'uM /umol mol-1')
  call restart_pft_memory(ncid,flag,'LeafO2Solubility_pft',LeafO2Solubility_pft,'uM /umol mol-1')
  call restart_pft_memory(ncid,flag,'Km4LeafaqCO2_pft',Km4LeafaqCO2_pft,'uM')
  call restart_pft_memory(ncid,flag,'Km4RubiscoCarboxy_pft',Km4RubiscoCarboxy_pft,'uM')
  call restart_pft_memory(ncid,flag,'RadNet2Canopy_pft',RadNet2Canopy_pft,'MJ d-2 h-1')
  call restart_pft_memory(ncid,flag,'LWRadCanopy_pft',LWRadCanopy_pft,'MJ d-2 h-1')
  call restart_pft_memory(ncid,flag,'CanopyEvapTransLHeat_pft',CanopyEvapTransLHeat_pft,'MJ d-2 h-1')
  call restart_pft_memory(ncid,flag,'HeatXAir2PCan_pft',HeatXAir2PCan_pft,'MJ d-2 h-1')
  call restart_pft_memory(ncid,flag,'HeatStorCanopy_pft',HeatStorCanopy_pft,'MJ d-2 h-1')
  ! Read the former water-resistance name from checkpoints predating the split.
  has_water_resistance=.true.
  has_heat_resistance=.true.
  if(flag=='read')then
    call ncd_inqvid(ncid,'RawH2OCanopy2Atm_pft',root_conc_varid,root_conc_vardesc,readvar=has_water_resistance)
    call ncd_inqvid(ncid,'RawHeatCanopy2Atm_pft',root_conc_varid,root_conc_vardesc,readvar=has_heat_resistance)
  endif
  if(has_water_resistance)then
    call restart_pft_memory(ncid,flag,'RawH2OCanopy2Atm_pft',RawH2OCanopy2Atm_pft,'h m-1')
  else
    call restart_pft_memory(ncid,flag,'RawCanopy2Atm_pft',RawH2OCanopy2Atm_pft,'h m-1')
  endif
  if(has_heat_resistance)then
    call restart_pft_memory(ncid,flag,'RawHeatCanopy2Atm_pft',RawHeatCanopy2Atm_pft,'h m-1')
  else
    ! Heat resistance is derived; preserve inactive/missing entries unchanged.
    RawHeatCanopy2Atm_pft=RawH2OCanopy2Atm_pft
    where(RawH2OCanopy2Atm_pft/=spval) &
      RawHeatCanopy2Atm_pft=RawH2OCanopy2Atm_pft/BndlResistH2OToHeatRatio
  endif
  call restart_pft_memory(ncid,flag,'Transpiration_pft',Transpiration_pft,'m3 d-2 h-1')
  call restart_pft_memory(ncid,flag,'VapXAir2Canopy_pft',VapXAir2Canopy_pft,'m3 d-2 h-1')
  call restart_pft_memory(ncid,flag,'CO2NetFix_pft',CO2NetFix_pft,'gC d-2 h-1')
  call restart_pft_memory(ncid,flag,'CanopyNoduleNonstCConc_pft',CanopyNoduleNonstCConc_pft,'gC gC-1')
  call restart_pft_memory(ncid,flag,'TFN_Carboxy_pft',TFN_Carboxy_pft,'-')
  call restart_pft_memory(ncid,flag,'TFN_Oxygen_pft',TFN_Oxygen_pft,'-')
  call restart_pft_memory(ncid,flag,'TFN_eTranspt_pft',TFN_eTranspt_pft,'-')
  call restart_pft_memory(ncid,flag,'CdH2ORootxSoil_pft',CdH2ORootxSoil_pft,'m3 h-1 MPa-1')
  call restart_plant_profile(ncid,flag,'Root2ndAxialResist_pft',Root2ndAxialResist_pft,'MPa h m-4')
  call restart_plant_profile(ncid,flag,'RootSinkWeight_pvr',RootSinkWeight_pvr,'d-2','levsoi')
  call restart_root_layer_memory(ncid,flag,'RootRadialKond2H2O_pvr',RootRadialKond2H2O_pvr,'m3 h-1 MPa-1')
  call restart_root_layer_memory(ncid,flag,'RootAxialKond2H2O_pvr',RootAxialKond2H2O_pvr,'m3 h-1 MPa-1')
  DO ig=idg_beg,idg_NH3
    root_cache_name='RootAtmGasConductance_'//trim(int2str(ig))
    call restart_root_layer_memory(ncid,flag,trim(root_cache_name),RootAtmGasConductance_rpvr(ig,:,:,:,:,:), &
      'm3 h-1')
  ENDDO
  if(flag=='read' .and. .not.retained_inactive_plant)then
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        DO NZ=1,NP0_col(NY,NX)
          if(iPlantStateLive_pft(NZ,NY,NX)==iTrue)cycle
          ! Legacy inactive slots did not save population. Discard the cold
          ! initialization population so the upcoming planting record is read.
          PlantPopuLive_pft(NZ,NY,NX)=0._r8
          loc=DATAPI(NZ,NY,NX)
          if(loc<=0)cycle
          ! Recover the missing hydraulic cache from the residue PFT, before
          ! the annual reader replaces its trait index with the next crop.
          RootSingleVesselArea_pft(NZ,NY,NX)=PICON*RootVesselRadius_tab(loc)**2
          RootSingleVesselRstaxial_pft(NZ,NY,NX)= &
            AlphaVesselAxialResist_tab(loc)*Rax_ref*(1.e-6_r8/RootVesselRadius_tab(loc))**4
        ENDDO
      ENDDO
    ENDDO
  endif
  if(flag=='read')then
    ! Recover the aggregate from its mode-specific source, without a second
    ! checkpoint copy. Daily SP setup subsequently refreshes its root profile.
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        DO NZ=1,NP0_col(NY,NX)
          if(ldo_sp_mode)then
            RootNonstOMConc_pvr(:,:,NZ,NY,NX)=PrescribedRootNonstOMConc_pft(NZ,NY,NX)
            RootMyco2ndStrutElms_pvr(:,:,:,NZ,NY,NX)=0._r8
            if(NZ<=NP_col(NY,NX))then
              ! SP PopuRootMycoC stores the prescribed stock per square meter.
              RootMyco2ndStrutElms_pvr(ielmc,ipltroot,:,NZ,NY,NX)= &
                PopuRootMycoC_pvr(ipltroot,:,NZ,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
            endif
          else
            RootNonstOMConc_pvr(:,:,NZ,NY,NX)=sum(RootNonstructElmConc_pvr(:,:,:,NZ,NY,NX),dim=1)
            RootMyco2ndStrutElms_pvr(:,:,:,NZ,NY,NX)=sum( &
              RootMyco2ndStrutElms_rpvr(:,:,:,1:NumStructuralRootAxes_pft(NZ,NY,NX),NZ,NY,NX),dim=4)
          endif
        ENDDO
      ENDDO
    ENDDO
  endif
  END SUBROUTINE restartnc_plant
!------------------------------------------------------------------------------------------
  subroutine restart_plant_profile(ncid,flag,name,values,units,dimension_name)
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  character(len=*), intent(in) :: flag,name,units
  real(r8), intent(inout) :: values(:,:,:,:)
  character(len=*), optional, intent(in) :: dimension_name
  character(len=16) :: second_dimension
  character(len=*), parameter :: subname='restart_plant_profile'
  real(r8), target :: buffer(bounds%begp:bounds%endp,size(values,1))
  real(r8), pointer :: data(:,:)
  type(var_desc_t) :: vardesc
  integer :: varid
  logical :: present_cache

  associate(NHW => bounds%NHW, NHE => bounds%NHE, NVN => bounds%NVN, NVS => bounds%NVS)
  call PrintInfo('beg '//subname)
  second_dimension='rootyps'
  if(present(dimension_name))second_dimension=dimension_name
  present_cache=.true.
  if(flag=='read')call ncd_inqvid(ncid,name,varid,vardesc,readvar=present_cache)
  if(present_cache)then
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,values,buffer)
    data=>buffer
    call restartvar(ncid,flag,varname=name,dim1name='pft',dim2name=trim(second_dimension), &
      long_name='persistent plant property: '//name,units=units,interpinic_flag='skip', &
      data=data,missing_value=spval,fill_value=spval)
    if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,values,buffer)
  else
    write(iulog,*) 'Restart lacks ',name,'; exact continuation is not guaranteed'
  endif
  call PrintInfo('end '//subname)
  end associate
  end subroutine restart_plant_profile


!------------------------------------------------------------------------------------------

  subroutine restart_root_layer_memory(ncid,flag,name,values,units)
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  character(len=*), intent(in) :: flag,name,units
  real(r8), intent(inout) :: values(:,:,:,:,:)
  character(len=*), parameter :: subname='restart_root_layer_memory'
  real(r8), target :: buffer(bounds%begp:bounds%endp,size(values,1),size(values,2))
  real(r8), pointer :: data(:,:,:)
  type(var_desc_t) :: vardesc
  integer :: varid
  logical :: present_cache

  associate(NHW => bounds%NHW, NHE => bounds%NHE, NVN => bounds%NVN, NVS => bounds%NVS)
  call PrintInfo('beg '//subname)
  present_cache=.true.
  if(flag=='read')call ncd_inqvid(ncid,name,varid,vardesc,readvar=present_cache)
  if(present_cache)then
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,values,buffer)
    data=>buffer
    call restartvar(ncid,flag,varname=name,dim1name='pft',dim2name='rootyps',dim3name='levsoi', &
      long_name='persistent root layer memory: '//name,units=units,interpinic_flag='skip', &
      data=data,missing_value=spval,fill_value=spval)
    if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,values,buffer)
  else
    write(iulog,*) 'Restart lacks ',name,'; exact continuation is not guaranteed'
  endif
  call PrintInfo('end '//subname)
  end associate
  end subroutine restart_root_layer_memory
!------------------------------------------------------------------------------------------
  subroutine restartnc(ncid,flag)
  implicit none
  character(len=*) , intent(in) :: flag
  type(file_desc_t), intent(inout) :: ncid ! netcdf id

  call restartnc_plant(ncid,flag)

  call restartnc_soil(ncid,flag)
  call restartnc_soil_caches(ncid,flag)
  call restartnc_accumulators(ncid,flag)
  call restartnc_daily_irrigation(ncid,flag)

  end subroutine restartnc
!------------------------------------------------------------------------------------------
  subroutine restartnc_soil_hydraulic_scheme(ncid,flag)
  ! Description:
  ! Write scheme identity/revision and separate drainage policy as global attributes.
  ! Validate them before restoring caches. Attribute-free legacy files are EcoSIM;
  ! cross-scheme or cross-policy cache migration is not performed implicitly.
  use netcdf, only: nf90_global,nf90_noerr,nf90_enotatt,nf90_char,nf90_int, &
    nf90_inquire_attribute,nf90_get_att
  implicit none
  type(file_desc_t), intent(inout) :: ncid ! opened checkpoint descriptor [-]
  character(len=*), intent(in) :: flag    ! metadata operation: define or read [-]
  character(len=32) :: selected_name,saved_name
  integer :: selected_revision,saved_revision,saved_policy,status,attribute_type,attribute_length
  logical :: has_name,has_revision,has_policy
  character(len=*), parameter :: subname='restartnc_soil_hydraulic_scheme'

  call PrintInfo('beg '//subname)
  selected_name=GetSoilHydraulicSchemeName()
  selected_revision=GetSoilHydraulicSchemeRevision()
  associate(handle => ncid%fh)
    if(flag=='define')then
      call check_ret(ncd_putatt(ncid,nf90_global,'soil_hydraulic_scheme',trim(selected_name)),subname)
      call check_ret(ncd_putatt(ncid,nf90_global,'soil_hydraulic_scheme_revision',selected_revision),subname)
      call check_ret(ncd_putatt(ncid,nf90_global,'soil_hydraulic_ecosys_drainage',merge(1,0,use_ecosys_drainage)),subname)
    elseif(flag=='read')then
      status=nf90_inquire_attribute(handle,nf90_global,'soil_hydraulic_scheme')
      has_name=status==nf90_noerr
      if(status/=nf90_enotatt)call check_ret(status,subname)
      status=nf90_inquire_attribute(handle,nf90_global,'soil_hydraulic_scheme_revision')
      has_revision=status==nf90_noerr
      if(status/=nf90_enotatt)call check_ret(status,subname)
      status=nf90_inquire_attribute(handle,nf90_global,'soil_hydraulic_ecosys_drainage')
      has_policy=status==nf90_noerr
      if(status/=nf90_enotatt)call check_ret(status,subname)
      if(.not.has_name.and..not.has_revision.and..not.has_policy)then
        if(selected_name/=DefaultSoilHydraulicScheme)then
          write(iulog,*) 'Legacy checkpoint has no scheme identity; requested ',trim(selected_name)
          call endrun('Soil hydraulic restart scheme mismatch')
        endif
        write(iulog,*) 'Restart lacks soil hydraulic scheme metadata; assuming ecosim and retaining configured drainage policy.'
      elseif(.not.(has_name.and.has_revision.and.has_policy))then
        write(iulog,*) 'Checkpoint must contain all three soil hydraulic scheme metadata attributes.'
        call endrun('Incomplete soil hydraulic restart metadata')
      else
        call check_ret(nf90_inquire_attribute(handle,nf90_global,'soil_hydraulic_scheme', &
          xtype=attribute_type,len=attribute_length),subname)
        if(attribute_type/=nf90_char.or.attribute_length<1.or.attribute_length>len(saved_name))then
          call endrun('Invalid soil hydraulic restart scheme attribute')
        endif
        saved_name=''
        call check_ret(nf90_get_att(handle,nf90_global,'soil_hydraulic_scheme',saved_name),subname)
        call check_ret(nf90_inquire_attribute(handle,nf90_global,'soil_hydraulic_scheme_revision', &
          xtype=attribute_type,len=attribute_length),subname)
        if(attribute_type/=nf90_int.or.attribute_length/=1)then
          call endrun('Invalid soil hydraulic restart revision attribute')
        endif
        call check_ret(nf90_get_att(handle,nf90_global,'soil_hydraulic_scheme_revision',saved_revision),subname)
        call check_ret(nf90_inquire_attribute(handle,nf90_global,'soil_hydraulic_ecosys_drainage', &
          xtype=attribute_type,len=attribute_length),subname)
        if(attribute_type/=nf90_int.or.attribute_length/=1)then
          call endrun('Invalid soil hydraulic restart drainage-policy attribute')
        endif
        call check_ret(nf90_get_att(handle,nf90_global,'soil_hydraulic_ecosys_drainage',saved_policy),subname)
        if(saved_name/=selected_name)then
          write(iulog,*) 'Checkpoint scheme=',trim(saved_name),' requested=',trim(selected_name)
          call endrun('Soil hydraulic restart scheme mismatch')
        endif
        if(saved_revision/=selected_revision)then
          write(iulog,*) 'Checkpoint scheme revision=',saved_revision,' requested=',selected_revision
          call endrun('Soil hydraulic restart revision mismatch')
        endif
        if(saved_policy/=merge(1,0,use_ecosys_drainage))then
          write(iulog,*) 'Checkpoint drainage policy=',saved_policy,' requested=',merge(1,0,use_ecosys_drainage)
          call endrun('Soil hydraulic restart drainage-policy mismatch')
        endif
      endif
    endif
  end associate
  call PrintInfo('end '//subname)
  end subroutine restartnc_soil_hydraulic_scheme
!------------------------------------------------------------------------------------------
  subroutine restartnc_accumulators(ncid,flag)
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  character(len=*), intent(in) :: flag
  character(len=*), parameter :: subname='restartnc_accumulators'
  logical :: complete

  call PrintInfo('beg '//subname)
  call restart_column_memory(ncid,flag,'SnowOnCanopy_col',SnowOnCanopy_col,'m3')
  call restart_column_memory(ncid,flag,'fSnowCanopy_col',fSnowCanopy_col,'none')
  call restart_column_memory(ncid,flag,'BulkFactor4Snow_col',BulkFactor4Snow_col,'m2 ton-1')
  call restart_column_memory(ncid,flag,'NWetDep_col',NWetDep_col,'gN')
  call restart_column_memory(ncid,flag,'CumDryDepoOM_col',CumDryDepoOM_col,'gC')
  call restart_column_memory(ncid,flag,'GDD_col',GDD_col,'degC day')
  call restart_column_memory(ncid,flag,'TRAD_col',TRAD_col,'MJ m-2 h-1')
  call restart_column_memory(ncid,flag,'EvapoTransp_col',EvapoTransp_col,'m3')
  call restart_column_memory(ncid,flag,'QDrain_cum_col',QDrain_cum_col,'m3')
  ! Radiation and snow unloading consume these before the current plant step.
  call restart_pft_memory(ncid,flag,'fSnowCanopy_pft',fSnowCanopy_pft,'none')
  call restart_pft_memory(ncid,flag,'BulkFactor4Snow_pft',BulkFactor4Snow_pft,'m2 ton-1')
  call restart_pft_memory(ncid,flag,'LeafAreaSunlit_pft',LeafAreaSunlit_pft,'m2')
  call restart_pft_memory(ncid,flag,'StomatalStress_pft',StomatalStress_pft,'none')
  call restart_pft_memory(ncid,flag,'TKCanopy_pft',TKCanopy_pft,'K')
  call restart_pft_memory(ncid,flag,'GrossRespC_CumYr_pft',GrossRespC_CumYr_pft,'gC')
  call restart_pft_memory(ncid,flag,'RootUptk_Nmin_cumYr_pft',RootUptk_Nmin_cumYr_pft,'gN')
  call restart_pft_memory(ncid,flag,'RootUptk_Pmin_cumYr_pft',RootUptk_Pmin_cumYr_pft,'gP')
  complete=.true.
  call restart_soil_profile_cache(ncid,flag,'trcg_mass_cumerr_col',trcg_mass_cumerr_col, &
    'cumulative gas tracer balance error','g tracer',complete,dimension_name='gastrcs')
  if(flag=='read' .and. .not.complete) &
    write(iulog,*) 'Restart lacks accumulated gas balance errors; history continuation is not exact'
  call PrintInfo('end '//subname)
  end subroutine restartnc_accumulators
!------------------------------------------------------------------------------------------
  subroutine restart_column_memory(ncid,flag,name,values,units)
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  character(len=*), intent(in) :: flag,name,units
  real(r8), intent(inout) :: values(:,:)
  character(len=*), parameter :: subname='restart_column_memory'
  real(r8), target :: buffer(bounds%begc:bounds%endc)
  real(r8), pointer :: data(:)
  type(var_desc_t) :: vardesc
  integer :: varid
  logical :: present_value

  associate(NHW => bounds%NHW, NHE => bounds%NHE, NVN => bounds%NVN, NVS => bounds%NVS)
  call PrintInfo('beg '//subname)
  present_value=.true.
  if(flag=='read')call ncd_inqvid(ncid,name,varid,vardesc,readvar=present_value)
  if(present_value)then
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,values,buffer)
    data=>buffer
    call restartvar(ncid,flag,varname=name,dim1name='column', &
      long_name='persistent column memory: '//name,units=units,interpinic_flag='skip', &
      data=data,missing_value=spval,fill_value=spval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,values,buffer)
  else
    write(iulog,*) 'Restart lacks ',name,'; retaining initialization; exact continuation is not guaranteed'
  endif
  call PrintInfo('end '//subname)
  end associate
  end subroutine restart_column_memory
!------------------------------------------------------------------------------------------
  subroutine restart_pft_memory(ncid,flag,name,values,units)
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  character(len=*), intent(in) :: flag,name,units
  real(r8), intent(inout) :: values(:,:,:)
  character(len=*), parameter :: subname='restart_pft_memory'
  real(r8), target :: buffer(bounds%begp:bounds%endp)
  real(r8), pointer :: data(:)
  type(var_desc_t) :: vardesc
  integer :: varid
  logical :: present_value

  associate(NHW => bounds%NHW, NHE => bounds%NHE, NVN => bounds%NVN, NVS => bounds%NVS)
  call PrintInfo('beg '//subname)
  present_value=.true.
  if(flag=='read')call ncd_inqvid(ncid,name,varid,vardesc,readvar=present_value)
  if(present_value)then
    if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,values,buffer)
    data=>buffer
    call restartvar(ncid,flag,varname=name,dim1name='pft', &
      long_name='persistent plant memory: '//name,units=units,interpinic_flag='skip', &
      data=data,missing_value=spval,fill_value=spval)
    if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP0_col,values,buffer)
  else
    write(iulog,*) 'Restart lacks ',name,'; retaining initialization; exact continuation is not guaranteed'
  endif
  call PrintInfo('end '//subname)
  end associate
  end subroutine restart_pft_memory
!------------------------------------------------------------------------------------------
  subroutine restartnc_soil_caches(ncid,flag)
  ! These coefficients are updated only after disturbances. Re-running that
  ! update changes geometry; leaving them at allocation defaults is unsafe.
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  character(len=*), intent(in) :: flag
  character(len=*), parameter :: subname='restartnc_soil_caches'
  logical :: complete,present_cache
  integer :: varid,L
  type(var_desc_t) :: vardesc
  real(r8), target :: conductivity(bounds%begc:bounds%endc,3,100,JZ+1)
  real(r8), target :: root_mass(bounds%begc:bounds%endc,NumPlantChemElms,jroots,JZ)
  real(r8), target :: irrigation_nutrients(bounds%begc:bounds%endc,size(trcn_irrig_vr,1),JZ)
  real(r8), pointer :: irrigation_data(:,:,:)
  real(r8), pointer :: conductivity_data(:,:,:,:)
  integer, target :: pore_numbers(bounds%begc:bounds%endc,JZ)
  integer, pointer :: pore_number_data(:,:)

  associate(NHW => bounds%NHW, NHE => bounds%NHE, NVN => bounds%NVN, NVS => bounds%NVS)
  call PrintInfo('beg '//subname)
  complete=.true.
  if(allocated(VGThetaR_vr))then
    call restart_soil_profile_cache(ncid,flag,'VG_THETA_R',VGThetaR_vr, &
      'van Genuchten residual liquid content per micropore soil volume','m3 m-3',complete)
    call restart_soil_profile_cache(ncid,flag,'VG_ALPHA',VGAlpha_vr, &
      'van Genuchten inverse pressure scale','MPa-1',complete)
    call restart_soil_profile_cache(ncid,flag,'VG_N',VGN_vr,'van Genuchten shape exponent','none',complete)
    call restart_soil_profile_cache(ncid,flag,'VG_L',VGL_vr,'Mualem pore-connectivity exponent','none',complete)
    if(.not.complete)call endrun('Van Genuchten checkpoint is missing parameter fields')
  endif
  ! Pending litter-property updates use the previous geometry before HOUR1
  ! rebuilds it. Other columns retain their geometry until a disturbance.
  call restart_soil_profile_cache(ncid,flag,'VGeomLayer_vr',VGeomLayer_vr, &
    'soil and surface litter geometric volume','m3',complete)
  call restart_soil_profile_cache(ncid,flag,'SoilDepthMidLay_vr',SoilDepthMidLay_vr, &
    'depth to soil layer midpoint','m',complete)
  call restart_soil_profile_cache(ncid,flag,'CumSoilThickMidL_vr',CumSoilThickMidL_vr, &
    'soil layer midpoint depth from grid cell surface','m',complete)
  call restart_soil_profile_cache(ncid,flag,'LOGPOROS_vr',LOGPOROS_vr, &
    'log soil porosity','none',complete)
  call restart_soil_profile_cache(ncid,flag,'LOGFldCapacity_vr',LOGFldCapacity_vr, &
    'log water content at field capacity','none',complete)
  call restart_soil_profile_cache(ncid,flag,'LOGWiltPoint_vr',LOGWiltPoint_vr, &
    'log water content at wilting point','none',complete)
  call restart_soil_profile_cache(ncid,flag,'PSD_vr',PSD_vr, &
    'log porosity to field-capacity ratio','none',complete)
  call restart_soil_profile_cache(ncid,flag,'FCD_vr',FCD_vr, &
    'log field-capacity to wilting-point ratio','none',complete)
  call restart_soil_profile_cache(ncid,flag,'SRP_vr',SRP_vr, &
    'water retention curve shape','none',complete)
  call restart_soil_profile_cache(ncid,flag,'CSAND_vr',CSAND_vr, &
    'soil sand mass fraction','none',complete)
  call restart_soil_profile_cache(ncid,flag,'CSILT_vr',CSILT_vr, &
    'soil silt mass fraction','none',complete)
  call restart_soil_profile_cache(ncid,flag,'CCLAY_vr',CCLAY_vr, &
    'soil clay mass fraction','none',complete)
  call restart_soil_profile_cache(ncid,flag,'EHUM_vr',EHUM_vr, &
    'microbial residue humification fraction','none',complete)
  call restart_soil_profile_cache(ncid,flag,'EPOC_vr',EPOC_vr, &
    'litter partition coefficient','none',complete)
  call restart_soil_profile_cache(ncid,flag,'NumerSolidThermCond_vr',NumerSolidThermCond_vr, &
    'solid thermal conductivity numerator','MJ m-1 h-1 K-1',complete)
  call restart_soil_profile_cache(ncid,flag,'DenomSolidThermCond_vr',DenomSolidThermCond_vr, &
    'solid thermal conductivity denominator','none',complete)
  call restart_soil_profile_cache(ncid,flag,'ThetaSat_vr',ThetaSat_vr, &
    'water content at conductivity threshold','m3 m-3',complete)
  if(use_ecosys_drainage)then
    ! This policy derives air-entry pressure from the conductivity table. It
    ! must travel with ThetaSat/table caches rather than be reset by STARTS.
    call restart_soil_profile_cache(ncid,flag,'PSISoilAirEntry_vr',PSISoilAirEntry_vr, &
      'matric pressure at drainage conductivity threshold','MPa',complete)
  endif
  call restart_soil_profile_cache(ncid,flag,'PathLenMacPore_vr',PathLenMacPore_vr, &
    'distance between macropores','m',complete)
  call restart_soil_profile_cache(ncid,flag,'MacPoreRadius_vr',MacPoreRadius_vr, &
    'macropore radius','m',complete)
  call restart_soil_profile_cache(ncid,flag,'HydroCondMacP_vr',HydroCondMacP_vr, &
    'macropore hydraulic conductivity','m h-1',complete)

  present_cache=.true.
  if(flag=='read')call ncd_inqvid(ncid,'HydroCond_3D',varid,vardesc,readvar=present_cache)
  if(present_cache)then
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroCond_3D,conductivity)
    conductivity_data=>conductivity
    call restartvar(ncid,flag,varname='HydroCond_3D',dim1name='column',dim2name='sdim', &
      dim3name='moisture_class',dim4name='levsoi1',long_name='micropore hydraulic conductivity table', &
      units='m h-1',interpinic_flag='skip',data=conductivity_data,missing_value=spval,fill_value=spval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,HydroCond_3D,conductivity)
  else
    complete=.false.
  endif
  present_cache=.true.
  if(flag=='read')call ncd_inqvid(ncid,'MacPoreNumbers_vr',varid,vardesc,readvar=present_cache)
  if(present_cache)then
    if(flag=='write')then
      do L=1,JZ
        call cpcol(flag,NHW,NHE,NVN,NVS,MacPoreNumbers_vr(L,:,:),pore_numbers(:,L))
      enddo
    endif
    pore_number_data=>pore_numbers
    call restartvar(ncid,flag,varname='MacPoreNumbers_vr',dim1name='column',dim2name='levsoi', &
      long_name='macropore number per cross-sectional area',units='m-2',interpinic_flag='skip', &
      data=pore_number_data,missing_value=ispval,fill_value=ispval)
    if(flag=='read')then
      do L=1,JZ
        call cpcol(flag,NHW,NHE,NVN,NVS,MacPoreNumbers_vr(L,:,:),pore_numbers(:,L))
      enddo
    endif
  else
    complete=.false.
  endif
  if(flag=='read' .and. .not.complete)then
    ! Legacy files did not save these coefficients. Retain the old recovery
    ! path, but do not describe that recovery as an exact restart.
    iResetSoilProf_col(NVN:NVS,NHW:NHE)=itrue
    write(iulog,*) 'Restart lacks soil property caches; rebuilding properties, exact continuation is not guaranteed'
  endif
  present_cache=.true.
  if(flag=='read')call ncd_inqvid(ncid,'trcn_irrig_vr',varid,vardesc,readvar=present_cache)
  if(present_cache)then
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcn_irrig_vr,irrigation_nutrients)
    irrigation_data=>irrigation_nutrients
    call restartvar(ncid,flag,varname='trcn_irrig_vr',dim1name='column',dim2name='nonband_nutrients', &
      dim3name='levsoi',long_name='nutrient concentration at subsurface inflow boundary',units='g m-3', &
      interpinic_flag='skip',data=irrigation_data,missing_value=spval,fill_value=spval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,trcn_irrig_vr,irrigation_nutrients)
  else
    write(iulog,*) 'Restart lacks subsurface nutrient boundary; exact continuation is not guaranteed'
  endif
  ! This aggregate is consumed by hydrology before plant pools are refreshed.
  present_cache=.true.
  if(flag=='read')call ncd_inqvid(ncid,'RootMycoMassElm_vr',varid,vardesc,readvar=present_cache)
  if(present_cache)then
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RootMycoMassElm_vr,root_mass)
    conductivity_data=>root_mass
    call restartvar(ncid,flag,varname='RootMycoMassElm_vr',dim1name='column',dim2name='elmnts', &
      dim3name='rootyps',dim4name='levsoi',long_name='root and mycorrhizal element mass by layer', &
      units='g',interpinic_flag='skip',data=conductivity_data,missing_value=spval,fill_value=spval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,RootMycoMassElm_vr,root_mass)
  else
    write(iulog,*) 'Restart lacks root biomass aggregate; exact continuation is not guaranteed'
  endif
  complete=.true.
  ! HOUR1 stages the completed-hour heat flux for hydrology before resetting it.
  ! The matching water flux TWaterPlantRoot2Soil_vr is already saved above.
  call restart_soil_profile_cache(ncid,flag,'THeatLossRoot2Soil_vr',THeatLossRoot2Soil_vr, &
    'previous-hour root heat uptake','MJ h-1',complete)
  call restart_soil_profile_cache(ncid,flag,'RCH4EcoDmnd_vr',RCH4EcoDmnd_vr, &
    'completed-hour ecosystem methane demand','gC h-1',complete)
  if(flag=='read' .and. .not.complete) &
    write(iulog,*) 'Restart lacks completed-hour root heat or methane demand; exact continuation is not guaranteed'
  call PrintInfo('end '//subname)
  end associate
  end subroutine restartnc_soil_caches

!------------------------------------------------------------------------------------------
  subroutine restart_soil_profile_cache(ncid,flag,name,values,long_name,units,complete,dimension_name)
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  character(len=*), intent(in) :: flag,name,long_name,units
  real(r8), intent(inout) :: values(:,:,:)
  logical, intent(inout) :: complete
  character(len=*), optional, intent(in) :: dimension_name
  character(len=*), parameter :: subname='restart_soil_profile_cache'
  real(r8), target :: buffer(bounds%begc:bounds%endc,size(values,1))
  real(r8), pointer :: data(:,:)
  type(var_desc_t) :: vardesc
  integer :: varid
  logical :: present_cache
  character(len=8) :: level_dim

  associate(NHW => bounds%NHW, NHE => bounds%NHE, NVN => bounds%NVN, NVS => bounds%NVS)
  call PrintInfo('beg '//subname)
  level_dim='levsoi'
  if(size(values,1)==JZ+1)level_dim='levsoi1'
  if(present(dimension_name))level_dim=dimension_name
  present_cache=.true.
  if(flag=='read')call ncd_inqvid(ncid,name,varid,vardesc,readvar=present_cache)
  if(present_cache)then
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,values,buffer)
    data=>buffer
    call restartvar(ncid,flag,varname=name,dim1name='column',dim2name=trim(level_dim), &
      long_name=long_name,units=units,interpinic_flag='skip',data=data,missing_value=spval,fill_value=spval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,values,buffer)
  else
    complete=.false.
  endif
  call PrintInfo('end '//subname)
  end associate
  end subroutine restart_soil_profile_cache

!------------------------------------------------------------------------------------------

  subroutine restartnc_daily_irrigation(ncid,flag)
  ! Retain the already-decided day's schedule for exact intraday continuation.
  use DebugToolMod, only: PrintInfo
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  character(len=*), intent(in) :: flag
  character(len=*), parameter :: subname='restartnc_daily_irrigation'
  real(r8), target :: hourly_irrigation(bounds%begc:bounds%endc,24)
  real(r8), target :: irrigation_depth(bounds%begc:bounds%endc)
  real(r8), pointer :: hourly_data(:,:),depth_data(:)
  type(var_desc_t) :: vardesc
  integer :: varid,yr,mon,dom,tod,schedule_day
  logical :: has_schedule,has_depth,transfer_schedule

  call PrintInfo('beg '//subname)
  associate(NHW => bounds%NHW,NHE => bounds%NHE,NVN => bounds%NVN,NVS => bounds%NVS)
  schedule_day=etimer%get_curr_doy()
  transfer_schedule=.true.
  if(flag=='read')then
    call etimer%get_curr_date(yr,mon,dom,tod)
    ! At midnight the new day is decided after loading, using restored water.
    transfer_schedule=tod/=0
    if(transfer_schedule)then
      call ncd_inqvid(ncid,'RRIG_day',varid,vardesc,readvar=has_schedule)
      call ncd_inqvid(ncid,'WDPTH_day',varid,vardesc,readvar=has_depth)
      transfer_schedule=has_schedule .and. has_depth
      if(.not.transfer_schedule .and. any(Lirri_auto_col(NVN:NVS,NHW:NHE)))then
        call endrun('Intraday automatic-irrigation restart lacks RRIG_day/WDPTH_day; use a midnight checkpoint',__LINE__)
      endif
      ! Older prescribed-irrigation runs retain the schedule reread from inputs.
    endif
  endif
  if(transfer_schedule)then
    if(flag=='write')then
      call cpcol(flag,NHW,NHE,NVN,NVS,RRIG(:,schedule_day,:,:),hourly_irrigation)
      call cpcol(flag,NHW,NHE,NVN,NVS,WDPTH(schedule_day,:,:),irrigation_depth)
    endif
    hourly_data => hourly_irrigation
    call restartvar(ncid,flag,varname='RRIG_day',dim1name='column',dim2name='hour',&
      long_name='Hourly irrigation schedule for the checkpoint day',units='m h-1',interpinic_flag='skip',&
      data=hourly_data,missing_value=spval,fill_value=spval)
    depth_data => irrigation_depth
    call restartvar(ncid,flag,varname='WDPTH_day',dim1name='column',&
      long_name='Irrigation application depth for the checkpoint day',units='m',interpinic_flag='skip',&
      data=depth_data,missing_value=spval,fill_value=spval)
    if(flag=='read')then
      call cpcol(flag,NHW,NHE,NVN,NVS,RRIG(:,schedule_day,:,:),hourly_irrigation)
      call cpcol(flag,NHW,NHE,NVN,NVS,WDPTH(schedule_day,:,:),irrigation_depth)
    endif
  endif
  end associate
  call PrintInfo('end '//subname)
  end subroutine restartnc_daily_irrigation
!------------------------------------------------------------------------------------------

  SUBROUTINE restartnc_soil(ncid,flag)
!
!     THIS SUBROUTINE WRITES OUT ALL SOIL VARIABLES TO
!     CHECKPOINT FILES AT THE FREQUENCY GIVEN IN THE OPTIONS
!     FILE SO THAT OUTPUTS FROM EARLIER MODEL RUNS CAN BE USED
!     TO INITIALIZE LATER MODEL RUNS
!

  implicit none
  character(len=*) , intent(in) :: flag
  type(file_desc_t), intent(inout) :: ncid ! netcdf id
  integer :: NHW,NHE,NVN,NVS

  integer :: K,M,NX,NY,N,L,NGL,idsa
!     execution begins here
  integer, target  :: datic_1d(bounds%begc:bounds%endc)
  real(r8), target :: datrc_1d(bounds%begc:bounds%endc)
  real(r8), target :: datrc_2d(bounds%begc:bounds%endc,1:JZ+1)
  real(r8), target, allocatable :: datrc_3d(:,:,:)
  real(r8), target, allocatable :: datrc_4d(:,:,:,:)
  real(r8), target, allocatable :: datrc_5d(:,:,:,:,:)
  integer, pointer :: dat1pr(:)
  real(r8),pointer :: datpr1(:),datpr2(:,:),datpr3(:,:,:)
  real(r8),pointer :: datpr4(:,:,:,:),datpr5(:,:,:,:,:)
  integer :: sz3,sz4,sz5
  integer :: ncols, npfts  
  integer :: litter_band_varid
  type(var_desc_t) :: litter_band_vardesc
  logical :: has_litter_band_demand
  integer :: gas_residual_varid
  type(Var_desc_t) :: gas_residual_vardesc
  logical :: has_gas_residual
  integer :: ndoms, feedback_varid
  type(var_desc_t) :: feedback_vardesc
  logical :: has_heat_deriv, has_heat_ref, has_vap_deriv, has_vap_ref
  integer :: clim_mean_varid
  type(var_desc_t) :: clim_mean_vardesc
  logical :: has_clim_mean_temp
  logical :: has_fertP_soil,has_fertP_band
  logical :: has_gpp_day_sum,has_gpp_prev_day,has_gpp_day_hours
  integer :: salt_boundary_varid
  type(var_desc_t) :: salt_boundary_vardesc
  logical :: has_salt_boundary_conc
  logical :: has_reset_soil_profile
  logical :: has_snow_layers
  ncols = bounds%ncols
  npfts = bounds%npfts

  ndoms =trc_confs%NDOMS
  sz3=max(trc_confs%NGasTracers,  &
          trc_confs%NSolutTracers,  &
          trc_confs%NSaltTracers, & 
          trc_confs%NPrecipTracers,  &
          trc_confs%nxtracers,  & 
          trc_confs%NnutrientTracers,  &
          NumMicrobAutoTrophCmplx,jcplx)
  sz4=NumLiveHeterBioms
  sz5=jsken
  allocate(datrc_3d(bounds%begc:bounds%endc,1:sz3,1:JZ+1))
  allocate(datrc_4d(bounds%begc:bounds%endc,1:sz4,1:sz3,1:JZ+1))    
  allocate(datrc_5d(bounds%begc:bounds%endc,1:sz5,1:sz4,1:sz3,1:JZ+1))    

  NHW = bounds%NHW;NVN = bounds%NVN
  NHE = bounds%NHE;NVS = bounds%NVS

  ! Preserve pending geometry updates: replaying one after loading changes
  ! litter density/porosity even when no new disturbance has occurred.
  if(flag=='read')then
    call ncd_inqvid(ncid,'iResetSoilProf_col',feedback_varid,feedback_vardesc, &
      readvar=has_reset_soil_profile)
    if(has_reset_soil_profile)then
      dat1pr => datic_1d
      call restartvar(ncid, flag, varname='iResetSoilProf_col', dim1name='column', &
        long_name='pending soil profile geometry update', units='none', &
        interpinic_flag='skip',data=dat1pr,missing_value=ispval,fill_value=ispval)
      call cpcol(flag,NHW,NHE,NVN,NVS,iResetSoilProf_col,datic_1d)
    else
      write(iulog,*) 'Restart lacks iResetSoilProf_col; exact continuation is not guaranteed'
    endif
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,iResetSoilProf_col,datic_1d)
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='iResetSoilProf_col', dim1name='column', &
      long_name='pending soil profile geometry update', units='none', &
      interpinic_flag='skip',data=dat1pr,missing_value=ispval,fill_value=ispval)
  endif

  call restartvar(ncid, flag, varname='CRAIN_lnd', &
       long_name='total precipitation', units='m3 d-2', &
       interpinic_flag='skip', data=CRAIN_lnd, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TSedmErossLoss_lnds', &
       long_name='total sediment subsurface flux', units='Mg d-2', &
       interpinic_flag='skip', data=TSedmErossLoss_lnds, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TAREA', &
       long_name='total surface grid area', units='m2', &
       interpinic_flag='skip', data=TAREA, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='HEATIN_lnd', &
       long_name='total surface heat flux', units='MJ d-2', &
       interpinic_flag='skip', data=HEATIN_lnd,missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='SurfGas_O2_lnd', &
       long_name='total surface O2 flux', units='g d-2', &
       interpinic_flag='skip', data=SurfGas_lnd(idg_O2), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='SurfGas_CO2_lnd', &
       long_name='total surface CO2 flux', units='g d-2', &
       interpinic_flag='skip', data=SurfGas_lnd(idg_CO2), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='SurfGas_Ar_lnd', &
       long_name='total surface Ar flux', units='g d-2', &
       interpinic_flag='skip', data=SurfGas_lnd(idg_Ar), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='SurfGas_H2_lnd', &
       long_name='total surface H2 flux', units='g d-2', &
       interpinic_flag='skip', data=SurfGas_lnd(idg_H2), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='SurfGas_CH4_lnd', &
       long_name='total surface CH4 flux', units='g d-2', &
       interpinic_flag='skip', data=SurfGas_lnd(idg_CH4), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='tAmendOrgC_lnd', &
       long_name='total organic C amendment', units='g d-2', &
       interpinic_flag='skip', data=tAmendOrgC_lnd, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TORGN', &
       long_name='total organic N amendment', units='g d-2', &
       interpinic_flag='skip', data=TORGN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TORGP', &
       long_name='total organic P amendment', units='g d-2', &
       interpinic_flag='skip', data=TORGP, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='QH2OLoss_lnds', &
       long_name='total subsurface water flux', units='m3 d-2', &
       interpinic_flag='skip', data=QH2OLoss_lnds, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='CEVAP', &
       long_name='total evaporation', units='m3 d-2', &
       interpinic_flag='skip', data=CEVAP, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='CRUN', &
       long_name='total surface runoff', units='m3 d-2', &
       interpinic_flag='skip', data=CRUN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='HeatOut_lnds', &
       long_name='total subsurface heat flux', units='MJ d-2', &
       interpinic_flag='skip', data=HeatOut_lnds, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='OXYGOU', &
       long_name='total subsurface O2 flux', units='g d-2', &
       interpinic_flag='skip', data=OXYGOU, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TCOU', &
       long_name='total subsurface C flux', units='g d-2', &
       interpinic_flag='skip', data=TOMOU_lnds(ielmc), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TZOU', &
       long_name='total subsurface N flux', units='g d-2', &
       interpinic_flag='skip', data=TOMOU_lnds(ielmn), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TPOU', &
       long_name='total subsurface P flux', units='g d-2', &
       interpinic_flag='skip', data=TOMOU_lnds(ielmp), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TZIN', &
       long_name='total surface N flux', units='g d-2', &
       interpinic_flag='skip', data=TZIN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='TPIN', &
       long_name='total surface P flux', units='g d-2', &
       interpinic_flag='skip', data=TPIN, missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='XCSN', &
       long_name='total LitrFall C', units='g d-2', &
       interpinic_flag='skip', data=Litrfall_lnds(ielmc), missing_value=spval, &
       fill_value=spval)

  call restartvar(ncid, flag, varname='XZSN', &
       long_name='total LitrFall N', units='g d-2', &
       interpinic_flag='skip', data=Litrfall_lnds(ielmn), missing_value=spval, &
       fill_value=spval)

  
  call restartvar(ncid, flag, varname='XPSN', &
       long_name='total LitrFall P', units='g d-2', &
       interpinic_flag='skip', data=Litrfall_lnds(ielmp), missing_value=spval, &
       fill_value=spval)
  
  
  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='iFertNH4Band_col', dim1name='column',&
       long_name='banded NH4 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,iFertNH4Band_col,datic_1d)
  else
    !print*,'iFertNH4Band_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,iFertNH4Band_col,datic_1d)
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='iFertNH4Band_col', dim1name='column',&
       long_name='banded NH4 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval,flag_values=(/0,1/))

  endif

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='NU', dim1name='column',&
       long_name='upper soil layer id', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,NU_col,datic_1d)
  else 
    !print*,'NP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NU_col,datic_1d)
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='NU', dim1name='column',&
       long_name='upper soil layer id', units='none', interpinic_flag='skip', &
       data=dat1pr, missing_value=ispval, fill_value=ispval)
  endif

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IDWaterTable_col', dim1name='column',&
       long_name='water table flag from site file', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,IDWaterTable_col,datic_1d)
  else
    !print*,'IDWaterTable'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,IDWaterTable_col,datic_1d)  
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='IDWaterTable_col', dim1name='column',&
       long_name='water table flag from site file', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
  endif  

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='iFertNO3Band_col', dim1name='column',&
       long_name='banded NO3 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,iFertNO3Band_col,datic_1d)
  else
    !print*,'iFertNO3Band_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,iFertNO3Band_col,datic_1d)  
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='iFertNO3Band_col', dim1name='column',&
       long_name='banded NO3 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
  endif  

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='iFertPO4Band_col', dim1name='column',&
       long_name='banded H2PO4 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,iFertPO4Band_col,datic_1d)
  else
    !print*,'iFertPO4Band_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,iFertPO4Band_col,datic_1d)  
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='iFertPO4Band_col', dim1name='column',&
       long_name='banded H2PO4 fertilizer flag', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
  endif  

  if(flag=='read')then
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='iUreaHydInhibitorType_col', dim1name='column',&
       long_name='urea hydrolysis inhibitor type (1=no,2=yes)', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
    call cpcol(flag,NHW,NHE,NVN,NVS,iUreaHydInhibitorType_col,datic_1d)
  else
    !print*,'iUreaHydInhibitorType_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,iUreaHydInhibitorType_col,datic_1d)  
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='iUreaHydInhibitorType_col', dim1name='column',&
       long_name='urea hydrolysis inhibitor type (1=no,2=yes)', units='none', &
       interpinic_flag='skip', data=dat1pr, missing_value=ispval, &
       fill_value=ispval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='CanopyHeight_col', dim1name='column',&
       long_name='canopy height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeight_col,datrc_1d)
  else
    !print*,'CanopyHeight_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeight_col,datrc_1d)  
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='CanopyHeight_col', dim1name='column',&
       long_name='canopy height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='ZeroPlaneDisplacem_col', dim1name='column',&
       long_name='threshold for plant shoot displacement', units='', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,ZeroPlaneDisplacem_col,datrc_1d)
  else
    !print*,'CanopyHeight_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZeroPlaneDisplacem_col,datrc_1d)  
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='ZeroPlaneDisplacem_col', dim1name='column',&
       long_name='threshold for plant shoot displacement', units='', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='RawIsoTAtm2CanopySinkZ_col', dim1name='column',&
       long_name='above canopy boundary layer resistance', units='h m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,RawIsoTAtm2CanopySinkZ_col,datrc_1d)
  else
    !print*,'CanopyHeight_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RawIsoTAtm2CanopySinkZ_col,datrc_1d)  
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='RawIsoTAtm2CanopySinkZ_col', dim1name='column',&
       long_name='above canopy boundary layer resistance', units='h m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='RainIntceptByCanopy_col', dim1name='column',&
       long_name='net water transfer to whole grid canopy', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,RainIntceptByCanopy_col,datrc_1d)
  else
    !print*,'RainIntceptByCanopy_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RainIntceptByCanopy_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='RainIntceptByCanopy_col', dim1name='column',&
       long_name='net water transfer to whole grid canopy', units='MJ d-2 t-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='TSED_col', dim1name='column',&
       long_name='Erosion rate', units='Mg d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,TSED_col,datrc_1d)
  else
    !print*,'TSED'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TSED_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='TSED_col', dim1name='column',&
       long_name='Erosion rate', units='Mg d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='SoilSurfRoughness_col', dim1name='column',&
       long_name='initial soil surface roughness height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilSurfRoughness_col,datrc_1d)
  else
    !print*,'SoilSurfRoughness_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilSurfRoughness_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='SoilSurfRoughness_col', dim1name='column',&
       long_name='initial soil surface roughness height', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='LWRadCanG_col', dim1name='column',&
       long_name='grid total canopy LW emission', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,LWRadCanG_col,datrc_1d)
  else
    !print*,'LWRadCanG_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LWRadCanG_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='LWRadCanG_col', dim1name='column',&
       long_name='grid total canopy LW emission', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='CanopyBiomWater_col', dim1name='column',&
       long_name='total water stored in dry matter', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyBiomWater_col,datrc_1d)
  else
    !print*,'LWRadCanG_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyBiomWater_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='CanopyBiomWater_col', dim1name='column',&
       long_name='total water stored in dry matter', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='CanopyHeatStor_col', dim1name='column',&
       long_name='total canopy heat content', units='MJ d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeatStor_col,datrc_1d)
  else
    !print*,'LWRadCanG_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeatStor_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='CanopyHeatStor_col', dim1name='column',&
       long_name='total canopy heat content', units='MJ d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='LWRadBySurf_col', dim1name='column',&
       long_name='longwave radiation emitted from ground surface (including snow and litter)', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
    call cpcol(flag,NHW,NHE,NVN,NVS,LWRadBySurf_col,datrc_1d)
  else
    !print*,'LWRadBySurf_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LWRadBySurf_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='LWRadBySurf_col', dim1name='column',&
       long_name='longwave radiation emitted from ground surface  (including snow and litter)', units='MJ d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='Eco_NEE_col', dim1name='column',&
       long_name='total canopy net CO2 exchange', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NEE_col,datrc_1d)
  else
    !print*,'Eco_NEE_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NEE_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='Eco_NEE_col', dim1name='column',&
       long_name='total canopy net CO2 exchange', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='RIB', dim1name='column',&
       long_name='Richardson number', units=' ', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RIB_col,datrc_1d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RIB_col,datrc_1d)  
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='RIB', dim1name='column',&
       long_name='Richardson number', units=' ', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  


  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='WatHeldOnCanopy_col', dim1name='column',&
       long_name='canopy surface water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,WatHeldOnCanopy_col,datrc_1d) 
  else
    !print*,'WatHeldOnCanopy_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,WatHeldOnCanopy_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='WatHeldOnCanopy_col', dim1name='column',&
       long_name='canopy surface water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='VWatStoreCapSurf_col', dim1name='column',&
       long_name='surface water storage capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VWatStoreCapSurf_col,datrc_1d) 
  else
    !print*,'VWatStoreCapSurf_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VWatStoreCapSurf_col,datrc_1d)   
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='VWatStoreCapSurf_col', dim1name='column',&
       long_name='surface water storage capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='QRain_CumYr_col', dim1name='column',&
       long_name='total precipitation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,QRain_CumYr_col,datrc_1d) 
  else
    !print*,'QRain_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,QRain_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='QRain_CumYr_col', dim1name='column',&
       long_name='total precipitation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='CanopyLeafArea_col', dim1name='column',&
       long_name='grid canopy leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafArea_col,datrc_1d) 
  else
    !print*,'CanopyLeafArea_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafArea_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='CanopyLeafArea_col', dim1name='column',&
       long_name='grid canopy leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='StemArea_col', dim1name='column',&
       long_name='grid canopy stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,StemArea_col,datrc_1d) 
  else
    !print*,'StemArea_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,StemArea_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='StemArea_col', dim1name='column',&
       long_name='total canopy stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d  
    call restartvar(ncid, flag, varname='PPT', dim1name='column',&
       long_name='total plant population', units='# d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PlantPopu_col,datrc_1d) 
  else
    !print*,'PPT'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PlantPopu_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='PPT', dim1name='column',&
       long_name='total plant population', units='# d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VLWatHeldCapSurf_col', dim1name='column',&
       long_name='soil surface water retention capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatHeldCapSurf_col,datrc_1d) 
  else
    !print*,'VLWatHeldCapSurf_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatHeldCapSurf_col,datrc_1d)   
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='VLWatHeldCapSurf_col', dim1name='column',&
       long_name='soil surface water retention capacity', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='SoiSurfRoughness', dim1name='column',&
       long_name='soil surface roughness height for calculating runoff velocity', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoiSurfRoughness,datrc_1d) 
  else
    !print*,'SoiSurfRoughness'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoiSurfRoughness,datrc_1d)   
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='SoiSurfRoughness', dim1name='column',&
       long_name='soil surface roughness height for calculating runoff velocity', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)    
  endif  

  
  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:trc_confs%NGasTracers)    
    call restartvar(ncid, flag, varname='GasHydroLoss_cumflx_col', dim1name='column',&
       dim2name='gastrcs',long_name='cumulative gas loss via hydrological pathway', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GasHydroLoss_cumflx_col,datrc_2d) 
  else
    !print*,'UCH4G'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GasHydroLoss_cumflx_col,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:trc_confs%NGasTracers)          
    call restartvar(ncid, flag, varname='GasHydroLoss_cumflx_col', dim1name='column',&
       dim2name='gastrcs',long_name='cumulative gas loss via hydrological pathway', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  ! EcoSIM_ZM: retain the cumulative gas-transport diagnostic through restart.
  ! This is a bookkeeping counter, not a gas pool or a new physical source.
  ! Older checkpoints cannot reconstruct it; preserve their legacy zero start
  ! and warn, rather than inventing cumulative history from current pools.
  if(flag=='read')then
    call ncd_inqvid(ncid,'Gas_Prod_TP_cumRes_col',gas_residual_varid,gas_residual_vardesc, &
      readvar=has_gas_residual)
    if(has_gas_residual)then
      datpr2 => datrc_2d(1:ncols,1:trc_confs%NGasTracers)
      call restartvar(ncid, flag, varname='Gas_Prod_TP_cumRes_col', dim1name='column', &
        dim2name='gastrcs',long_name='Cumulative gas production/transport residual', &
        units='g d-2',interpinic_flag='skip',data=datpr2,missing_value=spval,fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,Gas_Prod_TP_cumRes_col,datrc_2d)
    else
      Gas_Prod_TP_cumRes_col=0._r8
      write(iulog,*) 'Legacy checkpoint: cumulative gas transport residual starts at zero'
    endif
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Gas_Prod_TP_cumRes_col,datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:trc_confs%NGasTracers)
    call restartvar(ncid, flag, varname='Gas_Prod_TP_cumRes_col', dim1name='column', &
      dim2name='gastrcs',long_name='Cumulative gas production/transport residual', &
      units='g d-2',interpinic_flag='skip',data=datpr2,missing_value=spval,fill_value=spval)
  endif

  if(flag=='read')then
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='Hydroloss_NH4_cumflx_col', dim1name='column',&
       long_name='cumulative NH4 loss via hydrological pathway', units='gN d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Hydroloss_NH4_cumflx_col,datrc_1d) 
  else    
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Hydroloss_NH4_cumflx_col,datrc_1d) 
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='Hydroloss_NH4_cumflx_col', dim1name='column',&
       long_name='cumulative NH4 loss via hydrological pathway', units='gN d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='Hydroloss_NO3_cumflx_col', dim1name='column',&
       long_name='cumulative NO3 loss via hydrological pathway', units='gN d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Hydroloss_NO3_cumflx_col,datrc_1d) 
  else
    
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Hydroloss_NO3_cumflx_col,datrc_1d) 
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='Hydroloss_NO3_cumflx_col', dim1name='column',&
       long_name='cumulative NO3 loss via hydrological pathway', units='gN d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='CO2byFire_CumYr_col', dim1name='column',&
       long_name='total CO2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CO2byFire_CumYr_col,datrc_1d) 
  else
    !print*,'CO2byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CO2byFire_CumYr_col,datrc_1d) 
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='CO2byFire_CumYr_col', dim1name='column',&
       long_name='total CO2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='CH4byFire_CumYr_col', dim1name='column',&
       long_name='total CH4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CH4byFire_CumYr_col,datrc_1d) 
  else
    !print*,'CH4byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CH4byFire_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='CH4byFire_CumYr_col', dim1name='column',&
       long_name='total CH4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='O2byFire_CumYr_col', dim1name='column',&
       long_name='total O2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,O2byFire_CumYr_col,datrc_1d) 
  else
    !print*,'O2byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,O2byFire_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='O2byFire_CumYr_col', dim1name='column',&
       long_name='total O2 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='N2ObyFire_CumYr_col', dim1name='column',&
       long_name='total N2O flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,N2ObyFire_CumYr_col,datrc_1d) 
  else
    !print*,'N2ObyFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,N2ObyFire_CumYr_col,datrc_1d)
    datpr1 => datrc_1d             
    call restartvar(ncid, flag, varname='N2ObyFire_CumYr_col', dim1name='column',&
       long_name='total N2O flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='NH3byFire_CumYr_col', dim1name='column',&
       long_name='total NH3 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,NH3byFire_CumYr_col,datrc_1d) 
  else
    !print*,'NH3byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NH3byFire_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='NH3byFire_CumYr_col', dim1name='column',&
       long_name='total NH3 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='PO4byFire_CumYr_col', dim1name='column',&
       long_name='total PO4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PO4byFire_CumYr_col,datrc_1d) 
  else
    !print*,'PO4byFire_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PO4byFire_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='PO4byFire_CumYr_col', dim1name='column',&
       long_name='total PO4 flux from fire', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='AmendC_CumYr_flx_col', dim1name='column',&
       long_name='total C amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,AmendC_CumYr_flx_col,datrc_1d) 
  else
    !print*,'AmendC_CumYr_flx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,AmendC_CumYr_flx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='AmendC_CumYr_flx_col', dim1name='column',&
       long_name='total C amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='FertN_Flx_CumYr_col', dim1name='column',&
       long_name='total fertilizer N amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,FertN_Flx_CumYr_col,datrc_1d) 
  else
    !print*,'FertN_Flx_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertN_Flx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='FertN_Flx_CumYr_col', dim1name='column',&
       long_name='total fertilizer N amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='FerP_Flx_CumYr_col', dim1name='column',&
       long_name='total fertilizer P amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,FerP_Flx_CumYr_col,datrc_1d) 
  else
    !print*,'FerP_Flx_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FerP_Flx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='FerP_Flx_CumYr_col', dim1name='column',&
       long_name='total fertilizer P amendment', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='H2OLoss_CumYr_col', dim1name='column',&
       long_name='total subsurface water flux', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,H2OLoss_CumYr_col,datrc_1d) 
  else
    !print*,'H2OLoss_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,H2OLoss_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='H2OLoss_CumYr_col', dim1name='column',&
       long_name='total subsurface water flux', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='QEvap_CumYr_col', dim1name='column',&
       long_name='total evaporation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,QEvap_CumYr_col,datrc_1d) 
  else
    !print*,'QEvap_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,QEvap_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='QEvap_CumYr_col', dim1name='column',&
       long_name='total evaporation', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='Qrunoff_CumYr_col', dim1name='column',&
       long_name='total surface runoff', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Qrunoff_CumYr_col,datrc_1d) 
  else
    !print*,'Qrunoff_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Qrunoff_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='Qrunoff_CumYr_col', dim1name='column',&
       long_name='total surface runoff', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSufDOCFlx_col', dim1name='column',&
       long_name='total surface DOC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOCFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDOCFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOCFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDOCFlx_col', dim1name='column',&
       long_name='total surface DOC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSubsDOCFlx_col', dim1name='column',&
       long_name='total subsurface DOC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDOCFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDOCFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDOCFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSubsDOCFlx_col', dim1name='column',&
       long_name='total subsurface DOC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSufDONFlx_col', dim1name='column',&
       long_name='total surface DON flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDONFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDONFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDONFlx_col,datrc_1d) 
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDONFlx_col', dim1name='column',&
       long_name='total surface DON flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    

  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSubsDONFlx_col', dim1name='column',&
       long_name='total subsurface DON flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDONFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDONFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDONFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSubsDONFlx_col', dim1name='column',&
       long_name='total subsurface DON flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSufDOPFlx_col', dim1name='column',&
       long_name='total surface DOP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOPFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDOPFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDOPFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDOPFlx_col', dim1name='column',&
       long_name='total surface DOP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSubsDOPFlx_col', dim1name='column',&
       long_name='total subsurface DOP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDOPFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDOPFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDOPFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSubsDOPFlx_col', dim1name='column',&
       long_name='total subsurface DOP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='HydroSufDICFlx_col', dim1name='column',&
       long_name='total surface DIC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDICFlx_col,datrc_1d) 
  else
    !print*,'HydroSufDICFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDICFlx_col,datrc_1d)   
    datpr1 => datrc_1d          
    call restartvar(ncid, flag, varname='HydroSufDICFlx_col', dim1name='column',&
       long_name='total surface DIC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSubsDICFlx_col', dim1name='column',&
       long_name='total subsurface DIC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDICFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDICFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDICFlx_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSubsDICFlx_col', dim1name='column',&
       long_name='total subsurface DIC flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSufDINFlx_CumYr_col', dim1name='column',&
       long_name='total surface DIN flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDINFlx_CumYr_col,datrc_1d) 
  else
    !print*,'HydroSufDINFlx_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDINFlx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSufDINFlx_CumYr_col', dim1name='column',&
       long_name='total surface DIN flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSubsDINFlx_col', dim1name='column',&
       long_name='total subsurface DIN flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDINFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDINFlx_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDINFlx_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSubsDINFlx_col', dim1name='column',&
       long_name='total subsurface DIN flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSufDIPFlx_CumYr_col', dim1name='column',&
       long_name='total surface DIP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDIPFlx_CumYr_col,datrc_1d) 
  else
    !print*,'HydroSufDIPFlx_CumYr_col'
   if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,HydroSufDIPFlx_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSufDIPFlx_CumYr_col', dim1name='column',&
       long_name='total surface DIP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='HydroSubsDIPFlx_col', dim1name='column',&
       long_name='total subsurface DIP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDIPFlx_col,datrc_1d) 
  else
    !print*,'HydroSubsDIPFlx_col'
   if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroSubsDIPFlx_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='HydroSubsDIPFlx_col', dim1name='column',&
       long_name='total subsurface DIP flux', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumPlantChemElms)            
    call restartvar(ncid, flag, varname='LiterfalOrgM_col', dim1name='column',&
       dim2name='elmnts',long_name='total litrFall organic matter', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,LiterfalOrgM_col(:,1:JY0,1:JX0),datrc_2d) 
  else
!    print*,'LiterfalOrgC_col',size(LiterfalOrgM_col,1),size(LiterfalOrgM_col,2)
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LiterfalOrgM_col(:,1:JY0,1:JX0),datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='LiterfalOrgM_col', dim1name='column',&
       dim2name='elmnts',long_name='total litrFall organic matter', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumPlantChemElms)            
    call restartvar(ncid, flag, varname='EcoHavstElmnt_CumYr_col', dim1name='column',&
       dim2name='elmnts',long_name='total harvested organic matter', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,EcoHavstElmnt_CumYr_col(:,1:JY0,1:JX0),datrc_2d) 
  else
!    print*,'LiterfalOrgC_col',size(LiterfalOrgM_col,1),size(LiterfalOrgM_col,2)
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,EcoHavstElmnt_CumYr_col(:,1:JY0,1:JX0),datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:NumPlantChemElms)
    call restartvar(ncid, flag, varname='EcoHavstElmnt_CumYr_col', dim1name='column',&
       dim2name='elmnts',long_name='total harvested organic matter', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='RootResp_CumYr_col', dim1name='column',&
       long_name='total soil autotrophic respiration', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,RootResp_CumYr_col,datrc_1d) 
  else
    !print*,'RootResp_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RootResp_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='RootResp_CumYr_col', dim1name='column',&
       long_name='total soil autotrophic respiration', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='QDrain_col', dim1name='column',&
       long_name='total water drainage below root zone', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,QDrain_col,datrc_1d) 
  else
    !print*,'QDrain_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,QDrain_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='QDrain_col', dim1name='column',&
       long_name='total water drainage below root zone', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ZDRAIN', dim1name='column',&
       long_name='total N drainage below root zone', units='gN d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ZDRAIN_col,datrc_1d) 
  else
    !print*,'ZDRAIN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZDRAIN_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ZDRAIN', dim1name='column',&
       long_name='total N drainage below root zone', units='gN d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='PDRAIN', dim1name='column',&
       long_name='total P drainage below root zone', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PDRAIN_col,datrc_1d) 
  else
    !print*,'PDRAIN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PDRAIN_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='PDRAIN', dim1name='column',&
       long_name='total P drainage below root zone', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DPNH4', dim1name='column',&
       long_name='total depth of NH4 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthNH4_col,datrc_1d) 
  else
    !print*,'DPNH4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthNH4_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DPNH4', dim1name='column',&
       long_name='total depth of NH4 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DPNO3', dim1name='column',&
       long_name='total depth of NO3 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthNO3_col,datrc_1d) 
  else
    !print*,'DPNO3'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthNO3_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DPNO3', dim1name='column',&
       long_name='total depth of NO3 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DPPO4', dim1name='column',&
       long_name='total depth of PO4 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthPO4_col,datrc_1d) 
  else
    !print*,'DPPO4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandDepthPO4_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DPPO4', dim1name='column',&
       long_name='total depth of PO4 band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='SoilDetachability4Erosion1', dim1name='column',&
       long_name='soil detachment', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilDetachability4Erosion1,datrc_1d) 
  else
    !print*,'SoilDetachability4Erosion1'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilDetachability4Erosion1,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='SoilDetachability4Erosion1', dim1name='column',&
       long_name='soil detachability', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='SoilDetachability4Erosion2', dim1name='column',&
       long_name='soil detachment', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilDetachability4Erosion2,datrc_1d) 
  else
    !print*,'SoilDetachability4Erosion2'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilDetachability4Erosion2,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='SoilDetachability4Erosion2', dim1name='column',&
       long_name='soil detachment', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='CER_col', dim1name='column',&
       long_name='microbial C  erosion', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CER_col,datrc_1d) 
  else
    !print*,'CER'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CER_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='CER_col', dim1name='column',&
       long_name='microbial C  erosion', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='XER_col', dim1name='column',&
       long_name='soil detachment/deposition', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,XER_col,datrc_1d) 
  else
    !print*,'XER'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,XER_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='XER_col', dim1name='column',&
       long_name='soil detachment/deposition', units='none', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='SedmErossLoss_CumYr_col', dim1name='column',&
       long_name='total sediment subsurface flux', units='Mg d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SedmErossLoss_CumYr_col,datrc_1d) 
  else
    !print*,'SedmErossLoss_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SedmErossLoss_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='SedmErossLoss_CumYr_col', dim1name='column',&
       long_name='total sediment subsurface flux', units='Mg d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ROWSpaceNH4_col', dim1name='column',&
       long_name='row spacing of NH4 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ROWSpaceNH4_col,datrc_1d) 
  else
    !print*,'ROWSpaceNH4_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROWSpaceNH4_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ROWSpaceNH4_col', dim1name='column',&
       long_name='row spacing of NH4 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ROWSpaceNO3Band_col', dim1name='column',&
       long_name='row spacing of NO3 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ROWSpaceNO3Band_col,datrc_1d) 
  else
    !print*,'ROWSpaceNO3Band_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROWSpaceNO3Band_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ROWSpaceNO3Band_col', dim1name='column',&
       long_name='row spacing of NO3 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ROWSpacePO4_col', dim1name='column',&
       long_name='row spacing of PO4 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ROWSpacePO4_col,datrc_1d) 
  else
    !print*,'ROWSpacePO4_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ROWSpacePO4_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ROWSpacePO4_col', dim1name='column',&
       long_name='row spacing of PO4 fertilizer band', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='ExtWaterTableLocDepz_col', dim1name='column',&
       long_name='initial external water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTableLocDepz_col,datrc_1d) 
  else
    !print*,'ExtWaterTableLocDepz_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTableLocDepz_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='ExtWaterTableLocDepz_col', dim1name='column',&
       long_name='initial external water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DTBLD', dim1name='column',&
       long_name='depth of artificial water table adjusted for elevation', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TileWaterTableLocDepz_col,datrc_1d) 
  else
    !print*,'DTBLD'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TileWaterTableLocDepz_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DTBLD', dim1name='column',&
       long_name='depth of artificial water table adjusted for elevation', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_NBP_CumYr_col', dim1name='column',&
       long_name='total NBP', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NBP_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_NBP_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NBP_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_NBP_CumYr_col', dim1name='column',&
       long_name='total NBP', units='g d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='VLitR_col', dim1name='column',&
       long_name='surface litter volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLitR_col,datrc_1d) 
  else
    !print*,'VLitR_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLitR_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='VLitR_col', dim1name='column',&
       long_name='surface litter volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='SED_col', dim1name='column',&
       long_name='sediment transport', units='Mg d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SED_col,datrc_1d) 
  else
    !print*,'SED'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SED_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='SED_col', dim1name='column',&
       long_name='sediment transport', units='Mg d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_GPP_CumYr_col', dim1name='column',&
       long_name='ecosystem GPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_GPP_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_GPP_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_GPP_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_GPP_CumYr_col', dim1name='column',&
       long_name='ecosystem GPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  ! These three fields form one daily clock. Legacy checkpoints cannot
  ! reconstruct yesterday or a partial day's missing hours, so start all at zero.
  has_gpp_day_sum=.true.;has_gpp_prev_day=.true.;has_gpp_day_hours=.true.
  if(flag=='read')then
    call ncd_inqvid(ncid,'GPP_DaySum_col',feedback_varid,feedback_vardesc,readvar=has_gpp_day_sum)
    call ncd_inqvid(ncid,'GPP_PrevDay_col',feedback_varid,feedback_vardesc,readvar=has_gpp_prev_day)
    call ncd_inqvid(ncid,'GPP_DayHours_col',feedback_varid,feedback_vardesc,readvar=has_gpp_day_hours)
    GPP_DaySum_col=0._r8
    GPP_PrevDay_col=0._r8
    GPP_DayHours_col=0
  endif
  if(has_gpp_day_sum .and. has_gpp_prev_day .and. has_gpp_day_hours)then
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GPP_DaySum_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='GPP_DaySum_col', dim1name='column',&
      long_name='Accumulated GPP for the current day', units='g C d-2', interpinic_flag='skip',&
      data=datpr1, missing_value=spval, fill_value=spval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,GPP_DaySum_col,datrc_1d)

    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GPP_PrevDay_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='GPP_PrevDay_col', dim1name='column',&
      long_name='GPP for the previous completed day', units='g C d-2', interpinic_flag='skip',&
      data=datpr1, missing_value=spval, fill_value=spval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,GPP_PrevDay_col,datrc_1d)

    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GPP_DayHours_col,datic_1d)
    dat1pr => datic_1d
    call restartvar(ncid, flag, varname='GPP_DayHours_col', dim1name='column',&
      long_name='Hours accumulated in the current daily GPP sum', units='hours', interpinic_flag='skip',&
      data=dat1pr, missing_value=ispval, fill_value=ispval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,GPP_DayHours_col,datic_1d)
  endif

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_AutoR_CumYr_col', dim1name='column',&
       long_name='ecosystem autotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_AutoR_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_AutoR_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_AutoR_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_AutoR_CumYr_col', dim1name='column',&
       long_name='ecosystem autotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_NPP_CumYr_col', dim1name='column',&
       long_name='ecosystem NPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NPP_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_NPP_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_NPP_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_NPP_CumYr_col', dim1name='column',&
       long_name='ecosystem NPP', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Eco_HR_CumYr_col', dim1name='column',&
       long_name='ecosystem heterotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Eco_HR_CumYr_col,datrc_1d) 
  else
    !print*,'Eco_HR_CumYr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Eco_HR_CumYr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Eco_HR_CumYr_col', dim1name='column',&
       long_name='ecosystem heterotrophic respiration', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Air_Heat_Latent_store_col', dim1name='column',&
       long_name='total latent heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Air_Heat_Latent_store_col,datrc_1d) 
  else
    !print*,'Air_Heat_Latent_store_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Air_Heat_Latent_store_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Air_Heat_Latent_store_col', dim1name='column',&
       long_name='total latent heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='Air_Heat_Sens_store_col', dim1name='column',&
       long_name='total sensible heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,Air_Heat_Sens_store_col,datrc_1d) 
  else
    !print*,'Air_Heat_Sens_store_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,Air_Heat_Sens_store_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='Air_Heat_Sens_store_col', dim1name='column',&
       long_name='total sensible heat flux x boundary layer resistance', units='MJ m-1', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)      
  endif

  ! Preserve the climatological temperature used for atmospheric gas densities.
  ! Legacy checkpoints predate this field; recover the same site-input mean
  ! used by STARTS instead of retaining an uninitialized Kelvin temperature.
  has_clim_mean_temp=.true.
  if(flag=='read')then
    call ncd_inqvid(ncid,'TairKClimMean_col',clim_mean_varid,clim_mean_vardesc,readvar=has_clim_mean_temp)
    if(.not.has_clim_mean_temp)then
      TairKClimMean_col(NVN:NVS,NHW:NHE)=units%Celcius2Kelvin(ATCAI_col(NVN:NVS,NHW:NHE))
    endif

  endif
  if(has_clim_mean_temp)then
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TairKClimMean_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='TairKClimMean_col', dim1name='column',&
      long_name='Mean annual air temperature', units='K', interpinic_flag='skip',&
      data=datpr1, missing_value=spval, fill_value=spval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,TairKClimMean_col,datrc_1d)
  endif

  ! Preserve the paired hourly feedback history. Legacy restarts start with
  ! zero derivative and build a fresh resistance estimate after one hour.
  if(flag=='read')then
    call ncd_inqvid(ncid,'dSensHeatdTKQ_col',feedback_varid,feedback_vardesc,readvar=has_heat_deriv)
    call ncd_inqvid(ncid,'TKQHeatRef_col',feedback_varid,feedback_vardesc,readvar=has_heat_ref)
    dSensHeatdTKQ_col=0._r8
    TKQHeatRef_col=0._r8
    if(has_heat_deriv .and. has_heat_ref)then
      datpr1 => datrc_1d
      call restartvar(ncid, flag, varname='dSensHeatdTKQ_col', dim1name='column',&
        long_name='Resistance estimate of weighted sensible-heat derivative', units='MJ m-1 K-1',&
        interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,dSensHeatdTKQ_col,datrc_1d)
      datpr1 => datrc_1d
      call restartvar(ncid, flag, varname='TKQHeatRef_col', dim1name='column',&
        long_name='Canopy air temperature used for the last completed hour', units='K',&
        interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,TKQHeatRef_col,datrc_1d)
    endif
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,dSensHeatdTKQ_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='dSensHeatdTKQ_col', dim1name='column',&
      long_name='Resistance estimate of weighted sensible-heat derivative', units='MJ m-1 K-1',&
      interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TKQHeatRef_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='TKQHeatRef_col', dim1name='column',&
      long_name='Canopy air temperature used for the last completed hour', units='K',&
      interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
  endif

  ! Preserve the paired hourly feedback history. Legacy restarts start with
  ! zero derivative and build a fresh resistance estimate after one hour.
  if(flag=='read')then
    call ncd_inqvid(ncid,'dLatHeatdVPQ_col',feedback_varid,feedback_vardesc,readvar=has_vap_deriv)
    call ncd_inqvid(ncid,'VPQLatHeatRef_col',feedback_varid,feedback_vardesc,readvar=has_vap_ref)
    dLatHeatdVPQ_col=0._r8
    VPQLatHeatRef_col=-1._r8
    if(has_vap_deriv .and. has_vap_ref)then
      datpr1 => datrc_1d
      call restartvar(ncid, flag, varname='dLatHeatdVPQ_col', dim1name='column',&
        long_name='Resistance estimate of weighted latent-energy derivative', units='MJ m2 ton-1',&
        interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,dLatHeatdVPQ_col,datrc_1d)
      datpr1 => datrc_1d
      call restartvar(ncid, flag, varname='VPQLatHeatRef_col', dim1name='column',&
        long_name='Canopy vapor concentration used for the last completed hour', units='ton m-3',&
        interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,VPQLatHeatRef_col,datrc_1d)
    endif
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,dLatHeatdVPQ_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='dLatHeatdVPQ_col', dim1name='column',&
      long_name='Resistance estimate of weighted latent-energy derivative', units='MJ m2 ton-1',&
      interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VPQLatHeatRef_col,datrc_1d)
    datpr1 => datrc_1d
    call restartvar(ncid, flag, varname='VPQLatHeatRef_col', dim1name='column',&
      long_name='Canopy vapor concentration used for the last completed hour', units='ton m-3',&
      interpinic_flag='skip', data=datpr1, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr1 => datrc_1d            
    call restartvar(ncid, flag, varname='DayLensCurr_col', dim1name='column',&
       long_name='daylength', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DayLensCurr_col,datrc_1d) 
  else
    !print*,'DayLensCurr_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DayLensCurr_col,datrc_1d)   
    datpr1 => datrc_1d              
    call restartvar(ncid, flag, varname='DayLensCurr_col', dim1name='column',&
       long_name='daylength', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='DayLenthPrev_col', dim1name='column',&
       long_name='daylength of previous day', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DayLenthPrev_col,datrc_1d) 
  else
    !print*,'DayLenthPrev_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DayLenthPrev_col,datrc_1d)   
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='DayLenthPrev_col', dim1name='column',&
       long_name='daylength of previous day', units='h', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d      
    call restartvar(ncid, flag, varname='ExtWaterTable', dim1name='column',&
       long_name='current external water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTable_col,datrc_1d) 
  else
    !print*,'ExtWaterTable'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ExtWaterTable_col,datrc_1d)   
    datpr1 => datrc_1d        
    call restartvar(ncid, flag, varname='ExtWaterTable', dim1name='column',&
       long_name='current external water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='TileWaterTable_col', dim1name='column',&
       long_name='artificial water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TileWaterTable_col,datrc_1d) 
  else
    !print*,'TileWaterTable_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TileWaterTable_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='TileWaterTable_col', dim1name='column',&
       long_name='artificial water table depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:micpar%NumOfLitrCmplxs)    
    call restartvar(ncid, flag, varname='RC0', dim1name='column',dim2name='nlitromcomplx',&
       long_name='surface OM in each complex', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SurfLitterC_col(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'RC0',size(RC0,1),size(RC0,2),jcplx
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SurfLitterC_col(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:micpar%NumOfLitrCmplxs)    
    call restartvar(ncid, flag, varname='RC0', dim1name='column',dim2name='nlitromcomplx',&
       long_name='surface OM in each complex', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumDrySnoWE_col', dim1name='column',&
       long_name='snow volume in snowpack (water equivalent)', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumDrySnoWE_col,datrc_1d) 
  else
    !print*,'VcumDrySnoWE_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumDrySnoWE_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumDrySnoWE_col', dim1name='column',&
       long_name='snow volume in snowpack (water equivalent)', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumWatSnow_col', dim1name='column',&
       long_name='water volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumWatSnow_col,datrc_1d) 
  else
    !print*,'VcumWatSnow_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumWatSnow_col,datrc_1d) 
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumWatSnow_col', dim1name='column',&
       long_name='water volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumIceSnow_col', dim1name='column',&
       long_name='ice volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumIceSnow_col,datrc_1d) 
  else
    !print*,'VcumIceSnow_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumIceSnow_col,datrc_1d)   
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumIceSnow_col', dim1name='column',&
       long_name='ice volume in snowpack', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumSnoDWI_col', dim1name='column',&
       long_name='total snowpack volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VcumSnoDWI_col,datrc_1d) 
  else
    !print*,'VcumSnoDWI_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VcumSnoDWI_col,datrc_1d) 
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='VcumSnoDWI_col', dim1name='column',&
       long_name='total snowpack volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='SnowDepth', dim1name='column',&
       long_name='snowpack depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SnowDepth_col,datrc_1d) 
  else
    !print*,'SnowDepth'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnowDepth_col,datrc_1d) 
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='SnowDepth', dim1name='column',&
       long_name='snowpack depth', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='EnergyImpact4Erosion_col', dim1name='column',&
       long_name='cumulative rainfall energy impact on soil surface', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,EnergyImpact4Erosion_col,datrc_1d) 
  else
    !print*,'EnergyImpact4Erosion_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,EnergyImpact4Erosion_col,datrc_1d)    
    datpr1 => datrc_1d    
    call restartvar(ncid, flag, varname='EnergyImpact4Erosion_col', dim1name='column',&
       long_name='cumulative rainfall energy impact on soil surface', units='m', &
       interpinic_flag='skip', data=datpr1, missing_value=spval, &
       fill_value=spval)    
   endif 

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLDrySnoWE', dim1name='column',dim2name='levsno',&
       long_name='water equivalent dry snow in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLDrySnoWE_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLDrySnoWE'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLDrySnoWE_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)    
    call restartvar(ncid, flag, varname='VLDrySnoWE', dim1name='column',dim2name='levsno',&
       long_name='water equivalent dry snow in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLIceSnow', dim1name='column',dim2name='levsno',&
       long_name='snow ice volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLIceSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLIceSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLIceSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)    
    call restartvar(ncid, flag, varname='VLIceSnow', dim1name='column',dim2name='levsno',&
       long_name='snow ice volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLWatSnow', dim1name='column',dim2name='levsno',&
       long_name='snow water volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLWatSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLWatSnow', dim1name='column',dim2name='levsno',&
       long_name='snow water volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLSnoDWIprev_snvr', dim1name='column',dim2name='levsno',&
       long_name='total snow volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSnoDWIprev_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLSnoDWIprev_snvr
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSnoDWIprev_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLSnoDWIprev_snvr', dim1name='column',dim2name='levsno',&
       long_name='total snow volume in snowpack layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnoDensL', dim1name='column',dim2name='levsno',&
       long_name='snowpack density', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SnoDens_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'SnoDensL'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnoDens_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnoDensL', dim1name='column',dim2name='levsno',&
       long_name='snowpack density', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnowThickL_snvr', dim1name='column',dim2name='levsno',&
       long_name='snowpack layer depth', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SnowThickL_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'SnowThickL_snvr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnowThickL_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='SnowThickL_snvr', dim1name='column',dim2name='levsno',&
       long_name='snowpack layer depth', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLHeatCapSnow_snvr', dim1name='column',dim2name='levsno',&
       long_name='snowpack heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLHeatCapSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'VLHeatCapSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLHeatCapSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='VLHeatCapSnow_snvr', dim1name='column',dim2name='levsno',&
       long_name='snowpack heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  ! The initial gas balance is summed before snow physics runs. Preserve its
  ! layer count so restored snow tracers are included in that first sum.
  has_snow_layers=.true.
  if(flag=='read')call ncd_inqvid(ncid,'nsnol_col',feedback_varid,feedback_vardesc,readvar=has_snow_layers)
  if(has_snow_layers)then
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,nsnol_col,datic_1d)
    dat1pr=>datic_1d
    call restartvar(ncid,flag,varname='nsnol_col',dim1name='column', &
      long_name='number of active snow layers',units='none',interpinic_flag='skip', &
      data=dat1pr,missing_value=ispval,fill_value=ispval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,nsnol_col,datic_1d)
  else
    ! Match the completed-hour SnowBalance count without advancing physics.
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        nsnol_col(NY,NX)=count(VLHeatCapSnow_snvr(1:JS,NY,NX)>0._r8)
      ENDDO
    ENDDO
    write(iulog,*) 'Restart lacks nsnol_col; reconstructed from restored snow heat capacity'
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='TKSnow', dim1name='column',dim2name='levsno',&
       long_name='snow temperature', units='K', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TKSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'TKSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TKSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='TKSnow', dim1name='column',dim2name='levsno',&
       long_name='snow temperature', units='K', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='TCSnow', dim1name='column',dim2name='levsno',&
       long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TCSnow_snvr(:,1:JY0,1:JX0),datrc_2d) 
  else
    !print*,'TCSnow'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TCSnow_snvr(:,1:JY0,1:JX0),datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JS)        
    call restartvar(ncid, flag, varname='TCSnow', dim1name='column',dim2name='levsno',&
       long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JS)    
    call restartvar(ncid, flag, varname='trcg_solsml', dim1name='column',dim2name='gastrcs',&
       dim3name='levsno',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,trcg_solsml_snvr(:,:,1:JY0,1:JX0),datrc_3d) 
  else
    !print*,'trcg_solsml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_solsml_snvr(:,:,1:JY0,1:JX0),datrc_3d)   
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JS)        
    call restartvar(ncid, flag, varname='trcg_solsml', dim1name='column',dim2name='gastrcs',&
       dim3name='levsno',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif 

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NnutrientTracers,1:JS)            
    call restartvar(ncid, flag, varname='trcn_solsml', dim1name='column',dim2name='NnutrientTracers',&
       dim3name='levsno',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,trcn_solsml_snvr,datrc_3d) 
  else
    !print*,'trcn_solsml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcn_solsml_snvr,datrc_3d)   
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NnutrientTracers,1:JS)            
    call restartvar(ncid, flag, varname='trcn_solsml', dim1name='column',dim2name='NnutrientTracers',&
       dim3name='levsno',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SoilFracAsMacP_vr', dim1name='column',dim2name='levsoi',&
       long_name='macropore fraction', units='m3/m3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilFracAsMacP_vr,datrc_2d) 
  else
    !print*,'SoilFracAsMacP_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilFracAsMacP_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ)        
    call restartvar(ncid, flag, varname='SoilFracAsMacP_vr', dim1name='column',dim2name='levsoi',&
       long_name='thickness of soil layer', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:3,1:JZ+1)    
    call restartvar(ncid, flag, varname='DLYR_3D', dim1name='column',dim2name='sdim',&
       dim3name='levsoi1',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,DLYR_3D,datrc_3d) 
  else
    !print*,'DLYR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DLYR_3D,datrc_3d)   
    datpr3 => datrc_3d(1:ncols,1:3,1:JZ+1)            
    call restartvar(ncid, flag, varname='DLYR_3D', dim1name='column',dim2name='sdim',&
       dim3name='levsoi1',long_name='snow temperature', units='oC', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='CumDepz2LayBottom_vr', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CumDepz2LayBottom_vr,datrc_2d) 
  else
    !print*,'CumDepz2LayBottom_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CumDepz2LayBottom_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='CumDepz2LayBottom_vr', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='CumSoilThickness_vr', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer from grid surface', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CumSoilThickness_vr,datrc_2d) 
  else
    !print*,'CumSoilThickness_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CumSoilThickness_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='CumSoilThickness_vr', dim1name='column',dim2name='levsoi1',&
       long_name='depth to bottom of soil layer from grid surface', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)        
    call restartvar(ncid, flag, varname='SoiBulkDensityt0_vr', dim1name='column',dim2name='levsoi',&
       long_name='initial bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensityt0_vr,datrc_2d) 
  else
    !print*,'SoiBulkDensityt0_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoiBulkDensityt0_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='SoiBulkDensityt0_vr', dim1name='column',dim2name='levsoi',&
       long_name='initial bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SoiBulkDensity', dim1name='column',dim2name='levsoi1',&
       long_name='soil bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilBulkDensity_vr,datrc_2d) 
  else
    !print*,'SoiBulkDensity'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilBulkDensity_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SoiBulkDensity', dim1name='column',dim2name='levsoi1',&
       long_name='soil bulk density,0=water', units='Mg m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumPlantChemElms,1:JZ+1)    
    call restartvar(ncid, flag, varname='CSoilOrgM_vr', dim1name='column',dim2name='elmnts',dim3name='levsoi1',&
       long_name='soil organic matter content density', units='g kg-1 m-3', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,CSoilOrgM_vr,datrc_3d) 
  else
    !print*,'CORGC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CSoilOrgM_vr,datrc_3d)   
    datpr3 => datrc_3d(1:ncols,1:NumPlantChemElms,1:JZ+1)    
    call restartvar(ncid, flag, varname='CSoilOrgM_vr', dim1name='column',dim2name='elmnts',dim3name='levsoi1',&
       long_name='soil organic matter content density', units='g kg-1 m-3', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='POROS_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil porosity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,POROS_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,POROS_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='POROS_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil porosity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='PH', dim1name='column',dim2name='levsoi1',&
       long_name='soil pH', units='none', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PH_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PH_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='PH', dim1name='column',dim2name='levsoi1',&
       long_name='soil pH', units='none', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicPMass_vr', dim1name='column',dim2name='levsoi1',&
       long_name='mass of soil layer', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicPMass_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicPMass_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicPMass_vr', dim1name='column',dim2name='levsoi1',&
       long_name='mass of soil layer', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilPoreMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume of soil layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilPoreMicP_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilPoreMicP_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilPoreMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume of soil layer', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='PSISoilMatricP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume of soil layer', units='MPa', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,PSISoilMatricP_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PSISoilMatricP_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='PSISoilMatricP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil matric potential', units='MPa', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='THETW_vr', dim1name='column',dim2name='levsoi1',&
       long_name='volumetric water content', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,THETW_vr,datrc_2d) 
  else
    !print*,'POROS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,THETW_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='THETW_vr', dim1name='column',dim2name='levsoi1',&
       long_name='volumetric water content', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='FieldCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='Field capacity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,FieldCapacity_vr,datrc_2d) 
  else
    !print*,'FieldCapacity'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FieldCapacity_vr,datrc_2d)   
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='FieldCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='Field capacity', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='WiltPoint_vr', dim1name='column',dim2name='levsoi1',&
       long_name='Wilting point', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,WiltPoint_vr,datrc_2d)   
  else
    !print*,'WiltPoint_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,WiltPoint_vr,datrc_2d)    
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='WiltPoint_vr', dim1name='column',dim2name='levsoi1',&
       long_name='Wilting point', units='m3 m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SatHydroCondVert_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil vertical saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondVert_vr,datrc_2d)   
  else
    !print*,'SatHydroCondVert_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondVert_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SatHydroCondVert_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil vertical saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SatHydroCondHrzn_vr', dim1name='column',dim2name='levsoi',&
       long_name='soil horizontal saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondHrzn_vr,datrc_2d)   
  else
    !print*,'SatHydroCondHrzn_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SatHydroCondHrzn_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SatHydroCondHrzn_vr', dim1name='column',dim2name='levsoi',&
       long_name='soil horizontal saturated hydraulic conductivity', units='mm h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SAND', dim1name='column',dim2name='levsoi',&
       long_name='soil sand content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SAND_vr,datrc_2d)   
  else
    !print*,'SAND'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SAND_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SAND', dim1name='column',dim2name='levsoi',&
       long_name='soil sand content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SILT', dim1name='column',dim2name='levsoi',&
       long_name='soil silt content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SILT_vr,datrc_2d)   
  else
    !print*,'SILT'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SILT_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='SILT', dim1name='column',dim2name='levsoi',&
       long_name='soil silt content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CLAY', dim1name='column',dim2name='levsoi',&
       long_name='soil clay content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,CLAY_vr,datrc_2d)   
  else
    !print*,'CLAY'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CLAY_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CLAY', dim1name='column',dim2name='levsoi',&
       long_name='soil clay content', units='Mg d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLWatMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicP_vr,datrc_2d)   
  else
    !print*,'VLWatMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicP_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLWatMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLWatMicPX', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore water content before wetting front', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicPX_vr,datrc_2d)     
  else
    !print*,'VLWatMicPX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMicPX_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLWatMicPX', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore water content before wetting front', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLiceMicP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore ice content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMicP_vr,datrc_2d)     
  else
    !print*,'VLiceMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMicP_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLiceMicP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore ice content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLsoiAirP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore air content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLsoiAirP_vr,datrc_2d)     
  else
    !print*,'VLsoiAirP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLsoiAirP_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLsoiAirP', dim1name='column',dim2name='levsoi1',&
       long_name='soil micropore air content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total volume in micropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLMicP_vr,datrc_2d)     
  else
    !print*,'VLMicP_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLMicP_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total volume in micropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicP_vr,datrc_2d)     
  else
    !print*,'VLSoilMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLSoilMicP_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='VLSoilMicP_vr', dim1name='column',dim2name='levsoi1',&
       long_name='micropore volume', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SoilWatAirDry_vr', dim1name='column',dim2name='levsoi1',&
       long_name='air-dry water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SoilWatAirDry_vr,datrc_2d)     
  else
    !print*,'VLSoilMicP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilWatAirDry_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='SoilWatAirDry_vr', dim1name='column',dim2name='levsoi1',&
       long_name='air-dry water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%nxtracers, 1:JZ+1)    
    call restartvar(ncid, flag, varname='trcx_solml_vr', dim1name='column',dim2name='xtracers',&
       dim3name='levsoi1',long_name='exchangeable tracers', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,trcx_solml_vr,datrc_3d)     
  else
    !print*,'trcx_solml'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcx_solml_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%nxtracers, 1:JZ+1)        
    call restartvar(ncid, flag, varname='trcx_solml_vr', dim1name='column',dim2name='xtracers',&
       dim3name='levsoi1',long_name='exchangeable tracers', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='TCS', dim1name='column',dim2name='levsoi1',&
       long_name='soil temperature', units='oC', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,TCS_vr,datrc_2d)     
  else
    !print*,'TCS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TCS_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='TCS', dim1name='column',dim2name='levsoi1',&
       long_name='soil temperature', units='oC', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='TKS', dim1name='column',dim2name='levsoi1',&
       long_name='soil temperature', units='K', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,TKS_vr,datrc_2d)     
  else
    !print*,'TKS'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TKS_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='TKS', dim1name='column',dim2name='levsoi1',&
       long_name='soil temperature', units='K', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  


  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='soil heat capacity', units='MJ d-2 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacity_vr,datrc_2d)     
  else
    !print*,'VHeatCapacity'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapacity_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapacity', dim1name='column',dim2name='levsoi1',&
       long_name='soil heat capacity', units='MJ d-2 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapSolidSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil solid heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapSolidSoil_vr,datrc_2d)     
  else
    !print*,'VHeatCapSolidSoil_vr'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,VHeatCapSolidSoil_vr,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VHeatCapSolidSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='soil solid heat capacity', units='MJ m-3 K-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JZ)    
    call restartvar(ncid, flag, varname='trcg_gasml_vr', dim1name='column',dim2name='gastrcs',&
       dim3name='levsoi',long_name='layer mass of gases', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,trcg_gasml_vr,datrc_3d)     
  else
    !print*,'trcg_gasml_vr'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,trcg_gasml_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JZ)        
    call restartvar(ncid, flag, varname='trcg_gasml_vr', dim1name='column',dim2name='gastrcs',&
       dim3name='levsoi',long_name='layer mass of gases', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JZ)    
    call restartvar(ncid, flag, varname='trcg_root_vr', dim1name='column',dim2name='gastrcs',&
       dim3name='levsoi',long_name='layer mass of gases in roots', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,trcg_root_vr,datrc_3d)     
  else
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,trcg_root_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JZ)        
    call restartvar(ncid, flag, varname='trcg_root_vr', dim1name='column',dim2name='gastrcs',&
       dim3name='levsoi',long_name='layer mass of gases in roots', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)        
    call restartvar(ncid, flag, varname='trcs_solml_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,trcs_solml_vr,datrc_3d)     
  else
    !print*,'trcs_solml_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcs_solml_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)            
    call restartvar(ncid, flag, varname='trcs_solml_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)        
    call restartvar(ncid, flag, varname='trcs_solml_drib_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute dribble mass in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,trcs_solml_drib_vr,datrc_3d)     
  else
    !print*,'trcs_solml_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcs_solml_drib_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)            
    call restartvar(ncid, flag, varname='trcs_solml_drib_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute dribble mass in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)           
    call restartvar(ncid, flag, varname='trcs_soHml_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in macropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,trcs_soHml_vr,datrc_3d)     
  else
    !print*,'trcs_soHml_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcs_soHml_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)            
    call restartvar(ncid, flag, varname='trcs_soHml_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='solute mass in macropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)   

  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JZ+1)    
    call restartvar(ncid, flag, varname='RGasTranspFlxPrev_vr', dim1name='column',dim2name='gastrcs',&
       dim3name='levsoi1',long_name='net gaseous flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RGasTranspFlxPrev_vr,datrc_3d)     
  else
    !print*,'RO2GasXchangePrev_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RGasTranspFlxPrev_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NGasTracers,1:JZ+1)        
    call restartvar(ncid, flag, varname='RGasTranspFlxPrev_vr', dim1name='column',dim2name='gastrcs',&
       dim3name='levsoi1',long_name='net gaseous flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
  endif  


  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RO2AquaSourcePrev_vr', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous O2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RO2AquaSourcePrev_vr,datrc_2d)       
  else
    !print*,'RO2AquaSourcePrev_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RO2AquaSourcePrev_vr,datrc_2d)         
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RO2AquaSourcePrev_vr', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous O2 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCH4L', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous CH4 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RCH4PhysexchPrev_vr,datrc_2d)       
  else
    !print*,'RCH4L'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RCH4PhysexchPrev_vr,datrc_2d)         
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RCH4L', dim1name='column',dim2name='levsoi1',&
       long_name='net aqueous CH4 flux', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoO2DmndResp_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial O2 uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoO2DmndResp_vr,datrc_2d)      
  else
    !print*,'REcoO2DmndResp_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoO2DmndResp_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoO2DmndResp_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial O2 uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNH4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoNH4DmndSoil_vr,datrc_2d)      
  else
    !print*,'REcoNH4DmndSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoNH4DmndSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNH4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNO3DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoNO3DmndSoil_vr,datrc_2d)      
  else
    !print*,'REcoNO3DmndSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoNO3DmndSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNO3DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2EcoUptkSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2EcoUptkSoil_vr,datrc_2d)      
  else
    !print*,'RNO2EcoUptkSoil_vr'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,RNO2EcoUptkSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2EcoUptkSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN2OEcoUptkSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial N2O uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RN2OEcoUptkSoil_vr,datrc_2d)      
  else
    !print*,'RN2OEcoUptkSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RN2OEcoUptkSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RN2OEcoUptkSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial N2O uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH2PO4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoH2PO4DmndSoil_vr,datrc_2d)      
  else
    !print*,'REcoH2PO4DmndSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoH2PO4DmndSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH2PO4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake non-band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH1PO4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in non-band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoH1PO4DmndSoil_vr,datrc_2d)      
  else
    !print*,'REcoH1PO4DmndSoil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoH1PO4DmndSoil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH1PO4DmndSoil_vr', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in non-band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNH4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoNH4DmndBand_vr,datrc_2d)      
  else
    !print*,'REcoNH4DmndBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoNH4DmndBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNH4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NH4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNO3DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoNO3DmndBand_vr,datrc_2d)      
  else
    !print*,'REcoNO3DmndBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoNO3DmndBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoNO3DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO3 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2EcoUptkBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2EcoUptkBand_vr,datrc_2d)      
  else
    !print*,'RNO2EcoUptkBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2EcoUptkBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='RNO2EcoUptkBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial NO2 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH2PO4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoH2PO4DmndBand_vr,datrc_2d)      
  else
    !print*,'REcoH2PO4DmndBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoH2PO4DmndBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH2PO4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='total root + microbial PO4 uptake band', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH1PO4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,REcoH1PO4DmndBand_vr,datrc_2d)      
  else
    !print*,'REcoH1PO4DmndBand_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,REcoH1PO4DmndBand_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)        
    call restartvar(ncid, flag, varname='REcoH1PO4DmndBand_vr', dim1name='column',dim2name='levsoi1',&
       long_name='HPO4 demand in band by all microbial,root,myco populations', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)        
    call restartvar(ncid, flag, varname='RDOMEcoDmndK_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial DOC uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RDOMEcoDmndK_vr,datrc_3d)      
  else
    !print*,'RDOMEcoDmndK_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RDOMEcoDmndK_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='RDOMEcoDmndK_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial DOC uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)        
    call restartvar(ncid, flag, varname='FracBulkSOMC_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='mass fraction of complex OM in bulk OM', units='none', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,FracBulkSOMC_vr,datrc_3d)      
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FracBulkSOMC_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='FracBulkSOMC_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='mass fraction of complex OM in bulk OM', units='none', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='RAcetateEcoDmndK_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial acetate uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RAcetateEcoDmndK_vr,datrc_3d)      
  else
    !print*,'RAcetateEcoDmndK_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RAcetateEcoDmndK_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:jcplx, 1:JZ+1)            
    call restartvar(ncid, flag, varname='RAcetateEcoDmndK_vr', dim1name='column',dim2name='nomcomplx',&
       dim3name='levsoi1',long_name='total root + microbial acetate uptake', units='g d-2 h-1', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)            
    call restartvar(ncid, flag, varname='VLWatMacP', dim1name='column',dim2name='levsoi1',&
       long_name='soil macropore water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMacP_vr,datrc_2d)      
  else
    !print*,'VLWatMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLWatMacP_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)                
    call restartvar(ncid, flag, varname='VLWatMacP', dim1name='column',dim2name='levsoi1',&
       long_name='soil macropore water content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         

  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)                
    call restartvar(ncid, flag, varname='VLiceMacP', dim1name='column',dim2name='levsoi1',&
       long_name='soil macropore ice content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMacP_vr,datrc_2d)      
  else
    !print*,'VLiceMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLiceMacP_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ+1)                
    call restartvar(ncid, flag, varname='VLiceMacP', dim1name='column',dim2name='levsoi1',&
       long_name='soil macropore ice content', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)                
    call restartvar(ncid, flag, varname='TWaterPlantRoot2Soil_vr', dim1name='column',dim2name='levsoi',&
       long_name='water flux from plant roots to soil', units='m3 H2O d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,TWaterPlantRoot2Soil_vr,datrc_2d)      
  else
    !print*,'VLMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TWaterPlantRoot2Soil_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='TWaterPlantRoot2Soil_vr', dim1name='column',dim2name='levsoi',&
       long_name='water flux from plant roots to soil', units='m3 H2O d-2 h-1', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)                
    call restartvar(ncid, flag, varname='VLMacroPores_vr', dim1name='column',dim2name='levsoi',&
       long_name='total volume in macropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,VLMacroPores_vr,datrc_2d)      
  else
    !print*,'VLMacP'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,VLMacroPores_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='VLMacroPores_vr', dim1name='column',dim2name='levsoi',&
       long_name='total volume in macropores', units='m3 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='totRootLenDens_vr', dim1name='column',dim2name='levsoi',&
       long_name='total root length density', units='m m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,totRootLenDens_vr,datrc_2d)      
  else
    !print*,'totRootLenDens_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,totRootLenDens_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ)                    
    call restartvar(ncid, flag, varname='totRootLenDens_vr', dim1name='column',dim2name='levsoi',&
       long_name='total root length density', units='m m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers+1)                    
    call restartvar(ncid, flag, varname='CanopyHeightZ_col', dim1name='column',dim2name='levcan1',&
       long_name='canopy layer height', units='m m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeightZ_col,datrc_2d)      
  else
    !print*,'ZL'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyHeightZ_col,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers+1)                        
    call restartvar(ncid, flag, varname='CanopyHeightZ_col', dim1name='column',dim2name='levcan1',&
       long_name='canopy layer height', units='m', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers)                        
    call restartvar(ncid, flag, varname='CanopyLeafAareZ_col', dim1name='column',dim2name='levcan',&
       long_name='Grid total leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafAareZ_col,datrc_2d)      
  else
    !print*,'CanopyLeafAareZ_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyLeafAareZ_col,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers)    
    call restartvar(ncid, flag, varname='CanopyLeafAareZ_col', dim1name='column',dim2name='levcan',&
       long_name='Grid total leaf area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers+1)                        
    call restartvar(ncid, flag, varname='TAU_DirectSunLit', dim1name='column',dim2name='levcan1',&
       long_name='Direct radiation transmission coefficient', units='', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,TAU_DirectSunLit,datrc_2d)      
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TAU_DirectSunLit,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers+1)    
    call restartvar(ncid, flag, varname='TAU_DirectSunLit', dim1name='column',dim2name='levcan1',&
       long_name='Direct radiation transmission coefficient', units='', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers+1)                        
    call restartvar(ncid, flag, varname='TAU_DirectSunSha', dim1name='column',dim2name='levcan1',&
       long_name='Through canopy radiation coefficient', units='', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,TAU_DirectSunSha,datrc_2d)      
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TAU_DirectSunSha,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers+1)    
    call restartvar(ncid, flag, varname='TAU_DirectSunSha', dim1name='column',dim2name='levcan1',&
       long_name='Through canopy radiation coefficient', units='', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers)    
    call restartvar(ncid, flag, varname='CanopyStemAareZ_col', dim1name='column',dim2name='levcan',&
       long_name='total stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,CanopyStemAareZ_col,datrc_2d)      
  else
    !print*,'CanopyStemAareZ_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyStemAareZ_col,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers)    
    call restartvar(ncid, flag, varname='CanopyStemAareZ_col', dim1name='column',dim2name='levcan',&
       long_name='total stem area', units='m2 d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers)    
    call restartvar(ncid, flag, varname='tCanLeafC_clyr', dim1name='column',dim2name='levcan',&
       long_name='total leaf mass', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,tCanLeafC_clyr,datrc_2d)      
  else
    !print*,'tCanLeafC_clyr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,tCanLeafC_clyr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:NumCanopyLayers)    
    call restartvar(ncid, flag, varname='tCanLeafC_clyr', dim1name='column',dim2name='levcan',&
       long_name='total leaf mass', units='g d-2', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)    
    call restartvar(ncid, flag, varname='trcs_VLN_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='nutrient salt tracers', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,trcs_VLN_vr,datrc_3d)      
  else
    !print*,'trcs_VLN_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcs_VLN_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSolutTracers,1:JZ+1)    
    call restartvar(ncid, flag, varname='trcs_VLN_vr', dim1name='column',dim2name='soltrcs',&
       dim3name='levsoi1',long_name='nutrient salt tracers', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)         
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NPrecipTracers,1:JZ+1)    
    call restartvar(ncid, flag, varname='trcp_saltpml_vr', dim1name='column',dim2name='ptracers',&
       dim3name='levsoi1',long_name='salt precipitate in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,trcp_saltpml_vr,datrc_3d)      
  else
    !print*,'trcp_salml'
    if(flag=='write') call cpcol(flag,NHW,NHE,NVN,NVS,trcp_saltpml_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NPrecipTracers,1:JZ+1)        
    call restartvar(ncid, flag, varname='trcp_saltpml_vr', dim1name='column',dim2name='ptracers',&
       dim3name='levsoi1',long_name='salt precipitate in micropore', units='g d-2', &
       interpinic_flag='skip', data=datpr3, missing_value=spval, &
       fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)        
    call restartvar(ncid, flag, varname='SolutesIonConc_vr', dim1name='column',dim2name='levsoi',&
       long_name='soil solution ion concentration', units='mol m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
    call cpcol(flag,NHW,NHE,NVN,NVS,SolutesIonConc_vr,datrc_2d)      
  else
    !print*,'SolutesIonConc_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SolutesIonConc_vr,datrc_2d)        
    datpr2 => datrc_2d(1:ncols,1:JZ)           
    call restartvar(ncid, flag, varname='SolutesIonConc_vr', dim1name='column',dim2name='levsoi',&
       long_name='soil solution ion concentration', units='mol m-3', &
       interpinic_flag='skip', data=datpr2, missing_value=spval, &
       fill_value=spval)         
  endif  


  IF(salt_model)THEN

    ! Preserve prescribed boundary chemistry independently of irrigation fluxes.
    ! The shared salt dimension includes band entries; unused entries are zero.
    if(flag=='read')then
      call ncd_inqvid(ncid,'trcSalt_Irrig_vr',salt_boundary_varid,salt_boundary_vardesc, &
        readvar=has_salt_boundary_conc)
      if(.not.has_salt_boundary_conc)then
        write(iulog,*)'Salt restart lacks trcSalt_Irrig_vr boundary concentrations; regenerate the salt restart.'
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif
    endif
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ)
    if(flag=='write')then
      datpr3=0._r8
      call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_Irrig_vr,datrc_3d)
    endif
    call restartvar(ncid, flag, varname='trcSalt_Irrig_vr', dim1name='column',dim2name='satracers', &
      dim3name='levsoi',long_name='prescribed external-water salt concentration, unbanded',units='mol m-3', &
      interpinic_flag='skip',data=datpr3,missing_value=spval,fill_value=spval)
    if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_Irrig_vr,datrc_3d)

    if(flag=='read')then
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ+1)           
      call restartvar(ncid, flag, varname='trcSalt_solml_vr', dim1name='column',dim2name='satracers',&
        dim3name='levsoi1',long_name='soil aqueous salt content micropre', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)  
      call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_solml_vr,datrc_3d)      
    else
      !print*,'trcSalt_solml'
      if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_solml_vr,datrc_3d)        
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ+1)               
      call restartvar(ncid, flag, varname='trcSalt_solml_vr', dim1name='column',dim2name='satracers',&
        dim3name='levsoi1',long_name='soil aqueous salt content micropre', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)  
    endif  

    if(flag=='read')then
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JS)               
      call restartvar(ncid, flag, varname='trcs_solsml', dim1name='column',dim2name='satracers',&
         dim3name='levsno',long_name='snowpack salt dissolved tracers', units='mol d-2', &
         interpinic_flag='skip', data=datpr3, missing_value=spval, &
         fill_value=spval)  
        call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_ml_snvr,datrc_3d)      
    else
      !print*,'trcs_solsml'
      if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_ml_snvr,datrc_3d)        
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JS)    
      call restartvar(ncid, flag, varname='trcs_solsml', dim1name='column',dim2name='satracers',&
         dim3name='levsno',long_name='snowpack salt dissolved tracers', units='mol d-2', &
         interpinic_flag='skip', data=datpr3, missing_value=spval, &
         fill_value=spval)  

    endif  

    if(flag=='read')then
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ)    
      call restartvar(ncid, flag, varname='trcSalt_soHml_vr', dim1name='column',dim2name='satracers',&
        dim3name='levsoi',long_name='soil macropore aqueous salt dissolved tracers', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)  
      call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_soHml_vr,datrc_3d)      
    else
      !print*,'trcSalt_soHml'
      if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcSalt_soHml_vr,datrc_3d)          
      datpr3 => datrc_3d(1:ncols,1:trc_confs%NSaltTracers,1:JZ)    
      call restartvar(ncid, flag, varname='trcSalt_soHml_vr', dim1name='column',dim2name='satracers',&
        dim3name='levsoi',long_name='soil macropore aqueous salt dissolved tracers', units='mol d-2', &
        interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)
    endif  
  endif

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RO2DmndHetert_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial aqueous O2 demand', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RO2DmndHetert_vr,datrc_4d)      
  else
    !print*,'RO2DmndHetert'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RO2DmndHetert_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)        
    call restartvar(ncid, flag, varname='RO2DmndHetert_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial aqueous O2 demand', units='g d-2 h-1', &
      interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3ReduxDmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake non-band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3ReduxDmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RNO3ReduxDmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3ReduxDmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3ReduxDmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake non-band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndReduxSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndReduxSoilHeter_vr,datrc_4d)      
  else
    !print*,'RNO2DmndReduxSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndReduxSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndReduxSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RN2ODmndReduxHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total heterotrophic microbial N2O uptake unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RN2ODmndReduxHeter_vr,datrc_4d)      
  else
    !print*,'RN2ODmndReduxHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RN2ODmndReduxHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RN2ODmndReduxHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total heterotrophic microbial N2O uptake unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3ReduxDmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3ReduxDmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RNO3ReduxDmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3ReduxDmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3ReduxDmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO3 uptake band unconstrained by NO3', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndReduxBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndReduxBandHeter_vr,datrc_4d)      
  else
    !print*,'RNO2DmndReduxBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndReduxBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndReduxBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='total microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RNH4DmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RNO3DmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RH2PO4DmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial H1PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndSoilHeter_vr,datrc_4d)      
  else
    !print*,'RH2PO4DmndSoilHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndSoilHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4DmndSoilHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial H1PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RNH4DmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RNO3DmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='microbial NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H2PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RH2PO4DmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H2PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H1PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndBandHeter_vr,datrc_4d)      
  else
    !print*,'RH2PO4DmndBandHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndBandHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4DmndBandHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='substrate-unlimited H1PO4 mineraln-immobiln', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RDOCUptkHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial DOC flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
    call cpcol(flag,NHW,NHE,NVN,NVS,RDOCUptkHeter_vr,datrc_4d)      
  else
    !print*,'RDOCUptkHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RDOCUptkHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RDOCUptkHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial DOC flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)  
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RAcetateUptkHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial acetate flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RAcetateUptkHeter_vr,datrc_4d)      
  else
    !print*,'RAcetateUptkHeter_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RAcetateUptkHeter_vr,datrc_4d)        
    datpr4 => datrc_4d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RAcetateUptkHeter_vr', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='net microbial acetate flux', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RNH4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NH4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndLitrHeter_col,datrc_3d)      
  else
    !print*,'RNH4DmndLitrHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndLitrHeter_col,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RNH4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NH4 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif

  if(flag=='read')then
    call ncd_inqvid(ncid,'RNH4DmndLitrBandHeter_col',litter_band_varid,litter_band_vardesc, &
      readvar=has_litter_band_demand)
    if(has_litter_band_demand)then
      datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
      call restartvar(ncid, flag, varname='RNH4DmndLitrBandHeter_col', dim1name='column',dim2name='hetrmicb',&
        dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying band soil', &
        units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndLitrBandHeter_col,datrc_3d)
    else
      !Older restart files charged all supplemental uptake to nonband soil.
      RNH4DmndLitrBandHeter_col=0._r8
    endif
  else
    !print*,'RNH4DmndLitrBandHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4DmndLitrBandHeter_col,datrc_3d)
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
    call restartvar(ncid, flag, varname='RNH4DmndLitrBandHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RNO3DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NO3 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndLitrHeter_col,datrc_3d)      
  else
    !print*,'RNO3DmndLitrHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndLitrHeter_col,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RNO3DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial NO3 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif

  if(flag=='read')then
    call ncd_inqvid(ncid,'RNO3DmndLitrBandHeter_col',litter_band_varid,litter_band_vardesc, &
      readvar=has_litter_band_demand)
    if(has_litter_band_demand)then
      datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
      call restartvar(ncid, flag, varname='RNO3DmndLitrBandHeter_col', dim1name='column',dim2name='hetrmicb',&
        dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying band soil', &
        units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndLitrBandHeter_col,datrc_3d)
    else
      !Older restart files charged all supplemental uptake to nonband soil.
      RNO3DmndLitrBandHeter_col=0._r8
    endif
  else
    !print*,'RNO3DmndLitrBandHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3DmndLitrBandHeter_col,datrc_3d)
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
    call restartvar(ncid, flag, varname='RNO3DmndLitrBandHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RH2PO4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial PO4 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndLitrHeter_col,datrc_3d)      
  else
    !print*,'RH2PO4DmndLitrHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndLitrHeter_col,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)    
    call restartvar(ncid, flag, varname='RH2PO4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='microbial PO4 demand in surface litter', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif

  if(flag=='read')then
    call ncd_inqvid(ncid,'RH1PO4DmndLitrHeter_col',litter_band_varid,litter_band_vardesc, &
      readvar=has_litter_band_demand)
    if(has_litter_band_demand)then
      datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
      call restartvar(ncid, flag, varname='RH1PO4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
        dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying nonband soil', &
        units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndLitrHeter_col,datrc_3d)
    else
      !Older restart files charged all supplemental uptake to nonband soil.
      RH1PO4DmndLitrHeter_col=0._r8
    endif
  else
    !print*,'RH1PO4DmndLitrHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndLitrHeter_col,datrc_3d)
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
    call restartvar(ncid, flag, varname='RH1PO4DmndLitrHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying nonband soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)
  endif

  if(flag=='read')then
    call ncd_inqvid(ncid,'RH1PO4DmndLitrBandHeter_col',litter_band_varid,litter_band_vardesc, &
      readvar=has_litter_band_demand)
    if(has_litter_band_demand)then
      datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
      call restartvar(ncid, flag, varname='RH1PO4DmndLitrBandHeter_col', dim1name='column',dim2name='hetrmicb',&
        dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying band soil', &
        units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndLitrBandHeter_col,datrc_3d)
    else
      !Older restart files charged all supplemental uptake to nonband soil.
      RH1PO4DmndLitrBandHeter_col=0._r8
    endif
  else
    !print*,'RH1PO4DmndLitrBandHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4DmndLitrBandHeter_col,datrc_3d)
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
    call restartvar(ncid, flag, varname='RH1PO4DmndLitrBandHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)
  endif

  if(flag=='read')then
    call ncd_inqvid(ncid,'RH2PO4DmndLitrBandHeter_col',litter_band_varid,litter_band_vardesc, &
      readvar=has_litter_band_demand)
    if(has_litter_band_demand)then
      datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
      call restartvar(ncid, flag, varname='RH2PO4DmndLitrBandHeter_col', dim1name='column',dim2name='hetrmicb',&
        dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying band soil', &
        units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
        fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndLitrBandHeter_col,datrc_3d)
    else
      !Older restart files charged all supplemental uptake to nonband soil.
      RH2PO4DmndLitrBandHeter_col=0._r8
    endif
  else
    !print*,'RH2PO4DmndLitrBandHeter_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4DmndLitrBandHeter_col,datrc_3d)
    datpr3 => datrc_3d(1:ncols,1:NumHetetr1MicCmplx,1:jcplx)
    call restartvar(ncid, flag, varname='RH2PO4DmndLitrBandHeter_col', dim1name='column',dim2name='hetrmicb',&
      dim3name='nomcomplx',long_name='Litter microbial nutrient demand on underlying band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)
  endif  

  if(flag=='read')then
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='mBiomeHeter_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nlivehetermicb',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='heterotrophic microbial biomass element',units='g d-2', interpinic_flag='skip', &
      data=datpr5, missing_value=spval, fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,mBiomeHeter_vr,datrc_5d)      
  else
    !print*,'OMC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,mBiomeHeter_vr,datrc_5d)        
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='mBiomeHeter_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nlivehetermicb',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='heterotrophic microbial biomass element', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
  endif  

  
  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RO2MetaDmndAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='autotrophic aqueous O2 demand', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RO2MetaDmndAutor_vr,datrc_3d)      
  else
    !print*,'RO2MetaDmndAutor'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RO2MetaDmndAutor_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RO2MetaDmndAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='autotrophic aqueous O2 demand', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RCH4MetaDmndAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='autotrophic aqueous CH4 demand', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RCH4MetaDmndAutor_vr,datrc_3d)      
  else
    !print*,'RO2MetaDmndAutor'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RCH4MetaDmndAutor_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RCH4MetaDmndAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='autotrophic aqueous CH4 demand', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH3OxidAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH3OxidAutor_vr,datrc_3d)      
  else
    !print*,'RNH3OxidAutor'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH3OxidAutor_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH3OxidAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2XupAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2XupAutor_vr,datrc_3d)      
  else
    !print*,'RNO2XupAutor'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2XupAutor_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2XupAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3XupAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3XupAutor_vr,datrc_3d)      
  else
    !print*,'RNO2XupAutor'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3XupAutor_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3XupAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 uptake non-band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RN2ODmndReduxAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RN2ODmndReduxAutor_vr,datrc_3d)      
  else
    !print*,'RN2ODmndReduxAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RN2ODmndReduxAutor_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RN2ODmndReduxAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH3OxidAutorBand_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH3OxidAutorBand_vr,datrc_3d)      
  else
    !print*,'RNH3OxidAutorBand'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH3OxidAutorBand_vr,datrc_3d)        
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH3OxidAutorBand_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 uptake non-band unconstrained by NH4', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2XupAutorBand_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2XupAutorBand_vr,datrc_3d)     
  else
    !print*,'RNO2XupAutorBand'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2XupAutorBand_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2XupAutorBand_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO2 uptake band unconstrained by NO2', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)        
  endif

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkSoilAutor_vr,datrc_3d)     
  else
    !print*,'RNH4UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkSoilAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NH4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkSoilAutor_vr,datrc_3d)     
  else
    !print*,'RNO3UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkSoilAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial NO3 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkSoilAutor_vr,datrc_3d)     
  else
    !print*,'RH2PO4UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkSoilAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial H1PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkSoilAutor_vr,datrc_3d)     
  else
    !print*,'RH2PO4UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkSoilAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4UptkSoilAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial H1PO4 demand in soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial H1PO4 demand in band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkBandAutor_vr,datrc_3d)     
  else
    !print*,'RH2PO4UptkSoilAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkBandAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH1PO4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic microbial H1PO4 demand in band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)          
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkBandAutor_vr,datrc_3d)     
  else
    !print*,'RNH4UptkBandAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkBandAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNH4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NH4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkBandAutor_vr,datrc_3d)     
  else
    !print*,'RNO3UptkBandAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkBandAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO3UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic NO3 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic H2PO4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkBandAutor_vr,datrc_3d)     
  else
    !print*,'RH2PO4UptkBandAutor_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkBandAutor_vr,datrc_3d)       
    datpr3 => datrc_3d(1:ncols,1:NumMicrobAutoTrophCmplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='RH2PO4UptkBandAutor_vr', dim1name='column',dim2name='automicb',&
      dim3name='levsoi1',long_name='total autotrophic H2PO4 immobilization (+ve) - mineralization (-ve) band', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr3, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)    
    call restartvar(ncid, flag, varname='RNH4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NH4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkLitrAutor_col,datrc_2d)     
  else
    !print*,'RNH4UptkLitrAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkLitrAutor_col,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)       
    call restartvar(ncid, flag, varname='RNH4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NH4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)        
  endif

  if(flag=='read')then
    call ncd_inqvid(ncid,'RNH4UptkLitrBandAutor_col',litter_band_varid,litter_band_vardesc, &
      readvar=has_litter_band_demand)
    if(has_litter_band_demand)then
      datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)
      call restartvar(ncid, flag, varname='RNH4UptkLitrBandAutor_col', dim1name='column',dim2name='automicb',&
        long_name='Litter microbial nutrient demand on underlying band soil', &
        units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
        fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkLitrBandAutor_col,datrc_2d)
    else
      !Older restart files charged all supplemental uptake to nonband soil.
      RNH4UptkLitrBandAutor_col=0._r8
    endif
  else
    !print*,'RNH4UptkLitrBandAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNH4UptkLitrBandAutor_col,datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)
    call restartvar(ncid, flag, varname='RNH4UptkLitrBandAutor_col', dim1name='column',dim2name='automicb',&
      long_name='Litter microbial nutrient demand on underlying band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)       
    call restartvar(ncid, flag, varname='RNO3UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NO3 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkLitrAutor_col,datrc_2d)     
  else
    !print*,'RNO3UptkLitrAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkLitrAutor_col,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)       
    call restartvar(ncid, flag, varname='RNO3UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial NO3 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
  endif

  if(flag=='read')then
    call ncd_inqvid(ncid,'RNO3UptkLitrBandAutor_col',litter_band_varid,litter_band_vardesc, &
      readvar=has_litter_band_demand)
    if(has_litter_band_demand)then
      datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)
      call restartvar(ncid, flag, varname='RNO3UptkLitrBandAutor_col', dim1name='column',dim2name='automicb',&
        long_name='Litter microbial nutrient demand on underlying band soil', &
        units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
        fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkLitrBandAutor_col,datrc_2d)
    else
      !Older restart files charged all supplemental uptake to nonband soil.
      RNO3UptkLitrBandAutor_col=0._r8
    endif
  else
    !print*,'RNO3UptkLitrBandAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO3UptkLitrBandAutor_col,datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)
    call restartvar(ncid, flag, varname='RNO3UptkLitrBandAutor_col', dim1name='column',dim2name='automicb',&
      long_name='Litter microbial nutrient demand on underlying band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)       
    call restartvar(ncid, flag, varname='RH2PO4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial H2PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkLitrAutor_col,datrc_2d)     
  else
    !print*,'RH2PO4UptkLitrAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkLitrAutor_col,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)       
    call restartvar(ncid, flag, varname='RH2PO4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial H2PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
  endif

  if(flag=='read')then
    call ncd_inqvid(ncid,'RH2PO4UptkLitrBandAutor_col',litter_band_varid,litter_band_vardesc, &
      readvar=has_litter_band_demand)
    if(has_litter_band_demand)then
      datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)
      call restartvar(ncid, flag, varname='RH2PO4UptkLitrBandAutor_col', dim1name='column',dim2name='automicb',&
        long_name='Litter microbial nutrient demand on underlying band soil', &
        units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
        fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkLitrBandAutor_col,datrc_2d)
    else
      !Older restart files charged all supplemental uptake to nonband soil.
      RH2PO4UptkLitrBandAutor_col=0._r8
    endif
  else
    !print*,'RH2PO4UptkLitrBandAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH2PO4UptkLitrBandAutor_col,datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)
    call restartvar(ncid, flag, varname='RH2PO4UptkLitrBandAutor_col', dim1name='column',dim2name='automicb',&
      long_name='Litter microbial nutrient demand on underlying band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)       
    call restartvar(ncid, flag, varname='RH1PO4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial H1PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkLitrAutor_col,datrc_2d)     
  else
    !print*,'RH1PO4UptkLitrAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkLitrAutor_col,datrc_2d)       
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)       
    call restartvar(ncid, flag, varname='RH1PO4UptkLitrAutor_col', dim1name='column',dim2name='automicb',&
      long_name='total autotrophic microbial H1PO4 demand in surface litte', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)            
  endif

  if(flag=='read')then
    call ncd_inqvid(ncid,'RH1PO4UptkLitrBandAutor_col',litter_band_varid,litter_band_vardesc, &
      readvar=has_litter_band_demand)
    if(has_litter_band_demand)then
      datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)
      call restartvar(ncid, flag, varname='RH1PO4UptkLitrBandAutor_col', dim1name='column',dim2name='automicb',&
        long_name='Litter microbial nutrient demand on underlying band soil', &
        units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
        fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkLitrBandAutor_col,datrc_2d)
    else
      !Older restart files charged all supplemental uptake to nonband soil.
      RH1PO4UptkLitrBandAutor_col=0._r8
    endif
  else
    !print*,'RH1PO4UptkLitrBandAutor_col'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RH1PO4UptkLitrBandAutor_col,datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:NumMicrobAutoTrophCmplx)
    call restartvar(ncid, flag, varname='RH1PO4UptkLitrBandAutor_col', dim1name='column',dim2name='automicb',&
      long_name='Litter microbial nutrient demand on underlying band soil', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:NumPlantChemElms,1:NumLiveAutoBioms,1:JZ+1)       
    call restartvar(ncid, flag, varname='mBiomeAutor_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nliveautomicb',dim4name='levsoi1',long_name='autotrophic microbial biomass element', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)            
    call cpcol(flag,NHW,NHE,NVN,NVS,mBiomeAutor_vr,datrc_4d)     
  else
    !print*,'OMCff'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,mBiomeAutor_vr,datrc_4d)       
    datpr4 => datrc_4d(1:ncols,1:NumPlantChemElms,1:NumLiveAutoBioms,1:JZ+1)           
    call restartvar(ncid, flag, varname='mBiomeAutor_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nliveautomicb',dim4name='levsoi1',long_name='autotrophic microbial biomass element', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)            
  endif  

  if(flag=='read')then
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:nlbiomcp,1:jcplx,1:JZ+1)           
    call restartvar(ncid, flag, varname='OMBioResdu_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nlbiomcp',dim4name='nomcomplx',dim5name='levsoi1',long_name='microbial residue matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,OMBioResdu_vr,datrc_5d)     
  else
    !print*,'ORC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,OMBioResdu_vr,datrc_5d)       
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:nlbiomcp,1:jcplx,1:JZ+1)           
    call restartvar(ncid, flag, varname='OMBioResdu_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nlbiomcp',dim4name='nomcomplx',dim5name='levsoi1',long_name='microbial residue matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)               
    call restartvar(ncid, flag, varname='DOM_MicP_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic C micropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MicP_vr,datrc_4d)     
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MicP_vr,datrc_4d)       
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='DOM_MicP_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic matter micropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)               
    call restartvar(ncid, flag, varname='DOM_MicP_drib_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dribbling flux for dissolved organic C micropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MicP_drib_vr,datrc_4d)     
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MicP_drib_vr,datrc_4d)       
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='DOM_MicP_drib_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dribbling flux for dissolved organic matter micropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='DOM_MacP', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic matter macropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MacP_vr,datrc_4d)
  else
    !print*,'OQCH'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DOM_MacP_vr,datrc_4d)  
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='DOM_MacP', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='dissolved organic matter macropore', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='SorbedOM_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='adsorbed soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SorbedOM_vr,datrc_4d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SorbedOM_vr,datrc_4d)  
    datpr4 => datrc_4d(1:ncols,1:trc_confs%NDOMS,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='SorbedOM_vr', dim1name='column',dim2name='ndoms',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='adsorbed soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      

  endif  

  if(flag=='read')then
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:jsken,1:jcplx,1:JZ+1)                   
    call restartvar(ncid, flag, varname='SolidOM_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nkinecmp',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='solid soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SolidOM_vr,datrc_5d)
  else
    !print*,'OSC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SolidOM_vr,datrc_5d)  
    datpr5 => datrc_5d(1:ncols,1:NumPlantChemElms,1:jsken,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='SolidOM_vr', dim1name='column',dim2name='elmnts',&
      dim3name='nkinecmp',dim4name='nomcomplx',dim5name='levsoi1',&
      long_name='solid soil organic matter', &
      units='g d-2', interpinic_flag='skip', data=datpr5, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr4 => datrc_4d(1:ncols,1:jsken,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='SolidOMAct_vr', dim1name='column',dim2name='nkinecmp',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='colonized soil organic C', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,SolidOMAct_vr,datrc_4d)
  else
    !print*,'OSA'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SolidOMAct_vr,datrc_4d)  
    datpr4 => datrc_4d(1:ncols,1:jsken,1:jcplx,1:JZ+1)    
    call restartvar(ncid, flag, varname='SolidOMAct_vr', dim1name='column',dim2name='nkinecmp',&
      dim3name='nomcomplx',dim4name='levsoi1',long_name='colonized soil organic C', &
      units='g d-2', interpinic_flag='skip', data=datpr4, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndSoilChemo_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndSoilChemo_vr,datrc_2d)
  else
    !print*,'RNO2DmndSoilChemo_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndSoilChemo_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndSoilChemo_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake non-band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndBandChemo_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndBandChemo_vr,datrc_2d)
  else
    !print*,'RNO2DmndSoilChemo_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RNO2DmndBandChemo_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='RNO2DmndBandChemo_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total chemodenitrification N2O uptake band unconstrained by N2O', &
      units='g d-2 h-1', interpinic_flag='skip', data=datpr2, missing_value=spval, &
      fill_value=spval)      
  endif  
  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:NumPlantChemElms,1:JZ+1)    
    call restartvar(ncid, flag, varname='SoilOrgM_vr', dim1name='column',dim2name='elmnts',dim3name='levsoi1',&
      long_name='total soil organic matter', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      

    call cpcol(flag,NHW,NHE,NVN,NVS,SoilOrgM_vr,datrc_3d)

  else
    !print*,'ORGC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilOrgM_vr,datrc_3d)  
    datpr3 => datrc_3d(1:ncols,1:NumPlantChemElms,1:JZ+1)   
    call restartvar(ncid, flag, varname='SoilOrgM_vr', dim1name='column',dim2name='elmnts',dim3name='levsoi1',&
      long_name='total soil organic matter', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='OMLitrC_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total particulate organic C', units='g d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, &
      fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,OMLitrC_vr,datrc_2d)
  else
    !print*,'ORGR'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,OMLitrC_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='OMLitrC_vr', dim1name='column',dim2name='levsoi1',&
      long_name='total particulate organic C', units='g d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, &
      fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitro,1:JZ+1)    
    call restartvar(ncid, flag, varname='FertN_mole_soil_vr', dim1name='column',dim2name='fertN',&
      dim3name='levsoi1',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      
    call cpcol(flag,NHW,NHE,NVN,NVS,FertN_mole_soil_vr,datrc_3d)
  else
    !print*,'FertN_mole_soil_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertN_mole_soil_vr,datrc_3d)  
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitro,1:JZ+1)    
    call restartvar(ncid, flag, varname='FertN_mole_soil_vr', dim1name='column',dim2name='fertN', &
      dim3name='levsoi1',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)      
  endif  

  if(flag=='read')then
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitrob,1:JZ)    
    call restartvar(ncid, flag, varname='FertN_mole_Band_vr', dim1name='column',dim2name='fertNb',&
      dim3name='levsoi',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,FertN_mole_Band_vr,datrc_3d)
  else
    !print*,'FertN_mole_Band_vr'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertN_mole_Band_vr,datrc_3d)  
    datpr3 => datrc_3d(1:ncols,1:trc_confs%NFertNitrob,1:JZ)    
    call restartvar(ncid, flag, varname='FertN_mole_Band_vr', dim1name='column',dim2name='fertNb',&
      dim3name='levsoi',long_name='fertilizer application', units='g d-2', interpinic_flag='skip',&
      data=datpr3, missing_value=spval, fill_value=spval)   
  endif  

  ! Legacy restart files have no fertilizer-P fields; missing pools start at zero.
  if(flag=='read')then
    call ncd_inqvid(ncid,'FertP_mole_soil_vr',feedback_varid,feedback_vardesc,readvar=has_fertP_soil)
    FertP_mole_soil_vr=0._r8
    if(has_fertP_soil)then
      datpr2 => datrc_2d(1:ncols,1:JZ+1)
      call restartvar(ncid, flag, varname='FertP_mole_soil_vr', dim1name='column',dim2name='levsoi1',&
        long_name='Broadcast fertilizer P including surface litter', units='mol P d-2', interpinic_flag='skip',&
        data=datpr2, missing_value=spval, fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,FertP_mole_soil_vr,datrc_2d)
    endif
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertP_mole_soil_vr,datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:JZ+1)
    call restartvar(ncid, flag, varname='FertP_mole_soil_vr', dim1name='column',dim2name='levsoi1',&
      long_name='Broadcast fertilizer P including surface litter', units='mol P d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)
  endif

  ! Legacy restart files have no fertilizer-P fields; missing pools start at zero.
  if(flag=='read')then
    call ncd_inqvid(ncid,'FertP_mole_band_vr',feedback_varid,feedback_vardesc,readvar=has_fertP_band)
    FertP_mole_band_vr=0._r8
    if(has_fertP_band)then
      datpr2 => datrc_2d(1:ncols,1:JZ)
      call restartvar(ncid, flag, varname='FertP_mole_band_vr', dim1name='column',dim2name='levsoi',&
        long_name='Banded fertilizer P in soil', units='mol P d-2', interpinic_flag='skip',&
        data=datpr2, missing_value=spval, fill_value=spval)
      call cpcol(flag,NHW,NHE,NVN,NVS,FertP_mole_band_vr,datrc_2d)
    endif
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FertP_mole_band_vr,datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:JZ)
    call restartvar(ncid, flag, varname='FertP_mole_band_vr', dim1name='column',dim2name='levsoi',&
      long_name='Banded fertilizer P in soil', units='mol P d-2', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNHB', dim1name='column', &
      dim2name='levsoi',long_name='width of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthNH4_vr,datrc_2d)
  else
    !print*,'WDNHB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthNH4_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNHB', dim1name='column', &
      dim2name='levsoi',long_name='width of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='AEC_vr', dim1name='column', &
      dim2name='levsoi',long_name='anion exchange capacity', units='cmol kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,AEC_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,AEC_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='AEC_vr', dim1name='column', &
      dim2name='levsoi',long_name='anion exchange capacity', units='cmol kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CEC_vr', dim1name='column', &
      dim2name='levsoi',long_name='cation exchange capacity', units='cmol kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CEC_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CEC_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CEC_vr', dim1name='column', &
      dim2name='levsoi',long_name='cation exchange capacity', units='cmol kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CFE', dim1name='column', &
      dim2name='levsoi',long_name='soil Fe content', units='mg Fe kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CFE_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CFE_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CFE', dim1name='column', &
      dim2name='levsoi',long_name='soil Fe content', units='mg Fe kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CCA', dim1name='column', &
      dim2name='levsoi',long_name='soil Ca content', units='mg Ca kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CCA_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CCA_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CCA', dim1name='column', &
      dim2name='levsoi',long_name='soil Ca content', units='mg Ca kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CAL', dim1name='column', &
      dim2name='levsoi',long_name='soil Al content', units='mg Al kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CAL_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CAL_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CAL', dim1name='column', &
      dim2name='levsoi',long_name='soil Al content', units='mg Al kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CSO4', dim1name='column', &
      dim2name='levsoi',long_name='soil SO4 content', units='mg S kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CSO4_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CSO4_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CSO4', dim1name='column', &
      dim2name='levsoi',long_name='soil SO4 content', units='mg S kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CCL', dim1name='column', &
      dim2name='levsoi',long_name='soil Cl content', units='mg Cl kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CCL_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CCL_vr,datrc_2d)
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CCL', dim1name='column', &
      dim2name='levsoi',long_name='soil Cl content', units='mg Cl kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CALOH', dim1name='column', &
      dim2name='levsoi',long_name='soil AlOH content', units='mg Al kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,CALOH_vr,datrc_2d)
  else
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CALOH_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='CALOH', dim1name='column', &
      dim2name='levsoi',long_name='soil ALOH content', units='mg Al kg-1', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPNHB', dim1name='column', &
      dim2name='levsoi',long_name='depth of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessNH4_vr,datrc_2d)
  else
    !print*,'DPNHB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessNH4_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPNHB', dim1name='column', &
      dim2name='levsoi',long_name='depth of NH4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNOB', dim1name='column', &
      dim2name='levsoi',long_name='width of NO3 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
    call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthNO3_vr,datrc_2d)
  else
    !print*,'WDNOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthNO3_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDNOB', dim1name='column', &
      dim2name='levsoi',long_name='width of NO3 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPNOB', dim1name='column', &
      dim2name='levsoi',long_name='depth of NO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessNO3_vr,datrc_2d)
  else
    !print*,'DPNOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessNO3_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPNOB', dim1name='column', &
      dim2name='levsoi',long_name='depth of NO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDPOB', dim1name='column', &
      dim2name='levsoi',long_name='width of PO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthPO4_vr,datrc_2d)
  else
    !print*,'WDPOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandWidthPO4_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='WDPOB', dim1name='column', &
      dim2name='levsoi',long_name='width of PO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPPOB', dim1name='column', &
      dim2name='levsoi',long_name='depth of PO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessPO4_vr,datrc_2d)
  else
     !print*,'DPPOB'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BandThicknessPO4_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='DPPOB', dim1name='column', &
      dim2name='levsoi',long_name='depth of PO4 band', units='m', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNHUI', dim1name='column', &
      dim2name='levsoi1',long_name='current inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,ZNHUI_vr,datrc_2d)
  else
    !print*,'ZNHUI'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZNHUI_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNHUI', dim1name='column', &
      dim2name='levsoi1',long_name='current inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNHU0', dim1name='column', &
      dim2name='levsoi1',long_name='urea hydrolysis inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,ZNHU0_vr,datrc_2d)
  else
    !print*,'ZNHU0'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZNHU0_vr,datrc_2d)  
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNHU0', dim1name='column', &
      dim2name='levsoi1',long_name='urea hydrolysis inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNFNI', dim1name='column', &
      dim2name='levsoi1',long_name='current nitrification inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,ZNFNI_vr,datrc_2d)  
  else
    !print*,'ZNFNI'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZNFNI_vr,datrc_2d)    
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNFNI', dim1name='column', &
      dim2name='levsoi1',long_name='current nitrification inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNFN0', dim1name='column', &
      dim2name='levsoi1',long_name='initial nitrification inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,ZNFN0_vr,datrc_2d)  
  else
    !print*,'ZNFN0'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ZNFN0_vr,datrc_2d)    
    datpr2 => datrc_2d(1:ncols,1:JZ+1)    
    call restartvar(ncid, flag, varname='ZNFN0', dim1name='column', &
      dim2name='levsoi1',long_name='initial nitrification inhibition activity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDTPX', dim1name='column', &
      dim2name='month',long_name='accumulated change  for maximum temperature', &
      units='K', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDTPX,datrc_2d)    
  else
    !print*,'TDTPX'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDTPX,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDTPX', dim1name='column', &
      dim2name='month',long_name='accumulated change for maximum temperature', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDRAD', dim1name='column', &
      dim2name='month',long_name='accumulated change for radiation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDRAD,datrc_2d)    
  else
    !print*,'TDRAD'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDRAD,datrc_2d)      
    datpr2 => datrc_2d (1:ncols,1:12)   
    call restartvar(ncid, flag, varname='TDRAD', dim1name='column', &
      dim2name='month',long_name='accumulated change for radiation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDTPN', dim1name='column', &
      dim2name='month',long_name='accumulated change for minimum temperature', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDTPN,datrc_2d)    
  else
    !print*,'TDTPN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDTPN,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDTPN', dim1name='column', &
      dim2name='month',long_name='accumulated change for minimum temperature', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDWND', dim1name='column', &
      dim2name='month',long_name='accumulated change for wind speed', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDWND,datrc_2d)    
  else
    !print*,'TDWND'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDWND,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)   
    call restartvar(ncid, flag, varname='TDWND', dim1name='column', &
      dim2name='month',long_name='accumulated change for wind speed', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDHUM', dim1name='column', &
      dim2name='month',long_name='accumulated change for humidity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDHUM,datrc_2d)    
  else
    !print*,'TDHUM'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDHUM,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDHUM', dim1name='column', &
      dim2name='month',long_name='accumulated change for humidity', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDPRC', dim1name='column', &
      dim2name='month',long_name='accumulated change for precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
    call cpcol(flag,NHW,NHE,NVN,NVS,TDPRC,datrc_2d)    
  else
    !print*,'TDPRC'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDPRC,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDPRC', dim1name='column', &
      dim2name='month',long_name='accumulated change for precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)     
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDIRI', dim1name='column', &
      dim2name='month',long_name='accumulated change for irrigation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TDIRI,datrc_2d)    
  else
    !print*,'TDIRI'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDIRI,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDIRI', dim1name='column', &
      dim2name='month',long_name='accumulated change for irrigation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)     
    call restartvar(ncid, flag, varname='TDCN4', dim1name='column', &
      dim2name='month',long_name='accumulated change for NH4 in precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TDCN4,datrc_2d)    
  else
    !print*,'TDCN4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDCN4,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)     
    call restartvar(ncid, flag, varname='TDCN4', dim1name='column', &
      dim2name='month',long_name='accumulated change for NH4 in precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)   
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:12)     
    call restartvar(ncid, flag, varname='TDCNO', dim1name='column', &
      dim2name='month',long_name='accumulated change for NO3 in precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,TDCNO,datrc_2d)    
  else
    !print*,'TDCNO'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TDCNO,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:12)    
    call restartvar(ncid, flag, varname='TDCNO', dim1name='column', &
      dim2name='month',long_name='accumulated change for NO3 in precipitation', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKC4', dim1name='column', &
      dim2name='levsoi',long_name='Ca-NH4 Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKC4_vr,datrc_2d)    
  else
    !print*,'GKC4'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKC4_vr,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKC4', dim1name='column', &
      dim2name='levsoi',long_name='Ca-NH4 Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCH', dim1name='column', &
      dim2name='levsoi',long_name='Ca-H Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCH_vr,datrc_2d)    
  else
    !print*,'GKCH'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCH_vr,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCH', dim1name='column', &
      dim2name='levsoi',long_name='Ca-H Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCA', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Al Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCA_vr,datrc_2d)    
  else
    !print*,'GKCA'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCA_vr,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)    
    call restartvar(ncid, flag, varname='GKCA', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Al Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCM', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Mg Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCM_vr,datrc_2d)    
  else
    !print*,'GKCM'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCM_vr,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCM', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Mg Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCN', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Na Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCN_vr,datrc_2d)    
  else
    !print*,'GKCN'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCN_vr,datrc_2d)      
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCN', dim1name='column', &
      dim2name='levsoi',long_name='Ca-Na Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  if(flag=='read')then
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCK', dim1name='column', &
      dim2name='levsoi',long_name='Ca-K Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
    call cpcol(flag,NHW,NHE,NVN,NVS,GKCK_vr,datrc_2d)    
  else
    !print*,'GKCK'
    if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,GKCK_vr,datrc_2d)     
    datpr2 => datrc_2d(1:ncols,1:JZ)     
    call restartvar(ncid, flag, varname='GKCK', dim1name='column', &
      dim2name='levsoi',long_name='Ca-K Gapon selectivity coefficient', &
      units='none', interpinic_flag='skip',&
      data=datpr2, missing_value=spval, fill_value=spval)    
  endif  

  
  call destroy(datrc_3d)
  call destroy(datrc_4d)
  call destroy(datrc_5d)

  END SUBROUTINE restartnc_soil


  !-----------------------------------------------------------------------
  subroutine restFile_getfile( file, path )
    !
    ! !DESCRIPTION:
    ! Determine and obtain netcdf restart file
    !
    ! !USES:
    use EcoSIMConfig     , only : restartFileFullPath, nsrest, nsrStartup, nsrBranch 
    use EcoSIMConfig     , only : nsrContinue, case_name,  brnch_retain_casename
    use fileutil         , only : getfil
    !
    implicit none
    ! !ARGUMENTS:
    character(len=*), intent(out) :: file  ! name of netcdf restart file
    character(len=*), intent(out) :: path  ! full pathname of netcdf restart file
    !
    ! !LOCAL VARIABLES:
    integer :: status                      ! return status
    integer :: length                      ! temporary          
    character(len=256) :: ftest,ctest      ! temporaries
    !-----------------------------------------------------------------------

    ! Continue run:
    ! Restart file pathname is read restart pointer file 

    if (nsrest==nsrContinue) then
       call restFile_read_pfile( path )
       call getfil( path, file, 0 )
    
    ! Branch run: 
    ! Restart file pathname is obtained from namelist "restartFileFullPath"
    ! Check case name consistency (case name must be different for branch run, 
    ! unless namelist specification states otherwise)
    elseif (nsrest==nsrBranch) then
       length = len_trim(finidat)
       if (finidat(length-2:length) == '.nc') then
          path = trim(finidat) 
       else
          path = trim(finidat) // '.nc'
       end if
       call getfil( path, file, 0 )

       ! tcraig, adding xx. and .elm makes this more robust
       ctest = 'xx.'//trim(case_name)//'.ecosim'
       ftest = 'xx.'//trim(file)
       status = index(trim(ftest),trim(ctest))
       if (status /= 0 .and. .not.(brnch_retain_casename)) then
!          if (masterproc) then
             write(iulog,*) 'Must change case name on branch run if ',&
                  'brnch_retain_casename namelist is not set'
             write(iulog,*) 'previous case filename= ',trim(file),&
                  ' current case = ',trim(case_name), &
                  ' ctest = ',trim(ctest), &
                  ' ftest = ',trim(ftest)
!          end if
          call endrun(msg=errMsg(__FILE__, __LINE__)) 
       end if
    end if

  end subroutine restFile_getfile

!-----------------------------------------------------------------------
  subroutine restFile_open( flag, file, ncid )

  use ncdio_pio
  implicit none
  character(len=*),  intent(in) :: flag ! flag to specify read or write
  character(len=*),  intent(in) :: file ! filename
  type(file_desc_t), intent(out):: ncid ! netcdf id

    integer :: omode                              ! netCDF dummy variable
    character(len= 32) :: subname='restFile_open' ! subroutine name

    if (flag == 'write') then

       ! Create new netCDF file (in define mode) and set fill mode
       ! to "no fill" to optimize performance

!       if (masterproc) then	
          write(iulog,*)
          write(iulog,*)'restFile_open: writing restart dataset at ',&
               trim(file), ' at nstep = ',etimer%get_nstep()
          write(iulog,*)
!       end if+
       call ncd_pio_createfile(ncid, trim(file))

    else if (flag == 'read') then

       ! Open netcdf restart file

 !      if (masterproc) then
          write(iulog,*) 'Reading restart dataset'
 !      end if
       call ncd_pio_openfile (ncid, trim(file), 0)

    end if

  end subroutine restFile_open

  !-----------------------------------------------------------------------
  subroutine restFile_closeRestart( file )
    !
    ! !DESCRIPTION:
    ! Close restart file and write restart pointer file if
    ! in write mode, otherwise just close restart file if in read mode
    !
    ! !USES:
    implicit none
    !
    ! !ARGUMENTS:
    character(len=*) , intent(in) :: file  ! local output filename
    !
    ! !CALLED FROM:
    ! subroutine restart in this module
    !
    ! !REVISION HISTORY:
    ! Author: Mariana Vertenstein
    !
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                   !index
    !-----------------------------------------------------------------------

!    if (masterproc) then
       write(iulog,*) 'Successfully wrote local restart file ',trim(file)
       write(iulog,'(72a1)') ("-",i=1,60)
       write(iulog,*)
!    end if

  end subroutine restFile_closeRestart

  !-----------------------------------------------------------------------
  subroutine restFile_write_pfile( fnamer )
    !
    ! !DESCRIPTION:
    ! Open restart pointer file. Write names of current netcdf restart file.
    !
    ! !USES:
    use EcoSIMConfig, only : rpntdir, rpntfil
    use fileutil , only : relavu
    use fileutil , only : getavu, opnfil
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: fnamer
    !
    ! !LOCAL VARIABLES:
    integer :: m                    ! index
    integer :: nio                  ! restart pointer file
    character(len=256) :: filename  ! local file name
    character(len=datestrlen) :: curr_date
    !-----------------------------------------------------------------------

!    if (masterproc) then
       nio = getavu()
       filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
       call opnfil( filename, nio, 'f' )
       curr_date =etimer%get_calendar()
       write(nio,'(a)') fnamer
       call relavu( nio )
       write(iulog,*)'Successfully wrote local restart pointer file'
!    end if

  end subroutine restFile_write_pfile

  !-----------------------------------------------------------------------
  subroutine restFile_read_pfile( pnamer )
    !
    ! !DESCRIPTION:
    ! Setup restart file and perform necessary consistency checks
    !
    ! !USES:
    use fileutil    , only : opnfil, getavu, relavu
    use EcoSIMConfig, only : rpntfil, rpntdir, inst_suffix
    !
    ! !ARGUMENTS:
    character(len=*), intent(out) :: pnamer ! full path of restart file
    !
    ! !LOCAL VARIABLES:
    !EOP
    integer :: i                  ! indices
    integer :: nio                ! restart unit
    integer :: status             ! substring check status
    character(len=256) :: locfn   ! Restart pointer file name
    !-----------------------------------------------------------------------

    ! Obtain the restart file from the restart pointer file. 
    ! For restart runs, the restart pointer file contains the full pathname 
    ! of the restart file. For branch runs, the namelist variable 
    ! [restartFileFullPath] contains the full pathname of the restart file. 
    ! New history files are always created for branch runs.

!    if (masterproc) then
       write(iulog,*) 'Reading restart pointer file....'
!    endif

    nio = getavu()
    locfn = trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
    call opnfil (locfn, nio, 'f')
    read (nio,'(a256)') pnamer
    call relavu (nio)

!    if (masterproc) then
       write(iulog,*) 'Reading restart data.....'
       write(iulog,'(72a1)') ("-",i=1,60)
!    end if

  end subroutine restFile_read_pfile  

  !-----------------------------------------------------------------------
  subroutine restFile_read( bounds, file)

  implicit none
  ! !ARGUMENTS:
  character(len=*)               , intent(in)    :: file  ! output netcdf restart file
  type(bounds_type)              , intent(in)    :: bounds  

! !LOCAL VARIABLES:
  type(file_desc_t) :: ncid         ! netcdf id
  integer           :: nc
  integer :: i

  ! Open file
  print*,trim(file)
  call restFile_open( flag='read', file=file, ncid=ncid )
  call restartnc_soil_hydraulic_scheme(ncid,'read')
  print*,'restFile_open Successfully'
  ! Read file

  call restFile_dimcheck( ncid )

  print*,'hist_restart_ncd'
  call hist_restart_ncd (bounds, ncid, flag='read')

  print*,'Do error checking on file'
  
  call restFile_check_consistency(bounds, ncid)
  
  call timemgr_restart_io( ncid, flag='read')
  
  call restartnc(ncid,flag='read')
  
  ! Close file 
  call restFile_close( ncid )
!    if (masterproc) then
      write(iulog,'(72a1)') ("-",i=1,60)
      write(iulog,*) 'Successfully read restart data for restart run'
      write(iulog,*)
!    end if

  end subroutine restFile_read  

  !-----------------------------------------------------------------------
  subroutine restFile_check_consistency(bounds, ncid)
  implicit none
  ! !ARGUMENTS:
  type(bounds_type), intent(in)    :: bounds  ! bounds
  type(file_desc_t), intent(inout) :: ncid    ! netcdf id
  character(len=*), parameter :: subname=trim(mod_filename)//'::restFile_check_consistency'

  print*,subname
  end subroutine restFile_check_consistency

  !-----------------------------------------------------------------------
  subroutine restFile_dimcheck( ncid )
    !
    ! !DESCRIPTION:
  implicit none 
  ! !ARGUMENTS:
  type(file_desc_t), intent(inout) :: ncid

  if ( .not. is_cold_run() )then
!      call get_proc_global(ng=numg, nt=numt, nl=numl, nc=numc, np=nump, nCohorts=numCohort)
!      call check_dim(ncid, nameg, numg)
!      call check_dim(ncid, namet, numt)
!      call check_dim(ncid, namel, numl)
!      call check_dim(ncid, namec, numc)
!      call check_dim(ncid, namep, nump)
  end if

!  call check_dim(ncid, 'levsno'  , nlevsno)
!  call check_dim(ncid, 'levgrnd' , nlevgrnd)
!  call check_dim(ncid, 'levurb'  , nlevurb)
!  call check_dim(ncid, 'levlak'  , nlevlak) 

  call check_dim(ncid,'datestrlen',datestrlen)  
  call check_dim(ncid,'rootaxs',MaxNumRootAxes)
  call check_dim(ncid,'nodes',MaxNodesPerBranch)
  call check_dim(ncid,'nodes1',MaxNodesPerBranch+1)
  call check_dim(ncid,'cansecz',NumLeafInclinationClasses)
  call check_dim(ncid,'fertN',trc_confs%NFertNitro)
  call check_dim(ncid,'fertNb',trc_confs%NFertNitrob)
  call check_dim(ncid,'automicb',NumMicrobAutoTrophCmplx)
  call check_dim(ncid,'nlbiomcp',nlbiomcp)
  call check_dim(ncid, 'levsoi', JZ)
  call check_dim(ncid, 'levsoi1', JZ+1)
  call check_dim(ncid, 'levsno',  JS)
  call check_dim(ncid, 'levcan',NumCanopyLayers)
  call check_dim(ncid, 'levcan1',NumCanopyLayers+1)
  call check_dim(ncid, 'npfts',  JP)
  call check_dim(ncid, 'nbranches',MaxNumBranches)
  call check_dim(ncid, 'ngrstages',NumGrowthStages)
  call check_dim(ncid, 'elmnts',NumPlantChemElms)
  call check_dim(ncid, 'nkinecmp',jsken)
  call check_dim(ncid, 'nomcomplx',jcplx)
  call check_dim(ncid, 'sdim',3)   !grid dimension
  call check_dim(ncid,'hetrmicb',NumHetetr1MicCmplx)
  call check_dim(ncid,'rootyps'  , jroots)
  call check_dim(ncid,'xtracers' , trc_confs%nxtracers)

  if(salt_model)then
    call check_dim(ncid,'satracers', trc_confs%NSaltTracers)
  endif  
  call check_dim(ncid,'ptracers' , trc_confs%NPrecipTracers)
  call check_dim(ncid,'NnutrientTracers', trc_confs%NnutrientTracers)
  call check_dim(ncid,'gastrcs'  , trc_confs%NGasTracers)
  call check_dim(ncid,'soltrcs'  , trc_confs%NSolutTracers)

  call check_dim(ncid , 'string_length', 64        )
  call check_dim(ncid , 'month'   , 12        )
  end subroutine restFile_dimcheck
  !-----------------------------------------------------------------------
  character(len=256) function restFile_filename( rdate )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
  implicit none
  !
  ! !ARGUMENTS:
  character(len=*), intent(in) :: rdate   ! input date for restart file name 
  !-----------------------------------------------------------------------

  restFile_filename = "./"//trim(case_name)//".ecosim"//trim(inst_suffix)//&
        ".r."//trim(rdate)//".nc"
  call strip_null(restFile_filename)      
  call strip_space(restFile_filename)

!    if (masterproc) then
    write(iulog,*)'writing restart file ',trim(restFile_filename),' for model date = ',rdate
!    end if

  end function restFile_filename  

  !-----------------------------------------------------------------------
  subroutine restFile_write( bounds, file,  rdate, noptr)

  implicit none
  type(bounds_type), intent(in) :: bounds
  character(len=*)               , intent(in)    :: file             ! output netcdf restart file  
  character(len=*)               , intent(in), optional :: rdate
  logical                        , intent(in), optional :: noptr     ! if should NOT write to the restart pointer file  
   ! !LOCAL VARIABLES:
    type(file_desc_t) :: ncid ! netcdf id
    integer :: i       ! index
    logical :: ptrfile ! write out the restart pointer file
    !-----------------------------------------------------------------------

    if ( present(noptr) )then
       ptrfile = .not. noptr
    else
       ptrfile = .true.
    end if

    ! --------------------------------------------
    ! Open restart file
    ! --------------------------------------------

    call restFile_open( flag='write', file=file, ncid=ncid )

    ! --------------------------------------------
    ! Define dimensions and variables
    ! --------------------------------------------

    call restFile_dimset ( ncid )

    ! Define restart file variables

    call timemgr_restart_io(ncid, flag='define')
    
    if (present(rdate)) then       
      call hist_restart_ncd (bounds, ncid, flag='define', rdate=rdate )
    end if

    call restartnc(ncid,flag='define')

    call restFile_enddef( ncid )

    ! --------------------------------------------
    ! Write restart file variables
    ! --------------------------------------------
    
    call timemgr_restart_io( ncid, flag='write' )

    call hist_restart_ncd (bounds, ncid, flag='write' )

    call restartnc(ncid,flag='write')

    ! --------------------------------------------
    ! Close restart file and write restart pointer file
    ! --------------------------------------------
    
    call restFile_close( ncid )
    call restFile_closeRestart( file )

    ! Write restart pointer file
    
    if ( ptrfile ) call restFile_write_pfile( file )
    
    ! Write out diagnostic info

!    if (masterproc) then
       write(iulog,*) 'Successfully wrote out restart data at nstep = ',etimer%get_nstep()
       write(iulog,'(72a1)') ("-",i=1,60)
!    end if
    
  end subroutine restFile_write

  !=========================================================================================

  subroutine timemgr_restart_io( ncid, flag )

    !---------------------------------------------------------------------------------
    ! Read/Write information needed on restart to a netcdf file. 
    use ncdio_pio, only: ncd_int, var_desc_t, file_desc_t

    use EcosimConst, only : secspday
    !
    implicit none
    ! Arguments
    type(file_desc_t), intent(inout) :: ncid  ! netcdf id
    character(len=*), intent(in) :: flag  ! 'read' or 'write'
    !
    ! Local variables
    character(len=*), parameter :: sub = trim(mod_filename)//'::'//'timemgr_restart_io'
    integer :: rc                  ! return code
    logical :: readvar             ! determine if variable is on initial file
    character(len=datestrlen) :: curr_date   ! date of data in restart file, YYYYMMDDHHMMSS
    integer :: rst_caltype         ! calendar type
    integer :: isecspday
    integer :: varid
    integer, parameter :: noleap = 1
    integer, parameter :: gregorian = 2
    
    integer :: rst_start_ymd
    integer :: rst_ref_ymd
    integer :: rst_curr_ymd
    integer :: rst_start_tod
    integer :: rst_ref_tod
    integer :: rst_curr_tod
    integer :: rst_step_sec
    !---------------------------------------------------------------------------------

    if (flag == 'write') then
       rst_step_sec  = etimer%get_step_size()
       curr_date =etimer%get_calendar()
       read(start_date,'(I8)')rst_start_ymd       !the start time info of the simulation
       read(ref_date,'(I8)')rst_ref_ymd           !the referene time info of the simulation
       read(curr_date,'(I8)')rst_curr_ymd         !current date info 
       rst_start_tod=get_tod(start_date(9:14))
       rst_ref_tod  =get_tod(ref_date(9:14))
       rst_curr_tod =get_tod(curr_date(9:14))

       call ncd_putvar(ncid,'ref_date',ref_date)
       call ncd_putvar(ncid,'curr_date',curr_date)
       call ncd_putvar(ncid,'start_date',start_date)

    elseif(flag=='define')then

      call ncd_defvar(ncid,varname='ref_date',xtype=ncd_char, &
        dim1name='datestrlen', long_name='reference date of the simulation')
        
      call ncd_defvar(ncid,varname='start_date',xtype=ncd_char, &
        dim1name='datestrlen', long_name='starting date of the simulation')
        
      call ncd_defvar(ncid,varname='curr_date',xtype=ncd_char, &
        dim1name='datestrlen', long_name='current date of the simulation')

    elseif(flag=='read')then

    end if
    
    isecspday=secspday

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_step_sec', &
         long_name='seconds component of timestep size', units='sec',         &
         nvalid_range=(/0,isecspday/), fill_value=uninit_int,                &
         interpinic_flag='skip',  data=rst_step_sec)
    
    if ((flag == 'read') .and. ( rst_step_sec < 0 .or. rst_step_sec > isecspday )) then
       call endrun( sub//'ERROR: timemgr_rst_step_sec out of range',__LINE__)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_start_ymd', &
         long_name='start date', units='YYYYMMDD', fill_value=uninit_int,     &
         interpinic_flag='skip',  data=rst_start_ymd)

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_start_tod', &
         long_name='start time of day', units='sec',                           &
         nvalid_range=(/0,isecspday/), fill_value=uninit_int,                 &
         interpinic_flag='skip',  data=rst_start_tod)

    if ((flag == 'read') .and. ( rst_start_tod < 0 .or. rst_start_tod > isecspday )) then
       call endrun( sub//'ERROR: timemgr_rst_strart_tod out of range',__LINE__)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_ref_ymd',  &
         long_name='reference date', units='YYYYMMDD', fill_value=uninit_int, &
         interpinic_flag='skip',  data=rst_ref_ymd)

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_ref_tod',  &
         long_name='reference time of day', units='sec',                       &
         nvalid_range=(/0,isecspday/), fill_value=uninit_int,                 &
         interpinic_flag='skip',  data=rst_ref_tod)

    if ((flag == 'read') .and. ( rst_start_tod < 0 .or. rst_start_tod > isecspday )) then
       call endrun( sub//'ERROR: timemgr_rst_ref_tod out of range',__LINE__)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_curr_ymd', &
         long_name='current date', units='YYYYMMDD', fill_value=uninit_int,   &
         interpinic_flag='skip',  data=rst_curr_ymd)

    call restartvar(ncid=ncid, flag=flag, varname='timemgr_rst_curr_tod', &
         long_name='current time of day', units='sec',                         &
         nvalid_range=(/0,isecspday/), fill_value=uninit_int,                 &
         interpinic_flag='skip',  data=rst_curr_tod)

    if ((flag == 'read') .and. ( rst_curr_tod < 0 .or. rst_curr_tod > isecspday )) then
       call endrun( sub//'ERROR: timemgr_rst_ref_ymd out of range',__LINE__)
    end if

  end subroutine timemgr_restart_io
  !------------------------------------------------------------------------
  subroutine restFile_dimset( ncid )

  ! !ARGUMENTS:
  implicit none
  type(file_desc_t), intent(inout) :: ncid
  !
  ! !LOCAL VARIABLES:
  integer :: dimid               ! netCDF dimension id
  integer :: numg                ! total number of gridcells across all processors
  integer :: numt                ! total number of topounits across all processors
  integer :: numc                ! total number of columns across all processors
  integer :: nump                ! total number of pfts across all processors

  integer :: ier                 ! error status
  integer :: strlen_dimid        ! string dimension id
  character(len=  8) :: curdate  ! current date
  character(len=  8) :: curtime  ! current time
  character(len=256) :: str
  character(len= 32), parameter :: subname=trim(mod_filename)//'::restFile_dimset' ! subroutine name
  !------------------------------------------------------------------------

  numg=bounds%ngrid;numt=bounds%ntopou;numc=bounds%ncols;nump=bounds%npfts

    ! Define dimensions

  call ncd_defdim(ncid , nameg      , numg           ,  dimid)
  call ncd_defdim(ncid , namet      , numt           ,  dimid)
  call ncd_defdim(ncid , namec      , numc           ,  dimid)
  call ncd_defdim(ncid , namep      , nump           ,  dimid)
  
  ! "level" dimensions
  call ncd_defdim(ncid,'hour',24,dimid)
  call ncd_defdim(ncid,'datestrlen',datestrlen,dimid)  
  call ncd_defdim(ncid,'rootaxs',MaxNumRootAxes,dimid)
  call ncd_defdim(ncid,'nodes',MaxNodesPerBranch,dimid)
  call ncd_defdim(ncid,'nodes1',MaxNodesPerBranch+1,dimid)  
  call ncd_defdim(ncid,'cansecz',NumLeafInclinationClasses,dimid)
  call ncd_defdim(ncid,'fertN',trc_confs%NFertNitro, dimid)
  call ncd_defdim(ncid,'fertNb',trc_confs%NFertNitrob, dimid)
  call ncd_defdim(ncid,'automicb',NumMicrobAutoTrophCmplx, dimid)
  call ncd_defdim(ncid,'nlbiomcp',nlbiomcp,dimid)
  call ncd_defdim(ncid,'nlivehetermicb',NumLiveHeterBioms,dimid)
  call ncd_defdim(ncid,'nliveautomicb',NumLiveAutoBioms,dimid)
  call ncd_defdim(ncid, 'nonband_nutrients', size(trcn_irrig_vr,1), dimid)
  call ncd_defdim(ncid, 'moisture_class', 100, dimid)
  call ncd_defdim(ncid, 'levsoi', JZ, dimid)
  call ncd_defdim(ncid, 'levsoi1', JZ+1, dimid)
  call ncd_defdim(ncid, 'levsno',  JS,dimid)
  call ncd_defdim(ncid, 'levcan',NumCanopyLayers,dimid)
  call ncd_defdim(ncid, 'levcan1',NumCanopyLayers+1,dimid)
  call ncd_defdim(ncid, 'npfts',  JP,dimid)
  call ncd_defdim(ncid, 'nbranches',MaxNumBranches,dimid)
  call ncd_defdim(ncid, 'ngrstages',NumGrowthStages,dimid)
  call ncd_defdim(ncid, 'elmnts',NumPlantChemElms,dimid)
  call ncd_defdim(ncid, 'nkinecmp',jsken,dimid)
  call ncd_defdim(ncid, 'nomcomplx',jcplx,dimid)
  call ncd_defdim(ncid, 'nlitromcomplx',micpar%NumOfLitrCmplxs,dimid)
  call ncd_defdim(ncid, 'sdim',3,dimid)   !grid dimension
  call ncd_defdim(ncid,'hetrmicb',NumHetetr1MicCmplx,dimid)
  call ncd_defdim(ncid,'rootyps'  , jroots, dimid)
  call ncd_defdim(ncid,'xtracers' , trc_confs%nxtracers, dimid)
  call ncd_defdim(ncid,'ndoms',trc_confs%NDOMS,dimid)
  call ncd_defdim(ncid,'rootsegs',NMaxRootSegs,dimid)
  if(salt_model)call ncd_defdim(ncid,'satracers', trc_confs%NSaltTracers, dimid)
  call ncd_defdim(ncid,'ptracers' , trc_confs%NPrecipTracers, dimid)
  call ncd_defdim(ncid,'NnutrientTracers', trc_confs%NnutrientTracers, dimid)
  call ncd_defdim(ncid,'gastrcs'  , trc_confs%NGasTracers, dimid)
  call ncd_defdim(ncid,'soltrcs'  , trc_confs%NSolutTracers, dimid)
  call ncd_defdim(ncid , 'string_length', 64,  dimid)
  call ncd_defdim(ncid , 'month', 12,  dimid)


  if (do_budgets) then
  end if

    ! Define global attributes

!  call ncd_putatt(ncid, NCD_GLOBAL, 'Conventions', trim(conventions))
  call etimer%getdatetime(curdate, curtime)
  str = 'created on ' // curdate // ' ' // curtime
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'history' , trim(str)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'username', trim(username)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'host'    , trim(hostname)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'version' , trim(version)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'source'  , trim(source)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'case_title'     , trim(ctitle)),subname)
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'case_id'        , trim(case_name)),subname)
!  call ncd_putatt(ncid, NCD_GLOBAL, 'surface_dataset', trim(fsurdat))
!  call ncd_putatt(ncid, NCD_GLOBAL, 'flanduse_timeseries', trim(get_flanduse_timeseries()))
  call check_ret(ncd_putatt(ncid, NCD_GLOBAL, 'title', 'ECOSIM Restart information'),subname)
  call restartnc_soil_hydraulic_scheme(ncid,'define')

!  call restFile_add_ipft_metadata(ncid)
!  call restFile_add_icol_metadata(ncid)

  end subroutine restFile_dimset  
  !-----------------------------------------------------------------------
  subroutine restFile_close( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    !
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !
    ! !LOCAL VARIABLES:
    character(len=32) :: subname='restFile_close' ! subroutine name
    !-----------------------------------------------------------------------

    call ncd_pio_closefile(ncid)

  end subroutine restFile_close

!------------------------------------------------------------------------------------------

  integer function get_tod(hhmmss)
  !return the current time of day by seconds

  implicit none
  character(len=6) :: hhmmss
  integer :: hh,mm,ss

  read(hhmmss,'(I2,I2,I2)')hh,mm,ss
  get_tod=hh*3600+mm*60+ss

  end function get_tod

  !-----------------------------------------------------------------------
  subroutine restFile_enddef( ncid )
    !
    ! !DESCRIPTION:
    ! Read a CLM restart file.
    ! !ARGUMENTS:
    type(file_desc_t), intent(inout) :: ncid
    !-----------------------------------------------------------------------

    call ncd_enddef(ncid)

  end subroutine restFile_enddef

  !-----------------------------------------------------------------------
  subroutine get_restart_date(curr_date)
  !
  ! DESCRIPTION
  ! read the simulation date of restart file
  use EcoSIMConfig, only : rpntdir, rpntfil
  use fileutil , only : relavu
  use fileutil , only : getavu, opnfil
  implicit none
  character(len=datestrlen), intent(out) :: curr_date
  character(len=256) :: filename  ! local file name
  character(len=256) :: fnamer    ! restart file name
  integer :: nio
  integer :: ii,year,mon,day,tod,iostat,hh,mm,ss
  nio = getavu()
  filename= trim(rpntdir) //'/'// trim(rpntfil)//trim(inst_suffix)
  
  call opnfil( filename, nio, 'f' )  
  read(nio,'(A)',iostat=iostat) fnamer
  if (iostat == 0) then
    write(iulog,*) 'use restart file: ', fnamer
  else if (iostat > 0) then
    write(iulog,*) 'Error occurred during read operation: iostat = ', iostat
    call endrun(msg=errMsg(__FILE__, __LINE__)) 
  else
    write(iulog,*) 'End of file reached'
    call endrun(msg=errMsg(__FILE__, __LINE__)) 
  end if
  call relavu( nio )
  
  !get current date
  ii=2
  do
    ii=ii+1
    if(fnamer(ii-2:ii)=='.r.')exit
  enddo
  read(fnamer(ii+1:),'(I4)')year
  ii=ii+5
  read(fnamer(ii+1:),'(I2)')mon
  ii=ii+3
  read(fnamer(ii+1:),'(I2)')day
  ii=ii+3
  read(fnamer(ii+1:),'(I6)')tod
  ! Restart filenames use seconds since midnight; the driver matches HHMMSS.
  hh=tod/3600
  mm=mod(tod,3600)/60
  ss=mod(tod,60)
  write(curr_date,'(I4.4,5I2.2)')year,mon,day,hh,mm,ss
  end subroutine get_restart_date
end module restartMod
