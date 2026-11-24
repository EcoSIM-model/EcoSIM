module PlantSiteAPIData
  ! Description:
  ! Owns the site API type, instance and lifecycle procedures.
  use DebugToolMod, only: PrintInfo
  use PlantAPICommonData
  implicit none
  private
  public :: plt_site_Init, plt_site_destroy

  type, public :: plant_siteinfo_type
  integer  :: NY = 0
  integer  :: NX = 0
  real(r8) :: SoilSurfDepZ_col                     !soil surface depth, [m]
  real(r8) :: ALAT                                 !latitude,	[degrees north]
  real(r8) :: ATCA                                 !mean annual air temperature, [oC]
  real(r8) :: ALT                                  !altitude of the grid cell, [m]
  real(r8) :: CCO2EI_gperm3                        !initial atmospheric CO2 concentration, [g m-3]
  real(r8) :: CO2EI                                !initial atmospheric CO2 concentration, [umol mol-1]
  real(r8) :: COXYE                                !current atmospheric O2 concentration, [g m-3]
  real(r8) :: SolarNoonHour_col                    !time of solar noon, [h]
  real(r8) :: CO2E                                 !atmospheric CO2 concentration, [umol mol-1]
  real(r8) :: DayLenthPrev                         !daylength of previous day, [h]
  real(r8) :: DayLenthCurrent                      !current daylength of the grid, [h]
  real(r8) :: DayLenthMax_col                      !maximum daylength of the grid, [h]
  real(r8) :: SoilSurfRoughness_col              !initial soil surface roughness height, [m]
  real(r8) :: OXYE                                 !atmospheric O2 concentration, [umol mol-1]
  real(r8) :: PlantPopu_col                        !total plant population, [plants d-2]
  real(r8) :: POROS1                               !top layer soil porosity, [m3 m-3]
  real(r8) :: WindSpeedAtm_col                     !wind speed, [m h-1]
  real(r8) :: WindMesureHeight_col                 !wind speed measurement height, [m]
  real(r8) :: ZEROS2                               !threshold zero for numerical stability,[-]
  real(r8) :: ZEROS                                !threshold zero for numerical stability,[-]
  real(r8) :: ZERO                                 !threshold zero for numerical stability, [-]
  real(r8) :: ZERO2                                !threshold zero for numerical stability,[-]
  integer :: KoppenClimZone                        !Koppen climate zone for the grid,[-]
  integer :: NumActivePlants                       !number of active PFT in the grid, [-]
  integer :: NL                                    !lowest soil layer number,[-]
  integer :: NP0                                   !intitial number of plant species,[-]
  integer :: MaxNumRootLays                        !maximum root layer number,[-]
  integer :: NP                                    !current number of plant species,[-]
  integer :: NU                                    !current soil surface layer number, [-]
  integer :: NK                                    !current hydrologically active layer, [-]
  integer :: DazCurrYear                           !number of days in current year,[-]
  integer :: iYearCurrent                          !current year,[-]
  real(r8) :: QH2OLoss_lnds                        !total subsurface water loss flux over the landscape,	[m3 d-2]

  character(len=16), pointer :: DATAP(:)             => null()    !parameter file name,[-]
  CHARACTER(len=16), pointer :: DATA(:)              => null()    !pft file,[-]
  logical ,     pointer :: flag_active_pft(:)        => null()    !flag for active plants,  [-]
  real (r8), pointer :: PlantElemntStoreLandscape(:) => null()    !total plant element balance,                                              [g d-2]
  real (r8), pointer :: AtmGasc(:)                   => null()    !atmospheric gas concentrations,                                           [g m-3]
  real (r8), pointer :: AREA3(:)                     => null()    !soil cross section area (vertical plane defined by its normal direction), [m2]
  real (r8), pointer :: PlantElmBalCum_pft(:,:)       => null()    !cumulative plant element balance,                                         [g d-2]
  real (r8), pointer :: CumSoilThickness_vr(:)       => null()    !depth to bottom of soil layer from  surface of grid cell,                 [m]
  real (r8), pointer :: PPI_pft(:)                   => null()    !initial plant population,                                                 [plants d-2]
  real (r8), pointer :: PPatSeeding_pft(:)           => null()    !plant population at seeding,                                              [plants d-2]
  real (r8), pointer :: PPX_pft(:)                   => null()    !plant population,                                                         [plants m-2]
  real (r8), pointer :: PlantPopuLive_pft(:)         => null()    !plant population,                                         [d-2]
  real (r8), pointer :: PlantPopuDead_pft(:)         => null()     !live+standing dead plant population,                           [d-2]
  real (r8), pointer :: CumSoilThickMidL_vr(:)       => null()    !depth to middle of soil layer from  surface of grid cell, [m]
  real (r8), pointer :: FracSoiAsMicP_vr(:)          => null()    !micropore fraction,                                       [-]
  real (r8), pointer :: DLYR3(:)                     => null()    !vertical thickness of soil layer,                         [m]
  real (r8), pointer :: VLWatMicPM_vr(:,:)           => null()    !soil micropore water content,                             [m3 d-2]
  real (r8), pointer :: VLsoiAirPM_vr(:,:)           => null()    !soil air content,                                         [m3 d-2]
  real (r8), pointer :: TortMicPM_vr(:,:)            => null()    !micropore soil tortuosity,                                [m3 m-3]
  real (r8), pointer :: FILMM_vr(:,:)                => null()    !soil water film thickness,                                [m]
  real (r8), pointer :: SoilWeightStress_vr(:)         => null()    !soil bulk stress on root thickening, [MPa]
  real (r8), pointer :: SoilSuctStress_vr(:)         => null()    !soil suction stress on root thickening, [MPa]
  real (r8), pointer :: rSat_vr(:)                   => null()    !relative soil saturation, [-]
  contains
    procedure, public :: Init =>  plt_site_Init
    procedure, public :: Destroy => plt_site_destroy
  end type plant_siteinfo_type


  type(plant_siteinfo_type) , public, target :: plt_site      !site info

contains

  subroutine plt_site_Init(this)
  implicit none
  class(plant_siteinfo_type) :: this

  allocate(this%PlantElemntStoreLandscape(NumPlantChemElms));this%PlantElemntStoreLandscape=spval
  allocate(this%FracSoiAsMicP_vr(0:JZ1));this%FracSoiAsMicP_vr=spval
  allocate(this%AtmGasc(idg_beg:idg_NH3));this%AtmGasc=spval
  allocate(this%DATAP(JP1)); this%DATAP=''
  allocate(this%DATA(30)); this%DATA=''
  allocate(this%AREA3(0:JZ1));this%AREA3=spval
  allocate(this%DLYR3(0:JZ1)); this%DLYR3=spval
  allocate(this%PlantElmBalCum_pft(NumPlantChemElms,JP1));this%PlantElmBalCum_pft=spval
  allocate(this%CumSoilThickness_vr(0:JZ1));this%CumSoilThickness_vr=spval
  allocate(this%CumSoilThickMidL_vr(0:JZ1));this%CumSoilThickMidL_vr=spval
  allocate(this%PPI_pft(JP1));this%PPI_pft=spval
  allocate(this%PPatSeeding_pft(JP1));this%PPatSeeding_pft=spval
  allocate(this%PPX_pft(JP1));this%PPX_pft=spval
  allocate(this%PlantPopuLive_pft(JP1));this%PlantPopuLive_pft=spval
  allocate(this%PlantPopuDead_pft(JP1)); this%PlantPopuDead_pft=spval
  allocate(this%flag_active_pft(JP1));  this%flag_active_pft=.false.
  allocate(this%VLWatMicPM_vr(60,0:JZ1));this%VLWatMicPM_vr=spval
  allocate(this%VLsoiAirPM_vr(60,0:JZ1));this%VLsoiAirPM_vr=spval
  allocate(this%TortMicPM_vr(60,0:JZ1));this%TortMicPM_vr=spval
  allocate(this%FILMM_vr(60,0:JZ1)); this%FILMM_vr=spval
  allocate(this%SoilWeightStress_vr(JZ1));this%SoilWeightStress_vr=spval
  allocate(this%SoilSuctStress_vr(JZ1));this%SoilSuctStress_vr=spval
  allocate(this%rSat_vr(JZ1));this%rSat_vr=spval
  end subroutine plt_site_Init

  subroutine plt_site_destroy(this)
  ! Description:
  ! Release every array owned by this PlantSiteAPIData workspace.
  class(plant_siteinfo_type), intent(inout) :: this ! owned API state; component units and dimensions as declared above
  character(len=*), parameter :: subname='plt_site_destroy'

  call PrintInfo('beg '//subname)
  if(associated(this%DATAP))deallocate(this%DATAP);nullify(this%DATAP)
  if(associated(this%DATA))deallocate(this%DATA);nullify(this%DATA)
  if(associated(this%flag_active_pft))deallocate(this%flag_active_pft);nullify(this%flag_active_pft)
  if(associated(this%PlantElemntStoreLandscape))deallocate(this%PlantElemntStoreLandscape);nullify(this%PlantElemntStoreLandscape)
  if(associated(this%AtmGasc))deallocate(this%AtmGasc);nullify(this%AtmGasc)
  if(associated(this%AREA3))deallocate(this%AREA3);nullify(this%AREA3)
  if(associated(this%PlantElmBalCum_pft))deallocate(this%PlantElmBalCum_pft);nullify(this%PlantElmBalCum_pft)
  if(associated(this%CumSoilThickness_vr))deallocate(this%CumSoilThickness_vr);nullify(this%CumSoilThickness_vr)
  if(associated(this%PPI_pft))deallocate(this%PPI_pft);nullify(this%PPI_pft)
  if(associated(this%PPatSeeding_pft))deallocate(this%PPatSeeding_pft);nullify(this%PPatSeeding_pft)
  if(associated(this%PPX_pft))deallocate(this%PPX_pft);nullify(this%PPX_pft)
  if(associated(this%PlantPopuLive_pft))deallocate(this%PlantPopuLive_pft);nullify(this%PlantPopuLive_pft)
  if(associated(this%PlantPopuDead_pft))deallocate(this%PlantPopuDead_pft);nullify(this%PlantPopuDead_pft)
  if(associated(this%CumSoilThickMidL_vr))deallocate(this%CumSoilThickMidL_vr);nullify(this%CumSoilThickMidL_vr)
  if(associated(this%FracSoiAsMicP_vr))deallocate(this%FracSoiAsMicP_vr);nullify(this%FracSoiAsMicP_vr)
  if(associated(this%DLYR3))deallocate(this%DLYR3);nullify(this%DLYR3)
  if(associated(this%VLWatMicPM_vr))deallocate(this%VLWatMicPM_vr);nullify(this%VLWatMicPM_vr)
  if(associated(this%VLsoiAirPM_vr))deallocate(this%VLsoiAirPM_vr);nullify(this%VLsoiAirPM_vr)
  if(associated(this%TortMicPM_vr))deallocate(this%TortMicPM_vr);nullify(this%TortMicPM_vr)
  if(associated(this%FILMM_vr))deallocate(this%FILMM_vr);nullify(this%FILMM_vr)
  if(associated(this%SoilWeightStress_vr))deallocate(this%SoilWeightStress_vr);nullify(this%SoilWeightStress_vr)
  if(associated(this%SoilSuctStress_vr))deallocate(this%SoilSuctStress_vr);nullify(this%SoilSuctStress_vr)
  if(associated(this%rSat_vr))deallocate(this%rSat_vr);nullify(this%rSat_vr)
  call PrintInfo('end '//subname)
  end subroutine plt_site_destroy
end module PlantSiteAPIData
