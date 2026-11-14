module PlantPhaseWorkMod
  ! Description:
  ! Retain independent full-plant workspaces between hydraulic and BGC phases.
  ! Pointer-bearing API objects exchange ownership with a column workspace;
  ! assignments swap pointers rather than making shallow snapshot copies.
  ! The soil chemistry workspace is refreshed after soil/microbial updates.
  ! No phase workspace is needed at an hourly checkpoint boundary.
  use data_kind_mod, only: r8 => DAT_KIND_R8,yearIJ_type
  use PlantAPICommonData, only: JP1,JZ1,pltpar
  use abortutils, only: endrun
  use DebugToolMod, only: PrintInfo
  use PlantSiteAPIData, only: plant_siteinfo_type,plt_site
  use PlantPhotosynthesisAPIData, only: plant_photosyns_type,plt_photo
  use PlantRadiationAPIData, only: plant_radiation_type,plt_rad
  use PlantMorphologyAPIData, only: plant_morph_type,plt_morph
  use PlantPhenologyAPIData, only: plant_pheno_type,plt_pheno
  use PlantAllometryAPIData, only: plant_allometry_type,plt_allom
  use PlantBiomassAPIData, only: plant_biom_type,plt_biom
  use PlantEnergyWaterAPIData, only: plant_ew_type,plt_ew
  use PlantDisturbanceAPIData, only: plant_disturb_type,plt_distb
  use PlantBGCRatesAPIData, only: plant_bgcrate_type,plt_bgcr
  use PlantRootBGCAPIData, only: plant_rootbgc_type,plt_rbgc
  implicit none
  private
  public :: InitPlantPhaseWork,DestroyPlantPhaseWork,SwapPlantPhaseColumn
  public :: SavePlantHydraulicColumn,OpenPlantBGCColumn,ClosePlantBGCColumn
  public :: WaterLive_pft,WaterPassive_pft,WaterDryMatter_pft,WaterPath_pvr,WaterRadius_pvr
  public :: WaterRootShare_pvr,WaterRootDepth_pvr,WaterRootLength_pvr
  type :: PlantColumnWorkspace
    ! Description:
    ! Own one column's intermediate plant API state until its BGC phase finishes.
    logical :: Ready=.false.
    integer :: Year=0,Day=0,Hour=0
    type(plant_siteinfo_type) :: Site
    type(plant_photosyns_type) :: Photo
    type(plant_radiation_type) :: Radiation
    type(plant_morph_type) :: Morphology
    type(plant_pheno_type) :: Phenology
    type(plant_allometry_type) :: Allometry
    type(plant_biom_type) :: Biomass
    type(plant_ew_type) :: EnergyWater
    type(plant_disturb_type) :: Disturbance
    type(plant_bgcrate_type) :: BGCRates
    type(plant_rootbgc_type) :: RootBGC
  end type PlantColumnWorkspace
  type(PlantColumnWorkspace), allocatable, target :: Workspace(:,:)
  logical, allocatable :: WaterLive_pft(:,:,:),WaterPassive_pft(:,:,:) ! PFT/row/column [-]
  real(r8), allocatable :: WaterDryMatter_pft(:,:,:) ! accepted dry-matter fraction [-], PFT/row/column
  real(r8), allocatable :: WaterPath_pvr(:,:,:,:,:),WaterRadius_pvr(:,:,:,:,:) ! [m], root type/layer/PFT/row/column
  real(r8), allocatable :: WaterRootShare_pvr(:,:,:,:,:) ! root carbon share [-], root type/layer/PFT/row/column
  real(r8), allocatable :: WaterRootDepth_pvr(:,:,:,:) ! rooted layer fraction [-], layer/PFT/row/column
  real(r8), allocatable :: WaterRootLength_pvr(:,:,:,:,:) ! [m], root type/layer/PFT/row/column
contains
!------------------------------------------------------------------------------------------
  subroutine InitPlantPhaseWork(NHW,NHE,NVN,NVS)
    ! Description:
    ! Allocate column-owned plant and hydraulic work for the fixed active domain.
    integer, intent(in) :: NHW,NHE,NVN,NVS ! west/east/north/south column bounds [-]
    character(len=*), parameter :: subname='InitPlantPhaseWork'
    integer :: NY,NX

    call PrintInfo('beg '//subname)
    if(allocated(Workspace))then
      if(any(lbound(Workspace)/=[NVN,NHW]).or.any(ubound(Workspace)/=[NVS,NHE])) &
        call endrun('Plant phase workspace domain changed during the simulation',__LINE__)
    else
      allocate(Workspace(NVN:NVS,NHW:NHE))
      allocate(WaterLive_pft(JP1,NVN:NVS,NHW:NHE),WaterPassive_pft(JP1,NVN:NVS,NHW:NHE))
      allocate(WaterDryMatter_pft(JP1,NVN:NVS,NHW:NHE))
      allocate(WaterPath_pvr(pltpar%jroots,JZ1,JP1,NVN:NVS,NHW:NHE))
      allocate(WaterRadius_pvr(pltpar%jroots,JZ1,JP1,NVN:NVS,NHW:NHE))
      allocate(WaterRootShare_pvr(pltpar%jroots,JZ1,JP1,NVN:NVS,NHW:NHE))
      allocate(WaterRootDepth_pvr(JZ1,JP1,NVN:NVS,NHW:NHE))
      allocate(WaterRootLength_pvr(pltpar%jroots,JZ1,JP1,NVN:NVS,NHW:NHE))
      do NX=NHW,NHE
        do NY=NVN,NVS
          call InitPlantPhaseColumn(NY,NX)
        enddo
      enddo
    endif
    call PrintInfo('end '//subname)
  end subroutine InitPlantPhaseWork
!------------------------------------------------------------------------------------------
  subroutine InitPlantPhaseColumn(NY,NX)
    ! Description:
    ! Initialize an independent set of plant API arrays for one column.
    integer, intent(in) :: NY,NX ! column row and column indices [-]
    character(len=*), parameter :: subname='InitPlantPhaseColumn'
    call PrintInfo('beg '//subname)
    associate(                                         &
      Site           => Workspace(NY,NX)%Site        , & ! inout: owned Site API workspace; units/dimensions as in its type
      Photo          => Workspace(NY,NX)%Photo       , & ! inout: owned Photo API workspace; units/dimensions as in its type
      Radiation      => Workspace(NY,NX)%Radiation   , & ! inout: owned Radiation API workspace; units/dimensions as in its type
      Morphology     => Workspace(NY,NX)%Morphology  , & ! inout: owned Morphology API workspace; units/dimensions as in its type
      Phenology      => Workspace(NY,NX)%Phenology   , & ! inout: owned Phenology API workspace; units/dimensions as in its type
      Allometry      => Workspace(NY,NX)%Allometry   , & ! inout: owned Allometry API workspace; units/dimensions as in its type
      Biomass        => Workspace(NY,NX)%Biomass     , & ! inout: owned Biomass API workspace; units/dimensions as in its type
      EnergyWater    => Workspace(NY,NX)%EnergyWater , & ! inout: owned EnergyWater API workspace; units/dimensions as in its type
      Disturbance    => Workspace(NY,NX)%Disturbance , & ! inout: owned Disturbance API workspace; units/dimensions as in its type
      BGCRates       => Workspace(NY,NX)%BGCRates    , & ! inout: owned BGCRates API workspace; units/dimensions as in its type
      RootBGC        => Workspace(NY,NX)%RootBGC       & ! inout: owned RootBGC API workspace; units/dimensions as in its type
    )
    call Site%Init()
    call Photo%Init()
    call Radiation%Init()
    call Morphology%Init()
    call Phenology%Init()
    call Allometry%Init()
    call Biomass%Init()
    call EnergyWater%Init()
    call Disturbance%Init()
    call BGCRates%Init()
    call RootBGC%Init()
    end associate
    call PrintInfo('end '//subname)
  end subroutine InitPlantPhaseColumn
!------------------------------------------------------------------------------------------
  subroutine SwapPlantPhaseColumn(NY,NX)
    ! Description:
    ! Exchange column ownership of plant API pointers and scalar values. Each
    ! swap returns the previous shared workspace to the column, keeping targets
    ! distinct and retaining native API pointer association.
    integer, intent(in) :: NY,NX ! column row and column indices [-]
    character(len=*), parameter :: subname='SwapPlantPhaseColumn'
    type(plant_siteinfo_type) :: TempSite
    type(plant_photosyns_type) :: TempPhoto
    type(plant_radiation_type) :: TempRadiation
    type(plant_morph_type) :: TempMorphology
    type(plant_pheno_type) :: TempPhenology
    type(plant_allometry_type) :: TempAllometry
    type(plant_biom_type) :: TempBiomass
    type(plant_ew_type) :: TempEnergyWater
    type(plant_disturb_type) :: TempDisturbance
    type(plant_bgcrate_type) :: TempBGCRates
    type(plant_rootbgc_type) :: TempRootBGC
    call PrintInfo('beg '//subname)
    associate(                                      &
      Site         => Workspace(NY,NX)%Site       , & ! inout: owned Site API workspace; units/dimensions as in its type
      Photo        => Workspace(NY,NX)%Photo      , & ! inout: owned Photo API workspace; units/dimensions as in its type
      Radiation    => Workspace(NY,NX)%Radiation  , & ! inout: owned Radiation API workspace; units/dimensions as in its type
      Morphology   => Workspace(NY,NX)%Morphology , & ! inout: owned Morphology API workspace; units/dimensions as in its type
      Phenology    => Workspace(NY,NX)%Phenology  , & ! inout: owned Phenology API workspace; units/dimensions as in its type
      Allometry    => Workspace(NY,NX)%Allometry  , & ! inout: owned Allometry API workspace; units/dimensions as in its type
      Biomass      => Workspace(NY,NX)%Biomass    , & ! inout: owned Biomass API workspace; units/dimensions as in its type
      EnergyWater  => Workspace(NY,NX)%EnergyWater, & ! inout: owned EnergyWater API workspace; units/dimensions as in its type
      Disturbance  => Workspace(NY,NX)%Disturbance, & ! inout: owned Disturbance API workspace; units/dimensions as in its type
      BGCRates     => Workspace(NY,NX)%BGCRates   , & ! inout: owned BGCRates API workspace; units/dimensions as in its type
      RootBGC      => Workspace(NY,NX)%RootBGC      & ! inout: owned RootBGC API workspace; units/dimensions as in its type
    )
    TempSite        = plt_site
    plt_site        = Site
    Site            = TempSite
    TempPhoto       = plt_photo
    plt_photo       = Photo
    Photo           = TempPhoto
    TempRadiation   = plt_rad
    plt_rad         = Radiation
    Radiation       = TempRadiation
    TempMorphology  = plt_morph
    plt_morph       = Morphology
    Morphology      = TempMorphology
    TempPhenology   = plt_pheno
    plt_pheno       = Phenology
    Phenology       = TempPhenology
    TempAllometry   = plt_allom
    plt_allom       = Allometry
    Allometry       = TempAllometry
    TempBiomass     = plt_biom
    plt_biom        = Biomass
    Biomass         = TempBiomass
    TempEnergyWater = plt_ew
    plt_ew          = EnergyWater
    EnergyWater     = TempEnergyWater
    TempDisturbance = plt_distb
    plt_distb       = Disturbance
    Disturbance     = TempDisturbance
    TempBGCRates    = plt_bgcr
    plt_bgcr        = BGCRates
    BGCRates        = TempBGCRates
    TempRootBGC     = plt_rbgc
    plt_rbgc        = RootBGC
    RootBGC         = TempRootBGC
    end associate
    call PrintInfo('end '//subname)
  end subroutine SwapPlantPhaseColumn
!------------------------------------------------------------------------------------------
  subroutine SavePlantHydraulicColumn(yearIJ,NY,NX)
  !
  ! Description:
  ! Park accepted early plant state and record the within-hour handoff time.
  !
  type(yearIJ_type), intent(in) :: yearIJ ! current year/day/hour [-]
  integer, intent(in) :: NY,NX ! column row and column indices [-]
  character(len=*), parameter :: subname='SavePlantHydraulicColumn'

  call PrintInfo('beg '//subname)
  associate(                           &
    Ready   => Workspace(NY,NX)%Ready, & ! inout: hydraulic state awaits BGC [-]
    Year    => Workspace(NY,NX)%Year , & ! output: handoff year [-]
    Day     => Workspace(NY,NX)%Day  , & ! output: handoff day of year [-]
    Hour    => Workspace(NY,NX)%Hour   & ! output: handoff hour [-]
  )
  if(Ready)call endrun('Plant hydraulic phase repeated before BGC',__LINE__)
  call SwapPlantPhaseColumn(NY,NX)
  Ready=.true.
  Year=yearIJ%year;Day=yearIJ%I;Hour=yearIJ%J
  end associate
  call PrintInfo('end '//subname)
  end subroutine SavePlantHydraulicColumn
!------------------------------------------------------------------------------------------
  subroutine OpenPlantBGCColumn(yearIJ,NY,NX)
  !
  ! Description:
  ! Recover exactly this column's accepted hydraulic and prepared plant state.
  !
  type(yearIJ_type), intent(in) :: yearIJ ! current year/day/hour [-]
  integer, intent(in) :: NY,NX ! column row and column indices [-]
  character(len=*), parameter :: subname='OpenPlantBGCColumn'
  call PrintInfo('beg '//subname)
  if(.not.allocated(Workspace))call endrun('Plant BGC has no hydraulic workspace',__LINE__)
  associate( &
    Ready => Workspace(NY,NX)%Ready, & ! input: hydraulic phase completed [-]
    Year  => Workspace(NY,NX)%Year , & ! input: saved handoff year [-]
    Day   => Workspace(NY,NX)%Day  , & ! input: saved handoff day of year [-]
    Hour  => Workspace(NY,NX)%Hour   & ! input: saved handoff hour [-]
  )

  if(.not.Ready.or.Year/=yearIJ%year.or.Day/=yearIJ%I.or.Hour/=yearIJ%J) &
    call endrun('Plant BGC hydraulic handoff has the wrong column or time',__LINE__)
  call SwapPlantPhaseColumn(NY,NX)
  end associate
  call PrintInfo('end '//subname)
  end subroutine OpenPlantBGCColumn
!------------------------------------------------------------------------------------------
  subroutine ClosePlantBGCColumn(NY,NX)
  !
  ! Description:
  ! Return the shared API workspace and mark this hourly handoff consumed.
  !
  integer, intent(in) :: NY,NX ! column row and column indices [-]
  character(len=*), parameter :: subname='ClosePlantBGCColumn'

  call PrintInfo('beg '//subname)
  associate(Ready => Workspace(NY,NX)%Ready) ! output: no intermediate state remains [-]

  call SwapPlantPhaseColumn(NY,NX)
  Ready=.false.
  end associate
  call PrintInfo('end '//subname)
  end subroutine ClosePlantBGCColumn
!------------------------------------------------------------------------------------------
  subroutine DestroyPlantPhaseColumn(NY,NX)
    ! Description:
    ! Release the column-owned API arrays after all hydraulic handoffs finish.
    integer, intent(in) :: NY,NX ! column row and column indices [-]
    character(len=*), parameter :: subname='DestroyPlantPhaseColumn'
    call PrintInfo('beg '//subname)
    associate(                                     &
      Site        => Workspace(NY,NX)%Site        , & ! inout: owned Site API workspace; units/dimensions as in its type
      Photo       => Workspace(NY,NX)%Photo       , & ! inout: owned Photo API workspace; units/dimensions as in its type
      Radiation   => Workspace(NY,NX)%Radiation   , & ! inout: owned Radiation API workspace; units/dimensions as in its type
      Morphology  => Workspace(NY,NX)%Morphology  , & ! inout: owned Morphology API workspace; units/dimensions as in its type
      Phenology   => Workspace(NY,NX)%Phenology   , & ! inout: owned Phenology API workspace; units/dimensions as in its type
      Allometry   => Workspace(NY,NX)%Allometry   , & ! inout: owned Allometry API workspace; units/dimensions as in its type
      Biomass     => Workspace(NY,NX)%Biomass     , & ! inout: owned Biomass API workspace; units/dimensions as in its type
      EnergyWater => Workspace(NY,NX)%EnergyWater , & ! inout: owned EnergyWater API workspace; units/dimensions as in its type
      Disturbance => Workspace(NY,NX)%Disturbance , & ! inout: owned Disturbance API workspace; units/dimensions as in its type
      BGCRates    => Workspace(NY,NX)%BGCRates    , & ! inout: owned BGCRates API workspace; units/dimensions as in its type
      RootBGC     => Workspace(NY,NX)%RootBGC       & ! inout: owned RootBGC API workspace; units/dimensions as in its type
    )
  call Site%Destroy()
  call Photo%Destroy()
  call Radiation%Destroy()
  call Morphology%Destroy()
  call Phenology%Destroy()
  call Allometry%Destroy()
  call Biomass%Destroy()
  call EnergyWater%Destroy()
  call Disturbance%Destroy()
  call BGCRates%Destroy()
  call RootBGC%Destroy()
  end associate
  call PrintInfo('end '//subname)
  end subroutine DestroyPlantPhaseColumn
!------------------------------------------------------------------------------------------
  subroutine DestroyPlantPhaseWork
  !
  ! Description:
  ! Release within-hour plant work independently of native persistent state.
  character(len=*), parameter :: subname='DestroyPlantPhaseWork'
  integer :: NY,NX
  call PrintInfo('beg '//subname)
  if(allocated(Workspace))then
    if(any(Workspace%Ready))call endrun('Unconsumed plant hydraulic handoff at cleanup',__LINE__)
    do NX=lbound(Workspace,2),ubound(Workspace,2)
      do NY=lbound(Workspace,1),ubound(Workspace,1)
        call DestroyPlantPhaseColumn(NY,NX)
      enddo
    enddo
    deallocate(Workspace,WaterLive_pft,WaterPassive_pft,WaterDryMatter_pft,WaterPath_pvr,WaterRadius_pvr)
    deallocate(WaterRootShare_pvr,WaterRootDepth_pvr,WaterRootLength_pvr)
  endif
  call PrintInfo('end '//subname)
  end subroutine DestroyPlantPhaseWork
end module PlantPhaseWorkMod
