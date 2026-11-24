module PlantPhotosynthesisAPIData
  ! Description:
  ! Owns the photosynthesis API type, instance and lifecycle procedures.
  use DebugToolMod, only: PrintInfo
  use PlantAPICommonData
  implicit none
  private
  public :: plt_photo_init, plt_photo_destroy

  type, public :: plant_photosyns_type
  integer,  pointer :: iPlantPhotosynsType_pft(:)          => null()  !plant photosynthetic type (C3 or C4),[-]
  real(r8), pointer :: SpecLeafChlAct_pft(:)                => null()  !cholorophyll activity  at 25 oC,                                           [umol g-1 h-1]
  real(r8), pointer :: LeafProtein2Chl_pft(:)               => null()  !fraction of leaf protein that is chlorophyll-binded, [gC gC-1]
  real(r8), pointer :: LeafPEP2Protein_pft(:)               => null()  !leaf PEP carboxylase content,                                              [gC gC-1]
  real(r8), pointer :: fMesophyllChlProtein_pft(:)          => null()  !fraction of Chl-bound protein in mesophyll cell, [gC gC-1]
  real(r8), pointer :: LeafRubisco2Protein_pft(:)           => null()  !leaf rubisco content,                                                      [gC gC-1]
  real(r8), pointer :: VmaxPEPCarboxyRef_pft(:)             => null()  !PEP carboxylase activity at 25 oC                                          [umol g-1 h-1]
  real(r8), pointer :: VmaxRubOxyRef_pft(:)                 => null()  !rubisco oxygenase activity  at 25 oC,                                      [umol g-1 h-1]
  real(r8), pointer :: VmaxSpecRubCarboxyRef_pft(:)             => null()  !rubisco carboxylase activity  at 25 oC,                                    [umol g-1 h-1]
  real(r8), pointer :: XKCO2_pft(:)                         => null()  !Km for rubisco carboxylase activity,                                       [uM]
  real(r8), pointer :: XKO2_pft(:)                          => null()  !Km for rubisco oxygenase activity,                                         [uM]
  real(r8), pointer :: RubiscoActivity_brch(:,:)            => null()   !branch down-regulation of CO2 fixation,                                   [-]
  real(r8), pointer :: GrainFillDowreg_brch(:,:)           => null()   !down-regulation of C4 photosynthesis,                                     [-]
  real(r8), pointer :: aquCO2Intraleaf_pft(:)               => null()   !leaf aqueous CO2 concentration,                                           [uM]
  real(r8), pointer :: O2I_pft(:)                           => null()   !leaf gaseous O2 concentration,                                            [umol m-3]
  real(r8), pointer :: LeafIntracellularCO2_pft(:)          => null()   !leaf gaseous CO2 concentration,                                           [umol m-3]
  real(r8), pointer :: Km4RubiscoCarboxy_pft(:)             => null()   !leaf aqueous CO2 Km ambient O2,                                           [uM]
  real(r8), pointer :: Km4LeafaqCO2_pft(:)                  => null()   !leaf aqueous CO2 Km no O2,                                                [uM]
  real(r8), pointer :: RCS_pft(:)                           => null()   !shape parameter for calculating stomatal resistance from turgor pressure, [-]
  real(r8), pointer :: CanopyCi2CaRatio_pft(:)                    => null()   !Ci:Ca ratio,                                                              [-]
  real(r8), pointer :: H2OCuticleResist_pft(:)              => null()   !maximum stomatal resistance to vapor,                                     [s h-1]
  real(r8), pointer :: ChillHours_pft(:)                    => null()   !chilling effect on CO2 fixation,                                          [-]
  real(r8), pointer :: CO2Solubility_pft(:)                 => null()   !leaf CO2 solubility,                                                      [uM /umol mol-1]
  real(r8), pointer :: CanopyGasCO2_pft(:)                  => null()   !canopy gaesous CO2 concentration,                                         [umol mol-1]
  real(r8), pointer :: O2L_pft(:)                           => null()   !leaf aqueous O2 concentration,                                            [uM]
  real(r8), pointer :: Km4PEPCarboxy_pft(:)                 => null()   !Km for PEP carboxylase activity,                                          [uM]
  real(r8), pointer :: CanopyMinStomaResistH2O_pft(:)         => null()   !canopy minimum stomatal resistance,                          [s m-1]
  real(r8), pointer :: CuticleResist_pft(:)                 => null()   !maximum stomatal resistance to vapor,                        [s m-1]
  real(r8), pointer :: CanPStomaResistH2O_pft(:)            => null()   !canopy stomatal resistance,                                  [h m-1]
  real(r8), pointer :: RawH2OCanopy2Atm_pft(:)              => null()   !canopy boundary layer resistance,                            [h m-1]
  real(r8), pointer :: RawHeatCanopy2Atm_pft(:)              => null()   !canopy sensible-heat boundary layer resistance,                            [h m-1]
  real(r8), pointer :: LeafO2Solubility_pft(:)              => null()   !leaf O2 solubility,                                          [uM /umol mol-1]
  real(r8), pointer :: CO2CuticleResist_pft(:)              => null()   !maximum stomatal resistance to CO2,                          [s h-1]
  real(r8), pointer :: DiffCO2Atmos2Intracel_pft(:)         => null()   !gaesous CO2 concentration difference across stomates,        [umol m-3]
  real(r8), pointer :: LeafEffArea_zsec(:,:,:,:,:)        => null()  !leaf irradiated surface area in different leaf sector,       [m2 d-2]
  real(r8), pointer :: CPOOL3_node(:,:,:)                   => null()   !minimum sink strength for nonstructural C transfer,          [g d-2]
  real(r8), pointer :: CPOOL4_node(:,:,:)                   => null()   !leaf nonstructural C4 content in C4 photosynthesis,          [g d-2]
  real(r8), pointer :: CMassCO2BundleSheath_node(:,:,:)     => null()   !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: CO2CompenPoint_node(:,:,:)           => null()   !CO2 compensation point,                                      [uM]
  real(r8), pointer :: VoMaxRubiscoRef_brch(:,:)            => null()   !branch maximum rubisco oxygenation rate at reference temperature, [umol g-1 h-1]
  real(r8), pointer :: VcMaxRubiscoRef_brch(:,:)            => null()   !branch maximum rubisco carboxylation rate at reference temperature, [umol g-1 h-1]
  real(r8), pointer :: RubiscoCarboxyEff_node(:,:,:)        => null()   !carboxylation efficiency,                                    [umol umol-1]
  real(r8), pointer :: VcMaxPEPCarboxyRef_brch(:,:)         => null()   !branch reference maximum dark C4 carboxylation rate under saturating CO2, [umol s-1]
  real(r8), pointer :: ElectronTransptJmaxRef_brch(:,:)     => null()   !branch Jmax at reference temperature, [umol e- s-1]
  real(r8), pointer :: C4CarboxyEff_node(:,:,:)             => null()   !C4 carboxylation efficiency,                                 [umol umol-1]
  real(r8), pointer :: LigthSatCarboxyRate_node(:,:,:)      => null()   !maximum light carboxylation rate under saturating CO2,       [umol m-2 s-1]
  real(r8), pointer :: LigthSatC4CarboxyRate_node(:,:,:)    => null()   !maximum  light C4 carboxylation rate under saturating CO2,   [umol m-2 s-1]
  real(r8), pointer :: NutrientCtrlonC4Carboxy_node(:,:,:)  => null()   !down-regulation of C4 photosynthesis,                        [-]
  real(r8), pointer :: CMassHCO3BundleSheath_node(:,:,:)    => null()   !bundle sheath nonstructural C3 content in C4 photosynthesis, [g d-2]
  real(r8), pointer :: Vmax4RubiscoCarboxy_node(:,:,:)      => null()   !maximum dark carboxylation rate under saturating CO2,        [umol m-2 s-1]
  real(r8), pointer :: ProteinCperm2LeafArea_node(:,:,:)    => null()   !Protein C per m2 of leaf aera,                               [gC (leaf area m-2)]
  real(r8), pointer :: CO2lmtRubiscoCarboxyRate_node(:,:,:) => null()   !carboxylation rate,                                          [umol m-2 s-1]
  real(r8), pointer :: Vmax4PEPCarboxy_node(:,:,:)          => null()   !maximum dark C4 carboxylation rate under saturating CO2,     [umol m-2 s-1]
  real(r8), pointer :: CO2lmtPEPCarboxyRate_node(:,:,:)     => null()   !C4 carboxylation rate,                                       [umol m-2 s-1]
  real(r8), pointer :: AirConc_pft(:)                       => null()   !total gas concentration,                                     [mol m-3]
  real(r8), pointer :: CanopyVcMaxRubisco25C_pft(:)         => null()   !Canopy VcMax for rubisco carboxylation, [umol h-1 m-2]
  real(r8), pointer :: CanopyVoMaxRubisco25C_pft(:)         => null()   !Canopy VoMax for rubisco oxygenation, [umol h-1 m-2]
  real(r8), pointer :: CanopyVcMaxPEP25C_pft(:)             => null()   !Canopy VcMax in PEP C4 fixation, [umol h-1 m-2]
  real(r8), pointer :: ElectronTransptJmax25C_pft(:)        => null()   !Canopy Jmax at reference temperature, [umol e- s-1 m-2]
  real(r8), pointer :: TFN_Carboxy_pft(:)                   => null()   !temperature dependence of carboxylation, [-]
  real(r8), pointer :: TFN_Oxygen_pft(:)                    => null()   !temperature dependence of oxygenation, [-]
  real(r8), pointer :: TFN_eTranspt_pft(:)                  => null()   !temperature dependence of electron transport, [-]
  real(r8), pointer :: LeafAreaSunlit_pft(:)                => null()   !leaf irradiated surface area, [m2 d-2]
  real(r8), pointer :: PARSunlit_pft(:)                     => null()   !PAR absorbed by sunlit leaf, [umol m-2 s-1]
  real(r8), pointer :: PARSunsha_pft(:)                     => null()   !PAR absorbed by sun-shaded leaf, [umol m-2 s-1]
  real(r8), pointer :: DynCi2CaRatio_pft(:) => null() !dynamic intracellular-to-canopy CO2 ratio, [-]
  real(r8), pointer :: CO2Intra_pft(:) => null() !leaf-area-weighted intracellular CO2 sum, [umol mol-1 m2 d-2]
  real(r8), pointer :: CO2IntraScal_pft(:) => null() !intracellular CO2 leaf-area weight, [m2 d-2]
  real(r8), pointer :: CH2OSunlit_pft(:)                    => null()   !carbon fixation by sun-lit leaf, [gC d-2 h-1]
  real(r8), pointer :: CH2OSunsha_pft(:)                    => null()   !carbon fixation by sun-shaded leaf, [gC d-2 h-1]

  contains
    procedure, public :: Init    =>  plt_photo_init
    procedure, public :: Destroy => plt_photo_destroy
  end type plant_photosyns_type


  type(plant_photosyns_type), public, target :: plt_photo     !plant photosynthesis type

contains

  subroutine plt_photo_init(this)
  class(plant_photosyns_type) :: this

  allocate(this%CuticleResist_pft(JP1));this%CuticleResist_pft=spval
  allocate(this%CanopyMinStomaResistH2O_pft(JP1));this%CanopyMinStomaResistH2O_pft=spval
  allocate(this%LeafO2Solubility_pft(JP1));this%LeafO2Solubility_pft=spval
  allocate(this%RawH2OCanopy2Atm_pft(JP1));this%RawH2OCanopy2Atm_pft=spval
  allocate(this%RawHeatCanopy2Atm_pft(JP1));this%RawHeatCanopy2Atm_pft=spval
  allocate(this%CanPStomaResistH2O_pft(JP1));this%CanPStomaResistH2O_pft=spval
  allocate(this%DiffCO2Atmos2Intracel_pft(JP1));this%DiffCO2Atmos2Intracel_pft=spval
  allocate(this%AirConc_pft(JP1));this%AirConc_pft=spval
  allocate(this%CanopyVcMaxRubisco25C_pft(JP1));this%CanopyVcMaxRubisco25C_pft=0._r8
  allocate(this%CanopyVoMaxRubisco25C_pft(JP1));this%CanopyVoMaxRubisco25C_pft=0._r8
  allocate(this%CanopyVcMaxPEP25C_pft(JP1)); this%CanopyVcMaxPEP25C_pft=0._r8
  allocate(this%ElectronTransptJmax25C_pft(JP1));this%ElectronTransptJmax25C_pft=0._r8
  allocate(this%TFN_Carboxy_pft(JP1));this%TFN_Carboxy_pft=0._r8
  allocate(this%TFN_Oxygen_pft(JP1));this%TFN_Oxygen_pft=0._r8
  allocate(this%TFN_eTranspt_pft(JP1));this%TFN_eTranspt_pft=0._r8
  allocate(this%LeafAreaSunlit_pft(JP1)); this%LeafAreaSunlit_pft=0._r8
  allocate(this%PARSunlit_pft(JP1));this%PARSunlit_pft=0._r8
  allocate(this%PARSunsha_pft(JP1));this%PARSunsha_pft=0._r8
  allocate(this%DynCi2CaRatio_pft(JP1));this%DynCi2CaRatio_pft=0._r8
  allocate(this%CO2Intra_pft(JP1));this%CO2Intra_pft=0._r8
  allocate(this%CO2IntraScal_pft(JP1));this%CO2IntraScal_pft=0._r8
  allocate(this%CH2OSunlit_pft(JP1));this%CH2OSunlit_pft=0._r8
  allocate(this%CH2OSunsha_pft(JP1));this%CH2OSunsha_pft=0._r8
  allocate(this%CO2CuticleResist_pft(JP1));this%CO2CuticleResist_pft=spval
  allocate(this%LeafEffArea_zsec(NumLeafInclinationClasses1,NumCanopyLayers1,MaxNodesPerBranch1,MaxNumBranches,JP1));this%LeafEffArea_zsec=0._r8
  allocate(this%CPOOL3_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CPOOL3_node=spval
  allocate(this%CPOOL4_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CPOOL4_node=spval
  allocate(this%CMassCO2BundleSheath_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CMassCO2BundleSheath_node=spval
  allocate(this%CO2CompenPoint_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CO2CompenPoint_node=spval
  allocate(this%RubiscoCarboxyEff_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%RubiscoCarboxyEff_node=spval
  allocate(this%VcMaxPEPCarboxyRef_brch(MaxNumBranches,JP1));this%VcMaxPEPCarboxyRef_brch=0._r8
  allocate(this%ElectronTransptJmaxRef_brch(MaxNumBranches,JP1));this%ElectronTransptJmaxRef_brch=0._r8
  allocate(this%VcMaxRubiscoRef_brch(MaxNumBranches,JP1)); this%VcMaxRubiscoRef_brch=0._r8
  allocate(this%VoMaxRubiscoRef_brch(MaxNumBranches,JP1)); this%VoMaxRubiscoRef_brch=0._r8
  allocate(this%C4CarboxyEff_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%C4CarboxyEff_node=spval
  allocate(this%LigthSatCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%LigthSatCarboxyRate_node=spval
  allocate(this%LigthSatC4CarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%LigthSatC4CarboxyRate_node=spval
  allocate(this%NutrientCtrlonC4Carboxy_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%NutrientCtrlonC4Carboxy_node=spval
  allocate(this%CMassHCO3BundleSheath_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CMassHCO3BundleSheath_node=spval
  allocate(this%Vmax4RubiscoCarboxy_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%Vmax4RubiscoCarboxy_node=spval
  allocate(this%ProteinCperm2LeafArea_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%ProteinCperm2LeafArea_node=spval
  allocate(this%CO2lmtRubiscoCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CO2lmtRubiscoCarboxyRate_node=spval
  allocate(this%Vmax4PEPCarboxy_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%Vmax4PEPCarboxy_node=spval
  allocate(this%CO2lmtPEPCarboxyRate_node(MaxNodesPerBranch1,MaxNumBranches,JP1));this%CO2lmtPEPCarboxyRate_node=spval
  allocate(this%iPlantPhotosynsType_pft(JP1));this%iPlantPhotosynsType_pft=0
  allocate(this%Km4PEPCarboxy_pft(JP1));this%Km4PEPCarboxy_pft=spval
  allocate(this%O2L_pft(JP1));this%O2L_pft=spval
  allocate(this%CO2Solubility_pft(JP1));this%CO2Solubility_pft=spval
  allocate(this%CanopyGasCO2_pft(JP1));this%CanopyGasCO2_pft=spval
  allocate(this%ChillHours_pft(JP1));this%ChillHours_pft=spval
  allocate(this%SpecLeafChlAct_pft(JP1));this%SpecLeafChlAct_pft=spval
  allocate(this%LeafProtein2Chl_pft(JP1));this%LeafProtein2Chl_pft=spval
  allocate(this%LeafPEP2Protein_pft(JP1));this%LeafPEP2Protein_pft=spval
  allocate(this%fMesophyllChlProtein_pft(JP1));this%fMesophyllChlProtein_pft=spval
  allocate(this%LeafRubisco2Protein_pft(JP1));this%LeafRubisco2Protein_pft=spval
  allocate(this%VmaxPEPCarboxyRef_pft(JP1));this%VmaxPEPCarboxyRef_pft=spval
  allocate(this%VmaxRubOxyRef_pft(JP1));this%VmaxRubOxyRef_pft=spval
  allocate(this%VmaxSpecRubCarboxyRef_pft(JP1));this%VmaxSpecRubCarboxyRef_pft=spval
  allocate(this%XKCO2_pft(JP1));this%XKCO2_pft=spval
  allocate(this%XKO2_pft(JP1));this%XKO2_pft=spval
  allocate(this%RubiscoActivity_brch(MaxNumBranches,JP1));this%RubiscoActivity_brch=1._r8
  allocate(this%GrainFillDowreg_brch(MaxNumBranches,JP1));this%GrainFillDowreg_brch=spval
  allocate(this%aquCO2Intraleaf_pft(JP1));this%aquCO2Intraleaf_pft=spval
  allocate(this%Km4LeafaqCO2_pft(JP1));this%Km4LeafaqCO2_pft=spval
  allocate(this%Km4RubiscoCarboxy_pft(JP1));this%Km4RubiscoCarboxy_pft=spval
  allocate(this%LeafIntracellularCO2_pft(JP1));this%LeafIntracellularCO2_pft=spval
  allocate(this%O2I_pft(JP1));this%O2I_pft=spval
  allocate(this%RCS_pft(JP1));this%RCS_pft=spval
  allocate(this%CanopyCi2CaRatio_pft(JP1));this%CanopyCi2CaRatio_pft=spval
  allocate(this%H2OCuticleResist_pft(JP1));this%H2OCuticleResist_pft=spval

  end subroutine plt_photo_init

  subroutine plt_photo_destroy(this)
  ! Description:
  ! Release every array owned by this PlantPhotosynthesisAPIData workspace.
  class(plant_photosyns_type), intent(inout) :: this ! owned API state; component units and dimensions as declared above
  character(len=*), parameter :: subname='plt_photo_destroy'

  call PrintInfo('beg '//subname)
  if(associated(this%iPlantPhotosynsType_pft))deallocate(this%iPlantPhotosynsType_pft);nullify(this%iPlantPhotosynsType_pft)
  if(associated(this%SpecLeafChlAct_pft))deallocate(this%SpecLeafChlAct_pft);nullify(this%SpecLeafChlAct_pft)
  if(associated(this%LeafProtein2Chl_pft))deallocate(this%LeafProtein2Chl_pft);nullify(this%LeafProtein2Chl_pft)
  if(associated(this%LeafPEP2Protein_pft))deallocate(this%LeafPEP2Protein_pft);nullify(this%LeafPEP2Protein_pft)
  if(associated(this%fMesophyllChlProtein_pft))deallocate(this%fMesophyllChlProtein_pft);nullify(this%fMesophyllChlProtein_pft)
  if(associated(this%LeafRubisco2Protein_pft))deallocate(this%LeafRubisco2Protein_pft);nullify(this%LeafRubisco2Protein_pft)
  if(associated(this%VmaxPEPCarboxyRef_pft))deallocate(this%VmaxPEPCarboxyRef_pft);nullify(this%VmaxPEPCarboxyRef_pft)
  if(associated(this%VmaxRubOxyRef_pft))deallocate(this%VmaxRubOxyRef_pft);nullify(this%VmaxRubOxyRef_pft)
  if(associated(this%VmaxSpecRubCarboxyRef_pft))deallocate(this%VmaxSpecRubCarboxyRef_pft);nullify(this%VmaxSpecRubCarboxyRef_pft)
  if(associated(this%XKCO2_pft))deallocate(this%XKCO2_pft);nullify(this%XKCO2_pft)
  if(associated(this%XKO2_pft))deallocate(this%XKO2_pft);nullify(this%XKO2_pft)
  if(associated(this%RubiscoActivity_brch))deallocate(this%RubiscoActivity_brch);nullify(this%RubiscoActivity_brch)
  if(associated(this%GrainFillDowreg_brch))deallocate(this%GrainFillDowreg_brch);nullify(this%GrainFillDowreg_brch)
  if(associated(this%aquCO2Intraleaf_pft))deallocate(this%aquCO2Intraleaf_pft);nullify(this%aquCO2Intraleaf_pft)
  if(associated(this%O2I_pft))deallocate(this%O2I_pft);nullify(this%O2I_pft)
  if(associated(this%LeafIntracellularCO2_pft))deallocate(this%LeafIntracellularCO2_pft);nullify(this%LeafIntracellularCO2_pft)
  if(associated(this%Km4RubiscoCarboxy_pft))deallocate(this%Km4RubiscoCarboxy_pft);nullify(this%Km4RubiscoCarboxy_pft)
  if(associated(this%Km4LeafaqCO2_pft))deallocate(this%Km4LeafaqCO2_pft);nullify(this%Km4LeafaqCO2_pft)
  if(associated(this%RCS_pft))deallocate(this%RCS_pft);nullify(this%RCS_pft)
  if(associated(this%CanopyCi2CaRatio_pft))deallocate(this%CanopyCi2CaRatio_pft);nullify(this%CanopyCi2CaRatio_pft)
  if(associated(this%H2OCuticleResist_pft))deallocate(this%H2OCuticleResist_pft);nullify(this%H2OCuticleResist_pft)
  if(associated(this%ChillHours_pft))deallocate(this%ChillHours_pft);nullify(this%ChillHours_pft)
  if(associated(this%CO2Solubility_pft))deallocate(this%CO2Solubility_pft);nullify(this%CO2Solubility_pft)
  if(associated(this%CanopyGasCO2_pft))deallocate(this%CanopyGasCO2_pft);nullify(this%CanopyGasCO2_pft)
  if(associated(this%O2L_pft))deallocate(this%O2L_pft);nullify(this%O2L_pft)
  if(associated(this%Km4PEPCarboxy_pft))deallocate(this%Km4PEPCarboxy_pft);nullify(this%Km4PEPCarboxy_pft)
  if(associated(this%CanopyMinStomaResistH2O_pft))deallocate(this%CanopyMinStomaResistH2O_pft);nullify(this%CanopyMinStomaResistH2O_pft)
  if(associated(this%CuticleResist_pft))deallocate(this%CuticleResist_pft);nullify(this%CuticleResist_pft)
  if(associated(this%CanPStomaResistH2O_pft))deallocate(this%CanPStomaResistH2O_pft);nullify(this%CanPStomaResistH2O_pft)
  if(associated(this%RawH2OCanopy2Atm_pft))deallocate(this%RawH2OCanopy2Atm_pft);nullify(this%RawH2OCanopy2Atm_pft)
  if(associated(this%RawHeatCanopy2Atm_pft))deallocate(this%RawHeatCanopy2Atm_pft);nullify(this%RawHeatCanopy2Atm_pft)
  if(associated(this%LeafO2Solubility_pft))deallocate(this%LeafO2Solubility_pft);nullify(this%LeafO2Solubility_pft)
  if(associated(this%CO2CuticleResist_pft))deallocate(this%CO2CuticleResist_pft);nullify(this%CO2CuticleResist_pft)
  if(associated(this%DiffCO2Atmos2Intracel_pft))deallocate(this%DiffCO2Atmos2Intracel_pft);nullify(this%DiffCO2Atmos2Intracel_pft)
  if(associated(this%LeafEffArea_zsec))deallocate(this%LeafEffArea_zsec);nullify(this%LeafEffArea_zsec)
  if(associated(this%CPOOL3_node))deallocate(this%CPOOL3_node);nullify(this%CPOOL3_node)
  if(associated(this%CPOOL4_node))deallocate(this%CPOOL4_node);nullify(this%CPOOL4_node)
  if(associated(this%CMassCO2BundleSheath_node))deallocate(this%CMassCO2BundleSheath_node);nullify(this%CMassCO2BundleSheath_node)
  if(associated(this%CO2CompenPoint_node))deallocate(this%CO2CompenPoint_node);nullify(this%CO2CompenPoint_node)
  if(associated(this%VoMaxRubiscoRef_brch))deallocate(this%VoMaxRubiscoRef_brch);nullify(this%VoMaxRubiscoRef_brch)
  if(associated(this%VcMaxRubiscoRef_brch))deallocate(this%VcMaxRubiscoRef_brch);nullify(this%VcMaxRubiscoRef_brch)
  if(associated(this%RubiscoCarboxyEff_node))deallocate(this%RubiscoCarboxyEff_node);nullify(this%RubiscoCarboxyEff_node)
  if(associated(this%VcMaxPEPCarboxyRef_brch))deallocate(this%VcMaxPEPCarboxyRef_brch);nullify(this%VcMaxPEPCarboxyRef_brch)
  if(associated(this%ElectronTransptJmaxRef_brch))deallocate(this%ElectronTransptJmaxRef_brch);nullify(this%ElectronTransptJmaxRef_brch)
  if(associated(this%C4CarboxyEff_node))deallocate(this%C4CarboxyEff_node);nullify(this%C4CarboxyEff_node)
  if(associated(this%LigthSatCarboxyRate_node))deallocate(this%LigthSatCarboxyRate_node);nullify(this%LigthSatCarboxyRate_node)
  if(associated(this%LigthSatC4CarboxyRate_node))deallocate(this%LigthSatC4CarboxyRate_node);nullify(this%LigthSatC4CarboxyRate_node)
  if(associated(this%NutrientCtrlonC4Carboxy_node))deallocate(this%NutrientCtrlonC4Carboxy_node);nullify(this%NutrientCtrlonC4Carboxy_node)
  if(associated(this%CMassHCO3BundleSheath_node))deallocate(this%CMassHCO3BundleSheath_node);nullify(this%CMassHCO3BundleSheath_node)
  if(associated(this%Vmax4RubiscoCarboxy_node))deallocate(this%Vmax4RubiscoCarboxy_node);nullify(this%Vmax4RubiscoCarboxy_node)
  if(associated(this%ProteinCperm2LeafArea_node))deallocate(this%ProteinCperm2LeafArea_node);nullify(this%ProteinCperm2LeafArea_node)
  if(associated(this%CO2lmtRubiscoCarboxyRate_node))deallocate(this%CO2lmtRubiscoCarboxyRate_node);nullify(this%CO2lmtRubiscoCarboxyRate_node)
  if(associated(this%Vmax4PEPCarboxy_node))deallocate(this%Vmax4PEPCarboxy_node);nullify(this%Vmax4PEPCarboxy_node)
  if(associated(this%CO2lmtPEPCarboxyRate_node))deallocate(this%CO2lmtPEPCarboxyRate_node);nullify(this%CO2lmtPEPCarboxyRate_node)
  if(associated(this%AirConc_pft))deallocate(this%AirConc_pft);nullify(this%AirConc_pft)
  if(associated(this%CanopyVcMaxRubisco25C_pft))deallocate(this%CanopyVcMaxRubisco25C_pft);nullify(this%CanopyVcMaxRubisco25C_pft)
  if(associated(this%CanopyVoMaxRubisco25C_pft))deallocate(this%CanopyVoMaxRubisco25C_pft);nullify(this%CanopyVoMaxRubisco25C_pft)
  if(associated(this%CanopyVcMaxPEP25C_pft))deallocate(this%CanopyVcMaxPEP25C_pft);nullify(this%CanopyVcMaxPEP25C_pft)
  if(associated(this%ElectronTransptJmax25C_pft))deallocate(this%ElectronTransptJmax25C_pft);nullify(this%ElectronTransptJmax25C_pft)
  if(associated(this%TFN_Carboxy_pft))deallocate(this%TFN_Carboxy_pft);nullify(this%TFN_Carboxy_pft)
  if(associated(this%TFN_Oxygen_pft))deallocate(this%TFN_Oxygen_pft);nullify(this%TFN_Oxygen_pft)
  if(associated(this%TFN_eTranspt_pft))deallocate(this%TFN_eTranspt_pft);nullify(this%TFN_eTranspt_pft)
  if(associated(this%LeafAreaSunlit_pft))deallocate(this%LeafAreaSunlit_pft);nullify(this%LeafAreaSunlit_pft)
  if(associated(this%PARSunlit_pft))deallocate(this%PARSunlit_pft);nullify(this%PARSunlit_pft)
  if(associated(this%PARSunsha_pft))deallocate(this%PARSunsha_pft);nullify(this%PARSunsha_pft)
  if(associated(this%DynCi2CaRatio_pft))deallocate(this%DynCi2CaRatio_pft);nullify(this%DynCi2CaRatio_pft)
  if(associated(this%CO2Intra_pft))deallocate(this%CO2Intra_pft);nullify(this%CO2Intra_pft)
  if(associated(this%CO2IntraScal_pft))deallocate(this%CO2IntraScal_pft);nullify(this%CO2IntraScal_pft)
  if(associated(this%CH2OSunlit_pft))deallocate(this%CH2OSunlit_pft);nullify(this%CH2OSunlit_pft)
  if(associated(this%CH2OSunsha_pft))deallocate(this%CH2OSunsha_pft);nullify(this%CH2OSunsha_pft)
  call PrintInfo('end '//subname)
  end subroutine plt_photo_destroy
end module PlantPhotosynthesisAPIData
