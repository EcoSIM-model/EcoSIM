module PlantRootBGCAPIData
  ! Description:
  ! Owns the rootbgc API type, instance and lifecycle procedures.
  use DebugToolMod, only: PrintInfo
  use PlantAPICommonData
  implicit none
  private
  public :: plt_rootbgc_init, plt_rootbgc_destroy

  type, public :: plant_rootbgc_type
  real(r8), pointer :: canopy_growth_pft(:)              => null()  !canopy structural C growth rate,                                [gC d-2 h-1]
  real(r8), pointer :: RootSoilCO2Debt_vr(:)             => null()  !pending local soil CO2 debt, paired with mass on handoff, [gC d-2]
  real(r8), pointer :: TRootGasLossDisturb_col(:)        => null()  !total root gas content,                                         [g d-2]
  real(r8), pointer :: trcs_Soil2plant_uptake_vr(:,:)    => null()  !total root-soil solute flux non-band,                           [g d-2 h-1]
  real(r8), pointer :: trcs_Soil2plant_uptake_pvr(:,:,:) => null()  !plant root-soil solute flux non-band,                           [g d-2 h-1]
  real(r8), pointer :: Soil2RootMycoExudE_pft(:,:)         => null()  !total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
  real(r8), pointer :: CanopyN2Fix_pft(:)                => null()  !total canopy N2 fixation, [g d-2 h-1]
  real(r8), pointer :: RootN2Fix_pft(:)                  => null()  !total root N2 fixation,                                         [g d-2 h-1]
  real(r8), pointer :: RootNO3Uptake_pft(:)              => null()  !total root uptake of NO3,                                       [g d-2 h-1]
  real(r8), pointer :: RootNH4Uptake_pft(:)              => null()  !total root uptake of NH4,                                       [g d-2 h-1]
  real(r8), pointer :: RootHPO4Uptake_pft(:)             => null()  !total root uptake of HPO4,                                      [g d-2 h-1]
  real(r8), pointer :: RootH2PO4Uptake_pft(:)            => null()  !total root uptake of PO4,                                       [g d-2 h-1]
  real(r8), pointer :: Soil2RootMycoExudE_pvr(:,:,:,:,:)  => null()  !root uptake (+ve) - exudation (-ve) of DOE,                     [g d-2 h-1]
  real(r8), pointer :: PlantRootSoilElmNetX_pft(:,:)     => null()  !net root element uptake (+ve) - exudation (-ve),                [gC d-2 h-1]
  real(r8), pointer :: REcoUptkSoilO2M_vr(:,:)              => null()  !total O2 sink,                                                  [g d-2 t-1]
  real(r8), pointer :: ZERO4Uptk_pft(:)                  => null()  !threshold zero for uptake calculation,                          [-]
  real(r8), pointer :: CMinPO4Root_pft(:,:)              => null()  !minimum PO4 concentration for root NH4 uptake,                  [g m-3]
  real(r8), pointer :: VmaxPO4Root_pft(:,:)              => null()  !maximum root PO4 uptake rate,                                   [g m-2 h-1]
  real(r8), pointer :: KmPO4Root_pft(:,:)                => null()  !Km for root PO4 uptake,                                         [g m-3]
  real(r8), pointer :: CminNO3Root_pft(:,:)              => null()  !minimum NO3 concentration for root NH4 uptake,                  [g m-3]
  real(r8), pointer :: VmaxNO3Root_pft(:,:)              => null()  !maximum root NO3 uptake rate,                                   [g m-2 h-1]
  real(r8), pointer :: KmNO3Root_pft(:,:)                => null()  !Km for root NO3 uptake,                                         [g m-3]
  real(r8), pointer :: CMinNH4Root_pft(:,:)              => null()  !minimum NH4 concentration for root NH4 uptake,                  [g m-3]
  real(r8), pointer :: VmaxNH4Root_pft(:,:)              => null()  !maximum root NH4 uptake rate,                                   [g m-2 h-1]
  real(r8), pointer :: KmNH4Root_pft(:,:)                => null()  !Km for root NH4 uptake,                                         [g m-3]
  real(r8), pointer :: RCO2Emis2Root_rpvr(:,:,:)          => null()  !aqueous CO2 flux from roots to root water,                      [g d-2 h-1]
  real(r8), pointer :: RootO2Uptk_pvr(:,:,:)             => null()  !aqueous O2 flux from roots to root water,                       [g d-2 h-1]
  real(r8), pointer :: RootUptkSoiSol_pvr(:,:,:,:)       => null()  !aqueous CO2 flux from roots to soil water,                      [g d-2 h-1]
  real(r8), pointer :: trcg_air2root_flx_pvr(:,:,:,:)    => null()  !gaseous tracer flux through roots,                              [g d-2 h-1]
  real(r8), pointer :: trcg_Root_gas2aqu_flx_vr(:,:,:,:) => null()  !dissolution (+ve) - volatilization (-ve) gas flux in roots,     [g d-2 h-1]
  real(r8), pointer :: RootO2Dmnd4Resp_pvr(:,:,:)        => null()  !root  O2 demand from respiration,                               [g d-2 h-1]

  real(r8), pointer :: RootNH4DmndSoilPrev_pvr(:,:,:)    => null()  !previous root uptake of NH4 non-band unconstrained by NH4,      [g d-2 h-1]
  real(r8), pointer :: RootNH4DmndBandPrev_pvr(:,:,:)    => null()  !previous root uptake of NO3 band unconstrained by NO3,          [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndSoilPrev_pvr(:,:,:)    => null()  !previous root uptake of NH4 band unconstrained by NH4,          [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndBandPrev_pvr(:,:,:)    => null()  !previous root uptake of NO3 non-band unconstrained by NO3,      [g d-2 h-1]
  real(r8), pointer :: VmaxNH4Root_pvr(:,:,:)            => null()  !maximum NH4 uptake rate,                                        [gN h-1 (gC root)-1]
  real(r8), pointer :: VmaxNO3Root_pvr(:,:,:)            => null()  !maximum NO3 uptake rate,                                        [gN h-1 (gC root)-1]
  real(r8), pointer :: RootNH4DmndSoil_pvr(:,:,:)        => null()  !root uptake of NH4 non-band unconstrained by NH4,               [g d-2 h-1]
  real(r8), pointer :: RootNH4DmndBand_pvr(:,:,:)        => null()  !root uptake of NO3 band unconstrained by NO3,                   [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndSoil_pvr(:,:,:)        => null()  !root uptake of NH4 band unconstrained by NH4,                   [g d-2 h-1]
  real(r8), pointer :: RootNO3DmndBand_pvr(:,:,:)        => null()  !root uptake of NO3 non-band unconstrained by NO3,               [g d-2 h-1]
  real(r8), pointer :: RootH2PO4DmndSoil_pvr(:,:,:)      => null()  !root uptake of H2PO4 non-band,                                  [g d-2 h-1]
  real(r8), pointer :: RootH2PO4DmndBand_pvr(:,:,:)      => null()  !root uptake of H2PO4 band,                                      [g d-2 h-1]
  real(r8), pointer :: RootH1PO4DmndSoil_pvr(:,:,:)      => null()  !HPO4 demand in non-band by each root population,                [g d-2 h-1]
  real(r8), pointer :: RootH1PO4DmndBand_pvr(:,:,:)      => null()  !HPO4 demand in band by each root population,                    [g d-2 h-1]
  real(r8), pointer :: RootH2PO4DmndSoilPrev_pvr(:,:,:)  => null()  !previous step root uptake of H2PO4 non-band,                    [g d-2 h-1]
  real(r8), pointer :: RootH2PO4DmndBandPrev_pvr(:,:,:)  => null()  !previous step root uptake of H2PO4 band,                        [g d-2 h-1]
  real(r8), pointer :: RootH1PO4DmndSoilPrev_pvr(:,:,:)  => null()  !previous step HPO4 demand in non-band by each root population,  [g d-2 h-1]
  real(r8), pointer :: RootH1PO4DmndBandPrev_pvr(:,:,:)  => null()  !previous step HPO4 demand in band by each root population,      [g d-2 h-1]

  real(r8), pointer :: RAutoRootO2Limter_rpvr(:,:,:)     => null()  !O2 constraint to root respiration (0-1),                        [-]
  real(r8), pointer :: trcg_rootml_pvr(:,:,:,:)          => null() !root gas content,                                                [g d-2]
  real(r8), pointer :: trcs_rootml_pvr(:,:,:,:)          => null() !root aqueous content,                                            [g d-2]
  real(r8), pointer :: RootAtmGasConductance_rpvr(:,:,:,:)   => null()  !Conductance for gas diffusion                               [m3 d-2 h-1]
  real(r8), pointer :: RNodeInitiate_pft(:)               => null()    !node initiation rate, [h-1]
  real(r8), pointer :: RLeafAppear_pft(:)                 => null()    !leaf appearing rate, [h-1]
  real(r8), pointer :: NH3Dep2Can_brch(:,:)              => null()  !gaseous NH3 flux fron root disturbance band,                    [g d-2 h-1]
  real(r8), pointer :: GPP_brch(:,:)                    => null()  !dGPP (C4-C3 product) over branch, [gC d-2 h-1]
  real(r8), pointer :: CytokininMRConc_rpvr(:,:,:)      => null()  !cytokinin concentration in medium size roots, [gC m-3 H2O], [g d-2 h-1]
  real(r8), pointer :: Cytokinin1stConc_rpvr(:,:,:)     => null()   !cytokinin concentration in primary roots, [gC m-3 H2O]
  real(r8), pointer :: Cytokinin2ndConc_rpvr(:,:,:,:)    => null()  !cytokinin concentration in fine roots, [gC m-3 H2O]
  real(r8), pointer :: RootNutUptake_pvr(:,:,:,:)        => null()  !root uptake of Nutrient band,                                   [g d-2 h-1]
  real(r8), pointer :: RootOUlmNutUptake_pvr(:,:,:,:)    => null()  !root uptake of NH4 band unconstrained by O2,                    [g d-2 h-1]
  real(r8), pointer :: RootCUlmNutUptake_pvr(:,:,:,:)    => null()  !root uptake of NH4 band unconstrained by root nonstructural C,  [g d-2 h-1]
  real(r8), pointer :: RootRespPotent_pvr(:,:,:)         => null()  !root respiration unconstrained by O2,                           [g d-2 h-1]
  real(r8), pointer :: RootMRProdCytok_rpvr(:,:,:)       => null()  !cytokinin production rate due to medium root metabolism, [gC d-2 h-1]
  real(r8), pointer :: Root2ndProdCytok_rpvr(:,:,:,:)    => null()  !cytokinin production rate due to fine root/myco elongation, [gC d-2 h-1]
  real(r8), pointer :: RootMyco2ndSinkC_rpvr(:,:,:,:)    => null()  !fine root/myco carbon sink, [gC d-2 h-1]
  real(r8), pointer :: RootMyco1stSinkC_rpvr(:,:,:)      => null()  !primary root C sink, [gC d-2 h-1]
  real(r8), pointer :: RootCO2EmisPot_pvr(:,:,:)         => null()  !root CO2 efflux unconstrained by root nonstructural C,          [g d-2 h-1]
  real(r8), pointer :: RootCO2Autor_pvr(:,:,:)           => null()  !root respiration constrained by O2,                             [g d-2 h-1]
  real(r8), pointer :: TurgEff4CanopyResp_pft(:)         => null()  !Turgor pressure effect on canopy respiration, [-]
  real(r8), pointer :: RootCO2AutorX_pvr(:,:,:)          => null()  !root respiration from previous time step,                       [g d-2 h-1]
  real(r8), pointer :: PlantExudElm_CumYr_pft(:,:)       => null()  !total net root element uptake (+ve) - exudation (-ve),          [gC d-2 ]
  real(r8), pointer :: trcg_root_vr(:,:)                 => null()   !total root internal gas flux,                                  [g d-2 h-1]
  real(r8), pointer :: trcg_air2root_flx_vr(:,:)         => null()   !total internal root gas flux,                                  [gC d-2 h-1]
  real(r8), pointer :: CO2FixCL_pft(:)                   => null()   !Rubisco-limited CO2 fixation,                                  [gC d-2 h-1]
  real(r8), pointer :: RootNutUptakeN_pft(:)             => null()   !total N uptake by plant roots,                                 [gN d-h2 h-1]
  real(r8), pointer :: RootNutUptakeP_pft(:)             => null()   !total P uptake by plant roots,                                 [gP d-h2 h-1]
  real(r8), pointer :: CO2FixLL_pft(:)                   => null()   !Light-limited CO2 fixation,                                    [gC d-h2 h-1]
  real(r8), pointer :: RootUptk_N_CumYr_pft(:)           => null()  !cumulative plant N uptake,                                      [gN d-2]
  real(r8), pointer :: RootUptk_P_CumYr_pft(:)           => null()  !cumulative plant P uptake,                                      [gP d-2]
  real(r8), pointer :: RootCO2Ar2Soil_pvr(:,:)             => null()  !root respiration released to soil,                              [gC d-2 h-1]
  real(r8), pointer :: RootCO2Ar2RootX_pvr(:,:)          => null()  !root respiration released to root,                              [gC d-2 h-1]
  real(r8), pointer :: RootCO2Ar2RootX_rpvr(:,:,:)       => null()  !root/myco respiration released to root/myco,                    [gC d-2 h-1]
  real(r8), pointer :: trcs_deadroot2soil_pvr(:,:,:)     => null()  !gases released to soil upong dying roots,                       [g d-2 h-1]
  real(r8), pointer :: GroSrcRootStress_pvr(:,:)        => null()  !root growth stress due to nutrient and water, [-]
  contains
    procedure, public :: Init => plt_rootbgc_init
    procedure, public :: Destroy  => plt_rootbgc_destroy
  end type plant_rootbgc_type


  type(plant_rootbgc_type)  , public, target :: plt_rbgc      !root bgc

contains

  subroutine plt_rootbgc_init(this)

  implicit none
  class(plant_rootbgc_type) :: this
  allocate(this%trcs_Soil2plant_uptake_pvr(ids_beg:ids_end,JZ1,JP1)); this%trcs_Soil2plant_uptake_pvr=0._r8
  allocate(this%trcs_Soil2plant_uptake_vr(ids_beg:ids_end,JZ1)); this%trcs_Soil2plant_uptake_vr=0._r8
  allocate(this%trcg_rootml_pvr(idg_beg:idg_NH3,jroots,JZ1,JP1));this%trcg_rootml_pvr=spval
  allocate(this%trcs_rootml_pvr(idg_beg:idg_NH3,jroots,JZ1,JP1));this%trcs_rootml_pvr=spval
  allocate(this%RootAtmGasConductance_rpvr(idg_beg:idg_NH3,jroots,JZ1,JP1)); this%RootAtmGasConductance_rpvr=0._r8
  allocate(this%RootSoilCO2Debt_vr(0:JZ1));this%RootSoilCO2Debt_vr=0._r8
  allocate(this%TRootGasLossDisturb_col(idg_beg:idg_NH3));this%TRootGasLossDisturb_col=spval
  allocate(this%REcoUptkSoilO2M_vr(60,0:JZ1)); this%REcoUptkSoilO2M_vr=spval
  allocate(this%Soil2RootMycoExudE_pvr(NumPlantChemElms,jroots,1:jcplx,0:JZ1,JP1));this%Soil2RootMycoExudE_pvr=spval
  allocate(this%PlantRootSoilElmNetX_pft(NumPlantChemElms,JP1)); this%PlantRootSoilElmNetX_pft=spval
  allocate(this%PlantExudElm_CumYr_pft(NumPlantChemElms,JP1));this%PlantExudElm_CumYr_pft=spval
  allocate(this%RootUptk_N_CumYr_pft(JP1)); this%RootUptk_N_CumYr_pft=spval
  allocate(this%RootUptk_P_CumYr_pft(JP1)); this%RootUptk_P_CumYr_pft=spval
  allocate(this%Soil2RootMycoExudE_pft(NumPlantChemElms,JP1));this%Soil2RootMycoExudE_pft=spval
  allocate(this%RootN2Fix_pft(JP1)); this%RootN2Fix_pft=spval
  allocate(this%CanopyN2Fix_pft(JP1));this%CanopyN2Fix_pft=spval
  allocate(this%canopy_growth_pft(JP1)); this%canopy_growth_pft=spval
  allocate(this%RootNO3Uptake_pft(JP1)); this%RootNO3Uptake_pft=0._r8
  allocate(this%RootNH4Uptake_pft(JP1)); this%RootNH4Uptake_pft=0._r8
  allocate(this%RootHPO4Uptake_pft(JP1)); this%RootHPO4Uptake_pft=0._r8
  allocate(this%RootH2PO4Uptake_pft(JP1)); this%RootH2PO4Uptake_pft=0._r8

  allocate(this%ZERO4Uptk_pft(JP1)); this%ZERO4Uptk_pft=spval
  allocate(this%RootRespPotent_pvr(jroots,JZ1,JP1)); this%RootRespPotent_pvr=spval
  allocate(this%RootCO2EmisPot_pvr(jroots,JZ1,JP1)); this%RootCO2EmisPot_pvr=spval
  allocate(this%RootCO2Autor_pvr(jroots,JZ1,JP1)); this%RootCO2Autor_pvr=0._r8
  allocate(this%TurgEff4CanopyResp_pft(JP1)); this%TurgEff4CanopyResp_pft=0._r8
  allocate(this%trcs_deadroot2soil_pvr(idg_beg:idg_NH3,JZ1,JP1));this%trcs_deadroot2soil_pvr=0._r8
  allocate(this%RootCO2Ar2Soil_pvr(JZ1,JP1)); this%RootCO2Ar2Soil_pvr=0._r8
  allocate(this%RootCO2Ar2RootX_pvr(JZ1,JP1)); this%RootCO2Ar2RootX_pvr=0._r8
  allocate(this%RootCO2Ar2RootX_rpvr(jroots,JZ1,JP1));this%RootCO2Ar2RootX_rpvr=0._r8
  allocate(this%RootCO2AutorX_pvr(jroots,JZ1,JP1)); this%RootCO2AutorX_pvr=spval
  allocate(this%RootMyco2ndSinkC_rpvr(jroots,JZ1,MaxNumRootAxes,JP1)); this%RootMyco2ndSinkC_rpvr=0._r8
  allocate(this%RootMyco1stSinkC_rpvr(JZ1,MaxNumRootAxes,JP1)); this%RootMyco1stSinkC_rpvr=0._r8
  allocate(this%RootMRProdCytok_rpvr(JZ1,MaxNumRootAxes,JP1)); this%RootMRProdCytok_rpvr=0._r8
  allocate(this%Root2ndProdCytok_rpvr(jroots,JZ1,MaxNumRootAxes,JP1)); this%Root2ndProdCytok_rpvr=0._r8
  allocate(this%RootNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ1,JP1)); this%RootNutUptake_pvr=0._r8
  allocate(this%Cytokinin2ndConc_rpvr(jroots,JZ1,MaxNumRootAxes,JP1));this%Cytokinin2ndConc_rpvr=0._r8
  allocate(this%Cytokinin1stConc_rpvr(JZ1,MaxNumRootAxes,JP1)); this%Cytokinin1stConc_rpvr=0._r8
  allocate(this%CytokininMRConc_rpvr(JZ1,MaxNumRootAxes,JP1)); this%CytokininMRConc_rpvr=0._r8
  allocate(this%RootOUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ1,JP1));this%RootOUlmNutUptake_pvr=spval
  allocate(this%RootCUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ1,JP1));this%RootCUlmNutUptake_pvr=spval
  allocate(this%NH3Dep2Can_brch(MaxNumBranches,JP1));this%NH3Dep2Can_brch=0._r8
  allocate(this%GPP_brch(MaxNumBranches,JP1)); this%GPP_brch=spval
  allocate(this%RNodeInitiate_pft(JP1));this%RNodeInitiate_pft=spval
  allocate(this%RLeafAppear_pft(JP1));this%RLeafAppear_pft=spval
  allocate(this%CO2FixCL_pft(JP1)); this%CO2FixCL_pft=spval
  allocate(this%CO2FixLL_pft(JP1)); this%CO2FixLL_pft=spval
  allocate(this%GroSrcRootStress_pvr(JZ1,JP1));this%GroSrcRootStress_pvr=1._R8
  allocate(this%RootNutUptakeN_pft(JP1));this%RootNutUptakeN_pft=spval
  allocate(this%RootNutUptakeP_pft(JP1));this%RootNutUptakeP_pft=spval
  allocate(this%trcg_air2root_flx_vr(idg_beg:idg_NH3,JZ1));this%trcg_air2root_flx_vr=spval
  allocate(this%trcg_root_vr(idg_beg:idg_NH3,JZ1));this%trcg_root_vr=spval

  allocate(this%trcg_air2root_flx_pvr(idg_beg:idg_NH3,jroots,JZ1,JP1));this%trcg_air2root_flx_pvr=spval
  allocate(this%trcg_Root_gas2aqu_flx_vr(idg_beg:idg_NH3,jroots,JZ1,JP1));this%trcg_Root_gas2aqu_flx_vr=spval
  allocate(this%RootO2Dmnd4Resp_pvr(jroots,JZ1,JP1));this%RootO2Dmnd4Resp_pvr=spval
  allocate(this%RootNH4DmndSoil_pvr(jroots,JZ1,JP1));this%RootNH4DmndSoil_pvr=spval
  allocate(this%VmaxNH4Root_pvr(jroots,JZ1,JP1)); this%VmaxNH4Root_pvr=0._r8
  allocate(this%VmaxNO3Root_pvr(jroots,JZ1,JP1)); this%VmaxNO3Root_pvr=0._r8
  allocate(this%RootNH4DmndBand_pvr(jroots,JZ1,JP1));this%RootNH4DmndBand_pvr=spval
  allocate(this%RootNO3DmndSoil_pvr(jroots,JZ1,JP1));this%RootNO3DmndSoil_pvr=spval
  allocate(this%RootNO3DmndBand_pvr(jroots,JZ1,JP1));this%RootNO3DmndBand_pvr=spval
  allocate(this%RootH2PO4DmndSoil_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndSoil_pvr=spval
  allocate(this%RootH2PO4DmndBand_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndBand_pvr=spval
  allocate(this%RootH1PO4DmndSoil_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndSoil_pvr=spval
  allocate(this%RootH1PO4DmndBand_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndBand_pvr=spval

  allocate(this%RootNH4DmndSoilPrev_pvr(jroots,JZ1,JP1)); this%RootNH4DmndSoilPrev_pvr=0._r8
  allocate(this%RootNH4DmndBandPrev_pvr(jroots,JZ1,JP1)); this%RootNH4DmndBandPrev_pvr=0._r8
  allocate(this%RootNO3DmndSoilPrev_pvr(jroots,JZ1,JP1)); this%RootNO3DmndSoilPrev_pvr=0._r8
  allocate(this%RootNO3DmndBandPrev_pvr(jroots,JZ1,JP1)); this%RootNO3DmndBandPrev_pvr=0._r8
  allocate(this%RootH2PO4DmndSoilPrev_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndSoilPrev_pvr=0._r8
  allocate(this%RootH2PO4DmndBandPrev_pvr(jroots,JZ1,JP1));this%RootH2PO4DmndBandPrev_pvr=0._r8
  allocate(this%RootH1PO4DmndSoilPrev_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndSoilPrev_pvr=0._r8
  allocate(this%RootH1PO4DmndBandPrev_pvr(jroots,JZ1,JP1));this%RootH1PO4DmndBandPrev_pvr=0._r8

  allocate(this%RAutoRootO2Limter_rpvr(jroots,JZ1,JP1));this%RAutoRootO2Limter_rpvr=spval
  allocate(this%CMinPO4Root_pft(jroots,JP1));this%CMinPO4Root_pft=spval
  allocate(this%VmaxPO4Root_pft(jroots,JP1));this%VmaxPO4Root_pft=spval
  allocate(this%KmPO4Root_pft(jroots,JP1));this%KmPO4Root_pft=spval
  allocate(this%CminNO3Root_pft(jroots,JP1));this%CminNO3Root_pft=spval
  allocate(this%VmaxNO3Root_pft(jroots,JP1));this%VmaxNO3Root_pft=spval
  allocate(this%KmNO3Root_pft(jroots,JP1));this%KmNO3Root_pft=spval
  allocate(this%CMinNH4Root_pft(jroots,JP1));this%CMinNH4Root_pft=spval
  allocate(this%VmaxNH4Root_pft(jroots,JP1));this%VmaxNH4Root_pft=spval
  allocate(this%KmNH4Root_pft(jroots,JP1));this%KmNH4Root_pft=spval
  allocate(this%RCO2Emis2Root_rpvr(jroots,JZ1,JP1));this%RCO2Emis2Root_rpvr=spval
  allocate(this%RootO2Uptk_pvr(jroots,JZ1,JP1));this%RootO2Uptk_pvr=spval
  allocate(this%RootUptkSoiSol_pvr(idg_beg:idg_end,jroots,JZ1,JP1));this%RootUptkSoiSol_pvr=0._r8
  end subroutine plt_rootbgc_init

  subroutine plt_rootbgc_destroy(this)
  ! Description:
  ! Release every array owned by this PlantRootBGCAPIData workspace.
  class(plant_rootbgc_type), intent(inout) :: this ! owned API state; component units and dimensions as declared above
  character(len=*), parameter :: subname='plt_rootbgc_destroy'

  call PrintInfo('beg '//subname)
  if(associated(this%canopy_growth_pft))deallocate(this%canopy_growth_pft);nullify(this%canopy_growth_pft)
  if(associated(this%RootSoilCO2Debt_vr))deallocate(this%RootSoilCO2Debt_vr);nullify(this%RootSoilCO2Debt_vr)
  if(associated(this%TRootGasLossDisturb_col))deallocate(this%TRootGasLossDisturb_col);nullify(this%TRootGasLossDisturb_col)
  if(associated(this%trcs_Soil2plant_uptake_vr))deallocate(this%trcs_Soil2plant_uptake_vr);nullify(this%trcs_Soil2plant_uptake_vr)
  if(associated(this%trcs_Soil2plant_uptake_pvr))deallocate(this%trcs_Soil2plant_uptake_pvr);nullify(this%trcs_Soil2plant_uptake_pvr)
  if(associated(this%Soil2RootMycoExudE_pft))deallocate(this%Soil2RootMycoExudE_pft);nullify(this%Soil2RootMycoExudE_pft)
  if(associated(this%CanopyN2Fix_pft))deallocate(this%CanopyN2Fix_pft);nullify(this%CanopyN2Fix_pft)
  if(associated(this%RootN2Fix_pft))deallocate(this%RootN2Fix_pft);nullify(this%RootN2Fix_pft)
  if(associated(this%RootNO3Uptake_pft))deallocate(this%RootNO3Uptake_pft);nullify(this%RootNO3Uptake_pft)
  if(associated(this%RootNH4Uptake_pft))deallocate(this%RootNH4Uptake_pft);nullify(this%RootNH4Uptake_pft)
  if(associated(this%RootHPO4Uptake_pft))deallocate(this%RootHPO4Uptake_pft);nullify(this%RootHPO4Uptake_pft)
  if(associated(this%RootH2PO4Uptake_pft))deallocate(this%RootH2PO4Uptake_pft);nullify(this%RootH2PO4Uptake_pft)
  if(associated(this%Soil2RootMycoExudE_pvr))deallocate(this%Soil2RootMycoExudE_pvr);nullify(this%Soil2RootMycoExudE_pvr)
  if(associated(this%PlantRootSoilElmNetX_pft))deallocate(this%PlantRootSoilElmNetX_pft);nullify(this%PlantRootSoilElmNetX_pft)
  if(associated(this%REcoUptkSoilO2M_vr))deallocate(this%REcoUptkSoilO2M_vr);nullify(this%REcoUptkSoilO2M_vr)
  if(associated(this%ZERO4Uptk_pft))deallocate(this%ZERO4Uptk_pft);nullify(this%ZERO4Uptk_pft)
  if(associated(this%CMinPO4Root_pft))deallocate(this%CMinPO4Root_pft);nullify(this%CMinPO4Root_pft)
  if(associated(this%VmaxPO4Root_pft))deallocate(this%VmaxPO4Root_pft);nullify(this%VmaxPO4Root_pft)
  if(associated(this%KmPO4Root_pft))deallocate(this%KmPO4Root_pft);nullify(this%KmPO4Root_pft)
  if(associated(this%CminNO3Root_pft))deallocate(this%CminNO3Root_pft);nullify(this%CminNO3Root_pft)
  if(associated(this%VmaxNO3Root_pft))deallocate(this%VmaxNO3Root_pft);nullify(this%VmaxNO3Root_pft)
  if(associated(this%KmNO3Root_pft))deallocate(this%KmNO3Root_pft);nullify(this%KmNO3Root_pft)
  if(associated(this%CMinNH4Root_pft))deallocate(this%CMinNH4Root_pft);nullify(this%CMinNH4Root_pft)
  if(associated(this%VmaxNH4Root_pft))deallocate(this%VmaxNH4Root_pft);nullify(this%VmaxNH4Root_pft)
  if(associated(this%KmNH4Root_pft))deallocate(this%KmNH4Root_pft);nullify(this%KmNH4Root_pft)
  if(associated(this%RCO2Emis2Root_rpvr))deallocate(this%RCO2Emis2Root_rpvr);nullify(this%RCO2Emis2Root_rpvr)
  if(associated(this%RootO2Uptk_pvr))deallocate(this%RootO2Uptk_pvr);nullify(this%RootO2Uptk_pvr)
  if(associated(this%RootUptkSoiSol_pvr))deallocate(this%RootUptkSoiSol_pvr);nullify(this%RootUptkSoiSol_pvr)
  if(associated(this%trcg_air2root_flx_pvr))deallocate(this%trcg_air2root_flx_pvr);nullify(this%trcg_air2root_flx_pvr)
  if(associated(this%trcg_Root_gas2aqu_flx_vr))deallocate(this%trcg_Root_gas2aqu_flx_vr);nullify(this%trcg_Root_gas2aqu_flx_vr)
  if(associated(this%RootO2Dmnd4Resp_pvr))deallocate(this%RootO2Dmnd4Resp_pvr);nullify(this%RootO2Dmnd4Resp_pvr)
  if(associated(this%RootNH4DmndSoilPrev_pvr))deallocate(this%RootNH4DmndSoilPrev_pvr);nullify(this%RootNH4DmndSoilPrev_pvr)
  if(associated(this%RootNH4DmndBandPrev_pvr))deallocate(this%RootNH4DmndBandPrev_pvr);nullify(this%RootNH4DmndBandPrev_pvr)
  if(associated(this%RootNO3DmndSoilPrev_pvr))deallocate(this%RootNO3DmndSoilPrev_pvr);nullify(this%RootNO3DmndSoilPrev_pvr)
  if(associated(this%RootNO3DmndBandPrev_pvr))deallocate(this%RootNO3DmndBandPrev_pvr);nullify(this%RootNO3DmndBandPrev_pvr)
  if(associated(this%VmaxNH4Root_pvr))deallocate(this%VmaxNH4Root_pvr);nullify(this%VmaxNH4Root_pvr)
  if(associated(this%VmaxNO3Root_pvr))deallocate(this%VmaxNO3Root_pvr);nullify(this%VmaxNO3Root_pvr)
  if(associated(this%RootNH4DmndSoil_pvr))deallocate(this%RootNH4DmndSoil_pvr);nullify(this%RootNH4DmndSoil_pvr)
  if(associated(this%RootNH4DmndBand_pvr))deallocate(this%RootNH4DmndBand_pvr);nullify(this%RootNH4DmndBand_pvr)
  if(associated(this%RootNO3DmndSoil_pvr))deallocate(this%RootNO3DmndSoil_pvr);nullify(this%RootNO3DmndSoil_pvr)
  if(associated(this%RootNO3DmndBand_pvr))deallocate(this%RootNO3DmndBand_pvr);nullify(this%RootNO3DmndBand_pvr)
  if(associated(this%RootH2PO4DmndSoil_pvr))deallocate(this%RootH2PO4DmndSoil_pvr);nullify(this%RootH2PO4DmndSoil_pvr)
  if(associated(this%RootH2PO4DmndBand_pvr))deallocate(this%RootH2PO4DmndBand_pvr);nullify(this%RootH2PO4DmndBand_pvr)
  if(associated(this%RootH1PO4DmndSoil_pvr))deallocate(this%RootH1PO4DmndSoil_pvr);nullify(this%RootH1PO4DmndSoil_pvr)
  if(associated(this%RootH1PO4DmndBand_pvr))deallocate(this%RootH1PO4DmndBand_pvr);nullify(this%RootH1PO4DmndBand_pvr)
  if(associated(this%RootH2PO4DmndSoilPrev_pvr))deallocate(this%RootH2PO4DmndSoilPrev_pvr);nullify(this%RootH2PO4DmndSoilPrev_pvr)
  if(associated(this%RootH2PO4DmndBandPrev_pvr))deallocate(this%RootH2PO4DmndBandPrev_pvr);nullify(this%RootH2PO4DmndBandPrev_pvr)
  if(associated(this%RootH1PO4DmndSoilPrev_pvr))deallocate(this%RootH1PO4DmndSoilPrev_pvr);nullify(this%RootH1PO4DmndSoilPrev_pvr)
  if(associated(this%RootH1PO4DmndBandPrev_pvr))deallocate(this%RootH1PO4DmndBandPrev_pvr);nullify(this%RootH1PO4DmndBandPrev_pvr)
  if(associated(this%RAutoRootO2Limter_rpvr))deallocate(this%RAutoRootO2Limter_rpvr);nullify(this%RAutoRootO2Limter_rpvr)
  if(associated(this%trcg_rootml_pvr))deallocate(this%trcg_rootml_pvr);nullify(this%trcg_rootml_pvr)
  if(associated(this%trcs_rootml_pvr))deallocate(this%trcs_rootml_pvr);nullify(this%trcs_rootml_pvr)
  if(associated(this%RootAtmGasConductance_rpvr))deallocate(this%RootAtmGasConductance_rpvr);nullify(this%RootAtmGasConductance_rpvr)
  if(associated(this%RNodeInitiate_pft))deallocate(this%RNodeInitiate_pft);nullify(this%RNodeInitiate_pft)
  if(associated(this%RLeafAppear_pft))deallocate(this%RLeafAppear_pft);nullify(this%RLeafAppear_pft)
  if(associated(this%NH3Dep2Can_brch))deallocate(this%NH3Dep2Can_brch);nullify(this%NH3Dep2Can_brch)
  if(associated(this%GPP_brch))deallocate(this%GPP_brch);nullify(this%GPP_brch)
  if(associated(this%CytokininMRConc_rpvr))deallocate(this%CytokininMRConc_rpvr);nullify(this%CytokininMRConc_rpvr)
  if(associated(this%Cytokinin1stConc_rpvr))deallocate(this%Cytokinin1stConc_rpvr);nullify(this%Cytokinin1stConc_rpvr)
  if(associated(this%Cytokinin2ndConc_rpvr))deallocate(this%Cytokinin2ndConc_rpvr);nullify(this%Cytokinin2ndConc_rpvr)
  if(associated(this%RootNutUptake_pvr))deallocate(this%RootNutUptake_pvr);nullify(this%RootNutUptake_pvr)
  if(associated(this%RootOUlmNutUptake_pvr))deallocate(this%RootOUlmNutUptake_pvr);nullify(this%RootOUlmNutUptake_pvr)
  if(associated(this%RootCUlmNutUptake_pvr))deallocate(this%RootCUlmNutUptake_pvr);nullify(this%RootCUlmNutUptake_pvr)
  if(associated(this%RootRespPotent_pvr))deallocate(this%RootRespPotent_pvr);nullify(this%RootRespPotent_pvr)
  if(associated(this%RootMRProdCytok_rpvr))deallocate(this%RootMRProdCytok_rpvr);nullify(this%RootMRProdCytok_rpvr)
  if(associated(this%Root2ndProdCytok_rpvr))deallocate(this%Root2ndProdCytok_rpvr);nullify(this%Root2ndProdCytok_rpvr)
  if(associated(this%RootMyco2ndSinkC_rpvr))deallocate(this%RootMyco2ndSinkC_rpvr);nullify(this%RootMyco2ndSinkC_rpvr)
  if(associated(this%RootMyco1stSinkC_rpvr))deallocate(this%RootMyco1stSinkC_rpvr);nullify(this%RootMyco1stSinkC_rpvr)
  if(associated(this%RootCO2EmisPot_pvr))deallocate(this%RootCO2EmisPot_pvr);nullify(this%RootCO2EmisPot_pvr)
  if(associated(this%RootCO2Autor_pvr))deallocate(this%RootCO2Autor_pvr);nullify(this%RootCO2Autor_pvr)
  if(associated(this%TurgEff4CanopyResp_pft))deallocate(this%TurgEff4CanopyResp_pft);nullify(this%TurgEff4CanopyResp_pft)
  if(associated(this%RootCO2AutorX_pvr))deallocate(this%RootCO2AutorX_pvr);nullify(this%RootCO2AutorX_pvr)
  if(associated(this%PlantExudElm_CumYr_pft))deallocate(this%PlantExudElm_CumYr_pft);nullify(this%PlantExudElm_CumYr_pft)
  if(associated(this%trcg_root_vr))deallocate(this%trcg_root_vr);nullify(this%trcg_root_vr)
  if(associated(this%trcg_air2root_flx_vr))deallocate(this%trcg_air2root_flx_vr);nullify(this%trcg_air2root_flx_vr)
  if(associated(this%CO2FixCL_pft))deallocate(this%CO2FixCL_pft);nullify(this%CO2FixCL_pft)
  if(associated(this%RootNutUptakeN_pft))deallocate(this%RootNutUptakeN_pft);nullify(this%RootNutUptakeN_pft)
  if(associated(this%RootNutUptakeP_pft))deallocate(this%RootNutUptakeP_pft);nullify(this%RootNutUptakeP_pft)
  if(associated(this%CO2FixLL_pft))deallocate(this%CO2FixLL_pft);nullify(this%CO2FixLL_pft)
  if(associated(this%RootUptk_N_CumYr_pft))deallocate(this%RootUptk_N_CumYr_pft);nullify(this%RootUptk_N_CumYr_pft)
  if(associated(this%RootUptk_P_CumYr_pft))deallocate(this%RootUptk_P_CumYr_pft);nullify(this%RootUptk_P_CumYr_pft)
  if(associated(this%RootCO2Ar2Soil_pvr))deallocate(this%RootCO2Ar2Soil_pvr);nullify(this%RootCO2Ar2Soil_pvr)
  if(associated(this%RootCO2Ar2RootX_pvr))deallocate(this%RootCO2Ar2RootX_pvr);nullify(this%RootCO2Ar2RootX_pvr)
  if(associated(this%RootCO2Ar2RootX_rpvr))deallocate(this%RootCO2Ar2RootX_rpvr);nullify(this%RootCO2Ar2RootX_rpvr)
  if(associated(this%trcs_deadroot2soil_pvr))deallocate(this%trcs_deadroot2soil_pvr);nullify(this%trcs_deadroot2soil_pvr)
  if(associated(this%GroSrcRootStress_pvr))deallocate(this%GroSrcRootStress_pvr);nullify(this%GroSrcRootStress_pvr)
  call PrintInfo('end '//subname)
  end subroutine plt_rootbgc_destroy
end module PlantRootBGCAPIData
