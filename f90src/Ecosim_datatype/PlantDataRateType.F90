module PlantDataRateType
  ! Description:
  ! Store plant exchange rates, accumulated diagnostics and restart transfer state.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__


  real(r8),target,allocatable :: CanopyGrosRCO2_pft(:,:,:)                       !canopy autotrophic respiraiton, [gC d-2 h-1]
  real(r8),target,allocatable ::  Eco_NEE_col(:,:)                               !total canopy net CO2 exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  NH3Dep2Can_pft(:,:,:)                          !canopy NH3 flux, [g d-2 h-1]
  real(r8),target,allocatable ::  RNFixCO2_pft(:,:,:)                            !CO2 respired by nodules, [gC d-2]    
  real(r8),target,allocatable ::  NodulInfectElms_pft(:,:,:,:)                   !pft nodule infection [g d-2 h-1]
  real(r8),target,allocatable ::  NH3Emis_CumYr_pft(:,:,:)                       !total canopy NH3 flux, [g d-2 ]
  real(r8),target,allocatable ::  SurfLitrfalStrutElms_CumYr_pft(:,:,:,:)        !total surface LitrFall element, [g d-2]
  real(r8),target,allocatable ::  Soil2RootMycoExudE_pvr(:,:,:,:,:,:,:)           !root uptake (+ve) - exudation (-ve) of DOC, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNutUptake_pvr(:,:,:,:,:,:)                 !root uptake of NH4 non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RootN2Fix_pvr(:,:,:,:)                         !root N2 fixation, [gN d-2 h-1]
  real(r8),target,allocatable ::  RootN2Fix_vr(:,:,:)                            !vertical profile of root N2 fixation, [gN d-2 h-1]
  real(r8),target,allocatable ::  RootN2Fix_col(:,:)                             !root N2 fixation, [gN d-2 h-1]
  real(r8),target,allocatable ::  RootUptkSoiSol_pvr(:,:,:,:,:,:)                 !aqueous H2 flux from roots to soil water, [g d-2 h-1]
  real(r8),target,allocatable ::  RootH2PO4DmndSoil_pvr(:,:,:,:,:)               !root uptake of H2PO4 non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RootH2PO4DmndBand_pvr(:,:,:,:,:)               !root uptake of H2PO4 band, [g d-2 h-1]
  real(r8),target,allocatable ::  RootH1PO4DmndSoil_pvr(:,:,:,:,:)               !HPO4 demand in non-band by each root population, [g d-2 h-1]
  real(r8),target,allocatable ::  RootH1PO4DmndBand_pvr(:,:,:,:,:)               !HPO4 demand in band by each root population, [g d-2 h-1]
  real(r8),target,allocatable ::  LeafElmntRemobFlx_brch(:,:,:,:,:)              !cached remobilizable leaf element mass, [g d-2]
  real(r8),target,allocatable :: LeafSenescInitialElms_brch(:,:,:,:,:) !initial senescing leaf C/N/P mass; -1 means unset, [g d-2]
  real(r8),target,allocatable ::  PetolShethChemElmRemobFlx_brch(:,:,:,:,:)         !cached remobilizable sheath element mass, [g d-2]
  real(r8),target,allocatable :: PetolSenescInitialElms_brch(:,:,:,:,:) !initial senescing sheath/petiole C/N/P mass; -1 means unset, [g d-2]
  real(r8),target,allocatable ::  GrossCO2Fix_pft(:,:,:)                         !total gross CO2 fixation, [g d-2 h-1]
  real(r8),target,allocatable :: GPP_DaySum_col(:,:) !current day's accumulated GPP, [g C d-2]
  real(r8),target,allocatable :: GPP_PrevDay_col(:,:) !previous completed day's GPP; zero until available, [g C d-2]
  integer,target,allocatable :: GPP_DayHours_col(:,:) !hours accumulated; partial startup days are not promoted, [-]
  real(r8),target,allocatable ::  GrossCO2Fix_CumYr_pft(:,:,:)                   !cumulative total gross CO2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  LitrfallElms_pft(:,:,:,:)                  !total plant element LitrFall , [g d-2 ]
  real(r8),target,allocatable ::  PlantN2Fix_CumYr_pft(:,:,:)                    !total plant N2 fixation, [g d-2 ]
  real(r8),target,allocatable ::  GrossRespC_CumYr_pft(:,:,:)                    !cumulative total plant respiration, [g d-2 ]
  real(r8),target,allocatable ::  GrossResp_pft(:,:,:)                           !pft level plant respiraiton, [g d-2 h-1]
  real(r8),target,allocatable ::  PlantElmBalCum_pft(:,:,:,:)                     !plant element balance, [g d-2]
  real(r8),target,allocatable ::  LitrfalStrutElms_CumYr_pft(:,:,:,:)            !plant element LitrFall, [g d-2 h-1]
  real(r8),target,allocatable ::  LitrfallElms_pvr(:,:,:,:,:,:,:)            !plant LitrFall element, [g d-2 h-1]
  real(r8),target,allocatable ::  NetPrimProduct_pft(:,:,:)                      !Hourly net primary productivity, [g d-2 h-1]
  real(r8),target,allocatable :: PTSHTR_pft(:,:,:)                               !root-shoot coupling conductance, [h-1]
  real(r8),target,allocatable ::  cumNPP_pft(:,:,:)                              !cumulative net primary productivity, [g d-2]
  real(r8),target,allocatable ::  ETCanopy_CumYr_pft(:,:,:)                      !total transpiration  <0 into atmosphere, [m d-2]
  real(r8),target,allocatable ::  CanopyRespC_CumYr_pft(:,:,:)                   !total autotrophic respiration, [g d-2 ]
  real(r8),target,allocatable ::  EcoHavstElmnt_CumYr_pft(:,:,:,:)               !plant element harvest, [g d-2 ]
  real(r8),target,allocatable ::  EcoHavstElmntCum_pft(:,:,:,:)                  !total plant harvest, [g d-2 ]
  real(r8),target,allocatable ::  CO2ByFire_CumYr_pft(:,:,:)                     !plant CO2 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  FireLossE_pft(:,:,:,:)                         !plant element lost by fire, [g d-2 h-1]
  real(r8),target,allocatable ::  CH4ByFire_CumYr_pft(:,:,:)                     !plant CH4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  O2ByFire_CumYr_pft(:,:,:)                      !plant O2 uptake from fire, [g d-2 ]
  real(r8),target,allocatable ::  NH3byFire_CumYr_pft(:,:,:)                     !plant NH3 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  N2ObyFire_CumYr_pft(:,:,:)                     !plant N2O emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  PO4byFire_CumYr_pft(:,:,:)                     !plant PO4 emission from fire, [g d-2 ]
  real(r8),target,allocatable ::  RootO2Dmnd4Resp_pvr(:,:,:,:,:)                 !root  O2 demand from respiration, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_air2root_flx_pvr(:,:,:,:,:,:)             !gaseous tracer flux through roots, [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_Root_gas2aqu_flx_vr(:,:,:,:,:,:)          !dissolution (+ve) - volatilization (-ve) gas flux in roots, [g d-2 h-1]
  real(r8),target,allocatable ::  RCO2Emis2Root_rpvr(:,:,:,:,:)                     !aqueous CO2 flux from roots to root water , [g d-2 h-1]
  real(r8),target,allocatable ::  RootO2Uptk_pvr(:,:,:,:,:)                      !aqueous O2 flux from roots to root water , [g d-2 h-1]
  real(r8),target,allocatable ::  RootRespPotent_pvr(:,:,:,:,:)                  !root respiration unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCO2Autor_pvr(:,:,:,:,:)                    !root respiration constrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  Soil2RootMycoExudE_pft(:,:,:,:)                  !total root uptake (+ve) - exudation (-ve) of dissovled element, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNH4Uptake_pft(:,:,:)                       !total root uptake of NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNO3Uptake_pft(:,:,:)                       !total root uptake of NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RootH2PO4Uptake_pft(:,:,:)                     !total root uptake of PO4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootHPO4Uptake_pft(:,:,:)                      !total root uptake of HPO4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootN2Fix_pft(:,:,:)                           !total root N2 fixation, [g d-2 h-1]
  real(r8),target,allocatable ::  RootGasLossDisturb_pft(:,:,:,:)                !gas flux from root disturbance (<0 into atmosphere), [g d-2 h-1]
  real(r8),target,allocatable ::  RootOUlmNutUptake_pvr(:,:,:,:,:,:)             !root uptake of NH4 non-band unconstrained by O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCUlmNutUptake_pvr(:,:,:,:,:,:)             !root uptake of NH4 non-band unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCO2EmisPot_pvr(:,:,:,:,:)                  !root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
  real(r8),target,allocatable ::  VmaxNH4Root_pvr(:,:,:,:,:)                     !maximum NH4 uptake rate, [gN h-1 (gC root)-1]
  real(r8),target,allocatable ::  VmaxNO3Root_pvr(:,:,:,:,:)                     !maximum NO3 uptake rate, [gN h-1 (gC root)-1]
  real(r8),target,allocatable ::  RootNH4DmndSoil_pvr(:,:,:,:,:)                 !root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNO3DmndSoil_pvr(:,:,:,:,:)                 !root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNH4DmndBand_pvr(:,:,:,:,:)                 !root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  RootNO3DmndBand_pvr(:,:,:,:,:)                 !root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
  real(r8),target,allocatable ::  NH3Dep2Can_brch(:,:,:,:)                       !gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
  real(r8),target,allocatable ::  RAutoRootO2Limter_rpvr(:,:,:,:,:)              !O2 constraint to root respiration, [-]
  real(r8),target,allocatable ::  PlantRootSoilElmNetX_pft(:,:,:,:)              !net root element uptake (+ve) - exudation (-ve), [g d-2 h-1]
  real(r8),target,allocatable ::  PlantExudElm_CumYr_pft(:,:,:,:)                !total net root element uptake (+ve) - exudation (-ve), [g d-2 ]
  real(r8),target,allocatable ::  RootUptk_N_CumYr_pft(:,:,:)                    !pft cumulative N uptake, [gN d-2]
  real(r8),target,allocatable ::  RootUptk_P_CumYr_pft(:,:,:)                    !pft cumulative P uptake, [gP d-2]
  real(r8),target,allocatable ::  RootUptk_Nmin_cumYr_pft(:,:,:)                 !pft cumulative mineral N uptake, [gN d-2]
  real(r8),target,allocatable ::  RootUptk_Pmin_cumYr_pft(:,:,:)                 !pft cumulative mineral P uptake, [gP d-2]
  real(r8),target,allocatable ::  TPlantRootH2OUptake_col(:,:)                   !total root H2O uptake, [m3 d-2 h-1]
  ! Accepted hourly exchanges, signed positive from plant to soil. ATS reads
  ! these unchanged; the native soil solver apportions them over its substeps.
  real(r8),target,allocatable ::  TWaterPlantRoot2Soil_vr(:,:,:)                 !accepted root water profile [m3/cell/h]
  real(r8),target,allocatable ::  TWaterPlantRoot2SoilPrev_vr(:,:,:)             !root water actually applied to soil this hour [m3/cell/h]
  real(r8),target,allocatable ::  THeatPlantRoot2SoilPrev_vr(:,:,:)              !matching applied root donor enthalpy [MJ/cell/h]
  real(r8),target,allocatable ::  THeatLossRoot2Soil_vr(:,:,:)                   !matching accepted donor enthalpy profile [MJ/cell/h]
  real(r8),target,allocatable ::  THeatRootRelease_col(:,:)                      !total root heat relase, [MJ d-2 h-1]
  real(r8),target,allocatable ::  trcg_air2root_flx_vr(:,:,:,:)                  !total internal root gas flux , [g d-2 h-1]
  real(r8),target,allocatable ::  trcg_root_vr(:,:,:,:)                          !total root internal gas flux, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_Soil2plant_uptake_vr(:,:,:,:)                  !total root-soil solute flux, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_Soil2plant_uptake_col(:,:,:)                   !total root-soil solute flux, [g d-2 h-1]
  real(r8),target,allocatable ::  trcs_Soil2plant_uptakep_col(:,:,:)                  !total root-soil solute flux, [g d-2 h-1]  
  real(r8),target,allocatable ::  tRootMycoExud2Soil_vr(:,:,:,:,:)               !total root element exchange, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCO2Emis2Root_vr(:,:,:)                     !total root CO2 flux into roots, [g d-2 h-1]
  real(r8),target,allocatable ::  RootCO2Emis2Root_col(:,:)                     !total root CO2 flux into roots, [g d-2 h-1]
  real(r8),target,allocatable ::  RUptkRootO2_vr(:,:,:)                          !total root internal O2 flux taken away from root O2, [g d-2 h-1]
  real(r8),target,allocatable ::  RootO2_TotSink_vr(:,:,:)                          !root O2 consumption for autotrophic respiraiton, [gO d-2 h-1]
  real(r8),target,allocatable ::  RootO2_TotSink_col(:,:)                        !integrated root O2 uptake from soil and atmosphere for autotrophic respiraiton, [gO d-2 h-1]  
  real(r8),target,allocatable ::  RUptkRootO2_col(:,:)                           !total root internal O2 flux take away from root O2, [g d-2 h-1]  
  real(r8),target,allocatable ::  totRootLenDens_vr(:,:,:)                       !total root length density, [m m-3]
  real(r8),target,allocatable ::  RCH4EcoDmndPrev_vr(:,:,:)                      !total root + microbial CH4 uptake, [gCd-2 h-1]
  real(r8),target,allocatable ::  RCH4EcoDmnd_vr(:,:,:)                          !total root + microbial CH4 uptake, [gCd-2 h-1]
  real(r8),target,allocatable ::  REcoO2DmndResp_vr(:,:,:)                       !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RO2EcoDmndPrev_vr(:,:,:)                       !total root + microbial O2 uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoNH4DmndSoil_vr(:,:,:)                      !total root + microbial NH4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH4EcoDmndSoilPrev_vr(:,:,:)                  !total root + microbial NH4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoNO3DmndSoil_vr(:,:,:)                      !total root + microbial NO3 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO3EcoDmndSoilPrev_vr(:,:,:)                  !total root + microbial NO3 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2EcoUptkSoil_vr(:,:,:)                      !total root + microbial NO2 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2EcoUptkSoilPrev_vr(:,:,:)                  !total root + microbial NO2 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoH2PO4DmndSoil_vr(:,:,:)                    !total root + microbial PO4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH2PO4EcoDmndSoilPrev_vr(:,:,:)                !total root + microbial PO4 uptake non-band, [g d-2 h-1]
  real(r8),target,allocatable ::  RN2OEcoUptkSoil_vr(:,:,:)                      !total root + microbial N2O uptake , [g d-2 h-1]
  real(r8),target,allocatable ::  RN2OEcoUptkSoilPrev_vr(:,:,:)                  !total root + microbial N2O uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoNH4DmndBand_vr(:,:,:)                      !total root + microbial NH4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNH4EcoDmndBandPrev_vr(:,:,:)                  !total root + microbial NH4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoNO3DmndBand_vr(:,:,:)                      !total root + microbial NO3 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO3EcoDmndBandPrev_vr(:,:,:)                  !total root + microbial NO3 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2EcoUptkBand_vr(:,:,:)                      !total root + microbial NO2 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RNO2EcoUptkBandPrev_vr(:,:,:)                  !total root + microbial NO2 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  REcoH2PO4DmndBand_vr(:,:,:)                    !total root + microbial PO4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RH2PO4EcoDmndBandPrev_vr(:,:,:)                !total root + microbial PO4 uptake band, [g d-2 h-1]
  real(r8),target,allocatable ::  RDOMEcoDmndK_vr(:,:,:,:)                       !total root + microbial DOC uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RDOMEcoDmndPrev_vr(:,:,:,:)                    !total root + microbial DOC uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RAcetateEcoDmndK_vr(:,:,:,:)                   !total root + microbial acetate uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  RAcetateEcoDmndPrev_vr(:,:,:,:)                !total root + microbial acetate uptake, [g d-2 h-1]
  real(r8),target,allocatable ::  TRootH2Flx_col(:,:)                            !total root H2 flux, [g d-2]
  real(r8),target,allocatable ::  RootCO2Autor_vr(:,:,:)                         !root autotrophic respiraiton, [gC d-3 hr-1]
  real(r8),target,allocatable ::  RootCO2Autor_col(:,:)                          !current time step root autotrophic respiraiton, [gC d-2 h-1]  
  real(r8),target,allocatable ::  RootCO2AutorPrev_col(:,:)                      !previous time step root autotrophic respiraiton, [gC d-2 h-1]  
  real(r8),target,allocatable ::  fRootGrowPSISense_pvr(:,:,:,:,:)               !moisture dependence scalar for root growth, [-]
  real(r8),target,allocatable ::  RootCO2Ar2Soil_vr(:,:,:)                       !autotrophic root respiration released to soil, [gC d-2 h-1]
  real(r8),target,allocatable ::  RootCO2Ar2Soil_col(:,:)                        !total autotrophic root respiraiton released to soil, [gC d-2 h-1]
  real(r8),target,allocatable ::  RootCO2Ar2Root_vr(:,:,:)                       !autotrophic root respiration released to root, [gC d-2 h-1]
  real(r8),target,allocatable ::  RootCO2Ar2Root_col(:,:)                        !total autotrophic root respiration released to root, [gC d-2 h-1]
  real(r8),target,allocatable ::  trcs_deadroot2soil_vr(:,:,:,:)                 !gases released to soil due to dying roots, ([g d-2 h-1])
  real(r8),target,allocatable ::  trcs_deadroot2soil_col(:,:,:)                  !gas released to soil due to dying roots, [g d-2 h-1]
  real(r8),target,allocatable ::  TotEndVegE_pft(:,:,:,:)                        !total vegetation biomass at the end of time step, [g d-2]
  real(r8),target,allocatable ::  TotBegVegE_pft(:,:,:,:)                        !total vegetation biomass at the beginning of time step, [g d-2]
  real(r8),target,allocatable ::  CanopyN2Fix_pft(:,:,:)                         !total canopy N2 fixation, [g d-2 h-1]          
  real(r8),target,allocatable ::  ShootRootXferElm_pft(:,:,:,:)                  !shoot-root nonstructural element transfer, [ g d-2 h-1]  
  private :: InitAllocate
  contains

!----------------------------------------------------------------------
  subroutine InitPlantRates(NumOfPlantLitrCmplxs,jroots)
  implicit none
  integer, intent(in) :: NumOfPlantLitrCmplxs
  integer, intent(in) :: jroots

  call InitAllocate(NumOfPlantLitrCmplxs,jroots)

  end subroutine InitPlantRates
!----------------------------------------------------------------------

  subroutine InitAllocate(NumOfPlantLitrCmplxs,jroots)

  implicit none
  integer, intent(in) :: NumOfPlantLitrCmplxs
  integer, intent(in) :: jroots    !number of root types, root,mycos

  allocate(TotEndVegE_pft(NumPlantChemElms,JP,JY,JX)); TotEndVegE_pft=0._r8
  allocate(TotBegVegE_pft(NumPlantChemElms,JP,JY,JX)); TotBegVegE_pft=0._r8  
  allocate(trcs_deadroot2soil_vr(idg_beg:idg_NH3,JZ,JY,JX)); trcs_deadroot2soil_vr=0._r8
  allocate(trcs_deadroot2soil_col(idg_beg:idg_NH3,JY,JX)); trcs_deadroot2soil_col=0._r8
  allocate(TPlantRootH2OUptake_col(JY,JX)); TPlantRootH2OUptake_col=0._r8
  allocate(CanopyGrosRCO2_pft(JP,JY,JX)); CanopyGrosRCO2_pft=0._r8
  allocate(Eco_NEE_col(JY,JX));       Eco_NEE_col=0._r8
  allocate(NH3Dep2Can_pft(JP,JY,JX));    NH3Dep2Can_pft=0._r8
  allocate(RNFixCO2_pft(JP,JY,JX)); RNFixCO2_pft=0._r8
  allocate(NH3Emis_CumYr_pft(JP,JY,JX));    NH3Emis_CumYr_pft=0._r8  
  allocate(NodulInfectElms_pft(NumPlantChemElms,JP,JY,JX));NodulInfectElms_pft=0._r8
  allocate(SurfLitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP,JY,JX));    SurfLitrfalStrutElms_CumYr_pft=0._r8
  allocate(Soil2RootMycoExudE_pvr(NumPlantChemElms,jroots,1:jcplx,JZ,JP,JY,JX));Soil2RootMycoExudE_pvr=0._r8
  allocate(RootNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ,JP,JY,JX));RootNutUptake_pvr=0._r8
  allocate(RootN2Fix_pvr(JZ,JP,JY,JX)); RootN2Fix_pvr=0._r8
  allocate(RootN2Fix_vr(JZ,JY,JX)); RootN2Fix_vr=0._r8
  allocate(RootN2Fix_col(JY,JX)); RootN2Fix_col=0._r8  
  allocate(RootH2PO4DmndSoil_pvr(jroots,JZ,JP,JY,JX));RootH2PO4DmndSoil_pvr=0._r8
  allocate(RootH2PO4DmndBand_pvr(jroots,JZ,JP,JY,JX));RootH2PO4DmndBand_pvr=0._r8
  allocate(RootH1PO4DmndSoil_pvr(jroots,JZ,JP,JY,JX));RootH1PO4DmndSoil_pvr=0._r8
  allocate(RootH1PO4DmndBand_pvr(jroots,JZ,JP,JY,JX));RootH1PO4DmndBand_pvr=0._r8
  allocate(LeafElmntRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); LeafElmntRemobFlx_brch=0._r8
  allocate(LeafSenescInitialElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); LeafSenescInitialElms_brch=-1._r8
  allocate(PetolShethChemElmRemobFlx_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); PetolShethChemElmRemobFlx_brch=0._r8
  allocate(PetolSenescInitialElms_brch(NumPlantChemElms,MaxNumBranches,JP,JY,JX)); PetolSenescInitialElms_brch=-1._r8
  allocate(GrossCO2Fix_pft(JP,JY,JX));    GrossCO2Fix_pft=0._r8
  allocate(GPP_DaySum_col(JY,JX)); GPP_DaySum_col=0._r8
  allocate(GPP_PrevDay_col(JY,JX)); GPP_PrevDay_col=0._r8
  allocate(GPP_DayHours_col(JY,JX)); GPP_DayHours_col=0
  allocate(GrossCO2Fix_CumYr_pft(JP,JY,JX)); GrossCO2Fix_CumYr_pft=0._r8
  allocate(LitrfalStrutElms_CumYr_pft(NumPlantChemElms,JP,JY,JX));    LitrfalStrutElms_CumYr_pft=0._r8
  allocate(PlantN2Fix_CumYr_pft(JP,JY,JX));   PlantN2Fix_CumYr_pft=0._r8
  allocate(GrossResp_pft(JP,JY,JX));    GrossResp_pft=0._r8
  allocate(GrossRespC_CumYr_pft(JP,JY,JX)); GrossRespC_CumYr_pft=0._r8
  allocate(PlantElmBalCum_pft(NumPlantChemElms,JP,JY,JX));     PlantElmBalCum_pft=0._r8
  allocate(LitrfallElms_pft(NumPlantChemElms,JP,JY,JX));    LitrfallElms_pft=0._r8
  allocate(LitrfallElms_pvr(NumPlantChemElms,jsken,1:NumOfPlantLitrCmplxs,0:JZ,JP,JY,JX));LitrfallElms_pvr=0._r8
  allocate(NetPrimProduct_pft(JP,JY,JX));     NetPrimProduct_pft=0._r8
  allocate(cumNPP_pft(JP,JY,JX)); cumNPP_pft=0._r8
  allocate(PTSHTR_pft(JP,JY,JX)); PTSHTR_pft=0._r8
  allocate(ETCanopy_CumYr_pft(JP,JY,JX));    ETCanopy_CumYr_pft=0._r8
  allocate(CanopyRespC_CumYr_pft(JP,JY,JX));    CanopyRespC_CumYr_pft=0._r8
  allocate(EcoHavstElmnt_CumYr_pft(NumPlantChemElms,JP,JY,JX));    EcoHavstElmnt_CumYr_pft=0._r8
  allocate(EcoHavstElmntCum_pft(NumPlantChemElms,JP,JY,JX));   EcoHavstElmntCum_pft=0._r8
  allocate(CO2ByFire_CumYr_pft(JP,JY,JX));    CO2ByFire_CumYr_pft=0._r8
  allocate(FireLossE_pft(NumPlantChemElms,JP,JY,JX)); FireLossE_pft=0._r8
  allocate(CH4ByFire_CumYr_pft(JP,JY,JX));    CH4ByFire_CumYr_pft=0._r8
  allocate(O2ByFire_CumYr_pft(JP,JY,JX));    O2ByFire_CumYr_pft=0._r8
  allocate(NH3byFire_CumYr_pft(JP,JY,JX));    NH3byFire_CumYr_pft=0._r8
  allocate(N2ObyFire_CumYr_pft(JP,JY,JX));    N2ObyFire_CumYr_pft=0._r8
  allocate(PO4byFire_CumYr_pft(JP,JY,JX));    PO4byFire_CumYr_pft=0._r8
  allocate(RootO2Dmnd4Resp_pvr(jroots,JZ,JP,JY,JX));RootO2Dmnd4Resp_pvr=0._r8
  allocate(trcg_air2root_flx_pvr(idg_beg:idg_NH3,2,JZ,JP,JY,JX));trcg_air2root_flx_pvr=0._r8
  allocate(trcg_Root_gas2aqu_flx_vr(idg_beg:idg_NH3,2,JZ,JP,JY,JX));trcg_Root_gas2aqu_flx_vr=0._r8
  allocate(RootUptkSoiSol_pvr(idg_beg:idg_end,jroots,JZ,JP,JY,JX));RootUptkSoiSol_pvr=0._r8
  allocate(RCO2Emis2Root_rpvr(jroots,JZ,JP,JY,JX));RCO2Emis2Root_rpvr=0._r8
  allocate(RootO2Uptk_pvr(jroots,JZ,JP,JY,JX));RootO2Uptk_pvr=0._r8
  allocate(RootRespPotent_pvr(jroots,JZ,JP,JY,JX));RootRespPotent_pvr=0._r8
  allocate(RootCO2Autor_pvr(jroots,JZ,JP,JY,JX));RootCO2Autor_pvr=0._r8
  allocate(Soil2RootMycoExudE_pft(1:NumPlantChemElms,JP,JY,JX));    Soil2RootMycoExudE_pft=0._r8
  allocate(RootNH4Uptake_pft(JP,JY,JX));    RootNH4Uptake_pft=0._r8
  allocate(RootNO3Uptake_pft(JP,JY,JX));    RootNO3Uptake_pft=0._r8
  allocate(RootH2PO4Uptake_pft(JP,JY,JX));    RootH2PO4Uptake_pft=0._r8
  allocate(RootHPO4Uptake_pft(JP,JY,JX));    RootHPO4Uptake_pft=0._r8
  allocate(RootN2Fix_pft(JP,JY,JX));     RootN2Fix_pft=0._r8
  allocate(CanopyN2Fix_pft(JP,JY,JX)); CanopyN2Fix_pft=0._r8
  allocate(RootGasLossDisturb_pft(idg_beg:idg_NH3,JP,JY,JX)); RootGasLossDisturb_pft=0._r8
  allocate(RootCUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ,JP,JY,JX));RootCUlmNutUptake_pvr=0._r8
  allocate(RootOUlmNutUptake_pvr(ids_nutb_beg+1:ids_nuts_end,jroots,JZ,JP,JY,JX));RootOUlmNutUptake_pvr=0._r8
  allocate(RootCO2EmisPot_pvr(jroots,JZ,JP,JY,JX));RootCO2EmisPot_pvr=0._r8
  allocate(RootNH4DmndSoil_pvr(jroots,JZ,JP,JY,JX));RootNH4DmndSoil_pvr=0._r8
  allocate(VmaxNH4Root_pvr(jroots,JZ,JP,JY,JX)); VmaxNH4Root_pvr=0._r8
  allocate(VmaxNO3Root_pvr(jroots,JZ,JP,JY,JX)); VmaxNO3Root_pvr=0._r8
  allocate(RootNO3DmndSoil_pvr(jroots,JZ,JP,JY,JX));RootNO3DmndSoil_pvr=0._r8
  allocate(RootNH4DmndBand_pvr(jroots,JZ,JP,JY,JX));RootNH4DmndBand_pvr=0._r8
  allocate(RootNO3DmndBand_pvr(jroots,JZ,JP,JY,JX));RootNO3DmndBand_pvr=0._r8
  allocate(fRootGrowPSISense_pvr(jroots,JZ,JP,JY,JX)); fRootGrowPSISense_pvr=0._r8
  allocate(NH3Dep2Can_brch(MaxNumBranches,JP,JY,JX)); NH3Dep2Can_brch=0._r8
  allocate(RAutoRootO2Limter_rpvr(jroots,JZ,JP,JY,JX)); RAutoRootO2Limter_rpvr=0._r8
  allocate(PlantRootSoilElmNetX_pft(NumPlantChemElms,JP,JY,JX));   PlantRootSoilElmNetX_pft=0._r8
  allocate(PlantExudElm_CumYr_pft(NumPlantChemElms,JP,JY,JX));   PlantExudElm_CumYr_pft=0._r8
  allocate(RootUptk_N_CumYr_pft(JP,JY,JX)); RootUptk_N_CumYr_pft=0._r8
  allocate(RootUptk_P_CumYr_pft(JP,JY,JX)); RootUptk_P_CumYr_pft=0._r8
  allocate(RootUptk_Nmin_cumYr_pft(JP,JY,JX)); RootUptk_Nmin_cumYr_pft=0._r8
  allocate(RootUptk_Pmin_cumYr_pft(JP,JY,JX)); RootUptk_Pmin_cumYr_pft=0._r8
  allocate(TWaterPlantRoot2Soil_vr(JZ,JY,JX)); TWaterPlantRoot2Soil_vr=0._r8
  allocate(TWaterPlantRoot2SoilPrev_vr(JZ,JY,JX)); TWaterPlantRoot2SoilPrev_vr=0._r8  
  allocate(THeatPlantRoot2SoilPrev_vr(JZ,JY,JX)); THeatPlantRoot2SoilPrev_vr=0._r8
  allocate(THeatLossRoot2Soil_vr(0:JZ,JY,JX));  THeatLossRoot2Soil_vr=0._r8
  allocate(THeatRootRelease_col(JY,JX)); THeatRootRelease_col=0._r8
  allocate(trcg_air2root_flx_vr(idg_beg:idg_NH3,JZ,JY,JX));   trcg_air2root_flx_vr=0._r8
  allocate(trcg_root_vr(idg_beg:idg_NH3,JZ,JY,JX));   trcg_root_vr=0._r8
  allocate(trcs_Soil2plant_uptake_vr(ids_beg:ids_end,JZ,JY,JX));    trcs_Soil2plant_uptake_vr=0._r8
  allocate(trcs_Soil2plant_uptake_col(ids_beg:ids_end,JY,JX));    trcs_Soil2plant_uptake_col=0._r8
  allocate(trcs_Soil2plant_uptakep_col(ids_beg:ids_end,JY,JX));    trcs_Soil2plant_uptakep_col=0._r8  
  allocate(tRootMycoExud2Soil_vr(NumPlantChemElms,1:jcplx,JZ,JY,JX));tRootMycoExud2Soil_vr=0._r8
  allocate(RootCO2Emis2Root_vr(JZ,JY,JX));    RootCO2Emis2Root_vr=0._r8
  allocate(RootCO2Emis2Root_col(JY,JX));    RootCO2Emis2Root_col=0._r8  
  allocate(RUptkRootO2_vr(JZ,JY,JX));   RUptkRootO2_vr=0._r8
  allocate(RootO2_TotSink_vr(JZ,JY,JX));  RootO2_TotSink_vr=0._r8
  allocate(RootO2_TotSink_col(JY,JX));  RootO2_TotSink_col=0._r8  
  allocate(RUptkRootO2_col(JY,JX));   RUptkRootO2_col=0._r8
  allocate(totRootLenDens_vr(JZ,JY,JX));    totRootLenDens_vr=0._r8
  allocate(REcoO2DmndResp_vr(0:JZ,JY,JX));  REcoO2DmndResp_vr=0._r8
  allocate(RO2EcoDmndPrev_vr(0:JZ,JY,JX));  RO2EcoDmndPrev_vr=0._r8
  allocate(RCH4EcoDmndPrev_vr(0:JZ,JY,JX)); RCH4EcoDmndPrev_vr=0._r8
  allocate(RCH4EcoDmnd_vr(0:JZ,JY,JX)); RCH4EcoDmnd_vr=0._r8
  allocate(REcoNH4DmndSoil_vr(0:JZ,JY,JX));  REcoNH4DmndSoil_vr=0._r8
  allocate(RNH4EcoDmndSoilPrev_vr(0:JZ,JY,JX));  RNH4EcoDmndSoilPrev_vr=0._r8
  allocate(REcoNO3DmndSoil_vr(0:JZ,JY,JX));  REcoNO3DmndSoil_vr=0._r8
  allocate(RNO3EcoDmndSoilPrev_vr(0:JZ,JY,JX));  RNO3EcoDmndSoilPrev_vr=0._r8
  allocate(RNO2EcoUptkSoil_vr(0:JZ,JY,JX));  RNO2EcoUptkSoil_vr=0._r8
  allocate(RNO2EcoUptkSoilPrev_vr(0:JZ,JY,JX));  RNO2EcoUptkSoilPrev_vr=0._r8
  allocate(REcoH2PO4DmndSoil_vr(0:JZ,JY,JX));  REcoH2PO4DmndSoil_vr=0._r8
  allocate(RH2PO4EcoDmndSoilPrev_vr(0:JZ,JY,JX));  RH2PO4EcoDmndSoilPrev_vr=0._r8
  allocate(RN2OEcoUptkSoil_vr(0:JZ,JY,JX));  RN2OEcoUptkSoil_vr=0._r8
  allocate(RN2OEcoUptkSoilPrev_vr(0:JZ,JY,JX));  RN2OEcoUptkSoilPrev_vr=0._r8
  allocate(REcoNH4DmndBand_vr(0:JZ,JY,JX));  REcoNH4DmndBand_vr=0._r8
  allocate(RNH4EcoDmndBandPrev_vr(0:JZ,JY,JX));  RNH4EcoDmndBandPrev_vr=0._r8
  allocate(REcoNO3DmndBand_vr(0:JZ,JY,JX));  REcoNO3DmndBand_vr=0._r8
  allocate(RNO3EcoDmndBandPrev_vr(0:JZ,JY,JX));  RNO3EcoDmndBandPrev_vr=0._r8
  allocate(RNO2EcoUptkBand_vr(0:JZ,JY,JX));  RNO2EcoUptkBand_vr=0._r8
  allocate(RNO2EcoUptkBandPrev_vr(0:JZ,JY,JX));  RNO2EcoUptkBandPrev_vr=0._r8
  allocate(REcoH2PO4DmndBand_vr(0:JZ,JY,JX));  REcoH2PO4DmndBand_vr=0._r8
  allocate(RH2PO4EcoDmndBandPrev_vr(0:JZ,JY,JX));  RH2PO4EcoDmndBandPrev_vr=0._r8
  allocate(RDOMEcoDmndK_vr(1:jcplx,0:JZ,JY,JX));RDOMEcoDmndK_vr=0._r8
  allocate(RDOMEcoDmndPrev_vr(1:jcplx,0:JZ,JY,JX));RDOMEcoDmndPrev_vr=0._r8
  allocate(RAcetateEcoDmndK_vr(1:jcplx,0:JZ,JY,JX));RAcetateEcoDmndK_vr=0._r8
  allocate(RAcetateEcoDmndPrev_vr(1:jcplx,0:JZ,JY,JX));RAcetateEcoDmndPrev_vr=0._r8
  allocate(TRootH2Flx_col(JY,JX));       TRootH2Flx_col=0._r8
  allocate(RootCO2Autor_vr(JZ,JY,JX));   RootCO2Autor_vr=0._r8
  allocate(RootCO2Ar2Soil_vr(JZ,JY,JX)); RootCO2Ar2Soil_vr=0._r8
  allocate(RootCO2Ar2Soil_col(JY,JX)); RootCO2Ar2Soil_col=0._r8
  allocate(RootCO2Ar2Root_col(JY,JX)); RootCO2Ar2Root_col=0._r8
  allocate(RootCO2Ar2Root_vr(JZ,JY,JX)); RootCO2Ar2Root_vr=0._r8
  allocate(RootCO2Autor_col(JY,JX));     RootCO2Autor_col=0._r8
  allocate(RootCO2AutorPrev_col(JY,JX)); RootCO2AutorPrev_col=0._r8
  allocate(ShootRootXferElm_pft(NumPlantChemElms,JZ,JY,JX)); ShootRootXferElm_pft=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructPlantRates
  use abortutils, only : destroy
  implicit none
  call destroy(PTSHTR_pft)
  call destroy(TotBegVegE_pft)
  call destroy(TotEndVegE_pft)
  call destroy(trcs_deadroot2soil_vr)
  call destroy(RootCO2Ar2Soil_vr)
  call destroy(RootCO2Ar2Soil_col)
  call destroy(RootCO2Ar2Root_col)
  call destroy(RootCO2Ar2Root_vr)
  call destroy(RootCO2Autor_col)
  call destroy(RootCO2AutorPrev_col)
  call destroy(TPlantRootH2OUptake_col)
  call destroy(CanopyGrosRCO2_pft)
  call destroy(Eco_NEE_col)
  call destroy(RNFixCO2_pft)
  call destroy(NH3Dep2Can_pft)
  call destroy(NH3Emis_CumYr_pft)
  call destroy(SurfLitrfalStrutElms_CumYr_pft)
  call destroy(Soil2RootMycoExudE_pvr)
  call destroy(RootNutUptake_pvr)
  call destroy(RootUptkSoiSol_pvr)
  call destroy(RootN2Fix_pvr)
  call destroy(RootN2Fix_vr)
  call destroy(RootN2Fix_col)
  call destroy(RootH2PO4DmndSoil_pvr)
  call destroy(RootH2PO4DmndBand_pvr)
  call destroy(RootH1PO4DmndSoil_pvr)
  call destroy(RootH1PO4DmndBand_pvr)
  call destroy(LeafElmntRemobFlx_brch)
  call destroy(LeafSenescInitialElms_brch)
  call destroy(PetolShethChemElmRemobFlx_brch)
  call destroy(PetolSenescInitialElms_brch)
  call destroy(GrossCO2Fix_pft)
  call destroy(GPP_DaySum_col)
  call destroy(GPP_PrevDay_col)
  call destroy(GPP_DayHours_col)
  call destroy(GrossCO2Fix_CumYr_pft)
  call destroy(LitrfalStrutElms_CumYr_pft)
  call destroy(PlantN2Fix_CumYr_pft)
  call destroy(GrossResp_pft)
  call destroy(GrossRespC_CumYr_pft)
  call destroy(PlantElmBalCum_pft)
  call destroy(LitrfallElms_pft)
  call destroy(LitrfallElms_pvr)
  call destroy(NetPrimProduct_pft)
  call destroy(cumNPP_pft)
  call destroy(ETCanopy_CumYr_pft)
  call destroy(CanopyRespC_CumYr_pft)
  call destroy(EcoHavstElmnt_CumYr_pft)
  call destroy(EcoHavstElmntCum_pft)
  call destroy(NodulInfectElms_pft)
  call destroy(CO2ByFire_CumYr_pft)
  call destroy(FireLossE_pft)
  call destroy(CH4ByFire_CumYr_pft)
  call destroy(O2ByFire_CumYr_pft)
  call destroy(NH3byFire_CumYr_pft)
  call destroy(N2ObyFire_CumYr_pft)
  call destroy(PO4byFire_CumYr_pft)
  call destroy(RootO2Dmnd4Resp_pvr)  
  call destroy(RCO2Emis2Root_rpvr)
  call destroy(RootO2Uptk_pvr)
  call destroy(RootRespPotent_pvr)
  call destroy(RootCO2Autor_pvr)
  call destroy(Soil2RootMycoExudE_pft)
  call destroy(RootNH4Uptake_pft)
  call destroy(RootNO3Uptake_pft)
  call destroy(RootH2PO4Uptake_pft)
  call destroy(RootHPO4Uptake_pft)
  call destroy(RootN2Fix_pft)
  CALL Destroy(CanopyN2Fix_pft)
  call destroy(RootOUlmNutUptake_pvr)
  call destroy(RootCUlmNutUptake_pvr)
  call destroy(RootCO2EmisPot_pvr)
  call destroy(RootNH4DmndSoil_pvr)
  call destroy(VmaxNH4Root_pvr)
  call destroy(VmaxNO3Root_pvr)
  call destroy(RootNO3DmndSoil_pvr)
  call destroy(RootNH4DmndBand_pvr)
  call destroy(RootNO3DmndBand_pvr)
  call destroy(fRootGrowPSISense_pvr)
  call destroy(NH3Dep2Can_brch)
  call destroy(RAutoRootO2Limter_rpvr)
  call destroy(PlantRootSoilElmNetX_pft)
  call destroy(PlantExudElm_CumYr_pft)
  call destroy(RootUptk_N_CumYr_pft)
  call destroy(RootUptk_P_CumYr_pft)
  call destroy(RootUptk_Nmin_cumYr_pft)
  call destroy(RootUptk_Pmin_cumYr_pft)
  call destroy(TWaterPlantRoot2Soil_vr)
  call destroy(TWaterPlantRoot2SoilPrev_vr)  
  call destroy(THeatPlantRoot2SoilPrev_vr)
  call destroy(THeatLossRoot2Soil_vr)
  call destroy(THeatRootRelease_col)
  call destroy(tRootMycoExud2Soil_vr)
  call destroy(RootCO2Emis2Root_vr)
  call destroy(RootCO2Emis2Root_col)
  call destroy(RUptkRootO2_vr)
  call destroy(RootO2_TotSink_vr)
  call destroy(RootO2_TotSink_col)
  call destroy(RUptkRootO2_col)
  call destroy(totRootLenDens_vr)
  call destroy(REcoO2DmndResp_vr)
  call destroy(RO2EcoDmndPrev_vr)
  call destroy(RCH4EcoDmndPrev_vr)
  call destroy(RCH4EcoDmnd_vr)
  call destroy(REcoNH4DmndSoil_vr)
  call destroy(RNH4EcoDmndSoilPrev_vr)
  call destroy(REcoNO3DmndSoil_vr)
  call destroy(RNO3EcoDmndSoilPrev_vr)
  call destroy(RNO2EcoUptkSoil_vr)
  call destroy(RNO2EcoUptkSoilPrev_vr)
  call destroy(REcoH2PO4DmndSoil_vr)
  call destroy(RH2PO4EcoDmndSoilPrev_vr)
  call destroy(RN2OEcoUptkSoil_vr)
  call destroy(RN2OEcoUptkSoilPrev_vr)
  call destroy(REcoNH4DmndBand_vr)
  call destroy(RNH4EcoDmndBandPrev_vr)
  call destroy(REcoNO3DmndBand_vr)
  call destroy(RNO3EcoDmndBandPrev_vr)
  call destroy(RNO2EcoUptkBand_vr)
  call destroy(RNO2EcoUptkBandPrev_vr)
  call destroy(REcoH2PO4DmndBand_vr)
  call destroy(RH2PO4EcoDmndBandPrev_vr)
  call destroy(RDOMEcoDmndK_vr)
  call destroy(ShootRootXferElm_pft)
  call destroy(RDOMEcoDmndPrev_vr)
  call destroy(RAcetateEcoDmndK_vr)
  call destroy(RAcetateEcoDmndPrev_vr)
  call destroy(TRootH2Flx_col)
  call destroy(RootGasLossDisturb_pft)
  call destroy(trcs_Soil2plant_uptake_vr)
  call destroy(trcs_Soil2plant_uptake_col)
  call destroy(trcs_Soil2plant_uptakep_col)  
  call destroy(RootCO2Autor_vr)
  end subroutine DestructPlantRates

end module PlantDataRateType
