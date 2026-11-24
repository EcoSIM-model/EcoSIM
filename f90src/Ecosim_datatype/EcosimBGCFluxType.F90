module EcosimBGCFluxType

!
! Ecosystm fluxes for C, N, and P budget
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  Eco_NetRad_col(:,:)                               !ecosystem net radiation, [MJ d-2 h-1]
  real(r8),target,allocatable ::  Eco_Heat_Latent_col(:,:)                          !ecosystem latent heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  Eco_Heat_Sens_col(:,:)                            !ecosystem sensible heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  Eco_Heat_GrndSurf_col(:,:)                        !ecosystem storage heat flux, [MJ d-2 h-1]
  real(r8),target,allocatable ::  Eco_GPP_CumYr_col(:,:)                            !cumulative ecosystem GPP, [g d-2]
  real(r8),target,allocatable ::  Eco_AutoR_CumYr_col(:,:)                          !cumulative ecosystem autotrophic respiration, [gC d-2]
  real(r8),target,allocatable ::  Eco_NPP_CumYr_col(:,:)                            !cumulative ecosystem NPP, [gC d-2]
  real(r8),target,allocatable ::  Eco_HR_CumYr_col(:,:)                             !cumulative ecosystem heterotrophic respiration, [gC d-2]
  real(r8),target,allocatable ::  ECO_HR_CO2_col(:,:)                               !heterotrophic respiration as CO2, [gC d-2 hr-1]
  real(r8),target,allocatable ::  ECO_HR_CO2_vr(:,:,:)                              !Vertically resolved heterotrophic respiraiton as CO2 [gC/d2/hr]
  real(r8),target,allocatable ::  ECO_HR_CH4_vr(:,:,:)                              !heterotrophic respiration as CH4, [gC d-2]  
  real(r8),target,allocatable ::  ECO_HR_CH4_col(:,:)                               !heterotrophic respiration as CH4, [gC d-2]
  real(r8),target,allocatable ::  EcoHavstElmnt_CumYr_col(:,:,:)                    !ecosystem harvest , [gC d-2]
  real(r8),target,allocatable ::  NetNH4Mineralize_CumYr_col(:,:)                   !total NH4 net mineraln (-ve) or immobiln (+ve),[gN d-2]
  real(r8),target,allocatable ::  NetPO4Mineralize_CumYr_col(:,:)                   !total H2PO4 net mineraln (-ve) or immobiln (+ve), [gP d-2]
  real(r8),target,allocatable ::  Canopy_NEE_col(:,:)                               !total net CO2 fixation, [gC d-2 h-1]
  real(r8),target,allocatable ::  LitrFallStrutElms_col(:,:,:)                      !total LitrFall chemical elements, [g d-2 h-1]
  real(r8),target,allocatable ::  ECO_ER_col(:,:)                                   !ecosystem respiration, [gC d-2 h-1]
  real(r8),target,allocatable ::  Eco_NBP_CumYr_col(:,:)                            !cumulative NBP, [gC d-2]
  real(r8),target,allocatable ::  REcoH1PO4DmndSoil_vr(:,:,:)                       !Current HPO4 demand in non-band by all microbial,root,myco populations, [gN d-2 h-1]
  real(r8),target,allocatable ::  RH1PO4EcoDmndSoilPrev_vr(:,:,:)                   !previous time step HPO4 demand in non-band by all microbial,root,myco populations, [gN d-2 h-1]
  real(r8),target,allocatable ::  REcoH1PO4DmndBand_vr(:,:,:)                       !Current HPO4 demand in band by all microbial,root,myco populations, [gN d-2 h-1]
  real(r8),target,allocatable ::  RH1PO4EcoDmndBandPrev_vr(:,:,:)                   !Previous time step HPO4 demand in band by all microbial,root,myco populations, [gN d-2 h-1]
  real(r8),target,allocatable ::  CumDryDepoOM_col(:,:)                              !cumulative C from dry decomposition, [gC d-2]
!----------------------------------------------------------------------

contains
  subroutine InitEcosimBGCFluxData

  implicit none
  allocate(Eco_NetRad_col(JY,JX));         Eco_NetRad_col=0._r8
  allocate(Eco_Heat_Latent_col(JY,JX));         Eco_Heat_Latent_col=0._r8
  allocate(Eco_Heat_Sens_col(JY,JX));         Eco_Heat_Sens_col=0._r8
  allocate(Eco_Heat_GrndSurf_col(JY,JX));         Eco_Heat_GrndSurf_col=0._r8
  allocate(Eco_GPP_CumYr_col(JY,JX));        Eco_GPP_CumYr_col=0._r8
  allocate(Eco_AutoR_CumYr_col(JY,JX));        Eco_AutoR_CumYr_col=0._r8
  allocate(Eco_NPP_CumYr_col(JY,JX));        Eco_NPP_CumYr_col=0._r8
  allocate(Eco_HR_CumYr_col(JY,JX));        Eco_HR_CumYr_col=0._r8
  allocate(ECO_HR_CO2_col(JY,JX));          ECO_HR_CO2_col=0._r8
  allocate(ECO_HR_CO2_vr(0:JZ,JY,JX));     ECO_HR_CO2_vr=0._r8
  allocate(ECO_HR_CH4_vr(0:JZ,JY,JX));     ECO_HR_CH4_vr=0._r8  
  allocate(ECO_HR_CH4_col(JY,JX));          ECO_HR_CH4_col=0._r8
  allocate(EcoHavstElmnt_CumYr_col(NumPlantChemElms,JY,JX));      EcoHavstElmnt_CumYr_col=0._r8
  allocate(NetNH4Mineralize_CumYr_col(JY,JX));      NetNH4Mineralize_CumYr_col=0._r8
  allocate(NetPO4Mineralize_CumYr_col(JY,JX));      NetPO4Mineralize_CumYr_col=0._r8
  allocate(Canopy_NEE_col(JY,JX));       Canopy_NEE_col=0._r8
  allocate(LitrFallStrutElms_col(NumPlantChemElms,JY,JX));       LitrFallStrutElms_col=0._r8
  allocate(ECO_ER_col(JY,JX));        ECO_ER_col=0._r8
  allocate(Eco_NBP_CumYr_col(JY,JX));        Eco_NBP_CumYr_col=0._r8
  allocate(REcoH1PO4DmndSoil_vr(0:JZ,JY,JX));  REcoH1PO4DmndSoil_vr=0._r8
  allocate(RH1PO4EcoDmndSoilPrev_vr(0:JZ,JY,JX));  RH1PO4EcoDmndSoilPrev_vr=0._r8
  allocate(REcoH1PO4DmndBand_vr(0:JZ,JY,JX));  REcoH1PO4DmndBand_vr=0._r8
  allocate(RH1PO4EcoDmndBandPrev_vr(0:JZ,JY,JX));  RH1PO4EcoDmndBandPrev_vr=0._r8
  allocate(CumDryDepoOM_col(JY,JX)); CumDryDepoOM_col=0._r8
  end subroutine InitEcosimBGCFluxData

!----------------------------------------------------------------------
  subroutine DestructEcosimBGCFluxData
  use abortutils, only : destroy
  implicit none
  call destroy(Eco_NetRad_col)
  call destroy(Eco_Heat_Latent_col)
  call destroy(Eco_Heat_Sens_col)
  call destroy(Eco_Heat_GrndSurf_col)
  call destroy(Eco_GPP_CumYr_col)
  call destroy(Eco_AutoR_CumYr_col)
  call destroy(Eco_NPP_CumYr_col)
  call destroy(Eco_HR_CumYr_col)
  call destroy(ECO_HR_CO2_col)
  call destroy(ECO_HR_CO2_vr)
  call destroy(ECO_HR_CH4_vr)  
  call destroy(ECO_HR_CH4_col)
  call destroy(EcoHavstElmnt_CumYr_col)
  call destroy(NetNH4Mineralize_CumYr_col)
  call destroy(NetPO4Mineralize_CumYr_col)
  call destroy(Canopy_NEE_col)
  call destroy(LitrFallStrutElms_col)
  call destroy(ECO_ER_col)
  call destroy(Eco_NBP_CumYr_col)
  call destroy(REcoH1PO4DmndSoil_vr)
  call destroy(RH1PO4EcoDmndSoilPrev_vr)
  call destroy(REcoH1PO4DmndBand_vr)
  call destroy(RH1PO4EcoDmndBandPrev_vr)
  call destroy(CumDryDepoOM_col)
  end subroutine DestructEcosimBGCFluxData

end module EcosimBGCFluxType
