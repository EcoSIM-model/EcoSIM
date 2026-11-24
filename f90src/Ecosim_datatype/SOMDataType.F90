module SOMDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc, ndbiomcp=>NumDeadMicrbCompts
  use data_const_mod, only : spval => DAT_CONST_SPVAL  
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable :: RSC_vr(:,:,:,:)                       !initial surface litter C, [g m-2]
  real(r8),target,allocatable :: RSN_vr(:,:,:,:)                       !initial surface litter N, [g m-2]
  real(r8),target,allocatable :: RSP_vr(:,:,:,:)                       !initial surface litter P, [g m-2]
  real(r8),target,allocatable :: CFOSC_vr(:,:,:,:,:)                   !fraction of SOC in kinetic components,[-]
  real(r8),target,allocatable :: CNOSC_vr(:,:,:,:,:)                   !N:C,ratios of SOC kinetic components,[-]
  real(r8),target,allocatable :: CPOSC_vr(:,:,:,:,:)                   !P:C ratios of SOC kinetic components,[-]
  real(r8),target,allocatable :: SolidOM_vr(:,:,:,:,:,:)               !humus soil OM chemical element,	[g d-2]
  real(r8),target,allocatable :: TSolidOMActC_vr(:,:,:)                !total active solid organic C, [gC d-2]
  real(r8),target,allocatable :: TSolidOMC_vr(:,:,:)                  !total solid organic C, [gC d-2]
  real(r8),target,allocatable :: tOMActC_vr(:,:,:)                    !active heterotrophic microbial C in layer, [gC d-2]
  real(r8),target,allocatable :: SorbedOM_vr(:,:,:,:,:)               !adsorbed soil OM chemical element,	[g d-2]
  real(r8),target,allocatable :: OMBioResdu_vr(:,:,:,:,:,:)           !microbial residue chemical element, [C	g d-2]
  real(r8),target,allocatable :: DOM_MicP_vr(:,:,:,:,:)               !dissolved organic matter in micropore,	[g d-2]
  real(r8),target,allocatable :: DOM_MacP_vr(:,:,:,:,:)               !dissolved organic matter in macropore,	[g d-2]
  real(r8),target,allocatable :: SoilOrgM_vr(:,:,:,:)                 !total soil organic matter, [g d-2]
  real(r8),target,allocatable :: ORGCX_vr(:,:,:)                      !SOC concentration,	[gC Mg soil-1]
  real(r8),target,allocatable :: SolidOMAct_vr(:,:,:,:,:)             !colonized humus C in each complex, [g d-2]
  real(r8),target,allocatable :: OMLitrC_vr(:,:,:)                    !total particulate organic C, [g d-2]
  real(r8),target,allocatable :: CSoilOrgM_vr(:,:,:,:)                 !soil organic matter content, input as [kg Mg soil-1], used as [g Mg soil -1]
  real(r8),target,allocatable :: SolidOMPercent_vr(:,:,:)             !percentage of OM as soil mass, [%]
  real(r8),target,allocatable :: COMLitrC_vr(:,:,:)                    !soil litter particulate C content, [g kg-1]
  real(r8),target,allocatable :: ElmAllocmatMicrblitr2POM_vr(:,:,:,:)    !allocation coefficient to humus fractions,[-]
  real(r8),target,allocatable :: tMicBiome_col(:,:,:)                         !total micriobial biomass chemical element, [g d-2]
  real(r8),target,allocatable :: tSoilOrgM_col(:,:,:)                  !total soil organic matter, include everything organic (exclude live roots), [g d-2]
  real(r8),target,allocatable :: tLitrOM_col(:,:,:)                        !total litter chemical element, [g d-2]
  real(r8),target,allocatable :: litrOM_vr(:,:,:,:)                    !vertical layered litter chemical element, [g d-2]
  real(r8),target,allocatable :: tHumOM_col(:,:,:)                        !total humus chemical element, [g d-2]
  real(r8),target,allocatable :: EPOC_vr(:,:,:)                         !partitioning coefficient between POC and litter, [-]
  real(r8),target,allocatable :: EHUM_vr(:,:,:)                         !partitioning coefficient between humus and microbial residue, [-]
  real(r8),target,allocatable :: CDOM_vr(:,:,:,:,:)                     !DOC concentration, [g m-3]
  real(r8),target,allocatable :: FracBulkSOMC_vr(:,:,:,:)                    !fraction of total organic C in complex, [-]
  real(r8),target,allocatable :: DIC_mass_col(:,:)                        !total soil DIC, [g d-2]
  real(r8),target,allocatable :: tNH4_col(:,:)                         !total soil NH4 + NH3 content, [g d-2]
  real(r8),target,allocatable :: tNO3_col(:,:)                         !total soil NO3 + NO2 content, [g d-2]
  real(r8),target,allocatable :: tHxPO4_col(:,:)                         !total soil PO4 content, [g d-2]
  real(r8),target,allocatable :: FracLitrMix_vr(:,:,:)              !fraction of litter to be mixed downward,[-]
  real(r8),target,allocatable :: RHydlySOCK_vr(:,:,:,:)             !hydrolysis of soil organic C in each complex, [gC d-2 h-1]
  private :: InitAllocate
  contains

  subroutine InitSOMData(NumOfLitrCmplxs)

  implicit none
  integer, intent(in) :: NumOfLitrCmplxs

  call InitAllocate(NumOfLitrCmplxs)

  end subroutine InitSOMData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate(NumOfLitrCmplxs)
  implicit none
  integer, intent(in) :: NumOfLitrCmplxs

  allocate(RSC_vr(1:NumOfLitrCmplxs,0:JZ,JY,JX)); RSC_vr=0._r8
  allocate(RSN_vr(1:NumOfLitrCmplxs,0:JZ,JY,JX)); RSN_vr=0._r8
  allocate(RSP_vr(1:NumOfLitrCmplxs,0:JZ,JY,JX)); RSP_vr=0._r8
  allocate(CFOSC_vr(jsken,1:jcplx,0:JZ,JY,JX)); CFOSC_vr=0._r8
  allocate(CNOSC_vr(jsken,1:jcplx,0:JZ,JY,JX)); CNOSC_vr=0._r8
  allocate(CPOSC_vr(jsken,1:jcplx,0:JZ,JY,JX)); CPOSC_vr=0._r8
  allocate(SolidOM_vr(NumPlantChemElms,jsken,1:jcplx,0:JZ,JY,JX));SolidOM_vr=0._r8
  allocate(TSolidOMActC_vr(0:JZ,JY,JX));TSolidOMActC_vr=0._r8
  allocate(TSolidOMC_vr(0:JZ,JY,JX)); TSolidOMC_vr=0._r8
  allocate(tOMActC_vr(0:JZ,JY,JX)); tOMActC_vr=0._r8
  allocate(SorbedOM_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));SorbedOM_vr=0._r8
  allocate(OMBioResdu_vr(1:NumPlantChemElms,ndbiomcp,1:jcplx,0:JZ,JY,JX)); OMBioResdu_vr=0._r8
  allocate(DOM_MicP_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));DOM_MicP_vr=0._r8
  allocate(DOM_MacP_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));DOM_MacP_vr=0._r8
  allocate(SoilOrgM_vr(1:NumPlantChemElms,0:JZ,JY,JX));SoilOrgM_vr=0._r8
  allocate(ORGCX_vr(0:JZ,JY,JX)); ORGCX_vr=0._r8
  allocate(SolidOMAct_vr(jsken,1:jcplx,0:JZ,JY,JX));SolidOMAct_vr=0._r8
  allocate(OMLitrC_vr(0:JZ,JY,JX));OMLitrC_vr=0._r8
  allocate(CSoilOrgM_vr(1:NumPlantChemElms,0:JZ,JY,JX));CSoilOrgM_vr=0._r8
  allocate(SolidOMPercent_vr(JZ,JY,JX)); SolidOMPercent_vr=0._r8
  allocate(COMLitrC_vr(JZ,JY,JX));COMLitrC_vr=0._r8
  allocate(ElmAllocmatMicrblitr2POM_vr(2,JZ,JY,JX));; ElmAllocmatMicrblitr2POM_vr=0._r8
  allocate(tMicBiome_col(NumPlantChemElms,JY,JX));        tMicBiome_col=0._r8
  allocate(tSoilOrgM_col(NumPlantChemElms,JY,JX)); tSoilOrgM_col=0._r8
  allocate(tLitrOM_col(NumPlantChemElms,JY,JX));       tLitrOM_col=0._r8
  allocate(litrOM_vr(NumPlantChemElms,0:JZ,JY,JX)); litrOM_vr=0._r8
  allocate(tHumOM_col(NumPlantChemElms,JY,JX));       tHumOM_col=0._r8
  allocate(EPOC_vr(0:JZ,JY,JX));   EPOC_vr=0._r8
  allocate(EHUM_vr(0:JZ,JY,JX));   EHUM_vr=0._r8
  allocate(CDOM_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));CDOM_vr=0._r8
  allocate(FracBulkSOMC_vr(1:jcplx,0:JZ,JY,JX));FracBulkSOMC_vr=0._r8
  allocate(DIC_mass_col(JY,JX));       DIC_mass_col=0._r8
  allocate(tNH4_col(JY,JX));        tNH4_col=0._r8
  allocate(tNO3_col(JY,JX));        tNO3_col=0._r8
  allocate(tHxPO4_col(JY,JX));        tHxPO4_col=0._r8
  allocate(FracLitrMix_vr(0:JZ,JY,JX)); FracLitrMix_vr=0._r8
  allocate(RHydlySOCK_vr(1:jcplx,0:JZ,JY,JX)); RHydlySOCK_vr=0._r8
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------

  subroutine DestructSOMData

  use abortutils, only : destroy
  implicit none

  call destroy(FracLitrMix_vr)
  call destroy(TSolidOMActC_vr)
  call destroy(TSolidOMC_vr)
  call destroy(tOMActC_vr)
  call destroy(litrOM_vr)
  call destroy(RSC_vr)
  call destroy(RSN_vr)
  call destroy(RSP_vr)
  call destroy(CFOSC_vr)
  call destroy(CNOSC_vr)
  call destroy(CPOSC_vr)
  call destroy(SolidOM_vr)
  call destroy(SorbedOM_vr)
  call destroy(OMBioResdu_vr)
  call destroy(DOM_MicP_vr)
  call destroy(DOM_MacP_vr)
  call destroy(SoilOrgM_vr)
  call destroy(ORGCX_vr)
  call destroy(SolidOMAct_vr)
  call destroy(OMLitrC_vr)
  call destroy(CSoilOrgM_vr)
  call destroy(SolidOMPercent_vr)
  call destroy(COMLitrC_vr)
  call destroy(ElmAllocmatMicrblitr2POM_vr)
  call destroy(tMicBiome_col)
  call destroy(tSoilOrgM_col)
  call destroy(tLitrOM_col)
  call destroy(tHumOM_col)
  call destroy(EPOC_vr)
  call destroy(EHUM_vr)
  call destroy(RHydlySOCK_vr)
  call destroy(CDOM_vr)
  call destroy(FracBulkSOMC_vr)
  call destroy(DIC_mass_col)
  call destroy(tNH4_col)
  call destroy(tNO3_col)
  call destroy(tHxPO4_col)
  end subroutine DestructSOMData
end module SOMDataType
