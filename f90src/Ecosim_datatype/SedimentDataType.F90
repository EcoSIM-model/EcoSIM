module SedimentDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,ndbiomcp=>NumDeadMicrbCompts,jsken=>jskenc
implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  TSED_col(:,:)                     !erosion rate [Mg d-2 t-1]
  real(r8),target,allocatable ::  SoilDetachability4Erosion1(:,:)   !soil detachment, [-]
  real(r8),target,allocatable ::  SoilDetachability4Erosion2(:,:)   !soil detachability
  real(r8),target,allocatable ::  CER_col(:,:)                      !soil detachment/deposition, [h-1]
  real(r8),target,allocatable ::  XER_col(:,:)                      !soil detachment/deposition shape parameter, [-]
  real(r8),target,allocatable ::  PtclDensitySurfLay_col(:,:)     !particle density of surface layer,[-]
  real(r8),target,allocatable ::  VLS_col(:,:)                      !hourly sinking rate,[h-1]
  real(r8),target,allocatable ::  SED_col(:,:)                      !sediment transport, [Mg d-2 h-1]
  real(r8),target,allocatable ::  XSand_Eros_2D(:,:,:,:)            !total sand erosion , [Mg d-2 h-1]
  real(r8),target,allocatable ::  XSilt_Eros_2D(:,:,:,:)            !total silt erosion , [Mg d-2 h-1]
  real(r8),target,allocatable ::  XClay_Eros_2D(:,:,:,:)            !total clay erosion , [Mg d-2 h-1]
  real(r8),target,allocatable ::  XNH4Soil_Eros_molN_2D(:,:,:,:)    !total NH4 fertilizer erosion non-band , [molN d-2 h-1]
  real(r8),target,allocatable ::  XNH3Soil_Eros_molN_2D(:,:,:,:)    !total NH3 fertilizer erosion non-band, [molN d-2 h-1]
  real(r8),target,allocatable ::  XUreaSoil_Eros_molN_2D(:,:,:,:)   !total urea fertilizer erosion non-band , [molN d-2 h-1]
  real(r8),target,allocatable ::  XNO3Soil_Eros_molN_2D(:,:,:,:)    !total NO3 fertilizer erosion non-band , [molN d-2 h-1]
  real(r8),target,allocatable ::  XPO4Soil_Eros_molP_2D(:,:,:,:)    !total PO4 fertilizer erosion non-band, [molP d-2 h-1]
  real(r8),target,allocatable ::  XNH4Band_Eros_molN_2D(:,:,:,:)    !total NH4 fertilizer erosion band , [molN d-2 h-1]
  real(r8),target,allocatable ::  XNH3Band_Eros_molN_2D(:,:,:,:)    !total NH3 fertilizer erosion band, [molN d-2 h-1]
  real(r8),target,allocatable ::  XUreaBand_Eros_molN_2D(:,:,:,:)   !total urea fertilizer erosion band , [molN d-2 h-1]
  real(r8),target,allocatable ::  XNO3Band_Eros_molN_2D(:,:,:,:)    !total NO3 fertilizer erosion band , [molN d-2 h-1]
  real(r8),target,allocatable ::  XPO4Band_Eros_molP_2D(:,:,:,:)    !Total band PO4 fertilizer erosion loss, [molP d-2 h-1]
  real(r8),target,allocatable ::  trcx_Eros_mol_2D(:,:,:,:,:)       !total adsorbed sediment erosion non-band , [mol d-2 h-1]
  real(r8),target,allocatable ::  trcp_Eros_mol_2D(:,:,:,:,:)       !total adsorbed ALOH3  erosion , [mol d-2 h-1]
  real(r8),target,allocatable ::  cumSed_Eros_2D(:,:,:,:)           !sediment erosion, [Mg d-2 h-1]
  real(r8),target,allocatable ::  OMBioResdu_Eros_2D(:,:,:,:,:,:,:) !microbial residue C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  SorbedOM_Eros_2D(:,:,:,:,:,:)     !adsorbed C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  SolidOM_Eros_2D(:,:,:,:,:,:,:)    !humus C  erosion , [g d-2 h-1]
  real(r8),target,allocatable ::  SolidOMAct_Eros_2D(:,:,:,:,:,:)   !colonized humus C  erosion , [g d-2 h-1]
  private :: InitAllocate

  contains
!----------------------------------------------------------------------

  subroutine InitSedimentData

  implicit none

  call InitAllocate
  end subroutine InitSedimentData

!----------------------------------------------------------------------

  subroutine InitAllocate

  implicit none
  allocate(TSED_col(JY,JX));        TSED_col=0._r8
  allocate(SoilDetachability4Erosion1(JY,JX));        SoilDetachability4Erosion1=0._r8
  allocate(SoilDetachability4Erosion2(JY,JX));        SoilDetachability4Erosion2=0._r8
  allocate(CER_col(JY,JX));         CER_col=0._r8
  allocate(XER_col(JY,JX));         XER_col=0._r8
  allocate(PtclDensitySurfLay_col(JY,JX));      PtclDensitySurfLay_col=0._r8
  allocate(VLS_col(JY,JX));         VLS_col=0._r8
  allocate(SED_col(JY,JX));         SED_col=0._r8
  allocate(XSand_Eros_2D(2,2,JV,JH));  XSand_Eros_2D=0._r8
  allocate(XSilt_Eros_2D(2,2,JV,JH));  XSilt_Eros_2D=0._r8
  allocate(XClay_Eros_2D(2,2,JV,JH));  XClay_Eros_2D=0._r8
  allocate(XNH4Soil_Eros_molN_2D(2,2,JV,JH));  XNH4Soil_Eros_molN_2D=0._r8
  allocate(XNH3Soil_Eros_molN_2D(2,2,JV,JH));  XNH3Soil_Eros_molN_2D=0._r8
  allocate(XUreaSoil_Eros_molN_2D(2,2,JV,JH));  XUreaSoil_Eros_molN_2D=0._r8
  allocate(XNO3Soil_Eros_molN_2D(2,2,JV,JH));  XNO3Soil_Eros_molN_2D=0._r8
  allocate(XNH4Band_Eros_molN_2D(2,2,JV,JH));  XNH4Band_Eros_molN_2D=0._r8
  allocate(XNH3Band_Eros_molN_2D(2,2,JV,JH));  XNH3Band_Eros_molN_2D=0._r8
  allocate(XUreaBand_Eros_molN_2D(2,2,JV,JH));  XUreaBand_Eros_molN_2D=0._r8
  allocate(XNO3Band_Eros_molN_2D(2,2,JV,JH));  XNO3Band_Eros_molN_2D=0._r8
  allocate(trcx_Eros_mol_2D(idx_beg:idx_end,2,2,JV,JH));   trcx_Eros_mol_2D=0._r8
  allocate(trcp_Eros_mol_2D(idsp_beg:idsp_end,2,2,JV,JH));  trcp_Eros_mol_2D=0._r8
  allocate(XPO4Soil_Eros_molP_2D(2,2,JV,JH)); XPO4Soil_Eros_molP_2D = 0._r8
  allocate(XPO4Band_Eros_molP_2D(2,2,JV,JH)); XPO4Band_Eros_molP_2D=0._r8
  allocate(cumSed_Eros_2D(2,2,JV,JH));  cumSed_Eros_2D=0._r8
  allocate(OMBioResdu_Eros_2D(NumPlantChemElms,ndbiomcp,1:jcplx,2,2,JV,JH));OMBioResdu_Eros_2D=0._r8
  allocate(SorbedOM_Eros_2D(idom_beg:idom_end,1:jcplx,2,2,JV,JH));SorbedOM_Eros_2D=0._r8
  allocate(SolidOM_Eros_2D(NumPlantChemElms,jsken,1:jcplx,2,2,JV,JH));SolidOM_Eros_2D=0._r8
  allocate(SolidOMAct_Eros_2D(jsken,1:jcplx,2,2,JV,JH));SolidOMAct_Eros_2D=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructSedimentData

  use abortutils, only : destroy
  implicit none
  call destroy(XPO4Soil_Eros_molP_2D)
  call destroy(XPO4Band_Eros_molP_2D)
  call destroy(trcp_Eros_mol_2D)
  call destroy(trcx_Eros_mol_2D)
  call destroy(TSED_col)
  call destroy(SoilDetachability4Erosion1)
  call destroy(SoilDetachability4Erosion2)
  call destroy(CER_col)
  call destroy(XER_col)
  call destroy(PtclDensitySurfLay_col)
  call destroy(VLS_col)
  call destroy(SED_col)
  call destroy(XSand_Eros_2D)
  call destroy(XSilt_Eros_2D)
  call destroy(XClay_Eros_2D)
  call destroy(XNH4Soil_Eros_molN_2D)
  call destroy(XNH3Soil_Eros_molN_2D)
  call destroy(XUreaSoil_Eros_molN_2D)
  call destroy(XNO3Soil_Eros_molN_2D)
  call destroy(XNH4Band_Eros_molN_2D)
  call destroy(XNH3Band_Eros_molN_2D)
  call destroy(XUreaBand_Eros_molN_2D)
  call destroy(XNO3Band_Eros_molN_2D)
  call destroy(cumSed_Eros_2D)
  call destroy(OMBioResdu_Eros_2D)
  call destroy(SorbedOM_Eros_2D)
  call destroy(SolidOM_Eros_2D)
  call destroy(SolidOMAct_Eros_2D)
  end subroutine DestructSedimentData

end module SedimentDataType
