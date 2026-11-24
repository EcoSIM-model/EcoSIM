module ChemTranspDataType

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use EcoSIMConfig, only : jcplx => jcplxc
  use TracerIDMod
  implicit none

  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),target,allocatable ::  TScal4Difsvity_vr(:,:,:)                   !temperature effect on diffusivity [-]

  real(r8),target,allocatable ::  GasDifcT_vr(:,:,:,:)                        !temperature scaled gaseous diffusivity, [m2 h-1]
  real(r8),target,allocatable ::  SoluteDifusvtyT_vr(:,:,:,:)                 !temperature scaled solute diffusivity, [m2 h-1]

  real(r8),target,allocatable ::  O2AquaDiffusvity(:,:,:)                    !aqueous CO2 diffusivity,	[m2 h-1]

  real(r8),target,allocatable ::  DOMdiffusivity_vr(:,:,:,:)                 !aqueous DOC diffusivity in soil, [m2 h-1]
  real(r8),target,allocatable ::  WVapDifusvitySoil_vr(:,:,:)                !water vapor diffusivity in soil, [m2 h-1]
  real(r8),target,allocatable ::  H2OVapDifsc_snvr(:,:,:)                    !water vapor diffusivity in snow, [m2 h-1]
  real(r8),target,allocatable ::  VaporDiffusivityLitR_col(:,:)              !water vapor diffusivity in litter, [m2 h-1]
  real(r8),target,allocatable ::  WVapDifusvityAir_col(:,:)                  !water vapor diffusivity in air, [m2 h-1]

  real(r8),target,allocatable ::  GasSolbility_vr(:,:,:,:)                   !solubility of gases, [-]

  real(r8),target,allocatable ::  RGasTranspFlxPrev_vr(:,:,:,:)              !net gaseous flux from previous time step, [g d-2 h-1]
  real(r8),target,allocatable ::  RCH4PhysexchPrev_vr(:,:,:)                 !net aqueous CH4 flux from previous time step, [g d-2 h-1]
  real(r8),target,allocatable ::  RO2AquaSourcePrev_vr(:,:,:)                !net aqueous O2 flux from previous time step, [g d-2 h-1]
  real(r8),target,allocatable ::  AquaIonDifusivty_vr(:,:,:,:)
  real(r8),target,allocatable ::  trcSalt_FloXSurRunoff_2D(:,:,:,:,:)        !total Al in runoff, [mol d-2 h-1]
  real(r8),target,allocatable ::  trcg_FloXSurRunoff_2D(:,:,:,:,:)           !surface runoff gas flux, [g d-2 h-1]
  real(r8),target,allocatable ::  trcn_FloXSurRunoff_2D(:,:,:,:,:)           !surface runoff nutrient flux, [g d-2 h-1]
  real(r8),target,allocatable ::  DOM_FloXSurRunoff_2DH(:,:,:,:,:,:)         !surface runoff DOC flux, [g d-2 h-1]
  private :: InitAllocate

  contains


  subroutine InitChemTranspData(lsalt_model)

  implicit none
  logical, intent(in) :: lsalt_model

  call InitAllocate(lsalt_model)

  end subroutine InitChemTranspData
!------------------------------------------------------------------------------------------

  subroutine InitAllocate(lsalt_model)
  implicit none
  logical, intent(in) :: lsalt_model
  
  allocate(TScal4Difsvity_vr(0:JZ,JY,JX));   TScal4Difsvity_vr=0._r8
  allocate(GasDifcT_vr(idg_beg:idg_end,JZ,JY,JX));GasDifcT_vr=0._r8
  allocate(SoluteDifusvtyT_vr(ids_beg:ids_end,0:JZ,JY,JX));SoluteDifusvtyT_vr=0._r8

  allocate(O2AquaDiffusvity(0:JZ,JY,JX));  O2AquaDiffusvity=0._r8

  allocate(DOMdiffusivity_vr(idom_beg:idom_end,0:JZ,JY,JX));  DOMdiffusivity_vr=0._r8
  allocate(WVapDifusvitySoil_vr(JZ,JY,JX));    WVapDifusvitySoil_vr=0._r8
  allocate(H2OVapDifsc_snvr(JS,JY,JX));    H2OVapDifsc_snvr=0._r8
  allocate(VaporDiffusivityLitR_col(JY,JX));       VaporDiffusivityLitR_col=0._r8
  allocate(WVapDifusvityAir_col(JY,JX));       WVapDifusvityAir_col=0._r8
  allocate(GasSolbility_vr(idg_beg:idg_end,0:JZ,JY,JX)); GasSolbility_vr=0._r8
  
  allocate(RCH4PhysexchPrev_vr(0:JZ,JY,JX));  RCH4PhysexchPrev_vr=0._r8
  allocate(RGasTranspFlxPrev_vr(idg_beg:idg_end,0:JZ,JY,JX));  RGasTranspFlxPrev_vr=0._r8
  allocate(RO2AquaSourcePrev_vr(0:JZ,JY,JX));  RO2AquaSourcePrev_vr=0._r8

  if(lsalt_model)then
    allocate(AquaIonDifusivty_vr(idsalt_beg:idsalt_mend,JZ,JY,JX));AquaIonDifusivty_vr=0._r8
  endif


  allocate(trcSalt_FloXSurRunoff_2D(idsalt_beg:idsalt_end,2,2,JV,JH));   trcSalt_FloXSurRunoff_2D=0._r8
  allocate(trcg_FloXSurRunoff_2D(idg_beg:idg_NH3,2,2,JV,JH));  trcg_FloXSurRunoff_2D=0._r8
  allocate(trcn_FloXSurRunoff_2D(ids_nut_beg:ids_nuts_end,2,2,JV,JH));  trcn_FloXSurRunoff_2D=0._r8
  allocate(DOM_FloXSurRunoff_2DH(idom_beg:idom_end,1:jcplx,2,2,JV,JH));DOM_FloXSurRunoff_2DH=0._r8

  end subroutine InitAllocate
!------------------------------------------------------------------------------------------


  subroutine DestructChemTranspData
  use abortutils, only : destroy

  implicit none

  call destroy(trcSalt_FloXSurRunoff_2D)
  call destroy(GasDifcT_vr)
  call destroy(SoluteDifusvtyT_vr)
  call destroy(TScal4Difsvity_vr)

  call destroy(O2AquaDiffusvity)
  call destroy(DOMdiffusivity_vr)
  call destroy(WVapDifusvitySoil_vr)
  call destroy(H2OVapDifsc_snvr)
  call destroy(VaporDiffusivityLitR_col)
  call destroy(WVapDifusvityAir_col)

  call destroy(trcn_FloXSurRunoff_2D)
  call destroy(trcg_FloXSurRunoff_2D)
  call destroy(GasSolbility_vr)
  call destroy(AquaIonDifusivty_vr)
  call destroy(RGasTranspFlxPrev_vr)
  call destroy(RCH4PhysexchPrev_vr)
  call destroy(RO2AquaSourcePrev_vr)
  call destroy(DOM_FloXSurRunoff_2DH)
  end subroutine DestructChemTranspData
end module ChemTranspDataType
