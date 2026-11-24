module PlantSoilChemistryAPIData
  ! Owns the soilchemistry API type, instance and lifecycle procedures.
  use PlantAPICommonData
  implicit none
  private
  public :: plt_soilchem_init, plt_soilchem_destroy

  type, public :: plant_soilchem_type
  real(r8), pointer :: FracBulkSOMC_vr(:,:)           => null()  !fraction of total organic C in complex,       [-]
  real(r8), pointer :: PlantElmAllocMat4Litr(:,:,:,:)      => null() !litter kinetic fraction,                       [-]
  real(r8), pointer :: TScal4Difsvity_vr(:)           => null()  !temperature effect on diffusivity,[-]
  real(r8), pointer :: FracAirFilledSoilPoreM_vr(:,:) => null()  !soil air-filled porosity,                     [m3 m-3]
  real(r8), pointer :: DiffusivitySolutEffM_vr(:,:)   => null()  !coefficient for dissolution - volatilization, [-]
  real(r8), pointer :: SoilBulkModulus4RootPent_vr(:)   => null()  !elastic modulus of the undisturbed soil, [MPa]
  real(r8), pointer :: SoilModulus4RootRadialexp_vr(:) => null() ! soil modulus for root radial expansion, [MPa]
  real(r8), pointer :: SoilBulkDensity_vr(:)          => null()  !soil bulk density,                            [Mg m-3]
  real(r8), pointer :: trc_solcl_vr(:,:)              => null()  !aqueous tracer concentration, [g m-3]
  real(r8), pointer :: trcg_gascl_vr(:,:)             => null()  !gaseous tracer concentration, [g m-3]
  real(r8), pointer :: CSoilOrgM_vr(:,:)              => null()  !soil organic C content, [gC kg soil-1]
  real(r8), pointer :: HYCDMicP4RootUptake_vr(:) => null()  !soil micropore hydraulic conductivity for root water uptake, [m MPa-1 h-1]
  real(r8), pointer :: GasDifcT_vr(:,:)                => null()  !gaseous diffusivity, [m2 h-1]
  real(r8), pointer :: SoluteDifusvtyT_vr(:,:)         => null()  !aqueous diffusivity, [m2 h-1]
  real(r8), pointer :: trcg_gasml_vr(:,:)             => null()  !gas layer mass, [g d-2]
  real(r8), pointer :: GasSolbility_vr(:,:)           => null()  !gas solubility,                                [m3 m-3]
  real(r8), pointer :: THETW_vr(:)                    => null()  !volumetric water content, [m3 m-3]
  real(r8), pointer :: SoilWatAirDry_vr(:)            => null()  !air-dry water content,                        [m3 m-3]
  real(r8), pointer :: VLSoilPoreMicP_vr(:)           => null()  !volume of soil layer,	[m3 d-2]
  real(r8), pointer :: trcs_VLN_vr(:,:)               => null()  !effective relative tracer volume, [-]
  real(r8), pointer :: VLSoilMicP_vr(:)               => null()  !total micropore volume in layer, [m3 d-2]
  real(r8), pointer :: VLiceMicP_vr(:)                => null()  !soil micropore ice content,   [m3 d-2]
  real(r8), pointer :: VLWatMicP_vr(:)                => null()  !soil micropore water content, [m3 d-2]
  real(r8), pointer :: VLMicP_vr(:)                   => null()  !total volume in micropores, [m3 d-2]
  real(r8), pointer :: DOM_MicP_vr(:,:,:)             => null()  !dissolved organic matter in micropore,	[g d-2]
  real(r8), pointer :: DOM_MicP_drib_vr(:,:,:)        => null()  !dribbling flux for micropore dom,[g d-2]
  real(r8), pointer :: trcs_solml_vr(:,:)             => null() !aqueous tracer, [g d-2]
  contains
    procedure, public :: Init => plt_soilchem_init
    procedure, public :: Destroy => plt_soilchem_destroy
  end type plant_soilchem_type


  type(plant_soilchem_type) , public, target :: plt_soilchem  !soil bgc interface with plant root

contains

  subroutine  plt_soilchem_init(this)

  implicit none

  class(plant_soilchem_type) :: this

  allocate(this%FracBulkSOMC_vr(1:jcplx,0:JZ1));this%FracBulkSOMC_vr=spval
  allocate(this%PlantElmAllocMat4Litr(NumPlantChemElms,0:NumLitterGroups,jsken,JP1));this%PlantElmAllocMat4Litr=spval
  allocate(this%TScal4Difsvity_vr(0:JZ1));this%TScal4Difsvity_vr=spval
  allocate(this%FracAirFilledSoilPoreM_vr(60,0:JZ1));this%FracAirFilledSoilPoreM_vr=spval
  allocate(this%DiffusivitySolutEffM_vr(60,0:JZ1));this%DiffusivitySolutEffM_vr=spval
  allocate(this%VLSoilMicP_vr(0:JZ1));this%VLSoilMicP_vr=spval
  allocate(this%VLiceMicP_vr(0:JZ1));this%VLiceMicP_vr=spval
  allocate(this%VLWatMicP_vr(0:JZ1));this%VLWatMicP_vr=spval
  allocate(this%VLMicP_vr(0:JZ1));this%VLMicP_vr=spval
  allocate(this%trcs_VLN_vr(ids_nuts_beg:ids_nuts_end,0:JZ1));this%trcs_VLN_vr=spval
  allocate(this%DOM_MicP_vr(idom_beg:idom_end,1:jcplx,0:JZ1));this%DOM_MicP_vr=spval
  allocate(this%DOM_MicP_drib_vr(idom_beg:idom_end,1:jcplx,0:JZ1));this%DOM_MicP_drib_vr=spval
  allocate(this%trcs_solml_vr(ids_beg:ids_end,0:JZ1));this%trcs_solml_vr=spval
  allocate(this%trcg_gasml_vr(idg_beg:idg_NH3,0:JZ1));this%trcg_gasml_vr=spval
  allocate(this%trcg_gascl_vr(idg_beg:idg_NH3,0:JZ1));this%trcg_gascl_vr=spval
  allocate(this%CSoilOrgM_vr(1:NumPlantChemElms,0:JZ1));this%CSoilOrgM_vr=spval

  allocate(this%trc_solcl_vr(ids_beg:ids_end,0:jZ1));this%trc_solcl_vr=spval

  allocate(this%VLSoilPoreMicP_vr(0:JZ1));this%VLSoilPoreMicP_vr=spval
  allocate(this%THETW_vr(0:JZ1));this%THETW_vr=spval
  allocate(this%SoilWatAirDry_vr(0:JZ1));this%SoilWatAirDry_vr=spval

  allocate(this%GasSolbility_vr(idg_beg:idg_end,0:JZ1));this%GasSolbility_vr=spval
  allocate(this%GasDifcT_vr(idg_beg:idg_end,0:JZ1));this%GasDifcT_vr=spval
  allocate(this%SoluteDifusvtyT_vr(ids_beg:ids_end,0:JZ1));this%SoluteDifusvtyT_vr=spval
  allocate(this%SoilBulkModulus4RootPent_vr(JZ1));this%SoilBulkModulus4RootPent_vr=spval
  allocate(this%SoilModulus4RootRadialexp_vr(JZ1)); this%SoilModulus4RootRadialexp_vr=spval
  allocate(this%SoilBulkDensity_vr(0:JZ1));this%SoilBulkDensity_vr=spval
  allocate(this%HYCDMicP4RootUptake_vr(JZ1));this%HYCDMicP4RootUptake_vr=spval

  end subroutine plt_soilchem_init

  subroutine plt_soilchem_destroy(this)
  implicit none
  class(plant_soilchem_type) :: this

!  if(allocated(FracBulkSOMC_vr))deallocate(FracBulkSOMC_vr)

!  if(allocated(PlantElmAllocMat4Litr))deallocate(PlantElmAllocMat4Litr)
!  if(allocated(TScal4Difsvity_vr))deallocate(TScal4Difsvity_vr)
!  if(allocated(FracAirFilledSoilPoreM_vr))deallocate(FracAirFilledSoilPoreM_vr)
!  if(allocated(DiffusivitySolutEff))deallocate(DiffusivitySolutEff)
!  if(allocated(ZVSGL))deallocate(ZVSGL)
!  if(allocated(O2GSolubility))deallocate(O2GSolubility)

!  if(allocated(HYCDMicP4RootUptake_vr))deallocate(HYCDMicP4RootUptake_vr)
!  if(allocated(CGSGL))deallocate(CGSGL)
!  if(allocated(CHSGL))deallocate(CHSGL)
!  if(allocated(HGSGL))deallocate(HGSGL)
!  if(allocated(OGSGL))deallocate(OGSGL)
!  if(allocated(SoilBulkModulus4RootPent_vr))deallocate(SoilBulkModulus4RootPent_vr)

!   call destroy(this%GasSolbility_vr)
!  if(allocated(THETW))deallocate(THETW)
!  if(allocated(THETY))deallocate(THETY)
!  if(allocated(VLSoilPoreMicP_vr))deallocate(VLSoilPoreMicP_vr)
!  if(allocated(CCH4G))deallocate(CCH4G)
!  if(allocated(CZ2OG))deallocate(CZ2OG)
!  if(allocated(CNH3G))deallocate(CNH3G)
!  if(allocated(CH2GG))deallocate(CH2GG)
!  if(allocated(CORGC))deallocate(CORGC)
!  if(allocated(H2PO4))deallocate(H2PO4)
!  if(allocated(HLSGL))deallocate(HLSGL)

!   call destroy(this%trcs_VLN_vr)
!  if(allocated(O2AquaDiffusvity))deallocate(O2AquaDiffusvity)
!  if(allocated(POSGL))deallocate(POSGL)
!  if(allocated(VLSoilMicP))deallocate(VLSoilMicP)
!  if(allocated(VOLI))deallocate(VOLI)
!  if(allocated(VOLW))deallocate(VOLW)
!  if(allocated(VLMicP))deallocate(VLMicP)
!  if(allocated(ZOSGL))deallocate(ZOSGL)

!  if(allocated(OQC))deallocate(OQC)
!  if(allocated(OQN))deallocate(OQN)
!  if(allocated(OQP))deallocate(OQP)

!  if(allocated(ZNSGL))deallocate(ZNSGL)

!  if(allocated(Z2SGL))deallocate(Z2SGL)
!  if(allocated(ZHSGL))deallocate(ZHSGL)
!  if(allocated(SoilBulkDensity_vr))deallocate(SoilBulkDensity_vr)
!  if(allocated(CO2G))deallocate(CO2G)
!  if(allocated(CLSGL))deallocate(CLSGL)
!  if(allocated(CQSGL))deallocate(CQSGL)

  end subroutine plt_soilchem_destroy
end module PlantSoilChemistryAPIData
