module MicStateTraitTypeMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL     
  use EcoSiMParDataMod, only : micpar
  use TracerIDMod
  use abortutils, only : destroy
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  type, public :: micsttype
  real(r8) :: EPOC
  real(r8) :: EHUM
  real(r8) :: TMicHeterActivity    !heterotrophic micriobial activity
  real(r8) :: ZNH4B
  real(r8) :: ZNH4S
  real(r8) :: ZNO3B
  real(r8) :: ZNO3S
  real(r8) :: H1POB
  real(r8) :: H1PO4
  real(r8) :: ZNO2B
  real(r8) :: ZNO2S
  real(r8) :: H2POB
  real(r8) :: H2PO4
  real(r8) :: CCO2S
  real(r8) :: CNO2S
  real(r8) :: CNO2B
  real(r8) :: CNO3S
  real(r8) :: CNH3B       !NH3 concentation in band soil
  real(r8) :: CNH3S       !NH3 concentation in non-band soil
  real(r8) :: CNO3B
  real(r8) :: CNH4S
  real(r8) :: CNH4B
  real(r8) :: CZ2OS
  real(r8) :: Z2OS
  real(r8) :: OXYG
  real(r8) :: OXYS
  real(r8) :: COXYS
  real(r8) :: O2GSolubility
  real(r8) :: COXYG
  real(r8) :: CH2P4
  real(r8) :: CH2P4B
  real(r8) :: CH1P4
  real(r8) :: CH1P4B
  real(r8) :: CZ2GS
  real(r8) :: CH2GS
  real(r8) :: H2GS
  real(r8) :: CCH4G
  real(r8) :: CCH4S
  real(r8) :: CH4S
  real(r8) :: CH4AquaSolubility
  real(r8) :: ZNFN0
  real(r8) :: ZNFNI
  real(r8) ::ZNH4TU
  real(r8) ::ZNO3TU
  real(r8) ::H1P4TU
  real(r8) ::H2P4TU
  real(r8) ::CNH4BU
  real(r8) ::CNH4SU
  real(r8) ::CH2P4U
  real(r8) ::CH2P4BU
  real(r8) ::CH1P4U
  real(r8) ::CH1P4BU
  real(r8) ::CNO3SU
  real(r8) ::CNO3BU
  real(r8) :: TotActMicrobiom
  real(r8) :: TSens4MicbGrwoth
  real(r8) :: VWatMicrobAct
  real(r8) :: TSolidOMActC
  real(r8) :: TSolidOMC
  real(r8) :: tOMActC
  integer  :: Lay
  real(r8),allocatable :: FracBulkSOMC(:)
  real(r8),allocatable :: DOM(:,:)             !DOM mass, [g d-2]
  real(r8),allocatable :: DOM_MicP_drib(:,:)   !dribbling flux for DOM, [g d-2]
  real(r8),allocatable :: SorbedOM(:,:)        !adsorbed OM, [g d-2]
  real(r8),allocatable :: SolidOMAct(:,:)      !solid OM activity, [g d-2]
  real(r8),allocatable :: SolidOM(:,:,:)       !solid OM mass, [g d-2]
  real(r8),allocatable :: OMBioResdu(:,:,:)    !microbial residule, [g d-2]
  real(r8),allocatable :: CNOSC(:,:)           !CN ratio for solid OM, [gN gC-1]
  real(r8),allocatable :: CPOSC(:,:)           !CP ratio for solid OM, [gP gC-1]
  real(r8),allocatable :: mBiomeHeter(:,:,:)   !heterotrophic microbial biomass, [g d-2]
  real(r8),allocatable :: mBiomeAutor(:,:)     !autotrophic microbial biomass, [g d-2]
  real(r8),allocatable :: SOMPomProtein(:)     !Protein component in POM, [g d-2]
  real(r8),allocatable :: SOMHumProtein(:)     !Protein component in Humus, [g d-2]
  real(r8),allocatable :: SOMHumCarbohyd(:)    !carbonhydro component in humus, [g d-2]

  contains
  procedure, public :: Init
  procedure, public :: Destroy=> Destruct
  end type micsttype
  contains

!------------------------------------------------------------------------------------------
  subroutine Init(this)
  implicit none
  class(micsttype) :: this
  integer :: jcplx,jsken
  integer, pointer :: ndbiomcp, nlbiomcp
  integer, pointer :: NumMicrobAutoTrophCmplx, NumHetetr1MicCmplx
  integer, pointer :: NumLiveHeterBioms
  integer, pointer :: NumLiveAutoBioms


  jcplx=micpar%jcplx
  NumLiveAutoBioms => micpar%NumLiveAutoBioms
  jsken=micpar%jsken
  ndbiomcp =>micpar%ndbiomcp
  nlbiomcp =>micpar%nlbiomcp
  NumMicrobAutoTrophCmplx=>micpar%NumMicrobAutoTrophCmplx
  NumHetetr1MicCmplx=>micpar%NumHetetr1MicCmplx
  NumLiveHeterBioms => micpar%NumLiveHeterBioms

  allocate(this%FracBulkSOMC(1:jcplx));this%FracBulkSOMC=spval
  allocate(this%DOM(idom_beg:idom_end,1:jcplx));this%DOM=spval
  allocate(this%dom_micp_drib(idom_beg:idom_end,1:jcplx));this%dom_micp_drib=spval
  allocate(this%SorbedOM(idom_beg:idom_end,1:jcplx));this%SorbedOM=spval
  allocate(this%SolidOMAct(jsken,1:jcplx));this%SolidOMAct=spval
  allocate(this%SolidOM(NumPlantChemElms,jsken,1:jcplx));this%SolidOM=spval
  allocate(this%OMBioResdu(1:NumPlantChemElms,ndbiomcp,1:jcplx));this%OMBioResdu=spval
  allocate(this%CNOSC(jsken,1:jcplx));this%CNOSC=spval
  allocate(this%CPOSC(jsken,1:jcplx));this%CPOSC=spval
  allocate(this%mBiomeHeter(NumPlantChemElms,NumLiveHeterBioms,1:jcplx));this%mBiomeHeter=spval
  allocate(this%mBiomeAutor(NumPlantChemElms,NumLiveAutoBioms));this%mBiomeAutor=spval

  allocate(this%SOMPomProtein(1:NumPlantChemElms));
  allocate(this%SOMHumProtein(1:NumPlantChemElms));
  allocate(this%SOMHumCarbohyd(1:NumPlantChemElms));

  end subroutine Init

!------------------------------------------------------------------------------------------


  subroutine Destruct(this)

  implicit none
  class(micsttype) :: this

  call destroy(this%FracBulkSOMC)
  call destroy(this%DOM)
  call destroy(this%SorbedOM)
  call destroy(this%SolidOMAct)
  call destroy(this%SolidOM)
  call destroy(this%OMBioResdu)
  call destroy(this%CNOSC)
  call destroy(this%CPOSC)
  call destroy(this%mBiomeHeter)
  call destroy(this%mBiomeAutor)
  call destroy(this%SOMHumCarbohyd)
  call destroy(this%SOMPomProtein)
  call destroy(this%SOMHumProtein)
  call destroy(this%dom_micp_drib)
  end subroutine Destruct
end module MicStateTraitTypeMod
