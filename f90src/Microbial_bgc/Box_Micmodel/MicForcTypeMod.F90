module MicForcTypeMod

  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL   
  use EcoSiMParDataMod, only : micpar
  use abortutils, only : destroy
  implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  type, public :: micforctype
  integer  :: L                 !layer number
  real(r8) :: CCH4E
  real(r8) :: COXYE
  real(r8) :: O2_irrig_conc
  real(r8) :: O2_rain_conc
  real(r8) :: Irrig2LitRSurf_col
  real(r8) :: Rain2LitRSurf
  real(r8) :: TempOffset
  real(r8) :: VLitR
  real(r8) :: VWatLitRHoldCapcity
  real(r8) :: ZEROS2
  real(r8) :: ZEROS
  real(r8) :: VLWatMicP
  real(r8) :: VLsoiAirP
  real(r8) :: VOLW0
  real(r8) :: VLSoilMicP
  real(r8) :: THETY
  real(r8) :: POROS
  real(r8) :: TKS
  real(r8) :: FieldCapacity
  real(r8) :: THETW
  real(r8) :: pH
  real(r8) :: ZERO
  real(r8) :: PAR_rad              !PAR radiation [umol photon m-2 s-1]
  real(r8) :: SoilMicPMassLayer
  real(r8) :: VLSoilPoreMicP
  real(r8) :: TScal4Difsvity
  real(r8) :: VLNOB
  real(r8) :: VLNO3
  real(r8) :: VLNH4
  real(r8) :: VLNHB
  real(r8) :: VLPO4
  real(r8) :: VLPOB
  real(r8) :: AEC
  real(r8) :: PSISoilMatricP
  real(r8) :: O2AquaDiffusvity
  real(r8) :: ORGC
  real(r8) :: RNO2EcoUptkSoilPrev
  real(r8) :: RNO2EcoUptkBandPrev
  real(r8) :: RN2OEcoUptkSoilPrev
  real(r8) :: RO2EcoDmndPrev
  real(r8) :: RCH4EcoDmndPrev
  real(r8) :: RO2GasXchangePrev
  real(r8) :: RCH4PhysexchPrev
  real(r8) :: RCH4GasXchangePrev
  real(r8) :: RNH4EcoDmndSoilPrev
  real(r8) :: RNH4EcoDmndBandPrev
  real(r8) :: RNO3EcoDmndSoilPrev
  real(r8) :: RNO3EcoDmndBandPrev
  real(r8) :: RH2PO4EcoDmndSoilPrev
  real(r8) :: RH2PO4EcoDmndBandPrev
  real(r8) :: RH1PO4EcoDmndSoilPrev
  real(r8) :: RH1PO4EcoDmndBandPrev
  real(r8) :: SoilMicPMassLayer0
  logical  :: LitrM                     !true if it is the litter layer
  logical  :: Lsurf
  real(r8) :: RNH4EcoDmndLitrPrev
  real(r8) :: RNO3EcoDmndLitrPrev
  real(r8) :: RH2PO4EcoDmndLitrPrev
  real(r8) :: RH1PO4EcoDmndLitrPrev
  real(r8) :: RO2AquaXchangePrev
  real(r8) :: tRNH4MicrbImobilSoil
  real(r8) :: tRNH4MicrbImobilBand
  real(r8) :: tRNO3MicrbImobilSoil
  real(r8) :: tRNO3MicrbImobilBand
  real(r8) :: tRH2PO4MicrbImobilSoil
  real(r8) :: tRH2PO4MicrbImobilBand
  real(r8) :: tRH1PO4MicrbImobilSoil
  real(r8) :: tRH1PO4MicrbImobilBand

  !Underlying-soil compartments, independent of the litter layer's fractions.
  !Dimensions: (nonband/band, NH4/NO3/H2PO4/HPO4).
  real(r8) :: LitrSoilNutFrac(2,4)=0._r8
  real(r8) :: LitrSoilNutPool(2,4)=0._r8
  real(r8) :: LitrSoilNutConc(2,4)=0._r8
  real(r8) :: LitrSoilBandDemandPrev(4)=0._r8
  real(r8) :: VOLWU
  real(r8), allocatable :: ElmAllocmatMicrblitr2POM(:)
  real(r8), allocatable :: ElmAllocmatMicrblitr2POMU(:)
  real(r8), allocatable :: RDOMEcoDmndPrev(:)
  real(r8), allocatable :: RAcetateUptkHeterPrev(:,:)    
  real(r8), allocatable :: RAcetateEcoDmndPrev(:)
  real(r8), allocatable :: RDOCUptkHeterPrev(:,:)        
  real(r8), allocatable :: DiffusivitySolutEff(:)  !rate constant for air-water gas exchange
  real(r8), allocatable :: FILM(:)
  real(r8), allocatable :: THETPM(:)
  real(r8), allocatable :: VLWatMicPM(:)
  real(r8), allocatable :: TortMicPM(:)
  real(r8), allocatable :: VLsoiAirPM(:)
  real(r8),allocatable :: AttenfNH4HeterR(:,:)      !litter layer NH4 relative demand
  real(r8),allocatable :: AttenfNO3HeterR(:,:)      !litter layer NO3 relative demand
  real(r8),allocatable :: AttenfH2PO4HeterR(:,:)    !litter layer H2PO4 relative demand
  real(r8),allocatable :: AttenfH1PO4HeterR(:,:)    !litter layer H1PO4 relative demand
  real(r8),allocatable :: AttenfNH4AutorR(:)
  real(r8),allocatable :: AttenfNO3AutorR(:)
  real(r8),allocatable :: AttenfH2PO4AutorR(:)
  real(r8),allocatable :: AttenfH1PO4AutorR(:)
  contains
   procedure, public :: Init
   procedure, public :: destroy=>Destruct
  end type micforctype

  contains
!------------------------------------------------------------------------------------------

  subroutine Init(this)

  implicit none
  class(micforctype) :: this
  integer :: jcplx,NPH,NumHetetr1MicCmplx,NumMicrobAutoTrophCmplx

  jcplx=micpar%jcplx
  NPH=60
  NumMicrobAutoTrophCmplx=micpar%NumMicrobAutoTrophCmplx
  NumHetetr1MicCmplx=micpar%NumHetetr1MicCmplx
  allocate(this%ElmAllocmatMicrblitr2POM(1:micpar%ndbiomcp));this%ElmAllocmatMicrblitr2POM=spval
  allocate(this%ElmAllocmatMicrblitr2POMU(1:micpar%ndbiomcp));this%ElmAllocmatMicrblitr2POMU=spval
  allocate(this%RDOMEcoDmndPrev(1:jcplx));this%RDOMEcoDmndPrev=spval
  allocate(this%RAcetateEcoDmndPrev(1:jcplx));this%RAcetateEcoDmndPrev=spval
  allocate(this%RDOCUptkHeterPrev(1:micpar%NumHetetr1MicCmplx,1:jcplx))      ;this%RDOCUptkHeterPrev=0._r8  
  allocate(this%VLWatMicPM(NPH));this%VLWatMicPM=spval
  allocate(this%THETPM(NPH));this%THETPM=spval
  allocate(this%FILM(NPH));this%FILM=spval
  allocate(this%TortMicPM(NPH));this%TortMicPM=spval
  allocate(this%VLsoiAirPM(NPH));this%VLsoiAirPM=spval
  allocate(this%DiffusivitySolutEff(NPH));this%DiffusivitySolutEff=spval
  allocate(this%RAcetateUptkHeterPrev(1:NumHetetr1MicCmplx,1:jcplx))  ;this%RAcetateUptkHeterPrev=0._r8
  allocate(this%AttenfNH4HeterR(NumHetetr1MicCmplx,1:jcplx));this%AttenfNH4HeterR=spval
  allocate(this%AttenfNO3HeterR(NumHetetr1MicCmplx,1:jcplx));this%AttenfNO3HeterR=spval
  allocate(this%AttenfH2PO4HeterR(NumHetetr1MicCmplx,1:jcplx));this%AttenfH2PO4HeterR=spval
  allocate(this%AttenfH1PO4HeterR(NumHetetr1MicCmplx,1:jcplx));this%AttenfH1PO4HeterR=spval
  allocate(this%AttenfNH4AutorR(NumMicrobAutoTrophCmplx));this%AttenfNH4AutorR=spval
  allocate(this%AttenfNO3AutorR(NumMicrobAutoTrophCmplx));this%AttenfNO3AutorR=spval
  allocate(this%AttenfH2PO4AutorR(NumMicrobAutoTrophCmplx));this%AttenfH2PO4AutorR=spval
  allocate(this%AttenfH1PO4AutorR(NumMicrobAutoTrophCmplx));this%AttenfH1PO4AutorR=spval  
  end subroutine Init
!------------------------------------------------------------------------------------------


  subroutine Destruct(this)
  use abortutils, only : destroy
  implicit none
  class(micforctype) :: this

  call destroy(this%RAcetateUptkHeterPrev)
  call destroy(this%RDOCUptkHeterPrev)
  call destroy(this%VLWatMicPM)
  call destroy(this%THETPM)
  call destroy(this%FILM)
  call destroy(this%TortMicPM)
  call destroy(this%VLsoiAirPM)
  call destroy(this%DiffusivitySolutEff)
  call destroy(this%RDOMEcoDmndPrev)
  call destroy(this%RAcetateEcoDmndPrev)
  call destroy(this%ElmAllocmatMicrblitr2POM)
  call destroy(this%ElmAllocmatMicrblitr2POMU)
  call destroy(this%AttenfNH4HeterR)
  call destroy(this%AttenfNO3HeterR)
  call destroy(this%AttenfH2PO4HeterR)
  call destroy(this%AttenfH1PO4HeterR)
  call destroy(this%AttenfNH4AutorR)
  call destroy(this%AttenfNO3AutorR)
  call destroy(this%AttenfH2PO4AutorR)
  call destroy(this%AttenfH1PO4AutorR)  
  end subroutine Destruct

end module MicForcTypeMod
