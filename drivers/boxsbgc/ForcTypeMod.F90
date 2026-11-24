module ForcTypeMod
  use data_kind_mod     , only : r8 => DAT_KIND_R8
  use EcoSimConst
  use TracerPropMod
  use EcoSIMSolverPar
  use TracerIDMod
  use ChemTracerParsMod
  use MiniFuncMod
  use minimathmod, only : AZMAX1
implicit none

  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  logical :: first

  type, public :: forc_type
!===============================================================================
!primary variables
    real(r8) :: DLYR3       !soil layer thickness
    real(r8) :: AREA3       !cross section area
    real(r8) :: CEC         !
    real(r8) :: XCEC        !a variable derived from CEC and CNH4
    real(r8) :: AEC         !
    real(r8) :: XAEC        !a variable derived from AEC and CPO4
    real(r8) :: CFE         !
    real(r8) :: CCA
    real(r8) :: CMG
    real(r8) :: CNA
    real(r8) :: CKA
    real(r8) :: CSO4
    real(r8) :: CCL
    real(r8) :: CAL
    real(r8) :: ZMG
    real(r8) :: ZNA
    real(r8) :: ZKA
    real(r8) :: CALPO
    real(r8) :: CFEPO
    real(r8) :: CCAPD
    real(r8) :: CCAPH
    real(r8) :: CALOH
    real(r8) :: CFEOH
    real(r8) :: CCACO
    real(r8) :: CCASO
    real(r8) :: BKDS
    real(r8) :: ORGC        !total soil organic C [g d-2]
    real(r8) :: ATCS        !mean annual temperature for offset calculation
    real(r8) :: EHUM        !partitioning coefficient between humus and microbial residue, [], hour1.f
    real(r8) :: pH          !pH value
    real(r8) :: SoilMicPMassLayer        !mass of soil layer	Mg d-2
    real(r8) :: VLSoilPoreMicP        !volume of soil layer	m3 d-2
    real(r8) :: ZNH4S       !NH4 non-band micropore, [g d-2]
    real(r8) :: ZNO3S       !NO3 band micropore, [g d-2]
    real(r8) :: H2PO4       !PO4 non-band micropore, [g d-2]
    real(r8) :: H1PO4       !soil aqueous HPO4 content micropore non-band, [mol d-2]
    real(r8) :: ZNO2S       !NO2  non-band micropore, [g d-2]
    real(r8) :: VLSoilMicP      !micropore volume, [m3 d-2]
    real(r8) :: POROS     !soil porosity, [m3 m-3]
    real(r8) :: FieldCapacity        !water contents at field capacity, [m3 m-3]
    real(r8) :: SRP
    real(r8) :: WiltPoint
    real(r8) :: LOGPSIMX
    real(r8) :: LOGPSIMND
    real(r8) :: LOGPSIAtSat
    real(r8) :: PSISD
    real(r8) :: PSISE
    real(r8), allocatable :: ElmAllocmatMicrblitr2POM(:)  !allocation coefficient to humus fractions
    real(r8), allocatable :: CNOSC(:,:)
    real(r8), allocatable :: CPOSC(:,:)
    real(r8), allocatable :: SolidOM(:,:,:)
    real(r8), allocatable :: SolidOMAct(:,:)
    real(r8), allocatable :: OMBioResdu(:,:,:)
    real(r8), allocatable :: SorbedOM(:,:)
    real(r8), allocatable :: DOM(:,:)
    real(r8), allocatable :: mBiomeHeter(:,:,:)
    real(r8), allocatable :: mBiomeAutor(:,:)

!    real(r8) :: O2_irrig_conc        !surface irrigation  O2 concentration, [g m-3]
!    real(r8) :: O2_rain_conc        !precipitation  O2 concentration, [g m-3]
!    real(r8) :: Irrig2LitRSurf_col       !irrigation flux into surface litter, [m3 d-2 h-1]
!    real(r8) :: Rain2LitRSurf       !precipitation flux into surface litter, [m3 d-2 h-1]
!
    real(r8) :: TKS          !temperature in kelvin, [K]
    real(r8) :: THETW        !volumetric water content [m3 m-3]
!===============================================================================
!  derived variables
    real(r8) :: PARG
    real(r8) :: DCO2GQ
    real(r8) :: DCH4GQ
    real(r8) :: DOXYGQ
    real(r8) :: DZ2GGQ
    real(r8) :: DZ2OGQ
    real(r8) :: DNH3GQ
    real(r8) :: DH2GGQ

    real(r8) :: TempOffset      !offset for calculating temperature in Arrhenius curves, [oC]

    real(r8) :: THETY       !air-dry water content, [m3 m-3]
    real(r8) :: VLWatMicP        !soil micropore water content [m3 d-2]
    real(r8) :: VLsoiAirP        !soil air content, [m3 d-2]
    real(r8) :: VOLA        !total volume in micropores [m3 d-2]

    real(r8) :: VLNOB       !NO3 band volume fracrion, [-]
    real(r8) :: VLNO3       !NO3 non-band volume fraction,[-] := 1-VLNOB
    real(r8) :: VLNHB       !NH4 band volume fraction, [-]
    real(r8) :: VLNH4       !NH4 non-band volume fraction, [-]:=1-VLNHB
    real(r8) :: VLPOB       !PO4 band volume fracrion, [-]
    real(r8) :: VLPO4       !PO4 non-band volume fraction,[-]:=1-VLPOB

    real(r8) :: PSISoilMatricP       !soil micropore matric water potential, [MPa]
    real(r8) :: O2AquaDiffusvity       !aqueous O2 diffusivity, [m2 h-1], set in hour1
    real(r8) :: CLSGL       !aqueous CO2 diffusivity	[m2 h-1]
    real(r8) :: CQSGL       !aqueous CH4 diffusivity	m2 h-1
    real(r8) :: ZLSGL       !aqueous N2 diffusivity, [m2 h-1]
    real(r8) :: ZNSGL       !aqueous NH3 diffusivity, [m2 h-1]
    real(r8) :: ZVSGL       !aqueous N2O diffusivity, [m2 h-1]
    real(r8) :: HLSGL       !aqueous H2 diffusivity, [m2 h-1]

    real(r8) :: CGSGL       !gaseous CO2 diffusivity
    real(r8) :: CHSGL       !gaseous CH4 diffusivity
    real(r8) :: OGSGL       !gaseous O2 diffusivity
    real(r8) :: ZGSGL       !gaseous N2 diffusivity
    real(r8) :: Z2SGL       !gaseous N2O diffusivity
    real(r8) :: ZHSGL       !gaseous NH3 diffusivity
    real(r8) :: HGSGL       !gaseous H2 diffusivity

    real(r8) :: CO2E                      !initial atmospheric CO2 concentration, [umol mol-1]
    real(r8) :: OXYE                       !atmospheric O2 concentration, [umol mol-1]
    real(r8) :: Z2OE                       !atmospheric N2O concentration, [umol mol-1]
    real(r8) :: Z2GE                       !atmospheric N2 concentration, [umol mol-1]
    real(r8) :: ZNH3E                      !atmospheric NH3 concentration, [umol mol-1]
    real(r8) :: CH4E                       !atmospheric CH4 concentration, [umol mol-1]
    real(r8) :: H2GE                       !atmospheric H2 concentration, [umol mol-1]
    real(r8) :: CCO2E                     !initial atmospheric CO2 concentration, [g m-3]
    real(r8) :: COXYE                      !atmospheric O2 concentration, [g m-3]
    real(r8) :: CZ2OE                      !atmospheric N2O concentration, [g m-3]
    real(r8) :: CZ2GE                      !atmospheric N2 concentration, [g m-3]
    real(r8) :: CNH3E                      !atmospheric NH3 concentration, [g m-3]
    real(r8) :: CCH4E                      !atmospheric CH4 concentration, [g m-3]
    real(r8) :: CH2GE                      !atmospheric H2 concentration, [g m-3]


    real(r8) :: ZNH4B       !NH4 band micropore, [g d-2]
    real(r8) :: ZNO3B       !NO3 band micropore, [g d-2]
    real(r8) :: H2POB       !PO4 band micropore, [g d-2]
    real(r8) :: H1POB       !soil aqueous HPO4 content micropore band, [mol d-2]
    real(r8) :: ZNO2B       !NO2  band micropore, [g d-2]
    real(r8) :: CNH4B       !NH4 concentration band micropore	[g m-3], derived from ZNH4B
    real(r8) :: CNO3B       !NO3 concentration band micropore	[g m-3], derived from ZNO3B
    real(r8) :: CH2P4B      !aqueous PO4 concentration band	[g m-3], derived from H2POB
    real(r8) :: CH1P4B      !aqueous H1PO4 concentration band [g m-3], derived from H1POB
    real(r8) :: CNO2B       !aqueous HNO2 concentration band [g m-3], derived from ZNO2B

    real(r8) :: CNO2S       !NO2 concentration non-band micropore	[g m-3], derived from ZNO2S
    real(r8) :: CNH4S       !NH4 concentration non-band micropore	[g m-3], derived from ZNH4S
    real(r8) :: CNO3S       !NO3 concentration non-band micropore	[g m-3], derived from ZNO3S
    real(r8) :: CH2P4       !aqueous PO4 concentration non-band	[g m-3], derived from H2PO4
    real(r8) :: CH1P4       !aqueous H1PO4 concentration non-band [g m-3], derived from H1PO4
    real(r8) :: DiffusivitySolutEff        !coefficient for dissolution - volatilization, []
    real(r8) :: THETPM      !soil air-filled porosity, [m3 m-3]

!derived variables
    real(r8) :: FILM        !soil water film thickness , [m]
    real(r8) :: TortMicPM        !soil tortuosity, []

    real(r8) :: EPOC      !partitioning coefficient between POC and litter, [], hour1.f
    real(r8) :: CCO2S     !aqueous CO2 concentration micropore	[g m-3]
    real(r8) :: CZ2OS     !aqueous N2O concentration micropore	[g m-3]
    real(r8) :: COXYS     !aqueous O2 concentration micropore	[g m-3]
    real(r8) :: CZ2GS     !aqueous N2 concentration micropore	[g m-3]
    real(r8) :: COXYG     !gaseous O2 concentration	[g m-3]
    real(r8) :: CH2GS     !aqueous H2 concentration	[g m-3]
    real(r8) :: CCH4G     !gaseous CH4 concentration	[g m-3]
    real(r8) :: Z2OS      !aqueous N2O micropore, [g d-2]
    real(r8) :: OXYS      !aqueous O2  micropore	[g d-2]
    real(r8) :: H2GS      !aqueous H2 	[g d-2]
    real(r8) :: CH4S      !aqueous CO2  micropore	[g d-2]
    real(r8) :: CH4AquaSolubility     !solubility of CH4, [m3 m-3]
    real(r8) :: O2GSolubility     !solubility of O2, [m3 m-3]
    real(r8) :: SCO2L     !solubility of CO2, [m3 m-3]
    real(r8) :: SN2GL     !solubility of N2, [m3 m-3]
    real(r8) :: SN2OL     !solubility of N2O, [m3 m-3]
    real(r8) :: SNH3L     !solubility of NH3, [m3 m-3]
    real(r8) :: SH2GL     !solubility of H2, [m3 m-3]
    real(r8) :: ZNFN0     !initial nitrification inhibition activity
    real(r8) :: ZNFNI     !current nitrification inhibition activity

    real(r8) :: TScal4Difsvity      !temperature effect on aqueous diffusivity
    real(r8) :: ATKA      !mean annual air temperature [K]
 !litter layer
    real(r8) :: VLitR      !surface litter volume, [m3 d-2]
    real(r8) :: VWatLitRHoldCapcity    !surface litter water holding capacity, [m3 d-2]
 !non litter layer
    real(r8) :: RO2EcoDmndPrev       !total root + microbial O2 uptake from previous hour, [g d-2 h-1], updated in hour1
    real(r8) :: RN2OEcoUptkSoilPrev       !total root + microbial N2O uptake from previous hour, [g d-2 h-1]
    real(r8) :: RNO2EcoUptkSoilPrev       !total root + microbial NO2 uptake non-band from previous hour, [g d-2 h-1]
    real(r8) :: RNO2EcoUptkBandPrev       !total root + microbial NO2 uptake band from previous hour, [g d-2 h-1]
    real(r8) :: RNH4EcoDmndBandPrev       !total root + microbial NH4 uptake band from previous hour, [g d-2 h-1]
    real(r8) :: RNO3EcoDmndBandPrev       !total root + microbial NO3 uptake band from previous hour, [g d-2 h-1]
    real(r8) :: RH2PO4EcoDmndBandPrev       !total root + microbial PO4 uptake band from previous hour, [g d-2 h-1]
    real(r8) :: RH1PO4EcoDmndBandPrev       !HPO4 demand in band by all microbial, root, myco populations from previous hour
    real(r8) :: RNH4EcoDmndSoilPrev       !total root + microbial NH4 uptake non-band from previous hour, [g d-2 h-1]
    real(r8) :: RNO3EcoDmndSoilPrev       !total root + microbial NO3 uptake non-band from previous hour, [g d-2 h-1]
    real(r8) :: RH2PO4EcoDmndSoilPrev       !total root + microbial PO4 uptake non-band from previous hour, [g d-2 h-1]
    real(r8) :: RH1PO4EcoDmndSoilPrev       !HPO4 demand in non-band by all microbial, root, myco populations from previous hour
    real(r8) :: RO2GasXchangePrev       !net gaseous O2 flux, [g d-2 h-1], updated in redist.f
    real(r8) :: RCH4PhysexchPrev_vr       !net aqueous CH4 flux, [g d-2 h-1], updated in redist.f
    real(r8) :: RO2AquaXchangePrev       !net aqueous O2 flux from previous hour, [g d-2 h-1], updated in redist.f
    logical  :: disvolonly  !flag to only do dissolution/volatilization
  end type forc_type

  contains


!------------------------------------------------------------------------------------------

  subroutine ReadFORC(forc,fname)
  use ncdio_pio
  implicit none
  type(forc_type), intent(inout) :: forc
  character(len=*), intent(in) :: fname
  integer :: jcplx,ndbiomcp,nlbiomcp,NumMicbAFunGrupsPerCmplx
  integer :: NumMicbHFunGrupsPerCmplx,jsken,NumHetetr1MicCmplx,NumMicrobAutoTrophCmplx
  integer :: NumLiveHeterBioms,NumLiveAutoBioms
  integer :: NumPlantChemElms
  type(file_desc_t) :: ncf

  call ncd_pio_openfile(ncf, fname, ncd_nowrite)

  jcplx=get_dim_len(ncf,'jcplx')
  jsken =get_dim_len(ncf,'jsken')
  NumHetetr1MicCmplx=get_dim_len(ncf,'NumHetetr1MicCmplx')
  NumMicrobAutoTrophCmplx=get_dim_len(ncf,'NumMicrobAutoTrophCmplx')
  NumLiveHeterBioms = get_dim_len(ncf,'NumLiveHeterBioms')
  NumLiveAutoBioms = get_dim_len(ncf,'NumLiveAutoBioms')
  NumPlantChemElms=get_dim_len(ncf,'element')
  nlbiomcp=get_dim_len(ncf,'nlbiomcp')
  ndbiomcp=get_dim_len(ncf,'ndbiomcp')
  NumMicbAFunGrupsPerCmplx    =get_dim_len(ncf,'NumMicbAFunGrupsPerCmplx')
  NumMicbHFunGrupsPerCmplx    =get_dim_len(ncf,'NumMicbHFunGrupsPerCmplx')
  allocate(forc%mBiomeHeter(NumPlantChemElms,NumLiveHeterBioms,1:jcplx))
  allocate(forc%DOM(idom_beg:idom_end,1:jcplx))
  allocate(forc%SolidOM(1:NumPlantChemElms,jsken,1:jcplx))
  allocate(forc%SolidOMAct(jsken,1:jcplx))
  allocate(forc%OMBioResdu(1:NumPlantChemElms,ndbiomcp,1:jcplx))
  allocate(forc%mBiomeAutor(NumPlantChemElms,NumLiveAutoBioms))
  allocate(forc%SorbedOM(idom_beg:idom_end,1:jcplx))
  allocate(forc%ElmAllocmatMicrblitr2POM(ndbiomcp))
  allocate(forc%CNOSC(1:jsken,1:jcplx))
  allocate(forc%CPOSC(1:jsken,1:jcplx))
  call ncd_getvar(ncf,'ZNH4S',forc%ZNH4S)
  call ncd_getvar(ncf,'ZNO3S',forc%ZNO3S)
  call ncd_getvar(ncf,'ZNO2S',forc%ZNO2S)
  call ncd_getvar(ncf,'H2PO4',forc%H2PO4)
  call ncd_getvar(ncf,'H1PO4',forc%H1PO4)
  call ncd_getvar(ncf,'CNOSC',forc%CNOSC)
  call ncd_getvar(ncf,'CPOSC',forc%CPOSC)
  call ncd_getvar(ncf,'pH',forc%PH)
  call ncd_getvar(ncf,'VLSoilPoreMicP',forc%VLSoilPoreMicP)
  call ncd_getvar(ncf,'ORGC',forc%ORGC)
  call ncd_getvar(ncf,'ElmAllocmatMicrblitr2POM',forc%ElmAllocmatMicrblitr2POM)
  call ncd_getvar(ncf,'VLSoilMicP',forc%VLSoilMicP)
  call ncd_getvar(ncf,'BKVL',forc%SoilMicPMassLayer)
  call ncd_getvar(ncf,'POROS',forc%POROS)
  call ncd_getvar(ncf,'FC',forc%FieldCapacity)
  call ncd_getvar(ncf,'WP',forc%WiltPoint)
  call ncd_getvar(ncf,'SRP',forc%SRP)
  call ncd_getvar(ncf,'EHUM',forc%EHUM)
  call ncd_getvar(ncf,'EPOC',forc%EPOC)
  call ncd_getvar(ncf,'DLYR3',forc%DLYR3)
  call ncd_getvar(ncf,'CEC',forc%CEC)
  call ncd_getvar(ncf,'AEC',forc%AEC)
  call ncd_getvar(ncf,'CFE',forc%CFE)
  call ncd_getvar(ncf,'CCA',forc%CCA)
  call ncd_getvar(ncf,'CMG',forc%CMG)
  call ncd_getvar(ncf,'CNA',forc%CNA)
  call ncd_getvar(ncf,'CKA',forc%CKA)
  call ncd_getvar(ncf,'CSO4',forc%CSO4)
  call ncd_getvar(ncf,'CCL',forc%CCL)
  call ncd_getvar(ncf,'ZMG',forc%ZMG)
  call ncd_getvar(ncf,'ZNA',forc%ZNA)
  call ncd_getvar(ncf,'ZKA',forc%ZKA)
  call ncd_getvar(ncf,'CAL',forc%CAL)
  call ncd_getvar(ncf,'CALPO',forc%CALPO)
  call ncd_getvar(ncf,'CFEPO',forc%CFEPO)
  call ncd_getvar(ncf,'CCAPD',forc%CCAPD)
  call ncd_getvar(ncf,'CCAPH',forc%CCAPH)
  call ncd_getvar(ncf,'CALOH',forc%CALOH)
  call ncd_getvar(ncf,'CFEOH',forc%CFEOH)
  call ncd_getvar(ncf,'CCACO',forc%CCACO)
  call ncd_getvar(ncf,'CCASO',forc%CCASO)
  call ncd_getvar(ncf,'BKDS',forc%BKDS)
  call ncd_getvar(ncf,'ATCS',forc%ATCS)
  call ncd_getvar(ncf,'mBiomeHeter',forc%mBiomeHeter)
  call ncd_getvar(ncf,'mBiomeAutor',forc%mBiomeAutor)
  call ncd_getvar(ncf,'OSM',forc%SolidOM)
  call ncd_getvar(ncf,'OSA',forc%SolidOMAct)
  call ncd_getvar(ncf,'ORM',forc%OMBioResdu)
  call ncd_getvar(ncf,'OHM',forc%SorbedOM)
  call ncd_getvar(ncf,'DOM',forc%DOM)
  call ncd_getvar(ncf,'THETY',forc%THETY)
  call ncd_getvar(ncf,'PSIMX',forc%LOGPSIMX)
  call ncd_getvar(ncf,'PSIMD',forc%LOGPSIMND)
  call ncd_getvar(ncf,'PSIMS',forc%LOGPSIAtSat)
  call ncd_getvar(ncf,'PSISD',forc%PSISD)
  call ncd_getvar(ncf,'PSISE',forc%PSISE)
  call ncd_getvar(ncf,'ATKA', forc%ATKA)
  call ncd_pio_closefile(ncf)

  forc%TempOffset=f_microbe_TOffset(forc%ATCS)
  forc%VLNOB=0._r8
  forc%VLNO3=1._r8
  forc%VLNHB=0._r8
  forc%VLNH4=1._r8
  forc%VLPOB=0._r8
  forc%VLPO4=1._r8
  forc%ZNFN0=0._r8
  forc%ZNFNI=0._r8
  forc%VOLA=forc%POROS*forc%VLSoilMicP
  forc%AREA3=1._r8
  forc%CNH4B=0._r8       !NH4 concentration band micropore	[g m-3], derived from ZNH4B
  forc%ZNH4B=0._r8       !NH4 band micropore, [g d-2]
  forc%CNO3B=0._r8       !NO3 concentration band micropore	[g m-3], derived from ZNO3B
  forc%ZNO3B=0._r8       !NO3 band micropore, [g d-2]
  forc%CH2P4B=0._r8      !aqueous PO4 concentration band	[g m-3], derived from H2POB
  forc%H2POB =0._r8      !PO4 band micropore, [g d-2]
  forc%CH1P4B=0._r8      !aqueous H1PO4 concentration band [g m-3], derived from H1POB
  forc%H1POB =0._r8      !soil aqueous HPO4 content micropore band, [mol d-2]
  forc%CNO2B =0._r8      !aqueous HNO2 concentration band [g m-3], derived from ZNO2B
  forc%ZNO2B =0._r8      !NO2  band micropore, [g d-2]




  first=.true.
  end subroutine ReadForc
!------------------------------------------------------------------------------------------

  subroutine UpdateFORC(forc, forctype)
  use MiniFuncMod
  implicit none
  integer, intent(in) :: forctype  ! 0: (transient), 1: T const, 2: water const, 3: T and water const
  type(forc_type), target, intent(inout) :: forc

  real(r8) :: FCL, LOGWiltPoint,PSL,FCD,PSD
  real(r8) :: scalar, Z3S, TCS
  real(r8) :: TFACG
  real(r8) :: DFLG2
  real(r8) :: PARG
  real(r8) :: PARGCO
  real(r8) :: PARGCH
  real(r8) :: PARGOX
  real(r8) :: PARGNG
  real(r8) :: PARGN2
  real(r8) :: PARGN3
  real(r8) :: PARGH2

  real(r8) :: DCO2G
  real(r8) :: DCH4G
  real(r8) :: DOXYG
  real(r8) :: DZ2GG
  real(r8) :: DZ2OG
  real(r8) :: DNH3G
  real(r8) :: DH2GG
  real(r8) :: PARGM

  real(r8), parameter :: DTHETW=1.0E-06_r8
  real(r8) :: XNPD
  !since the model if configured for incubation
  !the primary variable forcing is temperature, and moisture
  associate(                     &
     LOGPSIMX    =>  forc%LOGPSIMX   , &
     LOGPSIMND    =>  forc%LOGPSIMND   , &
     PSISD    =>  forc%PSISD   , &
     LOGPSIAtSat    =>  forc%LOGPSIAtSat   , &
     PSISE    =>  forc%PSISE   , &
     SRP      =>  forc%SRP     , &
    THETW => forc%THETW        , &  !relative saturation
    TKS   => forc%TKS            &
  )

  TKS=298._r8
  THETW=0.65_r8
  if (first .or. forctype<=1)then
!  variable moisture
    FCL=log(forc%FieldCapacity)
    LOGWiltPoint=log(forc%WiltPoint)
    PSL=log(forc%POROS)
    PSD=PSL-FCL
    FCD=FCL-LOGWiltPoint
    IF(THETW.LT.forc%FieldCapacity)THEN
      forc%PSISoilMatricP=AMAX1(PSIHY,-EXP(LOGPSIMX+((FCL-LOG(THETW))/FCD*LOGPSIMND)))
    ELSEIF(THETW.LT.forc%POROS-DTHETW)THEN
      forc%PSISoilMatricP=-EXP(LOGPSIAtSat+(((PSL-LOG(THETW))/PSD)**SRP*PSISD))
    ELSE
      forc%PSISoilMatricP=PSISE
    ENDIF
    forc%FILM=FilmThickness(forc%PSISoilMatricP)
    forc%TortMicPM=TortMicporew(THETW)
    forc%THETPM=1._r8-THETW
    forc%VLWatMicP=forc%THETPM*forc%POROS*forc%VLSoilMicP
    forc%VLsoiAirP=forc%VOLA-forc%VLWatMicP
  endif

  if(first .or. forctype ==0 .or. forctype==2)then

    TCS=TKS-TC2K
    forc%CH4AquaSolubility=gas_solubility(idg_CH4,TCS)
    forc%O2GSolubility=gas_solubility(idg_O2, TCS)
    forc%SCO2L=gas_solubility(idg_CO2, TCS)     !solubility of CO2, [m3 m-3]
    forc%SN2GL=gas_solubility(idg_N2, TCS)      !solubility of N2, [m3 m-3]
    forc%SN2OL=gas_solubility(idg_N2O, TCS)     !solubility of N2O, [m3 m-3]
    forc%SNH3L=gas_solubility(idg_NH3, TCS)     !solubility of NH3, [m3 m-3]
    forc%SH2GL=gas_solubility(idg_H2, TCS)       !solubility of H2, [m3 m-3]


    forc%CCO2E=forc%CO2E*5.36E-04_r8*Tref/TKS    ![gC/m3]
    forc%CCH4E=forc%CH4E*5.36E-04_r8*Tref/TKS    ![gC/m3]
    forc%COXYE=forc%OXYE*1.43E-03_r8*Tref/TKS    ![gO/m3]
    forc%CZ2GE=forc%Z2GE*1.25E-03_r8*Tref/TKS    ![gN/m3]
    forc%CZ2OE=forc%Z2OE*1.25E-03_r8*Tref/TKS    ![gN/m3]
    forc%CNH3E=forc%ZNH3E*6.25E-04_r8*Tref/TKS   ![gN/m3]
    forc%CH2GE=forc%H2GE*8.92E-05_r8*Tref/TKS    ![gH/m3]

    !variable temperature
    forc%TScal4Difsvity=TEFAQUDIF(TKS)
    forc%O2AquaDiffusvity=OLSG*forc%TScal4Difsvity       !aqueous O2 diffusivity [m2 h-1]
    forc%CLSGL=CLSG*forc%TScal4Difsvity       !aqueous CO2 diffusivity	[m2 h-1]
    forc%CQSGL=CQSG*forc%TScal4Difsvity       !aqueous CH4 diffusivity	m2 h-1
    forc%ZLSGL=ZLSG*forc%TScal4Difsvity       !aqueous N2 diffusivity, [m2 h-1]
    forc%ZNSGL=ZNSG*forc%TScal4Difsvity       !aqueous NH3 diffusivity, [m2 h-1]
    forc%ZVSGL=ZVSG*forc%TScal4Difsvity       !aqueous N2O diffusivity, [m2 h-1]
    forc%HLSGL=HLSG*forc%TScal4Difsvity       !aqueous H2 diffusivity, [m2 h-1]

  endif

  if (first .or. forctype <=2)then
    Z3S=forc%FieldCapacity/forc%POROS
    XNPD=600.0_r8*dts_gas
    scalar=forc%TScal4Difsvity*XNPD
    forc%DiffusivitySolutEff=fDiffusivitySolutEff(scalar,THETW,Z3S)

    DFLG2=2.0_r8*AZMAX1(forc%THETPM)*POROQ &
      *forc%THETPM/forc%POROS*forc%AREA3/forc%DLYR3
    TFACG=TEFGASDIF(forc%TKS)
    forc%CGSGL=CGSG*TFACG*DFLG2
    forc%CHSGL=CHSG*TFACG*DFLG2
    forc%OGSGL=OGSG*TFACG*DFLG2
    forc%ZGSGL=ZGSG*TFACG*DFLG2
    forc%Z2SGL=Z2SG*TFACG*DFLG2
    forc%ZHSGL=ZHSG*TFACG*DFLG2
    forc%HGSGL=HGSG*TFACG*DFLG2

    forc%PARG=forc%AREA3*dts_HeatWatTP/(0.0139_r8+1.39E-03_r8)
    PARGM=forc%PARG*dt_GasCyc
    !GASEOUS BOUNDARY LAYER CONDUCTANCES
    PARGCO=PARGM*0.74_r8
    PARGCH=PARGM*1.04_r8
    PARGOX=PARGM*0.83_r8
    PARGNG=PARGM*0.86_r8
    PARGN2=PARGM*0.74_r8
    PARGN3=PARGM*1.02_r8
    PARGH2=PARGM*2.08_r8

    DCO2G=forc%CGSGL*dts_gas
    DCH4G=forc%CHSGL*dts_gas
    DOXYG=forc%OGSGL*dts_gas
    DZ2GG=forc%ZGSGL*dts_gas
    DZ2OG=forc%Z2SGL*dts_gas
    DNH3G=forc%ZHSGL*dts_gas
    DH2GG=forc%HGSGL*dts_gas


    forc%DCO2GQ=DCO2G*PARGCO/(DCO2G+PARGCO)
    forc%DCH4GQ=DCH4G*PARGCH/(DCH4G+PARGCH)
    forc%DOXYGQ=DOXYG*PARGOX/(DOXYG+PARGOX)
    forc%DZ2GGQ=DZ2GG*PARGNG/(DZ2GG+PARGNG)
    forc%DZ2OGQ=DZ2OG*PARGN2/(DZ2OG+PARGN2)
    forc%DNH3GQ=DNH3G*PARGN3/(DNH3G+PARGN3)
    forc%DH2GGQ=DH2GG*PARGH2/(DH2GG+PARGH2)

  endif

  forc%RO2EcoDmndPrev=0._r8
  forc%RN2OEcoUptkSoilPrev=0._r8
  forc%RNO2EcoUptkSoilPrev=0._r8
  forc%RNO2EcoUptkBandPrev=0._r8
  forc%RNH4EcoDmndBandPrev=0._r8
  forc%RNO3EcoDmndBandPrev=0._r8
  forc%RH2PO4EcoDmndBandPrev=0._r8
  forc%RH1PO4EcoDmndBandPrev=0._r8
  forc%RNH4EcoDmndSoilPrev=0._r8
  forc%RNO3EcoDmndSoilPrev=0._r8
  forc%RH2PO4EcoDmndSoilPrev=0._r8
  forc%RH1PO4EcoDmndSoilPrev=0._r8
  forc%RO2GasXchangePrev=0._r8
  forc%RCH4PhysexchPrev_vr=0._r8
  forc%RO2AquaXchangePrev=0._r8

  first=.false.
  end associate
  end subroutine UpdateForc

end module ForcTypeMod
