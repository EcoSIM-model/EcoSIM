module MicrobeDiagTypes
  ! !PUBLIC TYPES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL       
  USE abortutils, only : destroy
  use TracerIDMod
  use EcoSiMParDataMod, only : micpar
  implicit none
  save
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

! accumulative flux diagnostics
type, public :: Cumlate_Flux_Diag_type
! fraction diagnostics, used for substrate competition/uptake
    real(r8) :: TFNH4B           !total releative demand of NH4 in banded soil by microbes
    real(r8) :: TFNO3B           !total releative demand of NO3 in banded soil by microbes
    real(r8) :: TFNO2B           !total releative demand of NO2 in banded soil by microbes
    real(r8) :: TFP14B           !total releative demand of HPO4 in banded soil by microbes
    real(r8) :: TFPO4B           !total releative demand of H2PO4 in banded soil by microbes
    real(r8) :: tCH4ProdAceto    !acetoclastic CH4 production 
    real(r8) :: tCH4ProdH2       !Hydrogenotrophic methane production
    real(r8) :: TFOQC            !total relative DOC uptake by microbes
    real(r8) :: TFOQA            !total releative acetate uptake by microbes
    real(r8) :: TFOXYX           !total relative O2 demand by microbes
    real(r8) :: TFNH4X           !total relative demand of NH4 in non-banded soil by Microbes
    real(r8) :: TFNO3X           !total relative demand of NO3 in non-banded soil by Microbes
    real(r8) :: TFNO2X           !total relative demand of NO2 in non-banded soil by Microbes
    real(r8) :: TFN2OX           !total relative demand of N2O in non-banded soil by Microbes
    real(r8) :: TFP14X           !total relative demand of H1PO4 in non-banded soil by Microbes
    real(r8) :: TFPO4X           !total relative demand of H2PO4 in non-banded soil by Microbes    
    real(r8) :: tRespGrossHeter    !total gross respiration by heterotrophs
    real(r8) :: tRespGrossHeterUlm !total gross respiration by heterotrophs unlimited by O2
!fluxes
    real(r8) :: tCResp4H2Prod    !total respiration flux to support H2 production, [gC d-2 h-1]  
    real(r8) :: tRNH4MicrbImobilSoil
    real(r8) :: tRH2PO4MicrbImobilSoil
    real(r8) :: tRNO3MicrbImobilSoil
    real(r8) :: tRH1PO4MicrbImobilSoil
    real(r8) :: tRNH4MicrbImobilBand
    real(r8) :: tRNO3MicrbImobilBand
    real(r8) :: tRH2PO4MicrbImobilBand
    real(r8) :: tRH1PO4MicrbImobilBand
    real(r8) :: tNO2OxiAuto      !total NO2 oxidation into NO3 by nitrifiers
    real(r8) :: tRO2UptkHeterG   !total oxygen uptake by heterotrophs
    real(r8) :: tRO2DmndHeterG   !total oxygen demand by heterotrophs
    real(r8) :: tRCO2MicrbProd
    real(r8) :: tRCH4MicrbProd
    real(r8) :: tRNOxMicrbRedux
    real(r8) :: tRCO2GrothAutor       !CO2 taken up for autotrophic biomass growth
    real(r8) :: TProdH2               !H2 production by fermenters
    real(r8) :: tRO2MicrbUptk         !O2 uptake by microbes
    real(r8) :: TReduxNO3Soil         !NO3 reduction in non-band soil by denitrifiers
    real(r8) :: TReduxNO3Band         !NO3 reduction in banded soil by denitrifiers
    real(r8) :: TReduxNO2toN2OSoil         !NO2 reduction in non-band soil, by ammonia oxidizers, and denitrifiers
    real(r8) :: TReduxNO2toN2OBand         !NO2 reduction in banded soil, by ammonia oxidizers, and denitrifiers
    real(r8) :: TDeniReduxNO2Band
    real(r8) :: TDeniReduxNO2Soil
    real(r8) :: TNitNO2Redux2N2OBand
    real(r8) :: TNitNO2Redux2N2OSoil

    real(r8) :: TReduxN2OtoN2             !N2O reduction by detnitrifiers
    real(r8) :: TFixN2                !N2 fixation by aerobic and anaerobic N2 fixers
    real(r8) :: tCH4OxiAero           !aerobic CH4 oxidation
    real(r8) :: tCH4OxiANMO           !anaerobic CH4 oxidation
    real(r8) :: tRNH3Oxi              !NH3 oxidation by nitrifiers, soil+band

    real(r8) :: RNO2ReduxSoilChemo
    real(r8) :: RNO2ReduxBandChemo
    real(r8) :: RN2OProdSoilChemo
    real(r8) :: RN2OProdBandChemo
    real(r8) :: RNO3ProdSoilChemo
    real(r8) :: RNO3ProdBandChemo
    real(r8) :: RNO2ReduxChemo             !DON-N production from chemical nitrite reduction (band + nonband)

  contains
    procedure, public :: ZeroOut => nit_aqmf_diag    
  end type Cumlate_Flux_Diag_type

  type, public :: Microbe_State_type
  real(r8),allocatable :: rCNBiomeActHeter(:,:,:)  !Nutrient to carbon ratio for active heterotrophs
  real(r8),allocatable :: OMActHeter(:,:)

  real(r8),allocatable :: OMC2(:,:)
  real(r8),allocatable :: GrowthEnvScalHeter(:,:)
  real(r8),allocatable :: TempMaintRHeter(:,:)
  real(r8),allocatable :: OMN2(:,:)
  real(r8),allocatable :: FOM2(:,:)
  real(r8),allocatable :: OxyLimterHeter(:,:)

  real(r8),allocatable :: FCN(:,:)
  real(r8),allocatable :: FCP(:,:)
  real(r8),allocatable :: FBiomStoiScalarHeter(:,:)
  real(r8),allocatable :: FracOMActHeter(:,:)
  real(r8),allocatable :: FracNO2ReduxHeter(:,:)
  real(r8),allocatable :: FracHeterBiomOfActK(:,:)

  real(r8),allocatable :: rCNBiomeActAutor(:,:)   !nutrient to carbon ratio for active autotrophs
  real(r8),allocatable :: OMActAutor(:)
  real(r8),allocatable :: FracOMActAutor(:)
  real(r8),allocatable :: FracNO2XupAutor(:)
  real(r8),allocatable :: FracAutorBiomOfActK(:)
  real(r8),allocatable :: OMC2Autor(:)
  real(r8),allocatable :: GrowthEnvScalAutor(:)
  real(r8),allocatable :: TSensMaintRAutor(:)
  real(r8),allocatable :: TSensGroAutor(:)          !temperature sensitivity of autotrophic microbes
  real(r8),allocatable :: WSensGroAutor(:)          !moisture sensitivity of autotrophic microbes
  real(r8),allocatable :: TSensGroHeter(:,:)        !temperature sensitivity of microbial growth
  real(r8),allocatable :: WSensGroHeter(:,:)        !moisture sensitivity of microbial growth
  real(r8),allocatable :: OMN2Autor(:)
  real(r8),allocatable :: FOM2Autor(:)
  real(r8),allocatable :: fLimO2Autor(:)
  real(r8),allocatable :: FCNAutor(:)
  real(r8),allocatable :: FCPAutor(:)
  real(r8),allocatable :: FBiomNutStoiScalAutor(:)

  contains
   procedure, public :: Init => nit_mics_init
   procedure, public :: Destroy => nit_mics_destroy
  end type Microbe_State_type

  type, public :: Microbe_Flux_type
  !fluxes
  real(r8) :: RTotNH3OxidSoilAutor
  real(r8) :: RTotNH3OxidBandAutor  
! allocatable fluxes
  real(r8),allocatable :: DOMUptk4GrothHeter(:,:,:)
  real(r8),allocatable :: RMetabDOCUptkHeter(:,:)
  real(r8),allocatable :: RMetabAcetUptkHeter(:,:)
  real(r8),allocatable :: NonstX2stBiomHeter(:,:,:,:)  !nonstructural biomass export to labile and structural biomass, [g d-2 h-1]
  real(r8),allocatable :: RO2UptkHeter(:,:)
  real(r8),allocatable :: RO2UptkHeterG(:)      !complex summed O2 uptake 
  real(r8),allocatable :: Resp4NFixHeter(:,:)
  real(r8),allocatable :: RespGrossHeter(:,:)
  real(r8),allocatable :: RGrowthRespHeter(:,:)
  real(r8),allocatable :: RMaintRespHeter(:,:)
  real(r8),allocatable :: RO2Dmnd4RespHeter(:,:)
  real(r8),allocatable :: RO2DmndHeter(:,:)
  real(r8),allocatable :: RO2DmndHeterG(:)   !complex summed O2 demand
  real(r8),allocatable :: RO2Uptk4RespHeter(:,:)
  real(r8),allocatable :: RNO3ReduxHeterSoil(:,:)
  real(r8),allocatable :: RNO3ReduxHeterBand(:,:)
  real(r8),allocatable :: RNO2ReduxHeterSoil(:,:)
  real(r8),allocatable :: RNO2ReduxHeterBand(:,:)
  real(r8),allocatable :: RN2OReduxHeter(:,:)
  real(r8),allocatable :: RNOxDOCReduxRespDenitLim(:,:) !DOC respiration for denitrifcation
  real(r8),allocatable :: RNOxAcetReduxRespDenitLim(:,:) !Acetate respiration for denitrifcation  
  real(r8),allocatable :: RMaintDmndHeter(:,:,:)
  real(r8),allocatable :: RNH4imobilSoilHeter(:,:)
  real(r8),allocatable :: RNO3imobilSoilHeter(:,:)
  real(r8),allocatable :: RH2PO4imobilSoilHeter(:,:)
  real(r8),allocatable :: RNH4imobilBandHeter(:,:)
  real(r8),allocatable :: RNO3imobilBandHeter(:,:)
  real(r8),allocatable :: RH2PO4imobilBandHeter(:,:)
  real(r8),allocatable :: RkillLitfalOMHeter(:,:,:,:)
  real(r8),allocatable :: RkillLitrfal2HumOMHeter(:,:,:,:)
  real(r8),allocatable :: RkillLitrfal2ResduOMHeter(:,:,:,:)
  real(r8),allocatable :: RH2ProdHeter(:,:)
  real(r8),allocatable :: RMaintDefcitLitrfalOMHeter(:,:,:,:)
  real(r8),allocatable :: RMaintDefLitrfal2HumOMHeter(:,:,:,:)
  real(r8),allocatable :: RMaintDefLitrfal2ResduOMHeter(:,:,:,:)
  real(r8),allocatable :: RCCMEheter(:,:,:)
  real(r8),allocatable :: RCCMEAutor(:,:,:)
  real(r8),allocatable :: RN2FixHeter(:,:)
  real(r8),allocatable :: RKillOMHeter(:,:,:,:)
  real(r8),allocatable :: RkillRecycOMHeter(:,:,:,:)
  real(r8),allocatable :: RMaintDefcitKillOMHeter(:,:,:,:)
  real(r8),allocatable :: RMaintDefcitRecycOMHeter(:,:,:,:)
  real(r8),allocatable :: RNH4imobilLitrHeter(:,:)
  real(r8),allocatable :: RNO3imobilLitrHeter(:,:)
  real(r8),allocatable :: RH2PO4imobilLitrHeter(:,:)
  real(r8),allocatable :: RNOxReduxRespDenitUlm(:,:)
  real(r8),allocatable :: ROQC4HeterMicrobAct(:,:)
  real(r8),allocatable :: RCO2ProdHeter(:,:)
  real(r8),allocatable :: RCO2FixCyano(:,:)
  real(r8),allocatable :: RAcetateProdHeter(:,:)
  real(r8),allocatable :: RCH4ProdHeter(:,:)
  real(r8),allocatable :: RSMetaOxidSoilAutor(:)
  real(r8),allocatable :: RSMetaOxidBandAutor(:)
  real(r8),allocatable :: ECHZAutor(:)
  real(r8),allocatable :: ECHZHeter(:,:)
  real(r8),allocatable :: RGOCP(:,:)
  real(r8),allocatable :: RGOAP(:,:)
  real(r8),allocatable :: FOQC(:,:)
  real(r8),allocatable :: FOQA(:,:)
  REAL(R8),allocatable :: FGOCP(:,:)
  REAL(R8),allocatable :: FGOAP(:,:)  
  real(r8),allocatable :: fPhotoR(:,:)
  real(r8),allocatable :: XferBiomeHeterK(:,:,:,:)
  real(r8),allocatable :: RH1PO4imobilSoilHeter(:,:)
  real(r8),allocatable :: RH1PO4imobilBandHeter(:,:)
  real(r8),allocatable :: RH1PO4imobilLitrHeter(:,:)

  real(r8),allocatable :: RO2UptkAutor(:)
  real(r8),allocatable :: Resp4NFixAutor(:)
  real(r8),allocatable :: RespGrossAutor(:)
  real(r8),allocatable :: RO2Dmnd4GrossRespAutor(:)
  real(r8),allocatable :: RO2Uptk4RespAutor(:)
  real(r8),allocatable :: RNO3UptkAutor(:)
  real(r8),allocatable :: RNOxReduxAutorSoil(:)
  real(r8),allocatable :: RNOxReduxAutorBand(:)
  real(r8),allocatable :: RNOxReduxRespAutorLim(:)
  real(r8),allocatable :: RMaintDmndAutor(:,:)
  real(r8),allocatable :: RNH4TransfSoilAutor(:)
  real(r8),allocatable :: RNO3TransfSoilAutor(:)
  real(r8),allocatable :: RH2PO4TransfSoilAutor(:)
  real(r8),allocatable :: RNH4TransfBandAutor(:)
  real(r8),allocatable :: RNO3TransfBandAutor(:)
  real(r8),allocatable :: RH2PO4TransfBandAutor(:)
  real(r8),allocatable :: RkillLitfalOMAutor(:,:,:)
  real(r8),allocatable :: RkillLitrfal2HumOMAutor(:,:,:)
  real(r8),allocatable :: RkillLitrfal2ResduOMAutor(:,:,:)
  real(r8),allocatable :: DOMuptk4GrothAutor(:,:)            !C-uptake (CH4 or CO2) flux uptake by autotrophs, [g d-2 h-1]  
  real(r8),allocatable :: RMaintDefcitLitrfalOMAutor(:,:,:)
  real(r8),allocatable :: RMaintDefLitrfal2HumOMAutor(:,:,:)
  real(r8),allocatable :: RMaintDefLitrfal2ResduOMAutor(:,:,:)
  real(r8),allocatable :: RN2FixAutor(:)
  real(r8),allocatable :: RKillOMAutor(:,:,:)
  real(r8),allocatable :: RkillRecycOMAutor(:,:,:)
  real(r8),allocatable :: RMaintDefcitKillOMAutor(:,:,:)
  real(r8),allocatable :: RMaintDefcitRecycOMAutor(:,:,:)
  real(r8),allocatable :: RNH4TransfLitrAutor(:)
  real(r8),allocatable :: RNO3TransfLitrAutor(:)
  real(r8),allocatable :: RH2PO4TransfLitrAutor(:)
  real(r8),allocatable :: NonstX2stBiomAutor(:,:,:)
  real(r8),allocatable :: RCO2ProdAutor(:)    !gross CO2 production from autotrophs
  real(r8),allocatable :: RCO2XumpAutor(:)    !gross CO2 consumption from autotrophs
  real(r8),allocatable :: ROMProdCO2Autor(:)  !CO2 production from biomass C from autotrophs
  real(r8),allocatable :: RGrowthCAutor(:)    !gross C biomss growth flux from substrate uptake
  real(r8),allocatable :: RCH4ProdAutor(:)    !gross CH4 production
  real(r8),allocatable :: RH1PO4TransfSoilAutor(:)
  real(r8),allocatable :: RH1PO4TransfBandAutor(:)
  real(r8),allocatable :: RH1PO4TransfLitrAutor(:)

  contains
    procedure, public :: Init => nit_micf_init
    procedure, public :: ZeroOut => nit_micf_zero
    procedure, public :: destroy => nit_micf_destroy
  end type Microbe_Flux_type

  type, public :: OMCplx_Flux_type
    real(r8),allocatable :: RHydlysSolidOM(:,:,:)
    real(r8),allocatable :: RHumifySolidOM(:,:,:)
    real(r8),allocatable :: RDecmpProdDOM(:,:,:)
    real(r8),allocatable :: RHydlysBioResduOM(:,:,:)
    real(r8),allocatable :: RHydlysSorptOM(:,:)
    real(r8),allocatable :: RDOMSorp(:,:)                !rate of DOM adsorption
    real(r8),allocatable :: TDOMUptkHeter(:,:)
    real(r8),allocatable :: XferRespHeterK(:)
    real(r8),allocatable :: XferDOMK(:,:)
  contains
    procedure, public :: Init => nit_omcplxf_init
    procedure, public :: ZeroOut => nit_omcplxf_zero
    procedure, public :: Destroy => nit_omcplxf_destroy
  end type OMCplx_Flux_type

  type, public :: OMCplx_State_type
    real(r8),allocatable :: TOMEAutoK(:)
    real(r8),allocatable :: BulkSOMC(:)
    real(r8),allocatable :: FOCA(:)
    real(r8),allocatable :: FOAA(:)
    real(r8),allocatable :: rCNDOM(:)
    real(r8),allocatable :: rCPDOM(:)
    real(r8),allocatable :: OMBioResduK(:)
    real(r8),allocatable :: SolidOMCK(:)
    real(r8),allocatable :: SolidOMActK(:)
    real(r8),allocatable :: tMaxNActMicrbK(:)
    real(r8),allocatable :: tMaxPActMicrbK(:)
    real(r8),allocatable :: CDOM(:,:)            !dom concentration
  contains
    procedure, public :: Init => nit_omcplxs_init
    procedure, public :: ZeroOut => nit_omcplxs_zero
    procedure, public :: Destroy => nit_omcplxs_destroy
  end type OMCplx_State_type

  type, public :: Microbe_Diag_type
  real(r8) :: H1P4T
  real(r8) :: H2P4T
  real(r8) :: RH2UptkAutor
  real(r8) :: ThetaLitr
  real(r8) :: ThetaZ
  real(r8) :: TOMBioResdu
  real(r8) :: TotActMicrobiom
  real(r8) :: TotBiomNO2Consumers
  real(r8) :: TSensGrowth
  real(r8) :: TSensMaintR
  real(r8) :: VOLWZ
  real(r8) :: WatStressMicb
  real(r8) :: XCO2
  real(r8) :: ZNH4T
  real(r8) :: ZNO3T
  real(r8) :: ZNO2T
  real(r8),allocatable :: TOMEK(:,:)                     !total elemental biomass in complex K
  real(r8),allocatable :: FSBSTHeter(:,:)                !limitation of primary substrate for heterotrophs, [0->1, less limitation]              
  real(r8),allocatable :: FSBSTAutor(:)                  !limitation of primary substrate for autotrophs, [0->1, less limitation]              
  real(r8),allocatable :: ROQC4HeterMicActCmpK(:)        !microbial activity to trigger hydrolysis of organic complex
  real(r8),allocatable :: RHydrolysisScalCmpK(:)         !scalar for solid-organic matter hydrolysis, [0->1, faster]
  contains
    procedure, public :: Init => mic_diag_init
    procedure, public :: ZeroOut => mic_diag_zero
    procedure, public :: Destroy => mic_diag_destroy
    procedure, public :: Summary => mic_diag_summary
  end type Microbe_Diag_type
  contains
  
!------------------------------------------------------------------------------------------

  subroutine nit_aqmf_diag(this)
  implicit none
  class(Cumlate_Flux_Diag_type) :: this

  this%TFNH4B                 = 0._r8
  this%TFNO3B                 = 0._r8
  this%TFNO2B                 = 0._r8
  this%TFP14B                 = 0._r8
  this%TFPO4B                 = 0._r8
  this%tCH4ProdAceto          = 0._r8
  this%tCH4ProdH2             = 0._r8
  this%TFOQC                  = 0._r8
  this%TFOQA                  = 0._r8
  this%TFOXYX                 = 0._r8
  this%TFNH4X                 = 0._r8
  this%TFNO3X                 = 0._r8
  this%TFNO2X                 = 0._r8
  this%TFN2OX                 = 0._r8
  this%TFP14X                 = 0._r8
  this%TFPO4X                 = 0._r8
  this%tRespGrossHeterUlm     = 0._r8
  this%tCResp4H2Prod          = 0._r8
  this%tRNH4MicrbImobilSoil   = 0._r8
  this%tRH2PO4MicrbImobilSoil = 0._r8
  this%tRNO3MicrbImobilSoil   = 0.0_r8
  this%tRH1PO4MicrbImobilSoil = 0.0_r8
  this%tRNH4MicrbImobilBand   = 0.0_r8
  this%tRNO3MicrbImobilBand   = 0.0_r8
  this%tRH2PO4MicrbImobilBand = 0.0_r8
  this%tRH1PO4MicrbImobilBand = 0.0_r8
  this%tRCO2MicrbProd         = 0.0_r8
  this%tRCH4MicrbProd         = 0.0_r8
  this%tRNOxMicrbRedux        = 0.0_r8
  this%tRCO2GrothAutor        = 0.0_r8
  this%TProdH2                = 0.0_r8
  this%tRO2MicrbUptk          = 0.0_r8
  this%TReduxNO3Soil          = 0.0_r8
  this%TReduxNO3Band          = 0.0_r8
  this%TReduxNO2toN2OSoil          = 0.0_r8
  this%TReduxNO2toN2OBand     = 0.0_r8
  this%TReduxN2OtoN2              = 0.0_r8
  this%TFixN2                 = 0.0_r8
  this%tCH4OxiAero            = 0.0_r8
  this%tCH4OxiANMO            = 0._r8
  this%tRO2UptkHeterG         = 0._r8
  this%tNO2OxiAuto            = 0._r8
  this%tRO2DmndHeterG         = 0._r8
  this%TDeniReduxNO2Band      = 0._r8
  this%TDeniReduxNO2Soil      = 0._r8
  this%TNitNO2Redux2N2OBand   = 0._r8
  this%TNitNO2Redux2N2OSoil   = 0._r8
  this%tRNH3Oxi               = 0._r8
  this%RNO2ReduxSoilChemo     = 0._r8
  this%RNO2ReduxBandChemo     = 0._r8
  this%RN2OProdSoilChemo      = 0._r8
  this%RN2OProdBandChemo      = 0._r8
  this%RNO3ProdSoilChemo      = 0._r8
  this%RNO3ProdBandChemo      = 0._r8
  this%RNO2ReduxChemo         = 0._r8

  end subroutine nit_aqmf_diag

!------------------------------------------------------------------------------------------

  subroutine nit_micf_init(this,jcplx)
  implicit none
  class(Microbe_Flux_type) :: this
  integer, intent(in) :: jcplx
  integer :: ndbiomcp
  integer :: NumMicrobAutoTrophCmplx
  integer :: NumHetetr1MicCmplx
  ndbiomcp=micpar%ndbiomcp
  NumMicrobAutoTrophCmplx=micpar%NumMicrobAutoTrophCmplx
  NumHetetr1MicCmplx=micpar%NumHetetr1MicCmplx

  allocate(this%RO2UptkHeter(NumHetetr1MicCmplx,1:jcplx));this%RO2UptkHeter=spval
  allocate(this%RO2UptkHeterG(NumHetetr1MicCmplx));this%RO2UptkHeterG=spval

  allocate(this%Resp4NFixHeter(NumHetetr1MicCmplx,1:jcplx));this%Resp4NFixHeter=spval
  allocate(this%RespGrossHeter(NumHetetr1MicCmplx,1:jcplx));this%RespGrossHeter=spval
  allocate(this%RGrowthRespHeter(NumHetetr1MicCmplx,1:jcplx));this%RGrowthRespHeter=spval
  allocate(this%RMaintRespHeter(NumHetetr1MicCmplx,1:jcplx));this%RMaintRespHeter=spval
  allocate(this%RO2Dmnd4RespHeter(NumHetetr1MicCmplx,1:jcplx));this%RO2Dmnd4RespHeter=spval
  allocate(this%RO2DmndHeter(NumHetetr1MicCmplx,1:jcplx));this%RO2DmndHeter=spval
  allocate(this%RO2DmndHeterG(NumHetetr1MicCmplx));this%RO2DmndHeterG=spval
  allocate(this%RO2Uptk4RespHeter(NumHetetr1MicCmplx,1:jcplx));this%RO2Uptk4RespHeter=spval
  allocate(this%RNO3ReduxHeterSoil(NumHetetr1MicCmplx,1:jcplx));this%RNO3ReduxHeterSoil=spval
  allocate(this%RNO3ReduxHeterBand(NumHetetr1MicCmplx,1:jcplx));this%RNO3ReduxHeterBand=spval
  allocate(this%RNO2ReduxHeterSoil(NumHetetr1MicCmplx,1:jcplx));this%RNO2ReduxHeterSoil=spval
  allocate(this%RNO2ReduxHeterBand(NumHetetr1MicCmplx,1:jcplx));this%RNO2ReduxHeterBand=spval
  allocate(this%RN2OReduxHeter(NumHetetr1MicCmplx,1:jcplx));this%RN2OReduxHeter=spval
  allocate(this%RNOxDOCReduxRespDenitLim(NumHetetr1MicCmplx,1:jcplx));this%RNOxDOCReduxRespDenitLim=spval
  allocate(this%RNOxAcetReduxRespDenitLim(NumHetetr1MicCmplx,1:jcplx));this%RNOxAcetReduxRespDenitLim=spval
  allocate(this%RMaintDmndHeter(2,NumHetetr1MicCmplx,1:jcplx));this%RMaintDmndHeter=spval
  allocate(this%RNH4imobilSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RNH4imobilSoilHeter=spval
  allocate(this%RNO3imobilSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO3imobilSoilHeter=spval
  allocate(this%RH2PO4imobilSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RH2PO4imobilSoilHeter=spval
  allocate(this%RNH4imobilBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RNH4imobilBandHeter=spval
  allocate(this%RNO3imobilBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO3imobilBandHeter=spval
  allocate(this%RH2PO4imobilBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RH2PO4imobilBandHeter=spval
  allocate(this%RkillLitfalOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RkillLitfalOMHeter=spval
  allocate(this%RkillLitrfal2HumOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RkillLitrfal2HumOMHeter=spval
  allocate(this%RkillLitrfal2ResduOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RkillLitrfal2ResduOMHeter=spval
  allocate(this%DOMUptk4GrothHeter(NumPlantChemElms,NumHetetr1MicCmplx,1:jcplx));this%DOMUptk4GrothHeter=spval
  allocate(this%RH2ProdHeter(NumHetetr1MicCmplx,1:jcplx));this%RH2ProdHeter=spval
  allocate(this%RMaintDefcitLitrfalOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RMaintDefcitLitrfalOMHeter=spval
  allocate(this%RMaintDefLitrfal2HumOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RMaintDefLitrfal2HumOMHeter=spval
  allocate(this%RMaintDefLitrfal2ResduOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RMaintDefLitrfal2ResduOMHeter=spval
  allocate(this%RN2FixHeter(NumHetetr1MicCmplx,1:jcplx));this%RN2FixHeter=spval
  allocate(this%RKillOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RKillOMHeter=spval
  allocate(this%RkillRecycOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RkillRecycOMHeter=spval
  allocate(this%RMaintDefcitKillOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RMaintDefcitKillOMHeter=spval
  allocate(this%RMaintDefcitRecycOMHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%RMaintDefcitRecycOMHeter=spval
  allocate(this%RMetabDOCUptkHeter(NumHetetr1MicCmplx,1:jcplx));this%RMetabDOCUptkHeter=spval
  allocate(this%RMetabAcetUptkHeter(NumHetetr1MicCmplx,1:jcplx));this%RMetabAcetUptkHeter=spval
  allocate(this%RNH4imobilLitrHeter(NumHetetr1MicCmplx,1:jcplx));this%RNH4imobilLitrHeter=spval
  allocate(this%RNO3imobilLitrHeter(NumHetetr1MicCmplx,1:jcplx));this%RNO3imobilLitrHeter=spval
  allocate(this%RH2PO4imobilLitrHeter(NumHetetr1MicCmplx,1:jcplx));this%RH2PO4imobilLitrHeter=spval
  allocate(this%RNOxReduxRespDenitUlm(NumHetetr1MicCmplx,1:jcplx));this%RNOxReduxRespDenitUlm=spval
  allocate(this%NonstX2stBiomHeter(NumPlantChemElms,2,NumHetetr1MicCmplx,1:jcplx));this%NonstX2stBiomHeter=spval
  allocate(this%RCO2ProdHeter(NumHetetr1MicCmplx,1:jcplx));this%RCO2ProdHeter=spval
  allocate(this%RCO2FixCyano(micpar%JGniH(micpar%mid_HeterMixtCynoBacter):micpar%JGnfH(micpar%mid_HeterMixtCynoBacter),1:jcplx));this%RCO2FixCyano=spval
  allocate(this%RAcetateProdHeter(NumHetetr1MicCmplx,1:jcplx));this%RAcetateProdHeter=spval
  allocate(this%RCH4ProdHeter(NumHetetr1MicCmplx,1:jcplx));this%RCH4ProdHeter=spval
  allocate(this%RH1PO4imobilSoilHeter(NumHetetr1MicCmplx,1:jcplx));this%RH1PO4imobilSoilHeter=spval
  allocate(this%RH1PO4imobilBandHeter(NumHetetr1MicCmplx,1:jcplx));this%RH1PO4imobilBandHeter=spval
  allocate(this%RH1PO4imobilLitrHeter(NumHetetr1MicCmplx,1:jcplx));this%RH1PO4imobilLitrHeter=spval

  allocate(this%RSMetaOxidSoilAutor(NumMicrobAutoTrophCmplx));this%RSMetaOxidSoilAutor=spval
  allocate(this%RSMetaOxidBandAutor(NumMicrobAutoTrophCmplx));this%RSMetaOxidBandAutor=spval
  allocate(this%XferBiomeHeterK(1:NumPlantChemElms,3,NumHetetr1MicCmplx,1:jcplx));this%XferBiomeHeterK=spval
  allocate(this%ROQC4HeterMicrobAct(NumHetetr1MicCmplx,1:jcplx));this%ROQC4HeterMicrobAct=spval
  allocate(this%RCCMEheter(NumPlantChemElms,ndbiomcp,1:jcplx));this%RCCMEheter=spval
  allocate(this%RCCMEAutor(NumPlantChemElms,ndbiomcp,1:jcplx));this%RCCMEAutor=spval
  allocate(this%ECHZAutor(1:NumMicrobAutoTrophCmplx));this%ECHZAutor=spval
  allocate(this%ECHZHeter(1:NumHetetr1MicCmplx,1:jcplx));this%ECHZHeter=spval
  allocate(this%FOQC(1:NumHetetr1MicCmplx,1:jcplx));this%FOQC=spval
  allocate(this%FOQA(1:NumHetetr1MicCmplx,1:jcplx));this%FOQA=spval  
  allocate(this%RGOCP(1:NumHetetr1MicCmplx,1:jcplx));this%RGOCP=spval
  allocate(this%RGOAP(1:NumHetetr1MicCmplx,1:jcplx));this%RGOAP=spval
  allocate(this%FGOCP(1:NumHetetr1MicCmplx,1:jcplx));this%FGOCP=spval
  allocate(this%FGOAP(1:NumHetetr1MicCmplx,1:jcplx));this%FGOAP=spval  
  allocate(this%fPhotoR(micpar%JGniH(micpar%mid_HeterMixtCynoBacter):micpar%JGnfH(micpar%mid_HeterMixtCynoBacter),1:jcplx)); this%fPhotoR=spval
  allocate(this%RO2UptkAutor(NumMicrobAutoTrophCmplx));this%RO2UptkAutor=spval
  allocate(this%Resp4NFixAutor(NumMicrobAutoTrophCmplx));this%Resp4NFixAutor=spval
  allocate(this%RespGrossAutor(NumMicrobAutoTrophCmplx));this%RespGrossAutor=spval
  allocate(this%RO2Dmnd4GrossRespAutor(NumMicrobAutoTrophCmplx));this%RO2Dmnd4GrossRespAutor=spval
  allocate(this%RO2Uptk4RespAutor(NumMicrobAutoTrophCmplx));this%RO2Uptk4RespAutor=spval
  allocate(this%RNO3UptkAutor(NumMicrobAutoTrophCmplx));this%RNO3UptkAutor=spval
  allocate(this%RNOxReduxAutorSoil(NumMicrobAutoTrophCmplx));this%RNOxReduxAutorSoil=spval
  allocate(this%RNOxReduxAutorBand(NumMicrobAutoTrophCmplx));this%RNOxReduxAutorBand=spval
  allocate(this%RNOxReduxRespAutorLim(NumMicrobAutoTrophCmplx));this%RNOxReduxRespAutorLim=spval
  allocate(this%RMaintDmndAutor(2,NumMicrobAutoTrophCmplx));this%RMaintDmndAutor=spval
  allocate(this%RNH4TransfSoilAutor(NumMicrobAutoTrophCmplx));this%RNH4TransfSoilAutor=spval
  allocate(this%RNO3TransfSoilAutor(NumMicrobAutoTrophCmplx));this%RNO3TransfSoilAutor=spval
  allocate(this%RH2PO4TransfSoilAutor(NumMicrobAutoTrophCmplx));this%RH2PO4TransfSoilAutor=spval
  allocate(this%RNH4TransfBandAutor(NumMicrobAutoTrophCmplx));this%RNH4TransfBandAutor=spval
  allocate(this%RNO3TransfBandAutor(NumMicrobAutoTrophCmplx));this%RNO3TransfBandAutor=spval
  allocate(this%RH2PO4TransfBandAutor(NumMicrobAutoTrophCmplx));this%RH2PO4TransfBandAutor=spval
  allocate(this%RkillLitfalOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RkillLitfalOMAutor=spval
  allocate(this%RkillLitrfal2HumOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RkillLitrfal2HumOMAutor=spval
  allocate(this%RkillLitrfal2ResduOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RkillLitrfal2ResduOMAutor=spval
  allocate(this%DOMuptk4GrothAutor(idom_beg:idom_end,NumMicrobAutoTrophCmplx));this%DOMuptk4GrothAutor=spval
  allocate(this%RMaintDefcitLitrfalOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RMaintDefcitLitrfalOMAutor=spval
  allocate(this%RMaintDefLitrfal2HumOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RMaintDefLitrfal2HumOMAutor=spval
  allocate(this%RMaintDefLitrfal2ResduOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RMaintDefLitrfal2ResduOMAutor=spval
  allocate(this%RN2FixAutor(NumMicrobAutoTrophCmplx));this%RN2FixAutor=spval
  allocate(this%RKillOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RKillOMAutor=spval
  allocate(this%RkillRecycOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RkillRecycOMAutor=spval
  allocate(this%RMaintDefcitKillOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RMaintDefcitKillOMAutor=spval
  allocate(this%RMaintDefcitRecycOMAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%RMaintDefcitRecycOMAutor=spval
  allocate(this%RNH4TransfLitrAutor(NumMicrobAutoTrophCmplx));this%RNH4TransfLitrAutor=spval
  allocate(this%RNO3TransfLitrAutor(NumMicrobAutoTrophCmplx));this%RNO3TransfLitrAutor=spval
  allocate(this%RH2PO4TransfLitrAutor(NumMicrobAutoTrophCmplx));this%RH2PO4TransfLitrAutor=spval

  allocate(this%NonstX2stBiomAutor(NumPlantChemElms,2,NumMicrobAutoTrophCmplx));this%NonstX2stBiomAutor=spval
  allocate(this%RCO2ProdAutor(NumMicrobAutoTrophCmplx));this%RCO2ProdAutor=spval
  allocate(this%RCO2XumpAutor(NumMicrobAutoTrophCmplx));this%RCO2XumpAutor=spval
  allocate(this%ROMProdCO2Autor(NumMicrobAutoTrophCmplx));this%ROMProdCO2Autor=spval
  allocate(this%RGrowthCAutor(NumMicrobAutoTrophCmplx));this%RGrowthCAutor=spval
  allocate(this%RCH4ProdAutor(NumMicrobAutoTrophCmplx));this%RCH4ProdAutor=spval
  allocate(this%RH1PO4TransfSoilAutor(NumMicrobAutoTrophCmplx));this%RH1PO4TransfSoilAutor=spval
  allocate(this%RH1PO4TransfBandAutor(NumMicrobAutoTrophCmplx));this%RH1PO4TransfBandAutor=spval
  allocate(this%RH1PO4TransfLitrAutor(NumMicrobAutoTrophCmplx));this%RH1PO4TransfLitrAutor=spval

  call this%ZeroOut()
  end subroutine nit_micf_init

!------------------------------------------------------------------------------------------

  subroutine nit_mics_init(this, jcplx)

  implicit none
  class(Microbe_State_type) :: this
  integer, intent(in) :: jcplx
  integer :: NumMicrobAutoTrophCmplx,NumHetetr1MicCmplx
  NumMicrobAutoTrophCmplx=micpar%NumMicrobAutoTrophCmplx
  NumHetetr1MicCmplx=micpar%NumHetetr1MicCmplx

  allocate(this%rCNBiomeActHeter(NumPlantChemElms,NumHetetr1MicCmplx,1:jcplx));this%rCNBiomeActHeter=1._r8
  allocate(this%OMActHeter(NumHetetr1MicCmplx,1:jcplx));this%OMActHeter=spval
  allocate(this%FracOMActHeter(NumHetetr1MicCmplx,1:jcplx));this%FracOMActHeter=spval
  allocate(this%FracNO2ReduxHeter(NumHetetr1MicCmplx,1:jcplx));this%FracNO2ReduxHeter=spval
  allocate(this%FracHeterBiomOfActK(NumHetetr1MicCmplx,1:jcplx));this%FracHeterBiomOfActK=spval
  allocate(this%OMC2(NumHetetr1MicCmplx,1:jcplx));this%OMC2=spval
  allocate(this%GrowthEnvScalHeter(NumHetetr1MicCmplx,1:jcplx));this%GrowthEnvScalHeter=spval
  allocate(this%TempMaintRHeter(NumHetetr1MicCmplx,1:jcplx));this%TempMaintRHeter=spval
  allocate(this%OMN2(NumHetetr1MicCmplx,1:jcplx));this%OMN2=spval
  allocate(this%FOM2(NumHetetr1MicCmplx,1:jcplx));this%FOM2=spval
  allocate(this%OxyLimterHeter(NumHetetr1MicCmplx,1:jcplx));this%OxyLimterHeter=spval
  allocate(this%FCN(NumHetetr1MicCmplx,1:jcplx));this%FCN=spval
  allocate(this%FCP(NumHetetr1MicCmplx,1:jcplx));this%FCP=spval
  allocate(this%FBiomStoiScalarHeter(NumHetetr1MicCmplx,1:jcplx));this%FBiomStoiScalarHeter=spval

  allocate(this%rCNBiomeActAutor(NumPlantChemElms,NumMicrobAutoTrophCmplx));this%rCNBiomeActAutor=1.0_r8
  allocate(this%OMActAutor(NumMicrobAutoTrophCmplx));this%OMActAutor=spval
  allocate(this%FracOMActAutor(NumMicrobAutoTrophCmplx));this%FracOMActAutor=spval
  allocate(this%FracNO2XupAutor(NumMicrobAutoTrophCmplx));this%FracNO2XupAutor=spval
  allocate(this%FracAutorBiomOfActK(NumMicrobAutoTrophCmplx));this%FracAutorBiomOfActK=spval
  allocate(this%OMC2Autor(NumMicrobAutoTrophCmplx));this%OMC2Autor=spval
  allocate(this%GrowthEnvScalAutor(NumMicrobAutoTrophCmplx));this%GrowthEnvScalAutor=spval
  allocate(this%TSensMaintRAutor(NumMicrobAutoTrophCmplx));this%TSensMaintRAutor=spval
  allocate(this%TSensGroAutor(NumMicrobAutoTrophCmplx));this%TSensGroAutor=spval
  allocate(this%WSensGroAutor(NumMicrobAutoTrophCmplx));this%WSensGroAutor=spval
  allocate(this%WSensGroHeter(NumHetetr1MicCmplx,1:jcplx));this%WSensGroHeter=spval
  allocate(this%TSensGroHeter(NumHetetr1MicCmplx,1:jcplx));this%TSensGroHeter=spval
  allocate(this%OMN2Autor(NumMicrobAutoTrophCmplx));this%OMN2Autor=spval
  allocate(this%FOM2Autor(NumMicrobAutoTrophCmplx));this%FOM2Autor=spval
  allocate(this%fLimO2Autor(NumMicrobAutoTrophCmplx));this%fLimO2Autor=spval
  allocate(this%FCNAutor(NumMicrobAutoTrophCmplx));this%FCNAutor=spval
  allocate(this%FCPAutor(NumMicrobAutoTrophCmplx));this%FCPAutor=spval
  allocate(this%FBiomNutStoiScalAutor(NumMicrobAutoTrophCmplx));this%FBiomNutStoiScalAutor=spval
  end subroutine nit_mics_init
!------------------------------------------------------------------------------------------

  subroutine nit_micf_zero(this)
  implicit none
  class(Microbe_Flux_type) :: this

  this%RO2UptkHeter                     = 0._r8
  this%RO2UptkHeterG                    = 0._r8
  this%Resp4NFixHeter                   = 0._r8
  this%RespGrossHeter                   = 0._r8
  this%RGrowthRespHeter                 = 0._r8
  this%RMaintRespHeter                  = 0._r8
  this%RO2Dmnd4RespHeter                = 0._r8
  this%RO2DmndHeter                     = 0._r8
  this%RO2DmndHeterG                    = 0._r8
  this%RO2Uptk4RespHeter                = 0._r8
  this%RNO3ReduxHeterSoil               = 0._r8
  this%RNO3ReduxHeterBand               = 0._r8
  this%RNO2ReduxHeterSoil               = 0._r8
  this%RNO2ReduxHeterBand               = 0._r8
  this%RN2OReduxHeter                   = 0._r8
  this%RNOxDOCReduxRespDenitLim         = 0._r8
  this%RNOxAcetReduxRespDenitLim        = 0._r8
  this%RMaintDmndHeter                  = 0._r8
  this%RNH4imobilSoilHeter              = 0._r8
  this%RNO3imobilSoilHeter              = 0._r8
  this%RH2PO4imobilSoilHeter            = 0._r8
  this%RNH4imobilBandHeter              = 0._r8
  this%RNO3imobilBandHeter              = 0._r8
  this%RH2PO4imobilBandHeter            = 0._r8
  this%RkillLitfalOMHeter               = 0._r8
  this%RkillLitrfal2HumOMHeter          = 0._r8
  this%RkillLitrfal2ResduOMHeter        = 0._r8
  this%DOMUptk4GrothHeter               = 0._r8
  this%RH2ProdHeter                     = 0._r8
  this%RMaintDefLitrfal2HumOMHeter   = 0._r8
  this%RMaintDefLitrfal2ResduOMHeter = 0._r8
  this%RMaintDefcitLitrfalOMHeter       = 0._r8
  this%RCCMEheter                       = 0._r8
  this%RCCMEAutor                       = 0._r8
  this%RN2FixHeter                      = 0._r8

  this%RKillOMHeter             = 0._r8
  this%RkillRecycOMHeter        = 0._r8
  this%RMaintDefcitKillOMHeter  = 0._r8
  this%RMaintDefcitRecycOMHeter = 0._r8
  this%RMetabDOCUptkHeter      = 0._r8
  this%RMetabAcetUptkHeter     = 0._r8
  this%RNH4imobilLitrHeter      = 0._r8
  this%RNO3imobilLitrHeter      = 0._r8
  this%RH2PO4imobilLitrHeter    = 0._r8
  this%RNOxReduxRespDenitUlm    = 0._r8
  this%ROQC4HeterMicrobAct      = 0._r8
  this%NonstX2stBiomHeter       = 0._r8
  this%RCO2ProdHeter            = 0._r8
  this%RCO2FixCyano             = 0._r8
  this%RAcetateProdHeter          = 0._r8
  this%RCH4ProdHeter            = 0._r8
  this%RSMetaOxidSoilAutor          = 0._r8
  this%RSMetaOxidBandAutor          = 0._r8
  this%XferBiomeHeterK          = 0._r8
  this%RH1PO4imobilSoilHeter    = 0._r8
  this%RH1PO4imobilBandHeter    = 0._r8
  this%RH1PO4imobilLitrHeter    = 0._r8
  this%ECHZAutor                = 0._r8
  this%ECHZHeter                = 0._r8
  this%FOQC                     = 0._r8
  this%FOQA                     = 0._r8  
  this%FGOCP                    = 0._r8
  this%RGOCP                    = 0._r8
  this%RGOAP                    = 0._r8
  this%FGOAP                    = 0._r8
  this%fPhotoR                  = 0._r8
  this%RO2UptkAutor                     = 0._r8
  this%Resp4NFixAutor                   = 0._r8
  this%RespGrossAutor                   = 0._r8
  this%RO2Dmnd4GrossRespAutor                = 0._r8
  this%RO2Uptk4RespAutor                = 0._r8
  this%RNO3UptkAutor                    = 0._r8
  this%RNOxReduxAutorSoil               = 0._r8
  this%RNOxReduxAutorBand               = 0._r8
  this%RNOxReduxRespAutorLim            = 0._r8
  this%RMaintDmndAutor                  = 0._r8
  this%RNH4TransfSoilAutor              = 0._r8
  this%RNO3TransfSoilAutor              = 0._r8
  this%RH2PO4TransfSoilAutor            = 0._r8
  this%RNH4TransfBandAutor              = 0._r8
  this%RNO3TransfBandAutor              = 0._r8
  this%RH2PO4TransfBandAutor            = 0._r8
  this%RkillLitfalOMAutor               = 0._r8
  this%RkillLitrfal2HumOMAutor          = 0._r8
  this%RkillLitrfal2ResduOMAutor        = 0._r8
  this%DOMuptk4GrothAutor               = 0._r8
  this%RMaintDefcitLitrfalOMAutor       = 0._r8
  this%RMaintDefLitrfal2HumOMAutor   = 0._r8
  this%RMaintDefLitrfal2ResduOMAutor = 0._r8
  this%RN2FixAutor                      = 0._r8

  this%RKillOMAutor             = 0._r8
  this%RkillRecycOMAutor        = 0._r8
  this%RMaintDefcitKillOMAutor  = 0._r8
  this%RMaintDefcitRecycOMAutor = 0._r8
  this%RNH4TransfLitrAutor      = 0._r8
  this%RNO3TransfLitrAutor      = 0._r8
  this%RH2PO4TransfLitrAutor    = 0._r8  
  this%NonstX2stBiomAutor       = 0._r8
  this%RCO2ProdAutor            = 0._r8
  this%RCO2XumpAutor            = 0._r8
  this%ROMProdCO2Autor          = 0._r8
  this%RGrowthCAutor            = 0._r8 
  this%RCH4ProdAutor            = 0._r8
  this%RH1PO4TransfSoilAutor    = 0._r8
  this%RH1PO4TransfBandAutor    = 0._r8
  this%RH1PO4TransfLitrAutor    = 0._r8

  end subroutine nit_micf_zero


!------------------------------------------------------------------------------------------

  subroutine nit_micf_destroy(this)

  implicit none
  class(Microbe_Flux_type) :: this

  call destroy(this%RO2UptkHeter)
  call destroy(this%Resp4NFixHeter)
  call destroy(this%RespGrossHeter)
  call destroy(this%RGrowthRespHeter)
  call destroy(this%RMaintRespHeter)
  call destroy(this%RO2Dmnd4RespHeter)
  call destroy(this%RO2DmndHeter)
  call destroy(this%RO2DmndHeterG)
  call destroy(this%RO2Uptk4RespHeter)
  call destroy(this%RNO3ReduxHeterSoil)
  call destroy(this%RNO3ReduxHeterBand)
  call destroy(this%RNO2ReduxHeterSoil)
  call destroy(this%RNO2ReduxHeterBand)
  call destroy(this%RN2OReduxHeter)
  call destroy(this%RNOxDOCReduxRespDenitLim)
  call destroy(this%RNOxAcetReduxRespDenitLim)
  call destroy(this%RMaintDmndHeter)
  call destroy(this%RNH4imobilSoilHeter)
  call destroy(this%RNO3imobilSoilHeter)
  call destroy(this%RH2PO4imobilSoilHeter)
  call destroy(this%RNH4imobilBandHeter)
  call destroy(this%RNO3imobilBandHeter)
  call destroy(this%RH2PO4imobilBandHeter)
  call destroy(this%RkillLitfalOMHeter)
  call destroy(this%RkillLitrfal2HumOMHeter)
  call destroy(this%RkillLitrfal2ResduOMHeter)
  call destroy(this%DOMUptk4GrothHeter)
  call destroy(this%RH2ProdHeter)
  call destroy(this%RMaintDefLitrfal2HumOMHeter)
  call destroy(this%RMaintDefLitrfal2ResduOMHeter)
  call destroy(this%RMaintDefcitLitrfalOMHeter)
  call destroy(this%RCCMEheter)
  call destroy(this%RN2FixHeter)
  call destroy(this%RKillOMHeter)
  call destroy(this%RkillRecycOMHeter)
  call destroy(this%RMaintDefcitKillOMHeter)
  call destroy(this%RMaintDefcitRecycOMHeter)
  call destroy(this%RMetabDOCUptkHeter)
  call destroy(this%RMetabAcetUptkHeter)
  call destroy(this%RNH4imobilLitrHeter)
  call destroy(this%RNO3imobilLitrHeter)
  call destroy(this%RH2PO4imobilLitrHeter)
  call destroy(this%RNOxReduxRespDenitUlm)
  call destroy(this%ROQC4HeterMicrobAct)
  call destroy(this%NonstX2stBiomHeter)
  call destroy(this%RCO2ProdHeter)
  call destroy(this%RCO2FixCyano)
  call destroy(this%RAcetateProdHeter)
  call destroy(this%RCH4ProdHeter)
  call destroy(this%RSMetaOxidSoilAutor)
  call destroy(this%RSMetaOxidBandAutor)
  call destroy(this%XferBiomeHeterK)
  call destroy(this%RH1PO4imobilSoilHeter)
  call destroy(this%RH1PO4imobilBandHeter)
  call destroy(this%RH1PO4imobilLitrHeter)
  call destroy(this%RO2UptkHeterG)
  call destroy(this%RO2UptkAutor)
  call destroy(this%Resp4NFixAutor)
  call destroy(this%RespGrossAutor)
  call destroy(this%RO2Dmnd4GrossRespAutor)
  call destroy(this%RO2Uptk4RespAutor)
  call destroy(this%RNO3UptkAutor)
  call destroy(this%RNOxReduxAutorSoil)
  call destroy(this%RNOxReduxAutorBand)
  call destroy(this%RNOxReduxRespAutorLim)
  call destroy(this%RMaintDmndAutor)
  call destroy(this%RNH4TransfSoilAutor)
  call destroy(this%RNO3TransfSoilAutor)
  call destroy(this%RH2PO4TransfSoilAutor)
  call destroy(this%RNH4TransfBandAutor)
  call destroy(this%RNO3TransfBandAutor)
  call destroy(this%RH2PO4TransfBandAutor)
  call destroy(this%RkillLitfalOMAutor)
  call destroy(this%RkillLitrfal2HumOMAutor)
  call destroy(this%RkillLitrfal2ResduOMAutor)
  call destroy(this%DOMuptk4GrothAutor)
  call destroy(this%RMaintDefcitLitrfalOMAutor)
  call destroy(this%RMaintDefLitrfal2HumOMAutor)
  call destroy(this%RMaintDefLitrfal2ResduOMAutor)
  call destroy(this%RN2FixAutor)
  call destroy(this%RKillOMAutor)
  call destroy(this%RkillRecycOMAutor)
  call destroy(this%RMaintDefcitKillOMAutor)
  call destroy(this%RMaintDefcitRecycOMAutor)
  call destroy(this%RNH4TransfLitrAutor)
  call destroy(this%RNO3TransfLitrAutor)
  call destroy(this%RH2PO4TransfLitrAutor)
  call destroy(this%ROMProdCO2Autor)
  call destroy(this%NonstX2stBiomAutor)
  call destroy(this%RCO2ProdAutor)
  call destroy(this%RCO2XumpAutor)
  call destroy(this%RGrowthCAutor)
  call destroy(this%RCH4ProdAutor)
  call destroy(this%RH1PO4TransfSoilAutor)
  call destroy(this%RH1PO4TransfBandAutor)
  call destroy(this%RH1PO4TransfLitrAutor)

  end subroutine nit_micf_destroy
!------------------------------------------------------------------------------------------

  subroutine mic_diag_summary(this,groupid,FSubs_limiter)
  implicit none
  class(Microbe_Diag_type) :: this
  integer, intent(in) :: groupid
  real(r8),intent(out):: FSubs_limiter   !mean substrate limitation
  integer :: NGL,K

  FSubs_limiter=0._r8

  if(groupid==micpar%mid_HeterAerobBacter)then
  
    DO K=1,micpar%jcplx
      DO NGL=micpar%JGniH(groupid),micpar%JGnfH(groupid)
        FSubs_limiter=FSubs_limiter+this%FSBSTHeter(NGL,K)
      ENDDO
    ENDDO    
    FSubs_limiter=FSubs_limiter/real((micpar%JGnfH(groupid)-micpar%JGniH(groupid)+1)*micpar%jcplx,kind=r8)

  elseif(groupid==micpar%mid_Aerob_Fungi)then
    DO K=1,micpar%jcplx
      DO NGL=micpar%JGniH(groupid),micpar%JGnfH(groupid)
        FSubs_limiter=FSubs_limiter+this%FSBSTHeter(NGL,K)
      ENDDO
    ENDDO
    FSubs_limiter=FSubs_limiter/real((micpar%JGnfH(groupid)-micpar%JGniH(groupid)+1)*micpar%jcplx,kind=r8)
  endif

  end subroutine mic_diag_summary

!------------------------------------------------------------------------------------------  
  subroutine mic_diag_destroy(this)
  implicit none
  class(Microbe_Diag_type) :: this

  call destroy(this%FSBSTHeter)
  call destroy(this%FSBSTAutor)
  call destroy(this%ROQC4HeterMicActCmpK)
  call destroy(this%RHydrolysisScalCmpK)
  call destroy(this%TOMEK)
  end subroutine mic_diag_destroy
!------------------------------------------------------------------------------------------  
  subroutine nit_mics_destroy(this)
  implicit none
  class(Microbe_State_type) :: this

  call destroy(this%rCNBiomeActHeter)
  call destroy(this%OMActHeter)
  call destroy(this%FracOMActHeter)
  call destroy(this%FracNO2ReduxHeter)
  call destroy(this%FracHeterBiomOfActK)
  call destroy(this%OMC2)
  call destroy(this%GrowthEnvScalHeter)
  call destroy(this%TempMaintRHeter)
  call destroy(this%OMN2)
  call destroy(this%FOM2)
  call destroy(this%OxyLimterHeter)
  call destroy(this%FCN)
  call destroy(this%FCP)
  call destroy(this%FBiomStoiScalarHeter)
  call destroy(this%rCNBiomeActAutor)
  call destroy(this%OMActAutor)
  call destroy(this%FracOMActAutor)
  call destroy(this%FracNO2XupAutor)
  call destroy(this%FracAutorBiomOfActK)
  call destroy(this%OMC2Autor)
  call destroy(this%GrowthEnvScalAutor)
  call destroy(this%TSensMaintRAutor)
  call destroy(this%TSensGroAutor)
  call destroy(this%WSensGroAutor)  
  call destroy(this%TSensGroHeter)
  call destroy(this%WSensGroHeter)
  call destroy(this%OMN2Autor)
  call destroy(this%FOM2Autor)
  call destroy(this%fLimO2Autor)
  call destroy(this%FCNAutor)
  call destroy(this%FCPAutor)
  call destroy(this%FBiomNutStoiScalAutor)
  end subroutine nit_mics_destroy
!------------------------------------------------------------------------------------------
  subroutine nit_omcplxf_init(this)
  implicit none
  class(OMCplx_Flux_type) :: this
  integer :: nkinets
  integer :: ncplx
  integer :: ndbiomcp

  nkinets=micpar%jsken
  ncplx=micpar%jcplx
  ndbiomcp=micpar%ndbiomcp
  allocate(this%RHydlysSolidOM(NumPlantChemElms,nkinets,1:ncplx))
  allocate(this%RHumifySolidOM(1:NumPlantChemElms,nkinets,1:ncplx))
  allocate(this%RDecmpProdDOM(1:NumPlantChemElms,nkinets,1:ncplx))
  allocate(this%RHydlysBioResduOM(1:NumPlantChemElms,ndbiomcp,1:ncplx))
  allocate(this%RHydlysSorptOM(idom_beg:idom_end,1:ncplx))
  allocate(this%RDOMSorp(idom_beg:idom_end,1:ncplx))
  allocate(this%TDOMUptkHeter(idom_beg:idom_end,1:ncplx))
  allocate(this%XferRespHeterK(1:ncplx))
  allocate(this%XferDOMK(idom_beg:idom_end,1:ncplx))

  call this%ZeroOut()
  end subroutine nit_omcplxf_init
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxf_destroy(this)
  implicit none
  class(OMCplx_Flux_type) :: this

  call destroy(this%RHydlysSolidOM)
  call destroy(this%RHumifySolidOM)
  call destroy(this%RDecmpProdDOM)
  call destroy(this%RHydlysBioResduOM)
  call destroy(this%RHydlysSorptOM)
  call destroy(this%RDOMSorp)
  call destroy(this%XferRespHeterK)
  call destroy(this%XferDOMK)

  end subroutine nit_omcplxf_destroy
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxf_zero(this)
  implicit none
  class(OMCplx_Flux_type) :: this

  this%RHydlysSolidOM       = 0._r8
  this%RHumifySolidOM       = 0._r8
  this%RDecmpProdDOM         = 0._r8
  this%RHydlysBioResduOM    = 0._r8
  this%RHydlysSorptOM       = 0._r8
  this%RDOMSorp             = 0._r8
  this%XferRespHeterK       = 0._r8
  this%XferDOMK             = 0._r8
  end subroutine nit_omcplxf_zero
!------------------------------------------------------------------------------------------
  subroutine mic_diag_init(this)
  implicit none
  class(Microbe_Diag_type) :: this
  integer :: NumMicrobAutoTrophCmplx,NumHetetr1MicCmplx
  integer :: jcplx

  jcplx=micpar%jcplx
  NumMicrobAutoTrophCmplx = micpar%NumMicrobAutoTrophCmplx
  NumHetetr1MicCmplx    = micpar%NumHetetr1MicCmplx

  allocate(this%FSBSTHeter(1:NumHetetr1MicCmplx,1:jcplx));   this%FSBSTHeter=0._r8
  allocate(this%FSBSTAutor(1:NumMicrobAutoTrophCmplx)); this%FSBSTAutor=0._r8
  allocate(this%ROQC4HeterMicActCmpK(1:jcplx));this%ROQC4HeterMicActCmpK=0._r8
  allocate(this%RHydrolysisScalCmpK(1:jcplx));this%RHydrolysisScalCmpK=0._r8
  allocate(this%TOMEK(1:NumPlantChemElms,1:jcplx));this%TOMEK=spval
  end subroutine mic_diag_init

!------------------------------------------------------------------------------------------
  subroutine nit_omcplxs_init(this)
  implicit none
  class(OMCplx_State_type) :: this
  integer :: ncplx

  ncplx=micpar%jcplx
  allocate(this%BulkSOMC(1:ncplx));this%BulkSOMC=spval
  allocate(this%TOMEAutoK(1:NumPlantChemElms)); this%TOMEAutoK=spval
  allocate(this%FOCA(1:ncplx));this%FOCA=spval
  allocate(this%FOAA(1:ncplx));this%FOAA=spval
  allocate(this%rCNDOM(1:ncplx));this%rCNDOM=spval
  allocate(this%rCPDOM(1:ncplx));this%rCPDOM=spval
  allocate(this%OMBioResduK(1:ncplx));this%OMBioResduK=spval
  allocate(this%SolidOMCK(1:ncplx));this%SolidOMCK=spval
  allocate(this%SolidOMActK(1:ncplx));this%SolidOMActK=spval
  allocate(this%tMaxNActMicrbK(1:ncplx));this%tMaxNActMicrbK=spval
  allocate(this%tMaxPActMicrbK(1:ncplx));this%tMaxPActMicrbK=spval
  allocate(this%CDOM(idom_beg:idom_end,1:ncplx));this%CDOM=spval
  call this%ZeroOut()
  end subroutine nit_omcplxs_init
!------------------------------------------------------------------------------------------

  subroutine mic_diag_zero(this)
  implicit none
  class(Microbe_Diag_type) :: this

  this%ROQC4HeterMicActCmpK = 0._r8
  this%RHydrolysisScalCmpK =0._r8
  this%H1P4T               = 0._r8
  this%H2P4T               = 0._r8
  this%RH2UptkAutor        = 0._r8
  this%ThetaLitr           = 0._r8
  this%ThetaZ              = 0._r8
  this%TOMBioResdu         = 0._r8
  this%TotActMicrobiom     = 0._r8
  this%TotBiomNO2Consumers = 0._r8
  this%TSensGrowth         = 0._r8
  this%TSensMaintR         = 0._r8
  this%VOLWZ               = 0._r8
  this%WatStressMicb       = 0._r8
  this%XCO2                = 0._r8
  this%ZNH4T               = 0._r8
  this%ZNO3T               = 0._r8
  this%ZNO2T               = 0._r8
  this%TOMEK                = 0._r8
  end subroutine mic_diag_zero
!------------------------------------------------------------------------------------------

  subroutine nit_omcplxs_zero(this)

  implicit none
  class(OMCplx_State_type) :: this

  this%BulkSOMC=0._r8
  this%FOCA=0._r8
  this%FOAA=0._r8
  this%rCNDOM=0._r8
  this%rCPDOM=0._r8

  this%OMBioResduK=0._r8
  this%SolidOMCK=0._r8
  this%SolidOMActK=0._r8
  this%tMaxNActMicrbK=0._r8
  this%tMaxPActMicrbK=0._r8

  end subroutine nit_omcplxs_zero


!------------------------------------------------------------------------------------------
  subroutine nit_omcplxs_destroy(this)
  implicit none
  class(OMCplx_State_type) :: this

  call destroy(this%BulkSOMC)
  call destroy(this%FOCA)
  call destroy(this%FOAA)
  call destroy(this%rCNDOM)
  call destroy(this%rCPDOM)
  call destroy(this%OMBioResduK)
  call destroy(this%SolidOMCK)
  call destroy(this%SolidOMActK)
  call destroy(this%tMaxNActMicrbK)
  call destroy(this%tMaxPActMicrbK)
  call destroy(this%CDOM)

  end subroutine nit_omcplxs_destroy
end module MicrobeDiagTypes
