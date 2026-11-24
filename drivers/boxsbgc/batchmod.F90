module batchmod
!j6aG$1bZcd5xgo1eAroLF36nN^
  use abortutils     , only : endrun
  use minimathMod    , only : addone,pmod
  use data_kind_mod  , only : r8 => DAT_KIND_R8
  use ModelStatusType, only : model_status_type
  use EcoSiMParDataMod  , only : micpar
  use MicForcTypeMod    , only : micforctype
  use MicFLuxTypeMod    , only : micfluxtype
  use MicStateTraitTypeMod, only : micsttype
  use fileUtil
  use TracerIDMod
  use minimathmod    , only : AZMAX1,AZMIN1
  use EcoSIMSolverPar
  use MicIDMod
  use ChemIDMod
  use ForcTypeMod         , only : forc_type

implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  logical :: Litlayer


  public :: getvarllen, getvarlist,initmodel
  public :: BatchModelConfig
  public :: RunMicBGC
contains

  function getvarllen()result(nvars)
  use MicBGCMod, only : initNitro1Layer
  implicit none
  integer :: nvars
  integer :: nmicbguilds

  call micpar%Init()
  ! Load microbial defaults and optional input before SetPars copies SPOSC_in.
  ! Initialize once here so the box driver uses the same order as EcoSIM.
  call initNitro1Layer
  call micpar%SetPars()

  call Initboxbgc(nvars)

  end function getvarllen
! ----------------------------------------------------------------------
  subroutine initmodel(nvars, ystates0l, forc, err_status)
!
! set initial conditions for the boxsbgc
  implicit none
  integer, intent(in) :: nvars
  type(forc_type), intent(in) :: forc
  real(r8), intent(inout) :: ystates0l(nvars)
  type(model_status_type), intent(out) :: err_status

  call err_status%reset()


  associate(                                                   &
    nlbiomcp                => micpar%nlbiomcp,                &
    ndbiomcp                => micpar%ndbiomcp,                &
    jsken                   => micpar%jsken,                   &
    NumHetetr1MicCmplx      => micpar%NumHetetr1MicCmplx,      &
    NumMicrobAutoTrophCmplx => micpar%NumMicrobAutoTrophCmplx, &
    NumLiveHeterBioms       => micpar%NumLiveHeterBioms,       &
    NumLiveAutoBioms        => micpar%NumLiveAutoBioms,        &
    jcplx                   => micpar%jcplx                    &
  )

  ystates0l(cid_oqc_b:cid_oqc_e)=forc%DOM(idom_doc,1:jcplx)
  ystates0l(cid_oqn_b:cid_oqn_e)=forc%DOM(idom_don,1:jcplx)
  ystates0l(cid_oqp_b:cid_oqp_e)=forc%DOM(idom_dop,1:jcplx)
  ystates0l(cid_oqa_b:cid_oqa_e)=forc%DOM(idom_acetate,1:jcplx)
  ystates0l(cid_ohc_b:cid_ohc_e)=forc%SorbedOM(ielmc,1:jcplx)
  ystates0l(cid_ohn_b:cid_ohn_e)=forc%SorbedOM(ielmn,1:jcplx)
  ystates0l(cid_ohp_b:cid_ohp_e)=forc%SorbedOM(ielmp,1:jcplx)
  ystates0l(cid_oha_b:cid_oha_e)=forc%SorbedOM(idom_acetate,1:jcplx)
  ystates0l(cid_osc_b:cid_osc_e)=reshape(forc%SolidOM(ielmc,1:jsken,1:jcplx),(/jsken*jcplx/))
  ystates0l(cid_osa_b:cid_osa_e)=reshape(forc%SolidOMAct(1:jsken,1:jcplx),(/jsken*jcplx/))
  ystates0l(cid_osn_b:cid_osn_e)=reshape(forc%SolidOM(ielmn,1:jsken,1:jcplx),(/jsken*jcplx/))
  ystates0l(cid_osp_b:cid_osp_e)=reshape(forc%SolidOM(ielmp,1:jsken,1:jcplx),(/jsken*jcplx/))
  ystates0l(cid_orc_b:cid_orc_e)=reshape(forc%OMBioResdu(ielmc,1:ndbiomcp,1:jcplx),(/ndbiomcp*jcplx/))
  ystates0l(cid_orn_b:cid_orn_e)=reshape(forc%OMBioResdu(ielmn,1:ndbiomcp,1:jcplx),(/ndbiomcp*jcplx/))
  ystates0l(cid_orp_b:cid_orp_e)=reshape(forc%OMBioResdu(ielmp,1:ndbiomcp,1:jcplx),(/ndbiomcp*jcplx/))

  ystates0l(cid_mBiomeHeter_b:cid_mBiomeHeter_e)=reshape(forc%mBiomeHeter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx),&
    (/NumPlantChemElms*NumLiveHeterBioms*jcplx/))
  ystates0l(cid_mBiomeAutor_b:cid_mBiomeAutor_e)=reshape(forc%mBiomeAutor(1:NumPlantChemElms,1:NumLiveAutoBioms),&
    (/NumPlantChemElms*NumLiveAutoBioms/))

  end associate
  end subroutine initmodel

! ----------------------------------------------------------------------

  subroutine BatchModelConfig(nvars,ystates0l,forc,micfor,micstt,micflx,err_status)
!
! DESCRIPTION:
! configure the batch mode of the soil bgc
  use MicStateTraitTypeMod, only : micsttype
  use MicForcTypeMod      , only : micforctype
  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(in) :: ystates0l(nvars)
  type(forc_type), intent(in) :: forc
  type(micforctype), intent(inout) :: micfor
  type(micsttype), intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(model_status_type), intent(out) :: err_status

  real(r8), parameter :: ZERO=1.0E-15_r8
  real(r8), parameter :: ZEROS2=1.0E-08_r8
  integer :: kk,jcplx

  call err_status%reset()

  associate(                                                   &
    nlbiomcp                => micpar%nlbiomcp,                &
    ndbiomcp                => micpar%ndbiomcp,                &
    NumHetetr1MicCmplx      => micpar%NumHetetr1MicCmplx,      &
    NumMicrobAutoTrophCmplx => micpar%NumMicrobAutoTrophCmplx, &
    k_humus                 => micpar%k_humus,                 &
    k_POM                   => micpar%k_POM,                   &
    icarbhyro               => micpar%icarbhyro,               &
    iprotein                => micpar%iprotein,                &
    jsken                   => micpar%jsken,                   &
    NumLiveHeterBioms       => micpar%NumLiveHeterBioms,       &
    NumLiveAutoBioms        => micpar%NumLiveAutoBioms,        &
    jcplx                   => micpar%jcplx                    &
  )
  micfor%ZERO  =ZERO
  micfor%ZEROS2=ZEROS2
  micfor%ZEROS =ZERO

  micfor%CCH4E =forc%CCH4E
  micfor%COXYE =forc%COXYE
  micfor%O2_irrig_conc  =0._r8      !oxygen concentration in surface irrigation
  micfor%O2_rain_conc  =0._r8      !oxygen concentration in precipitation
  micfor%Irrig2LitRSurf_col =0._r8      !irrigation flux into surface litter
  micfor%Rain2LitRSurf =0._r8      !precipitation flux into surface litter
  micfor%TempOffset=forc%TempOffset
  micfor%VLitR  =forc%VLitR
  micfor%VWatLitRHoldCapcity=forc%VWatLitRHoldCapcity
  micfor%VLWatMicP  =forc%VLWatMicP
  micfor%VLSoilMicP  =forc%VLSoilMicP
  micfor%THETY =forc%THETY
  micfor%POROS =forc%POROS
  micfor%FieldCapacity    =forc%FieldCapacity
  micfor%TKS   =forc%TKS
  micfor%THETW =forc%THETW
  micfor%PH    =forc%PH
  micfor%SoilMicPMassLayer  =forc%SoilMicPMassLayer
  micfor%VLSoilPoreMicP  =forc%VLSoilPoreMicP
  micfor%TScal4Difsvity  =forc%TScal4Difsvity
  micfor%VLNOB =forc%VLNOB
  micfor%VLNO3 =forc%VLNO3
  micfor%VLNH4 =forc%VLNH4
  micfor%VLNHB =forc%VLNHB
  micfor%VLPO4 =forc%VLPO4
  micfor%VLPOB =forc%VLPOB
  micfor%PSISoilMatricP =forc%PSISoilMatricP
  micfor%O2AquaDiffusvity =forc%O2AquaDiffusvity
  micfor%ORGC                         = forc%ORGC
  micfor%RNO2EcoUptkSoilPrev          = ystates0l(fid_RNO2EcoUptkSoilPrev)
  micfor%RN2OEcoUptkSoilPrev          = ystates0l(fid_RN2OEcoUptkSoilPrev)
  micfor%RNO2EcoUptkBandPrev          = ystates0l(fid_RNO2EcoUptkBandPrev)
  micfor%RO2EcoDmndPrev               = ystates0l(fid_RO2EcoDmndPrev)
  micfor%RO2GasXchangePrev            = ystates0l(fid_RO2GasXchangePrev)
  micfor%RNH4EcoDmndBandPrev          = ystates0l(fid_RNH4EcoDmndBandPrev)
  micfor%RNO3EcoDmndBandPrev          = ystates0l(fid_RNO3EcoDmndBandPrev)
  micfor%RH2PO4EcoDmndBandPrev        = ystates0l(fid_RH2PO4EcoDmndBandPrev)
  micfor%RH1PO4EcoDmndBandPrev        = ystates0l(fid_RH1PO4EcoDmndBandPrev)
  micfor%RDOMEcoDmndPrev(1:jcplx)     = ystates0l(fid_RDOMEcoDmndPrev_b:fid_RDOMEcoDmndPrev_e)
  micfor%RAcetateEcoDmndPrev(1:jcplx) = ystates0l(fid_RAcetateEcoDmndPrev_b:fid_RAcetateEcoDmndPrev_e)
  micfor%RCH4PhysexchPrev             = 0._r8
  micfor%RO2AquaXchangePrev           = 0._r8
  micfor%ElmAllocmatMicrblitr2POM     = forc%ElmAllocmatMicrblitr2POM(1:ndbiomcp)
  micfor%litrm                        = .false.
  micfor%Lsurf=.True.
  if(micfor%litrm)then
!   the following hasn't been turned on yet
!
    micstt%ZNH4TU=AZMAX1(forc%ZNH4S)
    micstt%ZNO3TU=AZMAX1(forc%ZNO3S)
    micstt%H1P4TU=AZMAX1(forc%H1PO4)
    micstt%H2P4TU=AZMAX1(forc%H2PO4)
    micstt%CNH4BU=forc%CNH4B
    micstt%CNH4SU=forc%CNH4S
    micstt%CH2P4U=forc%CH2P4
    micstt%CH2P4BU=forc%CH2P4B
    micstt%CH1P4U=forc%CH1P4
    micstt%CH1P4BU=forc%CH1P4B
    micstt%CNO3SU=forc%CNO3S
    micstt%CNO3BU=forc%CNO3B
    micstt%SOMPomProtein(ielmc)=forc%SolidOM(ielmc,iprotein,k_POM)
    micstt%SOMPomProtein(ielmn)=forc%SolidOM(ielmn,iprotein,k_POM)
    micstt%SOMPomProtein(ielmp)=forc%SolidOM(ielmp,iprotein,k_POM)
    micstt%SOMHumProtein(ielmc)=forc%SolidOM(ielmc,iprotein,k_humus)
    micstt%SOMHumProtein(ielmn)=forc%SolidOM(ielmn,iprotein,k_humus)
    micstt%SOMHumProtein(ielmp)=forc%SolidOM(ielmp,iprotein,k_humus)
    micstt%SOMHumCarbohyd(ielmc)=forc%SolidOM(ielmc,icarbhyro,k_humus)
    micstt%SOMHumCarbohyd(ielmn)=forc%SolidOM(ielmn,icarbhyro,k_humus)
    micstt%SOMHumCarbohyd(ielmp)=forc%SolidOM(ielmp,icarbhyro,k_humus)
    micfor%RNH4EcoDmndLitrPrev=ystates0l(fid_RNH4EcoDmndSoilPrev)
    micfor%RNO3EcoDmndLitrPrev=ystates0l(fid_RNO3EcoDmndSoilPrev)
    micfor%RH2PO4EcoDmndLitrPrev=ystates0l(fid_RH2PO4EcoDmndSoilPrev)
    micfor%RH1PO4EcoDmndLitrPrev=ystates0l(fid_RH1PO4EcoDmndSoilPrev)
    micfor%VOLWU =forc%VLWatMicP
    micfor%ElmAllocmatMicrblitr2POMU=forc%ElmAllocmatMicrblitr2POM(1:ndbiomcp)
  else
    micfor%AEC=forc%AEC
    micstt%OXYG=ystates0l(cid_COXYG)*forc%VLsoiAirP
  endif
  micstt%CNH4B =forc%CNH4B
  micstt%CNO3B =forc%CNO3B
  micstt%CH2P4B=forc%CH2P4B

  micstt%CNH4S =ystates0l(cid_ZNH4S)/(forc%VLWatMicP*forc%VLNH4)
  micstt%CNO3S =ystates0l(cid_ZNO3S)/(forc%VLWatMicP*forc%VLNO3)
  micstt%CH2P4 =ystates0l(cid_H2PO4)/(forc%VLWatMicP*forc%VLPO4)
  micstt%CH1P4 =ystates0l(cid_H1PO4)/(forc%VLWatMicP*forc%VLPO4)
  micstt%CH1P4B=forc%CH1P4B
  micfor%RNH4EcoDmndSoilPrev =ystates0l(fid_RNH4EcoDmndSoilPrev)
  micfor%RNO3EcoDmndSoilPrev =ystates0l(fid_RNO3EcoDmndSoilPrev)
  micfor%RH2PO4EcoDmndSoilPrev =ystates0l(fid_RH2PO4EcoDmndSoilPrev)
  micfor%RH1PO4EcoDmndSoilPrev =ystates0l(fid_RH1PO4EcoDmndSoilPrev)
  micfor%VLWatMicP  =forc%VLWatMicP
  micfor%VLsoiAirP  =forc%VLsoiAirP
  if(micfor%Lsurf)then
    micfor%SoilMicPMassLayer0=forc%SoilMicPMassLayer
  endif
  micfor%DiffusivitySolutEff(1:NPH)  =forc%DiffusivitySolutEff
  micfor%FILM(1:NPH)  =forc%FILM
  micfor%THETPM(1:NPH)=forc%THETPM
  micfor%VLWatMicPM(1:NPH) =forc%VLWatMicP
  micfor%TortMicPM(1:NPH)  =forc%TortMicPM
  micfor%VLsoiAirPM(1:NPH) =forc%VLsoiAirP

  micstt%EPOC=forc%EPOC
  micstt%EHUM=forc%EHUM
  micstt%ZNH4B=forc%ZNH4B
  micstt%ZNH4S=ystates0l(cid_ZNH4S)
  micstt%ZNO3B=forc%ZNO3B
  micstt%ZNO3S=ystates0l(cid_ZNO3S)
  micstt%H1POB=forc%H1POB
  micstt%H1PO4=ystates0l(cid_H1PO4)
  micstt%ZNO2B=forc%ZNO2B
  micstt%ZNO2S=ystates0l(cid_ZNO2S)
  micstt%H2POB=forc%H2POB
  micstt%H2PO4=ystates0l(cid_H2PO4)
  micstt%CCO2S=ystates0l(cid_CO2S)/forc%VLWatMicP
  micstt%CNO2S=ystates0l(cid_ZNO2S)/(forc%VLWatMicP*forc%VLNO3)
  micstt%CNO2B=forc%CNO2B
  micstt%CZ2OS=ystates0l(cid_Z2OS)/forc%VLWatMicP
  micstt%Z2OS=ystates0l(cid_Z2OS)
  micstt%COXYS=ystates0l(cid_OXYS)/forc%VLWatMicP
  micstt%OXYS =ystates0l(cid_OXYS)
  micstt%O2GSolubility=forc%O2GSolubility
  micstt%COXYG=ystates0l(cid_COXYG)
  micstt%CZ2GS=ystates0l(cid_CZ2GS)
  micstt%CH2GS=ystates0l(cid_CH2GS)
  micstt%H2GS =ystates0l(cid_H2GS)
  micstt%CCH4G=ystates0l(cid_CCH4G)
  micstt%CH4S =ystates0l(cid_CH4S)
  micstt%CH4AquaSolubility=forc%CH4AquaSolubility
  micstt%ZNFN0=forc%ZNFN0
  micstt%ZNFNI=forc%ZNFNI

  micstt%DOM(idom_doc,1:jcplx)          = ystates0l(cid_oqc_b:cid_oqc_e)
  micstt%DOM(idom_don,1:jcplx)          = ystates0l(cid_oqn_b:cid_oqn_e)
  micstt%DOM(idom_dop,1:jcplx)          = ystates0l(cid_oqp_b:cid_oqp_e)
  micstt%DOM(idom_acetate,1:jcplx)      = ystates0l(cid_oqa_b:cid_oqa_e)
  micstt%SorbedOM(ielmc,1:jcplx)        = ystates0l(cid_ohc_b:cid_ohc_e)
  micstt%SorbedOM(ielmn,1:jcplx)        = ystates0l(cid_ohn_b:cid_ohn_e)
  micstt%SorbedOM(ielmp,1:jcplx)        = ystates0l(cid_ohp_b:cid_ohp_e)
  micstt%SorbedOM(idom_acetate,1:jcplx) = ystates0l(cid_oha_b:cid_oha_e)

  micstt%SolidOM(ielmc,1:jsken,1:jcplx)       = reshape(ystates0l(cid_osc_b:cid_osc_e),(/jsken,jcplx/))
  micstt%SolidOM(ielmn,1:jsken,1:jcplx)       = reshape(ystates0l(cid_osn_b:cid_osn_e),(/jsken,jcplx/))
  micstt%SolidOM(ielmp,1:jsken,1:jcplx)       = reshape(ystates0l(cid_osp_b:cid_osp_e),(/jsken,jcplx/))
  micstt%SolidOMAct(1:jsken,1:jcplx)          = reshape(ystates0l(cid_osa_b:cid_osa_e),(/jsken,jcplx/))
  micstt%OMBioResdu(ielmc,1:ndbiomcp,1:jcplx) = reshape(ystates0l(cid_orc_b:cid_orc_e),(/ndbiomcp,jcplx/))
  micstt%OMBioResdu(ielmn,1:ndbiomcp,1:jcplx) = reshape(ystates0l(cid_orn_b:cid_orn_e),(/ndbiomcp,jcplx/))
  micstt%OMBioResdu(ielmp,1:ndbiomcp,1:jcplx) = reshape(ystates0l(cid_orp_b:cid_orp_e),(/ndbiomcp,jcplx/))
  micstt%CNOSC(1:jsken,1:jcplx)               = forc%CNOSC(1:jsken,1:jcplx)
  micstt%CPOSC(1:jsken,1:jcplx)               = forc%CPOSC(1:jsken,1:jcplx)

  micstt%mBiomeHeter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx)=&
    reshape(ystates0l(cid_mBiomeHeter_b:cid_mBiomeHeter_e),(/NumPlantChemElms,NumLiveHeterBioms,jcplx/))
  micstt%mBiomeAutor(1:NumPlantChemElms,1:NumLiveAutoBioms)=reshape(ystates0l(cid_mBiomeAutor_b:cid_mBiomeAutor_e),&
    (/NumPlantChemElms,NumLiveAutoBioms/))
  
  micflx%RNH4DmndSoilHeter(1:NumHetetr1MicCmplx,1:jcplx)=reshape(ystates0l(fid_RNH4DmndSoilHeter_b:fid_RNH4DmndSoilHeter_e) &
    ,(/NumHetetr1MicCmplx,JCPLX/))
  micflx%RNH4DmndBandHeter(1:NumHetetr1MicCmplx,1:jcplx)=reshape(ystates0l(fid_RNH4DmndBandHeter_b:fid_RNH4DmndBandHeter_e),(/NumHetetr1MicCmplx,JCPLX/))
  micflx%RNO3DmndSoilHeter(1:NumHetetr1MicCmplx,1:jcplx)=reshape(ystates0l(fid_RNO3DmndSoilHeter_b:fid_RNO3DmndSoilHeter_e),(/NumHetetr1MicCmplx,JCPLX/))
  micflx%RNO3DmndBandHeter(1:NumHetetr1MicCmplx,1:jcplx)=reshape(ystates0l(fid_RNO3DmndBandHeter_b:fid_RNO3DmndBandHeter_e),(/NumHetetr1MicCmplx,JCPLX/))
  micflx%RH2PO4DmndSoilHeter(1:NumHetetr1MicCmplx,1:jcplx)=reshape(ystates0l(fid_RH2PO4DmndSoilHeter_b:fid_RH2PO4DmndSoilHeter_e),(/NumHetetr1MicCmplx,JCPLX/))
  micflx%RH2PO4DmndBandHeter(1:NumHetetr1MicCmplx,1:jcplx)=reshape(ystates0l(fid_RH2PO4DmndBandHeter_b:fid_RH2PO4DmndBandHeter_e),(/NumHetetr1MicCmplx,JCPLX/))
  micflx%RH1PO4DmndSoilHeter(1:NumHetetr1MicCmplx,1:jcplx)=reshape(ystates0l(fid_RH1PO4DmndSoilHeter_b:fid_RH1PO4DmndSoilHeter_e),(/NumHetetr1MicCmplx,JCPLX/))
  micflx%RH1PO4DmndBandHeter(1:NumHetetr1MicCmplx,1:jcplx)=reshape(ystates0l(fid_RH1PO4DmndBandHeter_b:fid_RH1PO4DmndBandHeter_e),(/NumHetetr1MicCmplx,JCPLX/))
  micflx%RO2DmndHetert(1:NumHetetr1MicCmplx,1:jcplx)=reshape(ystates0l(fid_RO2DmndHetert_b:fid_RO2DmndHetert_e),&
    (/NumHetetr1MicCmplx,JCPLX/))
  end associate
  end subroutine BatchModelConfig

! ----------------------------------------------------------------------


  subroutine Initboxbgc(nvars)
!
! DESCRIPTION
! Initialize the id of relevant model variables &
! obtain total number of variables
  implicit none

  integer, intent(out) :: nvars
  integer :: itemp

  associate(                                                   &
    jcplx                   => micpar%jcplx,                   &
    jsken                   => micpar%jsken,                   &
    NumMicrobAutoTrophCmplx   => micpar%NumMicrobAutoTrophCmplx,   &
    NumHetetr1MicCmplx       => micpar%NumHetetr1MicCmplx,       &
    NumLiveHeterBioms       => micpar%NumLiveHeterBioms,       &
    NumLiveAutoBioms        => micpar%NumLiveAutoBioms,        &
    ndbiomcp                => micpar%ndbiomcp,                &
    nlbiomcp                => micpar%nlbiomcp                 &
  )
  itemp=0
  cid_ZMG                     = addone(itemp)
  cid_ZNA                     = addone(itemp)
  cid_ZKA                     = addone(itemp)
  cid_CO2S                    = addone(itemp)
  cid_H1PO4_2e_aque_mole_conc           = addone(itemp)
  cid_H1PO4_2e_band_mole_conc      = addone(itemp)
  cid_H2PO4_1e_aque_mole_conc = addone(itemp)
  cid_H2PO4_1e_band_mole_conc      = addone(itemp)
  cid_NH3_aque_mole_conc      = addone(itemp)
  cid_NH3_aqu_band_mole_conc       = addone(itemp)
  cid_NH4_1p_aque_mole_conc   = addone(itemp)
  cid_NH4_1p_band_mole_conc        = addone(itemp)
  cid_XNH4_mole_conc               = addone(itemp)
  cid_XNH4_band_mole_conc          = addone(itemp)
  cid_XHPO4_band_mole_conc         = addone(itemp)
  cid_XH2PO4_band_mole_conc        = addone(itemp)
  cid_XROH_band_mole_conc           = addone(itemp)
  cid_XHPO4_mole_conc               = addone(itemp)
  cid_XROH2_band_mole_conc          = addone(itemp)
  cid_XH2PO4_mole_conc              = addone(itemp)
  cid_XROH1_mole_conc               = addone(itemp)
  cid_XROH2_mole_conc               = addone(itemp)
  cid_Precp_AlPO4_mole_conc         = addone(itemp)
  cid_PrecpB_AlPO4_mole_conc        = addone(itemp)
  cid_Precp_CaHPO4_mole_conc        = addone(itemp)
  cid_PrecpB_CaHPO4_mole_conc       = addone(itemp)
  cid_Precp_Ca5P3O12O3H3_mole_conc  = addone(itemp)
  cid_PrecpB_Ca5P3O12O3H3_mole_conc = addone(itemp)
  cid_Precp_CaH4P2O8_mole_conc      = addone(itemp)
  cid_PrecpB_CaH4P2O8_mole_conc     = addone(itemp)
  cid_Precp_FePO4_mole_conc         = addone(itemp)
  cid_PrecpB_FePO4_mole_conc         = addone(itemp)

  fid_TRChem_NH4_soil                  = addone(itemp)
  fid_TRChem_NH4_band_soil             = addone(itemp)
  fid_TRChem_NH3_soil_vr               = addone(itemp)
  fid_TRChem_NH3_band_soil             = addone(itemp)
  fid_TRChem_H1PO4_soil                = addone(itemp)
  fid_TRChem_H2PO4_soil                = addone(itemp)
  fid_TRChem_H1PO4_band_soil           = addone(itemp)
  fid_TRChem_H2PO4_band_soil           = addone(itemp)
  fid_TRChem_NH4_sorbed_soil           = addone(itemp)
  fid_TRChem_NH4_sorbed_band_soil      = addone(itemp)
  fid_TRChem_ROH_sorbed_soil           = addone(itemp)
  fid_TRChem_ROH2_sorbed_soil          = addone(itemp)
  fid_TRChem_RHPO4_sorbed_soil         = addone(itemp)
  fid_TRChem_RH2PO4_sorbed_soil        = addone(itemp)
  fid_TRChem_ROH_sorbed_band_soil      = addone(itemp)
  fid_TRChem_ROH2_sorbed_band_soil     = addone(itemp)
  fid_TRChem_RHPO4_sorbed_band_soil    = addone(itemp)
  fid_TRChem_RH2PO4_sorbed_band_soil   = addone(itemp)
  fid_TRChem_AlPO4_precip_soil         = addone(itemp)
  fid_TRChem_FePO4_precip_soil         = addone(itemp)
  fid_TRChem_CaHPO4_precip_soil        = addone(itemp)
  fid_TRChem_apatite_precip_soil       = addone(itemp)
  fid_TRChem_CaH4P2O8_precip_soil      = addone(itemp)
  fid_TRChem_AlPO4_precip_band_soil    = addone(itemp)
  fid_TRChem_FePO4_precip_band_soil    = addone(itemp)
  fid_TRChem_CaHPO4_precip_band_soil   = addone(itemp)
  fid_TRChem_apatite_precip_band_soil  = addone(itemp)
  fid_TRChem_CaH4P2O8_precip_band_soil = addone(itemp)
  fid_TRChem_Al_3p_soil                = addone(itemp)

  cid_oqc_b=addone(itemp);cid_oqc_e=cid_oqc_b+jcplx-1;itemp=cid_oqc_e
  cid_oqn_b=addone(itemp);cid_oqn_e=cid_oqn_b+jcplx-1;itemp=cid_oqn_e
  cid_oqp_b=addone(itemp);cid_oqp_e=cid_oqp_b+jcplx-1;itemp=cid_oqp_e
  cid_oqa_b=addone(itemp);cid_oqa_e=cid_oqa_b+jcplx-1;itemp=cid_oqa_e
  cid_ohc_b=addone(itemp);cid_ohc_e=cid_ohc_b+jcplx-1;itemp=cid_ohc_e
  cid_ohn_b=addone(itemp);cid_ohn_e=cid_ohn_b+jcplx-1;itemp=cid_ohn_e
  cid_ohp_b=addone(itemp);cid_ohp_e=cid_ohp_b+jcplx-1;itemp=cid_ohp_e
  cid_oha_b=addone(itemp);cid_oha_e=cid_oha_b+jcplx-1;itemp=cid_oha_e
  cid_osc_b=addone(itemp);cid_osc_e=cid_osc_b+jsken*jcplx-1;itemp=cid_osc_e
  cid_osa_b=addone(itemp);cid_osa_e=cid_osa_b+jsken*jcplx-1;itemp=cid_osa_e
  cid_osn_b=addone(itemp);cid_osn_e=cid_osn_b+jsken*jcplx-1;itemp=cid_osn_e
  cid_osp_b=addone(itemp);cid_osp_e=cid_osp_b+jsken*jcplx-1;itemp=cid_osp_e
  cid_orc_b=addone(itemp);cid_orc_e=cid_orc_b+ndbiomcp*jcplx-1;itemp=cid_orc_e
  cid_orn_b=addone(itemp);cid_orn_e=cid_orn_b+ndbiomcp*jcplx-1;itemp=cid_orn_e
  cid_orp_b=addone(itemp);cid_orp_e=cid_orp_b+ndbiomcp*jcplx-1;itemp=cid_orp_e

  cid_mBiomeHeter_b=addone(itemp);cid_mBiomeHeter_e=cid_mBiomeHeter_b+NumPlantChemElms*NumLiveHeterBioms*jcplx
  itemp=cid_mBiomeHeter_e
  cid_mBiomeAutor_b=addone(itemp);cid_mBiomeAutor_e=cid_mBiomeAutor_b+NumPlantChemElms*NumLiveAutoBioms
  itemp=cid_mBiomeAutor_e

  fid_RO2EcoDmndPrev=addone(itemp)
  fid_RO2GasXchangePrev=addone(itemp)
  fid_RNH4EcoDmndSoilPrev=addone(itemp)
  fid_RNO3EcoDmndSoilPrev=addone(itemp)
  fid_RNO2EcoUptkSoilPrev=addone(itemp)
  fid_RN2OEcoUptkSoilPrev=addone(itemp)
  fid_RH2PO4EcoDmndSoilPrev=addone(itemp)
  fid_RH1PO4EcoDmndSoilPrev=addone(itemp)
  fid_RNH4EcoDmndBandPrev=addone(itemp)
  fid_RNO3EcoDmndBandPrev=addone(itemp)
  fid_RNO2EcoUptkBandPrev=addone(itemp)
  fid_RH2PO4EcoDmndBandPrev=addone(itemp)
  fid_RH1PO4EcoDmndBandPrev=addone(itemp)
  fid_RDOMEcoDmndPrev_b=addone(itemp);fid_RDOMEcoDmndPrev_e=fid_RDOMEcoDmndPrev_b+jcplx;itemp=fid_RDOMEcoDmndPrev_e
  fid_RAcetateEcoDmndPrev_b=addone(itemp);fid_RAcetateEcoDmndPrev_e=fid_RAcetateEcoDmndPrev_b+jcplx;itemp=fid_RAcetateEcoDmndPrev_e
  fid_RNH4DmndSoilHeter_b=addone(itemp);fid_RNH4DmndSoilHeter_e=fid_RNH4DmndSoilHeter_b+NumHetetr1MicCmplx*jcplx;itemp=fid_RNH4DmndSoilHeter_e
  fid_RNH4DmndBandHeter_b=addone(itemp);fid_RNH4DmndBandHeter_e=fid_RNH4DmndBandHeter_b+NumHetetr1MicCmplx*jcplx;itemp=fid_RNH4DmndBandHeter_e
  fid_RNO3DmndSoilHeter_b=addone(itemp);fid_RNO3DmndSoilHeter_e=fid_RNO3DmndSoilHeter_b+NumHetetr1MicCmplx*jcplx;itemp=fid_RNO3DmndSoilHeter_e
  fid_RNO3DmndBandHeter_b=addone(itemp);fid_RNO3DmndBandHeter_e=fid_RNO3DmndBandHeter_b+NumHetetr1MicCmplx*jcplx;itemp=fid_RNO3DmndBandHeter_e
  fid_RH2PO4DmndSoilHeter_b=addone(itemp);fid_RH2PO4DmndSoilHeter_e=fid_RH2PO4DmndSoilHeter_b+NumHetetr1MicCmplx*jcplx;itemp=fid_RH2PO4DmndSoilHeter_e
  fid_RH2PO4DmndBandHeter_b=addone(itemp);fid_RH2PO4DmndBandHeter_e=fid_RH2PO4DmndBandHeter_b+NumHetetr1MicCmplx*jcplx;itemp=fid_RH2PO4DmndBandHeter_e
  fid_RH1PO4DmndSoilHeter_b=addone(itemp);fid_RH1PO4DmndSoilHeter_e=fid_RH1PO4DmndSoilHeter_b+NumHetetr1MicCmplx*jcplx;itemp=fid_RH1PO4DmndSoilHeter_e
  fid_RH1PO4DmndBandHeter_b=addone(itemp);fid_RH1PO4DmndBandHeter_e=fid_RH1PO4DmndBandHeter_b+NumHetetr1MicCmplx*jcplx;itemp=fid_RH1PO4DmndBandHeter_e
  fid_RO2DmndHetert_b=addone(itemp);fid_RO2DmndHetert_e=fid_RO2DmndHetert_b+NumHetetr1MicCmplx*jcplx;itemp=fid_RO2DmndHetert_e

  fid_XCODFS=addone(itemp)
  fid_XCHDFS=addone(itemp)
  fid_XOXDFS=addone(itemp)
  fid_XNGDFS=addone(itemp)
  fid_XN2DFS=addone(itemp)
  fid_XN3DFS=addone(itemp)
  fid_XNBDFS=addone(itemp)
  fid_XHGDFS=addone(itemp)
  fid_XCODFG=addone(itemp)
  fid_XCHDFG=addone(itemp)
  fid_XOXDFG=addone(itemp)
  fid_XNGDFG=addone(itemp)
  fid_XN2DFG=addone(itemp)
  fid_XN3DFG=addone(itemp)
  fid_XNBDFG=addone(itemp)
  fid_XHGDFG=addone(itemp)
  fid_XCOFLG=addone(itemp)
  fid_XCHFLG=addone(itemp)
  fid_XOXFLG=addone(itemp)
  fid_XNGFLG=addone(itemp)
  fid_XN2FLG=addone(itemp)
  fid_XN3FLG=addone(itemp)
  fid_XHGFLG=addone(itemp)

  cidg_CO2=addone(itemp)
  cidg_CH4=addone(itemp)
  cid_OXYG=addone(itemp)
  cid_Z2GG=addone(itemp)
  cid_Z2OG=addone(itemp)
  cid_H2GG=addone(itemp)
  cid_ZNH3G=addone(itemp)

  cid_ZNH4B=addone(itemp)
  cid_ZNH4S=addone(itemp)
  cid_ZNH3B=addone(itemp)
  cid_ZNH3S=addone(itemp)
  cid_ZNO3B=addone(itemp)
  cid_ZNO3S=addone(itemp)
  cid_Z2GS =addone(itemp)
  cid_H1POB=addone(itemp)
  cid_H1PO4=addone(itemp)
  cid_ZNO2B=addone(itemp)
  cid_ZNO2S=addone(itemp)
  cid_H2POB=addone(itemp)
  cid_H2PO4=addone(itemp)
  cid_CCO2S=addone(itemp)
  cid_CNO2S=addone(itemp)
  cid_CNO2B=addone(itemp)
  cid_CZ2OS=addone(itemp)
  cid_Z2OS =addone(itemp)
  cid_COXYS=addone(itemp)
  cid_OXYS =addone(itemp)
  cid_COXYG=addone(itemp)
  cid_CZ2GG=addone(itemp)
  cid_CZ2GS=addone(itemp)
  cid_CH2GS=addone(itemp)
  cid_H2GS =addone(itemp)
  cid_CCH4G=addone(itemp)
  cid_CH4S =addone(itemp)
  cid_ZNFN0=addone(itemp)
  cid_ZNFNI=addone(itemp)

  nvars=itemp
  end associate
  end subroutine Initboxbgc
! ----------------------------------------------------------------------
  subroutine UpdateStateVars(micfor, micstt,micflx,nvars,ystates0l,ystatesfl)
!
! DESCRIPTION
  implicit none
  type(micforctype), intent(in) :: micfor
  type(micfluxtype), intent(in) :: micflx
  type(micsttype)  , intent(in) :: micstt
  integer , intent(in) :: nvars
  real(r8), intent(in) :: ystates0l(nvars)
  real(r8), intent(inout) :: ystatesfl(nvars)

  integer :: K,N,NGL,M
  associate(                                                   &
    jcplx                   => micpar%jcplx,                   &
    NumMicbAFunGrupsPerCmplx => micpar%NumMicbAFunGrupsPerCmplx, &
    NumMicbHFunGrupsPerCmplx => micpar%NumMicbHFunGrupsPerCmplx, &    
    jsken                   => micpar%jsken,                   &
    k_humus                 => micpar%k_humus,                 &
    k_POM                   => micpar%k_POM,                   &
    iprotein                => micpar%iprotein,                &
    icarbhyro               => micpar%icarbhyro,               &
    nlbiomcp                => micpar%nlbiomcp,                &
    ndbiomcp                => micpar%ndbiomcp,                &
    is_litter               => micpar%is_litter,               &
    NumLiveAutoBioms        => micpar%NumLiveAutoBioms,        &
    NumLiveHeterBioms       => micpar%NumLiveHeterBioms,       &
    NumHetetr1MicCmplx       => micpar%NumHetetr1MicCmplx,       &
    NumMicrobAutoTrophCmplx   => micpar%NumMicrobAutoTrophCmplx,   &
    VLWatMicP               => micfor%VLWatMicP                &
  )
!atmospheric gaseous CO2,CH4,O2,NH3,N2,N2O,H2
!
  ystatesfl(cid_ZNH3B)=ystates0l(cid_ZNH3B)+ystatesfl(fid_TRChem_NH3_band_soil)+micflx%RNH4MicbReliz2Band
  ystatesfl(cid_ZNH3S)=ystates0l(cid_ZNH3S)+ystatesfl(fid_TRChem_NH3_soil_vr)+micflx%RNH4MicbReliz2Soil
  ystatesfl(cid_ZNH4B)=ystates0l(cid_ZNH4B)+ystatesfl(fid_TRChem_NH3_band_soil)+micflx%RNH4MicbReliz2Band
  ystatesfl(cid_ZNH4S)=ystates0l(cid_ZNH4S)+ystatesfl(fid_TRChem_NH4_soil)+micflx%RNH4MicbReliz2Soil
  ystatesfl(cid_H1POB)=ystates0l(cid_H1POB)+ystatesfl(fid_TRChem_H1PO4_band_soil)+micflx%RH1PO4MicbReliz2Band
  ystatesfl(cid_H1PO4)=ystates0l(cid_H1PO4)+ystatesfl(fid_TRChem_H1PO4_soil)+micflx%RH1PO4MicbReliz2Soil
  ystatesfl(cid_H2POB)=ystates0l(cid_H2POB)+ystatesfl(fid_TRChem_H2PO4_band_soil)+micflx%RH2PO4MicbReliz2Band
  ystatesfl(cid_H2PO4)=ystates0l(cid_H2PO4)+ystatesfl(fid_TRChem_H2PO4_soil)+micflx%RH2PO4MicbReliz2Soil
  ystatesfl(cid_ZNO3B)=ystates0l(cid_ZNO3B)+micflx%RNO3MicbReliz2Band
  ystatesfl(cid_ZNO3S)=ystates0l(cid_ZNO3S)+micflx%RNO3MicbReliz2Soil
  ystatesfl(cid_ZNO2B)=ystates0l(cid_ZNO2B)+micflx%RNO2MicbReliz2Band
  ystatesfl(cid_ZNO2S)=ystates0l(cid_ZNO2S)+micflx%RNO2MicbReliz2Soil

  ystatesfl(cid_CO2S) =ystates0l(cid_CO2S)-micflx%RCO2NetUptkMicb
  ystatesfl(cid_Z2OS) =ystates0l(cid_Z2OS)-micflx%RN2ONetUptkMicb
  ystatesfl(cid_OXYS) =ystates0l(cid_OXYS)-micflx%RO2UptkMicb
  ystatesfl(cid_H2GS) =ystates0l(cid_H2GS)-micflx%RH2NetUptkMicb
  ystatesfl(cid_CH4S) =ystates0l(cid_CH4S)-micflx%RCH4UptkAutor
  ystatesfl(cid_Z2GS) =ystates0l(cid_Z2GS)-micflx%RN2NetUptkMicb-micflx%MicrbN2Fix
  ystatesfl(cid_ZNFN0)=micstt%ZNFN0
  ystatesfl(cid_ZNFNI)=micstt%ZNFNI

  ystatesfl(cid_CO2S)=ystatesfl(cid_CO2S)+ystatesfl(fid_XCODFS)+ystatesfl(fid_XCODFG)
  ystatesfl(cid_CH4S)=ystatesfl(cid_CH4S)+ystatesfl(fid_XCHDFS)+ystatesfl(fid_XCHDFG)
  ystatesfl(cid_OXYS)=ystatesfl(cid_OXYS)+ystatesfl(fid_XOXDFS)+ystatesfl(fid_XOXDFG)
  ystatesfl(cid_Z2GS)=ystatesfl(cid_Z2GS)+ystatesfl(fid_XNGDFS)+ystatesfl(fid_XNGDFG)
  ystatesfl(cid_Z2OS)=ystatesfl(cid_Z2OS)+ystatesfl(fid_XN2DFS)+ystatesfl(fid_XN2DFG)
  ystatesfl(cid_ZNH3S)=ystatesfl(cid_ZNH3S)+ystatesfl(fid_XN3DFS)+ystatesfl(fid_XN3DFG)
  ystatesfl(cid_ZNH3B)=ystatesfl(cid_ZNH3B)+ystatesfl(fid_XNBDFS)+ystatesfl(fid_XNBDFG)
  ystatesfl(cid_H2GS)=ystatesfl(cid_H2GS)+ystatesfl(fid_XHGDFS)+ystatesfl(fid_XHGDFG)

  ystatesfl(cidg_CO2)=ystates0l(cidg_CO2)-ystatesfl(fid_XCODFG)+ystatesfl(fid_XHGDFG)
  ystatesfl(cidg_CH4)=ystates0l(cidg_CH4)-ystatesfl(fid_XCHDFG)+ystatesfl(fid_XCHFLG)
  ystatesfl(cid_OXYG)=ystates0l(cid_OXYG)-ystatesfl(fid_XOXDFG)+ystatesfl(fid_XOXFLG)
  ystatesfl(cid_Z2GG)=ystates0l(cid_Z2GG)-ystatesfl(fid_XNGDFG)+ystatesfl(fid_XNGFLG)
  ystatesfl(cid_Z2OG)=ystates0l(cid_Z2OG)-ystatesfl(fid_XN2DFG)+ystatesfl(fid_XN2FLG)
  ystatesfl(cid_ZNH3G)=ystates0l(cid_ZNH3G)-ystatesfl(fid_XN3DFG)-ystatesfl(fid_XNBDFG) &
    +ystatesfl(fid_XN3FLG)
  ystatesfl(cid_H2GG)=ystates0l(cid_H2GG)-ystatesfl(fid_XHGDFG)+ystatesfl(fid_XHGFLG)
  ystatesfl(fid_RO2GasXchangePrev)=ystatesfl(fid_XOXDFG)

  ystatesfl(cid_CNO2S)=ystatesfl(cid_ZNO2S)/(VLWatMicP*micfor%VLNO3)
  if(micfor%VLNOB>0._r8)ystatesfl(cid_CNO2B)=ystatesfl(cid_ZNO2B)/(VLWatMicP*micfor%VLNOB)
  ystatesfl(cid_CCO2S)=ystatesfl(cid_CO2S)/micfor%VLWatMicP
  ystatesfl(cid_CZ2OS)=ystatesfl(cid_Z2OS)/micfor%VLWatMicP
  ystatesfl(cid_CH2GS)=ystatesfl(cid_H2GS)/micfor%VLWatMicP
  ystatesfl(cid_COXYS)=ystatesfl(cid_OXYS)/micfor%VLWatMicP
  ystatesfl(cid_CZ2GS)=ystatesfl(cid_Z2GS)/micfor%VLWatMicP
  ystatesfl(cid_CZ2GG)=ystatesfl(cid_Z2GG)/micfor%VLsoiAirP
  ystatesfl(cid_COXYG)=ystatesfl(cid_OXYG)/micfor%VLsoiAirP
  ystatesfl(cid_CCH4G)=ystatesfl(cidg_CH4)/micfor%VLsoiAirP

! the following variables are updated in the microbial model
  ystatesfl(cid_oqc_b:cid_oqc_e)=micstt%DOM(idom_doc,1:jcplx)
  ystatesfl(cid_oqn_b:cid_oqn_e)=micstt%DOM(idom_don,1:jcplx)
  ystatesfl(cid_oqp_b:cid_oqp_e)=micstt%DOM(idom_dop,1:jcplx)
  ystatesfl(cid_oqa_b:cid_oqa_e)=micstt%DOM(idom_acetate,1:jcplx)
  ystatesfl(cid_ohc_b:cid_ohc_e)=micstt%SorbedOM(ielmc,1:jcplx)
  ystatesfl(cid_ohn_b:cid_ohn_e)=micstt%SorbedOM(ielmn,1:jcplx)
  ystatesfl(cid_ohp_b:cid_ohp_e)=micstt%SorbedOM(ielmp,1:jcplx)
  ystatesfl(cid_oha_b:cid_oha_e)=micstt%SorbedOM(idom_acetate,1:jcplx)
  ystatesfl(cid_osc_b:cid_osc_e)=reshape(micstt%SolidOM(ielmc,1:jsken,1:jcplx),(/jsken*jcplx/))
  ystatesfl(cid_osa_b:cid_osa_e)=reshape(micstt%SolidOMAct(1:jsken,1:jcplx),(/jsken*jcplx/))
  ystatesfl(cid_osn_b:cid_osn_e)=reshape(micstt%SolidOM(ielmn,1:jsken,1:jcplx),(/jsken*jcplx/))
  ystatesfl(cid_osp_b:cid_osp_e)=reshape(micstt%SolidOM(ielmp,1:jsken,1:jcplx),(/jsken*jcplx/))
  ystatesfl(cid_orc_b:cid_orc_e)=reshape(micstt%OMBioResdu(ielmc,1:ndbiomcp,1:jcplx),(/ndbiomcp*jcplx/))
  ystatesfl(cid_orn_b:cid_orn_e)=reshape(micstt%OMBioResdu(ielmn,1:ndbiomcp,1:jcplx),(/ndbiomcp*jcplx/))
  ystatesfl(cid_orp_b:cid_orp_e)=reshape(micstt%OMBioResdu(ielmp,1:ndbiomcp,1:jcplx),(/ndbiomcp*jcplx/))
  ystatesfl(cid_mBiomeHeter_b:cid_mBiomeHeter_e)=reshape(micstt%mBiomeHeter(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx),&
    (/NumPlantChemElms*NumLiveHeterBioms*jcplx/))
  ystatesfl(cid_mBiomeAutor_b:cid_mBiomeAutor_e)=reshape(micstt%mBiomeAutor(1:NumPlantChemElms,1:NumLiveAutoBioms),&
    (/NumPlantChemElms*NumLiveAutoBioms/))

! summarize diagnostic fluxes
  DO K=1,jcplx
    IF(.not.micfor%litrm.or.(micpar%is_litter(K)))THEN
      DO N=1,NumMicbHFunGrupsPerCmplx
        DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
          ystatesfl(fid_RO2EcoDmndPrev)=ystatesfl(fid_RO2EcoDmndPrev)+micflx%RO2DmndHetert(NGL,K)
          ystatesfl(fid_RNH4EcoDmndSoilPrev)=ystatesfl(fid_RNH4EcoDmndSoilPrev)+micflx%RNH4DmndSoilHeter(NGL,K)
          ystatesfl(fid_RNO3EcoDmndSoilPrev)=ystatesfl(fid_RNO3EcoDmndSoilPrev)+micflx%RNO3ReduxDmndSoilHeter(NGL,K)+micflx%RNO3DmndSoilHeter(NGL,K)
          ystatesfl(fid_RNO2EcoUptkSoilPrev)=ystatesfl(fid_RNO2EcoUptkSoilPrev)+micflx%RNO2DmndReduxSoilHeter(NGL,K)
          ystatesfl(fid_RN2OEcoUptkSoilPrev)=ystatesfl(fid_RN2OEcoUptkSoilPrev)+micflx%RN2ODmndReduxHeter(NGL,K)
          ystatesfl(fid_RH2PO4EcoDmndSoilPrev)=ystatesfl(fid_RH2PO4EcoDmndSoilPrev)+micflx%RH2PO4DmndSoilHeter(NGL,K)
          ystatesfl(fid_RH1PO4EcoDmndSoilPrev)=ystatesfl(fid_RH1PO4EcoDmndSoilPrev)+micflx%RH1PO4DmndSoilHeter(NGL,K)
          ystatesfl(fid_RNH4EcoDmndBandPrev)=ystatesfl(fid_RNH4EcoDmndBandPrev)+micflx%RNH4DmndBandHeter(NGL,K)
          ystatesfl(fid_RNO3EcoDmndBandPrev)=ystatesfl(fid_RNO3EcoDmndBandPrev)+micflx%RNO3ReduxDmndBandHeter(NGL,K)+micflx%RNO3DmndBandHeter(NGL,K)
          ystatesfl(fid_RNO2EcoUptkBandPrev)=ystatesfl(fid_RNO2EcoUptkBandPrev)+micflx%RNO2DmndReduxBandHeter(NGL,K)
          ystatesfl(fid_RH2PO4EcoDmndBandPrev)=ystatesfl(fid_RH2PO4EcoDmndBandPrev)+micflx%RH2PO4DmndBandHeter(NGL,K)
          ystatesfl(fid_RH1PO4EcoDmndBandPrev)=ystatesfl(fid_RH1PO4EcoDmndBandPrev)+micflx%RH1PO4DmndBandHeter(NGL,K)
          ystatesfl(fid_RDOMEcoDmndPrev_b+K)=ystatesfl(fid_RDOMEcoDmndPrev_b+K)+micflx%RDOCUptkHeter(NGL,K)
          ystatesfl(fid_RAcetateEcoDmndPrev_b+K)=ystatesfl(fid_RAcetateEcoDmndPrev_b+K)+micflx%RAcetateUptkHeter(NGL,K)
        enddo
      ENDDO
    ENDIF
  ENDDO

  DO  N=1,NumMicbAFunGrupsPerCmplx
    DO NGL=micpar%JGniA(N),micpar%JGnfA(N)
      ystatesfl(fid_RO2EcoDmndPrev)=ystatesfl(fid_RO2EcoDmndPrev)+micflx%RO2MetaDmndAutor(NGL)
      ystatesfl(fid_RNH4EcoDmndSoilPrev)=ystatesfl(fid_RNH4EcoDmndSoilPrev)+micflx%RNH3OxidAutor(NGL)+micflx%RNH4UptkSoilAutor(NGL)
      ystatesfl(fid_RNO3EcoDmndSoilPrev)=ystatesfl(fid_RNO3EcoDmndSoilPrev)+micflx%RNO3UptkSoilAutor(NGL)
      ystatesfl(fid_RNO2EcoUptkSoilPrev)=ystatesfl(fid_RNO2EcoUptkSoilPrev)+micflx%RNO2XupAutor(NGL)
      ystatesfl(fid_RH2PO4EcoDmndSoilPrev)=ystatesfl(fid_RH2PO4EcoDmndSoilPrev)+micflx%RH2PO4UptkSoilAutor(NGL)
      ystatesfl(fid_RH1PO4EcoDmndSoilPrev)=ystatesfl(fid_RH1PO4EcoDmndSoilPrev)+micflx%RH1PO4UptkSoilAutor(NGL)
      ystatesfl(fid_RNH4EcoDmndBandPrev)=ystatesfl(fid_RNH4EcoDmndBandPrev)+micflx%RNH3OxidAutorBand(NGL)+micflx%RNH4UptkBandAutor(NGL)
      ystatesfl(fid_RNO3EcoDmndBandPrev)=ystatesfl(fid_RNO3EcoDmndBandPrev)+micflx%RNO3UptkBandAutor(NGL)
      ystatesfl(fid_RNO2EcoUptkBandPrev)=ystatesfl(fid_RNO2EcoUptkBandPrev)+micflx%RNO2XupAutorBand(NGL)
      ystatesfl(fid_RH2PO4EcoDmndBandPrev)=ystatesfl(fid_RH2PO4EcoDmndBandPrev)+micflx%RH2PO4UptkBandAutor(NGL)
      ystatesfl(fid_RH1PO4EcoDmndBandPrev)=ystatesfl(fid_RH1PO4EcoDmndBandPrev)+micflx%RH1PO4UptkBandAutor(NGL)
    enddo
  ENDDO

  ystatesfl(fid_RNO2EcoUptkSoilPrev)=ystatesfl(fid_RNO2EcoUptkSoilPrev)+micflx%RNO2DmndSoilChemo
  ystatesfl(fid_RNO2EcoUptkBandPrev)=ystatesfl(fid_RNO2EcoUptkBandPrev)+micflx%RNO2DmndBandChemo

  end associate
  end subroutine UpdateStateVars

! ----------------------------------------------------------------------

  subroutine getvarlist(nvars, varl, varlnml, unitl, vartypes)

  use bhistMod, only : hist_var_str_len,hist_unit_str_len,hist_var_lon_str_len
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(nvars)            !variable name
  character(len=hist_var_lon_str_len), intent(out) :: varlnml(nvars)     !variable name
  character(len=hist_unit_str_len),intent(out) :: unitl(nvars)           !variable unit
  integer                         ,intent(out) :: vartypes(nvars)        !variable type, flx or state

  integer :: iknen,icplx
  integer :: jj,ll,k,m,n,ngl

  associate(                                                   &
    jcplx                   => micpar%jcplx,                   &
    jsken                   => micpar%jsken,                   &
    NumMicbAFunGrupsPerCmplx => micpar%NumMicbAFunGrupsPerCmplx, &
    NumMicbHFunGrupsPerCmplx => micpar%NumMicbHFunGrupsPerCmplx, &    
    nlbiomcp                => micpar%nlbiomcp,                &
    ndbiomcp                => micpar%ndbiomcp                 &
  )

  !configure variables
!  call getchemvarlist(nvars, varl, varlnml, unitl, vartypes)
  call getvarlist_nosalt(nvars, varl, varlnml, unitl, vartypes)

  varl(cid_ZMG) ='ZMG';varlnml(cid_ZMG)='soil aqueous Mg content micropore'
  unitl(cid_ZMG)='mol d-2';vartypes(cid_ZMG)=var_state_type

  varl(cid_ZNA) ='ZNA';varlnml(cid_ZNA)='soil aqueous Na content micropore'
  unitl(cid_ZNA) ='mol d-2';vartypes(cid_ZNA)=var_state_type

  varl(cid_ZKA) ='ZKA';varlnml(cid_ZKA)='soil aqueous K content micropore'
  unitl(cid_ZKA) ='mol d-2';vartypes(cid_ZKA)=var_state_type

  varl(cid_CO2S)='CO2S';varlnml(cid_CO2S)='aqueous CO2 concentration micropore';
  unitl(cid_CO2S)='gC d-2';vartypes(cid_CO2S)=var_state_type

  varl(cid_ZNH3B)='ZNH3B';varlnml(cid_ZNH3B)='band soil micropore NH3 mass'
  unitl(cid_ZNH3B)='gN d-2';vartypes(cid_ZNH3B)=var_state_type

  varl(cid_ZNH3S)='ZNH3S';varlnml(cid_ZNH3S)='non-band soil micropore NH3 mass'
  unitl(cid_ZNH3S)='gN d-2';vartypes(cid_ZNH3S)=var_state_type

  varl(cid_ZNH4B)='ZNH4B';varlnml(cid_ZNH4B)='band soil micropore NH4(+) mass'
  unitl(cid_ZNH4B)='gN d-2';vartypes(cid_ZNH4B)=var_state_type

  varl(cid_ZNH4S)='ZNH4S';varlnml(cid_ZNH4S)='non-band soil micropore NH4(+) mass'
  unitl(cid_ZNH4S)='gN d-2';vartypes(cid_ZNH4S)=var_state_type

  varl(cid_ZNO3B)='ZNO3B';varlnml(cid_ZNO3B)='band soil micropore NO3(-) mass'
  unitl(cid_ZNO3B)='gN d-2';vartypes(cid_ZNO3B)=var_state_type

  varl(cid_ZNO3S)='ZNO3S';varlnml(cid_ZNO3S)='non-band soil micropore NO3(-) mass'
  unitl(cid_ZNO3S)='gN d-2';vartypes(cid_ZNO3S)=var_state_type

  varl(cid_H1POB)='H1POB';varlnml(cid_H1POB)='band soil micropore aqueous HPO4 content'
  unitl(cid_H1POB)='gP m-2';vartypes(cid_H1POB)=var_state_type

  varl(cid_H1PO4)='H1PO4';varlnml(cid_H1PO4)='non-band soil micropore aqueous HPO4(--) content';
  unitl(cid_H1PO4)='gP m-2';vartypes(cid_H1PO4)=var_state_type

  varl(cid_ZNO2B)='ZNO2B';varlnml(cid_ZNO2B)='band soil micropore NO2(-) mass'
  unitl(cid_ZNO2B)='gN m-2';vartypes(cid_ZNO2B)=var_state_type

  varl(cid_ZNO2S)='ZNO2S';varlnml(cid_ZNO2S)='non-band soil micropore NO2(-) mass'
  unitl(cid_ZNO2S)='gN m-2';vartypes(cid_ZNO2S)=var_state_type

  varl(cid_H2POB)='H2POB';varlnml(cid_H2POB)='band soil micropore H2PO4 mass'
  unitl(cid_H2POB)='gP m-2';vartypes(cid_H2POB)=var_state_type

  varl(cid_H2PO4)='H2PO4';varlnml(cid_H2PO4)='non-band soil micropore H2PO4 mass'
  unitl(cid_H2PO4)='gP m-2';vartypes(cid_H2PO4)=var_state_type

  varl(cid_CCO2S)='CCO2S';varlnml(cid_CCO2S)='soil micropore aqueous CO2 concentration'
  unitl(cid_CCO2S)='gC m-3';vartypes(cid_CCO2S)=var_state_type

  varl(cid_CNO2S)='CNO2S';varlnml(cid_CNO2S)='non-band soil micropore NO2 concentration'
  unitl(cid_CNO2S)='gN m-3';vartypes(cid_CNO2S)=var_state_type

  varl(cid_CNO2B)='CNO2B';varlnml(cid_CNO2B)='band soil micropore NO2 concentration'
  unitl(cid_CNO2B)='gN m-3';vartypes(cid_CNO2B)=var_state_type

  varl(cid_CZ2OS)='CZ2OS';varlnml(cid_CZ2OS)='soil micropore aqueous N2O concentration'
  unitl(cid_CZ2OS)='gN m-3';vartypes(cid_CZ2OS)=var_state_type

  varl(cid_Z2OS) ='Z2OS';varlnml(cid_Z2OS)='soil micropore aqueous N2O mass'
  unitl(cid_Z2OS)='gN d-2';vartypes(cid_Z2OS)=var_state_type

  varl(cid_COXYS)='COXYS';varlnml(cid_COXYS)='soil micropore aqueous O2 concentration'
  unitl(cid_COXYS)='g m-3';vartypes(cid_COXYS)=var_state_type


  varl(cid_OXYS) ='OXYS';varlnml(cid_OXYS)='soil micropore aqueous O2 mass'
  unitl(cid_OXYS)='g d-2';vartypes(cid_OXYS)=var_state_type

  varl(cid_COXYG)='COXYG';varlnml(cid_COXYG)='soil micropore gaseous O2 concentration'
  unitl(cid_COXYG)='g m-3';vartypes(cid_COXYG)=var_state_type

  varl(cid_CZ2GG)='CZ2GG';varlnml(cid_CZ2GG)='soil micropore gaseous N2 concentration'
  unitl(cid_CZ2GG)='g m-3';vartypes(cid_CZ2GG)=var_state_type

  varl(cid_Z2GS)='Z2GS';varlnml(cid_Z2GS)='soil micropore aqueous N2 mass'
  unitl(cid_Z2GS)='gN d-2';vartypes(cid_Z2GS)=var_state_type

  varl(cid_CZ2GS)='CZ2GS';varlnml(cid_CZ2GS)='soil micropore aqueous N2 concentration'
  unitl(cid_CZ2GS)='gN m-3';vartypes(cid_CZ2GS)=var_state_type

  varl(cid_CH2GS)='CH2GS';varlnml(cid_CH2GS)='soil micropore aqueous H2 concentration'
  unitl(cid_CH2GS)='g m-3';vartypes(cid_CH2GS)=var_state_type

  varl(cid_H2GS) ='H2GS';varlnml(cid_H2GS)='soil micropore aqueous H2 mass'
  unitl(cid_H2GS)='g d-2';vartypes(cid_H2GS)=var_state_type

  varl(cid_CCH4G)='CCH4G';varlnml(cid_CCH4G)='soil micropore gaseous CH4 concentration'
  unitl(cid_CCH4G)='gC m-3';vartypes(cid_CCH4G)=var_state_type

  varl(cid_CH4S) ='CH4S';varlnml(cid_CH4S)='soil micropore aqueous CH4 mass'
  unitl(cid_CH4S)='gC d-2';vartypes(cid_CH4S)=var_state_type

  varl(cid_ZNFN0)='ZNFN0';varlnml(cid_ZNFN0)='initial nitrification inhibition activity'
  unitl(cid_ZNFN0)='none';vartypes(cid_ZNFN0)=var_state_type

  varl(cid_ZNFNI)='ZNFNI';varlnml(cid_ZNFNI)='current nitrification inhibition activity'
  unitl(cid_ZNFNI)='none';vartypes(cid_ZNFNI)=var_state_type

  varl(cidg_CO2)='CO2G';varlnml(cidg_CO2)='gaseous CO2 mass'
  unitl(cidg_CO2)='gC d-2';vartypes(cidg_CO2)=var_state_type

  varl(cidg_CH4)='CH4G';varlnml(cidg_CH4)='gaseous CH4 mass'
  unitl(cidg_CH4)='gC d-2';vartypes(cidg_CH4)=var_state_type

  varl(cid_OXYG)='OXYG';varlnml(cid_OXYG)='gaseous oxygen mass'
  unitl(cid_OXYG)='gO d-2';vartypes(cid_OXYG)=var_state_type

  varl(cid_Z2GG)='Z2GG';varlnml(cid_Z2GG)='gaseous N2 mass'
  unitl(cid_Z2GG)='gN d-2';vartypes(cid_Z2GG)=var_state_type

  varl(cid_Z2OG)='Z2OG';varlnml(cid_Z2OG)='gaseous N2O mass'
  unitl(cid_Z2OG)='gN d-2';vartypes(cid_Z2OG)=var_state_type

  varl(cid_H2GG)='H2GG';varlnml(cid_H2GG)='gaseous H2 mass'
  unitl(cid_H2GG)='gH d-2';vartypes(cid_H2GG)=var_state_type

  varl(cid_ZNH3G)='ZNH3G';varlnml(cid_ZNH3G)='gaseous NH3 mass'
  unitl(cid_ZNH3G)='gN d-2';vartypes(cid_ZNH3G)=var_state_type

  do jj=cid_oqc_b,cid_oqc_e
    write(varl(jj),'(A,I1)')'OQC',jj-cid_oqc_b
    varlnml(jj)='micropore dissolved organic C mass in complex '//trim(micpar%cplxname(jj-cid_oqc_b+1))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_oqn_b,cid_oqn_e
    write(varl(jj),'(A,I1)')'OQN',jj-cid_oqn_b
    varlnml(jj)='micropore dissolved N mass in complex '//trim(micpar%cplxname(jj-cid_oqn_b+1))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_oqp_b,cid_oqp_e
    write(varl(jj),'(A,I1)')'OQP',jj-cid_oqn_b
    varlnml(jj)='micropore dissolved N mass in complex '//trim(micpar%cplxname(jj-cid_oqp_b+1))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_oqa_b,cid_oqa_e
    write(varl(jj),'(A,I1)')'OQA',jj-cid_oqa_b
    varlnml(jj)='micropore dissolved acetate mass in complex '//trim(micpar%cplxname(jj-cid_oqa_b+1))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_ohc_b,cid_ohc_e
    write(varl(jj),'(A,I1)')'OHC',jj-cid_ohc_b
    varlnml(jj)='adsorbed soil C mass in complex'//trim(micpar%cplxname(jj-cid_ohc_b+1))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_ohn_b,cid_ohn_e
    write(varl(jj),'(A,I1)')'OHN',jj-cid_ohn_b
    varlnml(jj)='adsorbed soil N mass in complex'//trim(micpar%cplxname(jj-cid_ohn_b+1))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_ohp_b,cid_ohp_e
    write(varl(jj),'(A,I1)')'OHP',jj-cid_ohp_b
    varlnml(jj)='adsorbed soil P mass in complex'//trim(micpar%cplxname(jj-cid_ohp_b+1))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_oha_b,cid_oha_e
    write(varl(jj),'(A,I1)')'OHA',jj-cid_oha_b
    varlnml(jj)='adsorbed soil acetate mass in complex'//trim(micpar%cplxname(jj-cid_oha_b+1))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo
  do jj=cid_osc_b,cid_osc_e
    iknen=jj-cid_osc_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=mod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSC',iknen+1,icplx
    varlnml(jj)='humus soil C as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx+1))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_osn_b,cid_osn_e
    iknen=jj-cid_osn_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=pmod(iknen,jsken);
    write(varl(jj),'(A,I1,I1)')'OSN',iknen+1,icplx
    varlnml(jj)='humus soil N as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_osp_b,cid_osp_e
    iknen=jj-cid_osp_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=pmod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSP',iknen+1,icplx
    varlnml(jj)='humus soil P as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_osa_b,cid_osa_e
    iknen=jj-cid_osa_b
    icplx=floor((iknen+1-1.e-3_r8)/jsken)
    iknen=pmod(iknen,jsken)
    write(varl(jj),'(A,I1,I1)')'OSA',iknen+1,icplx
    varlnml(jj)='colonized humus soil C as '//trim(micpar%kiname(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_orc_b,cid_orc_e
    iknen=jj-cid_orc_b
    icplx=floor((iknen+1-1.e-3_r8)/ndbiomcp)
    iknen=pmod(iknen,ndbiomcp)
    write(varl(jj),'(A,I1,I1)')'ORC',iknen+1,icplx
    varlnml(jj)='microbial residue C as '//trim(micpar%micresb(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gC d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_orn_b,cid_orn_e
    iknen=jj-cid_orn_b
    icplx=floor((iknen+1-1.e-3_r8)/ndbiomcp)
    iknen=pmod(iknen,ndbiomcp)
    write(varl(jj),'(A,I1,I1)')'ORN',iknen+1,icplx
    varlnml(jj)='microbial residue N as '//trim(micpar%micresb(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gN d-2'
    vartypes(jj)=var_state_type
  enddo

  do jj=cid_orp_b,cid_orp_e
    iknen=jj-cid_orp_b
    icplx=floor((iknen+1-1.e-3_r8)/ndbiomcp)
    iknen=pmod(iknen,ndbiomcp)
    write(varl(jj),'(A,I1,I1)')'ORP',iknen+1,icplx
    varlnml(jj)='microbial residue P as '//trim(micpar%micresb(iknen))//' in complex '//trim(micpar%cplxname(icplx))
    unitl(jj)='gP d-2'
    vartypes(jj)=var_state_type
  enddo

  jj=0
  DO k=1,jcplx
  DO N=1,NumMicbHFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
  DO M=1,nlbiomcp
    ll=cid_mBiomeHeter_b+jj;jj=jj+1
    write(varl(ll),'(A,I2.2,A)')'OMC'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass C in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))//' in complex '//trim(micpar%cplxname(k))
    unitl(ll)='gC d-2'
    vartypes(ll)=var_state_type

    ll=cid_mBiomeHeter_b+jj;jj=jj+1
    write(varl(ll),'(A,I2.2,A)')'OMN'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass N in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))//' in complex '//trim(micpar%cplxname(k))
    unitl(ll)='gN d-2'
    vartypes(ll)=var_state_type

    ll=cid_mBiomeHeter_b+jj;jj=jj+1
    write(varl(ll),'(A,I2.2,A)')'OMP'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))//trim(micpar%cplxname(k))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass P in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))//' in complex '//trim(micpar%cplxname(k))
    unitl(ll)='gP d-2'
    vartypes(ll)=var_state_type    
  enddo
  ENDDO
  ENDDO
  ENDDO

  jj=0
  DO N=1,NumMicbAFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
  DO M=1,nlbiomcp
    ll=cid_mBiomeAutor_b+jj;jj=jj+1
    write(varl(ll),'(A,I2.2,A)')'OMC'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass C in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))
    unitl(ll)='gC d-2'
    vartypes(ll)=var_state_type

    ll=cid_mBiomeAutor_b+jj;jj=jj+1
    write(varl(ll),'(A,I2.2,A)')'OMN'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass N in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))
    unitl(ll)='gN d-2'
    vartypes(ll)=var_state_type

    ll=cid_mBiomeAutor_b+jj;jj=jj+1
    write(varl(ll),'(A,I2.2,A)')'OMP'//trim(micpar%micbiom(M))//'g',NGL,&
      trim(micpar%hmicname(N))
    write(varlnml(ll),'(A,I2.2,A)')trim(micpar%micbiom(M))//' microbial biomass P in guild ',NGL,&
      ' of '//trim(micpar%hmicname(N))
    unitl(ll)='gP d-2'
    vartypes(ll)=var_state_type

  enddo
  ENDDO
  ENDDO

  varl(fid_RO2EcoDmndPrev)='RO2EcoDmndPrev';varlnml(fid_RO2EcoDmndPrev)='total root + microbial O2 uptake potential'
  unitl(fid_RO2EcoDmndPrev)='g d-2 h-1'; vartypes(fid_RO2EcoDmndPrev)=var_flux_type

  varl(fid_RO2GasXchangePrev)='RO2GasXchangePrev';varlnml(fid_RO2GasXchangePrev)='net gaseous O2 flux from previous hour'
  unitl(fid_RO2GasXchangePrev)='g d-2 h-1'; vartypes(fid_RO2GasXchangePrev)=var_flux_type

  varl(fid_RNH4EcoDmndSoilPrev)='RNH4EcoDmndSoilPrev';varlnml(fid_RNH4EcoDmndSoilPrev)='total root + microbial NH4 uptake potential non-band soil'
  unitl(fid_RNH4EcoDmndSoilPrev)='gN d-2 h-1'; vartypes(fid_RNH4EcoDmndSoilPrev)=var_flux_type

  varl(fid_RNO3EcoDmndSoilPrev)='RNO3EcoDmndSoilPrev';varlnml(fid_RNO3EcoDmndSoilPrev)='total root + microbial NO3 uptake potential non-band soil'
  unitl(fid_RNO3EcoDmndSoilPrev)='gN d-2 h-1'; vartypes(fid_RNO3EcoDmndSoilPrev)=var_flux_type

  varl(fid_RNO2EcoUptkSoilPrev)='RNO2EcoUptkSoilPrev';varlnml(fid_RNO2EcoUptkSoilPrev)='total root + microbial NO2 uptake potential non-band soil'
  unitl(fid_RNO2EcoUptkSoilPrev)='gN d-2 h-1'; vartypes(fid_RNO2EcoUptkSoilPrev)=var_flux_type

  varl(fid_RN2OEcoUptkSoilPrev)='RN2OEcoUptkSoilPrev';varlnml(fid_RN2OEcoUptkSoilPrev)='total root + microbial N2O uptake potential';
  unitl(fid_RN2OEcoUptkSoilPrev)='gN d-2 h-1'; vartypes(fid_RN2OEcoUptkSoilPrev)=var_flux_type

  varl(fid_RH2PO4EcoDmndSoilPrev)='RH2PO4EcoDmndSoilPrev';varlnml(fid_RH2PO4EcoDmndSoilPrev)='total root + microbial PO4 uptake potential non-band soil'
  unitl(fid_RH2PO4EcoDmndSoilPrev)='gP d-2 h-1'; vartypes(fid_RH2PO4EcoDmndSoilPrev)=var_flux_type

  varl(fid_RH1PO4EcoDmndSoilPrev)='RH1PO4EcoDmndSoilPrev';varlnml(fid_RH1PO4EcoDmndSoilPrev)='total root + microbial HPO4 uptake non-band soil'
  unitl(fid_RH1PO4EcoDmndSoilPrev)='gP d-2 h-1'; vartypes(fid_RH1PO4EcoDmndSoilPrev)=var_flux_type

  varl(fid_RNH4EcoDmndBandPrev)='RNH4EcoDmndBandPrev';varlnml(fid_RNH4EcoDmndBandPrev)='total root + microbial NH4 uptake potential band soil'
  unitl(fid_RNH4EcoDmndBandPrev)='gN d-2 h-1'; vartypes(fid_RNH4EcoDmndBandPrev)=var_flux_type

  varl(fid_RNO3EcoDmndBandPrev)='RNO3EcoDmndBandPrev';varlnml(fid_RNO3EcoDmndBandPrev)='total root + microbial NO3 uptake potential band soil'
  unitl(fid_RNO3EcoDmndBandPrev)='gN d-2 h-1'; vartypes(fid_RNO3EcoDmndBandPrev)=var_flux_type

  varl(fid_RNO2EcoUptkBandPrev)='RNO2EcoUptkBandPrev';varlnml(fid_RNO2EcoUptkBandPrev)='total root + microbial NO2 uptake potential band soil'
  unitl(fid_RNO2EcoUptkBandPrev)='gN d-2 h-1'; vartypes(fid_RNO2EcoUptkBandPrev)=var_flux_type

  varl(fid_RH2PO4EcoDmndBandPrev)='RH2PO4EcoDmndBandPrev';varlnml(fid_RH2PO4EcoDmndBandPrev)='total root + microbial PO4 uptake potential band soil'
  unitl(fid_RH2PO4EcoDmndBandPrev)='gP d-2 h-1'; vartypes(fid_RH2PO4EcoDmndBandPrev)=var_flux_type

  varl(fid_RH1PO4EcoDmndBandPrev)='RH1PO4EcoDmndBandPrev';varlnml(fid_RH1PO4EcoDmndBandPrev)='total root + microbial HPO4 uptake potential band soil';
  unitl(fid_RH1PO4EcoDmndBandPrev)='gP d-2 h-1'; vartypes(fid_RH1PO4EcoDmndBandPrev)=var_flux_type

  varl(fid_XCODFS)='XCODFS';varlnml(fid_XCODFS)='CO2 dissolution (+)-volatiziation (-) with respect to atmosphere'
  unitl(fid_XCODFS)='gC d-2 h-1'; vartypes(fid_XCODFS)=var_flux_type

  varl(fid_XCHDFS)='XCHDFS';varlnml(fid_XCHDFS)='CH4 dissolution (+)-volatiziation (-) with respect to atmosphere'
  unitl(fid_XCHDFS)='gC d-2 h-1';vartypes(fid_XCHDFS)=var_flux_type

  varl(fid_XOXDFS)='XOXDFS';varlnml(fid_XOXDFS)='O2 dissolution (+)-volatiziation (-) with respect to atmosphere'
  unitl(fid_XOXDFS)='gO d-2 h-1';vartypes(fid_XOXDFS)=var_flux_type

  varl(fid_XNGDFS)='XNGDFS';varlnml(fid_XNGDFS)='N2 dissolution (+)-volatiziation (-) with respect to atmosphere'
  unitl(fid_XNGDFS)='gN d-2 h-1';vartypes(fid_XNGDFS)=var_flux_type

  varl(fid_XN2DFS)='XN2DFS';varlnml(fid_XN2DFS)='N2O dissolution (+)-volatiziation (-) with respect to atmosphere'
  unitl(fid_XN2DFS)='gN d-2 h-1';vartypes(fid_XN2DFS)=var_flux_type

  varl(fid_XN3DFS)='XN3DFS';varlnml(fid_XN3DFS)='NH3 dissolution (+)-volatiziation (-) with respect to atmosphere in non-band soil'
  unitl(fid_XN3DFS)='gN d-2 h-1';vartypes(fid_XN3DFS)=var_flux_type

  varl(fid_XNBDFS)='XNBDFS';varlnml(fid_XNBDFS)='NH3 dissolution (+)-volatiziation (-) with respect to atmosphere in band soil'
  unitl(fid_XNBDFS)='gN d-2 h-1';vartypes(fid_XNBDFS)=var_flux_type

  varl(fid_XHGDFS)='XHGDFS';varlnml(fid_XHGDFS)='H2 dissolution (+)-volatiziation (-) with respect to atmosphere'
  unitl(fid_XHGDFS)='gH d-2 h-1';vartypes(fid_XHGDFS)=var_flux_type

  varl(fid_XCODFG)='XCODFG';varlnml(fid_XCODFG)='CO2 dissolution (+)-volatiziation (-) in soil'
  unitl(fid_XCODFG)='gC d-2 h-1';vartypes(fid_XCODFG)=var_flux_type

  varl(fid_XCHDFG)='XCHDFG';varlnml(fid_XCHDFG)='CH4 dissolution (+)-volatiziation (-) in soil'
  unitl(fid_XCHDFG)='gC d-2 h-1';vartypes(fid_XCHDFG)=var_flux_type

  varl(fid_XOXDFG)='XOXDFG';varlnml(fid_XOXDFG)='O2 dissolution (+)-volatiziation (-) in soil'
  unitl(fid_XOXDFG)='gO d-2 h-1';vartypes(fid_XOXDFG)=var_flux_type

  varl(fid_XNGDFG)='XNGDFG';varlnml(fid_XNGDFG)='N2 dissolution (+)-volatiziation (-) in soil'
  unitl(fid_XNGDFG)='gN d-2 h-1';vartypes(fid_XNGDFG)=var_flux_type

  varl(fid_XN2DFG)='XN2DFG';varlnml(fid_XN2DFG)='N2O dissolution (+)-volatiziation (-) in soil'
  unitl(fid_XN2DFG)='gN d-2 h-1';vartypes(fid_XN2DFG)=var_flux_type

  varl(fid_XN3DFG)='XN3DFG';varlnml(fid_XN3DFG)='NH3 dissolution (+)-volatiziation (-) in soil'
  unitl(fid_XN3DFG)='gN d-2 h-1';vartypes(fid_XN3DFG)=var_flux_type

  varl(fid_XNBDFG)='XNBDFG';varlnml(fid_XNBDFG)='NH3 dissolution (+)-volatiziation (-) in band soil'
  unitl(fid_XNBDFG)='gN d-2 h-1';vartypes(fid_XNBDFG)=var_flux_type

  varl(fid_XHGDFG)='XHGDFG';varlnml(fid_XHGDFG)='H2 dissolution (+)-volatiziation (-) in soil'
  unitl(fid_XHGDFG)='gH d-2 h-1';vartypes(fid_XHGDFG)=var_flux_type

  varl(fid_XCOFLG)='XCOFLG';varlnml(fid_XCOFLG)='CO2 gaseous exchange with atmosphere (-) into atmosphere'
  unitl(fid_XCOFLG)='gC d-2 h-1';vartypes(fid_XCOFLG)=var_flux_type

  varl(fid_XCHFLG)='XCHFLG';varlnml(fid_XCHFLG)='CH4 gaseous exchange with atmosphere (-) into atmosphere'
  unitl(fid_XCHFLG)='gC d-2 h-1';vartypes(fid_XCHFLG)=var_flux_type

  varl(fid_XOXFLG)='XOXFLG';varlnml(fid_XOXFLG)='O2 gaseous exchange with atmosphere (-) into atmosphere'
  unitl(fid_XOXFLG)='gO d-2 h-1';vartypes(fid_XOXFLG)=var_flux_type

  varl(fid_XNGFLG)='XNGFLG';varlnml(fid_XNGFLG)='N2 gaseous exchange with atmosphere (-) into atmosphere'
  unitl(fid_XNGFLG)='gN d-2 h-1';vartypes(fid_XNGFLG)=var_flux_type

  varl(fid_XN2FLG)='XN2FLG';varlnml(fid_XN2FLG)='N2O gaseous exchange with atmosphere (-) into atmosphere'
  unitl(fid_XN2FLG)='gN d-2 h-1';vartypes(fid_XN2FLG)=var_flux_type

  varl(fid_XN3FLG)='XN3FLG';varlnml(fid_XN3FLG)='N3H gaseous exchange with atmosphere (-) into atmosphere'
  unitl(fid_XN3FLG)='gN d-2 h-1';vartypes(fid_XN3FLG)=var_flux_type

  varl(fid_XHGFLG)='XHGFLG';varlnml(fid_XHGFLG)='H2 gaseous exchange with atmosphere (-) into atmosphere'
  unitl(fid_XHGFLG)='gH d-2 h-1';vartypes(fid_XHGFLG)=var_flux_type

  do jj =fid_RDOMEcoDmndPrev_b,fid_RDOMEcoDmndPrev_e
    write(varl(jj),'(A,I2.2)')'RDOMEcoDmndPrev',jj-fid_RDOMEcoDmndPrev_b
    varlnml(jj)='total root + microbial DOC uptake in complex ' &
      //micpar%cplxname(jj-fid_RDOMEcoDmndPrev_b+1)
    vartypes(jj)=var_flux_type
    unitl(jj)='gC d-2 h-1'
  enddo
  do jj =fid_RAcetateEcoDmndPrev_b,fid_RAcetateEcoDmndPrev_e
    write(varl(jj),'(A,I2.2)')'RAcetateEcoDmndPrev',jj-fid_RAcetateEcoDmndPrev_b
    varlnml(jj)='total root + microbial acetate uptake in complex ' &
      //micpar%cplxname(jj-fid_RAcetateEcoDmndPrev_b+1)
    vartypes(jj)=var_flux_type
    unitl(jj)='gC d-2 h-1'
  enddo

  ll=0
  DO k=1,jcplx
  DO N=1,NumMicbAFunGrupsPerCmplx
  DO NGL=micpar%JGniA(N),micpar%JGnfA(N)
    jj=fid_RNH4DmndSoilHeter_b+ll
    write(varl(jj),'(A,I2.2)')'RNH4DmndSoilHeter',ll
    varlnml(jj)='microbial NH4 demand in soil' //micpar%cplxname(k)
    vartypes(jj)=var_flux_type
    unitl(jj)='gN d-2 h-1'
    ll=ll+1
  enddo
  enddo
  enddo

  ll=0
  DO k=1,jcplx
  DO N=1,NumMicbHFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    jj=fid_RNH4DmndBandHeter_b+ll
    write(varl(jj),'(A,I2.2)')'RNH4DmndBandHeter',ll
    varlnml(jj)='microbial NH4 immobilization (+ve) - mineralization (-ve) band' &
      //micpar%cplxname(k)
    vartypes(jj)=var_flux_type
    unitl(jj)='gN d-2 h-1'
    ll=ll+1
  enddo
  enddo
  enddo

  ll=0
  DO k=1,jcplx
  DO N=1,NumMicbHFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    jj=fid_RNO3DmndSoilHeter_b+ll
    write(varl(jj),'(A,I2.2)')'RNO3DmndSoilHeter',ll
    varlnml(jj)='microbial NO3 demand in soil'//micpar%cplxname(k)
    vartypes(jj)=var_flux_type
    unitl(jj)='gN d-2 h-1'
    ll=ll+1
  enddo
  enddo
  enddo

  ll=0
  DO k=1,jcplx
  DO N=1,NumMicbHFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    jj=fid_RNO3DmndBandHeter_b+ll
    write(varl(jj),'(A,I2.2)')'RNO3DmndBandHeter',ll
    varlnml(jj)='microbial NO3 immobilization (+ve) - mineralization (-ve) band' &
      //micpar%cplxname(k)
    vartypes(jj)=var_flux_type
    unitl(jj)='gN d-2 h-1'
    ll=ll+1
  enddo
  enddo
  enddo

  ll=0
  DO k=1,jcplx
  DO N=1,NumMicbHFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    jj=fid_RH2PO4DmndSoilHeter_b+ll
    write(varl(jj),'(A,I2.2)')'RH2PO4DmndSoilHeter',ll
    varlnml(jj)='microbial PO4 demand in soil'//micpar%cplxname(k)
    vartypes(jj)=var_flux_type
    unitl(jj)='gP d-2 h-1'
    ll=ll+1
  enddo
  enddo
  enddo

  ll=0
  DO k=1,jcplx
  DO N=1,NumMicbHFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    jj=fid_RH2PO4DmndBandHeter_b+ll
    write(varl(jj),'(A,I2.2)')'RH2PO4DmndBandHeter',ll
    varlnml(jj)='substrate-unlimited H2PO4 mineralization-immobilization'//micpar%cplxname(k)
    vartypes(jj)=var_flux_type
    unitl(jj)='gP d-2 h-1'
    ll=ll+1
  enddo
  enddo
  enddo

  ll=0
  DO k=1,jcplx
  DO N=1,NumMicbHFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    jj=fid_RH1PO4DmndSoilHeter_b+ll
    write(varl(jj),'(A,I2.2)')'RH1PO4DmndSoilHeter',ll
    varlnml(jj)='substrate-unlimited HPO4 immobilization'//micpar%cplxname(k)
    vartypes(jj)=var_flux_type
    unitl(jj)='gP d-2 h-1'
    ll=ll+1
  enddo
  enddo
  enddo

  ll=0
  DO k=1,jcplx
  DO N=1,NumMicbHFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    jj=fid_RH1PO4DmndBandHeter_b+ll
    write(varl(jj),'(A,I2.2)')'RH1PO4DmndBandHeter',ll
    varlnml(jj)='substrate-unlimited HPO4 mineralization-immobilization'//micpar%cplxname(k)
    vartypes(jj)=var_flux_type
    unitl(jj)='gP d-2 h-1'
    ll=ll+1
  enddo
  enddo
  enddo

  ll=0
  DO k=1,jcplx
  DO N=1,NumMicbHFunGrupsPerCmplx
  DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
    jj=fid_RO2DmndHetert_b+ll
    write(varl(jj),'(A,I2.2)')'RO2DmndHetert',ll
    varlnml(jj)='aqueous O2 demand'//micpar%cplxname(k)
    vartypes(jj)=var_flux_type
    unitl(jj)='gO d-2 h-1'
    ll=ll+1
  enddo
  enddo
  enddo

  end associate
  end subroutine getvarlist
! ----------------------------------------------------------------------
  subroutine RunMicBGC(nvars, ystates0l, ystatesfl, forc,micfor,micstt,micflx, err_status)
!
!
  use ChemMod
  use MicBGCMod           , only : SoilBGCOneLayer
  use MicrobeDiagTypes,     only: Cumlate_Flux_Diag_type,Microbe_Diag_type
  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(in) :: ystates0l(nvars)
  real(r8), intent(out) :: ystatesfl(nvars)
  type(forc_type), intent(in) :: forc
  type(micforctype), intent(inout)    :: micfor
  type(micsttype)  , intent(inout) :: micstt
  type(micfluxtype), intent(inout) :: micflx
  type(model_status_type), intent(out) :: err_status
  type(Cumlate_Flux_Diag_type)  :: naqfdiag
  type(Microbe_Diag_type)  :: nmicdiag

  integer :: I,J

  I=1;J=1
  call err_status%reset()

  ystatesfl=0._r8
  if(.not.forc%disvolonly)then
!    print*,'SoilBGCOneLayer'
    call SoilBGCOneLayer(I,J,micfor,micstt,micflx,naqfdiag,nmicdiag)
!    print*,'RunModel_nosalt'
    call RunModel_nosalt(forc,micfor,nvars,ystates0l, ystatesfl, err_status)
  endif
!  print*,'CalcSurflux'
  call CalcSurflux(forc,micfor, nvars, ystates0l,ystatesfl,err_status)
!  print*,'UpdateStateVars'
  call UpdateStateVars(micfor,micstt,micflx,nvars,ystates0l,ystatesfl)
!
  call UpdateSOMORGM(micfor,micstt)
!  print*,'RunMicBGC'
  end subroutine RunMicBGC
! ----------------------------------------------------------------------
  subroutine CalcSurflux(forc,micfor, nvars, ystates0l,ystatesfl,err_status)
!
!  DESCRIPTION
! calcualte fluxes due to gaseous exchange with respect to atmosphere
! dissolution-volatilization between gasesous and disolved phases in soil
! dissolution-volatiziation with respect to atmopshere
! no micropore-macropore exchange is considered.
  use EcoSimConst
!
  implicit none
  type(forc_type), intent(in) :: forc
  type(micforctype), intent(in)    :: micfor
  integer,  intent(in)  :: nvars
  real(r8), intent(in)  :: ystates0l(nvars)
  real(r8), intent(inout) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status

  real(r8) :: DFVCOG,DFVCHG,DFVOXG,DFVNGG
  real(r8) :: DFVN2G,DFVN3G,DFVHGG
  real(r8) :: VLWatMicPMCO,VLWatMicPMCH,VLWatMicPMOX
  real(r8) :: VLWatMicPMNG,VLWatMicPMN2,VLWatMicPMN3
  real(r8) :: VLWatMicPMNB,VLWatMicPMHG,VOLCOT
  real(r8) :: VOLCHT,VOLOXT,VOLNGT
  real(r8) :: VOLN2T,VOLN3T,VOLBranchNumber_pft
  real(r8) :: VOLHGT
  real(r8) :: VLsoiAirPMA,VLsoiAirPMB
  real(r8) :: CCO2G2,CCH4G2,COXYG2
  real(r8) :: CZ2GG2,CZ2OG2,CNH3G2
  real(r8) :: CH2GG2
  real(r8) :: VLWatMicPMA,VLWatMicPMB
  real(r8) :: DLYR1,TortMicPM1
  real(r8) :: DiffusivitySolutEffCO,DiffusivitySolutEffCH
  real(r8) :: DiffusivitySolutEffOX,DiffusivitySolutEffNG
  real(r8) :: DiffusivitySolutEffN2,DiffusivitySolutEffN3
  real(r8) :: DiffusivitySolutEffHL
  real(r8) :: RCODXS,RCODFS
  real(r8) :: RCHDXS,RCHDFS
  real(r8) :: ROXDXS,ROXDFS
  real(r8) :: RNGDXS,RNGDFS
  real(r8) :: RN2DXS,RN2DFS
  real(r8) :: RN3DXS,RN3DFS
  real(r8) :: RNBDXS,RNBDFS
  real(r8) :: RHGDXS,RHGDFS
  real(r8) :: CCO2GQ,CCO2S2
  real(r8) :: CCH4GQ,CCH4S2
  real(r8) :: COXYGQ,COXYS2
  real(r8) :: CZ2GGQ,CZ2GS2
  real(r8) :: CZ2OGQ,CZ2OS2
  real(r8) :: CZN3GQ,CNH3S2
  real(r8) :: CZN3BQ,CNH3B2
  real(r8) :: CH2GGQ,CH2GS2
  real(r8) :: CNH4B2,CNH4S2
  real(r8) :: CH4G2,CH4S2
  real(r8) :: CNH3B0,CNH3S0
  real(r8) :: CO2G2,CO2S2
  real(r8) :: CNH4B0,CNH4S0
  real(r8) :: H2GG2,H2GS2
  real(r8) :: OXYG2,OXYS2
  real(r8) :: Z2OG2,Z2GG2
  real(r8) :: Z2GS2,Z2OS2,ZNH4S2
  real(r8) :: ZN3G2,ZNH3B2,ZNH3S2,ZNH4B2
  real(r8) :: RCHDFG,RCODFG,RHGDFG
  real(r8) :: RN2DFG,RN3DFG,RNBDFG
  real(r8) :: RNGDFG,ROXDFG
  real(r8) :: VLWatMicPMXA,VLWatMicPMXB
  real(r8) :: XCODFS,XCHDFS
  real(r8) :: XOXDFS,XNGDFS
  real(r8) :: XN2DFS,XN3DFS
  real(r8) :: XNBDFS,XHGDFS
  real(r8) :: XCODFG,XCHDFG
  real(r8) :: XOXDFG,XNGDFG
  real(r8) :: XN2DFG,XN3DFG
  real(r8) :: XNBDFG,XHGDFG
  real(r8) :: XCOFLG,XCHFLG
  real(r8) :: XOXFLG,XNGFLG
  real(r8) :: XN2FLG,XN3FLG
  real(r8) :: XHGFLG

  integer  :: MM,M
  call err_status%reset()

  associate(                   &
    ZEROS2 =>  micfor%ZEROS2 , &
    ZEROS  =>  micfor%ZEROS  , &
    VLWatMicP  =>  forc%VLWatMicP     , &
    VLsoiAirP  =>  forc%VLsoiAirP     , &
    PARG   =>  forc%PARG       &
  )

  XCODFS=0._r8
  XCHDFS=0._r8
  XOXDFS=0._r8
  XNGDFS=0._r8
  XN2DFS=0._r8
  XN3DFS=0._r8
  XNBDFS=0._r8
  XHGDFS=0._r8

  XCODFG=0._r8
  XCHDFG=0._r8
  XOXDFG=0._r8
  XNGDFG=0._r8
  XN2DFG=0._r8
  XN3DFG=0._r8
  XNBDFG=0._r8
  XHGDFG=0._r8

  XCOFLG=0._r8
  XCHFLG=0._r8
  XOXFLG=0._r8
  XNGFLG=0._r8
  XN2FLG=0._r8
  XN3FLG=0._r8
  XHGFLG=0._r8
  if(VLsoiAirP>ZEROS2)then
      !gaseous flux between atmosphere and soil

    CO2G2=AZMAX1(ystates0l(cidg_CO2))
    CH4G2=AZMAX1(ystates0l(cidg_CH4))
    OXYG2=AZMAX1(ystates0l(cid_OXYG))
    Z2GG2=AZMAX1(ystates0l(cid_Z2GG))
    Z2OG2=AZMAX1(ystates0l(cid_Z2OG))
    H2GG2=AZMAX1(ystates0l(cid_H2GG))
    ZN3G2=AZMAX1(ystates0l(cid_ZNH3G))
  endif

  if(VLWatMicP > micfor%ZEROS2)then
    CO2S2=AZMAX1(ystates0l(cid_CO2S))
    CH4S2=AZMAX1(ystates0l(cid_CH4S))
    OXYS2=AZMAX1(ystates0l(cid_OXYS))
    Z2GS2=AZMAX1(ystates0l(cid_Z2GS))
    Z2OS2=AZMAX1(ystates0l(cid_Z2OS))
    H2GS2=AZMAX1(ystates0l(cid_H2GS))
    ZNH3S2=AZMAX1(ystates0l(cid_ZNH3S))
    ZNH4S2=AZMAX1(ystates0l(cid_ZNH4S))
    ZNH3B2=AZMAX1(ystates0l(cid_ZNH3B))
    ZNH4B2=AZMAX1(ystates0l(cid_ZNH4B))

  endif
  DFVCOG=0._r8
  DFVCHG=0._r8
  DFVOXG=0._r8
  DFVNGG=0._r8
  DFVN2G=0._r8
  DFVN3G=0._r8
  DFVHGG=0._r8

  DO MM=1,NPG
    M=MIN(NPH,INT((MM-1)*dt_GasCyc)+1)

    if(VLsoiAirP>ZEROS2)then
      !gaseous flux between atmosphere and soil

      CCO2G2=AZMAX1(CO2G2/VLsoiAirP)
      CCH4G2=AZMAX1(CH4G2/VLsoiAirP)
      COXYG2=AZMAX1(OXYG2/VLsoiAirP)
      CZ2GG2=AZMAX1(Z2GG2/VLsoiAirP)
      CZ2OG2=AZMAX1(Z2OG2/VLsoiAirP)
      CNH3G2=AZMAX1(ZN3G2/VLsoiAirP)
      CH2GG2=AZMAX1(H2GG2/VLsoiAirP)

      DFVCOG=forc%DCO2GQ*(forc%CCO2E-CCO2G2)
      DFVCHG=forc%DCH4GQ*(forc%CCH4E-CCH4G2)
      DFVOXG=forc%DOXYGQ*(forc%COXYE-COXYG2)
      DFVNGG=forc%DZ2GGQ*(forc%CZ2GE-CZ2GG2)
      DFVN2G=forc%DZ2OGQ*(forc%CZ2OE-CZ2OG2)
      DFVN3G=forc%DNH3GQ*(forc%CNH3E-CNH3G2)
      DFVHGG=forc%DH2GGQ*(forc%CH2GE-CH2GG2)
      XCOFLG=XCOFLG+DFVCOG
      XCHFLG=XCHFLG+DFVCHG
      XOXFLG=XOXFLG+DFVOXG
      XNGFLG=XNGFLG+DFVNGG
      XN2FLG=XN2FLG+DFVN2G
      XN3FLG=XN3FLG+DFVN3G
      XHGFLG=XHGFLG+DFVHGG
    endif

    if(VLWatMicP > micfor%ZEROS2)then
  ! dissolution/volatilization is computed as the difference between current dissolved concentration and
  ! the atmospheric equilibrium concentration, with some dissolution rate
    !dissolution between atmosphere and soil
      VLWatMicPMA=VLWatMicP*forc%VLNH4
      VLWatMicPMB=VLWatMicP*forc%VLNHB

      CCO2S2=AZMAX1(CO2S2/VLWatMicP)
      CCH4S2=AZMAX1(CH4S2/VLWatMicP)
      COXYS2=AZMAX1(OXYS2/VLWatMicP)
      CZ2GS2=AZMAX1(Z2GS2/VLWatMicP)
      CZ2OS2=AZMAX1(Z2OS2/VLWatMicP)
      CH2GS2=AZMAX1(H2GS2/VLWatMicP)
      IF(VLWatMicPMA.GT.ZEROS2)THEN
        CNH3S2=AZMAX1(ZNH3S2/VLWatMicPMA)
        CNH4S2=AZMAX1(ZNH4S2/VLWatMicPMA)
      ELSE
        CNH3S2=0.0_r8
        CNH4S2=0.0_r8
      ENDIF
      IF(VLWatMicPMB.GT.ZEROS2)THEN
        CNH3B2=AZMAX1(ZNH3B2/VLWatMicPMB)
        CNH4B2=AZMAX1(ZNH4B2/VLWatMicPMB)
      ELSE
        CNH3B2=CNH3S2
        CNH4B2=CNH4S2
      ENDIF

!   between atmosphere and topsoil
      DLYR1=forc%DLYR3
      TortMicPM1=forc%TortMicPM*forc%AREA3/(0.5_r8*DLYR1)
      DiffusivitySolutEffCO=forc%CLSGL*TortMicPM1*dts_HeatWatTP
      DiffusivitySolutEffCH=forc%CQSGL*TortMicPM1*dts_HeatWatTP
      DiffusivitySolutEffOX=forc%O2AquaDiffusvity*TortMicPM1*dts_HeatWatTP
      DiffusivitySolutEffNG=forc%ZLSGL*TortMicPM1*dts_HeatWatTP
      DiffusivitySolutEffN2=forc%ZNSGL*TortMicPM1*dts_HeatWatTP
      DiffusivitySolutEffN3=forc%ZVSGL*TortMicPM1*dts_HeatWatTP
      DiffusivitySolutEffHL=forc%HLSGL*TortMicPM1*dts_HeatWatTP

      CCO2GQ=(PARG*forc%CCO2E*forc%SCO2L+DiffusivitySolutEffCO*CCO2S2)/(DiffusivitySolutEffCO+PARG)
      CCH4GQ=(PARG*forc%CCH4E*forc%CH4AquaSolubility+DiffusivitySolutEffCH*CCH4S2)/(DiffusivitySolutEffCH+PARG)
      COXYGQ=(PARG*forc%COXYE*forc%O2GSolubility+DiffusivitySolutEffOX*COXYS2)/(DiffusivitySolutEffOX+PARG)
      CZ2GGQ=(PARG*forc%CZ2GE*forc%SN2GL+DiffusivitySolutEffNG*CZ2GS2)/(DiffusivitySolutEffNG+PARG)
      CZ2OGQ=(PARG*forc%CZ2OE*forc%SN2OL+DiffusivitySolutEffN2*CZ2OS2)/(DiffusivitySolutEffN2+PARG)
      CZN3GQ=(PARG*forc%CNH3E*forc%SNH3L+DiffusivitySolutEffN3*CNH3S2)/(DiffusivitySolutEffN3+PARG)
      CZN3BQ=(PARG*forc%CNH3E*forc%SNH3L+DiffusivitySolutEffN3*CNH3B2)/(DiffusivitySolutEffN3+PARG)
      CH2GGQ=(PARG*forc%CH2GE*forc%SH2GL+DiffusivitySolutEffHL*CH2GS2)/(DiffusivitySolutEffHL+PARG)

      RCODFS=(CCO2GQ-CCO2S2)*AMIN1(VLWatMicP,DiffusivitySolutEffCO)
      RCHDFS=(CCH4GQ-CCH4S2)*AMIN1(VLWatMicP,DiffusivitySolutEffCH)
      ROXDFS=(COXYGQ-COXYS2)*AMIN1(VLWatMicP,DiffusivitySolutEffOX)
      RNGDFS=(CZ2GGQ-CZ2GS2)*AMIN1(VLWatMicP,DiffusivitySolutEffNG)
      RN2DFS=(CZ2OGQ-CZ2OS2)*AMIN1(VLWatMicP,DiffusivitySolutEffN2)
      RN3DFS=(CZN3GQ-CNH3S2)*AMIN1(VLWatMicP*forc%VLNH4,DiffusivitySolutEffN3)
      RNBDFS=(CZN3BQ-CNH3B2)*AMIN1(VLWatMicP*forc%VLNHB,DiffusivitySolutEffN3)
      RHGDFS=(CH2GGQ-CH2GS2)*AMIN1(VLWatMicP,DiffusivitySolutEffHL)

!     accumulate atmospheric dissolution/volatilization
      XCODFS=XCODFS+RCODFS
      XCHDFS=XCHDFS+RCHDFS
      XOXDFS=XOXDFS+ROXDFS
      XNGDFS=XNGDFS+RNGDFS
      XN2DFS=XN2DFS+RN2DFS
      XN3DFS=XN3DFS+RN3DFS
      XNBDFS=XNBDFS+RNBDFS
      XHGDFS=XHGDFS+RHGDFS

      RCODXS=RCODFS*dt_GasCyc
      RCHDXS=RCHDFS*dt_GasCyc
      ROXDXS=ROXDFS*dt_GasCyc
      RNGDXS=RNGDFS*dt_GasCyc
      RN2DXS=RN2DFS*dt_GasCyc
      RN3DXS=RN3DFS*dt_GasCyc
      RNBDXS=RNBDFS*dt_GasCyc
      RHGDXS=RHGDFS*dt_GasCyc

  ! dissolution between gaseous and aqueous phases in soil
      VLWatMicPMXA=natomw*VLWatMicPMA
      VLWatMicPMXB=natomw*VLWatMicPMB
      VLsoiAirPMA=VLsoiAirP*forc%VLNH4
      VLsoiAirPMB=VLsoiAirP*forc%VLNHB
      VLWatMicPMCO=VLWatMicP*forc%SCO2L
      VLWatMicPMCH=VLWatMicP*forc%CH4AquaSolubility
      VLWatMicPMOX=VLWatMicP*forc%O2GSolubility
      VLWatMicPMNG=VLWatMicP*forc%SN2GL
      VLWatMicPMN2=VLWatMicP*forc%SN2OL
      VLWatMicPMN3=VLWatMicPMA*forc%SNH3L
      VLWatMicPMNB=VLWatMicPMB*forc%SNH3L
      VLWatMicPMHG=VLWatMicP*forc%SH2GL
      VOLCOT=VLWatMicPMCO+VLsoiAirP
      VOLCHT=VLWatMicPMCH+VLsoiAirP
      VOLOXT=VLWatMicPMOX+VLsoiAirP
      VOLNGT=VLWatMicPMNG+VLsoiAirP
      VOLN2T=VLWatMicPMN2+VLsoiAirP
      VOLN3T=VLWatMicPMN3+VLsoiAirPMA
      VOLBranchNumber_pft=VLWatMicPMNB+VLsoiAirPMB
      VOLHGT=VLWatMicPMHG+VLsoiAirP

      RCODFG=forc%DiffusivitySolutEff*(AMAX1(ZEROS,CO2G2)*VLWatMicPMCO-AMAX1(ZEROS,CO2S2+RCODXS)*VLsoiAirP)/VOLCOT
      RCHDFG=forc%DiffusivitySolutEff*(AMAX1(ZEROS,CH4G2)*VLWatMicPMCH-AMAX1(ZEROS,CH4S2+RCHDXS)*VLsoiAirP)/VOLCHT
      ROXDFG=forc%DiffusivitySolutEff*(AMAX1(ZEROS,OXYG2)*VLWatMicPMOX-AMAX1(ZEROS,OXYS2+ROXDXS)*VLsoiAirP)/VOLOXT
      RNGDFG=forc%DiffusivitySolutEff*(AMAX1(ZEROS,Z2GG2)*VLWatMicPMNG-AMAX1(ZEROS,Z2GS2+RNGDXS)*VLsoiAirP)/VOLNGT
      RN2DFG=forc%DiffusivitySolutEff*(AMAX1(ZEROS,Z2OG2)*VLWatMicPMN2-AMAX1(ZEROS,Z2OS2+RN2DXS)*VLsoiAirP)/VOLN2T
      RHGDFG=forc%DiffusivitySolutEff*(AMAX1(ZEROS,H2GG2)*VLWatMicPMHG-AMAX1(ZEROS,H2GS2+RHGDXS)*VLsoiAirP)/VOLHGT
      XCODFG=XCODFG+RCODFG
      XCHDFG=XCHDFG+RCHDFG
      XOXDFG=XOXDFG+ROXDFG
      XNGDFG=XNGDFG+RNGDFG
      XN2DFG=XN2DFG+RN2DFG
      XN3DFG=XN3DFG+RN3DFG
      XNBDFG=XNBDFG+RNBDFG
      XHGDFG=XHGDFG+RHGDFG

      IF(VOLN3T.GT.ZEROS2.AND.VLWatMicPMXA.GT.ZEROS2)THEN
        RN3DFG=forc%DiffusivitySolutEff*(AMAX1(ZEROS,ZN3G2)*VLWatMicPMN3-AMAX1(ZEROS,ZNH3S2+RN3DXS)*VLsoiAirPMA)/VOLN3T
        CNH3S0=AZMAX1((ZNH3S2+RN3DFG)/VLWatMicPMXA)
        CNH4S0=AZMAX1(ZNH4S2)/VLWatMicPMXA
      ELSE
        RN3DFG=0.0_r8
      ENDIF
      IF(VOLBranchNumber_pft.GT.ZEROS2.AND.VLWatMicPMXB.GT.ZEROS2)THEN
        RNBDFG=forc%DiffusivitySolutEff*(AMAX1(ZEROS,ZN3G2)*VLWatMicPMNB-AMAX1(ZEROS,ZNH3B2+RNBDXS)*VLsoiAirPMB)/VOLBranchNumber_pft
        CNH3B0=AZMAX1((ZNH3B2+RNBDFG)/VLWatMicPMXB)
        CNH4B0=AZMAX1(ZNH4B2)/VLWatMicPMXB
      ELSE
        RNBDFG=0.0_r8
      ENDIF

      CO2S2=CO2S2+RCODFS
      CH4S2=CH4S2+RCHDFS
      OXYS2=OXYS2+ROXDFS
      Z2GS2=Z2GS2+RNGDFS
      Z2OS2=Z2OS2+RN2DFS
      ZNH3S2=ZNH3S2+RN3DFS
      ZNH3B2=ZNH3B2+RNBDFS
      H2GS2=H2GS2+RHGDFS

      CO2S2=CO2S2+RCODFG
      CH4S2=CH4S2+RCHDFG
      OXYS2=OXYS2+ROXDFG
      Z2GS2=Z2GS2+RNGDFG
      Z2OS2=Z2OS2+RN2DFG
      ZNH3S2=ZNH3S2+RN3DFG
      ZNH3B2=ZNH3B2+RNBDFG
      H2GS2=H2GS2+RHGDFG

      CO2G2=CO2G2-RCODFG+DFVCOG
      CH4G2=CH4G2-RCHDFG+DFVCHG
      OXYG2=OXYG2-ROXDFG+DFVOXG
      Z2GG2=Z2GG2-RNGDFG+DFVNGG
      Z2OG2=Z2OG2-RN2DFG+DFVN2G
      ZN3G2=ZN3G2-RN3DFG-RNBDFG+DFVN3G
      H2GG2=H2GG2-RHGDFG+DFVHGG
    endif
  enddo

  ystatesfl(fid_XCODFS)=XCODFS
  ystatesfl(fid_XCHDFS)=XCHDFS
  ystatesfl(fid_XOXDFS)=XOXDFS
  ystatesfl(fid_XNGDFS)=XNGDFS
  ystatesfl(fid_XN2DFS)=XN2DFS
  ystatesfl(fid_XN3DFS)=XN3DFS
  ystatesfl(fid_XNBDFS)=XNBDFS
  ystatesfl(fid_XHGDFS)=XHGDFS

  ystatesfl(fid_XCODFG)=XCODFG
  ystatesfl(fid_XCHDFG)=XCHDFG
  ystatesfl(fid_XOXDFG)=XOXDFG
  ystatesfl(fid_XNGDFG)=XNGDFG
  ystatesfl(fid_XN2DFG)=XN2DFG
  ystatesfl(fid_XN3DFG)=XN3DFG
  ystatesfl(fid_XNBDFG)=XNBDFG
  ystatesfl(fid_XHGDFG)=XHGDFG

  ystatesfl(fid_XCOFLG)=XCOFLG
  ystatesfl(fid_XCHFLG)=XCHFLG
  ystatesfl(fid_XOXFLG)=XOXFLG
  ystatesfl(fid_XNGFLG)=XNGFLG
  ystatesfl(fid_XN2FLG)=XN2FLG
  ystatesfl(fid_XN3FLG)=XN3FLG
  ystatesfl(fid_XHGFLG)=XHGFLG



!  CO2S(NU_col(NY,NX),NY,NX)=CO2S(NU_col(NY,NX),NY,NX)+XCODFS(NY,NX)
!  CH4S(NU_col(NY,NX),NY,NX)=CH4S(NU_col(NY,NX),NY,NX)+XCHDFS(NY,NX)
!  OXYS(NU_col(NY,NX),NY,NX)=OXYS(NU_col(NY,NX),NY,NX)+XOXDFS(NY,NX)
!  Z2GS(NU_col(NY,NX),NY,NX)=Z2GS(NU_col(NY,NX),NY,NX)+XNGDFS(NY,NX)
!  Z2OS(NU_col(NY,NX),NY,NX)=Z2OS(NU_col(NY,NX),NY,NX)+XN2DFS(NY,NX)
!  ZNH3S(NU_col(NY,NX),NY,NX)=ZNH3S(NU_col(NY,NX),NY,NX)+XN3DFS(NY,NX)
!  ZNH3B(NU_col(NY,NX),NY,NX)=ZNH3B(NU_col(NY,NX),NY,NX)+XNBDFS(NY,NX)
!  H2GS(NU_col(NY,NX),NY,NX)=H2GS(NU_col(NY,NX),NY,NX)+XHGDFS(NY,NX)


!    CO2S(L,NY,NX)=CO2S(L,NY,NX)+XCODFG(L,NY,NX)-trcs_RMicbTransf_vr(idg_CO2,L,NY,NX)
!    CH4S(L,NY,NX)=CH4S(L,NY,NX)+XCHDFG(L,NY,NX)-trcs_RMicbTransf_vr(idg_CH4,L,NY,NX)
!    OXYS(L,NY,NX)=OXYS(L,NY,NX)+XOXDFG(L,NY,NX)-trcs_RMicbTransf_vr(idg_O2,L,NY,NX
!    Z2GS(L,NY,NX)=Z2GS(L,NY,NX)+XNGDFG(L,NY,NX)-trcs_RMicbTransf_vr(idg_N2,L,NY,NX)-Micb_N2Fixation_vr(L,NY,NX)
!    Z2OS(L,NY,NX)=Z2OS(L,NY,NX)+XN2DFG(L,NY,NX)-trcs_RMicbTransf_vr(idg_N2O,L,NY,NX)
!    H2GS(L,NY,NX)=H2GS(L,NY,NX)+XHGDFG(L,NY,NX)-trcs_RMicbTransf_vr(idg_H2,L,NY,NX)
!    ZNH3S(L,NY,NX)=ZNH3S(L,NY,NX)+XN3DFG(L,NY,NX)+TRProd_chem_sol_NH3_soil_vr(L,NY,NX)

!    CO2G(L,NY,NX)=CO2G(L,NY,NX)-XCODFG(L,NY,NX)+XCOFLG
!    CH4G(L,NY,NX)=CH4G(L,NY,NX)-XCHDFG(L,NY,NX)+XCHFLG
!    OXYG(L,NY,NX)=OXYG(L,NY,NX)-XOXDFG(L,NY,NX)+XOXFLG
!    Z2GG(L,NY,NX)=Z2GG(L,NY,NX)-XNGDFG(L,NY,NX)+XNGFLG
!    Z2OG(L,NY,NX)=Z2OG(L,NY,NX)-XN2DFG(L,NY,NX)+XN2FLG
!    ZNH3G(L,NY,NX)=ZNH3G(L,NY,NX)-XN3DFG(L,NY,NX)-XNBDFG(L,NY,NX)+XN3FLG
!    H2GG(L,NY,NX)=H2GG(L,NY,NX)-XHGDFG(L,NY,NX)+XHGFLG


  end associate
  end subroutine CalcSurflux
! ----------------------------------------------------------------------

  subroutine UpdateSOMORGM(micfor,micstt)
  implicit none
  type(micforctype), intent(inout)    :: micfor
  type(micsttype)  , intent(inout) :: micstt

  real(r8) :: DC,DN,DP        !litter
  real(r8) :: OC,ON,OP        !SOM
  real(r8) :: ORGC,ORGN,ORGR
  integer :: K,N,M,NGL,MID

  DC=0.0_r8
  DN=0.0_r8
  DP=0.0_r8
  OC=0.0_r8
  ON=0.0_r8
  OP=0.0_r8

  DO K=1,micpar%jcplx
    IF(micpar%is_litter(K))THEN
      DO N=1,micpar%NumMicbHFunGrupsPerCmplx
        DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
          DO  M=1,micpar%nlbiomcp
            MID=micpar%get_micb_id(M,NGL)
            DC=DC+micstt%mBiomeHeter(ielmc,MID,K)
            DN=DN+micstt%mBiomeHeter(ielmn,MID,K)
            DP=DP+micstt%mBiomeHeter(ielmp,MID,K)
          ENDDO
        enddo
      ENDDO
    ELSE
      DO N=1,micpar%NumMicbHFunGrupsPerCmplx
        DO NGL=micpar%JGniH(N),micpar%JGnfH(N)
          DO  M=1,micpar%nlbiomcp
            MID=micpar%get_micb_id(M,NGL)          
            OC=OC+micstt%mBiomeHeter(ielmc,MID,K)
            ON=ON+micstt%mBiomeHeter(ielmn,MID,K)
            OP=OP+micstt%mBiomeHeter(ielmp,MID,K)
          enddo
        enddo
      ENDDO
    ENDIF
  ENDDO
! abstract complex
  DO  N=1,micpar%NumMicbAFunGrupsPerCmplx
    DO NGL=micpar%JGniA(N),micpar%JGnfA(N)
      DO  M=1,micpar%nlbiomcp
        MID=micpar%get_micb_id(M,NGL)               
        OC=OC+micstt%mBiomeAutor(ielmc,MID)
        ON=ON+micstt%mBiomeAutor(ielmn,MID)
        OP=OP+micstt%mBiomeAutor(ielmp,MID)
      enddo
    enddo
  ENDDO
! microbial residue
  DO K=1,micpar%jcplx
    IF(micpar%is_litter(K))THEN
      DO M=1,micpar%ndbiomcp
        DC=DC+micstt%OMBioResdu(ielmc,M,K)
        DN=DN+micstt%OMBioResdu(ielmn,M,K)
        DP=DP+micstt%OMBioResdu(ielmp,M,K)
      ENDDO
!solutes in macropores are not subject to microbial attack
      DC=DC+micstt%DOM(idom_doc,K)+micstt%SorbedOM(ielmc,K)+micstt%DOM(idom_acetate,K)+micstt%SorbedOM(idom_acetate,K)  !+micstt%DOM_Macp(idom_doc,K)+micstt%DOM_Macp(idom_acetate,K)
      DN=DN+micstt%DOM(idom_don,K)+micstt%SorbedOM(ielmn,K) !+micstt%DOM_Macp(idom_don,K)
      DP=DP+micstt%DOM(idom_dop,K)+micstt%SorbedOM(ielmp,K) !+micstt%DOM_Macp(idom_dop,K)
      DO M=1,micpar%jsken
        DC=DC+micstt%SolidOM(ielmc,M,K)
        DN=DN+micstt%SolidOM(ielmn,M,K)
        DP=DP+micstt%SolidOM(ielmp,M,K)
      ENDDO
    ELSE
      DO M=1,micpar%ndbiomcp
        OC=OC+micstt%OMBioResdu(ielmc,M,K)
        ON=ON+micstt%OMBioResdu(ielmn,M,K)
        OP=OP+micstt%OMBioResdu(ielmp,M,K)
      ENDDO
      OC=OC+micstt%DOM(idom_doc,K)+micstt%SorbedOM(ielmc,K)+micstt%DOM(idom_acetate,K)+micstt%SorbedOM(idom_acetate,K)  !micstt%DOM_Macp(idom_acetate,K)++micstt%DOM_Macp(idom_doc,K)
      ON=ON+micstt%DOM(idom_don,K)+micstt%SorbedOM(ielmn,K)  !+micstt%DOM_Macp(idom_don,K)
      OP=OP+micstt%DOM(idom_dop,K)+micstt%SorbedOM(ielmp,K)  !+micstt%DOM_Macp(idom_dop,K)
      DO M=1,micpar%jsken
        OC=OC+micstt%SolidOM(ielmc,M,K)
        ON=ON+micstt%SolidOM(ielmn,M,K)
        OP=OP+micstt%SolidOM(ielmp,M,K)
      ENDDO
    ENDIF
  ENDDO
! DC is for litter complex, and OC is for POM and humus complex
  micfor%ORGC=DC+OC
  ORGN=DN+ON
  ORGR=DC
  end subroutine UpdateSOMORGM

end module batchmod
