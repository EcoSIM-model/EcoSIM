module SoilPhysParaMod
  ! Description:
  ! Evaluate soil/pond pressure and thermal properties. Hydraulic callers share
  ! the current retention curve and named auxiliary water-content calculations.
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only: micpar
  use DebugToolMod,     only: PrintInfo,DebugPrint
  use abortutils, only : endrun
  USE SoilWaterDataType
  use SoilPropertyDataType
  USE SoilPhysDataType
  use EcoSIMCtrlDataType
  use HydroThermData
  use PhysPars
  use EcoSimConst
  use FlagDataType
  USE AqueChemDatatype
  USE SOMDataType
  use SoilBGCDataType
  use MiniMathMod
  use SoilHeatDataType
  USE GridDataType
implicit none
  private
  character(len=*), parameter :: mod_filename=&
  __FILE__
  PUBLIC :: CalcSoilWatPotential
  public :: SetDeepSoil
  public :: CalcSoilThermConductivity
  public :: get_Tfrez
  public :: ComputePsiMCM  
  public :: ComputeSoilAirDryContent
  public :: ComputeKsatReferenceContent
  public :: ComputePSIPond
  public :: getMoistK
  public :: IsPondLayer

  contains
!------------------------------------------------------------------------------------------
  subroutine CalcSoilWatPotential(NY,NX,N1,N2,N3,PSISoilMatric,THETA1S)
  !
  !Description
  !Compute soil matric potential
  implicit none
  integer, intent(in) :: NY,NX,N1,N2,N3
  real(r8), intent(out) :: PSISoilMatric      ! soil matric potential
  real(r8), optional, intent(out) :: THETA1S  ! liquid-water filled porosity [m3 water m-3 micropore soil]
  character(len=*), parameter :: subname='CalcSoilWatPotential'

  real(r8) :: THETA1                          !volumetric soil moisture
  !     BKVL=soil mass
  !     FC,WP=water contents at field capacity,wilting point
  !     FCL,WPL=log FC,WP
  !     FCD,PSD=FCL-WPL,log(POROS)-FCL
  !     PSISoilMatric,PSIHY,PSISE=soil matric,hygroscopic,air entry potential
  !     PSIMX,PSIMD,LOGPSIAtSat=log water potential at FC,WP,saturation
  !     PSISD=PSIMX-LOGPSIAtSat
  !     SRP=parameter for deviation from linear log-log water retention
  !     PSISO=osmotic potential
  !     DTHETW=minimum water content for numerical purpose
  ! soil matric potential upper layer
  call PrintInfo('beg '//subname)

  THETA1=AMAX1(SoilWatAirDry_vr(N3,N2,N1),AMIN1(POROS_vr(N3,N2,N1),safe_adb(VLWatMicP1_vr(N3,N2,N1),VLSoilMicP_vr(N3,N2,N1))))

  IF(VLSoilMicPMass_vr(N3,N2,N1).GT.ZEROS(NY,NX))THEN
    !source layer is active soil
    call ComputePsiMCM(N3,N2,N1,THETA1,PSISoilMatric)    
    !
    !     SUBSURFCE UPPER WATER LAYER
    !
    !     FracSoiPAsIce,FracSoiPAsWat=ice,water concentration
    !     FCI,WPI=ice field capacity,wilting point
    !     PSISoilMatric=matric water potential
    !
  ELSEIF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    call ComputePSIPond(N3,N2,N1,FracSoiPAsIce_vr(N3,N2,N1),FracSoiPAsWat_vr(N3,N2,N1),&
      THETA1,PSISoilMatric,TKSoil1_vr(N3,N2,N1))
  ELSE
    THETA1        = POROS_vr(N3,N2,N1)
    PSISoilMatric = PSI_at_Sautration_vr(N3,N2,N1)
  ENDIF
  if(present(THETA1S))THETA1S=THETA1

  call PrintInfo('end '//subname)  
  end subroutine CalcSoilWatPotential
!------------------------------------------------------------------------------------------

  subroutine ComputePSIPond(L,NY,NX,ThetafI,ThetafW,ThetaW,PSI,TKS)
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: ThetafI
  real(r8), intent(inout) :: ThetafW
  real(r8), intent(inout) :: ThetaW
  real(r8), intent(out) :: PSI
  real(r8), optional, intent(in) :: TKS   !temperature [K]

  real(r8) :: FCX,WPX,FCLX,WPLX,PSDX,FCDX
  character(len=*), parameter :: subname='ComputePSIPond'

  call PrintInfo('beg '//subname)
  FCX  = FCI*ThetafI
  WPX  = WPI*ThetafI

  IF(ThetafW.LT.FCX)THEN
    !less than field capacity, but nonzero water
    FCLX = LOG(FCX)
    WPLX = LOG(WPX)    
    if(ThetafW.GT.ZEROS2(NY,NX))then
      FCDX = FCLX-WPLX
      PSI=AMAX1(PSIHY,-EXP(LOGPSIFLD_col(NY,NX)+(FCLX-LOG(ThetafW))*LOGPSIMND_col(NY,NX)/FCDX))
    else
      if(.not.present(TKS))then
        call endrun('TKS is not input to '//trim(subname)//' in '//trim(mod_filename)//' at line',__LINE__)        
      endif  
      !pure ice
      PSI=LtHeatIceMelt*(TKS-TFice)/TFice
    endif
  ELSE IF(ThetafW.LT.POROS_vr(L,NY,NX)-DTHETW .and. FCX > 0._r8 .and. WPX >0._r8)THEN  
    !more than field capacity
    FCLX = LOG(FCX)
    WPLX = LOG(WPX)
    PSDX = LOGPOROS_vr(L,NY,NX)-FCLX
    FCDX = FCLX-WPLX
    PSI=-EXP(LOGPSIAtSat(NY,NX)+(LOGPOROS_vr(L,NY,NX)-LOG(ThetafW))*LOGPSIMXD_col(NY,NX)/PSDX)
  ELSE
    ThetaW = POROS_vr(L,NY,NX)
    PSI    = PSI_at_Sautration_vr(L,NY,NX)
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine ComputePSIPond
!------------------------------------------------------------------------------------------

  subroutine SetDeepSoil(NY,NX,NM,JZ)

  implicit none
  integer, intent(in) :: NY,NX,NM,JZ
  integer :: L
  associate(                            &
  k_woody_comp => micpar%k_woody_comp , &
  k_fine_comp  => micpar%k_fine_comp  , &
  k_manure     => micpar%k_manure       &
  )

  DO L=NM+1,JZ
    CumDepz2LayBottom_vr(L,NY,NX) = 2.0_r8*CumDepz2LayBottom_vr(L-1,NY,NX)-1.0_r8*CumDepz2LayBottom_vr(L-2,NY,NX)
    SoiBulkDensityt0_vr(L,NY,NX)  = SoiBulkDensityt0_vr(L-1,NY,NX)
    FieldCapacity_vr(L,NY,NX)     = FieldCapacity_vr(L-1,NY,NX)
    WiltPoint_vr(L,NY,NX)         = WiltPoint_vr(L-1,NY,NX)
    SatHydroCondVert_vr(L,NY,NX)  = SatHydroCondVert_vr(L-1,NY,NX)
    SatHydroCondHrzn_vr(L,NY,NX)  = SatHydroCondHrzn_vr(L-1,NY,NX)
    CSAND_vr(L,NY,NX)             = CSAND_vr(L-1,NY,NX)
    CSILT_vr(L,NY,NX)             = CSILT_vr(L-1,NY,NX)
    CCLAY_vr(L,NY,NX)             = CCLAY_vr(L-1,NY,NX)
    SoilFracAsMacP_vr(L,NY,NX)    = SoilFracAsMacP_vr(L-1,NY,NX)
    ROCK_vr(L,NY,NX)              = ROCK_vr(L-1,NY,NX)
    PH_vr(L,NY,NX)                = PH_vr(L-1,NY,NX)
    CEC_vr(L,NY,NX)               = CEC_vr(L-1,NY,NX)
    AEC_vr(L,NY,NX)               = AEC_vr(L-1,NY,NX)

!     IF(IDWaterTable_col(NY,NX).EQ.0)THEN
!       0.25_r8 is the geometric decreasing ratio (tunable)
!  or different scheme can be used
    CSoilOrgM_vr(ielmc,L,NY,NX) = 0.25_r8*CSoilOrgM_vr(ielmc,L-1,NY,NX)
    COMLitrC_vr(L,NY,NX)        = 0.25_r8*COMLitrC_vr(L-1,NY,NX)
    CSoilOrgM_vr(ielmn,L,NY,NX) = 0.25_r8*CSoilOrgM_vr(ielmn,L-1,NY,NX)
    CSoilOrgM_vr(ielmp,L,NY,NX) = 0.25_r8*CSoilOrgM_vr(ielmp,L-1,NY,NX)
!     ELSE
!       CSoilOrgM_vr(ielmc,L,NY,NX)=CSoilOrgM_vr(ielmc,L-1,NY,NX)
!       COMLitrC_vr(L,NY,NX)=COMLitrC_vr(L-1,NY,NX)
!       CSoilOrgM_vr(ielmn,L,NY,NX)=CSoilOrgM_vr(ielmn,L-1,NY,NX)
!       CSoilOrgM_vr(ielmp,L,NY,NX)=CSoilOrgM_vr(ielmp,L-1,NY,NX)
!     ENDIF
    CNH4_vr(L,NY,NX)             = CNH4_vr(L-1,NY,NX)
    CNO3_vr(L,NY,NX)             = CNO3_vr(L-1,NY,NX)
    CPO4_vr(L,NY,NX)             = CPO4_vr(L-1,NY,NX)
    CAL_vr(L,NY,NX)              = CAL_vr(L-1,NY,NX)
    CFE_vr(L,NY,NX)              = CFE_vr(L-1,NY,NX)
    CCA_vr(L,NY,NX)              = CCA_vr(L-1,NY,NX)
    CMG_vr(L,NY,NX)              = CMG_vr(L-1,NY,NX)
    CNA_vr(L,NY,NX)              = CNA_vr(L-1,NY,NX)
    CKA_vr(L,NY,NX)              = CKA_vr(L-1,NY,NX)
    CSO4_vr(L,NY,NX)             = CSO4_vr(L-1,NY,NX)
    CCL_vr(L,NY,NX)              = CCL_vr(L-1,NY,NX)
    CALOH_vr(L,NY,NX)            = CALOH_vr(L-1,NY,NX)
    CFEOH_vr(L,NY,NX)            = CFEOH_vr(L-1,NY,NX)
    CCACO_vr(L,NY,NX)            = CCACO_vr(L-1,NY,NX)
    CCASO_vr(L,NY,NX)            = CCASO_vr(L-1,NY,NX)
    CALPO_vr(L,NY,NX)            = CALPO_vr(L-1,NY,NX)
    CFEPO_vr(L,NY,NX)            = CFEPO_vr(L-1,NY,NX)
    CCAPD_vr(L,NY,NX)            = CCAPD_vr(L-1,NY,NX)
    CCAPH_vr(L,NY,NX)            = CCAPH_vr(L-1,NY,NX)
    GKC4_vr(L,NY,NX)             = GKC4_vr(L-1,NY,NX)
    GKCH_vr(L,NY,NX)             = GKCH_vr(L-1,NY,NX)
    GKCA_vr(L,NY,NX)             = GKCA_vr(L-1,NY,NX)
    GKCM_vr(L,NY,NX)             = GKCM_vr(L-1,NY,NX)
    GKCN_vr(L,NY,NX)             = GKCN_vr(L-1,NY,NX)
    GKCK_vr(L,NY,NX)             = GKCK_vr(L-1,NY,NX)
    THW_vr(L,NY,NX)           = THW_vr(L-1,NY,NX)
    THI_vr(L,NY,NX)           = THI_vr(L-1,NY,NX)
    ISOIL_vr(1:4,L,NY,NX)        = ISOIL_vr(1:4,L-1,NY,NX)
    RSC_vr(k_fine_comp,L,NY,NX)  = 0.0_r8
    RSN_vr(k_fine_comp,L,NY,NX)  = 0.0_r8
    RSP_vr(k_fine_comp,L,NY,NX)  = 0.0_r8
    RSC_vr(k_woody_comp,L,NY,NX) = 0.0_r8
    RSN_vr(k_woody_comp,L,NY,NX) = 0.0_r8
    RSP_vr(k_woody_comp,L,NY,NX) = 0.0_r8
    RSC_vr(k_manure,L,NY,NX)     = 0.0_r8
    RSN_vr(k_manure,L,NY,NX)     = 0.0_r8
    RSP_vr(k_manure,L,NY,NX)     = 0.0_r8
  ENDDO
  END associate
  end subroutine SetDeepSoil

!------------------------------------------------------------------------------------------
  subroutine CalcSoilThermConductivity(N1,N2,N3,DTKX,TCND1)
  !
  !Description
  !Compute soil thermal conductivity
  !
  implicit none
  integer , intent(in) :: N1,N2,N3
  real(r8), intent(in) :: DTKX     !absolute temeprature gradient
  real(r8), intent(out):: TCND1

  real(r8) :: HeatDiffusByWat1,HeatDiffusByAir1,RYLXW1,RYLXA1,RYLNW1,RYLNA1
  REAL(R8) :: XNUSW1,XNUSA1,ThermalConducByWater,ThermalConducByAir,WTHET1

  IF(SoilBulkDensity_vr(N3,N2,N1).GT.ZERO .OR. FracSoiPAsWat_vr(N3,N2,N1)+FracSoiPAsIce_vr(N3,N2,N1).GT.ZERO)THEN
    !it is a soil layer or pure water layer
    HeatDiffusByWat1     = AZMAX1(FracSoiPAsWat_vr(N3,N2,N1)-TRBW)**3
    HeatDiffusByAir1     = AZMAX1(FracAirFilledSoilPore_vr(N3,N2,N1)-TRBA)**3
    RYLXW1               = DTKX*HeatDiffusByWat1*1.e-6_r8  !
    RYLXA1               = DTKX*HeatDiffusByAir1*1.e-6_r8  !convert into energy scale MJ
    RYLNW1               = AMIN1(1.0E+04_r8,RYLXW*RYLXW1)
    RYLNA1               = AMIN1(1.0E+04_r8,RYLXA*RYLXA1)
    XNUSW1               = AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNW1**0.25_r8/DNUSW)
    XNUSA1               = AMAX1(1.0_r8,0.68_r8+0.67_r8*RYLNA1**0.25_r8/DNUSA)
    ThermalConducByWater = 2.067E-03_r8*XNUSW1
    ThermalConducByAir   = 9.050E-05_r8*XNUSA1
    WTHET1               = 1.467_r8-0.467_r8*FracSoilAsAirt(N3,N2,N1)
    TCND1                = (NumerSolidThermCond_vr(N3,N2,N1)+FracSoiPAsWat_vr(N3,N2,N1)*ThermalConducByWater &
      +0.611_r8*FracSoiPAsIce_vr(N3,N2,N1)*7.844E-03_r8+WTHET1*FracAirFilledSoilPore_vr(N3,N2,N1)*ThermalConducByAir) &
      /(DenomSolidThermCond_vr(N3,N2,N1)+FracSoiPAsWat_vr(N3,N2,N1)+0.611_r8*FracSoiPAsIce_vr(N3,N2,N1) &
      +WTHET1*FracAirFilledSoilPore_vr(N3,N2,N1))
  ELSE
    TCND1=0.0_r8
  ENDIF
  end subroutine CalcSoilThermConductivity
!------------------------------------------------------------------------------------------

  pure function get_Tfrez(PSI)result(TFREEZ)
  !
  !computing freezing temeprature after considering depression.
  !apply the Clausis-Clapeyron equation for depressed freezing temperature, 
  !Cary and Mayland (1972), 9.0959E+04_r8~=273.15*LtHeatIceMelt  

  implicit none
  REAL(R8), intent(in) :: PSI   !matric pressure, [MPa]
  real(r8) :: TFREEZ

  !-5.5_r8 is used for numerical stability
  TFREEZ  = -9.0959E+04_r8/(AMAX1(PSI,-5.5_r8)-LtHeatIceMelt)

  END function get_Tfrez  

!------------------------------------------------------------------------------------------

  subroutine ComputePsiMCM(N6,N5,N4,THETWL,PSISM1)
  ! Description:
  !Compute soil matric pressure
  implicit none
  integer, intent(in) :: N6,N5,N4
  real(r8), intent(inout) :: THETWL  !liquid content [m3 m-3 soil]; set to porosity near saturation
  real(r8), intent(out) :: PSISM1    !Soil matric pressure [MPa]
  character(len=*), parameter :: subname='ComputePsiMCM'

  call PrintInfo('beg '//subname)
  IF(THETWL.LT.FieldCapacity_vr(N6,N5,N4))THEN
    !less than field capacity
    PSISM1=AMAX1(PSIHY,-EXP(LOGPSIFLD_col(N5,N4)+((LOGFldCapacity_vr(N6,N5,N4)-LOG(THETWL)) &
      /FCD_vr(N6,N5,N4)*LOGPSIMND_col(N5,N4))))
  ELSEIF(THETWL.LT.POROS_vr(N6,N5,N4)-DTHETW)THEN
    PSISM1=-EXP(LOGPSIAtSat(N5,N4)+(((LOGPOROS_vr(N6,N5,N4)-LOG(THETWL)) &
      /PSD_vr(N6,N5,N4))**SRP_vr(N6,N5,N4)*LOGPSIMXD_col(N5,N4)))
  ELSE
    !saturated
    THETWL = POROS_vr(N6,N5,N4)
    PSISM1 = PSI_at_Sautration_vr(N6,N5,N4)
  ENDIF
  call PrintInfo('end '//subname)
  END subroutine ComputePsiMCM

!------------------------------------------------------------------------------------------
  subroutine ComputeSoilAirDryContent(L,NY,NX,WaterContent)
  ! Description:
  ! Evaluate the current dry-branch water content at PSIHY. Callers retain
  ! responsibility for distinguishing ordinary soil, litter and pond layers.
  implicit none
  integer, intent(in) :: L       ! vertical layer index, including litter at zero [-]
  integer, intent(in) :: NY,NX   ! row and column indices [-]
  real(r8), intent(out) :: WaterContent ! dry-branch liquid content [m3 m-3 soil]
  character(len=*), parameter :: subname='ComputeSoilAirDryContent'

  call PrintInfo('beg '//subname)
  WaterContent=EXP((LOGPSIFLD_col(NY,NX)-LOG(-PSIHY))*FCD_vr(L,NY,NX)/LOGPSIMND_col(NY,NX)+LOGFldCapacity_vr(L,NY,NX))
  call PrintInfo('end '//subname)
  end subroutine ComputeSoilAirDryContent
!------------------------------------------------------------------------------------------
  function ComputeKsatReferenceContent(L,NY,NX) result(WaterContent)
  ! Description:
  ! Return the logarithmic -0.033 MPa reference content used to estimate Ksat.
  ! This historical expression omits SRP and is not a general retention inverse.
  implicit none
  integer, intent(in) :: L       ! ordinary-soil layer index [-]
  integer, intent(in) :: NY,NX   ! row and column indices [-]
  real(r8) :: WaterContent ! returned reference liquid content [m3 m-3 soil]
  character(len=*), parameter :: subname='ComputeKsatReferenceContent'

  call PrintInfo('beg '//subname)
  WaterContent=AMIN1(POROS_vr(L,NY,NX),EXP((LOGPSIAtSat(NY,NX)-LOG(0.033_r8)) &
    *(LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX))/LOGPSIMXD_col(NY,NX)+LOGPOROS_vr(L,NY,NX)))
  call PrintInfo('end '//subname)
  end function ComputeKsatReferenceContent
!------------------------------------------------------------------------------------------
  pure function getMoistK(thetawl,pores)result(K)
  implicit none
  real(r8), intent(in) :: thetawl  !liquid water [m3 H2O d-2]
  real(r8), intent(in) :: pores    !pore space [m3 pores d-2]
  integer :: K

  K=MAX(1,MIN(100,101-INT(100.0_r8*thetawl/pores)))
  end function getMoistK

!------------------------------------------------------------------------------------------
  pure function IsPondLayer(bulkDS)result(ans)
  implicit none
  real(r8), intent(in) :: bulkDS
  logical :: ans

  ans =bulkDS.LE.ZERO

  end function IsPondLayer

end  module SoilPhysParaMod
