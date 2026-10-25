module SoilHydroParaMod

  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only: micpar
  use minimathmod,      only: isclose, AZMAX1, AZMIN1, Viscosity_H2O
  Use SoilPhysParaMod,  only: ComputePsiMCM,ComputePSIPond,getMoistK
  use DebugToolMod
  use EcoSimConst
  use SoilWaterDataType
  use SoilPropertyDataType
  use SoilHeatDataType
  use GridDataType
  use SoilPhysDataType
  use AqueChemDatatype
  use EcoSIMCtrlDataType
  use SOMDataType
  use FlagDataType
  use LandSurfDataType
  use EcoSIMHistMod
  use EcoSIMConfig
  use SurfLitterDataType
  use EcoSIMCtrlMod
  use PhysPars

implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  real(r8), parameter :: FORGW=0.25E+06_r8   !threshold for  C concentration in organic soil 	gC Mg-1
  real(r8), parameter :: tiny_val=1.e-6_r8   !minimum value 
  ! ECOSIM_ZM: controlled legacy-ecosys hydraulic policy. These switches
  ! change constitutive properties, not water/heat accounting or boundary fluxes.
  ! Set false to recover upstream behavior. Not a universal soil calibration.
  public :: GetSoilHydraulicVars
  public :: SoilHydroProperty
  public :: LitterHydroproperty
  PUBLIC :: ComputeSoilHydroPars

  public :: SetColdRunSoilStates
contains
!------------------------------------------------------------------------------------------

  subroutine GetSoilHydraulicVars(I,J,NY,NX)
  !
  !DESCRIPTIONS
  !compute hydraulic properties
  !called in hour1.F90 before doing hydrology
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  character(len=*), parameter :: subname='GetSoilHydraulicVars'

  real(r8) :: THETW1           !water-filled soil micropores ![m3 water /(m3 micropore soil)]
  real(r8) :: WPX,WPLX,CC,DD,EE,CCLAY_pct,CORGCT
  integer :: K,L

  ! begin_execution
  call PrintInfo('beg '//subname)

  DO L=NUI_col(NY,NX),NLI_col(NY,NX)
    ! WATER POTENTIALS
    !
    ! FC,WP=water contents at field capacity,wilting point,saturation
    ! PSISM,PSISE=matric,saturation water potential
    ! SRP=parameter for deviation from linear log-log water retention
    ! FC,WP=water contents at field capacity,wilting point from soil file
    ! FCL,WPL=log FC,WP
    ! FCD,PSD=FCL-WPL,log(POROS)-FCL
    ! FCI,WPI=FC,WP of ice
    ! THETIX=ice concentration

    !significant soil mass
    IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX) .AND. VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
      THETW1=AZMAX1(AMIN1(POROS_vr(L,NY,NX),VLWatMicP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)),tiny_val)
      
      call ComputePsiMCM(L,NY,NX,THETW1,PSISoilMatricP_vr(L,NY,NX))

      !not significant soil mass, but have significant ice, for (partially) frozen ponding water 
    ELSE IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
      call ComputePSIPond(L,NY,NX,THETI_vr(L,NY,NX),THETW_vr(L,NY,NX),THETW1,PSISoilMatricP_vr(L,NY,NX),TKS_vr(L,NY,NX))

    ELSE
      PSISoilMatricP_vr(L,NY,NX)=PSI_at_Sautration_vr(L,NY,NX)
    ENDIF
    !
    !     SOIL OSMOTIC, GRAVIMETRIC AND MATRIC WATER POTENTIALS
    !
    !     PSISM,PSISO,PSIGrav_vr,PSIST=matric,osmotic,gravimetric,total water potential
    !
    PSISoilOsmotic_vr(L,NY,NX)          = -RGASC*TKS_vr(L,NY,NX)*SolutesIonConc_vr(L,NY,NX)*1.E-6_r8  !MPa
    PSIGrav_vr(L,NY,NX)                 = mGravAccelerat*(SurfElevation_col(NY,NX)-SoilDepthMidLay_vr(L,NY,NX))
    ElvAdjstedSoilH2OPSIMPa_vr(L,NY,NX) = AZMIN1(PSISoilMatricP_vr(L,NY,NX)+PSISoilOsmotic_vr(L,NY,NX)+PSIGrav_vr(L,NY,NX))
    !
    !     SoilBulkModulus4RootPent_vr=bulk modulus of the undisturbed soil [MPa]
    !
    IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      !Eq.(8) in Grant (1993), Simulation model of soil compaction and root growth I. Model structure.
      !fitted using data from Rickman et al. (1992)
      if(THETW_vr(L,NY,NX)>0.01_r8)THEN
        !a better algorithm may be developed.
        CCLAY_pct = CCLAY_vr(L,NY,NX)*1.0E+02_r8
        CORGCT    = CSoilOrgM_vr(ielmc,L,NY,NX)*1.0E-04_r8
        CC        = EXP(-3.6733_r8-0.1447_r8*CCLAY_pct+0.7653_r8*CORGCT)
        DD        = -0.4805_r8-0.1239_r8*CCLAY_pct+0.2080_r8*CORGCT
        EE        = 3.8521_r8+0.0963_r8*CCLAY_pct
        SoilBulkModulus4RootPent_vr(L,NY,NX) = CC*THETW_vr(L,NY,NX)**DD*SoilBulkDensity_vr(L,NY,NX)**EE
      ELSE
        SoilBulkModulus4RootPent_vr(L,NY,NX)=  300._r8
      ENDIF
      SoilWeightStress_vr(L,NY,NX)=SoilBulkDensity_vr(L,NY,NX)*mGravAccelerat*DLYR_3D(3,L,NY,NX) ![Mg m-3]*[km s-2]*[m]=[MPa]
    ELSE
       SoilBulkModulus4RootPent_vr(L,NY,NX) = 0.0_r8
       SoilWeightStress_vr(L,NY,NX)          = 0._r8
    ENDIF
    !
    !     SOIL HYDRAULIC CONDUCTIVITIES FROM AMBIENT SOIL WATER CONTENTS
    !
    !     HYCDMicP4RootUptake_vr=soil hydraulic conductivity for root uptake
    !
    K=getMoistK(THETW_vr(L,NY,NX),POROS_vr(L,NY,NX))
    HYCDMicP4RootUptake_vr(L,NY,NX)=0.5_r8*(HydroCond_3D(1,K,L,NY,NX)+HydroCond_3D(3,K,L,NY,NX))

  END DO
  call PrintInfo('end '//subname)
  end subroutine GetSoilHydraulicVars


!------------------------------------------------------------------------------------------

  subroutine SoilHydroProperty(I,J,L,NY,NX)
  !
  !Set up soil hydraulic property
  implicit none
  integer, intent(in) :: L,NY,NX
  integer, intent(in) :: I,J
  real(r8) :: THETF
  real(r8) :: H2OSOIK(100),PSISK(0:100)
  character(len=*), parameter :: subname='SoilHydroProperty'

  integer :: K,M,N,N100
  real(r8) :: XK,YK,SUM1,SUM2
  real(r8) :: VISCWL

  call PrintInfo('beg '//subname)
  !gC/Mg soil
  IF(CSoilOrgM_vr(ielmc,L,NY,NX).GT.FORGC)THEN
    !organic soil
    SRP_vr(L,NY,NX)=0.25_r8
  ELSE IF(CSoilOrgM_vr(ielmc,L,NY,NX).GT.0.5_r8*FORGC)THEN
    SRP_vr(L,NY,NX)=0.33_r8
  ELSE
    SRP_vr(L,NY,NX)=1.00_r8
  ENDIF
  ! ECOSIM_ZM: hour1.f uses SRP=0.67 for all mineral/organic soil branches.
  ! Rebuild K(theta) below from this same pressure curve; do not rescale fluxes.
  if (use_ecosys_drainage) SRP_vr(L,NY,NX)=0.67_r8
  
  ! double check is_cold_run() setup
  LOGPOROS_vr(L,NY,NX)=LOG(POROS_vr(L,NY,NX))

  ! Refresh estimated constitutive limits in every run mode. This changes
  ! neither water/ice volumes nor heat capacity; restricting it to cold runs
  ! freezes FC/WP after restart while continuous runs keep updating them.
  call EstimateSoilWaterLimits(L,NY,NX)
  LOGFldCapacity_vr(L,NY,NX) = LOG(FieldCapacity_vr(L,NY,NX))
  LOGWiltPoint_vr(L,NY,NX)   = LOG(WiltPoint_vr(L,NY,NX))
  PSD_vr(L,NY,NX)            = LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX)
  FCD_vr(L,NY,NX)            = LOGFldCapacity_vr(L,NY,NX)-LOGWiltPoint_vr(L,NY,NX)

  VLsoiAirP_vr(L,NY,NX)=AZMAX1(VLMicP_vr(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLiceMicP_vr(L,NY,NX)) &
    +AZMAX1(VLMacroPores_vr(L,NY,NX)-VLWatMacP_vr(L,NY,NX)-VLiceMacP_vr(L,NY,NX))

  IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
    !ratio of total air-filled pore to micropore
    ThetaAir_vr(L,NY,NX)=VLsoiAirP_vr(L,NY,NX)/VLSoilMicP_vr(L,NY,NX)
  ELSE
    ThetaAir_vr(L,NY,NX)=0.0_r8
  ENDIF

  IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
    !PSIHY: soil matric potential at dessication
    SoilWatAirDry_vr(L,NY,NX)=EXP((LOGPSIFLD_col(NY,NX)-LOG(-PSIHY))*FCD_vr(L,NY,NX)/LOGPSIMND_col(NY,NX)+LOGFldCapacity_vr(L,NY,NX))
  ELSE
    SoilWatAirDry_vr(L,NY,NX)=ZERO2
  ENDIF
      !
      !     SATURATED HYDRAULIC CONDUCTIVITY FROM SWC AT SATURATION VS.
      !     -0.033 MPA (MINERAL SOILS) IF NOT ENTERED IN SOIL FILE IN 'READS'
      !
      !     SCNV,SCNH=vertical,lateral saturated hydraulic conductivity
  !
  IF(ISOIL_vr(isoi_scnv,L,NY,NX).EQ.isoi_unset)THEN
    !computing vertical saturated hydraulic conductivity
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS_vr(L,NY,NX),EXP((LOGPSIAtSat(NY,NX)-LOG(0.033_r8)) &
        *(LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX))/LOGPSIMXD_col(NY,NX)+LOGPOROS_vr(L,NY,NX)))
      SatHydroCondVert_vr(L,NY,NX)=1.54_r8*((POROS_vr(L,NY,NX)-THETF)/THETF)**2
    ELSE
      SatHydroCondVert_vr(L,NY,NX) = 0.10_r8+75.0_r8*1.0E-15_r8**SoilBulkDensity_vr(L,NY,NX)
      SatHydroCondVert_vr(L,NY,NX) = SatHydroCondVert_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
    ENDIF
  ENDIF

  IF(ISOIL_vr(isoi_scnh,L,NY,NX).EQ.isoi_unset)THEN
    !computing horizontal saturated hydraulic conductivity
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      THETF=AMIN1(POROS_vr(L,NY,NX),EXP((LOGPSIAtSat(NY,NX)-LOG(0.033_r8)) &
        *(LOGPOROS_vr(L,NY,NX)-LOGFldCapacity_vr(L,NY,NX))/LOGPSIMXD_col(NY,NX)+LOGPOROS_vr(L,NY,NX)))
      SatHydroCondHrzn_vr(L,NY,NX)=1.54_r8*((POROS_vr(L,NY,NX)-THETF)/THETF)**2
    ELSE
      SatHydroCondHrzn_vr(L,NY,NX) = 0.10_r8+75.0_r8*1.0E-15_r8**SoilBulkDensity_vr(L,NY,NX)
      SatHydroCondHrzn_vr(L,NY,NX) = SatHydroCondHrzn_vr(L,NY,NX)*FracSoiAsMicP_vr(L,NY,NX)
    ENDIF
  ENDIF

  !
  !     HYDRAULIC CONDUCTIVITY FUNCTION FROM KSAT AND SOIL WATER RELEASE CURVE
  !
  !     THETK,PSISK=micropore class water content,potential
  !     HCND=lateral(1,2),vertical(3) micropore hydraulic conductivity
  !
  !     IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
  !using Eq.(1) and II in Table 1 from Green and Corey, 1971
  SUM2=0.0_r8;N100=100
  !integrate of 1/Psi(x) over [0,1] 
  DO  K=1,n100
    !air-filled relative pore volume
    XK  = (K-1._r8)/n100
    !moisture content
    H2OSOIK(K) = POROS_vr(L,NY,NX)*(1._r8-XK)
    call ComputePsiMCM(L,NY,NX,H2OSOIK(K),PSISK(K))
    SUM2=SUM2+(2._r8*K-1)/(PSISK(K)**2)
  ENDDO

  !progress from dry to saturation
  !completely dry at XK==1
  DO  K=1,N100
    
    !air-filled relative pore volume
    XK   = (K-1._r8)/N100
    !saturation at K
    YK   = (1._r8-XK)**1.33_r8
    SUM1 = 0.0_r8
    DO M    = K, N100
      SUM1=SUM1+(2*(M-K)+1._r8)/(PSISK(M)**2)
    ENDDO

    DO  N=1,3
      IF(N.EQ.3)THEN
        !vertical
        HydroCond_3D(N,K,L,NY,NX)=SatHydroCondVert_vr(L,NY,NX)*YK*SUM1/SUM2

        !obtain water content at saturation
        if (use_ecosys_drainage) then
          ! ECOSIM_ZM: legacy FSCNV=0.1 crossing, using preceding wetter
          ! entry for BOTH pressure and theta. Internal saturated/Green-Ampt
          ! switches and natural-outlet wetness now share the same definition.
          if (K.gt.1) then
            if (HydroCond_3D(3,K,L,NY,NX).lt.0.1_r8*SatHydroCondVert_vr(L,NY,NX) .and. &
                HydroCond_3D(3,K-1,L,NY,NX).ge.0.1_r8*SatHydroCondVert_vr(L,NY,NX)) then
              PSISoilAirEntry_vr(L,NY,NX)=PSISK(K-1)
              ThetaSat_vr(L,NY,NX)=H2OSOIK(K-1)
            endif
          endif
        else
          IF(K.GT.1 .AND. PSISK(K).LT.PSISoilAirEntry_vr(L,NY,NX))then
            if(PSISK(K-1).GE.PSISoilAirEntry_vr(L,NY,NX))THEN
              ThetaSat_vr(L,NY,NX)=H2OSOIK(K)
            endif
          ENDIF
        endif
      ELSE
        !horizontal
        HydroCond_3D(N,K,L,NY,NX)=SatHydroCondHrzn_vr(L,NY,NX)*YK*SUM1/SUM2
      ENDIF
    ENDDO
  ENDDO

  !     SOIL MACROPORE DIMENSIONS AND CONDUCTIVITY FROM MACROPORE FRACTION
  !     ENTERED IN 'READS'
  !
  !     PathLenMacP,MacPNumLayer,MacPRadius=path length between, number,radius of macropores
  !     CNDH=macropore hydraulic conductivity
  !VGeomLayert0_vr: geometric volume of the layer 
  !VLMacroPores_vr: volume of macropores in the layer
  !Assuming soil macropores are approximated as (cylydrical) pipes, and has layer thickness as pipe length
  !Eq.(3.3)-(3.4) from Dimitrov et al. (2010)
  
  MacPoreRadius_vr(L,NY,NX)   = 0.5E-03_r8
  MacPoreNumbers_vr(L,NY,NX) = INT(VLMacroPores_vr(L,NY,NX)/(PICON*MacPoreRadius_vr(L,NY,NX)**2*VGeomLayert0_vr(L,NY,NX)))

  IF(MacPoreNumbers_vr(L,NY,NX).GT.0.0_r8)THEN
    PathLenMacPore_vr(L,NY,NX)=1.0_r8/(SQRT(PICON*MacPoreNumbers_vr(L,NY,NX)))
  ELSE
    PathLenMacPore_vr(L,NY,NX)=1.0_r8
  ENDIF
  !water viscosity as a function of temperature
  VISCWL  = Viscosity_H2O(TCS_vr(L,NY,NX))
  !Apply the Poiseuille's Law
  HydroCondMacP_vr(L,NY,NX) = 3.6E+03_r8*PICON*MacPoreNumbers_vr(L,NY,NX)*MacPoreRadius_vr(L,NY,NX)**4/(8.0_r8*VISCWL)

  call PrintInfo('end '//subname)
  end subroutine SoilHydroProperty

!------------------------------------------------------------------------------------------
  subroutine EstimateSoilWaterLimits(L,NY,NX)
  ! Estimate both FC and WP when either was unspecified; retain the input flags.
  implicit none
  integer, intent(in) :: L,NY,NX
  character(len=*), parameter :: subname='EstimateSoilWaterLimits'

  call PrintInfo('beg '//subname)
  IF(ISOIL_vr(isoi_fc,L,NY,NX).EQ.isoi_unset .OR. ISOIL_vr(isoi_wp,L,NY,NX).EQ.isoi_unset)THEN
    !when either field capacity or wilting point is not provided, both are recomputed 
    !for consistency.
    !calculating FC or WP
    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      FieldCapacity_vr(L,NY,NX)=0.2576_r8-0.20_r8*CSAND_vr(L,NY,NX) &
            +0.36_r8*CCLAY_vr(L,NY,NX)+0.60E-06*CSoilOrgM_vr(ielmc,L,NY,NX)
    ELSE
      !organic soil
      IF(SoilBulkDensity_vr(L,NY,NX).LT.0.075_r8)THEN
        FieldCapacity_vr(L,NY,NX)=0.27_r8
      ELSEIF(SoilBulkDensity_vr(L,NY,NX).LT.0.195_r8)THEN
        FieldCapacity_vr(L,NY,NX)=0.62_r8
      ELSE
        FieldCapacity_vr(L,NY,NX)=0.71_r8
      ENDIF
    ENDIF
    FieldCapacity_vr(L,NY,NX) = FieldCapacity_vr(L,NY,NX)/(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
    FieldCapacity_vr(L,NY,NX) = AMIN1(0.75_r8*POROS_vr(L,NY,NX),FieldCapacity_vr(L,NY,NX))

    IF(CSoilOrgM_vr(ielmc,L,NY,NX).LT.FORGW)THEN
      WiltPoint_vr(L,NY,NX)=0.0260_r8+0.50_r8*CCLAY_vr(L,NY,NX)+0.32E-06_r8*CSoilOrgM_vr(ielmc,L,NY,NX)
    ELSE
      !organic soil
      IF(SoilBulkDensity_vr(L,NY,NX).LT.0.075_r8)THEN
        WiltPoint_vr(L,NY,NX)=0.04_r8
      ELSEIF(SoilBulkDensity_vr(L,NY,NX).LT.0.195_r8)THEN
        WiltPoint_vr(L,NY,NX)=0.15_r8
      ELSE
        WiltPoint_vr(L,NY,NX)=0.22_r8
      ENDIF
    ENDIF
    WiltPoint_vr(L,NY,NX) = WiltPoint_vr(L,NY,NX)/(1.0_r8-SoilFracAsMacP_vr(L,NY,NX))
    WiltPoint_vr(L,NY,NX) = AMIN1(0.75_r8*FieldCapacity_vr(L,NY,NX),WiltPoint_vr(L,NY,NX))
  ENDIF
  
  call PrintInfo('end '//subname)
  end subroutine EstimateSoilWaterLimits
!------------------------------------------------------------------------------------------

  subroutine SetColdRunSoilStates(L,NY,NX)
  ! Initialize volumes before chemistry. Supplied FC/WP retain the prescribed
  ! THW/THI interpretation; estimated limits retain the water-table saturation rule.
  implicit none
  integer, intent(in) :: L,NY,NX
  logical :: SaturateBelowWaterTable
  real(r8) :: WaterRelSaturation,IceRelSaturation
  character(len=*), parameter :: subname='SetColdRunSoilStates'

  call PrintInfo('beg '//subname)
  SaturateBelowWaterTable = (ISOIL_vr(isoi_fc,L,NY,NX).EQ.isoi_unset .OR. &
    ISOIL_vr(isoi_wp,L,NY,NX).EQ.isoi_unset) .AND. &
    SoilDepthMidLay_vr(L,NY,NX).GE.ExtWaterTableLocDepz_col(NY,NX)
  call EstimateSoilWaterLimits(L,NY,NX)

  IF(THW_vr(L,NY,NX).GT.1.0_r8 .OR. SaturateBelowWaterTable)THEN
    !below the water table, thus it is saturated
    THETW_vr(L,NY,NX)=POROS_vr(L,NY,NX)
  ELSEIF(isclose(THW_vr(L,NY,NX),1._r8))THEN
    !at field capacity
    THETW_vr(L,NY,NX)=FieldCapacity_vr(L,NY,NX)
  ELSEIF(isclose(THW_vr(L,NY,NX),0._r8))THEN
    !at wilting point
    THETW_vr(L,NY,NX)=WiltPoint_vr(L,NY,NX)
  ELSEIF(THW_vr(L,NY,NX).LT.0.0_r8)THEN
    !completely dry
    THETW_vr(L,NY,NX)=0.0_r8
  ELSE
    THETW_vr(L,NY,NX)=THW_vr(L,NY,NX)
  ENDIF

  IF(THI_vr(L,NY,NX).GT.1.0_r8.OR.SaturateBelowWaterTable)THEN
    THETI_vr(L,NY,NX)=AZMAX1(AMIN1(POROS_vr(L,NY,NX),POROS_vr(L,NY,NX)-THETW_vr(L,NY,NX)))
  ELSEIF(isclose(THI_vr(L,NY,NX),1._r8))THEN
    THETI_vr(L,NY,NX)=AZMAX1(AMIN1(FieldCapacity_vr(L,NY,NX),POROS_vr(L,NY,NX)-THETW_vr(L,NY,NX)))
  ELSEIF(isclose(THI_vr(L,NY,NX),0._r8))THEN
    THETI_vr(L,NY,NX)=AZMAX1(AMIN1(WiltPoint_vr(L,NY,NX),POROS_vr(L,NY,NX)-THETW_vr(L,NY,NX)))
  ELSEIF(THI_vr(L,NY,NX).LT.0.0_r8)THEN
    THETI_vr(L,NY,NX)=0.0_r8
  ELSE
    THETI_vr(L,NY,NX)=THI_vr(L,NY,NX)
  ENDIF

  ! THETW/THETI are volumes per unit micropore-domain soil volume. Macropore
  ! volume is already pore space, so initialize it with relative saturation.
  ! Give water priority and limit ice to the remaining macropore space.
  WaterRelSaturation=0._r8
  IceRelSaturation=0._r8
  IF(POROS_vr(L,NY,NX).GT.ZERO)THEN
    WaterRelSaturation=AZMAX1(AMIN1(1._r8,THETW_vr(L,NY,NX)/POROS_vr(L,NY,NX)))
    IceRelSaturation=AZMAX1(AMIN1(1._r8-WaterRelSaturation,THETI_vr(L,NY,NX)/POROS_vr(L,NY,NX)))
  ENDIF

  ! Initialize water/ice volumes and their contribution to heat capacity.
  VLWatMicP_vr(L,NY,NX)     = THETW_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
  VLWatMicPX_vr(L,NY,NX)    = VLWatMicP_vr(L,NY,NX)
  VLWatMacP_vr(L,NY,NX)     = WaterRelSaturation*VLMacroPores_vr(L,NY,NX)
  VLiceMicP_vr(L,NY,NX)     = THETI_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)
  VLiceMacP_vr(L,NY,NX)     = IceRelSaturation*VLMacroPores_vr(L,NY,NX)
  ! Match the relative pore occupancies reported after hourly redistribution.
  VOLTX_vr(L,NY,NX)=POROS_vr(L,NY,NX)*VLSoilPoreMicP_vr(L,NY,NX)+VLMacroPores_vr(L,NY,NX)
  IF(VOLTX_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
    ThetaH2OZ_vr(L,NY,NX)=AZMAX1(AMIN1(1._r8, &
      (VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX))/VOLTX_vr(L,NY,NX)))
    ThetaICEZ_vr(L,NY,NX)=AZMAX1(AMIN1(1._r8, &
      (VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))/VOLTX_vr(L,NY,NX)))
  ELSE
    ThetaH2OZ_vr(L,NY,NX)=0._r8
    ThetaICEZ_vr(L,NY,NX)=0._r8
  ENDIF
  VHeatCapacity_vr(L,NY,NX) = VHeatCapSolidSoil_vr(L,NY,NX)+Cpw*(VLWatMicP_vr(L,NY,NX) &
    +VLWatMacP_vr(L,NY,NX))+Cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))
  VLsoiAirP_vr(L,NY,NX)=AZMAX1(VLMicP_vr(L,NY,NX)-VLWatMicP_vr(L,NY,NX)-VLiceMicP_vr(L,NY,NX)) &
    +AZMAX1(VLMacroPores_vr(L,NY,NX)-VLWatMacP_vr(L,NY,NX)-VLiceMacP_vr(L,NY,NX))
  call PrintInfo('end '//subname)
  end subroutine SetColdRunSoilStates
!------------------------------------------------------------------------------------------
  subroutine LitterHydroproperty(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  integer  :: K,M
  real(r8) :: XK,YK
  real(r8) :: SUM1,SUM2
  integer, parameter :: n100=100
  real(r8) :: H2OSOIK(n100),PSISK(0:n100)

  IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX))THEN
    SoilBulkDensity_vr(0,NY,NX)=VLSoilMicPMass_vr(0,NY,NX)/VGeomLayer_vr(0,NY,NX)
  ELSE
    SoilBulkDensity_vr(0,NY,NX)=BulkDensLitR(micpar%k_fine_comp)
  ENDIF
  SoilWatAirDry_vr(0,NY,NX)=EXP((LOGPSIFLD_col(NY,NX)-LOG(-PSIHY))*FCD_vr(0,NY,NX)/LOGPSIMND_col(NY,NX)+LOGFldCapacity_vr(0,NY,NX))
  SUM2             =0.0_r8
  D1220: DO  K=1,n100
    XK         = (K-1._r8)/n100
    H2OSOIK(K) = POROS0_col(NY,NX)*(1._r8-XK)
    call ComputePsiMCM(0,NY,NX,H2OSOIK(K),PSISK(K))
    SUM2=SUM2+(2*K-1)/(PSISK(K)**2)
  ENDDO D1220

  D1235: DO  K=1,n100
    SUM1  = 0.0_r8
    XK    = (K-1._r8)/N100
    !air-filled pore
    YK    = (1._r8-XK)**1.33_r8
    D1230: DO M = K, n100
      SUM1=SUM1+(2*M+1-2*K)/(PSISK(M)**2)
    ENDDO D1230
    HydroCond_3D(3,K,0,NY,NX) = SatHydroCondVert_vr(0,NY,NX)*YK*SUM1/SUM2
    HydroCond_3D(1,K,0,NY,NX) = 0.0_r8
    HydroCond_3D(2,K,0,NY,NX) = 0.0_r8
    if (use_ecosys_drainage) then
      ! ECOSIM_ZM: apply the same conductivity-defined threshold to litter,
      ! matching legacy hour1.f. Preserve its existing retention curvature.
      if (K.gt.1) then
        if (HydroCond_3D(3,K,0,NY,NX).lt.0.1_r8*SatHydroCondVert_vr(0,NY,NX) .and. &
            HydroCond_3D(3,K-1,0,NY,NX).ge.0.1_r8*SatHydroCondVert_vr(0,NY,NX)) then
          PSISoilAirEntry_vr(0,NY,NX)=PSISK(K-1)
          ThetaSat_vr(0,NY,NX)=H2OSOIK(K-1)
        endif
      endif
    else
      IF(K.GT.1.AND.(PSISK(K).LT.PSISoilAirEntry_vr(0,NY,NX)))THEN
        IF(PSISK(K-1).GE.PSISoilAirEntry_vr(0,NY,NX))THEN
          ThetaSat_vr(0,NY,NX)=H2OSOIK(K)
        ENDIF
      ENDIF
    endif
  ENDDO D1235
  end subroutine LitterHydroproperty

!------------------------------------------------------------------------------------------

  subroutine ComputeSoilHydroPars(NY,NX,NU,NM)
  implicit none
  integer, intent(in) :: NY,NX,NU,NM
  integer :: L
  DO  L=NU,NM
    IF(FieldCapacity_vr(L,NY,NX).LT.0.0_r8)THEN
      !computing field capacity
      ISOIL_vr(isoi_fc,L,NY,NX)  = isoi_unset
      PSIAtFldCapacity_col(NY,NX) = -0.033_r8
    ELSE
      ISOIL_vr(isoi_fc,L,NY,NX)=isoi_set
    ENDIF

    IF(WiltPoint_vr(L,NY,NX).LT.0.0_r8)THEN
      !computing wilting point
      ISOIL_vr(isoi_wp,L,NY,NX) = isoi_unset
      PSIAtWiltPoint_col(NY,NX)  = -1.5_r8
    ELSE
      ISOIL_vr(isoi_wp,L,NY,NX)=isoi_set
    ENDIF
    IF(SatHydroCondVert_vr(L,NY,NX).LT.0.0_r8)THEN
      !soil vertical saturated hydraulic conductivity
      ISOIL_vr(isoi_scnv,L,NY,NX)=isoi_unset
    ELSE
      ISOIL_vr(isoi_scnv,L,NY,NX)=isoi_set
    ENDIF
    IF(SatHydroCondHrzn_vr(L,NY,NX).LT.0.0_r8)THEN
      !soil horizontal saturated hydraulic conductivity
      ISOIL_vr(isoi_scnh,L,NY,NX)=isoi_unset
    ELSE
      ISOIL_vr(isoi_scnh,L,NY,NX)=isoi_set
    ENDIF
  ENDDO
  END subroutine ComputeSoilHydroPars
!------------------------------------------------------------------------------------------

  subroutine calculate_soil_properties(sand_in, silt_in, clay_in, d10, d30, d50, d60, cu, cc)
  ! -------------------------------------------------------------------------
  ! Subroutine: calculate_soil_properties
  ! Purpose: Estimates D10, D30, D60, Cu, and Cc using piecewise log-linear 
  !          interpolation anchored to USDA boundaries.
  ! -------------------------------------------------------------------------
  implicit none
  real(r8), intent(in)  :: sand_in, silt_in, clay_in !texture fraction of mineral soil
  real(r8), intent(out) :: d10, d30, d50, d60 !particle size at 10%, 30%, 50% and 60% cdf [mm]
  real(r8), intent(out) :: cu, cc
  
  real(r8) :: diameters(4), percent_passing(4)
  real(r8) :: sand, clay, silt

  sand = sand_in * 100.0_r8
  silt = silt_in * 100.0_r8
  clay = clay_in * 100.0_r8
  

  ! 2. Define Boundaries (USDA Standards) & Cumulative Percent Passing
  ! Array index 1: Top (2.0mm) -> 100% Passing
  ! Array index 2: Sand/Silt Boundary (0.05mm) -> Silt + Clay Passing
  ! Array index 3: Silt/Clay Boundary (0.002mm) -> Clay Passing
  ! Array index 4: Minimum assumption (0.0001mm) -> 0% Passing
  
  diameters(1) = 2.0_r8
  percent_passing(1) = 100.0_r8
  
  diameters(2) = 0.05_r8
  percent_passing(2) = silt + clay
  
  diameters(3) = 0.002_r8
  percent_passing(3) = clay
  
  diameters(4) = 0.0001_r8
  percent_passing(4) = 0.0_r8

  ! 3. Interpolate D-Values
  d10 = get_diameter_at_percent(10.0_r8, diameters, percent_passing)
  d30 = get_diameter_at_percent(30.0_r8, diameters, percent_passing)
  d50 = get_diameter_at_percent(50.0_r8, diameters, percent_passing)
  d60 = get_diameter_at_percent(60.0_r8, diameters, percent_passing)

  ! 4. Calculate Coefficients (Cu, Cc)
  if (d10 > 1.0e-7_r8) then
    cu = d60 / d10
    cc = (d30**2) / (d10 * d60)
  else
    cu = 0.0_r8
    cc = 0.0_r8
  end if

  end subroutine calculate_soil_properties
!------------------------------------------------------------------------------------------

  function get_diameter_at_percent(target_pct, d_arr, p_arr) result(d_val)
  !
  !Description
  ! Function: get_diameter_at_percent
  ! Purpose: Helper function to perform the log-linear interpolation
  ! -------------------------------------------------------------------------
  implicit none
  real(r8), intent(in) :: target_pct
  real(r8), intent(in) :: d_arr(4), p_arr(4)
  real(r8) :: d_val
  
  integer :: i
  real(r8) :: d_upper, p_upper, d_lower, p_lower, ratio, log_d

  ! Handle Edge Cases
  if (target_pct >= 100.0_r8) then
    d_val = 2.0_r8
    return
  end if
  if (target_pct <= 0.0_r8) then
    d_val = 0.0001_r8
    return
  end if

  ! Loop through intervals to find where the target percent falls
  d_val = 0.0001_r8 ! Default fallback
  
  do i = 1, 3
    d_upper = d_arr(i)
    p_upper = p_arr(i)
    
    d_lower = d_arr(i+1)
    p_lower = p_arr(i+1)
    
    ! Check if target is in this interval (inclusive)
    ! Note: Arrays are descending in size/percent (100 -> 0)
    if (target_pct <= p_upper .and. target_pct >= p_lower) then        
      if (abs(p_upper - p_lower) < 1.0e-6_r8) then
        ! Flat shelf (avoid divide by zero)
        d_val = d_lower
      else
        ! Log-Linear Interpolation Formula
        ! ratio = (target - lower%) / (upper% - lower%)
        ratio = (target_pct - p_lower) / (p_upper - p_lower)
        
        ! log(D) = log(d_lower) + ratio * (log(d_upper) - log(d_lower))
        log_d = log10(d_lower) + ratio * (log10(d_upper) - log10(d_lower))
        
        d_val = 10.0_r8**log_d
      end if
      
      return ! Exit once found
    end if
  end do
        
  end function get_diameter_at_percent


! -------------------------------------------------------------------------
  subroutine estimate_natural_soil_limits(sand_in, silt_in, clay_in, om_in, D50_mm, e_min, e_max, status)
  !Description:    
  ! Estimates e_max (loosest) and e_min (densest) void ratios for 
  !          natural soils based on texture and organic matter.
  ! Inputs
  real(r8), intent(in)  :: sand_in   ! fraction as Sand (0-1)
  real(r8), intent(in)  :: silt_in   ! fraction as Silt (0-1)
  real(r8), intent(in)  :: clay_in   ! fraction as Clay (0-1)
  real(r8), intent(in)  :: om_in     ! fraction as Organic Matter (0-1)
  real(r8), intent(in)  :: D50_mm    ! 50% particle size, [mm]
  ! Outputs
  real(r8), intent(out) :: e_min     ! Estimated Minimum Void Ratio
  real(r8), intent(out) :: e_max     ! Estimated Maximum Void Ratio
  integer, intent(out)  :: status    ! 0 = Success, -1 = Error
  
  ! Local Variables
  real(r8) :: total_min, f_fines, ratio
  real(r8) :: e_min_base, e_max_base
  real(r8) :: om_factor_min, om_factor_max
  real(r8) :: opt_decimal,safe_D50

  ! Constants for Packing Limits (V-Shape Curve)
  ! Pure Sand (Uniform)
  real(r8), parameter :: ES_MIN_SAND = 0.55_r8
  real(r8), parameter :: ES_MAX_SAND = 0.85_r8
  ! Optimal Mix (Well Graded / Silty Sand ~25-30% fines)
  real(r8), parameter :: ES_MIN_OPT  = 0.35_r8
  real(r8), parameter :: ES_MAX_OPT  = 0.65_r8
  ! Pure Fines (Clay/Silt Matrix)
  real(r8), parameter :: ES_MIN_CLAY = 0.70_r8
  real(r8), parameter :: ES_MAX_CLAY = 1.30_r8
      
  ! 1. Validate and Normalize Inputs
  total_min = sand_in + silt_in + clay_in
  
  if (total_min < 0.1e-2_r8) then
    status = -1
    e_min = 0.0_r8
    e_max = 0.0_r8
    return
  end if

  safe_D50 = max(D50_mm, 0.001_r8)

  ! Calculate fraction of fines (Silt + Clay)
  f_fines = silt_in + clay_in 
  
  ! 2. Determine Base Mineral Limits (V-Shape Curve)
  ! Zone 1: Sand dominated (Fines fill voids) -> Void ratio decreases
  ! Boundary assumed at 30% fines

  ! Your Formula: Base 20% + scaling with size
  opt_decimal = 0.20_r8 + 0.15_r8 * log10(safe_D50)
    
  !SAFETY CLAMPS:
  ! 1. Min: Fine sands can usually hold at least 5-10% fines before expanding.
  ! 2. Max: Even coarse gravels usually lose grain-to-grain contact >40% fines.
    
  opt_decimal = max(0.05_r8, min(0.40_r8, opt_decimal))

  if (f_fines <= opt_decimal) then
    ratio = f_fines / opt_decimal
    
    ! Linear Interpolation between Pure Sand and Optimal Mix
    e_min_base = ES_MIN_SAND - ratio * (ES_MIN_SAND - ES_MIN_OPT)
    e_max_base = ES_MAX_SAND - ratio * (ES_MAX_SAND - ES_MAX_OPT)      
  else
    ! Zone 2: Fines dominated (Sand floats) -> Void ratio increases
    ! Interpolate between Optimal Mix and Pure Fines
    ratio = (f_fines - opt_decimal) / (1._r8-opt_decimal)
    
    e_min_base = ES_MIN_OPT + ratio * (ES_MIN_CLAY - ES_MIN_OPT)
    e_max_base = ES_MAX_OPT + ratio * (ES_MAX_CLAY - ES_MAX_OPT)
  end if
  
  ! 3. Apply Organic Matter (OM) Correction
  ! OM adds "fluffiness".
  ! Heuristic: e_min increases by ~0.03 per 1% OM
  !            e_max increases by ~0.05 per 1% OM
  
  om_factor_min = 3._r8 * om_in
  om_factor_max = 5._r8 * om_in
  
  e_min = e_min_base + om_factor_min
  e_max = e_max_base + om_factor_max
  
  status = 0
      
  end subroutine estimate_natural_soil_limits


! -------------------------------------------------------------------------
  subroutine estimate_friction_natural(e_in, e_min_in, e_max_in, &
    D50_mm_in, sigma_v_total_kPa_in, rel_sat_in, suction_kPa_in, &
    clay_in, phi_peak, phi_cv, dilatancy, p_prime, root_strain_pct, status)
  !Description:        
  ! Purpose: Estimates Peak Friction Angle using Stress-Dilatancy Theory 
  !          adapted for natural (unsaturated/aggregated) soils.
  ! -------------------------------------------------------------------------

  ! --- Inputs ---
  real(r8), intent(in)  :: e_in                 ! Current Void Ratio
  real(r8), intent(in)  :: e_min_in, e_max_in   ! Packing Limits
  real(r8), intent(in)  :: D50_mm_in            ! Median Grain Size [mm]
  real(r8), intent(in)  :: sigma_v_total_kPa_in ! Total Vertical Overburden [kPa]
  real(r8), intent(in)  :: clay_in              ! clay fraction [0-1]
  real(r8), intent(in)  :: rel_sat_in           ! Relative Saturation (0.0 - 1.0)
  real(r8), intent(in)  :: suction_kPa_in       ! Matric Suction Magnitude [kPa]
        
  ! --- Outputs ---
  real(r8), intent(out) :: phi_peak         ! Peak Friction Angle [deg]
  real(r8), intent(out) :: phi_cv           ! Constant Volume (Base) Friction [deg]
  real(r8), intent(out) :: dilatancy        ! Dilatancy Component [deg]
  real(r8), intent(out) :: p_prime          ! Mean Effective Stress [kPa]
  real(r8), intent(out) :: root_strain_pct  !Strain at Peak Failure of the soil [%]
  integer, intent(out)  :: status           ! 0 = Success, -1 = Error
        
  ! --- Local Variables ---
  real(r8) :: Dr                  ! Relative Density (0-1)
  real(r8) :: Sr_clamped          ! Saturation clamped
  real(r8) :: sigma_suction_eff   ! Stress boost from suction
  real(r8) :: sigma_v_eff         ! Total Vertical Effective Stress
  real(r8) :: sigma_h_eff         ! Horizontal Effective Stress
  real(r8) :: K0                  ! Coeff of Earth Pressure at Rest
  real(r8) :: reduction           ! Phi_cv reduction for organics
  real(r8) :: Q_natural           ! Crushing Parameter
  real(r8) :: suppression_term    ! Stress suppression factor
  
  ! Constants
  real(r8), parameter :: K0_DEFAULT = 0.5_r8
  real(r8), parameter :: K0_wetland = 0.65_r8
  real(r8), parameter :: PHI_CV_QUARTZ = 30.0_r8
  real(r8), parameter :: Q_NATURAL_VAL = 8.0_r8
  
  ! 1. Validate Inputs
  if (e_max_in <= e_min_in) then
      status = -1
      phi_peak = 0.0_r8
      return
  end if
  
  ! 2. Calculate Relative Density (Dr)
  if (e_max_in <= e_min_in) then
    Dr = 0.5_r8 ! Fallback
  else
    Dr = (e_max_in - e_in) / (e_max_in - e_min_in)
    Dr = max(0.0_r8, min(1.0_r8, Dr))
  end if  

  ! 2.1 Linear Heuristic: 
  ! Loose (Dr=0) -> 15% Strain
  ! Dense (Dr=1) -> 3% Strain
  root_strain_pct = 15.0_r8 - (12.0_r8 * Dr)
  
  ! 3. Clay Adjustment
  ! Clays are more ductile than sands, increasing failure strain
  if (clay_in >= 0.3_r8) then
    root_strain_pct = root_strain_pct + 3.0_r8
  end if
  
  ! 4. Safety Bounds
  ! Rarely fails below 1% or above 20% in natural conditions
  root_strain_pct = max(1.0_r8, min(20.0_r8, root_strain_pct))

  ! 3. Calculate Effective Stress (Bishop's Method)
  ! Clamp Saturation 0.0 - 1.0
  Sr_clamped = max(0.0_r8, min(1.0_r8, rel_sat_in))
  
  ! Suction "Invisible Weight" = Sr * Suction
  sigma_suction_eff = Sr_clamped * suction_kPa_in
  
  ! Total Vertical Effective Stress
  sigma_v_eff = sigma_v_total_kPa_in + sigma_suction_eff
  
  ! 4. Calculate Mean Effective Stress (p')
  ! Assume K0 condition (At Rest) for natural soil
  K0 = K0_DEFAULT
  sigma_h_eff = K0 * sigma_v_eff
  
  ! Average of principal stresses: (sv + 2*sh) / 3
  p_prime = (sigma_v_eff + 2.0_r8 * sigma_h_eff) / 3.0_r8
  
  ! Safety floor to prevent log(<=0) errors
  p_prime = max(1.0_r8, p_prime)
  
  ! 5. Determine Base Friction (Phi_cv)
  phi_cv = PHI_CV_QUARTZ
  
  ! Texture Adjustment (Gravels interlock better)
  if (D50_mm_in > 2.0_r8) then
      phi_cv = phi_cv + 2.0_r8
  end if
  
  ! Structure/Organic Adjustment
  ! High e_max (>1.0) implies organic coatings or soft aggregation
  if (e_max_in > 1.0_r8) then
      reduction = (e_max_in - 1.0_r8) * 5.0_r8
      ! Limit reduction to max 3 degrees
      phi_cv = phi_cv - min(3.0_r8, reduction)
  end if
  
  ! 6. Calculate Peak Friction (Bolton's Dilatancy)
  ! Q parameter for natural aggregates (weaker than clean sand)
  Q_natural = Q_NATURAL_VAL
  
  ! Stress Suppression Term: Q - ln(p')
  ! As p' increases (due to suction), suppression increases, dilatancy drops.
  suppression_term = Q_natural - log(p_prime)
  
  ! Dilatancy cannot be negative (soil doesn't weaken below critical state)
  suppression_term = max(0.0_r8, suppression_term)
  
  ! 3.0 factor for Plane Strain conditions
  dilatancy = 3.0_r8 * Dr * suppression_term
  
  ! Cap dilatancy for natural aggregates (rarely exceeds 12 deg)
  dilatancy = min(dilatancy, 12.0_r8)
  
  phi_peak = phi_cv + dilatancy
  status = 0
  
  end subroutine estimate_friction_natural


! -------------------------------------------------------------------------

  function estimate_mineral_cohesion(clay_fraction, om_fraction, e_current) result(c_prime)
  !Description:
  ! Function: estimate_mineral_cohesion
  ! Purpose: Estimates Effective Cohesion (c') based on soil texture/state.
  !          Note: Does NOT include suction (suction is in Effective Stress).
  ! -------------------------------------------------------------------------
  implicit none
  real(r8), intent(in) :: clay_fraction  ! 0.0 to 1.0
  real(r8), intent(in) :: om_fraction    ! 0.0 to 1.0
  real(r8), intent(in) :: e_current      ! Void Ratio
  real(r8)             :: c_prime        ! Result in kPa
  
  real(r8) :: base_c
  real(r8) :: om_bonus,om_pct

  ! 1. Baseline based on Texture
  !Holtz, R. D., & Kovacs, W. D. (1981). An Introduction to Geotechnical Engineering. Prentice-Hall

  if (clay_fraction < 0.10_r8) then
    ! SAND: No cohesion
    base_c = 0.0_r8
  else if (clay_fraction < 0.35_r8) then
    ! LOAM/SILT: Weak cohesion (1-5 kPa)
    base_c = 2.0_r8 + (clay_fraction * 10.0_r8)
  else
    ! CLAY: Moderate cohesion (5-15 kPa)
    base_c = 5.0_r8 + (clay_fraction * 20.0_r8)
  end if
  
  ! 2. Density Adjustment
  ! Lower void ratio (denser soil) = Higher cohesion (particles closer)
  ! Simple multiplier: stronger if e is low (e.g. e=0.4), weaker if e is high (e=1.0)
  
  if (base_c > 0.0_r8) then
    c_prime = base_c * (0.8_r8 / max(0.3_r8, e_current))
  else
    c_prime = 0.0_r8
  end if
  
  ! --- 3. Organic Matter "Friability" Effect (Dilution) ---
  ! Zhang, B., Horn, R., & Hallett, P. D. (2005). "Mechanical resilience of degraded soil amended with organic matter." 
  ! Soil Science Society of America Journal, 69(3), 864-871.
  ! In heavy clays (>40%), high OM (>3%) prevents particles from locking 
  ! tightly, effectively reducing the mineral cohesion term (making soil crumbly).
  if (clay_fraction > 0.40_r8 .and. om_fraction > 0.03_r8) then
    base_c = base_c * 0.85_r8
  end if
        
  ! --- 4. Organic Matter "Glue" Effect (Binding) ---
  ! Ekwue, E. I. (1990). "Organic-matter effects on soil strength properties." Soil and Tillage Research, 16(3), 289-297.        
  ! OM acts as a binder, adding cohesion especially in sandy/loamy soils.
  ! Heuristic: ~0.8 kPa per 1% OM.
  ! Saturation: This effect caps around 6% OM.
  
  om_pct = om_fraction * 100.0_r8
  om_bonus = min(6.0_r8, om_pct) * 0.8_r8
  
  ! --- 5. Total and Safety Bounds ---
  c_prime = base_c + om_bonus
  
  ! Cap max cohesion for natural agricultural soils (rarely > 35 kPa)
  c_prime = min(35.0_r8, c_prime)      
  end function estimate_mineral_cohesion

  ! -------------------------------------------------------------------------
  subroutine estimate_root_stiffness(e_in, A_factor_in, &
    sigma_v_total_in, suction_in, rel_sat_in, friction_deg_in, &
    root_strain_pct_in, clay_fraction_in, om_fraction_in, &
    G_max, G_root, deg_factor, p_prime, status)
  !  
  ! Description:
  ! Subroutine: estimate_root_stiffness
  ! Purpose: Calculates G_max and the degraded G_root for biological modeling.
  ! -------------------------------------------------------------------------
  implicit none      
  ! --- Inputs ---
  real(r8), intent(in)  :: e_in               ! Void Ratio
  real(r8), intent(in)  :: A_factor_in        ! Stiffness Coefficient [sqrt(kPa)]
  real(r8), intent(in)  :: sigma_v_total_in   ! Total Vertical Stress [kPa]
  real(r8), intent(in)  :: suction_in         ! Matric Suction [kPa]
  real(r8), intent(in)  :: rel_sat_in         ! Relative Saturation (0.0 - 1.0)
  real(r8), intent(in)  :: friction_deg_in    ! Friction Angle [deg]
  real(r8), intent(in)  :: root_strain_pct_in ! Strain caused by root (e.g. 1.0)
  real(r8), intent(in)  :: clay_fraction_in   ! 0.0 to 1.0
  real(r8), intent(in)  :: om_fraction_in     ! 0.0 to 1.0
  
  ! --- Outputs ---
  real(r8), intent(out) :: G_max       ! Small Strain Modulus [kPa]
  real(r8), intent(out) :: G_root      ! EFFECTIVE ROOT STIFFNESS [kPa]
  real(r8), intent(out) :: deg_factor  ! Degradation Ratio (G_root / G_max)
  real(r8), intent(out) :: p_prime     ! Mean Effective Stress [kPa]
  integer, intent(out)  :: status      ! 0 = Success, -1 = Error
        
  real(r8) :: GAMMA_REF


  ! --- Local Variables ---

  real(r8) :: Sr_clamped, strain_clamped
  real(r8) :: sigma_v_eff, sigma_h_eff, K0
  real(r8) :: F_e, phi_rad
  
  ! 1. Validate Inputs
  if (e_in < 0.0_r8) then
    status = -1
    return
  end if

  ! Reference Strain for Natural Soils (0.07%)
  ! This defines the curve where stiffness drops to 50%
  GAMMA_REF = get_gamma_ref(clay_fraction_in, om_fraction_in)

  ! 2. Calculate Effective Stress State (Bishop's Method)
  Sr_clamped = max(0.0_r8, min(1.0_r8, rel_sat_in))
  
  ! Vertical Effective Stress (Gravity + Suction)
  sigma_v_eff = sigma_v_total_in + (Sr_clamped * suction_in)
  
  ! Horizontal Effective Stress (K0)
  phi_rad = friction_deg_in * RadianPerDegree
  K0 = 1.0_r8 - sin(phi_rad)
  K0 = max(0.1_r8, K0) ! Safety floor
  
  sigma_h_eff = K0 * sigma_v_eff
  
  ! Mean Effective Stress (p')
  p_prime = (sigma_v_eff + 2.0_r8 * sigma_h_eff) / 3.0_r8
  p_prime = max(1.0_r8, p_prime)
  
  ! 3. Calculate G_max (Small Strain - The "Seismic" Stiffness)
  ! Void Ratio Function
  if (e_in > 2.0_r8) then
    F_e = 1.0_r8 / (e_in**1.3_r8)
  else
    F_e = ((2.97_r8 - e_in)**2) / (1.0_r8 + e_in)
  end if
  
  ! Hardin & Black Equation
  G_max = A_factor_in * F_e * sqrt(p_prime)
  
  ! 4. Apply Hyperbolic Degradation for Root Strain
  ! Roots impose large strains (Plastic deformation zone).
  ! Model: G / G_max = 1 / (1 + (strain / gamma_ref))
  
  strain_clamped = max(0.0001_r8, root_strain_pct_in)
        
  deg_factor = 1.0_r8 / (1.0_r8 + (strain_clamped / GAMMA_REF))
  
  ! 5. Calculate Effective Root Stiffness
  G_root = G_max * deg_factor
  
  status = 0
        
  end subroutine estimate_root_stiffness


! -------------------------------------------------------------------------
  function get_gamma_ref(clay_fraction, om_fraction) result(gamma_ref)
  ! Function: get_gamma_ref
  ! Purpose: Calculates the Reference Strain (gamma_ref) for modulus degradation.
  !          Based on Vucetic & Dobry (1991) trends for PI.
  !          Adjusted for Organic Matter (OM increases linear elastic range).
  ! -------------------------------------------------------------------------

  implicit none
  real(r8), intent(in) :: clay_fraction ! 0.0 to 1.0
  real(r8), intent(in) :: om_fraction   ! 0.0 to 1.0
  real(r8)             :: gamma_ref     ! Result in % (e.g. 0.05)

  real(r8) :: PI_mineral, PI_effective

  ! 1. Estimate Mineral Plasticity Index (PI) from Clay
  ! Assumption: Skempton's Activity ~ 0.7 (Typical for mixed clays like Illite/Kaolinite)
  ! PI = Activity * %Clay
  ! Example: 20% Clay -> PI = 14
  PI_mineral = (clay_fraction * 100.0_r8) * 0.7_r8

  ! 2. Account for Organic Matter
  ! OM acts as a "ductility booster." 
  ! Heuristic: 1% OM behaves mechanically like adding ~2% Clay to the plasticity.
  ! It extends the linear range (makes soil spongier).
  PI_effective = PI_mineral + (om_fraction * 100.0_r8 * 2.0_r8)

  ! 3. Vucetic & Dobry (1991) Correlation Curve
  ! PI = 0  (Sand) -> gamma_ref ~ 0.03%
  ! PI = 15 (Loam) -> gamma_ref ~ 0.06%
  ! PI = 50 (Clay) -> gamma_ref ~ 0.15% - 0.20%
  
  ! Linear fit: 0.03 + (0.003 * PI)
  gamma_ref = 0.03_r8 + (0.003_r8 * PI_effective)

  ! 4. Bounds
  ! Rarely drops below 0.02% (Clean Quartz Sand)
  ! Rarely exceeds 0.25% (High Plasticity Clay / Peat)
  gamma_ref = max(0.02_r8, min(0.25_r8, gamma_ref))

  end function get_gamma_ref
end module SoilHydroParaMod
