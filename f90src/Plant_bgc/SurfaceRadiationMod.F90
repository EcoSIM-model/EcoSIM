module SurfaceRadiationMod
  ! Description:
  ! Prepare canopy geometry, optical properties and radiation, and diagnose
  ! exposed canopy and aerodynamic paths for the current hour.
  use data_kind_mod,      only: r8 => DAT_KIND_R8
  use minimathmod,        only: AZMAX1
  use PlantBGCPars,       only: iforward, ibackward
  use PrescribePhenolMod, only: SetCanopyProfile
  use EcoSIMCtrlMod,      only: ldo_sp_mode,ldo_radiation_test,etimer
  use RadiationDataMod
  use DebugToolMod
  use EcoSimConst
  use PlantMathFuncMod, only: CanopySnowOptics,CanopyRadiationFractions
  use EcoSIMConfig
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantRadiationAPIData, only : plt_rad
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantSoilChemistryAPIData, only : plt_soilchem
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: CanopyConditionModel
  public :: CalcBoundaryLayerProperties

  real(r8), parameter :: RAM                    = 2.78E-03_r8                               !minimum boundary layer resistance (h m-1)
  real(r8), parameter :: RadSWStemAlbedo(0:1)  = (/0.1_r8,0.1_r8/)                         !stalk albedo for shortwave
  real(r8), parameter :: RadPARStemAlbedo(0:1) = (/0.1_r8,0.1_r8/)                         !stalk albedo for PAR
  real(r8), parameter :: RadSWStemTransmitance(0:1) = (/0.12_r8,0.001_r8/)
  real(r8), parameter :: RadPARStemTransmitance(0:1)= (/0.12_r8,0.001_r8/)
  real(r8), parameter :: StemClumpFactor(0:1) = (/0.8_r8,0.5_r8/)          !stem clumping factor, FORGW = minimum SOC or organic soil (g Mg-1)
  real(r8), parameter :: SnowSWAlbedo         = 0.8_r8
  real(r8), parameter :: SnowPARAlbedo        = 0.95_r8
  real(r8), parameter :: SnowSWTransmitance   = 0._r8
  real(r8), parameter :: SnowPARTransmitance  = 0._r8
  real(r8), parameter :: StdeadSWAlbedo(0:1)  = (/0.25_r8,0.1_r8/)         !shortwave albedo for standing dead (grass, trees)
  real(r8), parameter :: StdeadPARAlbedo(0:1) = (/0.2_r8,0.14_r8/)         !PAR albedo for standing dead (grass, trees) 
  real(r8), parameter :: StdeadSWTransmitance(0:1)=(/0.1_r8,0.1_r8/)       !shortwave transmittance for standing dead (grass, trees)
  real(r8), parameter :: StdeadPARTransmitance(0:1)=(/0.4_r8,0.4_r8/)      !PAR transmittance for standing dead (grass, trees)
  
  integer :: iyrc
  logical :: dowrite

  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine CanopyConditionModel(I,J,DepthSurfWatIce)
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DepthSurfWatIce !water+ice depth at surface
  character(len=*), parameter :: subname='CanopyConditionModel'

  real(r8) :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  real(r8) :: StemAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  real(r8) :: SurfAreaZsecDead_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  call PrintInfo('beg '//subname)

  iyrc=etimer%get_curr_yearAD()

  if(ldo_sp_mode)then
    !do prescribed phenolgoy mode
    call SetCanopyProfile(I,J,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,SurfAreaZsecDead_lpft)
  else
    call DeriveCanopyHeightProfile(I,J,SurfAreaZsecDead_lpft)

    call SummarizeLeafStemAreaProfile(I,J,DepthSurfWatIce,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft)
  endif

  call UpdateOpticalProperties(I,J,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,SurfAreaZsecDead_lpft)

  call SurfaceRadiation(I,J,DepthSurfWatIce,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,SurfAreaZsecDead_lpft)

  call CalcBoundaryLayerProperties(DepthSurfWatIce)
  call PrintInfo('end '//subname)
  end subroutine CanopyConditionModel


!----------------------------------------------------------------------------------------------------
  subroutine CalcBoundaryLayerProperties(DepthSurfWatIce)
  ! Description:
  ! Diagnose canopy exposure above snow and compute neutral atmospheric and
  ! within-canopy transport properties using that shared exposure criterion.
  implicit none
  real(r8), intent(in) :: DepthSurfWatIce ! surface water/ice depth [m]; currently excluded from exposure test
  character(len=*), parameter :: subname='CalcBoundaryLayerProperties'
  real(r8) :: ARLSC
  real(r8) :: ARLSG
  real(r8) :: ZX,ZY,ZE
  REAL(R8) :: ZWind,HZnorm
  real(r8) :: WindH !wind speed at canopy height [m/h]
!     begin_execution
  call PrintInfo('beg '//subname)
  associate(                                                          &
    WindSpeedAtm_col            => plt_site%WindSpeedAtm_col         ,& !input  :wind speed, [m h-1]
    WindMesureHeight_col        => plt_site%WindMesureHeight_col     ,& !input  :wind speed measurement height, [m]
    ZERO                        => plt_site%ZERO                     ,& !input  :threshold zero for numerical stability, [-]
    AREA3                       => plt_site%AREA3                    ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    KoppenClimZone              => plt_site%KoppenClimZone           ,& !input  :Koppen climate zone for the grid,[-]
    SoilSurfRoughness_col       => plt_site%SoilSurfRoughness_col    ,& !input  :initial soil surface roughness height, [m]
    ZEROS                       => plt_site%ZEROS                    ,& !input  :threshold zero for numerical stability,[-]
    NU                          => plt_site%NU                       ,& !input  :current soil surface layer number, [-]
    TairK                       => plt_ew%TairK                      ,& !input  :air temperature, [K]
    canopy_exposed_col          => plt_ew%canopy_exposed_col         ,& !output :canopy area protrudes above snow, [-]
    VLHeatCapSnowMin_col        => plt_ew%VLHeatCapSnowMin_col       ,& !input  :minimum snowpack heat capacity, [MJ d-2 K-1]
    SnowDepth                   => plt_ew%SnowDepth                  ,& !input  :snowpack depth, [m]
    VLHeatCapSurfSnow_col       => plt_ew%VLHeatCapSurfSnow_col      ,& !input  :snowpack heat capacity, [MJ m-3 K-1]
    CanopyHeight_col            => plt_morph%CanopyHeight_col        ,& !input  :canopy height , [m]
    StemArea_col                => plt_morph%StemArea_col            ,& !input  :grid canopy stem area, [m2 d-2]
    CanopyLeafArea_col          => plt_morph%CanopyLeafArea_col      ,& !input  :grid canopy leaf area, [m2 d-2]
    StandDeadSurfArea_col       => plt_morph%StandDeadSurfArea_col   ,& !input  :grid canopy standing-dead area, [m2 d-2]
    RawIsoTSurf2CanopyHScal_col => plt_ew%RawIsoTSurf2CanopyHScal_col,& !output :scalar for isothermal aerodynamic resistance between zero-sink height and ground surface, [h m-1]
    RawIsoTAtm2CanopySinkZ_col  => plt_ew%RawIsoTAtm2CanopySinkZ_col ,& !output :isothermal aerodynamic resistance between zero-sink height and wind ref height in atmosphere, [h m-1]
    RawCanopyH2SinkZ_col        => plt_ew%RawCanopyH2SinkZ_col       ,& !output :isothermal aerodynamic resistance bewtween canopy height and zero sink height, [h m-1]
    ZeroPlaneDisplacem_col      => plt_ew%ZeroPlaneDisplacem_col     ,& !output :zero plane displacement height, [m]
    RoughnessLength             => plt_ew%RoughnessLength            ,& !output :canopy surface roughness height, [m]
    RIB                         => plt_ew%RIB                         & !output :Richardson number for calculating boundary layer resistance, [-]
  )
  !     CANOPY ZERO PLANE AND ROUGHNESS HEIGHTS
  !
  !     CanopyLeafArea_col,StemArea_col=leaf,stalk area of combined canopy
  !     SnowDepth,DepthSurfWatIce=snowpack,surface water depths
  !     ZWind=reference height for wind speed
  !
  ARLSC=CanopyLeafArea_col+StemArea_col+StandDeadSurfArea_col
!  IF(ARLSC.GT.ZEROS .AND. CanopyHeight_col.GE.SnowDepth-ZERO .AND. CanopyHeight_col.GE.DepthSurfWatIce-ZERO)THEN
  canopy_exposed_col=ARLSC.GT.ZEROS .AND. CanopyHeight_col.GE.SnowDepth-ZERO
  IF(canopy_exposed_col)then
    !
    ARLSG = ARLSC/AREA3(NU)
    ZX    = EXP(-0.5_r8*ARLSG)
    ZY    = 1.0_r8-ZX
    ZE    = CanopyHeight_col*AMAX1(0.05_r8,ZX*ZY)
    !
    ZeroPlaneDisplacem_col = CanopyHeight_col*AZMAX1(1.0_r8-2.0_r8/ARLSG*ZY)
  ELSE
    ZE                     = 0.0_r8
    ZeroPlaneDisplacem_col = 0.0_r8
  ENDIF
  !
  IF(iFlagRaiseZ0GbyVeg.EQ.iTrue)THEN !when wind measurement height in available.
    ZWind=WindMesureHeight_col+CanopyHeight_col
  ELSE
    ZWind=AMAX1(WindMesureHeight_col,ZeroPlaneDisplacem_col+2.0_r8)
  ENDIF

  IF(KoppenClimZone.GE.0)THEN
    !significant snowcover on ground
    IF(VLHeatCapSurfSnow_col.GT.VLHeatCapSnowMin_col)THEN 
      RoughnessLength=AMAX1(0.001_r8,ZE,ZW)
    ELSE
      RoughnessLength=AMAX1(0.001_r8,ZE,SoilSurfRoughness_col)
    ENDIF
    !
    !     CANOPY ISOTHERMAL BOUNDARY LAYER RESISTANCE
    !
    !   RawIsoTAtm2CanopySinkZ_col,RAM=biome canopy,minimum isothermal boundary layer resistance
    !   WindSpeedAtm_col=wind speed
    !   RIB=canopy isothermal Richardson number
    !   here 0.168 = cvonkarman**2
    !eq.(17) from Choudhury and Monteith (1988)
    RawIsoTAtm2CanopySinkZ_col = AMAX1(RAM,(LOG((ZWind-ZeroPlaneDisplacem_col)/RoughnessLength))**2/(0.168_r8*WindSpeedAtm_col))
    HZnorm                     = (CanopyHeight_col-ZeroPlaneDisplacem_col)/RoughnessLength

    if(canopy_exposed_col .and. HZnorm.GT.1._r8)then
      WindH=WindSpeedAtm_col*log(HZnorm)/(log(ZWind-ZeroPlaneDisplacem_col)-log(RoughnessLength))
      !eq. (24) from Choudhury and Monteith (1988)  
      RawIsoTSurf2CanopyHScal_col = log((CanopyHeight_col-ZeroPlaneDisplacem_col)/RoughnessLength) &
        /(0.168*(CanopyHeight_col-ZeroPlaneDisplacem_col)*WindH)
      RawCanopyH2SinkZ_col = AMAX1(RAM,(LOG(HZnorm))**2/(0.168_r8*WindH))
    else
      RawIsoTSurf2CanopyHScal_col = 0._r8
      WindH                       = WindSpeedAtm_col
      RawCanopyH2SinkZ_col        = 0._r8
    endif  
    !
    !1.27E+08_r8=g*(3600^2), from eq. (19) in Choudhury and Monteith (1988)
    RIB  = 1.27E+08_r8*(ZWind-RoughnessLength)/(WindSpeedAtm_col**2*TairK)
  ELSE
    RawIsoTAtm2CanopySinkZ_col  = RAM
    RawIsoTSurf2CanopyHScal_col = 0._r8
    RawCanopyH2SinkZ_col        = 0._r8
    RIB                         = 0.0_r8
    RoughnessLength             = 0._r8
  ENDIF
  end associate
  call PrintInfo('end '//subname)
  end subroutine CalcBoundaryLayerProperties

!----------------------------------------------------------------------------------------------------
  subroutine DeriveCanopyHeightProfile(I,J,SurfAreaZsecDead_lpft)
  implicit none
  integer, intent(in) :: I,J
  character(len=*), parameter :: subname='DeriveCanopyHeightProfile'
  real(r8), intent(out) :: SurfAreaZsecDead_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)

  real(r8) :: SLAProfileZ(NumCanopyLayers1),SLA
  real(r8) :: RSTK !stalk radius
  real(r8) :: EffHeightDead,FARSTD
  real(r8) :: CanopySurfAreaProfDeadZ_col(NumCanopyLayers1)
  integer :: NZ,L
  !     begin_execution
  associate(                                                            &
    NP                         => plt_site%NP                          ,& !input  :current number of plant species,[-]
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]        
    CanopyHeightLive_pft       => plt_morph%CanopyHeightLive_pft       ,& !input  :live plant canopy height, [m]
    CanopyHeightDead_pft       => plt_morph%CanopyHeightDead_pft       ,& !input  :canopy height for standing dead, [m]
    CanopyStemAareZ_col        => plt_morph%CanopyStemAareZ_col        ,& !input  :total stem area, [m2 d-2]
    CanopyLeafAareZ_col        => plt_morph%CanopyLeafAareZ_col        ,& !input  :total leaf area, [m2 d-2]
    StemArea_col               => plt_morph%StemArea_col               ,& !input  :grid canopy stem area, [m2 d-2]
    CanopyLeafArea_col         => plt_morph%CanopyLeafArea_col         ,& !input  :grid canopy leaf area, [m2 d-2]
    StemSpecVolume_pft         => plt_morph%StemSpecVolume_pft         ,& !input  :stalk specific volume, [m3 gC-1]        
    PlantPopuDead_pft          => plt_site%PlantPopuDead_pft           ,& !inoput :live+standing dead plant population, [d-2]            
    CanopyHeight_col           => plt_morph%CanopyHeight_col           ,& !inoput :canopy height , [m]
    CanopyHeightZ_col          => plt_morph%CanopyHeightZ_col          ,& !output :canopy layer height, [m]
    ZERO4Groth_pft             => plt_biom%ZERO4Groth_pft              ,& !input  :threshold zero for plang growth calculation, [-]        
    StandDeadSurfArea_col      => plt_morph%StandDeadSurfArea_col      ,& !output :grid canopy standing-dead area, [m2 d-2]
    StandDeadSurfArea_pft      => plt_morph%StandDeadSurfArea_pft      ,& !output :standing dead surface area, [m]    
    CanopySurfAreaProfDead_pft => plt_morph%CanopySurfAreaProfDead_pft ,& !output  :standing dead canopy surface area profile, [m2 d-2]    
    StandDeadStrutElms_pft     => plt_biom%StandDeadStrutElms_pft       & !output :standing dead element, [g d-2]    
  )
  call PrintInfo('beg '//subname)
  !
  !     DIVISION OF CANOPY INTO NumCanopyLayers LAYERS WITH EQUAL LAI
  !
  ! Standing-dead area and canopy-layer boundaries are mutually dependent.
  ! The dead-area profile is evaluated using the layer boundaries from the
  ! preceding calculation and is then used to update the current boundaries.
  ! This is one lagged fixed-point iteration, consistent with the legacy
  ! subhourly treatment. It avoids an iterative solve and assumes that canopy
  ! geometry changes gradually between calls.
  !
  SurfAreaZsecDead_lpft = 0._r8
  CanopySurfAreaProfDeadZ_col = 0._r8
  D9670: DO NZ=1,NP
    if(StandDeadStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ) .AND. PlantPopuDead_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
      CanopyHeightDead_pft(NZ) = AMAX1(1.e-2_r8,CanopyHeightDead_pft(NZ),CanopyHeightLive_pft(NZ))

      !stalk radius, assuming cylinderical shape. 
      RSTK                      = SQRT(StemSpecVolume_pft(NZ)*(StandDeadStrutElms_pft(ielmc,NZ)/PlantPopuDead_pft(NZ))/(PICON*CanopyHeightDead_pft(NZ)))
      StandDeadSurfArea_pft(NZ) = 2._R8*PICON*RSTK*CanopyHeightDead_pft(NZ)*PlantPopuDead_pft(NZ)
      EffHeightDead             = AMIN1(CanopyHeightDead_pft(NZ),CanopyHeightZ_col(NumCanopyLayers1))
      
      DO L=1,NumCanopyLayers1
        if(CanopyHeightDead_pft(NZ).GT.ZERO .and. CanopyHeightZ_col(L-1).LT.CanopyHeightDead_pft(NZ) &
          .and. CanopyHeightZ_col(L) .GT. CanopyHeightZ_col(L-1)) then
          FARSTD=AMIN1(1.0_r8,(CanopyHeightDead_pft(NZ)-CanopyHeightZ_col(L-1))/(CanopyHeightZ_col(L)-CanopyHeightZ_col(L-1)))
          CanopySurfAreaProfDead_pft(L,NZ)=FARSTD*StandDeadSurfArea_pft(NZ)*(CanopyHeightZ_col(L)-CanopyHeightZ_col(L-1))/EffHeightDead
          CanopySurfAreaProfDeadZ_col(L) = CanopySurfAreaProfDeadZ_col(L)+CanopySurfAreaProfDead_pft(L,NZ)
        else
          FARSTD=0._r8
          CanopySurfAreaProfDead_pft(L,NZ)=0._r8
        endif  
        SurfAreaZsecDead_lpft(:,L,NZ) = CanopySurfAreaProfDead_pft(L,NZ) &
          / real(NumLeafInclinationClasses1*NumOfLeafAzimuthSectors1,r8)
      ENDDO
    ELSE
      CanopyHeightDead_pft(NZ) = 0._r8
      StandDeadSurfArea_pft(NZ)= 0._r8  
      CanopySurfAreaProfDead_pft(1:NumCanopyLayers1,NZ) = 0._r8    
    ENDIF
  ENDDO D9670
  !
  !grid-scale canopy height
  CanopyHeight_col=0.0_r8
  D9685: DO NZ=1,NP
    CanopyHeight_col=AMAX1(CanopyHeight_col,CanopyHeightLive_pft(NZ),CanopyHeightDead_pft(NZ))
  ENDDO D9685

  !derive the grid level canopy profile (leaf+stem)
  DO L= 1, NumCanopyLayers1
    SLAProfileZ(L)=CanopyLeafAareZ_col(L)+CanopyStemAareZ_col(L)+CanopySurfAreaProfDeadZ_col(L) 
  ENDDO
  StandDeadSurfArea_col=SUM(StandDeadSurfArea_pft(1:NP))
  SLA = CanopyLeafArea_col+StemArea_col+StandDeadSurfArea_col

  call DeriveCanopyHeightZ(SLA, CanopyHeight_col,SLAProfileZ,CanopyHeightZ_col)

  call PrintInfo('end '//subname)
  end associate
  end subroutine DeriveCanopyHeightProfile

!----------------------------------------------------------------------------------------------------
  subroutine DeriveCanopyHeightZ(SLA_col, CanopyHeight_col,SLAProfileZ,CanopyHeightZ_col)
  implicit none
  real(r8), intent(in)   :: SLA_col
  real(r8), intent(in)   :: CanopyHeight_col
  real(r8), intent(in)   :: SLAProfileZ(NumCanopyLayers1) !stem+leaf Area in each layer
  real(r8), intent(inout):: CanopyHeightZ_col(0:NumCanopyLayers1)

  character(len=*), parameter :: subname='DeriveCanopyHeightZ'
  real(r8) :: CHeightZ1(0:NumCanopyLayers1)    !temporary copy of canopy height profile
  real(r8) :: AreaInterval,AreaL
  real(r8) :: ARX  !interval canopy area: leaf+stem
  real(r8) :: DZL  !canopy interval height
  integer  :: L
  associate(                                               &  
    ZEROS               => plt_site%ZEROS                  & !input  :threshold zero for numerical stability,[-]
  )

  call PrintInfo('beg '//subname)
  CanopyHeightZ_col(NumCanopyLayers1) = CanopyHeight_col+0.01_r8
  CHeightZ1(NumCanopyLayers1)         = CanopyHeightZ_col(NumCanopyLayers1) !top
  CHeightZ1(0)                        = 0.0_r8                              !bottom

  !divide total are into NumCanopyLayers1, from top to bottom
  AreaInterval=SLA_col/NumCanopyLayers1

  !the division assumes within each layer, total leaf area is uniformly distributed with a relative error of 0.01
  IF(AreaInterval.GT.ZEROS)THEN
    D2765: DO L=NumCanopyLayers1,2,-1
      !leaf+stem area in layer L
      AreaL=SLAProfileZ(L)

      !greater than the mean leaf area or area-interval
      IF(AreaL.GT.1.01_r8*AreaInterval)THEN
        !interval length for AreaL
        DZL  = CanopyHeightZ_col(L)-CanopyHeightZ_col(L-1)
        !
        CHeightZ1(L-1) = CanopyHeightZ_col(L-1)+0.5_r8*AMIN1(1.0_r8,(AreaL-AreaInterval)/AreaL)*DZL
      ELSEIF(AreaL.LT.0.99_r8*AreaInterval)THEN
        !leaf+stem area in layer L-1
        ARX = SLAProfileZ(L-1)

        !interval length for ARX
        DZL = CanopyHeightZ_col(L-1)-CanopyHeightZ_col(L-2)

        !layer L-1 has significant leaf+stem (canopy) area
        IF(ARX.GT.ZEROS)THEN
          CHeightZ1(L-1)=CanopyHeightZ_col(L-1)-0.5_r8*AMIN1(1.0_r8,(AreaInterval-AreaL)/ARX)*DZL
        ELSE
          CHeightZ1(L-1)=CanopyHeightZ_col(L-1)
        ENDIF
      ELSE
        CHeightZ1(L-1)=CanopyHeightZ_col(L-1)
      ENDIF
    ENDDO D2765

    D2770: DO L=NumCanopyLayers1,2,-1
      CanopyHeightZ_col(L-1)=CHeightZ1(L-1)
    ENDDO D2770
  ELSE
    DO L=NumCanopyLayers1,2,-1
      CanopyHeightZ_col(L-1)=0._r8
    ENDDO
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine DeriveCanopyHeightZ
!----------------------------------------------------------------------------------------------------
  subroutine SummarizeLeafStemAreaProfile(I,J,DepthSurfWatIce,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft)
  !
  !Description:
  !Summarize canopy leaf and stem area
  implicit none
  integer , intent(in) :: I,J
  real(r8), intent(in) :: DepthSurfWatIce   !surface water/ice thickness above soil surface
  character(len=*), parameter :: subname='SummarizeLeafStemAreaProfile'
  integer :: NZ,NB,L,K,N
  real(r8), intent(out) :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)    !leaf area in different angle sector
  real(r8), intent(out) :: StemAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)    !stem area in different angle sector
  
  associate(                                                            &
    NP                         => plt_site%NP                          ,& !input  :current number of plant species,[-]
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    SnowDepth                  => plt_ew%SnowDepth                     ,& !input  :snowpack depth, [m]
    NumOfBranches_pft          => plt_morph%NumOfBranches_pft          ,& !input  :number of branches,[-]
    LeafAreaZsec_brch          => plt_morph%LeafAreaZsec_brch          ,& !input  :leaf surface area, [m2 d-2]
    CanopyLeafArea_lnode       => plt_morph%CanopyLeafArea_lnode       ,& !input  :layer/node/branch leaf area, [m2 d-2]
    StemAreaZsec_brch          => plt_morph%StemAreaZsec_brch          ,& !input  :stem surface area, [m2 d-2]
    CanopyStalkSurfArea_lbrch  => plt_morph%CanopyStalkSurfArea_lbrch  ,& !input  :plant canopy layer branch stem area, [m2 d-2]
    CanopySurfAreaProfDead_pft => plt_morph%CanopySurfAreaProfDead_pft ,& !input  :standing dead canopy surface area profile, [m2 d-2]
    CanopyHeightZ_col          => plt_morph%CanopyHeightZ_col          ,& !input  :canopy layer height, [m]
    LeafStalkAreaAct_pft       => plt_morph%LeafStalkAreaAct_pft       ,& !output :radiation-active plant leaf+stem/stalk area, [m2 d-2]
    StandDeadSurfAreaAct_pft   => plt_morph%StandDeadSurfAreaAct_pft   ,& !output :radiation-active standing dead surface area, [m2 d-2]
    LeafStalkAreaAll_col       => plt_morph%LeafStalkAreaAll_col        & !output :live+dead leaf+stalk area of combined PFT canopy,[m^2 d-2]
  )
  call PrintInfo('beg '//subname)

  LeafStalkAreaAll_col=0.0_r8  
  
  D1135: DO NZ=1,NP
    StandDeadSurfAreaAct_pft(NZ) = 0.0_r8
    LeafStalkAreaAct_pft(NZ)     = 0.0_r8
    !move upward
    DO  L=1,NumCanopyLayers1    
        !above snow and water    
!        IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO .AND. CanopyHeightZ_col(L-1).GE.DepthSurfWatIce-ZERO)THEN
      IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO)then
        !above snow depth and above water/ice surface
        StandDeadSurfAreaAct_pft(NZ) = StandDeadSurfAreaAct_pft(NZ)+CanopySurfAreaProfDead_pft(L,NZ)
        LeafStalkAreaAll_col         = LeafStalkAreaAll_col+CanopySurfAreaProfDead_pft(L,NZ)
        DO  NB=1,NumOfBranches_pft(NZ)
          !
          !add all nodes over a branch
          D1130: DO K=1,MaxNodesPerBranch1
            LeafStalkAreaAct_pft(NZ) = LeafStalkAreaAct_pft(NZ)+CanopyLeafArea_lnode(L,K,NB,NZ)
            LeafStalkAreaAll_col     = LeafStalkAreaAll_col+CanopyLeafArea_lnode(L,K,NB,NZ)
          ENDDO D1130
          !add stem/stalk area
          LeafStalkAreaAct_pft(NZ) = LeafStalkAreaAct_pft(NZ)+CanopyStalkSurfArea_lbrch(L,NB,NZ)
          LeafStalkAreaAll_col     = LeafStalkAreaAll_col+CanopyStalkSurfArea_lbrch(L,NB,NZ)
        ENDDO  
      ENDIF      
    enddo    
  ENDDO D1135

  !summarize branches into different angle/depth classes
  D1150: DO NZ=1,NP
    DO  L=1,NumCanopyLayers1
      DO  N=1,NumLeafInclinationClasses1
        LeafAreaZsecLive_lpft(N,L,NZ)=0.0_r8
        StemAreaZsecLive_lpft(N,L,NZ)=0.0_r8
      enddo
    enddo
  ENDDO D1150

  D1200: DO NZ=1,NP
    DO  L=1,NumCanopyLayers1
!        IF(CanopyHeightZ_col(L-1).GT.SnowDepth-ZERO .AND. CanopyHeightZ_col(L-1).GT.DepthSurfWatIce-ZERO)THEN
      IF(CanopyHeightZ_col(L-1).GT.SnowDepth-ZERO)then

        DO  NB=1,NumOfBranches_pft(NZ)
          D1205: DO N=1,NumLeafInclinationClasses1
            D1210: DO K=1,MaxNodesPerBranch1
              LeafAreaZsecLive_lpft(N,L,NZ)=LeafAreaZsecLive_lpft(N,L,NZ)+LeafAreaZsec_brch(N,L,K,NB,NZ)
            ENDDO D1210
            StemAreaZsecLive_lpft(N,L,NZ)=StemAreaZsecLive_lpft(N,L,NZ)+StemAreaZsec_brch(N,L,NB,NZ)
          ENDDO D1205
        ENDDO
      ENDIF      
    enddo
  ENDDO D1200
  call PrintInfo('end '//subname)
  end associate
  end subroutine SummarizeLeafStemAreaProfile

!----------------------------------------------------------------------------------------------------
  subroutine SurfaceRadiation(I,J,DepthSurfWatIce,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,SurfAreaZsecDead_lpft)
  ! Description:
  ! Calculate multilayer shortwave absorption and the canopy/ground longwave
  ! partition using physical plus retained-snow area for both radiation paths.
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DepthSurfWatIce   !surface water/ice thickness above soil surface
  real(r8), intent(in) :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  real(r8), intent(in) :: StemAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  real(r8), intent(in) :: SurfAreaZsecDead_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)    
  character(len=*), parameter :: subname='SurfaceRadiation'

  real(r8) :: DGAZI
  real(r8) :: SolarAzimuthAngle,CosineSunInclAngle,cosGroundIncidentSolarAngle
  integer  :: NZ,N
  real(r8) :: GroundRadiationFraction
  real(r8) :: RadSW2Ground,RadPAR2Ground
  !     begin_execution
  call PrintInfo('beg '//subname)
  associate(                                                               &
    ZEROS                        => plt_site%ZEROS                        ,& !input  :threshold zero for numerical stability,[-]
    ZERO                         => plt_site%ZERO                         ,& !input  :threshold zero for numerical stability, [-]
    NP                           => plt_site%NP                           ,& !input  :current number of plant species,[-]
    NU                           => plt_site%NU                           ,& !input  :current soil surface layer number, [-]
    AREA3                        => plt_site%AREA3                        ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    CosineGrndSlope_col          => plt_rad%CosineGrndSlope_col           ,& !input  :cosine of slope, [-]
    SineGrndSlope_col            => plt_rad%SineGrndSlope_col             ,& !input  :sine of slope, [-]
    GroundSurfaceAzimuth_col     => plt_rad%GroundSurfaceAzimuth_col      ,& !input  :azimuth angle of the slope, [-]
    TotSineSkyAngles_grd         => plt_rad%TotSineSkyAngles_grd          ,& !input  :sine of sky angles,[-]
    OMEGA2Ground                 => plt_rad%OMEGA2Ground                  ,& !input  :sine of solar beam on leaf surface, [-]
    SolarNoonHour_col            => plt_site%SolarNoonHour_col            ,& !input  :time of solar noon, [h]
    SineSunInclinationAngle_col  => plt_rad%SineSunInclinationAngle_col   ,& !input  :sine of solar inclination angle, [-]
    LeafStalkAreaAct_pft         => plt_morph%LeafStalkAreaAct_pft        ,& !input  :plant leaf+stem/stalk area, [m2 d-2]    
    StandDeadSurfAreaAct_pft     => plt_morph%StandDeadSurfAreaAct_pft    ,& !input  :standing-dead surface area, [m2 d-2]
    LeafStalkAreaAll_col         => plt_morph%LeafStalkAreaAll_col        ,& !input  :stalk area of combined, each PFT canopy,[m^2 d-2]
    SnowOnCanopy_pft             => plt_ew%SnowOnCanopy_pft               ,& ! input: retained intercepted SWE [m3/cell, PFT]
    SnowType                    => plt_pheno%iPlantSnowIntercepType_pft ,& ! input: snow interception class [0:4, PFT]
    RadSWDirect_col              => plt_rad%RadSWDirect_col               ,& !inoput :direct shortwave radiation, [W m-2]
    RadSWDiffus_col              => plt_rad%RadSWDiffus_col               ,& !inoput :diffuse shortwave radiation, [W m-2]
    RadDirectPAR_col             => plt_rad%RadDirectPAR_col              ,& !inoput :direct PAR, [umol m-2 s-1]
    RadPARDiffus_col             => plt_rad%RadPARDiffus_col              ,& !inoput :diffuse PAR, [umol m-2 s-1]
    FracSWRad2Grnd_col           => plt_rad%FracSWRad2Grnd_col            ,& !inoput :fraction of radiation directed at ground surface, [-]
    RadSWCanopyAbsorption_pft    => plt_rad%RadSWCanopyAbsorption_pft     ,& !output :canopy absorbed shortwave radiation, [MJ d-2 h-1]
    RadSWCanopyLAbsroption_pft   => plt_rad%RadSWCanopyLAbsroption_pft    ,& !output :profile of canopy absorbed shortwave radiation, [MJ d-2 h-1]     
    RadPARCanopyLAbsorption_pft  => plt_rad%RadPARCanopyLAbsorption_pft   ,& !output :profile of canopy absorbed PAR, [MJ d-2 h-1]
    RadSWSolarBeam_col           => plt_rad%RadSWSolarBeam_col            ,& !output :shortwave radiation in solar beam, [MJ m-2 h-1]
    RadPARSolarBeam_col          => plt_rad%RadPARSolarBeam_col           ,& !output :PAR radiation in solar beam, [umol m-2 s-1]
    RadPARCanopyAbsorption_pft   => plt_rad%RadPARCanopyAbsorption_pft    ,& !output :canopy absorbed PAR, [umol m-2 s-1]
    FracPARads2Canopy_pft        => plt_rad%FracPARads2Canopy_pft         ,& !output :fraction of incoming PAR absorbed by total canopy, [-]
    FracPARads2LiveCanopy_pft    => plt_rad%FracPARads2LiveCanopy_pft     ,& !output :fraction of incoming PAR absorbed by live canopy, [-]
    RadSWGrnd_col                => plt_rad%RadSWGrnd_col                 ,& !output :radiation intercepted by ground surface, [MJ m-2 h-1]
    RadPARGrnd_col               => plt_rad%RadPARGrnd_col                ,& !output :PAR radiation reaching the ground, [umol m-2 s-1]
    ClumpFactorNow_pft           => plt_morph%ClumpFactorNow_pft           & !output :clumping factor for self-shading in canopy layer at current LAI, [-]
  )

  !
  !     MULTILAYER CANOPY INTERECEPTION OF DIRECT AND DIFFUSE RADIATION
  !     IN SW AND VISIBLE BANDS BY INCLINATION N, AZIMUTH M, LAYER L,
  !     NODE K, BRANCH NB, PFT NZ
  !
  !     LeafStalkAreaAct_pft,LeafStalkAreaAll_col=leaf+stalk area of combined,each PFT canopy
  !     ZL=height to bottom of canopy layer
  !     SnowDepth,DepthSurfWatIce=snowpack,surface water depths
  !     CanopyLeafArea_lnode,CanopyStalkSurfArea_lbrch=leaf,stalk areas of PFT
  !     RAD,RadPARSolarBeam_col=vertical direct+diffuse SW,PAR
  !     RADS,RADY,RAPS,RadPARDiffus_col=solar beam direct,diffuse SW,PAR
  !     SineSunInclinationAngle_col,TotSineSkyAngles_grd=sine of solar,sky angles
  !     ClumpFactorNow_pft=clumping factor for self-shading
  !
  IF(SineSunInclinationAngle_col.GT.ZERO)THEN
    RadSWSolarBeam_col  = RadSWDirect_col*SineSunInclinationAngle_col+RadSWDiffus_col*TotSineSkyAngles_grd
    RadPARSolarBeam_col = RadDirectPAR_col*SineSunInclinationAngle_col+RadPARDiffus_col*TotSineSkyAngles_grd
  ELSE
    RadSWDirect_col     = 0.0_r8
    RadSWDiffus_col     = 0.0_r8
    RadDirectPAR_col    = 0.0_r8
    RadPARDiffus_col    = 0.0_r8
    RadSWSolarBeam_col  = 0.0_r8
    RadPARSolarBeam_col = 0.0_r8
  ENDIF

  D1025: DO NZ=1,NP
    RadSWCanopyAbsorption_pft(NZ)     = 0.0_r8
    RadPARCanopyAbsorption_pft(NZ)    = 0.0_r8
    RadSWCanopyLAbsroption_pft(:,NZ)  = 0.0_r8
    RadPARCanopyLAbsorption_pft(:,NZ) = 0.0_r8
  ENDDO D1025
  call DebugPrint('SineSunInclinationAngle_col',SineSunInclinationAngle_col)
  !
  !     ANGLE BETWEEN SUN AND GROUND SURFACE
  !
  !     CosineSunInclAngle=solar azimuth,cosine of solar angle, radian, east be zero
  !     SolarAzimuthAngle uses the mathematical definition from Campbell and Norman, 1998
  !     cosGroundIncidentSolarAngle=cosine of incident solar angle at ground surface
  !     CosineGrndSlope_col,SineGrndSlope_col=cos,sin of ground surface
  !     SolarNoonHour_col=hour of solar noon from weather file
  !     0.2618=pi/12 (hrs), 4.7124=1.5*pi=270 degree
  IF(SineSunInclinationAngle_col.GT.ZERO)THEN
    !the solar azimuth angle is computed according to north hemisphere,
    SolarAzimuthAngle   = PICON*(SolarNoonHour_col-J)/12._r8+1.5_r8*PICON
    CosineSunInclAngle  = SQRT(1.0_r8-SineSunInclinationAngle_col**2)
    DGAZI               = COS(GroundSurfaceAzimuth_col-SolarAzimuthAngle)

    cosGroundIncidentSolarAngle = AZMAX1(AMIN1(1.0_r8,CosineGrndSlope_col*SineSunInclinationAngle_col+&
      SineGrndSlope_col*CosineSunInclAngle*DGAZI))

    !when there is canopy
    IF(LeafStalkAreaAll_col.GT.0.0_r8)THEN

      if(ldo_radiation_test)then
        RadSW2Ground=ABS(cosGroundIncidentSolarAngle)*RadSWDirect_col
        RadPAR2Ground=ABS(cosGroundIncidentSolarAngle)*RadDirectPAR_col
        D121: DO N=1,NumOfSkyAzimuthSects1
          RadSW2Ground=RadSW2Ground+ABS(OMEGA2Ground(N))*RadSWDiffus_col
          RadPAR2Ground=RadPAR2Ground+ABS(OMEGA2Ground(N))*RadPARDiffus_col
        ENDDO D121
      else
        call MultiCanopyLayerRadiation(I,J,DepthSurfWatIce,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,SurfAreaZsecDead_lpft,&
          SolarAzimuthAngle,CosineSunInclAngle,cosGroundIncidentSolarAngle,RadSW2Ground,RadPAR2Ground)
      endif
      !     RADIATION AT GROUND SURFACE IF NO CANOPY
    ELSE
      !plug in lake radiation below, direct beam
      RadSW2Ground=ABS(cosGroundIncidentSolarAngle)*RadSWDirect_col
      RadPAR2Ground=ABS(cosGroundIncidentSolarAngle)*RadDirectPAR_col
      D120: DO N=1,NumOfSkyAzimuthSects1
        RadSW2Ground=RadSW2Ground+ABS(OMEGA2Ground(N))*RadSWDiffus_col
        RadPAR2Ground=RadPAR2Ground+ABS(OMEGA2Ground(N))*RadPARDiffus_col
      ENDDO D120

    ENDIF
    RadSWGrnd_col=RadSW2Ground*AREA3(NU)
    RadPARGrnd_col=RadPAR2Ground*AREA3(NU)
    !
    !     IF NO RADIATION
    !
  ELSE
    RadSWGrnd_col=0.0_r8
    RadPARGrnd_col=0._r8
  ENDIF
  !
  !     CANOPY AND GROUND SKY FRACTIONS USED FOR BOUNDARY LAYER CALCULNS
  !
  !     FracSWRad2Grnd_col=fraction of radiation received by ground surface
  !     FRADP=fraction of radiation received by each PFT canopy
  !     LeafStalkAreaAll_col,LeafStalkAreaAct_pft=leaf+stalk area of all PFTs,each PFT
  !
  FracSWRad2Grnd_col=1.0_r8
  
  ! Use retained PFT snow rather than the column output sum cleared by HOUR1.
  ! PFT radiative weights and the column optical depth share the same added areas.
  call CanopyRadiationFractions(LeafStalkAreaAct_pft(1:NP),StandDeadSurfAreaAct_pft(1:NP), &
    SnowOnCanopy_pft(1:NP),SnowType(1:NP),AREA3(NU),FracPARads2Canopy_pft(1:NP), &
    FracPARads2LiveCanopy_pft(1:NP),GroundRadiationFraction)
  if(.not.ldo_radiation_test)FracSWRad2Grnd_col=GroundRadiationFraction
  end associate
  call PrintInfo('end '//subname)
  end subroutine SurfaceRadiation

!----------------------------------------------------------------------------------------------------
  subroutine MultiCanopyLayerRadiation(I,J,DepthSurfWatIce,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,SurfAreaZsecDead_lpft,&
    SolarAzimuthAngle,CosineSunInclAngle,cosGroundIncidentSolarAngle,RadSW2Ground,RadPAR2Ground)
  !
  !Description:
  ! Model multiple canopy layer radiation using bidirectional reflectance distribution function.
  !Ref: Grant et al., (1989), AFM, SIMULATION OF CANOPY PHOTOSYNTHESIS IN MAIZE AND SOYBEAN.
  !
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(in) :: DepthSurfWatIce   !surface water/ice thickness above soil surface
  real(r8), intent(in) :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  real(r8), intent(in) :: StemAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  real(r8), intent(in) :: SurfAreaZsecDead_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)      
  real(r8), intent(in) :: SolarAzimuthAngle
  real(r8), intent(in) :: CosineSunInclAngle
  real(r8), intent(in) :: cosGroundIncidentSolarAngle
  real(r8), intent(out):: RadSW2Ground                !radiation reaching the ground
  real(r8), intent(out):: RadPAR2Ground               !PAR radiation reaching the ground
  character(len=*), parameter :: subname='MultiCanopyLayerRadiation'

  integer :: NB,NZ,L,K,M,N,NN
  integer :: iScatteringDirect(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1)
  real(r8) :: TAU_RDiffusTransmitance(0:NumCanopyLayers1+1)
  real(r8) :: RadDirPAR2Leaf_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)
  real(r8) :: RadDirPAR2Stem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)
  real(r8) :: RadSWFwdScat2NextL(0:NumCanopyLayers1+1)   !total fwd scattered SW to next layer  
  real(r8) :: RadSWBakScat2NextL(0:NumCanopyLayers1+1)
  real(r8) :: RadPARBakScat2NextL(0:NumCanopyLayers1+1)
  real(r8) :: RadPARFwdScat2NextL(0:NumCanopyLayers1+1)
  
  real(r8) :: BETX(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1)                         !sine of direct solar radiation on leaf surface/sine of direct solar radiation, [-]
  REAL(R8) :: RadDirSW2Leaf_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)    !direct shortwave radiation on leave in each canopy sector
  real(r8) :: RadDirSW2Stem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)    !direct shortwave radiation on stalk in each canopy sector
  real(r8) :: RadDirSW2DStem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)   !direct shortwave radiation on dead stem in each canopy sector
  real(r8) :: RadDirPAR2DStem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)  !PAR radiation on dead stem in each canopy sector

  real(r8) :: RadSWDiffusL
  real(r8) :: RadPARDiffusL
  real(r8) :: XAREA
  real(r8) :: YAREA


  associate(                                                              &
    SineSunInclinationAngle_col  => plt_rad%SineSunInclinationAngle_col  ,& !input  :sine of solar angle, [-]
    NumOfBranches_pft            => plt_morph%NumOfBranches_pft          ,& !input  :number of branches,[-]    
    LeafAreaZsec_brch            => plt_morph%LeafAreaZsec_brch          ,& !input  :leaf surface area, [m2 d-2]    
    NU                           => plt_site%NU                          ,& !input  :current soil surface layer number, [-]
    AREA3                        => plt_site%AREA3                       ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    NP                           => plt_site%NP                          ,& !input  :current number of plant species,[-]
    ClumpFactorNow_pft           => plt_morph%ClumpFactorNow_pft         ,& !input  :clumping factor for self-shading in canopy layer at current LAI, [-]
    LeafEffArea_zsec             => plt_photo%LeafEffArea_zsec           ,& !inoput :leaf irradiated surface area, [m2 d-2]    
    RadTotPARAbsorption_zsec     => plt_rad%RadTotPARAbsorption_zsec     ,& !inoput :total incoming PAR, [umol m-2 s-1]
    RadDifPARAbsorption_zsec     => plt_rad%RadDifPARAbsorption_zsec     ,& !inoput :diffuse incoming PAR absorbed, [umol m-2 s-1]
    TAU_DirectSunLit             => plt_rad%TAU_DirectSunLit             ,& !inoput :fraction of radiation intercepted by canopy layer, [-]
    PARSunlit_pft                => plt_photo%PARSunlit_pft              ,& !inoput :PAR absorbed by sunlit leaf, [umol m-2 s-1]
    PARSunsha_pft                => plt_photo%PARSunsha_pft              ,& !inoput :PAR absorbed by sun-shaded leaf, [umol m-2 s-1]    
    LeafAreaSunlit_pft           => plt_photo%LeafAreaSunlit_pft         ,& !output :leaf irradiated surface area, [m2 d-2]        
    TAU_DirectSunSha             => plt_rad%TAU_DirectSunSha              & !output :fraction of radiation transmitted by canopy layer, [-]
  )
  
  call PrintInfo('beg '//subname)

  !distribute radiation into different leaf/canopy sector
  CALL AngularDistributeRadiation(SolarAzimuthAngle,CosineSunInclAngle,BETX,iScatteringDirect, &
    RadDirSW2Leaf_zsec,RadDirSW2Stem_zsec,RadDirPAR2Leaf_zsec,RadDirPAR2Stem_zsec)
  RadDirSW2DStem_zsec  = RadDirSW2Stem_zsec
  RadDirPAR2DStem_zsec = RadDirPAR2Stem_zsec

  XAREA = 1.0_r8/AREA3(NU)                      !average over the grid
  YAREA = XAREA/REAL(NumOfSkyAzimuthSects1,R8)  !area for one azimuthal zone
  !
  !distribute radiation from top of canopy to ground
  !it is only applicable to upland vegetation
  !
  !     RESET ARRAYS OF SUNLIT AND SHADED LEAF AREAS IN DIFFERENT
  !     LAYERS AND ANGLE CLASSES
  !
  !     TSURF,StemAreaZsecLive_lpft,SURF,StemAreaZsec_brch=leaf,stalk total,PFT surface area
  !
  !
  !     CALCULATE ABSORPTION, REFLECTION AND TRANSMISSION OF DIRECT AND
  !     DIFFUSE DOWNWARD TOTAL AND VISIBLE RADIATION BY EACH SPECIES
  !     NZ IN EACH LAYER L
  !  
  !direct radiation only applies to sunlit leaf/stalk, diffusive radiation applies to all.
  !
  call DownwellSweep(DepthSurfWatIce,XAREA,YAREA,BETX,iScatteringDirect,         &
    LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,SurfAreaZsecDead_lpft,           &
    RadDirSW2Leaf_zsec,RadDirSW2Stem_zsec,RadDirSW2DStem_zsec,                   &
    RadDirPAR2Leaf_zsec,RadDirPAR2Stem_zsec,RadDirPAR2DStem_zsec,                &
    RadSWDiffusL,RadPARDiffusL,                                                   &
    TAU_RDiffusTransmitance,RadSWFwdScat2NextL,RadSWBakScat2NextL,               &
    RadPARFwdScat2NextL,RadPARBakScat2NextL)

  !set lower boundary condition
  call RadiationAtGround(cosGroundIncidentSolarAngle, RadSWDiffusL,RadPARDiffusL, &
    TAU_DirectSunLit(1),TAU_RDiffusTransmitance(1),                               &
    RadSWFwdScat2NextL(1),RadPARFwdScat2NextL(1),                                 &
    RadSWBakScat2NextL(0),RadPARBakScat2NextL(0),RadSW2Ground,RadPAR2Ground)
  !
  !     ADD RADIATION FROM SCATTERING THROUGH CANOPY LAYERS
  !
  call UpwellSweep(YAREA,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,  &
    SurfAreaZsecDead_lpft, TAU_RDiffusTransmitance,RadSWBakScat2NextL, &
    RadPARBakScat2NextL,RadSWFwdScat2NextL,RadPARFwdScat2NextL)
  !
  !summarize sunlit/shaded PAR
  DO NZ=1,NP
    LeafAreaSunlit_pft(NZ) = 0._r8;PARSunlit_pft(NZ) = 0._r8;ParSunsha_pft(NZ) = 0._r8

    D500: DO NB=1,NumOfBranches_pft(NZ)
      D550: DO K=1,MaxNodesPerBranch1
        D600: DO L=NumCanopyLayers1,1,-1
          D650: DO N=1,NumLeafInclinationClasses1
            LeafEffArea_zsec(N,L,K,NB,NZ) = LeafAreaZsec_brch(N,L,K,NB,NZ)*ClumpFactorNow_pft(NZ)            
            DO M=1,NumOfSkyAzimuthSects1
              LeafAreaSunlit_pft(NZ) = LeafAreaSunlit_pft(NZ)+LeafEffArea_zsec(N,L,K,NB,NZ)*TAU_DirectSunLit(L+1)
              PARSunlit_pft(NZ)      = PARSunlit_pft(NZ)+RadTotPARAbsorption_zsec(N,M,L,NZ)*LeafEffArea_zsec(N,L,K,NB,NZ)*TAU_DirectSunLit(L+1)
              ParSunsha_pft(NZ)      = ParSunsha_pft(NZ)+RadDifPARAbsorption_zsec(N,M,L,NZ)*LeafEffArea_zsec(N,L,K,NB,NZ)*TAU_DirectSunSha(L+1)
            ENDDO
          ENDDO D650
        ENDDO D600
      ENDDO D550 
    ENDDO D500
  ENDDO
  call PrintInfo('end '//subname)
  end associate
  end subroutine MultiCanopyLayerRadiation
!----------------------------------------------------------------------------------------------------

  subroutine DownwellSweep(DepthSurfWatIce,XAREA,YAREA,BETX,                   &
    iScatteringDirect, LeafAreaZsecLive_lpft, StemAreaZsecLive_lpft,           &
    SurfAreaZsecDead_lpft, RadDirSW2Leaf_zsec,RadDirSW2Stem_zsec,              &
    RadDirSW2DStem_zsec,RadDirPAR2Leaf_zsec,RadDirPAR2Stem_zsec,               &
    RadDirPAR2DStem_zsec,RadSWDiffusL,RadPARDiffusL,                           &
    TAU_RDiffusTransmitance,RadSWFwdScat2NextL,RadSWBakScat2NextL,             &
    RadPARFwdScat2NextL,RadPARBakScat2NextL)

  implicit none
  real(r8), intent(in)    :: DepthSurfWatIce
  real(r8), intent(in)    :: XAREA,YAREA
  real(r8), intent(in)    :: BETX(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1)                         !sine of direct solar radiation on leaf surface/sine of direct solar radiation, [-]  
  integer,  INTENT(IN)    :: iScatteringDirect(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1)    
  real(r8), intent(in)    :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  real(r8), intent(in)    :: StemAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  real(r8), intent(in)    :: SurfAreaZsecDead_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)        
  REAL(R8), intent(in)    :: RadDirSW2Leaf_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)     !direct shortwave radiation on leave in each canopy sector
  real(r8), intent(in)    :: RadDirSW2Stem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)    !direct shortwave radiation on stalk in each canopy sector
  real(r8), intent(in)    :: RadDirSW2DStem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)   !direct shortwave radiation on dead stem in each canopy sector
  real(r8), intent(in)    :: RadDirPAR2Leaf_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)    !PAR radiation on leave in each canopy sector
  real(r8), intent(in)    :: RadDirPAR2Stem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)   !PAR radiation on stalk in each canopy sector
  real(r8), intent(in)    :: RadDirPAR2DStem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)  !PAR radiation on dead stem in each canopy sector
  real(r8), intent(out)   :: RadSWDiffusL,RadPARDiffusL
  real(r8), intent(out)   :: TAU_RDiffusTransmitance(0:NumCanopyLayers1+1)
  real(r8), intent(out)   :: RadSWFwdScat2NextL(0:NumCanopyLayers1+1)   !total fwd scattered SW to next layer
  real(r8), intent(out)   :: RadSWBakScat2NextL(0:NumCanopyLayers1+1)   !total backscattered SW to next layer
  real(r8), intent(out)   :: RadPARFwdScat2NextL(0:NumCanopyLayers1+1)  !total fwd scattered PAR to next layer
  real(r8), intent(out)   :: RadPARBakScat2NextL(0:NumCanopyLayers1+1)  !total backscattered PAR to next layer
  character(len=*), parameter :: subname='DownwellSweep'
  real(r8) :: RadDirSWLeafLAbsorption_pft(JP1)
  real(r8) :: RadDirPARLeafLAbsorption_pft(JP1)
  real(r8) :: RadDifSWStemAbsorption_pft(JP1)
  real(r8) :: RadDifPARStemAbsorption_pft(JP1)
  real(r8) :: RadDirSWStemLAbsorption_pft(JP1)
  real(r8) :: RadDirPARStemLAbsorption_pft(JP1)
  real(r8) :: RadDirSWDStemLAbsorption_pft(JP1)
  real(r8) :: RadDirPARDStemLAbsorption_pft(JP1)
  real(r8) :: RadDifSWLeafAbsorption_pft(JP1),RadDifPARLeafAbsorption_pft(JP1)  
  real(r8) :: RadDifSWDStemAbsorption_pft(JP1),RadDifPARDStemAbsorption_pft(JP1)
  real(r8) :: backScatDirSWLeafLAbsorpt_pft(JP1),FwdScatDirSWLeafLAbsorpt_pft(JP1)
  real(r8) :: backScatDirPARLeafLAbsorpt_pft(JP1),FwdScatDirPARLeafLAbsorpt_pft(JP1)
  real(r8) :: backScatDirPARStemLAbsorpt_pft(JP1),FwdScatDirPARStemLAbsorpt_pft(JP1)
  real(r8) :: backScatDirSWStemLAbsorpt_pft(JP1),FwdScatDirSWStemLAbsorpt_pft(JP1)
  real(r8) :: backScatDirSWDStemLAbsorpt_pft(JP1),FwdScatDirSWDStemLAbsorpt_pft(JP1)
  real(r8) :: backScatDirPARDStemLAbsorpt_pft(JP1),FwdScatDirPARDStemLAbsorpt_pft(JP1)
  real(r8) :: BackScatDifSWLeafLAbsorpt_pft(JP1),FwdScatDifSWLeafLAbsorpt_pft(JP1)
  real(r8) :: BackScatDifPARLeafLAbsorpt_pft(JP1),FwdScatDifPARLeafLAbsorpt_pft(JP1)
  real(r8) :: BackScatDifSWStemLAbsorpt_pft(JP1),FwdScatDifSWStemLAbsorpt_pft(JP1)
  real(r8) :: BackScatDifPARStemLAbsorpt_pft(JP1),FwdScatDifPARStemLAbsorpt_pft(JP1)
  real(r8) :: BackScatDifSWDStemLAbsorpt_pft(JP1),FwdScatDifSWDStemLAbsorpt_pft(JP1)
  real(r8) :: BackScatDifPARDStemLAbsorpt_pft(JP1),FwdScatDifPARDStemLAbsorpt_pft(JP1)

  integer :: L, NZ
  real(r8) :: FractionDirRadAbsorbtCum, FractionDifRadAbsorbtCum
  real(r8) :: FractionDirRadAbsorbtL, FractionDifRadAbsorbtL

  associate(                                                        &
    SnowDepth                  => plt_ew%SnowDepth                 ,& !input  :snowpack depth, [m]  
    NP                         => plt_site%NP                      ,& !input  :current number of plant species,[-]    
    ZERO                       => plt_site%ZERO                    ,& !input  :threshold zero for numerical stability, [-]    
    RadSWDiffus_col            => plt_rad%RadSWDiffus_col          ,& !input  :diffuse shortwave radiation, [W m-2]    
    RadPARDiffus_col           => plt_rad%RadPARDiffus_col         ,& !input  :diffuse PAR, [umol m-2 s-1]    
    CanopyHeightZ_col          => plt_morph%CanopyHeightZ_col      ,& !input  :canopy layer height, [m]      
    TAU_DirectSunSha           => plt_rad%TAU_DirectSunSha         ,& !output :fraction of radiation transmitted by canopy layer, [-]    
    TAU_DirectSunLit           => plt_rad%TAU_DirectSunLit          & !inoput :fraction of radiation intercepted by canopy layer, [-]
  )
  call PrintInfo('beg '//subname)

  !from the sky to canopy top
  !incoming SW and PAR
  RadSWDiffusL  = RadSWDiffus_col           
  RadPARDiffusL = RadPARDiffus_col          

  FractionDirRadAbsorbtCum                    = 0.0_r8
  TAU_DirectSunLit(NumCanopyLayers1+1)        = 1.0_r8
  TAU_RDiffusTransmitance(NumCanopyLayers1+1) = 1.0_r8
  RadSWFwdScat2NextL(NumCanopyLayers1+1)      = 0.0_r8
  RadPARFwdScat2NextL(NumCanopyLayers1+1)     = 0.0_r8

  D1800: DO L=NumCanopyLayers1,1,-1 !top -> down
    !only compute radiation for canopy layers above snow and ponding water
    !IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO .AND. CanopyHeightZ_col(L-1).GE.DepthSurfWatIce-ZERO)THEN
    IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO)THEN
      !
      !incoming radiation, shortwave /par, due to transmission and downward scattering from L+1 to L
      RadSWDiffusL           = RadSWDiffusL *TAU_RDiffusTransmitance(L+1)+RadSWFwdScat2NextL(L+1)
      RadPARDiffusL          = RadPARDiffusL*TAU_RDiffusTransmitance(L+1)+RadPARFwdScat2NextL(L+1)
      !
      !initialize to zero 
      RadSWFwdScat2NextL(L)    = 0.0_r8
      RadPARFwdScat2NextL(L)   = 0.0_r8
      RadSWBakScat2NextL(L)    = 0.0_r8
      RadPARBakScat2NextL(L)   = 0.0_r8
      FractionDifRadAbsorbtCum = 0.0_r8
      FractionDirRadAbsorbtL   = 0.0_r8
      FractionDifRadAbsorbtL   = 0.0_r8
      !
      !  RESET ACCUMULATORS OB ABSORBED, REFLECTED AND TRANSMITTED RADIATION
      !
      !
      D1500: DO NZ=1,NP
        !plant radiation in canopy layer L
        CALL DownwellRadiationCanopyZL(NZ,L,XAREA,YAREA,BETX,iScatteringDirect,                                       &
            LeafAreaZsecLive_lpft(:,L,NZ),StemAreaZsecLive_lpft(:,L,NZ), SurfAreaZsecDead_lpft(:,L,NZ),               &
            RadSWDiffusL,RadPARDiffusL,FractionDirRadAbsorbtL,FractionDifRadAbsorbtL,                                 &
            RadDirSW2Leaf_zsec,RadDirSW2Stem_zsec,RadDirSW2DStem_zsec,                                                &
            RadDirPAR2Leaf_zsec, RadDirPAR2Stem_zsec,RadDirPAR2DStem_zsec,                                            &
            RadDirSWLeafLAbsorption_pft(NZ),backScatDirSWLeafLAbsorpt_pft(NZ),FwdScatDirSWLeafLAbsorpt_pft(NZ),       &
            RadDirSWStemLAbsorption_pft(NZ),backScatDirSWStemLAbsorpt_pft(NZ),FwdScatDirSWStemLAbsorpt_pft(NZ),       &
            RadDirSWDStemLAbsorption_pft(NZ),backScatDirSWDStemLAbsorpt_pft(NZ),FwdScatDirSWDStemLAbsorpt_pft(NZ),    &    
            RadDirPARLeafLAbsorption_pft(NZ),backScatDirPARLeafLAbsorpt_pft(NZ),FwdScatDirPARLeafLAbsorpt_pft(NZ),    &
            RadDirPARStemLAbsorption_pft(NZ),backScatDirPARStemLAbsorpt_pft(NZ),FwdScatDirPARStemLAbsorpt_pft(NZ),    &
            RadDirPARDStemLAbsorption_pft(NZ),backScatDirPARDStemLAbsorpt_pft(NZ),FwdScatDirPARDStemLAbsorpt_pft(NZ), &
            RadDifSWLeafAbsorption_pft(NZ),BackScatDifSWLeafLAbsorpt_pft(NZ),FwdScatDifSWLeafLAbsorpt_pft(NZ),        &
            RadDifSWStemAbsorption_pft(NZ),BackScatDifSWStemLAbsorpt_pft(NZ),FwdScatDifSWStemLAbsorpt_pft(NZ),        &
            RadDifSWDStemAbsorption_pft(NZ),BackScatDifSWDStemLAbsorpt_pft(NZ),FwdScatDifSWDStemLAbsorpt_pft(NZ),     &
            RadDifPARLeafAbsorption_pft(NZ),BackScatDifPARLeafLAbsorpt_pft(NZ),FwdScatDifPARLeafLAbsorpt_pft(NZ),     &
            RadDifPARStemAbsorption_pft(NZ),BackScatDifPARStemLAbsorpt_pft(NZ),FwdScatDifPARStemLAbsorpt_pft(NZ),     &
            RadDifPARDStemAbsorption_pft(NZ),BackScatDifPARDStemLAbsorpt_pft(NZ),FwdScatDifPARDStemLAbsorpt_pft(NZ))
   
      ENDDO D1500
      !
      !     ACCUMULATED INTERCEPTION BY CANOPY LAYER
      !
      call AccumulateRadation4CanopyL(L,YAREA,RadDirPAR2Leaf_zsec,                                                           &
        FractionDirRadAbsorbtCum,FractionDifRadAbsorbtCum,                                                                   &
        FractionDirRadAbsorbtL,FractionDifRadAbsorbtL,RadDirSWLeafLAbsorption_pft,RadDirPARLeafLAbsorption_pft,              &
        RadDirSWStemLAbsorption_pft,RadDirPARStemLAbsorption_pft,RadDirSWDStemLAbsorption_pft,RadDirPARDStemLAbsorption_pft, &
        RadDifSWLeafAbsorption_pft,RadDifPARLeafAbsorption_pft,RadDifSWStemAbsorption_pft,RadDifPARStemAbsorption_pft,       &
        RadDifSWDStemAbsorption_pft,RadDifPARDStemAbsorption_pft,                                                            &
        backScatDirSWLeafLAbsorpt_pft,FwdScatDirSWLeafLAbsorpt_pft,backScatDirSWStemLAbsorpt_pft,FwdScatDirSWStemLAbsorpt_pft, &
        backScatDirSWDStemLAbsorpt_pft,FwdScatDirSWDStemLAbsorpt_pft,                                                        &
        backScatDirPARLeafLAbsorpt_pft,FwdScatDirPARLeafLAbsorpt_pft,backScatDirPARStemLAbsorpt_pft,FwdScatDirPARStemLAbsorpt_pft, &
        backScatDirPARDStemLAbsorpt_pft,FwdScatDirPARDStemLAbsorpt_pft,                                                      &
        BackScatDifSWLeafLAbsorpt_pft,FwdScatDifSWLeafLAbsorpt_pft,BackScatDifSWStemLAbsorpt_pft,FwdScatDifSWStemLAbsorpt_pft, &
        BackScatDifSWDStemLAbsorpt_pft,FwdScatDifSWDStemLAbsorpt_pft,                                                        &
        BackScatDifPARLeafLAbsorpt_pft,FwdScatDifPARLeafLAbsorpt_pft,BackScatDifPARStemLAbsorpt_pft,FwdScatDifPARStemLAbsorpt_pft, &
        BackScatDifPARDStemLAbsorpt_pft,FwdScatDifPARDStemLAbsorpt_pft,                                                      &
        RadSWBakScat2NextL,RadPARBakScat2NextL,RadSWFwdScat2NextL,RadPARFwdScat2NextL,TAU_RDiffusTransmitance)

      FractionDirRadAbsorbtCum   = FractionDirRadAbsorbtCum+FractionDirRadAbsorbtL
      FractionDifRadAbsorbtCum   = FractionDifRadAbsorbtCum+FractionDifRadAbsorbtL
      TAU_DirectSunLit(L)        = 1.0_r8-FractionDirRadAbsorbtCum
      TAU_DirectSunSha(L)        = 1.0_r8-TAU_DirectSunLit(L)
      TAU_RDiffusTransmitance(L) = 1.0_r8-FractionDifRadAbsorbtCum
      !
    ELSE
      !layer L is not in active canopy, so just make a copy from layer L+1
      RadSWFwdScat2NextL(L)      = RadSWFwdScat2NextL(L+1)
      RadPARFwdScat2NextL(L)     = RadPARFwdScat2NextL(L+1)
      TAU_DirectSunLit(L)        = TAU_DirectSunLit(L+1)
      TAU_DirectSunSha(L)        = 1.0_r8-TAU_DirectSunLit(L)
      TAU_RDiffusTransmitance(L) = TAU_RDiffusTransmitance(L+1)
    ENDIF
  ENDDO D1800
  call PrintInfo('end '//subname)
  end associate
  end subroutine DownwellSweep

!----------------------------------------------------------------------------------------------------

  subroutine UpwellSweep(YAREA,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,&
    SurfAreaZsecDead_lpft, TAU_RDiffusTransmitance,RadSWBakScat2NextL,&
    RadPARBakScat2NextL,RadSWFwdScat2NextL,RadPARFwdScat2NextL)

  implicit none
  real(r8), intent(in)     :: YAREA
  real(r8), intent(in)     :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)
  real(r8), intent(in)     :: StemAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)  
  real(r8), intent(in)     :: SurfAreaZsecDead_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)          
  real(r8), intent(inout)  :: TAU_RDiffusTransmitance(0:NumCanopyLayers1+1)
  real(r8), intent(inout)  :: RadSWBakScat2NextL(0:NumCanopyLayers1+1)
  real(r8), intent(inout)  :: RadPARBakScat2NextL(0:NumCanopyLayers1+1)
  real(r8), intent(inout)  :: RadSWFwdScat2NextL(0:NumCanopyLayers1+1)
  real(r8), intent(inout)  :: RadPARFwdScat2NextL(0:NumCanopyLayers1+1)
  character(len=*), parameter :: subname='UpwellSweep'

  REAL(R8) :: EffLeafArea,EffStemArea,EffDeadStemArea  
  real(r8) :: RadSWDiffusL,DiffSWLeafAbsorptAzclass,DiffSWStemAbsorptAzclass,DiffSWDStemAbsorptAzclass
  real(r8) :: RadPARDiffusL,DiffPARLeafAbsorptAzclass,DiffPARStemAbsorptAzclass,DiffPARDStemAbsorptAzclass
  real(r8) :: RadDifSWLeafAbsorption_pft(JP1)
  real(r8) :: RadDifPARLeafAbsorption_pft(JP1) 
  real(r8) :: RadDifSWStemAbsorption_pft(JP1)
  real(r8) :: RadDifPARStemAbsorption_pft(JP1)  
  real(r8) :: RadDifSWDStemAbsorption_pft(JP1)
  real(r8) :: RadDifPARDStemAbsorption_pft(JP1)  
  integer  :: L,NZ,M,N,NN

  associate(                                                            &
    OMEGA2Leaf                 => plt_rad%OMEGA2Leaf                   ,& !input  :sine of indirect sky radiation on leaf surface,[-]  
    NP                         => plt_site%NP                          ,& !input  :current number of plant species,[-]    
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]
    iPlant2ndGrothPattern_pft  => plt_pheno%iPlant2ndGrothPattern_pft  ,& !input  :plant expression of secondary growth, [-]            
    ClumpFactorNow_pft         => plt_morph%ClumpFactorNow_pft         ,& !input  :clumping factor for self-shading in canopy layer at current LAI, [-]    
    RadSWCanopyAbsorption_pft  => plt_rad%RadSWCanopyAbsorption_pft    ,& !inoput :canopy absorbed shortwave radiation, [MJ d-2 h-1]
    RadPARCanopyAbsorption_pft => plt_rad%RadPARCanopyAbsorption_pft   ,& !inoput :canopy absorbed PAR, [umol m-2 s-1]        
    RadSWCanopyLAbsroption_pft => plt_rad%RadSWCanopyLAbsroption_pft   ,& !output :profile of canopy absorbed shortwave radiation, [MJ d-2 h-1]         
    RadPARCanopyLAbsorption_pft=> plt_rad%RadPARCanopyLAbsorption_pft  ,& !output :profile of canopy absorbed PAR, [umol m-2 s-1]
    RadTotPARAbsorption_zsec   => plt_rad%RadTotPARAbsorption_zsec     ,& !inoput :total incoming PAR, [umol m-2 s-1]
    RadDifPARAbsorption_zsec   => plt_rad%RadDifPARAbsorption_zsec     ,& !inoput :diffuse incoming PAR, [umol m-2 s-1]
    LeafPARabsorptivity_pft    => plt_rad%LeafPARabsorptivity_pft      ,& !input  :canopy PAR absorptivity,[-]
    LeafSWabsorptivity_pft     => plt_rad%LeafSWabsorptivity_pft       ,& !input  :canopy shortwave absorptivity, [-]    
    SnowDepth                  => plt_ew%SnowDepth                     ,& !input  :snowpack depth, [m]    
    CanopyHeightZ_col          => plt_morph%CanopyHeightZ_col           & !input  :canopy layer height, [m]    
  )
  call PrintInfo('beg '//subname)
  !
  TAU_RDiffusTransmitance(0) = 1.0_r8
  RadSWFwdScat2NextL(0)  = 0.0_r8
  RadPARFwdScat2NextL(0) = 0.0_r8  
  RadSWDiffusL           = 0.0_r8
  RadPARDiffusL          = 0.0_r8
  D2800: DO L=1,NumCanopyLayers1 !bottom->up
    !the following line shuts off radiation when it is below water 
!    IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO .AND. CanopyHeightZ_col(L-1).GE.DepthSurfWatIce-ZERO)THEN
    IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO)THEN

      RadSWDiffusL  = RadSWDiffusL*TAU_RDiffusTransmitance(L-1)+RadSWFwdScat2NextL(L-1)+RadSWBakScat2NextL(L-1)
      RadPARDiffusL = RadPARDiffusL*TAU_RDiffusTransmitance(L-1)+RadPARFwdScat2NextL(L-1)+RadPARBakScat2NextL(L-1)

      RadSWFwdScat2NextL(L)  = 0.0
      RadPARFwdScat2NextL(L) = 0.0_r8

      D2500: DO NZ=1,NP
        RadDifSWLeafAbsorption_pft(NZ)  = 0.0_r8
        RadDifSWStemAbsorption_pft(NZ)  = 0.0_r8
        RadDifSWDStemAbsorption_pft(NZ) = 0.0_r8

        RadDifPARLeafAbsorption_pft(NZ)  = 0.0_r8
        RadDifPARStemAbsorption_pft(NZ)  = 0.0_r8
        RadDifPARDStemAbsorption_pft(NZ) = 0.0_R8

        D2600: DO N=1,NumLeafInclinationClasses1
          EffLeafArea     = LeafAreaZsecLive_lpft(N,L,NZ)*ClumpFactorNow_pft(NZ)
          EffStemArea     = StemAreaZsecLive_lpft(N,L,NZ)*StemClumpFactor_pft(NZ)
          EffDeadStemArea = SurfAreaZsecDead_lpft(N,L,NZ)*StemClumpFactor_pft(NZ)

          D2700: DO M = 1, NumOfSkyAzimuthSects1
            D2750: DO NN=1,NumOfLeafAzimuthSectors1
              DiffSWLeafAbsorptAzclass   = RadSWDiffusL*OMEGA2Leaf(M,N,NN)*LeafSWabsorptivityL_pft(L,NZ)
              DiffSWStemAbsorptAzclass  = RadSWDiffusL*OMEGA2Leaf(M,N,NN)*StemSWAbsorptivityL_pft(L,NZ)
              DiffSWDStemAbsorptAzclass = RadSWDiffusL*OMEGA2Leaf(M,N,NN)*DStemSWabsorptivityL_pft(L,NZ)

              DiffPARLeafAbsorptAzclass  = RadPARDiffusL*OMEGA2Leaf(M,N,NN)*LeafPARabsorptivityL_pft(L,NZ)
              DiffPARStemAbsorptAzclass = RadPARDiffusL*OMEGA2Leaf(M,N,NN)*StemPARAbsorptivityL_pft(L,NZ)
              DiffPARDStemAbsorptAzclass= RadPARDiffusL*OMEGA2Leaf(M,N,NN)*DStemPARabsorptivityL_pft(L,NZ)

              RadDifPARAbsorption_zsec(N,M,L,NZ) = RadDifPARAbsorption_zsec(N,M,L,NZ)+DiffPARLeafAbsorptAzclass
              RadTotPARAbsorption_zsec(N,M,L,NZ) = RadTotPARAbsorption_zsec(N,M,L,NZ)+DiffPARLeafAbsorptAzclass
              

              RadDifSWLeafAbsorption_pft(NZ) = RadDifSWLeafAbsorption_pft(NZ)+EffLeafArea*DiffSWLeafAbsorptAzclass
              RadDifSWStemAbsorption_pft(NZ) = RadDifSWStemAbsorption_pft(NZ)+EffStemArea*DiffSWStemAbsorptAzclass
              RadDifSWDStemAbsorption_pft(NZ)= RadDifSWDStemAbsorption_pft(NZ)+EffDeadStemArea*DiffSWDStemAbsorptAzclass

              RadDifPARLeafAbsorption_pft(NZ)  = RadDifPARLeafAbsorption_pft(NZ)+EffLeafArea*DiffPARLeafAbsorptAzclass
              RadDifPARStemAbsorption_pft(NZ) = RadDifPARStemAbsorption_pft(NZ)+EffStemArea*DiffPARStemAbsorptAzclass
              RadDifPARDStemAbsorption_pft(NZ)= RadDifPARDStemAbsorption_pft(NZ)+EffDeadStemArea*DiffPARDStemAbsorptAzclass
            ENDDO D2750
          ENDDO D2700
        ENDDO D2600

        RadSWFwdScat2NextL(L)  = RadSWFwdScat2NextL(L)+( &          
          RadDifSWLeafAbsorption_pft(NZ)*RadSWLeafTransmitanceL_pft(L,NZ) + &
          RadDifSWStemAbsorption_pft(NZ)*RadSWStemTransmitanceL_pft(L,NZ) + &
          RadDifSWDStemAbsorption_pft(NZ)*RadSWDstemTransmitanceL_pft(L,NZ))*YAREA

        RadPARFwdScat2NextL(L) = RadPARFwdScat2NextL(L)+( &
         RadDifPARLeafAbsorption_pft(NZ)*RadPARLeafTransmitanceL_pft(L,NZ) + &
         RadDifPARStemAbsorption_pft(NZ)*RadPARStemTransmitanceL_pft(L,NZ) + &
         RadDifPARDStemAbsorption_pft(NZ)*RadPARDstemTransmitanceL_pft(L,NZ))*YAREA

        RadSWCanopyAbsorption_pft(NZ)  = RadSWCanopyAbsorption_pft(NZ)+RadDifSWLeafAbsorption_pft(NZ) &
          +RadDifSWStemAbsorption_pft(NZ)+RadDifSWDStemAbsorption_pft(NZ)
        RadPARCanopyAbsorption_pft(NZ) = RadPARCanopyAbsorption_pft(NZ)+RadDifPARLeafAbsorption_pft(NZ) &
          +RadDifPARStemAbsorption_pft(NZ)+RadDifPARDStemAbsorption_pft(NZ)
        RadSWCanopyLAbsroption_pft(L,NZ)  = RadSWCanopyLAbsroption_pft(L,NZ)+RadDifSWLeafAbsorption_pft(NZ) &
          +RadDifSWStemAbsorption_pft(NZ)+RadDifSWDStemAbsorption_pft(NZ)
        RadPARCanopyLAbsorption_pft(L,NZ) =  RadPARCanopyLAbsorption_pft(L,NZ) +RadDifPARLeafAbsorption_pft(NZ) &
          +RadDifPARStemAbsorption_pft(NZ) + RadDifPARDStemAbsorption_pft(NZ)
      ENDDO D2500
    ELSE
      RadSWFwdScat2NextL(L)  = RadSWFwdScat2NextL(L-1)
      RadPARFwdScat2NextL(L) = RadPARFwdScat2NextL(L-1)
      RadSWBakScat2NextL(L)  = RadSWBakScat2NextL(L-1)
      RadPARBakScat2NextL(L) = RadPARBakScat2NextL(L-1)
    ENDIF
  ENDDO D2800  
  call PrintInfo('end '//subname)     
  end associate
  end subroutine UpwellSweep

!----------------------------------------------------------------------------------------------------
  subroutine AccumulateRadation4CanopyL(L,YAREA,RadDirPAR2Leaf_zsec,                                                    &
    FractionDirRadAbsorbtCum,FractionDifRadAbsorbtCum,FractionDirRadAbsorbtL,                                           &
    FractionDifRadAbsorbtL,RadDirSWLeafLAbsorption_pft,RadDirPARLeafLAbsorption_pft,                                    &
    RadDirSWStemLAbsorption_pft,RadDirPARStemLAbsorption_pft,RadDirSWDStemLAbsorption_pft,RadDirPARDStemLAbsorption_pft,&
    RadDifSWLeafAbsorption_pft,RadDifPARLeafAbsorption_pft,RadDifSWStemAbsorption_pft,RadDifPARStemAbsorption_pft,        &
    RadDifSWDStemAbsorption_pft,RadDifPARDStemAbsorption_pft,                                                             &
    backScatDirSWLeafLAbsorpt_pft,FwdScatDirSWLeafLAbsorpt_pft,backScatDirSWStemLAbsorpt_pft,FwdScatDirSWStemLAbsorpt_pft,&
    backScatDirSWDStemLAbsorpt_pft,FwdScatDirSWDStemLAbsorpt_pft,                                                         &
    backScatDirPARLeafLAbsorpt_pft,FwdScatDirPARLeafLAbsorpt_pft,backScatDirPARStemLAbsorpt_pft,FwdScatDirPARStemLAbsorpt_pft,&
    backScatDirPARDStemLAbsorpt_pft,FwdScatDirPARDStemLAbsorpt_pft,                                                     &
    BackScatDifSWLeafLAbsorpt_pft,FwdScatDifSWLeafLAbsorpt_pft,BackScatDifSWStemLAbsorpt_pft,FwdScatDifSWStemLAbsorpt_pft,&
    BackScatDifSWDStemLAbsorpt_pft,FwdScatDifSWDStemLAbsorpt_pft,                                                       &
    BackScatDifPARLeafLAbsorpt_pft,FwdScatDifPARLeafLAbsorpt_pft,BackScatDifPARStemLAbsorpt_pft,FwdScatDifPARStemLAbsorpt_pft,&
    BackScatDifPARDStemLAbsorpt_pft,FwdScatDifPARDStemLAbsorpt_pft,                                                     &
    RadSWBakScat2NextL,RadPARBakScat2NextL,RadSWFwdScat2NextL,RadPARFwdScat2NextL,                                     &
    TAU_RDiffusTransmitance)

  implicit none
  integer, intent(in)     :: L
  real(r8), intent(in)    :: YAREA 
  real(r8), intent(in)    :: RadDirPAR2Leaf_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)
  real(r8), intent(in)    :: FractionDirRadAbsorbtCum  
  REAL(R8), intent(in)    :: FractionDifRadAbsorbtCum  
  real(r8), intent(inout) :: FractionDirRadAbsorbtL
  REAL(R8), intent(inout) :: FractionDifRadAbsorbtL

  real(r8), intent(inout) :: RadDirSWLeafLAbsorption_pft(JP1), RadDirPARLeafLAbsorption_pft(JP1)
  real(r8), intent(inout) :: RadDirSWStemLAbsorption_pft(JP1), RadDirPARStemLAbsorption_pft(JP1)
  real(r8), intent(inout) :: RadDirSWDStemLAbsorption_pft(JP1),RadDirPARDStemLAbsorption_pft(JP1)

  real(r8), intent(inout) :: RadDifSWLeafAbsorption_pft(JP1),RadDifPARLeafAbsorption_pft(JP1)
  real(r8), intent(inout) :: RadDifSWStemAbsorption_pft(JP1),RadDifPARStemAbsorption_pft(JP1)
  real(r8), intent(inout) :: RadDifSWDStemAbsorption_pft(JP1),RadDifPARDStemAbsorption_pft(JP1)

  real(r8), intent(inout) :: backScatDirSWLeafLAbsorpt_pft(JP1),FwdScatDirSWLeafLAbsorpt_pft(JP1)
  real(r8), intent(inout) :: backScatDirSWStemLAbsorpt_pft(JP1),FwdScatDirSWStemLAbsorpt_pft(JP1)
  real(r8), intent(inout) :: backScatDirSWDStemLAbsorpt_pft(JP1),FwdScatDirSWDStemLAbsorpt_pft(JP1)

  real(r8), intent(inout) :: backScatDirPARLeafLAbsorpt_pft(JP1),FwdScatDirPARLeafLAbsorpt_pft(JP1)  
  real(r8), intent(inout) :: backScatDirPARStemLAbsorpt_pft(JP1),FwdScatDirPARStemLAbsorpt_pft(JP1)
  real(r8), intent(inout) :: backScatDirPARDStemLAbsorpt_pft(JP1),FwdScatDirPARDStemLAbsorpt_pft(JP1)

  real(r8), intent(inout) :: BackScatDifSWLeafLAbsorpt_pft(JP1),FwdScatDifSWLeafLAbsorpt_pft(JP1)
  real(r8), intent(inout) :: BackScatDifPARLeafLAbsorpt_pft(JP1),FwdScatDifPARLeafLAbsorpt_pft(JP1)

  real(r8), intent(inout) :: BackScatDifSWStemLAbsorpt_pft(JP1),FwdScatDifSWStemLAbsorpt_pft(JP1)
  real(r8), intent(inout) :: BackScatDifPARStemLAbsorpt_pft(JP1),FwdScatDifPARStemLAbsorpt_pft(JP1)

  real(r8), intent(inout) :: BackScatDifSWDStemLAbsorpt_pft(JP1),FwdScatDifSWDStemLAbsorpt_pft(JP1)
  real(r8), intent(inout) :: BackScatDifPARDStemLAbsorpt_pft(JP1),FwdScatDifPARDStemLAbsorpt_pft(JP1)

  real(r8), intent(inout) :: RadSWBakScat2NextL(0:NumCanopyLayers1+1)
  real(r8), intent(inout) :: RadPARBakScat2NextL(0:NumCanopyLayers1+1)  
  real(r8), intent(inout) :: RadSWFwdScat2NextL(0:NumCanopyLayers1+1)
  real(r8), intent(inout) :: RadPARFwdScat2NextL(0:NumCanopyLayers1+1)
  real(r8), intent(inout) :: TAU_RDiffusTransmitance(0:NumCanopyLayers1+1)

  character(len=*), parameter :: subname='AccumulateRadation4CanopyL'
  real(r8) :: RadDirSWLeafAbsorbT,RadDirSWStemAbsorbT,RadDirSWDStemAbsorbT
  real(r8) :: RadPARbyLeafT,RadPARbyStemT,RadPARbyDStemT
  real(r8) :: bakScatRadDirSWLeafAbsorbT,bakScatRadDirSWStemAbsorbT,bakScatRadDirSWDStemAbsorbT
  real(r8) :: bakScatRadPARbyLeafT,bakScatRadPARbyStemT,bakScatRadPARbyDStemT
  real(r8) :: fwdScatRadDirSWLeafAbsorbT,fwdScatRadDirSWStemAbsorbT,fwdScatRadDirSWDStemAbsorbT
  real(r8) :: fwdScatRadPARbyLeafT,fwdScatRadPARbyStemT,fwdScatRadPARbyDStemT
  real(r8) :: XTAUS,XTAUY
  integer  :: N,M,NZ

  associate(                                                            &
    NP                         => plt_site%NP                          ,& !input  :current number of plant species,[-]      
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]  
    TAU_DirectSunLit           => plt_rad%TAU_DirectSunLit             ,& !inoput :fraction of radiation intercepted by canopy layer, [-]  
    RadSWCanopyAbsorption_pft  => plt_rad%RadSWCanopyAbsorption_pft    ,& !inoput :canopy absorbed shortwave radiation, [MJ d-2 h-1]
    RadPARCanopyAbsorption_pft => plt_rad%RadPARCanopyAbsorption_pft   ,& !inoput :canopy absorbed PAR, [umol m-2 s-1]    
    RadSWCanopyLAbsroption_pft => plt_rad%RadSWCanopyLAbsroption_pft   ,& !output :profile of canopy absorbed shortwave radiation, [MJ d-2 h-1]         
    RadPARCanopyLAbsorption_pft=> plt_rad%RadPARCanopyLAbsorption_pft  ,& !output :profile of canopy absorbed PAR, [umol m-2 s-1]
    RadTotPARAbsorption_zsec   => plt_rad%RadTotPARAbsorption_zsec     ,& !inoput :total incoming PAR absorbed, [umol m-2 s-1]    
    RadDifPARAbsorption_zsec   => plt_rad%RadDifPARAbsorption_zsec      & !inoput :diffuse incoming PAR absorbed, [umol m-2 s-1]
  )
  !
  !     XTAUS=interception of direct radiation in current layer
  !     STOPZ=accumulated interception of direct radiation from topmost layer
  !     TAU_DirectSunLit=transmission of direct radiation to next lower layer
  !
  call PrintInfo('beg '//subname)

  !Eq. (14) in Grant et al. 1989
  IF(FractionDirRadAbsorbtCum+FractionDirRadAbsorbtL.GT.1.0_r8)THEN
    IF(FractionDirRadAbsorbtL.GT.ZERO)THEN
      !attenuate to deal with the over absorption situation in p81, Grant et al., 1989
      XTAUS=(1.0_r8-FractionDirRadAbsorbtCum)/FractionDirRadAbsorbtL
    ELSE
      XTAUS=0.0_r8
    ENDIF
    TAU_DirectSunLit(L+1)  = TAU_DirectSunLit(L+1)*XTAUS
    FractionDirRadAbsorbtL = FractionDirRadAbsorbtL*XTAUS
    
    D1510: DO NZ=1,NP
      RadDirSWLeafLAbsorption_pft(NZ)  = RadDirSWLeafLAbsorption_pft(NZ)*XTAUS
      RadDirSWStemLAbsorption_pft(NZ)  = RadDirSWStemLAbsorption_pft(NZ)*XTAUS
      RadDirSWDStemLAbsorption_pft(NZ)= RadDirSWDStemLAbsorption_pft(NZ)*XTAUS

      RadDirPARLeafLAbsorption_pft(NZ) = RadDirPARLeafLAbsorption_pft(NZ)*XTAUS
      RadDirPARStemLAbsorption_pft(NZ) = RadDirPARStemLAbsorption_pft(NZ)*XTAUS
      RadDirPARDStemLAbsorption_pft(NZ)= RadDirPARDStemLAbsorption_pft(NZ)*XTAUS

      backScatDirSWLeafLAbsorpt_pft(NZ)  = backScatDirSWLeafLAbsorpt_pft(NZ)*XTAUS
      backScatDirSWStemLAbsorpt_pft(NZ)  = backScatDirSWStemLAbsorpt_pft(NZ)*XTAUS
      backScatDirSWDStemLAbsorpt_pft(NZ) = backScatDirSWDStemLAbsorpt_pft(NZ)*XTAUS

      backScatDirPARLeafLAbsorpt_pft(NZ) = backScatDirPARLeafLAbsorpt_pft(NZ)*XTAUS
      backScatDirPARStemLAbsorpt_pft(NZ) = backScatDirPARStemLAbsorpt_pft(NZ)*XTAUS
      backScatDirPARDStemLAbsorpt_pft(NZ)= backScatDirPARDStemLAbsorpt_pft(NZ)*XTAUS

      FwdScatDirSWLeafLAbsorpt_pft(NZ)  = FwdScatDirSWLeafLAbsorpt_pft(NZ)*XTAUS
      FwdScatDirSWStemLAbsorpt_pft(NZ)  = FwdScatDirSWStemLAbsorpt_pft(NZ)*XTAUS
      FwdScatDirSWDStemLAbsorpt_pft(NZ) = FwdScatDirSWDStemLAbsorpt_pft(NZ)*XTAUS

      FwdScatDirPARLeafLAbsorpt_pft(NZ) = FwdScatDirPARLeafLAbsorpt_pft(NZ)*XTAUS
      FwdScatDirPARStemLAbsorpt_pft(NZ) = FwdScatDirPARStemLAbsorpt_pft(NZ)*XTAUS
      FwdScatDirPARDStemLAbsorpt_pft(NZ)= FwdScatDirPARDStemLAbsorpt_pft(NZ)*XTAUS
    ENDDO D1510
  ENDIF      
  !
  !     XTAUY=interception of diffuse radiation in current layer
  !     FracDifRadAbsorbt=accumulated interception of diffuse radiation from topmost layer
  !     TAU_RDiffusTransmitance=transmission of diffuse radiation to next lower layer
  ! Eq. (15) in Grant et al. 1989
  IF(FractionDifRadAbsorbtCum+FractionDifRadAbsorbtL.GT.1.0_r8)THEN
    XTAUY=(1.0_r8-FractionDifRadAbsorbtCum)/FractionDifRadAbsorbtL

    TAU_RDiffusTransmitance(L+1) = TAU_RDiffusTransmitance(L+1)*XTAUY
    FractionDifRadAbsorbtL       = FractionDifRadAbsorbtL*XTAUY

    D1520: DO NZ=1,NP
      RadDifSWLeafAbsorption_pft(NZ)  = RadDifSWLeafAbsorption_pft(NZ)*XTAUY
      RadDifSWStemAbsorption_pft(NZ)  = RadDifSWStemAbsorption_pft(NZ)*XTAUY
      RadDifSWDStemAbsorption_pft(NZ) = RadDifSWDStemAbsorption_pft(NZ)*XTAUY

      RadDifPARLeafAbsorption_pft(NZ) = RadDifPARLeafAbsorption_pft(NZ)*XTAUY
      RadDifPARStemAbsorption_pft(NZ) = RadDifPARStemAbsorption_pft(NZ)*XTAUY
      RadDifPARDStemAbsorption_pft(NZ)= RadDifPARDStemAbsorption_pft(NZ)*XTAUY

      BackScatDifSWLeafLAbsorpt_pft(NZ)  = BackScatDifSWLeafLAbsorpt_pft(NZ)*XTAUY
      BackScatDifSWStemLAbsorpt_pft(NZ)  = BackScatDifSWStemLAbsorpt_pft(NZ)*XTAUY
      BackScatDifSWDStemLAbsorpt_pft(NZ) = BackScatDifSWDStemLAbsorpt_pft(NZ)*XTAUY

      BackScatDifPARLeafLAbsorpt_pft(NZ) = BackScatDifPARLeafLAbsorpt_pft(NZ)*XTAUY
      BackScatDifPARStemLAbsorpt_pft(NZ) = BackScatDifPARStemLAbsorpt_pft(NZ)*XTAUY
      BackScatDifPARDStemLAbsorpt_pft(NZ)= BackScatDifPARDStemLAbsorpt_pft(NZ)*XTAUY

      FwdScatDifSWLeafLAbsorpt_pft(NZ)  = FwdScatDifSWLeafLAbsorpt_pft(NZ)*XTAUY
      FwdScatDifSWStemLAbsorpt_pft(NZ)  = FwdScatDifSWStemLAbsorpt_pft(NZ)*XTAUY
      FwdScatDifSWDStemLAbsorpt_pft(NZ) = FwdScatDifSWDStemLAbsorpt_pft(NZ)*XTAUY 

      FwdScatDifPARLeafLAbsorpt_pft(NZ)   = FwdScatDifPARLeafLAbsorpt_pft(NZ)*XTAUY
      FwdScatDifPARStemLAbsorpt_pft(NZ) = FwdScatDifPARStemLAbsorpt_pft(NZ)*XTAUY
      FwdScatDifPARDStemLAbsorpt_pft(NZ) = FwdScatDifPARDStemLAbsorpt_pft(NZ)*XTAUY 

      D1730: DO N=1,NumLeafInclinationClasses1
        DO  M=1,NumOfSkyAzimuthSects1
          RadDifPARAbsorption_zsec(N,M,L,NZ) = RadDifPARAbsorption_zsec(N,M,L,NZ)*XTAUY
          RadTotPARAbsorption_zsec(N,M,L,NZ) = RadDirPAR2Leaf_zsec(N,M,NZ)*LeafPARabsorptivityL_pft(L,NZ)+RadDifPARAbsorption_zsec(N,M,L,NZ)
        enddo
      ENDDO D1730
    ENDDO D1520
  ENDIF

  !
  ! Summarize radiation in layer L
  !     TOTAL RADIATION ABSORBED, REFLECTED AND TRANSMITTED BY ALL PFTs
  !
  !     FractionDirRadAbsorbtCum,FractionDifRadAbsorbtCum=accumulated interception of direct,diffuse radiation
  !     TAU_DirectSunLit,TAU_RDiffusTransmitance=transmission of direct,diffuse radiation to next lower layer
  !
  D1530: DO NZ=1,NP

    RadDirSWLeafAbsorbT = RadDirSWLeafLAbsorption_pft(NZ)+RadDifSWLeafAbsorption_pft(NZ)
    RadDirSWStemAbsorbT = RadDirSWStemLAbsorption_pft(NZ)+RadDifSWStemAbsorption_pft(NZ)
    RadDirSWDStemAbsorbT= RadDirSWDStemLAbsorption_pft(NZ)+RadDifSWDStemAbsorption_pft(NZ)

    RadPARbyLeafT  = RadDirPARLeafLAbsorption_pft(NZ)+RadDifPARLeafAbsorption_pft(NZ)
    RadPARbyStemT  = RadDirPARStemLAbsorption_pft(NZ)+RadDifPARStemAbsorption_pft(NZ)
    RadPARbyDStemT = RadDirPARDStemLAbsorption_pft(NZ)+RadDifPARDStemAbsorption_pft(NZ)

    bakScatRadDirSWLeafAbsorbT   = backScatDirSWLeafLAbsorpt_pft(NZ)+BackScatDifSWLeafLAbsorpt_pft(NZ)
    bakScatRadDirSWStemAbsorbT   = backScatDirSWStemLAbsorpt_pft(NZ)+BackScatDifSWStemLAbsorpt_pft(NZ)
    bakScatRadDirSWDstemAbsorbT  = backScatDirSWDStemLAbsorpt_pft(NZ)+BackScatDifSWDStemLAbsorpt_pft(NZ)

    bakScatRadPARbyLeafT  = backScatDirPARLeafLAbsorpt_pft(NZ)+BackScatDifPARLeafLAbsorpt_pft(NZ)
    bakScatRadPARbyStemT  = backScatDirPARStemLAbsorpt_pft(NZ)+BackScatDifPARStemLAbsorpt_pft(NZ)
    bakScatRadPARbyDStemT = backScatDirPARDStemLAbsorpt_pft(NZ)+BackScatDifPARDStemLAbsorpt_pft(NZ)

    fwdScatRadDirSWLeafAbsorbT  = FwdScatDirSWLeafLAbsorpt_pft(NZ)+FwdScatDifSWLeafLAbsorpt_pft(NZ)
    fwdScatRadDirSWStemAbsorbT  = FwdScatDirSWStemLAbsorpt_pft(NZ)+FwdScatDifSWStemLAbsorpt_pft(NZ)
    fwdScatRadDirSWDStemAbsorbT = FwdScatDirSWDStemLAbsorpt_pft(NZ)+FwdScatDifSWDStemLAbsorpt_pft(NZ)

    fwdScatRadPARbyLeafT  = FwdScatDirPARLeafLAbsorpt_pft(NZ)+FwdScatDifPARLeafLAbsorpt_pft(NZ)
    fwdScatRadPARbyStemT  = FwdScatDirPARStemLAbsorpt_pft(NZ)+FwdScatDifPARStemLAbsorpt_pft(NZ)
    fwdScatRadPARbyDStemT = FwdScatDirPARDStemLAbsorpt_pft(NZ)+FwdScatDifPARDStemLAbsorpt_pft(NZ)

    RadSWFwdScat2NextL(L)=RadSWFwdScat2NextL(L)+ ( &
       RadDirSWLeafAbsorbT*RadSWLeafTransmitanceL_pft(L,NZ) + &
       RadDirSWStemAbsorbT*RadSWStemTransmitanceL_pft(L,NZ) + &
       RadDirSWDStemAbsorbT*RadSWDstemTransmitanceL_pft(L,NZ) + &
       fwdScatRadDirSWLeafAbsorbT*RadSWLeafAlbedoL_pft(L,NZ) + &
       fwdScatRadDirSWStemAbsorbT*RadSWStemAlbedoL_pft(L,NZ) + &
       fwdScatRadDirSWDStemAbsorbT*RadSWDStemAlbedoL_pft(L,NZ))*YAREA

    RadPARFwdScat2NextL(L)=RadPARFwdScat2NextL(L)+( &
      RadPARbyLeafT*RadPARLeafTransmitanceL_pft(L,NZ) +&
      RadPARbyStemT*RadPARStemTransmitanceL_pft(L,NZ) + &
      RadPARbyDStemT*RadPARDstemTransmitanceL_pft(L,NZ) + &
      fwdScatRadPARbyLeafT*RadPARLeafAlbedoL_pft(L,NZ) + &
      fwdScatRadPARbyStemT*RadPARStemAlbedoL_pft(L,NZ) + &
      fwdScatRadPARbyDStemT*RadPARDstemAlbedoL_pft(L,NZ))*YAREA

    RadSWBakScat2NextL(L)=RadSWBakScat2NextL(L)+ &
      (bakScatRadDirSWLeafAbsorbT*RadSWLeafAlbedoL_pft(L,NZ)+ &
       bakScatRadDirSWStemAbsorbT*RadSWStemAlbedoL_pft(L,NZ) + &
       bakScatRadDirSWDstemAbsorbT*RadSWDStemAlbedoL_pft(L,NZ))*YAREA

    RadPARBakScat2NextL(L)=RadPARBakScat2NextL(L)+&
      (bakScatRadPARbyLeafT*RadPARLeafAlbedoL_pft(L,NZ)+&
       bakScatRadPARbyStemT*RadPARStemAlbedoL_pft(L,NZ)+&
       bakScatRadPARbyDStemT*RadPARDstemAlbedoL_pft(L,NZ))*YAREA

    !accumulate shortwave radiation on canopy

    RadSWCanopyAbsorption_pft(NZ)    = RadSWCanopyAbsorption_pft(NZ)+RadDirSWLeafAbsorbT+RadDirSWStemAbsorbT+RadDirSWDStemAbsorbT
    RadSWCanopyLAbsroption_pft(L,NZ) = RadDirSWLeafAbsorbT+RadDirSWStemAbsorbT+RadDirSWDStemAbsorbT
    RadPARCanopyAbsorption_pft(NZ)   = RadPARCanopyAbsorption_pft(NZ)+RadPARbyLeafT+RadPARbyStemT+RadPARbyDStemT
    RadPARCanopyLAbsorption_pft(L,NZ) = RadPARbyLeafT+RadPARbyStemT+RadPARbyDStemT
  ENDDO D1530

  call PrintInfo('end '//subname)
  end associate    
  end subroutine AccumulateRadation4CanopyL

!----------------------------------------------------------------------------------------------------
  subroutine RadiationAtGround(cosGroundIncidentSolarAngle,       &
    RadSWDiffusL,RadPARDiffusL,                                   &
    TAU_DirectSunLit1,TAU_RDiffusTransmitance1,                   &
    RadSWFwdScat2NextL1,RadPARFwdScat2NextL1                     ,&
    RadSWBakScat2NextL0,RadPARBakScat2NextL0,                     &
    RadSW2Ground,RadPAR2Ground)
  implicit none
  real(r8), intent(in)  :: cosGroundIncidentSolarAngle          !
  real(r8), intent(in)  :: RadSWDiffusL,RadPARDiffusL        !diffuse radiation and PAR radiation off the canopy
  real(r8), intent(in)  :: TAU_DirectSunLit1,TAU_RDiffusTransmitance1 !cumulative transmittance off the canopy
  real(r8), intent(in)  :: RadSWFwdScat2NextL1,RadPARFwdScat2NextL1          !incoming radiation off the canopy
  real(r8), intent(out) :: RadSWBakScat2NextL0,RadPARBakScat2NextL0          !reflected radiation at the ground
  real(r8), intent(out) :: RadSW2Ground          !incident radiation onto ground
  real(r8), intent(out) :: RadPAR2Ground         !incident PAR radiation onto ground

  character(len=*), parameter :: subname='RadiationAtGround'
  real(r8) :: RADSG,RADYG  
  real(r8) :: RAPSG,RAPYG
  real(r8) :: THETW1
  real(r8) :: SnowpackAlbedo,GrndAlbedo
  real(r8) :: FracGrndBySnow
  integer  :: N

  associate(                                                    &
    SnowDepth              => plt_ew%SnowDepth                 ,& !input  :snowpack depth, [m]  
    SurfAlbedo_col         => plt_rad%SurfAlbedo_col           ,& !input  :Surface albedo,[-]    
    SoilAlbedo             => plt_rad%SoilAlbedo               ,& !input  :soil albedo,[-]    
    OMEGA2Ground           => plt_rad%OMEGA2Ground             ,& !input  :cosine of solar beam onto ground, [-]  
    VcumIceSnow_col        => plt_ew%VcumIceSnow_col           ,& !input  :ice volume in snowpack, [m3 d-2]    
    VcumWatSnow_col        => plt_ew%VcumWatSnow_col           ,& !input  :water volume in snowpack, [m3 d-2]    
    RadDirectPAR_col       => plt_rad%RadDirectPAR_col         ,& !input  :direct PAR, [umol m-2 s-1]        
    RadSWDirect_col        => plt_rad%RadSWDirect_col          ,& !input  :direct shortwave radiation, [W m-2]
    NU                     => plt_site%NU                      ,& !input  :current soil surface layer number, [-]    
    VLHeatCapSnowMin_col   => plt_ew%VLHeatCapSnowMin_col      ,& !input  :minimum snowpack heat capacity, [MJ d-2 K-1]    
    VcumDrySnoWE_col       => plt_ew%VcumDrySnoWE_col          ,& !input  :snow volume in snowpack (water equivalent), [m3 d-2]    
    VLHeatCapSurfSnow_col  => plt_ew%VLHeatCapSurfSnow_col     ,& !input  :snowpack heat capacity, [MJ m-3 K-1]    
    ZEROS2                 => plt_site%ZEROS2                  ,& !input  :threshold zero for numerical stability,[-]    
    VLSoilMicP_vr          => plt_soilchem%VLSoilMicP_vr       ,& !input  :total micropore volume in layer, [m3 d-2]
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr        ,& !input  :soil micropore water content, [m3 d-2]
    POROS1                 => plt_site%POROS1                  ,& !input  :top layer soil porosity, [m3 m-3]    
    VLSoilPoreMicP_vr      => plt_soilchem%VLSoilPoreMicP_vr    & !input  :volume of soil layer, [m3 d-2]    
  )
  call PrintInfo('beg '//subname)
  !
  !     DIRECT AND DIFFUSE RADIATION ABSORBED AT GROUND SURFACE
  !
  !     RADSG,RADYG,RAPSG,RAPYG=direct,diffuse SW,PAR at horiCanopyHeightZ_col ground surface
  !     RADS,RadDirectPAR_col=solar beam direct SW,PAR flux
  !     TAU_DirectSunLit,TAU_RDiffusTransmitance=transmission of direct,diffuse radiation below canopy
  !     RadSWDiffusL,RadPARDiffusL=solar beam diffuse SW,PAR flux
  !     RadSW2Ground,RadPAR2Ground=total SW,PAR at ground surface
  !     cosGroundIncidentSolarAngle,OMEGA2Ground=incident solar,sky angle at ground surface
  !
  RADSG = RadSWDirect_col*TAU_DirectSunLit1
  RADYG = RadSWDiffusL*TAU_RDiffusTransmitance1+RadSWFwdScat2NextL1
  RAPSG = RadDirectPAR_col*TAU_DirectSunLit1
  RAPYG = RadPARDiffusL*TAU_RDiffusTransmitance1+RadPARFwdScat2NextL1

  RadSW2Ground  = ABS(cosGroundIncidentSolarAngle)*RADSG
  RadPAR2Ground = ABS(cosGroundIncidentSolarAngle)*RAPSG

  !add diffuse radiation from different sky azimuthal region
  D20: DO N=1,NumOfSkyAzimuthSects1
    RadSW2Ground  = RadSW2Ground+ABS(OMEGA2Ground(N))*RADYG
    RadPAR2Ground = RadPAR2Ground+ABS(OMEGA2Ground(N))*RAPYG
  ENDDO D20
  !
  !     RADIATION REFLECTED FROM GROUND SURFACE
  !
  !     VHCPW,VLHeatCapSnowMin_col=current,minimum snowpack heat capacity
  !     SnowpackAlbedo,VcumDrySnoWE_col,VcumWatSnow_col,VcumIceSnow_col=snowpack albedo,snow,water,ice volume
  !     GrndAlbedo,SoilAlbedo,FracGrndBySnow=ground,soil albedo,snow cover fraction
  !     THETW1=soil surface water content
  !     RadSWBakScat2NextL,RadDirPARLeafLAbsorption_pft=SW,PAR backscatter from ground surface
  !     TRADG,TRAPG=SW,PAR absorbed by ground surface
  !
  IF(VLHeatCapSurfSnow_col.GT.VLHeatCapSnowMin_col)THEN
    SnowpackAlbedo=(0.85_r8*VcumDrySnoWE_col+0.30_r8*VcumIceSnow_col+0.06_r8*VcumWatSnow_col) &
      /(VcumDrySnoWE_col+VcumIceSnow_col+VcumWatSnow_col)
    !the following partition differs from that used in the surface physics module
    FracGrndBySnow = AMIN1((SnowDepth/0.07_r8)**2._r8,1.0_r8)
    GrndAlbedo     = FracGrndBySnow*SnowpackAlbedo+(1.0_r8-FracGrndBySnow)*SoilAlbedo
  ELSE
    IF(VLSoilPoreMicP_vr(NU).GT.ZEROS2)THEN
      THETW1=AMIN1(POROS1,VLWatMicP_vr(NU)/VLSoilMicP_vr(NU))
    ELSE
      THETW1=0.0_r8
    ENDIF
    GrndAlbedo=AMIN1(SurfAlbedo_col,SoilAlbedo+AZMAX1(SurfAlbedo_col-THETW1))
  ENDIF
  !reflection on ground
  RadSWBakScat2NextL0  = RadSW2Ground*GrndAlbedo/real(NumOfSkyAzimuthSects1,r8)
  RadPARBakScat2NextL0 = RadPAR2Ground*GrndAlbedo/real(NumOfSkyAzimuthSects1,r8)

  call PrintInfo('end '//subname)
  end associate
  end subroutine RadiationAtGround
!----------------------------------------------------------------------------------------------------
  subroutine DownwellRadiationCanopyZL(NZ,L,XAREA,YAREA,BETX,iScatteringDirect,                   &
    LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,SurfAreaZsecDead_lpft,                            &
    RadSWDiffusL,RadPARDiffusL,FractionDirRadAbsorbtL,FractionDifRadAbsorbtL,                     &
    RadDirSW2Leaf_zsec,RadDirSW2Stem_zsec,RadDirSW2DStem_zsec,                                    &
    RadDirPAR2Leaf_zsec, RadDirPAR2Stem_zsec,RadDirPAR2DStem_zsec,                                &
    RadDirSWLeafLAbsorption_pft,backScatDirSWLeafLAbsorpt_pft,FwdScatDirSWLeafLAbsorpt_pft,       &
    RadDirSWStemLAbsorption_pft,backScatDirSWStemLAbsorpt_pft,FwdScatDirSWStemLAbsorpt_pft,       &
    RadDirSWDStemLAbsorption_pft,backScatDirSWDStemLAbsorpt_pft,FwdScatDirSWDStemLAbsorpt_pft,    &    
    RadDirPARLeafLAbsorption_pft,backScatDirPARLeafLAbsorpt_pft,FwdScatDirPARLeafLAbsorpt_pft,    &
    RadDirPARStemLAbsorption_pft,backScatDirPARStemLAbsorpt_pft,FwdScatDirPARStemLAbsorpt_pft,    &
    RadDirPARDStemLAbsorption_pft,backScatDirPARDStemLAbsorpt_pft,FwdScatDirPARDStemLAbsorpt_pft, &
    RadDifSWLeafAbsorption_pft,BackScatDifSWLeafLAbsorpt_pft,FwdScatDifSWLeafLAbsorpt_pft,        &
    RadDifSWStemAbsorption_pft,BackScatDifSWStemLAbsorpt_pft,FwdScatDifSWStemLAbsorpt_pft,        &
    RadDifSWDStemAbsorption_pft,BackScatDifSWDStemLAbsorpt_pft,FwdScatDifSWDStemLAbsorpt_pft,     &
    RadDifPARLeafAbsorption_pft,BackScatDifPARLeafLAbsorpt_pft,FwdScatDifPARLeafLAbsorpt_pft,     &
    RadDifPARStemAbsorption_pft,BackScatDifPARStemLAbsorpt_pft,FwdScatDifPARStemLAbsorpt_pft,     &
    RadDifPARDStemAbsorption_pft,BackScatDifPARDStemLAbsorpt_pft,FwdScatDifPARDStemLAbsorpt_pft)

  implicit none
  integer, intent(in) :: NZ,L
  real(r8), intent(in) :: XAREA,YAREA
  real(r8), INTENT(IN) :: BETX(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1)                         !sine of direct solar radiation on leaf surface/sine of direct solar radiation, [-]  
  integer,  INTENT(IN) :: iScatteringDirect(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1)  
  real(r8), intent(in) :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses1)
  real(r8), intent(in) :: StemAreaZsecLive_lpft(NumLeafInclinationClasses1)  
  real(r8), intent(in) :: SurfAreaZsecDead_lpft(NumLeafInclinationClasses1)          
  real(r8), intent(in) :: RadSWDiffusL,RadPARDiffusL  !incoming radiation to layer L
  REAL(R8), intent(in) :: RadDirSW2Leaf_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)    !direct shortwave radiation on leave in each canopy sector
  real(r8), intent(in) :: RadDirSW2Stem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)    !direct shortwave radiation on stalk in each canopy sector
  real(r8), intent(in) :: RadDirSW2DStem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)   !direct shortwave radiation on dead stem in each canopy sector  
  real(r8), intent(in) :: RadDirPAR2Leaf_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)   !PAR radiation on leave in each canopy sector
  real(r8), intent(in) :: RadDirPAR2Stem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)   !PAR radiation on stalk in each canopy sector
  real(r8), intent(in) :: RadDirPAR2DStem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)  !PAR radiation on dead stem in each canopy sector
  real(r8), intent(inout):: FractionDirRadAbsorbtL,FractionDifRadAbsorbtL
  real(r8), intent(out) :: RadDirSWLeafLAbsorption_pft,backScatDirSWLeafLAbsorpt_pft,FwdScatDirSWLeafLAbsorpt_pft
  real(r8), intent(out) :: RadDirSWStemLAbsorption_pft,backScatDirSWStemLAbsorpt_pft,FwdScatDirSWStemLAbsorpt_pft
  real(r8), intent(out) :: RadDirPARLeafLAbsorption_pft,backScatDirPARLeafLAbsorpt_pft,FwdScatDirPARLeafLAbsorpt_pft
  real(r8), intent(out) :: RadDirPARStemLAbsorption_pft,backScatDirPARStemLAbsorpt_pft,FwdScatDirPARStemLAbsorpt_pft
  real(r8), intent(out) :: RadDifSWLeafAbsorption_pft,BackScatDifSWLeafLAbsorpt_pft,FwdScatDifSWLeafLAbsorpt_pft
  real(r8), intent(out) :: RadDifSWStemAbsorption_pft,BackScatDifSWStemLAbsorpt_pft,FwdScatDifSWStemLAbsorpt_pft
  real(r8), intent(out) :: RadDifPARLeafAbsorption_pft,BackScatDifPARLeafLAbsorpt_pft,FwdScatDifPARLeafLAbsorpt_pft
  real(r8), intent(out) :: RadDifPARStemAbsorption_pft,BackScatDifPARStemLAbsorpt_pft,FwdScatDifPARStemLAbsorpt_pft
  real(r8), intent(out) :: RadDirSWDStemLAbsorption_pft,backScatDirSWDStemLAbsorpt_pft,FwdScatDirSWDStemLAbsorpt_pft
  real(r8), intent(out) :: RadDirPARDStemLAbsorption_pft,backScatDirPARDStemLAbsorpt_pft,FwdScatDirPARDStemLAbsorpt_pft
  real(r8), intent(out) :: RadDifSWDStemAbsorption_pft,BackScatDifSWDStemLAbsorpt_pft,FwdScatDifSWDStemLAbsorpt_pft
  real(r8), intent(out) :: RadDifPARDStemAbsorption_pft,BackScatDifPARDStemLAbsorpt_pft,FwdScatDifPARDStemLAbsorpt_pft

  CHARACTER(LEN=*), parameter :: subname='DownwellRadiationCanopyZL'
  integer :: N,M,NN
  REAL(R8) :: EffLeafArea,EffLeafAreaAzclass
  real(r8) :: LeafArea4Interception,SulitLeafEffArea
  real(r8) :: EffStemArea,EffStemAreaAzclass,SunlitStemEffArea,StemArea4Intception
  real(r8) :: DiffSWStemAbsorptAzclass,DiffSWLeafAbsorptAzclass,DiffSWDStemAbsorptAzclass
  real(r8) :: DiffPARLeafAbsorptAzclass,DiffPARStemAbsorptAzclass,DiffPARDStemAbsorptAzclass
  real(r8) :: RadDiffSWLeafAbsorbAzc,RadDiffSWStemAbsorbAzc,RadDiffSWDStemAbsorbAzc
  real(r8) :: RadDiffPARLeafAbsorbAzc,RadDiffPARStemAbsorbAzc,RadDiffPARDStemAbsorbAzc
  real(r8) :: DirSWLeafAbsorptionAzclass,DirSWStemAbsorptionAzclass
  real(r8) :: DirPARLeafAbsorptionAzclass,DirPARStemAbsorptionAzclass  
  real(r8) :: EffDeadStemArea,EffDeadStemAreaAzclass,SunlitDeadStemEffArea,DeadStemArea4Intception
  real(r8) :: DirSWDeadStemAbsorpAzclass,DirPARDeadStemAbsorpAzclass
  real(r8) :: RadSWDiffsLMNN,RadPARDiffsLMNN

  associate(                                                                &
    OMEGX                          => plt_rad%OMEGX                        ,& !input  :cosine of indirect sky radiation on leaf surface/sine of indirect sky radiation,[-]  
    OMEGA2Leaf                     => plt_rad%OMEGA2Leaf                   ,& !input  :cosine of indirect sky radiation on leaf surface,[-]  
    iScatteringDiffus              => plt_rad%iScatteringDiffus            ,& !input  :flag for calculating backscattering of radiation in canopy,[-]    
    LeafPARabsorptivity_pft        => plt_rad%LeafPARabsorptivity_pft      ,& !input  :canopy PAR absorptivity,[-]        
    LeafSWabsorptivity_pft         => plt_rad%LeafSWabsorptivity_pft       ,& !input  :canopy shortwave absorptivity, [-]    
    iPlant2ndGrothPattern_pft      => plt_pheno%iPlant2ndGrothPattern_pft  ,& !input  :plant expression of secondary growth, [-]    
    RadTotPARAbsorption_zsec       => plt_rad%RadTotPARAbsorption_zsec     ,& !inoput :total incoming PAR, [umol m-2 s-1]    
    RadDifPARAbsorption_zsec       => plt_rad%RadDifPARAbsorption_zsec     ,& !inoput :diffuse PAR absorption, [umol m-2 s-1]    
    ClumpFactorNow_pft             => plt_morph%ClumpFactorNow_pft         ,& !input  :clumping factor for self-shading in canopy layer at current LAI, [-]
    TAU_DirectSunLit               => plt_rad%TAU_DirectSunLit              & !inoput :fraction of radiation intercepted by canopy layer, [-]    
  )
  call PrintInfo('beg '//subname)

  RadDirSWLeafLAbsorption_pft  = 0._r8;BackScatDirSWLeafLAbsorpt_pft  = 0._r8;FwdScatDirSWLeafLAbsorpt_pft  = 0._r8
  RadDirSWStemLAbsorption_pft  = 0._r8;BackScatDirSWStemLAbsorpt_pft  = 0._r8;FwdScatDirSWStemLAbsorpt_pft  = 0._r8
  RadDirSWDStemLAbsorption_pft = 0._r8;BackScatDirSWDStemLAbsorpt_pft = 0._r8;FwdScatDirSWDStemLAbsorpt_pft = 0._r8

  RadDirPARLeafLAbsorption_pft = 0._r8;BackScatDirPARLeafLAbsorpt_pft = 0._r8;FwdScatDirPARLeafLAbsorpt_pft = 0._r8
  RadDirPARStemLAbsorption_pft = 0._r8;BackScatDirPARStemLAbsorpt_pft = 0._r8;FwdScatDirPARStemLAbsorpt_pft = 0._r8
  RadDirPARDStemLAbsorption_pft= 0._r8;BackScatDirPARDStemLAbsorpt_pft= 0._r8;FwdScatDirPARDStemLAbsorpt_pft= 0._r8

  RadDifSWLeafAbsorption_pft   = 0._r8;BackScatDifSWLeafLAbsorpt_pft  = 0._r8;FwdScatDifSWLeafLAbsorpt_pft  = 0._r8
  RadDifSWStemAbsorption_pft   = 0._r8;BackScatDifSWStemLAbsorpt_pft  = 0._r8;FwdScatDifSWStemLAbsorpt_pft  = 0._r8
  RadDifSWDStemAbsorption_pft  = 0._r8;BackScatDifSWDStemLAbsorpt_pft = 0._r8;FwdScatDifSWDStemLAbsorpt_pft = 0._r8

  RadDifPARLeafAbsorption_pft  = 0._r8;BackScatDifPARLeafLAbsorpt_pft = 0._r8;FwdScatDifPARLeafLAbsorpt_pft = 0._r8
  RadDifPARStemAbsorption_pft  = 0._r8;BackScatDifPARStemLAbsorpt_pft = 0._r8;FwdScatDifPARStemLAbsorpt_pft = 0._r8
  RadDifPARDStemAbsorption_pft = 0._r8;BackScatDifPARDStemLAbsorpt_pft= 0._r8;FwdScatDifPARDStemLAbsorpt_pft= 0._r8 
  !
  !     LEAF SURFACE AREA IN EACH INCLINATION CLASS N, AZIMUTH CLASS M,
  !     LAYER L AND SPECIES NZ
  !
  !     EffLeafArea=unself-shaded leaf area
  !     EffLeafAreaAzclass=unself-shaded leaf area m-2 in each azimuth class
  !     SulitLeafEffArea=EffLeafArea with shading from canopy layers above
  !     LeafArea4Interception=SulitLeafEffArea m-2
  !     EffStemArea=unself-shaded stalk area
  !     EffStemAreaAzclass=unself-shaded stalk area m-2 in each azimuth class
  !     SunlitStemEffArea=EffStemArea with shading from canopy layers above
  !     StemArea4Intception=SunlitStemEffArea m-2
  !     XAREA = 1.0_r8/AREA3(NU)                      !average over the grid
  !     YAREA = XAREA/REAL(NumOfSkyAzimuthSects1,R8)  !area for one azimuthal zone  
  !
  D1600: DO N=1,NumLeafInclinationClasses1
    !          
    EffLeafArea           = LeafAreaZsecLive_lpft(N)*ClumpFactorNow_pft(NZ)
    EffLeafAreaAzclass    = EffLeafArea*YAREA
    SulitLeafEffArea      = EffLeafArea*TAU_DirectSunLit(L+1)
    LeafArea4Interception = SulitLeafEffArea*XAREA

    EffStemArea         = StemAreaZsecLive_lpft(N)*StemClumpFactor_pft(NZ)
    EffStemAreaAzclass  = EffStemArea*YAREA    
    SunlitStemEffArea   = EffStemArea*TAU_DirectSunLit(L+1)
    StemArea4Intception = SunlitStemEffArea*XAREA

    EffDeadStemArea         = SurfAreaZsecDead_lpft(N)*StemClumpFactor_pft(NZ)
    EffDeadStemAreaAzclass  = EffDeadStemArea*YAREA
    SunlitDeadStemEffArea   = EffDeadStemArea*TAU_DirectSunLit(L+1)
    DeadStemArea4Intception = SunlitDeadStemEffArea*XAREA
    !
    !     ABSORPTION OF DIRECT RADIATION BY SUNLIT LEAF SURFACES
    !
    !     STOPZ=accumulated horizontal area of intercepted direct radiation
    ! assuming azimuthally uniform distribution of leaf and stalk distribution
    D1700: DO M=1,NumOfSkyAzimuthSects1
      
      DirSWLeafAbsorptionAzclass  = SulitLeafEffArea*RadDirSW2Leaf_zsec(N,M,NZ)*LeafSWabsorptivityL_pft(L,NZ)
      DirSWStemAbsorptionAzclass  = SunlitStemEffArea*RadDirSW2Stem_zsec(N,M,NZ)*StemSWAbsorptivityL_pft(L,NZ)
      DirSWDeadStemAbsorpAzclass  = SunlitDeadStemEffArea*RadDirSW2DStem_zsec(N,M,NZ)*DStemSWabsorptivityL_pft(L,NZ)

      RadDirSWLeafLAbsorption_pft   = RadDirSWLeafLAbsorption_pft+DirSWLeafAbsorptionAzclass
      RadDirSWStemLAbsorption_pft   = RadDirSWStemLAbsorption_pft+DirSWStemAbsorptionAzclass
      RadDirSWDStemLAbsorption_pft = RadDirSWDStemLAbsorption_pft+DirSWDeadStemAbsorpAzclass

      DirPARLeafAbsorptionAzclass  = SulitLeafEffArea*RadDirPAR2Leaf_zsec(N,M,NZ)*LeafPARabsorptivityL_pft(L,NZ)
      DirPARStemAbsorptionAzclass  = SunlitStemEffArea*RadDirPAR2Stem_zsec(N,M,NZ)*StemPARAbsorptivityL_pft(L,NZ)
      DirPARDeadStemAbsorpAzclass  = SunlitDeadStemEffArea*RadDirPAR2DStem_zsec(N,M,NZ)*DStemPARabsorptivityL_pft(L,NZ)

      RadDirPARLeafLAbsorption_pft = RadDirPARLeafLAbsorption_pft+DirPARLeafAbsorptionAzclass
      RadDirPARStemLAbsorption_pft = RadDirPARStemLAbsorption_pft+DirPARStemAbsorptionAzclass
      RadDirPARDStemLAbsorption_pft = RadDirPARDStemLAbsorption_pft+DirPARDeadStemAbsorpAzclass
      !eq. (14)
      FractionDirRadAbsorbtL    = FractionDirRadAbsorbtL+(LeafArea4Interception+StemArea4Intception+DeadStemArea4Intception)*BETX(N,M)
      !
      !     BACKSCATTERING OF REFLECTED DIRECT RADIATION
      !
      IF(iScatteringDirect(N,M).EQ.ibackward)THEN
        backScatDirSWLeafLAbsorpt_pft = backScatDirSWLeafLAbsorpt_pft+DirSWLeafAbsorptionAzclass
        backScatDirSWStemLAbsorpt_pft = backScatDirSWStemLAbsorpt_pft+DirSWStemAbsorptionAzclass
        backScatDirSWDStemLAbsorpt_pft= backScatDirSWDStemLAbsorpt_pft+DirSWDeadStemAbsorpAzclass

        backScatDirPARLeafLAbsorpt_pft = backScatDirPARLeafLAbsorpt_pft+DirPARLeafAbsorptionAzclass
        backScatDirPARStemLAbsorpt_pft = backScatDirPARStemLAbsorpt_pft+DirPARStemAbsorptionAzclass
        backScatDirPARDStemLAbsorpt_pft= backScatDirPARDStemLAbsorpt_pft+DirPARDeadStemAbsorpAzclass
        !
        ! FORWARD SCATTERING OF REFLECTED DIRECT RADIATION
        !
      ELSE
        FwdScatDirSWLeafLAbsorpt_pft = FwdScatDirSWLeafLAbsorpt_pft+DirSWLeafAbsorptionAzclass
        FwdScatDirSWStemLAbsorpt_pft = FwdScatDirSWStemLAbsorpt_pft+DirSWStemAbsorptionAzclass
        FwdScatDirSWDStemLAbsorpt_pft= FwdScatDirSWDStemLAbsorpt_pft+DirSWDeadStemAbsorpAzclass

        FwdScatDirPARLeafLAbsorpt_pft = FwdScatDirPARLeafLAbsorpt_pft+DirPARLeafAbsorptionAzclass
        FwdScatDirPARStemLAbsorpt_pft = FwdScatDirPARStemLAbsorpt_pft+DirPARStemAbsorptionAzclass
        FwdScatDirPARDStemLAbsorpt_pft= FwdScatDirPARDStemLAbsorpt_pft+DirPARDeadStemAbsorpAzclass
      ENDIF
      !
      !     INTENSITY OF ABSORBED DIFFUSE RADIATION AT LEAF SURFACES
      !
      !     DiffSWLeafAbsorptAzclass,DiffSWStemAbsorptAzclass,DiffPARLeafAbsorptAzclass,
      !     DiffPARStemAbsorptAzclass=diffuse SW,PAR flux absorbed by leaf,stalk surf
      !     OMEGA,OMEGX=incident angle of diffuse radn at leaf,horizontal surface
      !
      D1750: DO NN=1,NumOfLeafAzimuthSectors1
        RadSWDiffsLMNN              = RadSWDiffusL*OMEGA2Leaf(M,N,NN)
        DiffSWLeafAbsorptAzclass  = RadSWDiffsLMNN*LeafSWabsorptivityL_pft(L,NZ)
        DiffSWStemAbsorptAzclass = RadSWDiffsLMNN*StemSWAbsorptivityL_pft(L,NZ)
        DiffSWDStemAbsorptAzclass = RadSWDiffsLMNN*DStemSWabsorptivityL_pft(L,NZ)

        RadPARDiffsLMNN= RadPARDiffusL*OMEGA2Leaf(M,N,NN)
        DiffPARLeafAbsorptAzclass  = RadPARDiffsLMNN*LeafPARabsorptivityL_pft(L,NZ)
        DiffPARStemAbsorptAzclass = RadPARDiffsLMNN*StemPARAbsorptivityL_pft(L,NZ)
        DiffPARDStemAbsorptAzclass = RadPARDiffsLMNN*DStemPARabsorptivityL_pft(L,NZ)

        RadDifPARAbsorption_zsec(N,M,L,NZ) = RadDifPARAbsorption_zsec(N,M,L,NZ)+DiffPARLeafAbsorptAzclass
        RadTotPARAbsorption_zsec(N,M,L,NZ) = RadTotPARAbsorption_zsec(N,M,L,NZ)+DiffPARLeafAbsorptAzclass
        

        !     ABSORPTION OF DIFFUSE RADIATION BY SHADED LEAF SURFACES
        !
        !     FracDifRadAbsorbt=accumulated horizontal area of intercepted diffuse radiation
        !
        RadDiffSWLeafAbsorbAzc     = EffLeafArea*DiffSWLeafAbsorptAzclass
        RadDiffSWStemAbsorbAzc     = EffStemArea*DiffSWStemAbsorptAzclass
        RadDiffSWDStemAbsorbAzc    = EffDeadStemArea*DiffSWDStemAbsorptAzclass

        RadDifSWLeafAbsorption_pft = RadDifSWLeafAbsorption_pft+RadDiffSWLeafAbsorbAzc
        RadDifSWStemAbsorption_pft = RadDifSWStemAbsorption_pft+RadDiffSWStemAbsorbAzc
        RadDifSWDStemAbsorption_pft= RadDifSWDStemAbsorption_pft+RadDiffSWDStemAbsorbAzc

        RadDiffPARLeafAbsorbAzc     = EffLeafArea*DiffPARLeafAbsorptAzclass
        RadDiffPARStemAbsorbAzc     = EffStemArea*DiffPARStemAbsorptAzclass
        RadDiffPARDStemAbsorbAzc    = EffDeadStemArea*DiffPARDStemAbsorptAzclass

        RadDifPARLeafAbsorption_pft = RadDifPARLeafAbsorption_pft+RadDiffPARLeafAbsorbAzc
        RadDifPARStemAbsorption_pft = RadDifPARStemAbsorption_pft+RadDiffPARStemAbsorbAzc
        RadDifPARDStemAbsorption_pft= RadDifPARDStemAbsorption_pft+RadDiffPARDStemAbsorbAzc

        FractionDifRadAbsorbtL   = FractionDifRadAbsorbtL+(EffLeafAreaAzclass+EffStemAreaAzclass+EffDeadStemAreaAzclass)*OMEGX(M,N,NN)

        !     BACKSCATTERING OF REFLECTED DIFFUSE RADIATION
        !
        IF(iScatteringDiffus(M,N,NN).EQ.ibackward)THEN
          BackScatDifSWLeafLAbsorpt_pft = BackScatDifSWLeafLAbsorpt_pft+RadDiffSWLeafAbsorbAzc
          BackScatDifSWStemLAbsorpt_pft = BackScatDifSWStemLAbsorpt_pft+RadDiffSWStemAbsorbAzc
          BackScatDifSWDStemLAbsorpt_pft= BackScatDifSWDStemLAbsorpt_pft+RadDiffSWDStemAbsorbAzc

          BackScatDifPARLeafLAbsorpt_pft  = BackScatDifPARLeafLAbsorpt_pft+RadDiffPARLeafAbsorbAzc
          BackScatDifPARStemLAbsorpt_pft  = BackScatDifPARStemLAbsorpt_pft+RadDiffPARStemAbsorbAzc
          BackScatDifPARDStemLAbsorpt_pft = BackScatDifPARDStemLAbsorpt_pft+RadDiffPARDStemAbsorbAzc

          !
          !     FORWARD SCATTERING OF REFLECTED DIFFUSE RADIATION
          !
        ELSE
          FwdScatDifSWLeafLAbsorpt_pft = FwdScatDifSWLeafLAbsorpt_pft+RadDiffSWLeafAbsorbAzc
          FwdScatDifSWStemLAbsorpt_pft = FwdScatDifSWStemLAbsorpt_pft+RadDiffSWStemAbsorbAzc
          FwdScatDifSWDStemLAbsorpt_pft=FwdScatDifSWDStemLAbsorpt_pft+RadDiffSWDStemAbsorbAzc

          FwdScatDifPARLeafLAbsorpt_pft = FwdScatDifPARLeafLAbsorpt_pft+RadDiffPARLeafAbsorbAzc
          FwdScatDifPARStemLAbsorpt_pft = FwdScatDifPARStemLAbsorpt_pft+RadDiffPARStemAbsorbAzc
          FwdScatDifPARDStemLAbsorpt_pft=FwdScatDifPARDStemLAbsorpt_pft+RadDiffPARDStemAbsorbAzc
        ENDIF
      ENDDO D1750
    ENDDO D1700
  ENDDO D1600
  call PrintInfo('end '//subname)
  end associate
  end subroutine DownwellRadiationCanopyZL

!----------------------------------------------------------------------------------------------------
  subroutine UpdateOpticalProperties(I,J,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft,SurfAreaZsecDead_lpft)
  ! Description:
  ! Derive snow optics from retained SWE and current exposed area, then distribute
  ! physical plus snow-added area over the exposed leaf/live-stem/dead profile.
  implicit none
  integer, intent(in) :: I,J
  real(r8), intent(inout) :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)    !leaf area in different angle sector
  real(r8), intent(inout) :: StemAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)    !stem area in different angle sector
  real(r8), intent(inout) :: SurfAreaZsecDead_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)  
  character(len=*), parameter :: subname='UpdateOpticalProperties'
  integer :: NZ,L,N
  real(r8), parameter :: dOMEGA_SNOW=0.2_R8 !assuming snow cause 20% reduction of clumping factor
  real(r8) :: StemSWAlbedoEff,StemSWTransmEff
  real(r8) :: StemPARAlbedoEff,StemPARTransmEff
  real(r8) :: LeafSWAlbedoEff,LeafSWTransmEff
  real(r8) :: LeafPARAlbedoEff,leafPARTransmEff
  real(r8) :: DstemSWAlbedoEff,DstemPARAlbedoEff,DstemSWTransmEff,DstemPARTransmEff  
  real(r8) :: PhysicalArea,EffectiveArea,AreaMultiplier

  call PrintInfo('beg '//subname)
  associate(                                                            &
    ZEROS                      => plt_site%ZEROS                       ,& !input  :threshold zero for numerical stability,[-]
    NP                         => plt_site%NP                          ,& !input  :current number of plant species,[-]    
    NU                         => plt_site%NU                          ,& !input  :current soil surface layer number, [-]    
    SnowOnCanopy_pft           => plt_ew%SnowOnCanopy_pft              ,& !input  :snow water equivalent held on canopy, [m3 d-2]    
    ZERO                       => plt_site%ZERO                        ,& !input  :threshold zero for numerical stability, [-]    
    iPlant2ndGrothPattern_pft  => plt_pheno%iPlant2ndGrothPattern_pft  ,& !input  :Plant expression of secondary growth, [0/1] 
    fSnowCanopy_pft            => plt_ew%fSnowCanopy_pft               ,& ! output: retained snow cover on current exposed area [-, PFT]
    SnowType                  => plt_pheno%iPlantSnowIntercepType_pft,& ! input: snow interception class [0:4, PFT]
    LiveArea                  => plt_morph%LeafStalkAreaAct_pft       ,& ! input: exposed physical live area [m2/cell, PFT]
    DeadArea                  => plt_morph%StandDeadSurfAreaAct_pft   ,& ! input: exposed physical dead area [m2/cell, PFT]
    ClumpFactor_pft            => plt_morph%ClumpFactor_pft            ,& !input  :clumping factor for self-shading in canopy layer, [-]    
    SnowDepth                  => plt_ew%SnowDepth                     ,& !input  :snowpack depth, [m]    
    RadSWLeafAlbedo_pft        => plt_rad%RadSWLeafAlbedo_pft          ,& !input  :canopy shortwave albedo, [-]
    RadPARLeafAlbedo_pft       => plt_rad%RadPARLeafAlbedo_pft         ,& !input  :canopy PAR albedo, [-]
    RadPARLeafTransmitance_pft => plt_rad%RadPARLeafTransmitance_pft   ,& !input  :canopy PAR transmissivity, [-]
    RadSWLeafTransmitance_pft  => plt_rad%RadSWLeafTransmitance_pft    ,& !input  :canopy shortwave transmissivity, [-]    
    CanopyLeafArea_pft         => plt_morph%CanopyLeafArea_pft         ,& !input  :plant canopy leaf area, [m2 d-2]    
    AREA3                      => plt_site%AREA3                       ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]    
    CanopyHeightZ_col          => plt_morph%CanopyHeightZ_col          ,& !input  :canopy layer height, [m]    
    ClumpFactorNow_pft         => plt_morph%ClumpFactorNow_pft          & !input  :clumping factor for self-shading in canopy layer at current LAI, [-]    
  )

  DO NZ=1,NP
    PhysicalArea=LiveArea(NZ)+DeadArea(NZ)
    call CanopySnowOptics(SnowType(NZ),PhysicalArea,SnowOnCanopy_pft(NZ),fSnowCanopy_pft(NZ),EffectiveArea)
    AreaMultiplier=1._r8
    if(PhysicalArea>0._r8)AreaMultiplier=EffectiveArea/PhysicalArea
    ClumpFactorNow_pft(NZ)  = ClumpFactor_pft(NZ)*(1.0_r8-0.025_r8*CanopyLeafArea_pft(NZ)/AREA3(NU))
    StemClumpFactor_pft(NZ) = StemClumpFactor(iPlant2ndGrothPattern_pft(NZ))

    if(SnowOnCanopy_pft(NZ).GT.ZERO)then

      !assuming snow cause 20% reduction of clumping factor
      ClumpFactorNow_pft(NZ) = ClumpFactorNow_pft(NZ) *(1._r8-dOMEGA_SNOW*fSnowCanopy_pft(NZ))

      LeafSWAlbedoEff  = (1._r8-fSnowCanopy_pft(NZ))*RadSWLeafAlbedo_pft(NZ)+fSnowCanopy_pft(NZ)*SnowSWAlbedo
      LeafPARAlbedoEff = (1._r8-fSnowCanopy_pft(NZ))*RadPARLeafAlbedo_pft(NZ)+fSnowCanopy_pft(NZ)*SnowPARAlbedo
      LeafSWTransmEff  = (1._r8-fSnowCanopy_pft(NZ))*RadSWLeafTransmitance_pft(NZ)+fSnowCanopy_pft(NZ)*SnowSWTransmitance
      leafPARTransmEff = (1._r8-fSnowCanopy_pft(NZ))*RadPARLeafTransmitance_pft(NZ)+fSnowCanopy_pft(NZ)*SnowPARTransmitance

      StemSWAlbedoEff  = (1._r8-fSnowCanopy_pft(NZ))*RadSWStemAlbedo(iPlant2ndGrothPattern_pft(NZ))+fSnowCanopy_pft(NZ)*SnowSWAlbedo
      StemSWTransmEff  = (1._r8-fSnowCanopy_pft(NZ))*RadSWStemTransmitance(iPlant2ndGrothPattern_pft(NZ))+fSnowCanopy_pft(NZ)*SnowSWTransmitance
      StemPARAlbedoEff = (1._r8-fSnowCanopy_pft(NZ))*RadPARStemAlbedo(iPlant2ndGrothPattern_pft(NZ))+fSnowCanopy_pft(NZ)*SnowPARAlbedo
      StemPARTransmEff = (1._r8-fSnowCanopy_pft(NZ))*RadPARStemTransmitance(iPlant2ndGrothPattern_pft(NZ))+fSnowCanopy_pft(NZ)*SnowPARTransmitance

      DstemSWAlbedoEff = (1._r8-fSnowCanopy_pft(NZ))*StdeadSWAlbedo(iPlant2ndGrothPattern_pft(NZ))+fSnowCanopy_pft(NZ)*SnowSWAlbedo
      DstemPARAlbedoEff= (1._r8-fSnowCanopy_pft(NZ))*StdeadPARAlbedo(iPlant2ndGrothPattern_pft(NZ))+fSnowCanopy_pft(NZ)*SnowPARAlbedo
      DstemSWTransmEff = (1._r8-fSnowCanopy_pft(NZ))*StdeadSWTransmitance(iPlant2ndGrothPattern_pft(NZ))+fSnowCanopy_pft(NZ)*SnowSWTransmitance
      DstemPARTransmEff= (1._r8-fSnowCanopy_pft(NZ))*StdeadPARTransmitance(iPlant2ndGrothPattern_pft(NZ))+fSnowCanopy_pft(NZ)*SnowPARTransmitance

      !move donward
      DO  L=NumCanopyLayers1,1,-1
        IF(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO)THEN      
          DO  N=1,NumLeafInclinationClasses1
            LeafAreaZsecLive_lpft(N,L,NZ) = LeafAreaZsecLive_lpft(N,L,NZ)*AreaMultiplier
            StemAreaZsecLive_lpft(N,L,NZ) = StemAreaZsecLive_lpft(N,L,NZ)*AreaMultiplier
            SurfAreaZsecDead_lpft(N,L,NZ) = SurfAreaZsecDead_lpft(N,L,NZ)*AreaMultiplier
          enddo

          if(L.eq.NumCanopyLayers1)then

            RadSWLeafAlbedoL_pft(L,NZ)        = LeafSWAlbedoEff
            RadPARLeafAlbedoL_pft(L,NZ)       = LeafPARAlbedoEff
            RadSWLeafTransmitanceL_pft(L,NZ)  = LeafSWTransmEff
            RadPARLeafTransmitanceL_pft(L,NZ) = leafPARTransmEff

            RadSWStemAlbedoL_pft(L,NZ)        = StemSWAlbedoEff
            RadPARStemAlbedoL_pft(L,NZ)       = StemPARAlbedoEff
            RadSWStemTransmitanceL_pft(L,NZ)  = StemSWTransmEff
            RadPARStemTransmitanceL_pft(L,NZ) = StemPARTransmEff

            RadSWDStemAlbedoL_pft(L,NZ)        = DstemSWAlbedoEff
            RadPARDstemAlbedoL_pft(L,NZ)       = DstemPARAlbedoEff
            RadSWDstemTransmitanceL_pft(L,NZ)  = DstemSWTransmEff
            RadPARDstemTransmitanceL_pft(L,NZ) = DstemPARTransmEff

            LeafSWabsorptivityL_pft(L,NZ)  = 1._r8-RadSWLeafAlbedoL_pft(L,NZ)-RadSWLeafTransmitanceL_pft(L,NZ)
            LeafPARabsorptivityL_pft(L,NZ) = 1._r8-RadPARLeafAlbedoL_pft(L,NZ)-RadPARLeafTransmitanceL_pft(L,NZ)

            StemSWabsorptivityL_pft(L,NZ)  = 1._r8-RadSWStemAlbedoL_pft(L,NZ)-RadSWStemTransmitanceL_pft(L,NZ)
            StemPARabsorptivityL_pft(L,NZ) = 1._r8-RadPARStemAlbedoL_pft(L,NZ)-RadPARStemTransmitanceL_pft(L,NZ)

            DStemSWabsorptivityL_pft(L,NZ) = 1._r8-RadSWDStemAlbedoL_pft(L,NZ)-RadSWDstemTransmitanceL_pft(L,NZ)
            DStemPARabsorptivityL_pft(L,NZ)= 1._r8-RadPARDstemAlbedoL_pft(L,NZ)-RadPARDstemTransmitanceL_pft(L,NZ)

            RadSWLeafAlbedoL_pft(L,NZ)        = RadSWLeafAlbedoL_pft(L,NZ)/LeafSWabsorptivityL_pft(L,NZ)
            RadPARLeafAlbedoL_pft(L,NZ)       = RadPARLeafAlbedoL_pft(L,NZ)/LeafPARabsorptivityL_pft(L,NZ)
            RadSWLeafTransmitanceL_pft(L,NZ)  = RadSWLeafTransmitanceL_pft(L,NZ)/LeafSWabsorptivityL_pft(L,NZ)
            RadPARLeafTransmitanceL_pft(L,NZ) = RadPARLeafTransmitanceL_pft(L,NZ)/LeafPARabsorptivityL_pft(L,NZ)

            RadSWStemAlbedoL_pft(L,NZ)        = RadSWStemAlbedoL_pft(L,NZ)/StemSWabsorptivityL_pft(L,NZ)
            RadPARStemAlbedoL_pft(L,NZ)       = RadPARStemAlbedoL_pft(L,NZ)/StemPARabsorptivityL_pft(L,NZ)
            RadSWStemTransmitanceL_pft(L,NZ)  = RadSWStemTransmitanceL_pft(L,NZ)/StemSWabsorptivityL_pft(L,NZ)
            RadPARStemTransmitanceL_pft(L,NZ) = RadPARStemTransmitanceL_pft(L,NZ)/StemPARabsorptivityL_pft(L,NZ)

            RadSWDStemAlbedoL_pft(L,NZ)        = RadSWDStemAlbedoL_pft(L,NZ) /DStemSWabsorptivityL_pft(L,NZ)
            RadPARDstemAlbedoL_pft(L,NZ)       = RadPARDstemAlbedoL_pft(L,NZ) /DStemPARabsorptivityL_pft(L,NZ)
            RadSWDstemTransmitanceL_pft(L,NZ)  = RadSWDstemTransmitanceL_pft(L,NZ)  /DStemSWabsorptivityL_pft(L,NZ)
            RadPARDstemTransmitanceL_pft(L,NZ) = RadPARDstemTransmitanceL_pft(L,NZ) /DStemPARabsorptivityL_pft(L,NZ)

          else
            LeafSWabsorptivityL_pft(L,NZ)  = LeafSWabsorptivityL_pft(L+1,NZ)
            LeafPARabsorptivityL_pft(L,NZ) = LeafPARabsorptivityL_pft(L+1,NZ)

            StemSWabsorptivityL_pft(L,NZ)  = StemSWabsorptivityL_pft(L+1,NZ)
            StemPARabsorptivityL_pft(L,NZ) = StemPARabsorptivityL_pft(L+1,NZ)

            DStemSWabsorptivityL_pft(L,NZ) = DStemSWabsorptivityL_pft(L+1,NZ) 
            DStemPARabsorptivityL_pft(L,NZ)= DStemPARabsorptivityL_pft(L+1,NZ)

            RadSWLeafAlbedoL_pft(L,NZ)        = RadSWLeafAlbedoL_pft(L+1,NZ)
            RadPARLeafAlbedoL_pft(L,NZ)       = RadPARLeafAlbedoL_pft(L+1,NZ)
            RadSWLeafTransmitanceL_pft(L,NZ)  = RadSWLeafTransmitanceL_pft(L+1,NZ)
            RadPARLeafTransmitanceL_pft(L,NZ) = RadPARLeafTransmitanceL_pft(L+1,NZ)

            RadSWStemAlbedoL_pft(L,NZ)        = RadSWStemAlbedoL_pft(L+1,NZ)
            RadPARStemAlbedoL_pft(L,NZ)       = RadPARStemAlbedoL_pft(L+1,NZ)
            RadSWStemTransmitanceL_pft(L,NZ)  = RadSWStemTransmitanceL_pft(L+1,NZ)
            RadPARStemTransmitanceL_pft(L,NZ) = RadPARStemTransmitanceL_pft(L+1,NZ)

            RadSWDStemAlbedoL_pft(L,NZ)        = RadSWDStemAlbedoL_pft(L+1,NZ)
            RadPARDstemAlbedoL_pft(L,NZ)       = RadPARDstemAlbedoL_pft(L+1,NZ) 
            RadSWDstemTransmitanceL_pft(L,NZ)  = RadSWDstemTransmitanceL_pft(L+1,NZ) 
            RadPARDstemTransmitanceL_pft(L,NZ) = RadPARDstemTransmitanceL_pft(L+1,NZ)
          endif
        endif
      enddo            
    else
      
      DO  L=NumCanopyLayers1,1,-1
        if(L.eq.NumCanopyLayers1)THEN
          RadSWLeafAlbedoL_pft(L,NZ)        = RadSWLeafAlbedo_pft(NZ)
          RadPARLeafAlbedoL_pft(L,NZ)       = RadPARLeafAlbedo_pft(NZ)
          RadSWLeafTransmitanceL_pft(L,NZ)  = RadSWLeafTransmitance_pft(NZ)
          RadPARLeafTransmitanceL_pft(L,NZ) = RadPARLeafTransmitance_pft(NZ)

          RadSWStemAlbedoL_pft(L,NZ)        = RadSWStemAlbedo(iPlant2ndGrothPattern_pft(NZ))
          RadPARStemAlbedoL_pft(L,NZ)       = RadPARStemAlbedo(iPlant2ndGrothPattern_pft(NZ))
          RadSWStemTransmitanceL_pft(L,NZ)  = RadSWStemTransmitance(iPlant2ndGrothPattern_pft(NZ))
          RadPARStemTransmitanceL_pft(L,NZ) = RadPARStemTransmitance(iPlant2ndGrothPattern_pft(NZ))

          RadSWDStemAlbedoL_pft(L,NZ)        = StdeadSWAlbedo(iPlant2ndGrothPattern_pft(NZ))
          RadPARDstemAlbedoL_pft(L,NZ)       = StdeadPARAlbedo(iPlant2ndGrothPattern_pft(NZ))
          RadSWDstemTransmitanceL_pft(L,NZ)  = StdeadSWTransmitance(iPlant2ndGrothPattern_pft(NZ))
          RadPARDstemTransmitanceL_pft(L,NZ) = StdeadPARTransmitance(iPlant2ndGrothPattern_pft(NZ))

          LeafSWabsorptivityL_pft(L,NZ)  = 1._r8-RadSWLeafAlbedoL_pft(L,NZ)-RadSWLeafTransmitanceL_pft(L,NZ)
          LeafPARabsorptivityL_pft(L,NZ) = 1._r8-RadPARLeafAlbedoL_pft(L,NZ)-RadPARLeafTransmitanceL_pft(L,NZ)

          StemSWabsorptivityL_pft(L,NZ)  = 1._r8-RadSWStemAlbedoL_pft(L,NZ)-RadSWStemTransmitanceL_pft(L,NZ)
          StemPARabsorptivityL_pft(L,NZ) = 1._r8-RadPARStemAlbedoL_pft(L,NZ)-RadPARStemTransmitanceL_pft(L,NZ)

          DStemSWabsorptivityL_pft(L,NZ) = 1._r8-RadSWDStemAlbedoL_pft(L,NZ)-RadSWDstemTransmitanceL_pft(L,NZ)
          DStemPARabsorptivityL_pft(L,NZ)= 1._r8-RadPARDstemAlbedoL_pft(L,NZ)-RadPARDstemTransmitanceL_pft(L,NZ)

          RadSWLeafAlbedoL_pft(L,NZ)        = RadSWLeafAlbedoL_pft(L,NZ)/LeafSWabsorptivityL_pft(L,NZ)
          RadPARLeafAlbedoL_pft(L,NZ)       = RadPARLeafAlbedoL_pft(L,NZ) /LeafPARabsorptivityL_pft(L,NZ)
          RadSWLeafTransmitanceL_pft(L,NZ)  = RadSWLeafTransmitanceL_pft(L,NZ) /LeafSWabsorptivityL_pft(L,NZ)
          RadPARLeafTransmitanceL_pft(L,NZ) = RadPARLeafTransmitanceL_pft(L,NZ) /LeafPARabsorptivityL_pft(L,NZ)

          RadSWStemAlbedoL_pft(L,NZ)        = RadSWStemAlbedoL_pft(L,NZ) /StemSWabsorptivityL_pft(L,NZ)
          RadPARStemAlbedoL_pft(L,NZ)       = RadPARStemAlbedoL_pft(L,NZ)/StemPARabsorptivityL_pft(L,NZ)
          RadSWStemTransmitanceL_pft(L,NZ)  = RadSWStemTransmitanceL_pft(L,NZ)/StemSWabsorptivityL_pft(L,NZ)
          RadPARStemTransmitanceL_pft(L,NZ) = RadPARStemTransmitanceL_pft(L,NZ)/StemPARabsorptivityL_pft(L,NZ)

          RadSWDStemAlbedoL_pft(L,NZ)        = RadSWDStemAlbedoL_pft(L,NZ) /DStemSWabsorptivityL_pft(L,NZ)
          RadPARDstemAlbedoL_pft(L,NZ)       = RadPARDstemAlbedoL_pft(L,NZ) /DStemPARabsorptivityL_pft(L,NZ)
          RadSWDstemTransmitanceL_pft(L,NZ)  = RadSWDstemTransmitanceL_pft(L,NZ)  /DStemSWabsorptivityL_pft(L,NZ)
          RadPARDstemTransmitanceL_pft(L,NZ) = RadPARDstemTransmitanceL_pft(L,NZ) /DStemPARabsorptivityL_pft(L,NZ)          

        ELSE
          LeafSWabsorptivityL_pft(L,NZ)  = LeafSWabsorptivityL_pft(L+1,NZ)
          LeafPARabsorptivityL_pft(L,NZ) = LeafPARabsorptivityL_pft(L+1,NZ)
          
          StemSWabsorptivityL_pft(L,NZ)  = StemSWabsorptivityL_pft(L+1,NZ)
          StemPARabsorptivityL_pft(L,NZ) = StemPARabsorptivityL_pft(L+1,NZ)

          DStemSWabsorptivityL_pft(L,NZ) = DStemSWabsorptivityL_pft(L+1,NZ) 
          DStemPARabsorptivityL_pft(L,NZ)= DStemPARabsorptivityL_pft(L+1,NZ)

          RadSWLeafAlbedoL_pft(L,NZ)        = RadSWLeafAlbedoL_pft(L+1,NZ)
          RadPARLeafAlbedoL_pft(L,NZ)       = RadPARLeafAlbedoL_pft(L+1,NZ)
          RadSWLeafTransmitanceL_pft(L,NZ)  = RadSWLeafTransmitanceL_pft(L+1,NZ)
          RadPARLeafTransmitanceL_pft(L,NZ) = RadPARLeafTransmitanceL_pft(L+1,NZ)

          RadSWStemAlbedoL_pft(L,NZ)        = RadSWStemAlbedoL_pft(L+1,NZ)
          RadPARStemAlbedoL_pft(L,NZ)       = RadPARStemAlbedoL_pft(L+1,NZ)
          RadSWStemTransmitanceL_pft(L,NZ)  = RadSWStemTransmitanceL_pft(L+1,NZ)
          RadPARStemTransmitanceL_pft(L,NZ) = RadPARStemTransmitanceL_pft(L+1,NZ)

          RadSWDStemAlbedoL_pft(L,NZ)        = RadSWDStemAlbedoL_pft(L+1,NZ)
          RadPARDstemAlbedoL_pft(L,NZ)       = RadPARDstemAlbedoL_pft(L+1,NZ) 
          RadSWDstemTransmitanceL_pft(L,NZ)  = RadSWDstemTransmitanceL_pft(L+1,NZ) 
          RadPARDstemTransmitanceL_pft(L,NZ) = RadPARDstemTransmitanceL_pft(L+1,NZ)          
        ENDIF
      ENDDO
    endif
  ENDDO
  end associate
  call PrintInfo('end '//subname)
  end subroutine UpdateOpticalProperties
!----------------------------------------------------------------------------------------------------

  subroutine AngularDistributeRadiation(SolarAzimuthAngle,CosineSunInclAngle,BETX,iScatteringDirect,&
    RadDirSW2Leaf_zsec,RadDirSW2Stem_zsec,RadDirPAR2Leaf_zsec,RadDirPAR2Stem_zsec)
  implicit none
  real(r8), intent(in) :: SolarAzimuthAngle  
  real(r8), intent(in) :: CosineSunInclAngle  
  integer , INTENT(OUT) :: iScatteringDirect(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1)  
  real(r8), intent(out) :: BETX(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1)                          !sine of direct solar radiation on leaf surface/sine
  REAL(R8), intent(out) :: RadDirSW2Leaf_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)     !direct shortwave radiation on leave in each canopy sector
  real(r8), intent(out) :: RadDirSW2Stem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)    !direct shortwave radiation on stalk in each canopy sector
  real(r8), intent(out) :: RadDirPAR2Leaf_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)    !PAR radiation on leave in each canopy sector
  real(r8), intent(out) :: RadDirPAR2Stem_zsec(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1,JP1)   !PAR radiation on stalk in each canopy sector
  real(r8)  :: BETA(NumLeafInclinationClasses1,NumOfSkyAzimuthSects1)                          !sine of direct solar radiation on leaf surface, [-]
  REAL(R8) :: SolarAngle  
  REAL(R8) :: LeafAzimuthAngle,ZAGL,DAZI
  real(r8) :: BETY,BETZ,DAzimuthAngle  
  integer :: M,N,NZ,L  
  character(len=*), parameter :: subname='subroutine AngularDistributeRadiation'

  associate(                                                              &
    CanopyHeightZ_col          => plt_morph%CanopyHeightZ_col         ,&
    SnowDepth                  => plt_ew%SnowDepth                     ,&
    ZERO                       => plt_site%ZERO                       ,&
    NP                           => plt_site%NP                          ,& !input  :current number of plant species,[-]  
    LeafPARabsorptivity_pft      => plt_rad%LeafPARabsorptivity_pft      ,& !input  :canopy PAR absorptivity,[-]
    LeafSWabsorptivity_pft       => plt_rad%LeafSWabsorptivity_pft       ,& !input  :canopy shortwave absorptivity, [-]  
    RadSWDirect_col              => plt_rad%RadSWDirect_col              ,& !input  :direct shortwave radiation, [W m-2]
    RadDirectPAR_col             => plt_rad%RadDirectPAR_col             ,& !input  :direct PAR, [umol m-2 s-1]      
    iPlant2ndGrothPattern_pft    => plt_pheno%iPlant2ndGrothPattern_pft  ,& !input  :plant expression of secondary growth, [-]    
    SineSunInclinationAngle_col  => plt_rad%SineSunInclinationAngle_col  ,& !input  :sine of solar angle, [-]  
    CosineLeafAngle              => plt_rad%CosineLeafAngle              ,& !input  :cosine of leaf angle,[-]    
    SineLeafAngle                => plt_rad%SineLeafAngle                ,& !input  :sine of leaf angle,[-]  
    RadTotPARAbsorption_zsec     => plt_rad%RadTotPARAbsorption_zsec     ,& !inoput :total incoming PAR, [umol m-2 s-1]
    RadDifPARAbsorption_zsec     => plt_rad%RadDifPARAbsorption_zsec      & !inoput :diffuse incoming PAR absorbed, [umol m-2 s-1]
  )
  call PrintInfo('beg '//subname)
  SolarAngle=ASIN(SineSunInclinationAngle_col)
  !
  !     ABSORBED RADIATION FROM OPTICAL PROPERTIES ENTERED IN 'READS'
  !
  ! First separate into different sectors, then attenuate through canopy layers
  !
  !     ANGLES BETWEEN SUN OR SKY ZONES AND FOLIAR SURFACES
  !
  !     LeafAzimuthAngle=leaf azimuth
  !     BETA,BETX=incident angle of direct radiation at leaf,horizontal surface
  !     ZAGL=determines forward vs backscattering
  !     iScatteringDirect=flag for forward vs backscattering
  !
  DAzimuthAngle=PICON/real(NumOfSkyAzimuthSects1,r8)
  D1100: DO M=1,NumOfSkyAzimuthSects1
    LeafAzimuthAngle = SolarAzimuthAngle+(M-0.5_r8)*DAzimuthAngle
    DAZI             = COS(LeafAzimuthAngle-SolarAzimuthAngle)
    DO N = 1, NumLeafInclinationClasses1
      !BETY=cos(\theta_{L,M}), incident angle on leaf (leafAzimuthAngle,leafInclinationAngle), in eq. (12) of Grant et al. 1989
      BETY      = CosineLeafAngle(N)*SineSunInclinationAngle_col+SineLeafAngle(N)*CosineSunInclAngle*DAZI  
      BETA(N,M) = ABS(BETY)
      BETX(N,M) = BETA(N,M)/SineSunInclinationAngle_col

      !compute incident angle BETZ, eq. (12)/(13)
      IF(CosineLeafAngle(N).GT.SineSunInclinationAngle_col)THEN
        BETZ=ACOS(BETY)
      ELSE
        BETZ=-ACOS(BETY)
      ENDIF

      IF(BETZ.GT.-PICON2h)THEN
        ZAGL=SolarAngle+2.0_r8*BETZ
      ELSE
        ZAGL=SolarAngle-2.0_r8*(PICON+BETZ)
      ENDIF
      !
      IF(ZAGL.GT.0.0_r8 .AND. ZAGL.LT.PICON)THEN
        iScatteringDirect(N,M)=ibackward
      ELSE
        iScatteringDirect(N,M)=iforward
      ENDIF
      !
      !     INTENSITY OF ABSORBED DIRECT RADIATION AT LEAF SURFACES, direct radiation absorption in each leaf class
      !
      DO  NZ=1,NP
        !Whole canopy incident radiation in sector (N,M), specified by azimuth (M) and zenith angles (N), eq. (16).
        RadDirSW2Leaf_zsec(N,M,NZ)  = RadSWDirect_col*BETA(N,M)
        RadDirSW2Stem_zsec(N,M,NZ)  = RadSWDirect_col*BETA(N,M)
        RadDirPAR2Leaf_zsec(N,M,NZ) = RadDirectPAR_col*BETA(N,M)
        RadDirPAR2Stem_zsec(N,M,NZ) = RadDirectPAR_col*BETA(N,M)

        DO L=1,NumCanopyLayers1
          RadDifPARAbsorption_zsec(N,M,L,NZ) = 0.0_r8
          ! Optical properties are defined only above the snow surface.
          ! Buried layers must not reuse an earlier column/hour's coefficients.
          if(CanopyHeightZ_col(L-1).GE.SnowDepth-ZERO)then
            RadTotPARAbsorption_zsec(N,M,L,NZ) = RadDirPAR2Leaf_zsec(N,M,NZ)*LeafPARabsorptivityL_pft(L,NZ)
          else
            RadTotPARAbsorption_zsec(N,M,L,NZ) = 0._r8
          endif
        enddo
      enddo
    enddo
  ENDDO D1100
  call PrintInfo('end '//subname)
  end associate
  end subroutine AngularDistributeRadiation

  ![tail]
end module SurfaceRadiationMod
