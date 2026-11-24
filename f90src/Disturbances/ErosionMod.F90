module ErosionMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : isclose,AZMAX1
  use DebugToolMod, only : PrintInfo
  use EcoSIMConfig, only : nlbiomcp => NumLiveMicrbCompts
  use EcoSIMConfig, only : ndbiomcp=> NumDeadMicrbCompts
  use EcoSIMConfig, only : jcplx1=> jcplxcm1, NumMicbAFunGrupsPerCmplx, NumMicbHFunGrupsPerCmplx,jcplx=>jcplxc
  use EcoSIMConfig, only : column_mode
  use EcoSiMParDataMod    , only : micpar  
  use EcoSIMCtrlMod, only : iErosionMode
  use MicrobialDataType
  use SOMDataType
  use EcoSIMSolverPar
  use FertilizerDataType
  use GridConsts
  use FlagDataType
  use SoilPhysDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use LandSurfDataType
  use SurfSoilDataType
  use ChemTranspDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  USE SedimentDataType
  use GridDataType
  implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  real(r8), PARAMETER :: FSINK=0.01_r8

  real(r8), allocatable :: SedErosionM(:,:,:,:)
  real(r8), allocatable :: TERSED_col(:,:)
  real(r8), allocatable :: RDTSED_col(:,:)
  real(r8), allocatable :: FVOLIM(:,:)
  real(r8), allocatable :: FVLWatMicPM(:,:)
  real(r8), allocatable :: FracVol4Erosion(:,:)
  real(r8), allocatable :: BaseErosionRate(:,:)

  public :: erosion
  public :: InitErosion
  public :: DestructErosion
  contains

  subroutine InitErosion()

  implicit none

  character(len=*), parameter :: subname='InitErosion'

  call PrintInfo('beg '//subname)
  allocate(SedErosionM(2,2,JV,JH))
  allocate(TERSED_col(JY,JX))
  allocate(RDTSED_col(JY,JX))
  allocate(FVOLIM(JY,JX))
  allocate(FVLWatMicPM(JY,JX))
  allocate(FracVol4Erosion(JY,JX))
  allocate(BaseErosionRate(JY,JX))

  call PrintInfo('end '//subname)
  end subroutine InitErosion
!------------------------------------------------------------------------------------------

  SUBROUTINE erosion(I,J,NHW,NHE,NVN,NVS)
      !
!     THIS SUBROUTINE CALCULATES DETACHMENT AND OVERLAND TRANSPORT
!     OF SURFACE SEDIMENT FROM PRECIPITATION IN WEATHER FILE AND
!     FROM RUNOFF IN 'WATSUB'
  implicit none

  character(len=*), parameter :: subname='erosion'

  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: M

!     execution begins here
!
  call PrintInfo('beg '//subname)
  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    DO M=1,NPH
      call SedimentDetachmentM(M,NHW,NHE,NVN,NVS)

      call SedimentTransportM(M,NHW,NHE,NVN,NVS)
    ENDDO
  !
  !     INTERNAL SEDIMENT FLUXES
  !
    call LateralXGridSedmentFlux(NHW, NHE,NVN,NVS)
  !
  !     EXTERNAL BOUNDARY SEDIMENT FLUXES
  !
    if(.not.column_mode)call XBoundErosionFlux(NHW,NHE,NVN,NVS)
  ENDIF
  call PrintInfo('end '//subname)
  END subroutine erosion
!------------------------------------------------------------------------------------------

  subroutine SedimentDetachmentM(M,NHW,NHE,NVN,NVS)
!     INTERNAL TIME STEP AT WHICH SEDIMENT DETACHMENT AND TRANSPORT
!     IS CALCULATED. DETACHMENT IS THE SUM OF THAT BY RAINFALL AND
!     OVERLAND FLOW
!
  implicit none

  character(len=*), parameter :: subname='SedimentDetachmentM'

  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8) :: DETW,SEDX,CSEDD,CSEDX
  real(r8) :: DEPI,DETI,SoilDetachRate,STPR
  integer :: NY,NX


  call PrintInfo('beg '//subname)
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      TERSED_col(NY,NX)=0._r8
      RDTSED_col(NY,NX)=0._r8
      FVOLIM(NY,NX)=AMIN1(1.0_r8,AZMAX1(XVLiceMicPM(M,NY,NX)/VWatStoreCapSurf_col(NY,NX)))
      FVLWatMicPM(NY,NX)=AMIN1(1.0_r8,AZMAX1(XVLMobileWatMicPM(M,NY,NX)/VWatStoreCapSurf_col(NY,NX)))
      FracVol4Erosion(NY,NX)=(1.0_r8-FVOLIM(NY,NX))*FVLWatMicPM(NY,NX)
!
!     DETACHMENT BY RAINFALL WHEN SURFACE WATER IS PRESENT
!
      IF(SoilBulkDensity_vr(NU_col(NY,NX),NY,NX).GT.ZERO.AND.EnergyImpact4Erosion_colM(M,NY,NX).GT.0.0_r8 &
        .AND.XVLMobileWatMicPM(M,NY,NX).GT.ZEROS(NY,NX))THEN
!
!     DETACHMENT OF SEDIMENT FROM SURFACE SOIL DEPENDS ON RAINFALL
!     KINETIC ENERGY AND FROM DETACHMENT COEFFICIENT IN 'HOUR1'
!     ATTENUATED BY DEPTH OF SURFACE WATER
!
        DETW=SoilDetachability4Erosion1(NY,NX)*(1.0_r8+2.0_r8*VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)/VLMicP_vr(NU_col(NY,NX),NY,NX))
        SoilDetachRate=AMIN1(VLSoilMicPMass_vr(NU_col(NY,NX),NY,NX)*dts_wat &
          ,DETW*EnergyImpact4Erosion_colM(M,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX) &
          *FracSoiAsMicP_vr(NU_col(NY,NX),NY,NX)*FracSurfSnoFree_col(NY,NX)*(1.0-FVOLIM(NY,NX)))
        RDTSED_col(NY,NX)=RDTSED_col(NY,NX)+SoilDetachRate
      ENDIF
!
!     DEPOSITION OF SEDIMENT TO SOIL SURFACE FROM IMMOBILE SURFACE WATER
!

      IF(SoilBulkDensity_vr(NU_col(NY,NX),NY,NX).GT.ZERO.AND.FracVol4Erosion(NY,NX).GT.ZERO)THEN

        SEDX=SED_col(NY,NX)+RDTSED_col(NY,NX)
        IF(XVLMobileWaterLitRM(M,NY,NX).LE.VWatStoreCapSurf_col(NY,NX))THEN
          IF(SEDX.GT.ZEROS(NY,NX))THEN
            CSEDD=AZMAX1(SEDX/XVLMobileWatMicPM(M,NY,NX))
        
            DEPI=AMAX1(-SEDX,VLS_col(NY,NX)*(0.0_r8-CSEDD)*AREA_3D(3,NU_col(NY,NX),NY,NX) &
              *FracVol4Erosion(NY,NX)*FracSoiAsMicP_vr(NU_col(NY,NX),NY,NX)*dts_HeatWatTP)
            RDTSED_col(NY,NX)=RDTSED_col(NY,NX)+DEPI
          ENDIF
        ELSE
          !     DETACHMENT IN SURFACE WATER FROM OVERLAND WATER
          !     VELOCITY FROM 'WATSUB' USED TO CALCULATE STREAM POWER,
          !     AND FROM SEDIMENT TRANSPORT CAPACITY VS. CURRENT SEDIMENT
          !     CONCENTRATION IN SURFACE WATER, MODIFIED BY SOIL COHESION
          !     FROM 'HOUR1'
          !
          STPR=1.0E+02_r8*RunoffVelocityM_col(M,NY,NX)*ABS(SLOPE_col(0,NY,NX))
          CSEDX=PtclDensitySurfLay_col(NY,NX)*CER_col(NY,NX)*AZMAX1(STPR-0.4_r8)**XER_col(NY,NX)
          CSEDD=AZMAX1(SEDX/XVLMobileWatMicPM(M,NY,NX))

          IF(CSEDX.GT.CSEDD)THEN
            DETI=AMIN1(VLSoilMicPMass_vr(NU_col(NY,NX),NY,NX)*dts_wat &
              ,SoilDetachability4Erosion2(NY,NX)*(CSEDX-CSEDD)*AREA_3D(3,NU_col(NY,NX),NY,NX) &
              *FracVol4Erosion(NY,NX)*FracSoiAsMicP_vr(NU_col(NY,NX),NY,NX)*dts_HeatWatTP)
          ELSE
            IF(SEDX.GT.ZEROS(NY,NX))THEN
              DETI=AMAX1(-SEDX,VLS_col(NY,NX)*(CSEDX-CSEDD)*AREA_3D(3,NU_col(NY,NX),NY,NX) &
                *FracVol4Erosion(NY,NX)*FracSoiAsMicP_vr(NU_col(NY,NX),NY,NX)*dts_HeatWatTP)
            ELSE
              DETI=0._r8
            ENDIF
          ENDIF
          RDTSED_col(NY,NX)=RDTSED_col(NY,NX)+DETI
        ENDIF
      ENDIF
!
!
!     TRANSPORT OF SEDIMENT IN OVERLAND FLOW FROM SEDIMENT
!     CONCENTRATION TIMES OVERLAND WATER FLUX FROM 'WATSUB'
!

    ENDDO
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine SedimentDetachmentM
!------------------------------------------------------------------------------------------
  subroutine OverLandFlowSedTransp(M,NY,NX,NHW,NHE,NVN,NVS)

  implicit none

  character(len=*), parameter :: subname='OverLandFlowSedTransp'
  integer, intent(in) :: M,NY,NX,NHW,NHE,NVN,NVS
  integer :: N1,N2,N,NN,N4,N5,N4B,N5B
  real(r8) :: SEDX,CSEDE,FERM

  call PrintInfo('beg '//subname)
  N1=NX
  N2=NY
!     SEDX=SED_col(N2,N1)
  IF(SurfRunoffPotentM_col(M,N2,N1).LE.0.0_r8.OR.SoilBulkDensity_vr(NU_col(N2,N1),N2,N1).LE.ZERO)THEN
    BaseErosionRate(N2,N1)=0._r8
  ELSE
    IF(XVLMobileWatMicPM(M,N2,N1).GT.ZEROS2(N2,N1))THEN
      SEDX=SED_col(N2,N1)+RDTSED_col(N2,N1)
      CSEDE=AZMAX1(SEDX/XVLMobileWatMicPM(M,N2,N1))
       BaseErosionRate(N2,N1)=AMIN1(SEDX,CSEDE*SurfRunoffPotentM_col(M,N2,N1)*(1.0_r8-FVOLIM(N2,N1)))
    ELSE
      BaseErosionRate(N2,N1)=0._r8
    ENDIF
  ENDIF
!
!     LOCATE INTERNAL BOUNDARIES
!
  IF(BaseErosionRate(N2,N1).GT.0.0_r8)THEN
    DO  N=1,2
      DO  NN=1,2
        IF(N.EQ.iWestEastDirection)THEN
          !skip boundariers        
          IF(NX.EQ.NHE .AND. NN.EQ.iFront .OR. NX.EQ.NHW .AND. NN.EQ.iBehind)THEN
            cycle
          !internal grids  
          ELSE            
            N4=NX+1;N5=NY    !dest in the east        
            N4B=NX-1;N5B=NY  !dest in the west
          ENDIF
        ELSEIF(N.EQ.iNorthSouthDirection)THEN
          !skip boundariers
          IF(NY.EQ.NVS .AND. NN.EQ.iFront .OR. NY.EQ.NVN .AND. NN.EQ.iBehind)THEN
            cycle
          ELSE            
            N4=NX;N5=NY+1   !dest in the south             
            N4B=NX;N5B=NY-1 !dest in the north
          ENDIF
        ENDIF

        IF(BaseErosionRate(N2,N1).GT.ZEROS(N2,N1))THEN
          IF(NN.EQ.iFront)THEN  !            
            FERM                      = QflxSurfRunoffM_2DH(M,N,iBehind,N5,N4)/SurfRunoffPotentM_col(M,N2,N1)
            SedErosionM(N,iBehind,N5,N4)    = BaseErosionRate(N2,N1)*FERM
          ENDIF

          IF(NN.EQ.iBehind)THEN
            IF(N4B.GT.0.AND.N5B.GT.0)THEN
              !well-defined dest grid
              FERM                             = QflxSurfRunoffM_2DH(M,N,iFront,N5B,N4B)/SurfRunoffPotentM_col(M,N2,N1)
              SedErosionM(N,iFront,N5B,N4B)    = BaseErosionRate(N2,N1)*FERM
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine OverLandFlowSedTransp
!------------------------------------------------------------------------------------------
  subroutine SedimentTransportM(M,NHW,NHE,NVN,NVS)
!     INTERNAL TIME STEP AT WHICH SEDIMENT DETACHMENT AND TRANSPORT
!     IS CALCULATED. DETACHMENT IS THE SUM OF THAT BY RAINFALL AND
!     OVERLAND FLOW
!
  implicit none

  integer, intent(in) :: M,NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='SedimentTransportM'
  real(r8) :: CSEDE,SEDX,AvailableSediment,ProposedOutflow,SedimentScale,AcceptedOutflow
  real(r8) :: OutwardSign(2,2)
  integer :: N1,N2,NY,NX,N,Side,FY,FX,NN
  integer :: FaceNN(2,2),FaceNY(2,2),FaceNX(2,2)
  logical :: InternalFace(2,2)

  call PrintInfo('beg '//subname)
  ! Reset instantaneous faces once, including quiet cells and halo slots.
  ! Keep all donor states unchanged until every face has been finalized.
  SedErosionM=0._r8
  TERSED_col(NVN:NVS,NHW:NHE)=0._r8
!
!     BOUNDARY SEDIMENT FLUXES
!
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      N1=NX;N2=NY
      call OverLandFlowSedTransp(M,NY,NX,NHW,NHE,NVN,NVS)
      
      IF(SurfRunoffPotentM_col(M,N2,N1).LE.0.0.OR.SoilBulkDensity_vr(NU_col(N2,N1),N2,N1).LE.ZERO)THEN
        BaseErosionRate(N2,N1)=0._r8
      ELSE
        IF(XVLMobileWatMicPM(M,N2,N1).GT.ZEROS2(N2,N1))THEN
          SEDX  = SED_col(NY,NX)+RDTSED_col(NY,NX)
          CSEDE = AZMAX1(SEDX/XVLMobileWatMicPM(M,N2,N1))
          BaseErosionRate(N2,N1)=AMIN1(SEDX,CSEDE*SurfRunoffPotentM_col(M,N2,N1))
        ELSE
          BaseErosionRate(N2,N1)=0._r8
        ENDIF
      ENDIF
!
      if(.not. column_mode)call XBoundSedTranspM(M,NY,NX,NHW,NHE,NVN,NVS,N1,N2)
    ENDDO
  ENDDO

  ! Accept all outward faces against one donor snapshot. Incoming sediment
  ! does not fund outgoing transport within this explicit substep.
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO N=1,2
        FaceNY(N,:)=NY;FaceNX(N,:)=NX
        if(N==iWestEastDirection)then
          FaceNX(N,iFront)=NX+1;FaceNX(N,iBehind)=NX-1
          InternalFace(N,iFront)=NX<NHE
          InternalFace(N,iBehind)=NX>NHW
        else
          FaceNY(N,iFront)=NY+1;FaceNY(N,iBehind)=NY-1
          InternalFace(N,iFront)=NY<NVS
          InternalFace(N,iBehind)=NY>NVN
        endif
        ! Internal outflows are positive in receiving slots. External E/S
        ! outflows are positive outside; external W/N outflows are negative
        ! in this donor's slot, matching the water and material flux ledgers.
        FaceNN(N,iFront)=iBehind;FaceNN(N,iBehind)=iFront
        OutwardSign(N,:)=1._r8
        if(.not.InternalFace(N,iFront))FaceNN(N,iFront)=iFront
        if(.not.InternalFace(N,iBehind))then
          FaceNN(N,iBehind)=iBehind
          FaceNY(N,iBehind)=NY;FaceNX(N,iBehind)=NX
          OutwardSign(N,iBehind)=-1._r8
        endif
      ENDDO

      ProposedOutflow=0._r8
      DO N=1,2
        DO Side=iFront,iBehind
          ProposedOutflow=ProposedOutflow+AZMAX1(OutwardSign(N,Side) &
            *SedErosionM(N,FaceNN(N,Side),FaceNY(N,Side),FaceNX(N,Side)))
        ENDDO
      ENDDO
      AvailableSediment=AZMAX1(SED_col(NY,NX)+RDTSED_col(NY,NX))
      SedimentScale=1._r8
      if(ProposedOutflow>AvailableSediment)SedimentScale=AvailableSediment/ProposedOutflow
      DO N=1,2
        DO Side=iFront,iBehind
          NN=FaceNN(N,Side);FY=FaceNY(N,Side);FX=FaceNX(N,Side)
          AcceptedOutflow=AZMAX1(OutwardSign(N,Side)*SedErosionM(N,NN,FY,FX))*SedimentScale
          SedErosionM(N,NN,FY,FX)=OutwardSign(N,Side)*AcceptedOutflow
          TERSED_col(NY,NX)=TERSED_col(NY,NX)-AcceptedOutflow
          if(InternalFace(N,Side))then
            ! Credit the identical accepted transfer to its receiving cell.
            TERSED_col(FY,FX)=TERSED_col(FY,FX)+AcceptedOutflow
          endif
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  ! Only accepted faces enter the hourly sediment and material-transfer ledger.
  cumSed_Eros_2D=cumSed_Eros_2D+SedErosionM
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      SED_col(NY,NX)=SED_col(NY,NX)+RDTSED_col(NY,NX)+TERSED_col(NY,NX)
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SedimentTransportM
!------------------------------------------------------------------------------------------
  subroutine XBoundSedTranspM(M,NY,NX,NHW,NHE,NVN,NVS,N1,N2)

  implicit none

  character(len=*), parameter :: subname='XBoundSedTranspM'
  integer, intent(in) :: M,N1,N2,NY,NX,NHW,NHE,NVN,NVS
!     LOCATE EXTERNAL BOUNDARIES
!
  integer :: N,NN,N4,N5,N4B,N5B
  integer :: M1,M2,M4,M5
  real(r8) :: RCHQF,FERM,XN

  call PrintInfo('beg '//subname)
  D9580: DO  N=1,2
    D9575: DO  NN=1,2
      IF(N.EQ.iWestEastDirection)THEN
        N4=NX+1;N5=NY
        N4B=NX-1;N5B=NY
        IF(NN.EQ.iFront)THEN 
          IF(NX.EQ.NHE)THEN  !eastern boundary
            M1=NX;M2=NY
            M4=NX+1;M5=NY
            XN    = -1.0_r8
            RCHQF = RechargEastSurf_col(M2,M1)
          ELSE
            cycle
          ENDIF
        ELSEIF(NN.EQ.iBehind)THEN
          IF(NX.EQ.NHW)THEN  !western boundary
            M1=NX;M2=NY
            M4=NX;M5=NY
            XN    = 1.0_r8
            RCHQF = RechargWestSurf_col(M5,M4)
          ELSE
            cycle
          ENDIF
        ENDIF
      ELSEIF(N.EQ.iNorthSouthDirection)THEN
        N4  = NX;N5  = NY+1
        N4B = NX;N5B = NY-1
        IF(NN.EQ.iFront)THEN          
          IF(NY.EQ.NVS)THEN !southern boundary
            M1 = NX;M2 = NY
            M4 = NX;M5 = NY+1
            XN = -1.0_r8
            RCHQF=RechargSouthSurf_col(M2,M1)
          ELSE
            cycle
          ENDIF
        ELSEIF(NN.EQ.iBehind)THEN          
          IF(NY.EQ.NVN)THEN !north boundary
            M1 = NX;M2 = NY
            M4 = NX;M5 = NY
            XN = 1.0_r8
            RCHQF=RechargNorthSurf_col(M5,M4)
          ELSE
            cycle
          ENDIF
        ENDIF
      ENDIF
!
!     SEDIMENT TRANSPORT ACROSS BOUNDARY FROM BOUNDARY RUNOFF
!     IN 'WATSUB' TIMES BOUNDARY SEDIMENT CONCENTRATION IN
!     SURFACE WATER
!
      IF(.not.XGridRunoffFlag_2DH(NN,N,N2,N1) .OR. isclose(RCHQF,0._r8) &
        .OR. BaseErosionRate(N2,N1).LE.ZEROS(N2,N1))THEN
        SedErosionM(N,NN,M5,M4)=0._r8
      ELSE
        IF(SurfRunoffPotentM_col(M,N2,N1).GT.ZEROS(N2,N1))THEN
          IF((NN.EQ.iFront .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
            .OR. (NN.EQ.iBehind .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
            FERM                       = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)/SurfRunoffPotentM_col(M,N2,N1)
            SedErosionM(N,NN,M5,M4)    = BaseErosionRate(N2,N1)*FERM
          ELSEIF((NN.EQ.iBehind.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
            .OR.(NN.EQ.iFront.AND.QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
            SedErosionM(N,NN,M5,M4)=0._r8
          ELSE
            SedErosionM(N,NN,M5,M4)=0._r8
          ENDIF
        ENDIF
      ENDIF
    ENDDO D9575
  ENDDO D9580
  call PrintInfo('end '//subname)
  end subroutine XBoundSedTranspM
!------------------------------------------------------------------------------------------

  subroutine CalcErosionMaterialDenom(NY,NX,NHW,NHE,NVN,NVS,MaterialDenom)
  implicit none
  integer, intent(in) :: NY,NX,NHW,NHE,NVN,NVS
  real(r8), intent(out) :: MaterialDenom
  character(len=*), parameter :: subname='CalcErosionMaterialDenom'
  real(r8) :: OutgoingSediment

  call PrintInfo('beg '//subname)
  ! Sum this donor's hourly outgoing sediment magnitudes across all faces.
  ! Internal transfers use receiving slots; external W/N exports are negative
  ! in the donor slot. Incoming sediment does not offset the material budget.
  OutgoingSediment=0._r8
  if(NX<NHE)then
    OutgoingSediment=OutgoingSediment+AZMAX1(cumSed_Eros_2D(iWestEastDirection,iBehind,NY,NX+1))
  elseif(.not.column_mode)then
    OutgoingSediment=OutgoingSediment+AZMAX1(cumSed_Eros_2D(iWestEastDirection,iFront,NY,NX+1))
  endif
  if(NX>NHW)then
    OutgoingSediment=OutgoingSediment+AZMAX1(cumSed_Eros_2D(iWestEastDirection,iFront,NY,NX-1))
  elseif(.not.column_mode)then
    OutgoingSediment=OutgoingSediment+AZMAX1(-cumSed_Eros_2D(iWestEastDirection,iBehind,NY,NX))
  endif
  if(NY<NVS)then
    OutgoingSediment=OutgoingSediment+AZMAX1(cumSed_Eros_2D(iNorthSouthDirection,iBehind,NY+1,NX))
  elseif(.not.column_mode)then
    OutgoingSediment=OutgoingSediment+AZMAX1(cumSed_Eros_2D(iNorthSouthDirection,iFront,NY+1,NX))
  endif
  if(NY>NVN)then
    OutgoingSediment=OutgoingSediment+AZMAX1(cumSed_Eros_2D(iNorthSouthDirection,iFront,NY-1,NX))
  elseif(.not.column_mode)then
    OutgoingSediment=OutgoingSediment+AZMAX1(-cumSed_Eros_2D(iNorthSouthDirection,iBehind,NY,NX))
  endif
  ! A common denominator keeps the sum of exported material fractions <= 1,
  ! preserving each face's sign and relative share. Protect empty soil too.
  ! This limits soil-pool withdrawals; it does not track suspended composition
  ! or alter the separately accepted suspended-sediment transport ledger.
  MaterialDenom=max(ZERO,SoilMicPMassLayerMX(NY,NX),OutgoingSediment)
  call PrintInfo('end '//subname)
  end subroutine CalcErosionMaterialDenom
!------------------------------------------------------------------------------------------
  subroutine LateralXGridSedmentFlux(NHW, NHE,NVN,NVS)
  implicit none

  character(len=*), parameter :: subname='LateralXGridSedmentFlux'
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NGL,NTX,NTP,NE,idom
  integer :: K,NN,N,NO,M,NY,NX
  integer :: N1,N2,N4,N5,N4B,N5B
  real(r8) :: FSEDER,MaterialDenom

  call PrintInfo('beg '//subname)
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
        call CalcErosionMaterialDenom(NY,NX,NHW,NHE,NVN,NVS,MaterialDenom)
        N1=NX
        N2=NY
        DO  N=1,2
          DO  NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              IF(NX.EQ.NHE .AND. NN.EQ.iFront .OR. NX.EQ.NHW .AND. NN.EQ.iBehind)THEN
                cycle
              ELSE
                N4  = NX+1
                N5  = NY
                N4B = NX-1
                N5B = NY
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              IF(NY.EQ.NVS .AND. NN.EQ.iFront .OR. NY.EQ.NVN .AND. NN.EQ.iBehind)THEN
                cycle
              ELSE
                N4=NX
                N5=NY+1
                N4B=NX
                N5B=NY-1
              ENDIF
            ENDIF
!
!     FLUXES OF ALL SOLID MATERIALS IN SEDIMENT ARE CALCULATED
!     FROM VALUES OF THEIR CURRENT STATE VARIABLES MULTIPLIED
!     BY THE FRACTION OF THE TOTAL SURFACE LAYER MASS THAT IS
!     TRANSPORTED IN SEDIMENT
!
!     SOIL MINERALS
!
            IF(NN.EQ.iFront)THEN
              FSEDER = cumSed_Eros_2D(N,iBehind,N5,N4)/MaterialDenom

              XSand_Eros_2D(N,iBehind,N5,N4) = FSEDER*SAND_vr(NU_col(N2,N1),N2,N1)
              XSilt_Eros_2D(N,iBehind,N5,N4) = FSEDER*SILT_vr(NU_col(N2,N1),N2,N1)
              XClay_Eros_2D(N,iBehind,N5,N4) = FSEDER*CLAY_vr(NU_col(N2,N1),N2,N1)
!
!     FERTILIZER POOLS
!
!     *ER=sediment flux from erosion.f
!
              XNH4Soil_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_soil_vr(ifert_N_nh4,NU_col(N2,N1),N2,N1)
              XNH3Soil_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_soil_vr(ifert_N_nh3,NU_col(N2,N1),N2,N1)
              XUreaSoil_Eros_molN_2D(N,iBehind,N5,N4) = FSEDER*FertN_mole_soil_vr(ifert_N_urea,NU_col(N2,N1),N2,N1)
              XNO3Soil_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_soil_vr(ifert_N_no3,NU_col(N2,N1),N2,N1)
              XPO4Soil_Eros_molP_2D(N,iBehind,N5,N4) =  FSEDER*FertP_mole_soil_vr(NU_col(N2,N1),N2,N1)

              XNH4Band_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_Band_vr(ifert_N_nh4_band,NU_col(N2,N1),N2,N1)
              XNH3Band_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_Band_vr(ifert_N_nh3_band,NU_col(N2,N1),N2,N1)
              XUreaBand_Eros_molN_2D(N,iBehind,N5,N4) = FSEDER*FertN_mole_Band_vr(ifert_N_urea_band,NU_col(N2,N1),N2,N1)
              XNO3Band_Eros_molN_2D(N,iBehind,N5,N4)  = FSEDER*FertN_mole_Band_vr(ifert_N_no3_band,NU_col(N2,N1),N2,N1)
              XPO4Band_Eros_molP_2D(N,iBehind,N5,N4)  = FSEDER*FertP_mole_band_vr(NU_col(N2,N1),N2,N1)
!
!     EXCHANGEABLE CATIONS AND ANIONS
!     sediment code
!
              DO NTX=idx_beg,idx_end
                trcx_Eros_mol_2D(NTX,N,iBehind,N5,N4)=FSEDER*trcx_solml_vr(NTX,NU_col(N2,N1),N2,N1)
              ENDDO
!
!     PRECIPITATES
!
!     sediment code
!
              DO NTP=idsp_beg,idsp_end
                trcp_Eros_mol_2D(NTP,N,iBehind,N5,N4)   =FSEDER*trcp_saltpml_vr(NTP,NU_col(N2,N1),N2,N1)
              ENDDO
!
!     ORGANIC MATTER
!
              DO  K=1,jcplx
                DO NO=1,NumMicbHFunGrupsPerCmplx
                  DO NGL=JGniH(NO),JGnfH(NO)
                    DO M=1,nlbiomcp
                      associate(MID => micpar%get_micb_id(M,NGL))
                        DO NE=1,NumPlantChemElms
                          OMEERhetr_2D(NE,MID,K,N,iBehind,N5,N4)=FSEDER*mBiomeHeter_vr(NE,MID,K,NU_col(N2,N1),N2,N1)
                        ENDDO
                      end associate
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO

              DO NO=1,NumMicbAFunGrupsPerCmplx
                DO NGL=JGniA(NO),JGnfA(NO)
                  DO M=1,nlbiomcp
                    associate(MID => micpar%get_micb_id(M,NGL))
                      DO NE=1,NumPlantChemElms
                        OMEERauto_2D(NE,MID,N,iBehind,N5,N4)=FSEDER*mBiomeAutor_vr(NE,MID,NU_col(N2,N1),N2,N1)
                      ENDDO
                    end associate
                  ENDDO
                ENDDO
              ENDDO

              DO  K=1,jcplx
                DO  M=1,ndbiomcp
                  DO NE=1,NumPlantChemElms
                    OMBioResdu_Eros_2D(NE,M,K,N,iBehind,N5,N4)=FSEDER*OMBioResdu_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                  ENDDO
                ENDDO
                DO idom=idom_beg,idom_end
                  SorbedOM_Eros_2D(idom,K,N,iBehind,N5,N4)=FSEDER*SorbedOM_vr(idom,K,NU_col(N2,N1),N2,N1)
                ENDDO
                DO  M=1,jsken
                  DO NE=1,NumPlantChemElms
                    SolidOM_Eros_2D(NE,M,K,N,iBehind,N5,N4)=FSEDER*SolidOM_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                  ENDDO
                  SolidOMAct_Eros_2D(M,K,N,iBehind,N5,N4)=FSEDER*SolidOMAct_vr(M,K,NU_col(N2,N1),N2,N1)
                ENDDO
              ENDDO
            ENDIF
            IF(NN.EQ.iBehind)THEN
              IF(N4B.GT.0.AND.N5B.GT.0)THEN
                FSEDER=cumSed_Eros_2D(N,iFront,N5B,N4B)/MaterialDenom
                XSand_Eros_2D(N,iFront,N5B,N4B) = FSEDER*SAND_vr(NU_col(N2,N1),N2,N1)
                XSilt_Eros_2D(N,iFront,N5B,N4B) = FSEDER*SILT_vr(NU_col(N2,N1),N2,N1)
                XClay_Eros_2D(N,iFront,N5B,N4B) = FSEDER*CLAY_vr(NU_col(N2,N1),N2,N1)
!
!     FERTILIZER POOLS
!
                XNH4Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_soil_vr(ifert_N_nh4,NU_col(N2,N1),N2,N1)
                XNH3Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_soil_vr(ifert_N_nh3,NU_col(N2,N1),N2,N1)
                XUreaSoil_Eros_molN_2D(N,iFront,N5B,N4B) = FSEDER*FertN_mole_soil_vr(ifert_N_urea,NU_col(N2,N1),N2,N1)
                XNO3Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_soil_vr(ifert_N_no3,NU_col(N2,N1),N2,N1)
                XPO4Soil_Eros_molP_2D(N,iFront,N5B,N4B)  = FSEDER*FertP_mole_soil_vr(NU_col(N2,N1),N2,N1)
                XNH4Band_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_Band_vr(ifert_N_nh4_band,NU_col(N2,N1),N2,N1)
                XNH3Band_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_Band_vr(ifert_N_nh3_band,NU_col(N2,N1),N2,N1)
                XUreaBand_Eros_molN_2D(N,iFront,N5B,N4B) = FSEDER*FertN_mole_Band_vr(ifert_N_urea_band,NU_col(N2,N1),N2,N1)
                XNO3Band_Eros_molN_2D(N,iFront,N5B,N4B)  = FSEDER*FertN_mole_Band_vr(ifert_N_no3_band,NU_col(N2,N1),N2,N1)
                XPO4Band_Eros_molP_2D(N,iFront,N5B,N4B)  = FSEDER*FertP_mole_band_vr(NU_col(N2,N1),N2,N1)

!
!     EXCHANGEABLE CATIONS AND ANIONS
!
!     sediment code
!
                DO NTX=idx_beg,idx_end
                  trcx_Eros_mol_2D(NTX,N,iFront,N5B,N4B)=FSEDER*trcx_solml_vr(NTX,NU_col(N2,N1),N2,N1)
                ENDDO
!
!     PRECIPITATES
!
!     sediment code
                DO NTP=idsp_beg,idsp_end
                  trcp_Eros_mol_2D(NTP,N,iFront,N5B,N4B)=FSEDER*trcp_saltpml_vr(NTP,NU_col(N2,N1),N2,N1)
                ENDDO
!
!     ORGANIC MATTER
!
                DO  K=1,jcplx
                  DO  NO=1,NumMicbHFunGrupsPerCmplx
                    DO NGL=JGniH(NO),JGnfH(NO)
                      DO  M=1,nlbiomcp
                        associate(MID => micpar%get_micb_id(M,NGL))
                          DO NE=1,NumPlantChemElms
                            OMEERhetr_2D(NE,MID,K,N,iFront,N5B,N4B)=FSEDER*mBiomeHeter_vr(NE,MID,K,NU_col(N2,N1),N2,N1)
                          ENDDO
                        end associate
                      enddo
                    enddo
                  ENDDO
                ENDDO
                DO  NO=1,NumMicbAFunGrupsPerCmplx
                  DO NGL=JGniA(NO),JGnfA(NO)
                    DO  M=1,nlbiomcp
                      associate(MID => micpar%get_micb_id(M,NGL))
                        DO NE=1,NumPlantChemElms
                          OMEERauto_2D(NE,MID,N,iFront,N5B,N4B)=FSEDER*mBiomeAutor_vr(NE,MID,NU_col(N2,N1),N2,N1)
                        ENDDO
                      end associate
                    enddo
                  enddo
                ENDDO

                DO  K=1,jcplx
                  DO  M=1,ndbiomcp
                    DO NE=1,NumPlantChemElms
                      OMBioResdu_Eros_2D(NE,M,K,N,iFront,N5B,N4B)=FSEDER*OMBioResdu_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                    ENDDO
                  ENDDO
                  DO idom=idom_beg,idom_end
                    SorbedOM_Eros_2D(idom,K,N,iFront,N5B,N4B)=FSEDER*SorbedOM_vr(idom,K,NU_col(N2,N1),N2,N1)
                  ENDDO
                  DO  M=1,jsken
                    SolidOMAct_Eros_2D(M,K,N,iFront,N5B,N4B)=FSEDER*SolidOMAct_vr(M,K,NU_col(N2,N1),N2,N1)
                    DO NE=1,NumPlantChemElms
                      SolidOM_Eros_2D(NE,M,K,N,iFront,N5B,N4B)=FSEDER*SolidOM_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                    ENDDO  
                  ENDDO
                ENDDO
              ELSE
                XSand_Eros_2D(N,iFront,N5B,N4B) = 0._r8
                XSilt_Eros_2D(N,iFront,N5B,N4B) = 0._r8
                XClay_Eros_2D(N,iFront,N5B,N4B) = 0._r8
!
!     FERTILIZER POOLS
!
                XNH4Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XNH3Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XUreaSoil_Eros_molN_2D(N,iFront,N5B,N4B) = 0._r8
                XNO3Soil_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XPO4Soil_Eros_molP_2D(N,iFront,N5B,N4B)  = 0._r8
                XNH4Band_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XNH3Band_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XUreaBand_Eros_molN_2D(N,iFront,N5B,N4B) = 0._r8
                XNO3Band_Eros_molN_2D(N,iFront,N5B,N4B)  = 0._r8
                XPO4Band_Eros_molP_2D(N,iFront,N5B,N4B)= 0._r8
!
!     EXCHANGEABLE CATIONS AND ANIONS
!
                trcx_Eros_mol_2D(idx_beg:idx_end,N,iFront,N5B,N4B)=0._r8
!
!     PRECIPITATES
!
                trcp_Eros_mol_2D(idsp_beg:idsp_end,N,iFront,N5B,N4B)=0._r8
!
!     ORGANIC MATTER
!
                DO  K=1,jcplx
                  DO  NO=1,NumMicbHFunGrupsPerCmplx
                    DO NGL=JGniH(NO),JGnfH(NO)
                      DO  M=1,nlbiomcp
                        associate(MID => micpar%get_micb_id(M,NGL))
                          DO NE=1,NumPlantChemElms
                            OMEERhetr_2D(NE,MID,K,N,iFront,N5B,N4B)=0._r8
                          ENDDO
                        end associate
                      enddo
                    ENDDO
                  enddo
                ENDDO

                DO  NO=1,NumMicbAFunGrupsPerCmplx
                  DO NGL=JGniA(NO),JGnfA(NO)
                    DO  M=1,nlbiomcp
                      associate(MID => micpar%get_micb_id(M,NGL))
                        DO NE=1,NumPlantChemElms
                          OMEERauto_2D(NE,MID,N,iFront,N5B,N4B)=0._r8
                        ENDDO
                      end associate
                    enddo
                  ENDDO
                enddo

                DO  K=1,jcplx
                  DO  M=1,ndbiomcp
                    DO NE=1,NumPlantChemElms
                      OMBioResdu_Eros_2D(NE,M,K,N,iFront,N5B,N4B)=0._r8
                    ENDDO  
                  ENDDO
                  DO idom=idom_beg,idom_end
                    SorbedOM_Eros_2D(idom,K,N,iFront,N5B,N4B)=0._r8
                  ENDDO
                  DO  M=1,jsken
                    SolidOMAct_Eros_2D(M,K,N,iFront,N5B,N4B)=0._r8
                    DO NE=1,NumPlantChemElms
                      SolidOM_Eros_2D(NE,M,K,N,iFront,N5B,N4B)=0._r8
                    ENDDO  
                  ENDDO
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine LateralXGridSedmentFlux
!------------------------------------------------------------------------------------------

  subroutine XBoundErosionFlux(NHW,NHE,NVN,NVS)
  implicit none

  character(len=*), parameter :: subname='XBoundErosionFlux'

  integer, intent(in) :: NHW,NHE,NVN,NVS

  real(r8) :: RCHQF,FSEDER,MaterialDenom
  integer :: NGL
  integer :: NY,NX,NE,idom
  integer :: N,NN,K,NO,M,NTX,NTP
  integer :: N1,N2,N4,N5
  integer :: M1,M2,M4,M5
  real(r8) :: XN

  call PrintInfo('beg '//subname)
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      IF(SoilBulkDensity_vr(NU_col(NY,NX),NY,NX).GT.ZERO)THEN
        call CalcErosionMaterialDenom(NY,NX,NHW,NHE,NVN,NVS,MaterialDenom)
        N1=NX
        N2=NY
        D8980: DO  N=1,2
          D8975: DO  NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              N4=NX+1
              N5=NY
              IF(NN.EQ.iFront)THEN
                IF(NX.EQ.NHE)THEN
                  M1    = NX
                  M2    = NY
                  M4    = NX+1
                  M5    = NY
                  XN    = -1.0_r8
                  RCHQF = RechargEastSurf_col(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN
                IF(NX.EQ.NHW)THEN
                  M1 = NX
                  M2 = NY
                  M4 = NX
                  M5 = NY
                  XN = 1.0_r8
                  RCHQF=RechargWestSurf_col(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              N4=NX
              N5=NY+1
              IF(NN.EQ.iFront)THEN
                IF(NY.EQ.NVS)THEN
                  M1 = NX
                  M2 = NY
                  M4 = NX
                  M5 = NY+1
                  XN = -1.0_r8
                  RCHQF=RechargSouthSurf_col(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN
                IF(NY.EQ.NVN)THEN
                  M1=NX
                  M2=NY
                  M4=NX
                  M5=NY
                  XN=1.0_r8
                  RCHQF=RechargNorthSurf_col(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ENDIF
            IF(.not.XGridRunoffFlag_2DH(NN,N,N2,N1).OR.isclose(RCHQF,0._r8) &
              .OR. ABS(cumSed_Eros_2D(N,NN,M5,M4)).LE.ZEROS(N2,N1))THEN
              XSand_Eros_2D(N,NN,M5,M4)=0._r8
              XSilt_Eros_2D(N,NN,M5,M4)=0._r8
              XClay_Eros_2D(N,NN,M5,M4)=0._r8
!
!     FERTILIZER POOLS
!
              XNH4Soil_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XNH3Soil_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XUreaSoil_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XNO3Soil_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XPO4Soil_Eros_molP_2D(N,NN,M5,M4)=0._r8
              XNH4Band_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XNH3Band_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XUreaBand_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XNO3Band_Eros_molN_2D(N,NN,M5,M4)=0._r8
              XPO4Band_Eros_molP_2D(N,NN,M5,M4)=0._r8
!
!     EXCHANGEABLE CATIONS AND ANIONS
!
              trcx_Eros_mol_2D(idx_beg:idx_end,N,NN,M5,M4)=0._r8
!
!     PRECIPITATES
!
              trcp_Eros_mol_2D(idsp_beg:idsp_end,N,NN,M5,M4)=0._r8
!
!     ORGANIC MATTER
!
              DO  K=1,jcplx
                DO  NO=1,NumMicbHFunGrupsPerCmplx
                  DO NGL=JGniH(NO),JGnfH(NO)
                    DO  M=1,nlbiomcp
                      associate(MID => micpar%get_micb_id(M,NGL))
                        DO NE=1,NumPlantChemElms
                          OMEERhetr_2D(NE,MID,K,N,NN,M5,M4)=0._r8
                        ENDDO
                      end associate
                    enddo
                  ENDDO
                enddo
              enddo
              DO  NO=1,NumMicbAFunGrupsPerCmplx
                DO NGL=JGniA(NO),JGnfA(NO)
                  DO  M=1,nlbiomcp
                    associate(MID => micpar%get_micb_id(M,NGL))
                      DO NE=1,NumPlantChemElms
                        OMEERauto_2D(NE,MID,N,NN,M5,M4)=0._r8
                      ENDDO
                    end associate
                  enddo
                ENDDO
              enddo
              DO  K=1,jcplx
                DO  M=1,ndbiomcp
                  DO NE=1,NumPlantChemElms
                    OMBioResdu_Eros_2D(NE,M,K,N,NN,M5,M4)=0._r8
                  ENDDO  
                ENDDO
                DO idom=idom_beg,idom_end
                  SorbedOM_Eros_2D(idom,K,N,NN,M5,M4)=0._r8
                ENDDO
                DO  M=1,jsken
                  SolidOMAct_Eros_2D(M,K,N,NN,M5,M4)=0._r8
                  DO NE=1,NumPlantChemElms
                    SolidOM_Eros_2D(NE,M,K,N,NN,M5,M4)=0._r8
                  ENDDO
                ENDDO
              ENDDO
!
!     CALCULATE FRACTION OF SURACE MATERIAL ERODED
!
            ELSE
              FSEDER=cumSed_Eros_2D(N,NN,M5,M4)/MaterialDenom
!
!     SOIL MINERALS
!
              XSand_Eros_2D(N,NN,M5,M4)=FSEDER*SAND_vr(NU_col(N2,N1),N2,N1)
              XSilt_Eros_2D(N,NN,M5,M4)=FSEDER*SILT_vr(NU_col(N2,N1),N2,N1)
              XClay_Eros_2D(N,NN,M5,M4)=FSEDER*CLAY_vr(NU_col(N2,N1),N2,N1)
!
!     FERTILIZER POOLS
!
              XNH4Soil_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_soil_vr(ifert_N_nh4,NU_col(N2,N1),N2,N1)
              XNH3Soil_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_soil_vr(ifert_N_nh3,NU_col(N2,N1),N2,N1)
              XUreaSoil_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_soil_vr(ifert_N_urea,NU_col(N2,N1),N2,N1)
              XNO3Soil_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_soil_vr(ifert_N_no3,NU_col(N2,N1),N2,N1)
              XPO4Soil_Eros_molP_2D(N,NN,M5,M4)=FSEDER*FertP_mole_soil_vr(NU_col(N2,N1),N2,N1)

              XNH4Band_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_Band_vr(ifert_N_nh4_band,NU_col(N2,N1),N2,N1)
              XNH3Band_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_Band_vr(ifert_N_nh3_band,NU_col(N2,N1),N2,N1)
              XUreaBand_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_Band_vr(ifert_N_urea_band,NU_col(N2,N1),N2,N1)
              XNO3Band_Eros_molN_2D(N,NN,M5,M4)=FSEDER*FertN_mole_Band_vr(ifert_N_no3_band,NU_col(N2,N1),N2,N1)
              XPO4Band_Eros_molP_2D(N,NN,M5,M4)=FSEDER*FertP_mole_band_vr(NU_col(N2,N1),N2,N1)
!
!     EXCHANGEABLE CATIONS AND ANIONS
!
              DO NTX=idx_beg,idx_end
                trcx_Eros_mol_2D(NTX,N,NN,M5,M4)=FSEDER*trcx_solml_vr(NTX,NU_col(N2,N1),N2,N1)
              ENDDO
!
!     PRECIPITATES
!
              DO NTP=idsp_beg,idsp_end
                trcp_Eros_mol_2D(NTP,N,NN,M5,M4)=FSEDER*trcp_saltpml_vr(NTP,NU_col(N2,N1),N2,N1)
              ENDDO
!
!     ORGANIC MATTER
!
              DO  K=1,jcplx
                DO NO=1,NumMicbHFunGrupsPerCmplx
                  DO NGL=JGniH(NO),JGnfH(NO)
                    DO M=1,nlbiomcp
                      associate(MID => micpar%get_micb_id(M,NGL))
                        DO NE=1,NumPlantChemElms
                          OMEERhetr_2D(NE,MID,K,N,NN,M5,M4)=FSEDER*mBiomeHeter_vr(NE,MID,K,NU_col(N2,N1),N2,N1)
                        ENDDO
                      end associate
                    ENDDO
                  ENDDO
                ENDDO
              ENDDO
              DO NO=1,NumMicbAFunGrupsPerCmplx
                DO NGL=JGniA(NO),JGnfH(NO)
                  DO M=1,nlbiomcp
                    associate(MID => micpar%get_micb_id(M,NGL))
                      DO NE=1,NumPlantChemElms
                        OMEERauto_2D(NE,MID,N,NN,M5,M4)=FSEDER*mBiomeAutor_vr(NE,MID,NU_col(N2,N1),N2,N1)
                      ENDDO
                    end associate
                  ENDDO
                ENDDO
              ENDDO

              DO  K=1,jcplx
                DO  M=1,ndbiomcp
                  DO NE=1,NumPlantChemElms
                    OMBioResdu_Eros_2D(NE,M,K,N,NN,M5,M4)=FSEDER*OMBioResdu_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                  ENDDO  
                ENDDO
                DO idom=idom_beg,idom_end
                  SorbedOM_Eros_2D(idom,K,N,NN,M5,M4)=FSEDER*SorbedOM_vr(idom,K,NU_col(N2,N1),N2,N1)
                ENDDO
                DO  M=1,jsken
                  SolidOMAct_Eros_2D(M,K,N,NN,M5,M4)=FSEDER*SolidOMAct_vr(M,K,NU_col(N2,N1),N2,N1)
                  DO NE=1,NumPlantChemElms
                    SolidOM_Eros_2D(NE,M,K,N,NN,M5,M4)=FSEDER*SolidOM_vr(NE,M,K,NU_col(N2,N1),N2,N1)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
          ENDDO D8975
        ENDDO D8980
      ENDIF
    ENDDO
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine XBoundErosionFlux

!------------------------------------------------------------------------------------------
  subroutine DestructErosion()
  use abortutils, only : destroy

  implicit none

  character(len=*), parameter :: subname='DestructErosion'

  call PrintInfo('beg '//subname)
  call destroy(SedErosionM)
  call destroy(TERSED_col)
  call destroy(RDTSED_col)
  call destroy(FVOLIM)
  call destroy(FVLWatMicPM)
  call destroy(FracVol4Erosion)
  call destroy(BaseErosionRate)

  call PrintInfo('end '//subname)
  end subroutine DestructErosion
  end module ErosionMod
