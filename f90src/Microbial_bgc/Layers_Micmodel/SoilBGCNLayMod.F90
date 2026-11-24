module SoilBGCNLayMod
!!
! DESCRIPTION:
! codes to do soil biological transformations
!
! USES:
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,       only: endrun
  use minimathmod,      only: safe_adb, AZMAX1,AZERO
  use EcoSiMParDataMod, only: micpar
  use SoilHeatDataType, only: TCS_vr
  use PlantMgmtDataType,only: iDayPlanting_pft
  use DebugToolMod
  use SoilWaterDataType
  use SurfLitterDataType
  use MicrobialDataType
  use NitroPars
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use MicrobeDiagTypes
  use SoilDisturbMod
  use GridConsts
  use SoilBGCDataType
  use EcoSIMCtrlDataType
  use SoilPhysDataType
  use SoilPropertyDataType
  use GridDataType
  use MicBGCMod
  use EcoSIMConfig , only : ndbiomcp => NumDeadMicrbCompts  
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

!
  public :: initNitro
  public :: DownwardMixOM
  contains

!------------------------------------------------------------------------------------------

  subroutine initNitro

  implicit none

  call initNitro1Layer

  end subroutine initNitro

!------------------------------------------------------------------------------------------

  subroutine DownwardMixOM(I,J,L,NY,NX,FracLitrMix)

  implicit none
  integer, intent(in) :: I,J,L,NY,NX
  real(r8), intent(out) :: FracLitrMix
  real(r8) :: FOSCXD
  real(r8) :: ORGRL,ORGRLL
  real(r8) :: OSCXD
  integer :: LL,LN,K
  real(r8) :: DOC_s,DOC_u,ActL,ActLL,ActD
  character(len=*), parameter :: subname='DownwardMixOM'

!     begin_execution
  call PrintInfo('beg '//subname)

  IF(FOSCZ0.GT.ZERO)THEN
!     OMLitrC_vr=total litter C
!     FOSCZ0=rate constant for mixing surface litter
!     FracLitrMix=mixing fraction for surface litter
!     TMicHeterActivity_vr=total active biomass respiration activity
!     VLSoilPoreMicP_vr=soil layer volume
!     OSCXD=mixing required for equilibrating litter concentration
!     FOSCXD=mixing fraction for equilibrating subsurface litter
!     FracLitrMix=mixing fraction for subsurface litter
!
    IF(L.LT.NL_col(NY,NX))THEN
!     get mixing rate
      IF(L.EQ.0)THEN
        LL=NU_col(NY,NX)
        DOC_s=0._r8;DOC_u=0._r8
        DO K=1,micpar%NumOfLitrCmplxs
          if(.not.micpar%is_finelitter(K))cycle
          DOC_s=DOC_s+DOM_MicP_vr(idom_doc,K,L,NY,NX)
          DOC_u=DOC_u+DOM_MicP_vr(idom_doc,K,LL,NY,NX)
        ENDDO
        ActL  = DOC_s
        ActLL = DOC_u
        ActD  = (ActL*OMLitrC_vr(LL,NY,NX)-ActLL*OMLitrC_vr(L,NY,NX))/(OMLitrC_vr(L,NY,NX)+OMLitrC_vr(LL,NY,NX))

        IF(OMLitrC_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
          FracLitrMix=FOSCZ0/OMLitrC_vr(L,NY,NX)*TMicHeterActivity_vr(L,NY,NX)
        ELSE
          FracLitrMix=0.0_r8
        ENDIF

      ELSE
        D1100: DO LN=L+1,NL_col(NY,NX)
          IF(VLSoilPoreMicP_vr(LN,NY,NX).GT.ZEROS2(NY,NX))THEN
            LL=LN
            exit
          ENDIF
        ENDDO D1100
        ORGRL  = AZMAX1(OMLitrC_vr(L,NY,NX))
        ORGRLL = AZMAX1(OMLitrC_vr(LL,NY,NX))
        OSCXD  = (ORGRL*VGeomLayer_vr(LL,NY,NX)-ORGRLL*VGeomLayer_vr(L,NY,NX))/(VGeomLayer_vr(L,NY,NX)+VGeomLayer_vr(LL,NY,NX))

        IF(OSCXD.GT.0.0_r8 .AND. OMLitrC_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
         !mass gradient pointing from LL to L
          FOSCXD=OSCXD/OMLitrC_vr(L,NY,NX)
        ELSEIF(OSCXD.LT.0.0_r8 .AND. OMLitrC_vr(LL,NY,NX).GT.ZEROS(NY,NX))THEN
          !mass gradient pointing from layer L to LL
          FOSCXD=OSCXD/OMLitrC_vr(LL,NY,NX)
        ELSE
          FOSCXD=0.0_r8
        ENDIF

        IF(VGeomLayer_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          FracLitrMix=FOSCZL*FOSCXD*TMicHeterActivity_vr(L,NY,NX)/VGeomLayer_vr(L,NY,NX)  
        ELSE
          FracLitrMix=0.0_r8
        ENDIF
      ENDIF

!     apply mixing
      call ApplyVerticalMix(FracLitrMix,L,LL,NY,NX)
    ENDIF

  ENDIF
  call PrintInfo('end '//subname)

  end subroutine DownwardMixOM
!------------------------------------------------------------------------------------------

  subroutine ApplyVerticalMix(FracLitrMix,L,LL,NY,NX)
  use EcoSIMConfig, only :   NumMicbHFunGrupsPerCmplx  
  implicit none
  real(r8), intent(in) :: FracLitrMix
  integer, intent(in) :: NY,NX        !horizontal location of the grid
  integer, intent(in) :: L            !source grid
  integer, intent(in) :: LL           !destination grid

  real(r8) :: OMEXS,ORMXS
  real(r8) :: OQMXS,OQMHXS,OHMXS
  real(r8) :: OSMXS,OSAXS
  integer :: K,M,N,NGL,MID,NE,L1
  character(len=*), parameter :: subname='ApplyVerticalMix'
!     begin_execution
  call PrintInfo('beg '//subname)

  IF(FracLitrMix.GT.0.0_r8)THEN
    L1=L
  ELSE
    L1=LL
  ENDIF
  !how about heat, water flux and capacity?
! only downward mixing is considered
! upward mixing is yet to be considered, even though FracLitrMix could be negative

  IF(FracLitrMix.GT.ZERO)THEN
    !mix microbial biomass
    D7971: DO K=1,micpar%NumOfLitrCmplxs
      if(.not.micpar%is_finelitter(K))cycle
      D7961: DO N=1,NumMicbHFunGrupsPerCmplx
        DO NGL=JGniH(N),JGnfH(N)
          D7962: DO M=1,micpar%nlbiomcp
            MID=micpar%get_micb_id(M,NGL)            
            DO NE=1,NumPlantChemElms            
              OMEXS                             = FracLitrMix*AZMAX1(mBiomeHeter_vr(NE,MID,K,L1,NY,NX))
              mBiomeHeter_vr(NE,MID,K,L,NY,NX)  = mBiomeHeter_vr(NE,MID,K,L,NY,NX)-OMEXS
              mBiomeHeter_vr(NE,MID,K,LL,NY,NX) = mBiomeHeter_vr(NE,MID,K,LL,NY,NX)+OMEXS
            ENDDO
          ENDDO D7962
        ENDDO
      ENDDO D7961
    ENDDO D7971

    !mix microbial residual
    D7901: DO K=1,micpar%NumOfLitrCmplxs
      if(.not.micpar%is_finelitter(K))cycle
      !mix fine litter
      D7941: DO M=1,micpar%ndbiomcp

        DO NE=1,NumPlantChemElms
          ORMXS                          = FracLitrMix*AZMAX1(OMBioResdu_vr(NE,M,K,L1,NY,NX))
          OMBioResdu_vr(NE,M,K,L,NY,NX)  = OMBioResdu_vr(NE,M,K,L,NY,NX)-ORMXS
          OMBioResdu_vr(NE,M,K,LL,NY,NX) = OMBioResdu_vr(NE,M,K,LL,NY,NX)+ORMXS
        ENDDO
      ENDDO D7941

      !mix dissolved organic matter

      DO NE=idom_beg,idom_end
        OQMXS  = FracLitrMix*AZMAX1(DOM_MicP_vr(NE,K,L1,NY,NX))
        OQMHXS = FracLitrMix*AZMAX1(DOM_MacP_vr(NE,K,L1,NY,NX))
        OHMXS  = FracLitrMix*AZMAX1(SorbedOM_vr(NE,K,L1,NY,NX))

        DOM_MicP_vr(NE,K,L,NY,NX) = DOM_MicP_vr(NE,K,L,NY,NX)-OQMXS
        DOM_MacP_vr(NE,K,L,NY,NX) = DOM_MacP_vr(NE,K,L,NY,NX)-OQMHXS
        SorbedOM_vr(NE,K,L,NY,NX) = SorbedOM_vr(NE,K,L,NY,NX)-OHMXS

        DOM_MicP_vr(NE,K,LL,NY,NX) = DOM_MicP_vr(NE,K,LL,NY,NX)+OQMXS
        DOM_MacP_vr(NE,K,LL,NY,NX) = DOM_MacP_vr(NE,K,LL,NY,NX)+OQMHXS
        SorbedOM_vr(NE,K,LL,NY,NX) = SorbedOM_vr(NE,K,LL,NY,NX)+OHMXS

      ENDDO

      !mix solid organic matter
      D7931: DO M=1,jsken
        DO NE=1,NumPlantChemElms
          OSMXS                       = FracLitrMix*AZMAX1(SolidOM_vr(NE,M,K,L1,NY,NX))
          SolidOM_vr(NE,M,K,L,NY,NX)  = SolidOM_vr(NE,M,K,L,NY,NX)-OSMXS
          SolidOM_vr(NE,M,K,LL,NY,NX) = SolidOM_vr(NE,M,K,LL,NY,NX)+OSMXS
        ENDDO
        OSAXS                       = FracLitrMix*AZMAX1(SolidOMAct_vr(M,K,L1,NY,NX))
        SolidOMAct_vr(M,K,L,NY,NX)  = SolidOMAct_vr(M,K,L,NY,NX)-OSAXS
        SolidOMAct_vr(M,K,LL,NY,NX) = SolidOMAct_vr(M,K,LL,NY,NX)+OSAXS
      ENDDO D7931
    ENDDO D7901
  ENDIF
  call PrintInfo('end '//subname)

  end subroutine ApplyVerticalMix


end module SoilBGCNLayMod
