module InitSOMBGCMOD
  use data_kind_mod,    only : r8 => DAT_KIND_R8
  use EcoSIMConfig,     only : nlbiomcp => NumLiveMicrbCompts, ndbiomcp=> NumDeadMicrbCompts, NumMicbAFunGrupsPerCmplx,NumMicbHFunGrupsPerCmplx
  use MicrobialDiagMod, only : sumorgmlayl,sumLitrOMLayL, sumMicBiomLayL
  use minimathmod,      only : AZMAX1,safe_adb
  use EcoSiMParDataMod, only : micpar
  use EcosimBGCFluxType, only : CumDryDepoOM_col
  use DebugToolMod
  use MicrobialDataType
  use SOMDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use EcosimConst
  use SurfLitterDataType
  use SoilPropertyDataType
  use GridDataType
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__
  real(r8), allocatable :: CORGCX(:)  !C concentations from OM complexes, gC (Mg soil)^-1
  real(r8), allocatable :: CORGNX(:)  !N concentations from OM complexes, gN (Mg soil)-1
  real(r8), allocatable :: CORGPX(:)  !P concentations from OM complexes, gP (Mg soil)^-1

  public :: InitSOMVars
  public :: InitSOMProfile
  public :: InitSOMConsts
  public :: InitSOMBGC
  public :: DestructSOMBGC
  public :: MicrobeByLitterFall
  public :: InoculateCyanoBacter
  contains
!------------------------------------------------------------------------------------------

  subroutine InitSOMBGC()
  implicit none

  call MicPar%Init()

  JGniH  => micpar%JGniH
  JGnfH  => micpar%JGnfH
  JGniA  => micpar%JGniA
  JGnfA  => micpar%JGnfA
  NumMicrobAutoTrophCmplx = micpar%NumMicrobAutoTrophCmplx
  NumHetetr1MicCmplx      = micpar%NumHetetr1MicCmplx
  NumLiveHeterBioms       = micpar%NumLiveHeterBioms
  NumLiveAutoBioms        = micpar%NumLiveAutoBioms
  allocate(CORGCX(1:jcplx))
  allocate(CORGNX(1:jcplx))
  allocate(CORGPX(1:jcplx))
  end subroutine InitSOMBGC

!------------------------------------------------------------------------------------------

  subroutine InitSOMConsts
  implicit none

  call MicPar%SetPars

  end subroutine InitSOMConsts
!------------------------------------------------------------------------------------------

  subroutine InitSOMVars(L,NY,NX,FCX)
  use NitroPars, only : CyanoInocC
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: FCX

  real(r8), parameter :: DCKR=0.25_r8
  real(r8), parameter :: DCKM=2.5E+04_r8

  integer :: K,M,KK,N,NGL,NN,KL
  real(r8) :: OC,ON,OP,X
  real(r8) :: OME1(1:NumPlantChemElms)
  real(r8) :: FOSCI,FOSNI,FOSPI
  real(r8) :: RNT,RPT
  real(r8) :: FRNT,FRPT
  real(r8) :: TOSNI,TOSPI,TOSCI
  real(r8) :: OSCI(1:jcplx),OSNI(1:jcplx),OSPI(1:jcplx)
  real(r8) :: OSCM(1:jcplx)
  real(r8) :: TOSCK(1:jcplx),TOSNK(1:jcplx),TOSPK(1:jcplx)
  real(r8) :: CNOSCT(1:jcplx), CPOSCT(1:jcplx)
  real(r8) :: OSCX(1:jcplx)
  real(r8) :: OSNX(1:jcplx)
  real(r8) :: OSPX(1:jcplx)
  real(r8) :: litrOM(NumPlantChemElms)
  real(r8) :: ORGM(NumPlantChemElms)
  real(r8) :: tglds,tCyanoInocC
  integer  :: MID,NE
  ! begin_execution

  associate(                                                   &
    rNCOMC_ave              => micpar%rNCOMC_ave,              &
    rPCOMC_ave              => micpar%rPCOMC_ave,              &
    nlbiomcp                => micpar%nlbiomcp,                &
    k_humus                 => micpar%k_humus,                 &
    mid_HeterAerobBacter    => micpar%mid_HeterAerobBacter,    &
    mid_HeterMixtCynoBacter => micpar%mid_HeterMixtCynoBacter, &
    OHCK                    => micpar%OHCK,                    &
    OMCK                    => micpar%OMCK,                    &
    OQCK                    => micpar%OQCK,                    &
    ORCK                    => micpar%ORCK,                    &
    ORCI                    => micpar%ORCI,                    &
    OMCI                    => micpar%OMCI,                    &
    CNRH                    => micpar%CNRH,                    &
    CPRH                    => micpar%CPRH,                    &
    CNOFC                   => micpar%CNOFC,                   &
    CPOFC                   => micpar%CPOFC,                   &
    OMCF                    => micpar%OMCF,                    &
    OMCA                    => micpar%OMCA                     &
  )

  D975: DO K=1,micpar%NumOfLitrCmplxs
    CNOSCT(K)=0.0_r8
    CPOSCT(K)=0.0_r8
    IF(RSC_vr(K,L,NY,NX).GT.ZEROS(NY,NX))THEN
      RNT=0.0_r8
      RPT=0.0_r8
      D970: DO M=1,jsken
        RNT = RNT+RSC_vr(K,L,NY,NX)*CFOSC_vr(M,K,L,NY,NX)*CNOFC(M,K)
        RPT = RPT+RSC_vr(K,L,NY,NX)*CFOSC_vr(M,K,L,NY,NX)*CPOFC(M,K)
      ENDDO D970
      FRNT=RSN_vr(K,L,NY,NX)/RNT
      FRPT=RSP_vr(K,L,NY,NX)/RPT
      D960: DO M=1,jsken
        CNOSC_vr(M,K,L,NY,NX) = CNOFC(M,K)*FRNT
        CPOSC_vr(M,K,L,NY,NX) = CPOFC(M,K)*FRPT
        CNOSCT(K)          = CNOSCT(K)+CFOSC_vr(M,K,L,NY,NX)*CNOSC_vr(M,K,L,NY,NX)
        CPOSCT(K)          = CPOSCT(K)+CFOSC_vr(M,K,L,NY,NX)*CPOSC_vr(M,K,L,NY,NX)
      ENDDO D960
    ELSE
      D965: DO M=1,jsken
        CNOSC_vr(M,K,L,NY,NX)=CNRH(K)
        CPOSC_vr(M,K,L,NY,NX)=CPRH(K)
      ENDDO D965
      CNOSCT(K)=CNRH(K)
      CPOSCT(K)=CPRH(K)
    ENDIF
  ENDDO D975

  D990: DO K=micpar%NumOfLitrCmplxs+1,jcplx
    CNOSCT(K)=0.0_r8
    CPOSCT(K)=0.0_r8
    IF(CORGCX(K).GT.ZERO)THEN
      D985: DO M=1,jsken
        CNOSC_vr(M,K,L,NY,NX) = CORGNX(K)/CORGCX(K)
        CPOSC_vr(M,K,L,NY,NX) = CORGPX(K)/CORGCX(K)
        CNOSCT(K)          = CNOSCT(K)+CFOSC_vr(M,K,L,NY,NX)*CNOSC_vr(M,K,L,NY,NX)
        CPOSCT(K)          = CPOSCT(K)+CFOSC_vr(M,K,L,NY,NX)*CPOSC_vr(M,K,L,NY,NX)
      ENDDO D985
    ELSE
      D980: DO M=1,jsken
        CNOSC_vr(M,K,L,NY,NX)=CNRH(K)
        CPOSC_vr(M,K,L,NY,NX)=CPRH(K)
      ENDDO D980
      CNOSCT(K)=CNRH(K)
      CPOSCT(K)=CPRH(K)
    ENDIF
!    write(*,*)K,L,CNOSCT(K),CPOSCT(K),jcplx
  ENDDO D990
!
!     MICROBIAL BIOMASS,RESIDUE, DOC, ADSORBED
!
!     OSCI,OSNI,OSPI=initial SOC,SON,SOP mass in each complex (g)
!     OMCK,ORCK,OQCK,OHCK=fractions of SOC in biomass,litter,DOC adsorbed C
!     OSCM=total biomass in each complex (g)
!     DCKR,DCKM=parameters to initialize microbial biomass from SOC
!
  TOSCI=0.0_r8
  TOSNI=0.0_r8
  TOSPI=0.0_r8
  D995: DO K=1,jcplx
    IF(L.EQ.0)THEN
      KK=K
    ELSE
      !humus complex
      KK=micpar%k_humus
    ENDIF
    IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
      OSCI(K) = CORGCX(K)*VLSoilMicPMass_vr(L,NY,NX)
      OSNI(K) = CORGNX(K)*VLSoilMicPMass_vr(L,NY,NX)
      OSPI(K) = CORGPX(K)*VLSoilMicPMass_vr(L,NY,NX)
    ELSE
      OSCI(K) = CORGCX(K)*VGeomLayer_vr(L,NY,NX)
      OSNI(K) = CORGNX(K)*VGeomLayer_vr(L,NY,NX)
      OSPI(K) = CORGPX(K)*VGeomLayer_vr(L,NY,NX)
    ENDIF

    TOSCK(K) = OMCK(K)+ORCK(K)+OQCK(K)+OHCK(K)
    TOSNK(K) = ORCK(K)*CNRH(K)+OQCK(K)*CNOSCT(KK)+OHCK(K)*CNOSCT(KK)
    TOSPK(K) = ORCK(K)*CPRH(K)+OQCK(K)*CPOSCT(KK)+OHCK(K)*CPOSCT(KK)
!   based on aerobic heterotrophs

    TOSNK(K)=TOSNK(K)+OMCI(iLbiom_kinetic,K)*rNCOMC_ave(iLbiom_kinetic,mid_HeterAerobBacter,K)+OMCI(iLbiom_struct,K)*rNCOMC_ave(iLbiom_struct,mid_HeterAerobBacter,K)
    TOSPK(K)=TOSPK(K)+OMCI(iLbiom_kinetic,K)*rPCOMC_ave(iLbiom_kinetic,mid_HeterAerobBacter,K)+OMCI(iLbiom_struct,K)*rPCOMC_ave(iLbiom_struct,mid_HeterAerobBacter,K)

    TOSCI   = TOSCI+OSCI(K)*TOSCK(K)
    TOSNI   = TOSNI+OSCI(K)*TOSNK(K)
    TOSPI   = TOSPI+OSCI(K)*TOSPK(K)
    OSCX(K) = 0.0_r8
    OSNX(K) = 0.0_r8
    OSPX(K) = 0.0_r8
  ENDDO D995

  mBiomeAutor_vr(1:NumPlantChemElms,1:NumLiveAutoBioms,L,NY,NX)=0._r8

  D8995: DO K=1,jcplx
    IF(L.EQ.0)THEN
      OSCM(K) = AMIN1(DCKR,1._r8)*CORGCX(K)*VLSoilMicPMass_vr(L,NY,NX)
      X       = 0.0_r8
      KK      = K
      FOSCI   = 1.0_r8
      FOSNI   = 1.0_r8
      FOSPI   = 1.0_r8
    ELSE
      IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
        IF(K.LE.micpar%NumOfLitrCmplxs)THEN
          OSCM(K)=AMIN1(DCKR,1._r8)*CORGCX(K)*VLSoilMicPMass_vr(L,NY,NX)
        ELSE
          OSCM(K)=AMIN1(FCX,1._r8)*CORGCX(K)*VLSoilMicPMass_vr(L,NY,NX)*DCKM/(CORGCX(k_humus)+DCKM)
        ENDIF
      ELSE
        IF(K.LE.micpar%NumOfLitrCmplxs)THEN
          OSCM(K)=AMIN1(DCKR,1._r8)*CORGCX(K)*VGeomLayer_vr(L,NY,NX)
        ELSE
          OSCM(K)=AMIN1(FCX,1._r8)*CORGCX(K)*VGeomLayer_vr(L,NY,NX)*DCKM/(CORGCX(k_humus)+DCKM)
        ENDIF
      ENDIF

      X  = 1.0_r8
      KK = micpar%k_humus
      IF(TOSCI.GT.ZEROS(NY,NX))THEN
        FOSCI = AMIN1(1.0_r8,OSCI(KK)/TOSCI)
        FOSNI = AMIN1(1.0_r8,OSCI(KK)*CNOSCT(KK)/TOSNI)
        FOSPI = AMIN1(1.0_r8,OSCI(KK)*CPOSCT(KK)/TOSPI)
      ELSE
        FOSCI = 0.0_r8
        FOSNI = 0.0_r8
        FOSPI = 0.0_r8
      ENDIF
    ENDIF
    !
    !     MICROBIAL C, N AND P
    !
    !     OMC,OMN,OMP=microbial C,N,P
    !     OMCI=microbial biomass content in litter
    !     OMCF,OMCA=hetero,autotrophic biomass composition in litter
    !     rNCOMC,rPCOMC=maximum N:C and P:C ratios in microbial biomass
    !     OSCX,OSNX,OSPX=remaining unallocated SOC,SON,SOP
    !     The reason that initialization of complex-5 microbes is repated for each
    !     complex is because complex 5 is shared by all the other complexes

    D8990: DO N=1,NumMicbHFunGrupsPerCmplx
      
      D8991: DO M=1,nlbiomcp
        OME1(ielmc) = AZMAX1(OSCM(K)*OMCI(M,K)*OMCF(N)*FOSCI)
        OME1(ielmn) = AZMAX1(OME1(ielmc)*rNCOMC_ave(M,N,K)*FOSNI)
        OME1(ielmp) = AZMAX1(OME1(ielmc)*rPCOMC_ave(M,N,K)*FOSPI)

        tglds=JGnfH(N)-JGniH(N)+1._r8
        DO NGL=JGniH(N),JGnfH(N)
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeHeter_vr(NE,MID,K,L,NY,NX)=OME1(NE)/tglds
          ENDDO                
        ENDDO

        OSCX (KK)    = OSCX(KK)+OME1(ielmc)
        OSNX (KK)    = OSNX(KK)+OME1(ielmn)
        OSPX (KK)    = OSPX(KK)+OME1(ielmp)

        D8992: DO NN = 1, NumMicbAFunGrupsPerCmplx
          tglds=JGnfA(NN)-JGniA(NN)+1._r8
          do NGL=JGniA(NN),JGnfA(NN)
            MID=micpar%get_micb_id(M,NGL)
            DO NE=1,NumPlantChemElms
              mBiomeAutor_vr(NE,MID,L,NY,NX)=mBiomeAutor_vr(NE,MID,L,NY,NX)+OME1(NE)*OMCA(NN)/tglds
            ENDDO
          ENDDO
          OSCX(KK) = OSCX(KK)+OME1(ielmc)*OMCA(NN)
          OSNX(KK) = OSNX(KK)+OME1(ielmn)*OMCA(NN)
          OSPX(KK) = OSPX(KK)+OME1(ielmp)*OMCA(NN)

        ENDDO D8992
      ENDDO D8991
    ENDDO D8990

    !for cyanobacteria
    tCyanoInocC=CyanoInocC*AREA_3D(3,NU_col(NY,NX),NY,NX)    
    if(L.eq.0 .and. OSCM(K).GT.0._r8)then
      KL=micpar%NumOfLitrCmplxs

      call InoculateCyanoBacter(K,L,NY,NX,KL,tCyanoInocC)

    elseif(L.eq.NU_col(NY,NX) .and. OSCM(K).GT.0._r8)then
      KL=jcplx
      call InoculateCyanoBacter(K,L,NY,NX,KL,tCyanoInocC)
    endif
    !
    !     MICROBIAL RESIDUE C, N AND P
    !
    !     ORC,ORN,ORP=microbial residue C,N,P
    !     ORCI=allocation of microbial residue to kinetic components
    !  X is an indicator of surface residual layer, 1 for subsurface layer
    D8985: DO M=1,ndbiomcp
      OMBioResdu_vr(ielmc,M,K,L,NY,NX) = X*AZMAX1(OSCM(K)*ORCI(M,K)*FOSCI)
      OMBioResdu_vr(ielmn,M,K,L,NY,NX) = AZMAX1(OMBioResdu_vr(ielmc,M,K,L,NY,NX)*rNCOMC_ave(M,mid_HeterAerobBacter,K)*FOSNI)
      OMBioResdu_vr(ielmp,M,K,L,NY,NX) = AZMAX1(OMBioResdu_vr(ielmc,M,K,L,NY,NX)*rPCOMC_ave(M,mid_HeterAerobBacter,K)*FOSPI)
      OSCX(KK)                         = OSCX(KK)+OMBioResdu_vr(ielmc,M,K,L,NY,NX)
      OSNX(KK)                         = OSNX(KK)+OMBioResdu_vr(ielmn,M,K,L,NY,NX)
      OSPX(KK)                         = OSPX(KK)+OMBioResdu_vr(ielmp,M,K,L,NY,NX)
    ENDDO D8985
    !
    !     DOC, DON AND DOP
    !
    !     OQC,OQN,OQP,OQA=DOC,DON,DOP,acetate in micropores (g)
    !     OQCH,OQNH,OQPH,OQAH=DOC,DON,DOP,acetate in macropores (g)
    !
    DOM_MicP_vr(idom_doc,K,L,NY,NX)               = X*AZMAX1(OSCM(K)*OQCK(K)*FOSCI)
    DOM_MicP_vr(idom_don,K,L,NY,NX)               = AZMAX1(DOM_MicP_vr(idom_doc,K,L,NY,NX)*CNOSCT(KK)*FOSNI)
    DOM_MicP_vr(idom_dop,K,L,NY,NX)               = AZMAX1(DOM_MicP_vr(idom_doc,K,L,NY,NX)*CPOSCT(KK)*FOSPI)
    DOM_MicP_vr(idom_acetate,K,L,NY,NX)           = 0.0_r8
    DOM_MacP_vr(idom_beg:idom_end,K,L,NY,NX) = 0.0_r8
    OSCX(KK)                                 = OSCX(KK)+DOM_MicP_vr(idom_doc,K,L,NY,NX)
    OSNX(KK)                                 = OSNX(KK)+DOM_MicP_vr(idom_don,K,L,NY,NX)
    OSPX(KK)                                 = OSPX(KK)+DOM_MicP_vr(idom_dop,K,L,NY,NX)
    !
    !     ADSORBED C, N AND P
    !
    !     OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
    !
    SorbedOM_vr(ielmc,K,L,NY,NX)        = X*AZMAX1(OSCM(K)*OHCK(K)*FOSCI)
    SorbedOM_vr(ielmn,K,L,NY,NX)        = AZMAX1(SorbedOM_vr(ielmc,K,L,NY,NX)*CNOSCT(KK)*FOSNI)
    SorbedOM_vr(ielmp,K,L,NY,NX)        = AZMAX1(SorbedOM_vr(ielmc,K,L,NY,NX)*CPOSCT(KK)*FOSPI)
    SorbedOM_vr(idom_acetate,K,L,NY,NX) = 0.0_r8
    OSCX(KK)                            = OSCX(KK)+SorbedOM_vr(ielmc,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)
    OSNX(KK)                            = OSNX(KK)+SorbedOM_vr(ielmn,K,L,NY,NX)
    OSPX(KK)                            = OSPX(KK)+SorbedOM_vr(ielmp,K,L,NY,NX)
    !
    !     HUMUS C, N AND P
    !
    !     OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP

    D8980: DO M=1,jsken

      SolidOM_vr(ielmc,M,K,L,NY,NX)=AZMAX1(CFOSC_vr(M,K,L,NY,NX)*(OSCI(K)-OSCX(K)))
      IF(CNOSCT(K).GT.ZERO)THEN
        SolidOM_vr(ielmn,M,K,L,NY,NX)=AZMAX1(CFOSC_vr(M,K,L,NY,NX)*CNOSC_vr(M,K,L,NY,NX) &
          /CNOSCT(K)*(OSNI(K)-OSNX(K)))
      ELSE
        SolidOM_vr(ielmn,M,K,L,NY,NX)=0.0_r8
      ENDIF
      IF(CPOSCT(K).GT.ZERO)THEN
        SolidOM_vr(ielmp,M,K,L,NY,NX)=AZMAX1(CFOSC_vr(M,K,L,NY,NX)*CPOSC_vr(M,K,L,NY,NX) &
          /CPOSCT(K)*(OSPI(K)-OSPX(K)))
      ELSE
        SolidOM_vr(ielmp,M,K,L,NY,NX)=0.0_r8
      ENDIF

      IF(K.EQ.micpar%k_woody_comp)THEN
        SolidOMAct_vr(M,K,L,NY,NX)=SolidOM_vr(ielmc,M,K,L,NY,NX)*OMCI(1,K)
      ELSE
        SolidOMAct_vr(M,K,L,NY,NX)=SolidOM_vr(ielmc,M,K,L,NY,NX)
      ENDIF
    ENDDO D8980
  ENDDO D8995

  !
  !     ADD ALL LITTER,POC,HUMUS COMPONENTS TO GET TOTAL SOC
  !

  RO2DmndHetert_vr(:,:,L,NY,NX)             = 0.0_r8
  RNO3ReduxDmndSoilHeter_vr(:,:,L,NY,NX) = 0.0_r8
  RNO2DmndReduxSoilHeter_vr(:,:,L,NY,NX) = 0.0_r8
  RN2ODmndReduxHeter_vr(:,:,L,NY,NX)     = 0.0_r8
  RNH4DmndSoilHeter_vr(:,:,L,NY,NX)      = 0.0_r8
  RNH4DmndBandHeter_vr(:,:,L,NY,NX)      = 0.0_r8
  RNO3DmndSoilHeter_vr(:,:,L,NY,NX)      = 0.0_r8
  RH2PO4DmndSoilHeter_vr(:,:,L,NY,NX)    = 0.0_r8
  RNO3DmndBandHeter_vr(:,:,L,NY,NX)      = 0.0_r8
  RH2PO4DmndBandHeter_vr(:,:,L,NY,NX)    = 0.0_r8
  RH1PO4DmndSoilHeter_vr(:,:,L,NY,NX)    = 0.0_r8
  RH1PO4DmndBandHeter_vr(:,:,L,NY,NX)    = 0.0_r8
  IF(L.EQ.0)THEN
    RNH4DmndLitrHeter_col(:,:,NY,NX)   = 0.0_r8
    RNH4DmndLitrBandHeter_col(:,:,NY,NX)   = 0.0_r8
    RNO3DmndLitrHeter_col(:,:,NY,NX)   = 0.0_r8
    RNO3DmndLitrBandHeter_col(:,:,NY,NX)   = 0.0_r8
    RH2PO4DmndLitrHeter_col(:,:,NY,NX) = 0.0_r8
    RH2PO4DmndLitrBandHeter_col(:,:,NY,NX) = 0.0_r8
    RH1PO4DmndLitrHeter_col(:,:,NY,NX) = 0.0_r8
    RH1PO4DmndLitrBandHeter_col(:,:,NY,NX) = 0.0_r8
  ENDIF

  RO2MetaDmndAutor_vr(:,L,NY,NX)       = 0.0_r8
  RNH3OxidAutor_vr(:,L,NY,NX)          = 0.0_r8
  RNO2XupAutor_vr(:,L,NY,NX)          = 0.0_r8
  RNO3XupAutor_vr(:,L,NY,NX)          = 0.0_r8
  RN2ODmndReduxAutor_vr(:,L,NY,NX)  = 0.0_r8
  RNH4UptkSoilAutor_vr(:,L,NY,NX)   = 0.0_r8
  RNH4UptkBandAutor_vr(:,L,NY,NX)   = 0.0_r8
  RNO3UptkSoilAutor_vr(:,L,NY,NX)   = 0.0_r8
  RNO3UptkBandAutor_vr(:,L,NY,NX)   = 0.0_r8
  RH2PO4UptkSoilAutor_vr(:,L,NY,NX) = 0.0_r8
  RH2PO4UptkBandAutor_vr(:,L,NY,NX) = 0.0_r8
  RH1PO4UptkSoilAutor_vr(:,L,NY,NX) = 0.0_r8
  RH1PO4UptkBandAutor_vr(:,L,NY,NX) = 0.0_r8

  IF(L.EQ.0)THEN
    RNH4UptkLitrAutor_col(:,NY,NX)   = 0.0_r8
    RNH4UptkLitrBandAutor_col(:,NY,NX)   = 0.0_r8
    RNO3UptkLitrAutor_col(:,NY,NX)   = 0.0_r8
    RNO3UptkLitrBandAutor_col(:,NY,NX)   = 0.0_r8
    RH2PO4UptkLitrAutor_col(:,NY,NX) = 0.0_r8
    RH2PO4UptkLitrBandAutor_col(:,NY,NX) = 0.0_r8
    RH1PO4UptkLitrAutor_col(:,NY,NX) = 0.0_r8
    RH1PO4UptkLitrBandAutor_col(:,NY,NX) = 0.0_r8
  ENDIF
  
  call sumORGMLayL(L,NY,NX,ORGM)

  SoilOrgM_vr(1:NumPlantChemElms,L,NY,NX)=ORGM(1:NumPlantChemElms)
  
  ORGCX_vr(L,NY,NX)=SoilOrgM_vr(ielmc,L,NY,NX)

  call sumLitrOMLayL(L,NY,NX,litrOM)

  OMLitrC_vr(L,NY,NX)=litrOM(ielmc)

  end associate
  end subroutine InitSOMVars
!------------------------------------------------------------------------------------------
  subroutine InoculateCyanoBacter(K,L,NY,NX,KL,CyanoInocC,add_to_existing)
  implicit none
  integer, intent(in) :: K,L,NY,NX,KL
  real(r8), intent(in) :: CyanoInocC
  logical, optional, intent(in) :: add_to_existing
  integer :: N,M,NGL,MID,NE
  real(r8) :: OME1(1:NumPlantChemElms)
  real(r8) :: tglds,tOMCI
  logical :: additive

  additive=.false.
  if(present(add_to_existing))additive=add_to_existing

  N = micpar%mid_HeterMixtCynoBacter
  tglds = JGnfH(N)-JGniH(N)+1._r8
  tOMCI=sum(micpar%OMCI(:,K))
  DO M=1,micpar%nlbiomcp
    OME1(ielmc) = CyanoInocC*micpar%OMCI(M,K)/(tOMCI*KL)
    OME1(ielmn) = AZMAX1(OME1(ielmc)*micpar%rNCOMC_ave(M,N,K))
    OME1(ielmp) = AZMAX1(OME1(ielmc)*micpar%rPCOMC_ave(M,N,K))

    DO NGL=JGniH(N),JGnfH(N)
      MID=micpar%get_micb_id(M,NGL)
      DO NE=1,NumPlantChemElms
        if(additive)then
          mBiomeHeter_vr(NE,MID,K,L,NY,NX)=mBiomeHeter_vr(NE,MID,K,L,NY,NX)+OME1(NE)/tglds
        else
          mBiomeHeter_vr(NE,MID,K,L,NY,NX)=OME1(NE)/tglds
        endif
      ENDDO
    ENDDO
  ENDDO
  end subroutine InoculateCyanoBacter

!------------------------------------------------------------------------------------------
  subroutine InitSurfResiduKinetiComponent(L,NY,NX)
  implicit none
  integer, intent(in) :: L, NY,NX
!     begin_execution
  associate(                               &
   k_woody_comp  => micpar%k_woody_comp  , &
   k_fine_comp   => micpar%k_fine_comp   , &
   k_manure      => micpar%k_manure        &
  )
  IF(L.EQ.0)THEN
    !
    !     CFOSC=fraction of litter in protein(1),nonstructural(2)
    !     cellulose(3) and lignin(4)
    !
    !     PREVIOUS COARSE WOODY RESIDUE
    !
    CFOSC_vr(1:jsken,k_woody_comp,L,NY,NX)=real((/0.00,0.045,0.660,0.295/),r8)

    !
    !     MAIZE
    !
    IF(iLitrType_col(1,NY,NX).EQ.1)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.080,0.245,0.613,0.062/),r8)
      !
      !     WHEAT
      !
    ELSEIF(iLitrType_col(1,NY,NX).EQ.2)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.125,0.171,0.560,0.144/),r8)
      !
      !     SOYBEAN
      !
    ELSEIF(iLitrType_col(1,NY,NX).EQ.3)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.138,0.426,0.316,0.120/),r8)

      !     NEW STRAW
      !
    ELSEIF(iLitrType_col(1,NY,NX).EQ.4)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.036,0.044,0.767,0.153/),r8)

      !     OLD STRAW
      !
    ELSEIF(iLitrType_col(1,NY,NX).EQ.5)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.075,0.125,0.550,0.250/),r8)

      !
      !     COMPOST
      !
    ELSEIF(iLitrType_col(1,NY,NX).EQ.6)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.143,0.015,0.640,0.202/),r8)
      !
      !     GREEN MANURE
      !
    ELSEIF(iLitrType_col(1,NY,NX).EQ.7)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.202,0.013,0.560,0.225/),r8)

      !     NEW DECIDUOUS FOREST
!
    ELSEIF(iLitrType_col(1,NY,NX).EQ.8)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.070,0.41,0.36,0.16/),r8)

!
      !     NEW CONIFEROUS FOREST
!
    ELSEIF(iLitrType_col(1,NY,NX).EQ.9)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.07,0.25,0.38,0.30/),r8)

      !     OLD DECIDUOUS FOREST
!
    ELSEIF(iLitrType_col(1,NY,NX).EQ.10)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.02,0.06,0.34,0.58/),r8)
      !
      !     OLD CONIFEROUS FOREST
      !
    ELSEIF(iLitrType_col(1,NY,NX).EQ.11)THEN
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.02,0.06,0.34,0.58/),r8)
      !     DEFAULT
      !
    ELSE
      CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)=real((/0.075,0.125,0.550,0.250/),r8)

    ENDIF
    !
    !     PREVIOUS COARSE (K=0) AND FINE (K=1) ROOTS
    !
  ELSE
    CFOSC_vr(1:jsken,k_woody_comp,L,NY,NX) = real((/0.00,0.00,0.20,0.80/),r8)
    CFOSC_vr(1:jsken,k_fine_comp,L,NY,NX)  = real((/0.02,0.06,0.34,0.58/),r8)
  ENDIF
  end associate
  end subroutine InitSurfResiduKinetiComponent
!------------------------------------------------------------------------------------------
  subroutine InitManureKinetiComponent(L,NY,NX)

  implicit none
  integer, intent(in) :: L, NY, NX
  !     begin_execution
  associate(                      &
    k_manure  => micpar%k_manure, &
    jsken     => micpar%jsken     &
  )
  !
  !     RUMINANT
!
  IF(iLitrType_col(2,NY,NX).EQ.1)THEN
    CFOSC_vr(1:jsken,k_manure,L,NY,NX)=real((/0.036,0.044,0.630,0.290/),r8)

!
    !     NON-RUMINANT
!
  ELSEIF(iLitrType_col(2,NY,NX).EQ.2)THEN
    CFOSC_vr(1:jsken,k_manure,L,NY,NX)=real((/0.138,0.401,0.316,0.145/),r8)
!
!     OTHER
!
  ELSE
    CFOSC_vr(1:jsken,k_manure,L,NY,NX)=real((/0.138,0.401,0.316,0.145/),r8)
  ENDIF
  end associate
  end subroutine InitManureKinetiComponent
!------------------------------------------------------------------------------------------
  subroutine InitPOMKinetiComponent(L,NY,NX,HCX,TORGL,LandScape1stSoiLayDepth,FCX)
  implicit none
  integer, intent(in)  :: L, NY, NX
  real(r8), intent(in) :: HCX
  real(r8), intent(in) :: TORGL
  real(r8), intent(in) :: LandScape1stSoiLayDepth
  real(r8), intent(out):: FCX
  real(r8) :: FCY,FC0,FC1

! begin_execution
  associate(                           &
    k_POM       => micpar%k_POM,       &
    k_humus     => micpar%k_humus,     &
    k_fine_comp => micpar%k_fine_comp, &
    k_manure    => micpar%k_manure,    &
    iprotein    => micpar%iprotein,    &
    icarbhyro   => micpar%icarbhyro,   &
    icellulos   => micpar%icellulos,   &
    ilignin     => micpar%ilignin,     &
    jsken       => micpar%jsken        &
  )
!
! CFOSC=siNGLe kinetic fraction in POM
!
  IF(L.NE.0)THEN
    CFOSC_vr(1:jsken,k_POM,L,NY,NX)=real((/1.0,0.0,0.0,0.0/),r8)

!
!  HUMUS PARTITIONED TO DIFFERENT FRACTIONS
!  BASED ON SOC ACCUMULATION ABOVE EACH LAYER
!
!  NATURAL SOILS
!
    IF(ISOILR_col(NY,NX).EQ.0)THEN
!
      !     DRYLAND SOIL
      !
      !     CORGC,FORGC=SOC,minimum SOC for organic soil(g Mg-1)
      !     DPTH,ExtWaterTableLocDepz_col=depth to layer midpoint,external water table(m)
      !     FC0=partitioning to less resistant component at DPTH=0
      !     FCX=reduction in FC0 at DPTH
      !     CORGCX,CORGNX,CORGPX=C,N,P concentations in humus
!
      IF(CSoilOrgM_vr(ielmc,L,NY,NX).LE.FORGC .OR. SoilDepthMidLay_vr(L,NY,NX).LE.ExtWaterTableLocDepz_col(NY,NX) &
        +CumDepz2LayBottom_vr(NU_col(NY,NX),NY,NX)-LandScape1stSoiLayDepth)THEN
        FCY=0.60_r8
        IF(CORGCX(k_humus).GT.1.0E-32_r8)THEN
          FC0=FCY*EXP(-5.0_r8*(AMIN1(CORGNX(k_humus), &
            10.0_r8*CORGPX(k_humus))/CORGCX(k_humus)))
        ELSE
          FC0=FCY
        ENDIF
        FCX=EXP(HCX*TORGL)
        !     WETLAND
        !
      ELSE
        FCY=0.60_r8
        IF(CORGCX(k_humus).GT.1.0E-32_r8)THEN
          FC0=FCY*EXP(-5.0_r8*(AMIN1(CORGNX(k_humus),10.0_r8*CORGPX(k_humus))/CORGCX(k_humus)))
        ELSE
          FC0=FCY
        ENDIF
!     FCX=(EXP(HCX*TORGL))**0.5_r8
        FCX=EXP(HCX*TORGL)
      ENDIF
!
      !     RECONSTRUCTED SOILS
!
    ELSE
      FCY=0.30_r8
      IF(CORGCX(k_humus).GT.1.0E-32_r8)THEN
        FC0=FCY*EXP(-5.0_r8*(AMIN1(CORGNX(k_humus), 10.0_r8*CORGPX(k_humus))/CORGCX(k_humus)))
      ELSE
        FC0=FCY
      ENDIF
      FCX=1.0_r8
    ENDIF
!
!   PARTITION HUMUS
!
!   CFOSC=fraction of humus in less(1),more(2) resistant component
!
    FC1=FC0*FCX
    CFOSC_vr(iprotein,k_humus,L,NY,NX)  = FC1
    CFOSC_vr(icarbhyro,k_humus,L,NY,NX) = 1.0_r8-FC1
    CFOSC_vr(icellulos,k_humus,L,NY,NX) = 0.00_r8
    CFOSC_vr(ilignin,k_humus,L,NY,NX)   = 0.00_r8
!
!   MICROBIAL DETRITUS ALLOCATED TO HUMUS MAINTAINS
!   HUMUS PARTITIONING TO COMPONENTS
!
!   ElmAllocmatMicrblitr2POM_vr=fraction of microbial litter allocated to humus components
!
    ElmAllocmatMicrblitr2POM_vr(1,L,NY,NX)=3.0_r8*FC1/(2.0_r8*FC1+1.0_r8)
    ElmAllocmatMicrblitr2POM_vr(2,L,NY,NX)=1.0_r8-ElmAllocmatMicrblitr2POM_vr(1,L,NY,NX)
  ENDIF

  IF(L.GT.0)THEN
    IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
      !assume 0.55 of organic matter is carbon
      SolidOMPercent_vr(L,NY,NX)=(CORGCX(k_fine_comp)+CORGCX(k_manure)+CORGCX(k_POM)+CORGCX(k_humus))*1.e-4_r8/OMCMassFrac ![g Mg-1] -> [g /100g-1]
    endif
  endif
  end associate
  end subroutine InitPOMKinetiComponent

!------------------------------------------------------------------------------------------

  subroutine InitSOMProfile(L,NY,NX,HCX,TORGL,LandScape1stSoiLayDepth,FCX)

  implicit none
  integer,  intent(in)  :: L,NY,NX
  real(r8), intent(in)  :: HCX
  real(r8), intent(in)  :: TORGL
  real(r8), intent(in)  :: LandScape1stSoiLayDepth
  real(r8), intent(out) :: FCX

  call PrintInfo('beg InitSOMProfile')

  call InitLitterProfile(L,NY,NX)

  !     SURFACE RESIDUE KINETIC COMPONENTS
  call InitSurfResiduKinetiComponent(L,NY,NX)
  !
  !     ANIMAL MANURE
  call InitManureKinetiComponent(L,NY,NX)
  !
  !     POM
  call InitPOMKinetiComponent(L,NY,NX,HCX,TORGL,LandScape1stSoiLayDepth,FCX)
  
  call PrintInfo('end InitSOMProfile')
  end subroutine InitSOMProfile

!------------------------------------------------------------------------------------------

  subroutine InitLitterProfile(L,NY,NX)

  implicit none
  integer, intent(in) :: L,NY,NX

  real(r8) :: CORGCZ,CORGNZ,CORGPZ,CORGRZ
  real(r8) :: scal
! begin_execution
  associate(                                   &
    NumOfLitrCmplxs => micpar%NumOfLitrCmplxs, &
    k_POM           => micpar%k_POM,           &
    k_humus         => micpar%k_humus,         &
    CNRH            => micpar%CNRH,            &
    CPRH            => micpar%CPRH             &
  )

  call PrintInfo('beg InitLitterProfile')

  IF(VLSoilMicPMass_vr(L,NY,NX).GT.ZEROS(NY,NX))THEN
    scal                      = AREA_3D(3,L,NY,NX)/VLSoilMicPMass_vr(L,NY,NX) ![m2 d-2] /[Mg soil d-2]=[m2/Mg soil]
    CORGCX(1:NumOfLitrCmplxs) = RSC_vr(1:NumOfLitrCmplxs,L,NY,NX)*scal ![gC m-2]*[m2/Mg soil]=[gC/Mg soil]
    CORGNX(1:NumOfLitrCmplxs) = RSN_vr(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGPX(1:NumOfLitrCmplxs) = RSP_vr(1:NumOfLitrCmplxs,L,NY,NX)*scal
  ELSE
    scal                      = safe_adb(AREA_3D(3,L,NY,NX),VGeomLayer_vr(L,NY,NX))
    CORGCX(1:NumOfLitrCmplxs) = RSC_vr(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGNX(1:NumOfLitrCmplxs) = RSN_vr(1:NumOfLitrCmplxs,L,NY,NX)*scal
    CORGPX(1:NumOfLitrCmplxs) = RSP_vr(1:NumOfLitrCmplxs,L,NY,NX)*scal
  ENDIF
    !
    !     ALLOCATE SOC TO POC(3) AND HUMUS(4)
    !
  IF(L.GT.0)THEN
    CORGCZ = CSoilOrgM_vr(ielmc,L,NY,NX) !total SOC, [gC/Mg soil]
    CORGRZ = COMLitrC_vr(L,NY,NX)        !POC, [gC/Mg soil]
    CORGNZ = CSoilOrgM_vr(ielmn,L,NY,NX)
    CORGPZ = CSoilOrgM_vr(ielmp,L,NY,NX)
    IF(CORGCZ.GT.ZERO)THEN
      CORGCX(k_POM)   = CORGRZ
      CORGCX(k_humus) = AZMAX1(CORGCZ-CORGCX(k_POM))
      CORGNX(k_POM)   = AMIN1(CNRH(k_POM)*CORGCX(k_POM),CORGNZ)
      CORGNX(k_humus) = AZMAX1(CORGNZ-CORGNX(k_POM))
      CORGPX(k_POM)   = AMIN1(CPRH(k_POM)*CORGCX(k_POM),CORGPZ)
      CORGPX(k_humus) = AZMAX1(CORGPZ-CORGPX(k_POM))      
    ELSE
      CORGCX(k_POM)   = 0.0_r8
      CORGCX(k_humus) = 0.0_r8
      CORGNX(k_POM)   = 0.0_r8
      CORGNX(k_humus) = 0.0_r8
      CORGPX(k_POM)   = 0.0_r8
      CORGPX(k_humus) = 0.0_r8
    ENDIF
  ELSE
    CORGCX(k_POM)   = 0.0_r8
    CORGCX(k_humus) = 0.0_r8
    CORGNX(k_POM)   = 0.0_r8
    CORGNX(k_humus) = 0.0_r8
    CORGPX(k_POM)   = 0.0_r8
    CORGPX(k_humus) = 0.0_r8
  ENDIF
  
  call PrintInfo('end InitLitterProfile')
  end associate
  end subroutine InitLitterProfile
!------------------------------------------------------------------------------------------

  subroutine MicrobeByLitterFall(I,J,K,NY,NX,OSCMK,OME_in,mscal)
  !
  !seeding microbes added through litterfall
  implicit none
  integer, intent(in) :: I,J,K
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: OSCMK                       !input organic carbon
  real(r8), intent(out) :: OME_in(NumPlantChemElms)  !summary of total microbial input
  real(r8),optional,intent(in) :: mscal
  integer :: M,N,NGL,MID,NE,NN
  real(r8) :: FOSCI,FOSNI,FOSPI,tglds,tOMCI,fOMCH,tOMCA
  real(r8) :: OME1(1:NumPlantChemElms),OSCMKK
  
  real(r8) :: scal    !scalar for incoming microbial biomass associated with litterfall, [1%].
  associate(                         &
    rNCOMC_ave => micpar%rNCOMC_ave, &
    rPCOMC_ave => micpar%rPCOMC_ave, &
    nlbiomcp   => micpar%nlbiomcp,   &
    OMCI       => micpar%OMCI,       &
    OMCF       => micpar%OMCF,       &
    OMCA       => micpar%OMCA        &
  )

  if(present(mscal))then
    scal=mscal
  else
    scal=0.01_r8 !microbial biomass fraction of the incoming organic matter, [1%] 
  endif
  FOSCI=1._r8; FOSNI=1._r8; FOSPI=1._r8
  OME_in=0._r8; tOMCI=sum(OMCI(:,K))
  tOMCA=sum(OMCA(:));fOMCH=1._r8/(1._r8+tOMCA)
  OSCMKK=OSCMK*scal*fOMCH/tOMCI

  DO N=1,NumMicbHFunGrupsPerCmplx    
    DO M=1,nlbiomcp
      OME1(ielmc) = AZMAX1(OSCMKK*OMCI(M,K)*OMCF(N)*FOSCI)
      OME1(ielmn) = AZMAX1(OME1(ielmc)*rNCOMC_ave(M,N,K)*FOSNI)
      OME1(ielmp) = AZMAX1(OME1(ielmc)*rPCOMC_ave(M,N,K)*FOSPI)
      OME_in=OME_in+OME1
      tglds=JGnfH(N)-JGniH(N)+1._r8
      do NGL=JGniH(N),JGnfH(N)
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          mBiomeHeter_vr(NE,MID,K,0,NY,NX)=mBiomeHeter_vr(NE,MID,K,0,NY,NX)+OME1(NE)/tglds
        ENDDO
      ENDDO

      DO NN = 1, NumMicbAFunGrupsPerCmplx
        tglds=JGnfA(NN)-JGniA(NN)+1._r8
        do NGL=JGniA(NN),JGnfA(NN)
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeAutor_vr(NE,MID,0,NY,NX)=mBiomeAutor_vr(NE,MID,0,NY,NX)+OME1(NE)*OMCA(NN)/tglds
          ENDDO
        ENDDO
        OME_in=OME_in+OME1(:)*OMCA(NN)
      ENDDO
    ENDDO
  ENDDO
  
  end associate
  end subroutine MicrobeByLitterFall

!------------------------------------------------------------------------------------------

  subroutine DestructSOMBGC
  use abortutils, only : destroy
  implicit none

  call destroy(CORGCX)
  call destroy(CORGNX)
  call destroy(CORGPX)

  end subroutine DestructSOMBGC
end module InitSOMBGCMOD
