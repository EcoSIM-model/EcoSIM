module IngridTranspMod
  ! Route dissolved salts through surface, snow and soil water pathways.
  use SoilFaceMappingMod, only: SoilFaceSlot,SoilFaceIncomingRange,GetSoilFaceGeometry
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : AZMAX1,AZMIN1,isclose
  use DebugToolMod
  use TranspSaltDataMod
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use GridDataType
  use SurfLitterDataType
  use SoilBGCDataType
  use SurfSoilDataType
  use SnowDataType
  use ChemTranspDataType
  use EcoSIMSolverPar
  use IrrigationDataType
  implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  real(r8), PARAMETER :: XFRS=0.05_r8
  public :: GetSaltTranspFlxM
  public :: SurfaceSaltFluxM
  public :: SaltFall2Snowpack
  contains

!------------------------------------------------------------------------------------------

  subroutine GetSaltTranspFlxM(I,M,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,M,NHW,NHE,NVN,NVS
  
!     begin_execution

  call SurfaceSaltFluxM(M,NHW,NHE,NVN,NVS)
  !
  CALL SaltSurfRofM(M,NHW,NHE,NVN,NVS)

  call SaltFlowThruSnowDriftM(M,NHW,NHE,NVN,NVS)

  call SaltTranspXInterGridsM(M,NHW,NHE,NVN,NVS)

  call SaltXGridBoundTransptM(M,NHW,NHE,NVN,NVS)

  end subroutine GetSaltTranspFlxM

!------------------------------------------------------------------------------------------
  subroutine SurfaceSaltFluxM(M,NHW,NHE,NVN,NVS,band_transfer_only)

  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  ! True only for the once-per-substep litter/snow-to-band transfer pass.
  logical, optional, intent(in) :: band_transfer_only

  integer :: NY,NX,LS
  logical :: band_pass
  real(r8) :: snow_soil_water,snow_mac_frac
  real(r8) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_Adv2MacP_flx(idsalt_beg:idsaltb_end)
  
  band_pass=.false.
  if(present(band_transfer_only))band_pass=band_transfer_only
  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      if(band_pass)then
        call SnowSaltTranptM(M,NY,NX,LS)
      else
        call SnowSaltTranptM(M,NY,NX)
      endif

      call SoilLitrSaltAdvExchange(M,NY,NX,trcSalt_Adv2MicP_flx,trcSalt_Adv2MacP_flx)

      call Litr2SoilDifusTransptM(M,NY,NX,trcSalt_flx_diffusM)

      ! Split accepted soil-directed snow solutes with the liquid-water flow.
      ! Zero flow retains the existing micropore numerical fallback.
      snow_soil_water=WatFlowSno2MicPM_col(M,NY,NX)+WatFlowSno2MacPM_col(M,NY,NX)
      snow_mac_frac=0._r8
      if(snow_soil_water>0._r8)snow_mac_frac=WatFlowSno2MacPM_col(M,NY,NX)/snow_soil_water
      if(band_pass)then
        call AcceptSurfaceSaltBandTransfersM(NY,NX,LS,trcSalt_Adv2MicP_flx+trcSalt_flx_diffusM,trcSalt_Adv2MacP_flx,snow_mac_frac)
        cycle
      endif
      ! These ID-changing transfers have already been accepted locally.
      trcSalt_Adv2MicP_flx(idsalt_pband_beg:idsalt_pband_end)=0._r8
      trcSalt_flx_diffusM(idsalt_pband_beg:idsalt_pband_end)=0._r8
      trcSalt_Adv2MacP_flx(idsalt_pband_beg:idsalt_pband_end)=0._r8

      call LitrSoilAdvDifusTransptM(NY,NX,trcSalt_Adv2MicP_flx,trcSalt_flx_diffusM,trcSalt_Adv2MacP_flx,snow_mac_frac)

      !
      !     MACROPORE TO MICROPORE TRANSFER
      !
      IF(FWatExMacP2MicPM_vr(M,NU_col(NY,NX),NY,NX).GT.0.0)THEN
        call MacPoreSaltAdvExM(M,NY,NX)
        !
        !     MICROPORE TO MACROPORE TRANSFER
        !
      ELSEIF(FWatExMacP2MicPM_vr(M,NU_col(NY,NX),NY,NX).LT.0.0)THEN
        call MicPoreSaltAdvExM(M,NY,NX)
        !
        !     NO MACROPORE TO MICROPORE TRANSFER
        !
      ELSE
        trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)=0.0_r8
      ENDIF
      !      
      call MacMicPoreSaltDifusExM(M,NY,NX)
      !
      call SaltFlowThruLitrRofM(M,NX,NY)
    ENDDO
  ENDDO
  end subroutine SurfaceSaltFluxM
!------------------------------------------------------------------------------------------

  subroutine AcceptSurfaceSaltBandTransfersM(NY,NX,LS,litter_to_micp,litter_to_macp,snow_mac_frac)
  ! Accept ID-changing phosphate transfers without coupling ordinary soil pools.
  implicit none
  integer, intent(in) :: NY,NX,LS ! LS=0 means no bottom snow donor
  real(r8), intent(in) :: snow_mac_frac ! fraction of soil-directed snow liquid entering macropores
  real(r8), intent(in) :: litter_to_micp(idsalt_beg:idsaltb_end) ! signed litter-to-band proposal [mol d-2]
  real(r8), intent(in) :: litter_to_macp(idsalt_beg:idsaltb_end) ! nonnegative litter-to-band proposal [mol d-2]
  character(len=*), parameter :: subname='AcceptSurfaceSaltBandTransfersM'
  integer :: idsalt,idsband,L
  real(r8) :: litter_mass,soil_mass,snow_mass,outflow,scale
  real(r8) :: micp_flux,macp_flux,snow_flux

  call PrintInfo('beg '//subname)
  L=NU_col(NY,NX)
  DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
    idsband=idsalt-idsalt_psoil_beg+idsalt_pband_beg
    litter_mass=AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX))
    soil_mass=AZMAX1(trcSalt_solml2_vr(idsband,L,NY,NX))
    ! Share the litter donor snapshot between its two outgoing destinations.
    ! Incoming soil/snow transfers do not fund outgoing transfers in this pass.
    outflow=AZMAX1(litter_to_micp(idsband))+litter_to_macp(idsband)
    scale=1._r8
    if(outflow>0._r8)scale=MIN(1._r8,litter_mass/outflow)
    if(litter_to_micp(idsband)>=0._r8)then
      micp_flux=MIN(litter_mass,litter_to_micp(idsband)*scale)
    else
      micp_flux=-MIN(soil_mass,-litter_to_micp(idsband))
    endif
    macp_flux=MIN(litter_to_macp(idsband)*scale,AZMAX1(litter_mass-AZMAX1(micp_flux)))

    snow_flux=0._r8
    if(LS>0)then
      snow_mass=AZMAX1(trcSalt_ml2_snvr(idsalt,LS,NY,NX))
      snow_flux=MIN(snow_mass,AZMAX1(trcSalt_AquaADV_Snow2Soil_flxM(idsband,NY,NX)))
      trcSalt_ml2_snvr(idsalt,LS,NY,NX)=snow_mass-snow_flux
      trcSalt_ml_snvr(idsalt,LS,NY,NX)=trcSalt_ml2_snvr(idsalt,LS,NY,NX)
    endif

    trcSalt_solml2_vr(idsalt,0,NY,NX)=litter_mass-micp_flux-macp_flux
    trcSalt_solml2_vr(idsband,L,NY,NX)=soil_mass+micp_flux+snow_flux*(1._r8-snow_mac_frac)
    trcSalt_soHml2_vr(idsband,L,NY,NX)=trcSalt_soHml2_vr(idsband,L,NY,NX)+macp_flux+snow_flux*snow_mac_frac
    trcSalt_solml_vr(idsalt,0,NY,NX)=trcSalt_solml2_vr(idsalt,0,NY,NX)
    trcSalt_solml_vr(idsband,L,NY,NX)=trcSalt_solml2_vr(idsband,L,NY,NX)
    trcSalt_soHml_vr(idsband,L,NY,NX)=trcSalt_soHml2_vr(idsband,L,NY,NX)
    ! Accumulate this accepted snow transfer once, outside ordinary retries.
    trcSalt_AquaADV_Snow2Soil_flx(idsband,NY,NX)=trcSalt_AquaADV_Snow2Soil_flx(idsband,NY,NX)+snow_flux
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine AcceptSurfaceSaltBandTransfersM
!------------------------------------------------------------------------------------------

  subroutine SoilLitrSaltAdvExchange(M,NY,NX,trcSalt_Adv2MicP_flx,trcSalt_Adv2MacP_flx)
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(out) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)  
  real(r8),intent(out) :: trcSalt_Adv2MacP_flx(idsalt_beg:idsaltb_end)
  integer :: ids,idsalt
  real(r8) :: FLWRM1,VFLW,FracMic,FracMac

  FLWRM1=WatFLoLitr2SoilM_col(M,NY,NX)
  trcSalt_Adv2MicP_flx=0._r8; trcSalt_Adv2MacP_flx=0._r8
  call LitterInfiltrationFractions(FLWRM1,WatFLoLitr2MacPM_col(M,NY,NX),SurfRunoffPotentM_col(M,NY,NX),&
    VLWatMicPM_vr(M,0,NY,NX),ZEROS2(NY,NX),VFLWX,FracMic,FracMac)
!
  IF(FLWRM1.GT.0.0_r8)THEN
    !litter 2 soil
    VFLW=FracMic
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX))
    ENDDO

    DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
      ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B  
      trcSalt_Adv2MicP_flx(idsalt) = VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX))*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
      trcSalt_Adv2MicP_flx(ids)    = VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX))*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)
    ENDDO
   !soil advection to litter
  ELSEIF(FLWRM1.LT.0._r8)THEN
    IF(VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX))
    ENDDO  
  ENDIF

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_Adv2MacP_flx(idsalt)=FracMac*AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX))
  ENDDO
  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
    trcSalt_Adv2MacP_flx(idsalt)=FracMac*AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX))*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
    trcSalt_Adv2MacP_flx(ids)=FracMac*AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX))*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)
  ENDDO
  end subroutine SoilLitrSaltAdvExchange
!------------------------------------------------------------------------------------------

  subroutine SaltXGridBoundTransptM(M,NHW,NHE,NVN,NVS)
  !
  ! Description:
  ! X-boundary tracer transport
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  integer :: NY,NX,L,N1,N2,N3,N4,N5,N6,N4B,N5B,N,NN
  integer :: M1,M2,M3,M4,M5,M6,idsalt
  real(r8) :: RCHQF,RCHGFU,RCHGFT
  real(r8) :: FQRM

!     begin_execution

  D9595: DO NX=NHW,NHE
    D9590: DO NY=NVN,NVS
      trcSalt_FloXSurRof_flxM_col(idsalt_beg:idsalt_end,NY,NX) = 0._r8
      trcSalt_SnowDrift_flxM_col(idsalt_beg:idsalt_end,NY,NX)  = 0.0
      N1=NX;N2=NY

      D9585: DO L=NU_col(NY,NX),NL_col(NY,NX)
        N3=L
        M1 = NX;M2 = NY;M3 = L        
        !locate dest grid
        D9580: DO N=FlowDirIndicator_col(NY,NX),3
          D9575: DO NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              N4  = NX+1;N5  = NY
              N4B = NX-1;N5B = NY;N6  = L
              !configure boundaries
              IF(NN.EQ.iFront)THEN
                IF(NX.EQ.NHE)THEN               !on the eastern boundary
                  M4 = NX+1;M5 = NY;M6 = L
                  RCHQF  = RechargEastSurf_col(M2,M1)
                  !not boundary
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN
                IF(NX.EQ.NHW)THEN             !on the western boundary
                  M4 = NX;M5 = NY;M6 = L
                  RCHQF  = RechargWestSurf_col(M5,M4)
                ELSE
                  CYCLE
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              N4  = NX;N5  = NY+1
              N4B = NX;N5B = NY-1; N6=L
              IF(NN.EQ.iFront)THEN          
                IF(NY.EQ.NVS)THEN             !on the southern boundary
                  M4 = NX;M5 = NY+1;M6 = L
                  RCHQF  = RechargSouthSurf_col(M2,M1)
                ELSE
                  CYCLE
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN
                IF(NY.EQ.NVN)THEN      !on the northern boundary
                  M4=NX;M5=NY;M6=L
                  RCHQF  = RechargNorthSurf_col(M5,M4)
                ELSE
                  CYCLE
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN
              N4 = NX;N5 = NY;N6 = L+1
              IF(NN.EQ.iFront)THEN
                IF(L.EQ.NL_col(NY,NX))THEN        !at the bottom boundary
                  M4 = NX;M5 = NY;M6 = L+1
                ELSE
                  CYCLE
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN
                CYCLE
              ENDIF
            ENDIF
            !
            !surface horizontal flow
            IF(L.EQ.NUM_col(M2,M1) .AND. N.NE.iVerticalDirection)THEN
              !No surface runoff
              IF(.not.XGridRunoffFlag_2DH(NN,N,N2,N1) .OR. isclose(RCHQF,0.0_r8) &
                .OR. SurfRunoffPotentM_col(M,N2,N1).LE.ZEROS(N2,N1))THEN
                trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,N,NN,M5,M4)=0.0_r8
                !
              ELSE
                IF((NN.EQ.iFront .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
                  .OR.(NN.EQ.iBehind .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
                  !
                  FQRM=QflxSurfRunoffM_2DH(M,N,NN,M5,M4)/SurfRunoffPotentM_col(M,N2,N1)
                  DO idsalt=idsalt_beg,idsalt_end
                    trcSalt_FloXSurRof_flxM_2DH(idsalt,N,NN,M5,M4)=trcSalt_FloXSurRof_flxM(idsalt,N2,N1)*FQRM
                  ENDDO    
                  !
                ELSEIF((NN.EQ.iBehind .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
                  .OR. (NN.EQ.iFront .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN
                  trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,N,NN,M5,M4)=0.0_r8
                ELSE
                  trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,N,NN,M5,M4)=0.0_r8
                ENDIF
              ENDIF              
            ENDIF

            ! SOLUTE LOSS WITH SUBSURFACE MICROPORE WATER LOSS
            !
            IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS(NY,NX))THEN  !legitimate grid
              IF(FlowDirIndicator_col(M2,M1).NE.3 .OR. N.EQ.iVerticalDirection)THEN
                IF(NN.EQ.iFront .AND. WaterFlow2MicPM_3D(M,N,M6,M5,M4).GT.0.0_r8 &          !flow out at outgoing boundary
                  .OR. NN.EQ.iBehind .AND. WaterFlow2MicPM_3D(M,N,M6,M5,M4).LT.0.0_r8)THEN  !flow in at incoming boundary
                  !upstream advection
                  call SoluteMicporeLossXBoundaryM(M,N,M1,M2,M3,M4,M5,M6)

                  !     SOLUTE GAIN WITH SUBSURFACE MICROPORE WATER GAIN
                ELSE
                  !
                  call SoluteMicporeGainXBoundaryM(M,N,M1,M2,M3,M4,M5,M6)
                ENDIF
                !
                !     SOLUTE LOSS WITH SUBSURFACE MACROPORE WATER LOSS
                !
                IF(NN.EQ.iFront .AND. WaterFlow2MacPM_3D(M,N,M6,M5,M4).GT.0.0_r8 &
                  .OR. NN.EQ.iBehind .AND. WaterFlow2MacPM_3D(M,N,M6,M5,M4).LT.0.0_r8)THEN

                  call SoluteMacporeLossXBoundaryM(M,N,M1,M2,M3,M4,M5,M6)

                ELSE
                  !     NO SOLUTE GAIN IN SUBSURFACE MACROPORES
                  trcSalt_MacpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,M6,M5,M4)=0.0_r8
                ENDIF
              ENDIF
            ENDIF
          ENDDO D9575
          !
          !     TOTAL SOLUTE FLUXES IN EACH GRID CELL
          !
          IF(L.EQ.NUM_col(N2,N1))THEN    !topsoil layer
            IF(N.NE.iVerticalDirection)THEN

              !     NET OVERLAND SOLUTE FLUX IN WATER
              call AccumRofSaltFlxM(M,N,N1,N2,N4,N5,N4B,N5B,trcSalt_FloXSurRof_flxM_col)

              !     NET OVERLAND SOLUTE FLUX IN SNOW
              call AccumSnowDriftSaltFlxM(N,N1,N2,N4,N5,N4B,N5B, trcSalt_SnowDrift_flxM_col)
            ELSEIF(N.EQ.iVerticalDirection)THEN
              !     NET SOLUTE FLUX IN SNOWPACK
              call AccumSaltTranspInSnowM(M,NY,NX)
            ENDIF
          ENDIF

          !     TOTAL SOLUTE FLUX IN MICROPORES AND MACROPORES
          call AccumSoilPoreSaltTranspM(N,N1,N2,N3,N4,N5,N6)

        ENDDO D9580
      ENDDO D9585
    ENDDO D9590
  ENDDO D9595
  end subroutine SaltXGridBoundTransptM

!------------------------------------------------------------------------------------------

  subroutine AccumSoilPoreSaltTranspM(N,N1,N2,N3,N4,N5,N6)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,N1,N2,N3,N4,N5
  integer, intent(inout) :: N6
  integer :: LL,idsalt
  integer :: first_face,last_face

  !
  !     begin_execution
  !
  IF(FlowDirIndicator_col(N2,N1).NE.3 .OR. N.EQ.iVerticalDirection)THEN
    ! Match the destination mask used by the face producer. Keep the
    ! supplied index when no active layer exists (external/ghost face).
    D1200: DO LL=N6,NL_col(N5,N4)
      IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
        N6=LL
        exit
      ENDIF
    ENDDO D1200

  call SoilFaceIncomingRange(N,N3,N2,N1,first_face,last_face)
    IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      DO idsalt=idsalt_beg,idsaltb_end
        trcSalt_Transp2Micp_flxM_vr(idsalt,N3,N2,N1)=trcSalt_Transp2Micp_flxM_vr(idsalt,N3,N2,N1) &
          +sum(trcSalt_MicpTranspFlxM_3D(idsalt,N,first_face:last_face,N2,N1))-trcSalt_MicpTranspFlxM_3D(idsalt,N,SoilFaceSlot(N,N3,N6),N5,N4)
        trcSalt_Transp2Macp_flxM_vr(idsalt,N3,N2,N1)=trcSalt_Transp2Macp_flxM_vr(idsalt,N3,N2,N1) &
          +sum(trcSalt_MacpTranspFlxM_3D(idsalt,N,first_face:last_face,N2,N1))-trcSalt_MacpTranspFlxM_3D(idsalt,N,SoilFaceSlot(N,N3,N6),N5,N4)
      ENDDO

    ELSE
      trcSalt_Transp2Micp_flxM_vr(idsalt_beg:idsaltb_end,N3,N2,N1)=0.0
      trcSalt_Transp2Macp_flxM_vr(idsalt_beg:idsaltb_end,N3,N2,N1)=0.0
    ENDIF
  ENDIF
  end subroutine AccumSoilPoreSaltTranspM


!------------------------------------------------------------------------------------------

  subroutine SoluteMacporeLossXBoundaryM(M,N,M1,M2,M3,M4,M5,M6)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N,M1,M2,M3,M4,M5,M6
  real(r8) :: VFLW
  integer :: idsalt,ids
!     begin_execution

  IF(VLWatMacPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
    VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,M6,M5,M4)/VLWatMacPM_vr(M,M3,M2,M1)))

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_MacpTranspFlxM_3D(idsalt,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,M3,M2,M1))
    ENDDO

    DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
      ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
      trcSalt_MacpTranspFlxM_3D(idsalt,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,M3,M2,M1)) &
        *trcs_VLN_vr(ids_H1PO4,M3,M2,M1)

      trcSalt_MacpTranspFlxM_3D(ids,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_soHml2_vr(ids,M3,M2,M1)) &
        *trcs_VLN_vr(ids_H1PO4B,M3,M2,M1)
    ENDDO
  ELSE
    trcSalt_MacpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,M6,M5,M4)=0._r8
  ENDIF

  end subroutine SoluteMacporeLossXBoundaryM


!------------------------------------------------------------------------------------------

  subroutine SoluteMicporeGainXBoundaryM(M,N,M1,M2,M3,M4,M5,M6)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: M3,M2,M1
  integer, intent(in) :: M6,M5,M4
  integer :: idsalt, ids
  real(r8) :: VFLOW

  ! begin_execution

  VFLOW=WaterFlow2MicPM_3D(M,N,M6,M5,M4)
  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_MicpTranspFlxM_3D(idsalt,N,M6,M5,M4)=VFLOW*trcSalt_Irrig_vr(idsalt,M3,M2,M1)
  ENDDO

  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
    trcSalt_MicpTranspFlxM_3D(idsalt,N,M6,M5,M4) = VFLOW*trcSalt_Irrig_vr(idsalt,M3,M2,M1)*trcs_VLN_vr(ids_H1PO4,M3,M2,M1)
    trcSalt_MicpTranspFlxM_3D(ids,N,M6,M5,M4)    = VFLOW*trcSalt_Irrig_vr(idsalt,M3,M2,M1)*trcs_VLN_vr(ids_H1PO4B,M3,M2,M1)
  ENDDO
  end subroutine SoluteMicporeGainXBoundaryM  
!------------------------------------------------------------------------------------------

  subroutine AccumSnowDriftSaltFlxM(N,N1,N2,N4,N5,N4B,N5B,trcSalt_SnowDrift_flxM_col)
!
!     Description:
!
  implicit none
  integer, intent(in) :: N,N1,N2,N4,N5,N4B,N5B
  real(r8), intent(inout) :: trcSalt_SnowDrift_flxM_col(idsalt_beg:idsalt_end,JY,JX)
  integer :: idsalt,NN

!     begin_execution
!
  DO NN=1,2
    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_SnowDrift_flxM_col(idsalt,N2,N1)=trcSalt_SnowDrift_flxM_col(idsalt,N2,N1)+ &
        trcSalt_SnowDrift_flxM_2D(idsalt,N,NN,N2,N1)
    ENDDO
  ENDDO

  NN=iBehind
  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_SnowDrift_flxM_col(idsalt,N2,N1)=trcSalt_SnowDrift_flxM_col(idsalt,N2,N1)- &
      trcSalt_SnowDrift_flxM_2D(idsalt,N,NN,N5,N4)
  ENDDO  

  if(N4B.GT.0 .and. N5B.GT.0)then
    NN=iFront
    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_SnowDrift_flxM_col(idsalt,N2,N1)=trcSalt_SnowDrift_flxM_col(idsalt,N2,N1)- &
        trcSalt_SnowDrift_flxM_2D(idsalt,N,NN,N5B,N4B)
    ENDDO  
  endif
  
  end subroutine AccumSnowDriftSaltFlxM
!------------------------------------------------------------------------------------------

  subroutine SoluteMicporeLossXBoundaryM(M,N,M1,M2,M3,M4,M5,M6)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: M1,M2,M3
  integer, intent(in) :: M4,M5,M6
  real(r8) :: VFLW
  integer :: idsalt,ids
!     begin_execution

  IF(VLWatMicPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
    VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,M6,M5,M4)/VLWatMicPM_vr(M,M3,M2,M1)))
  ELSE
    VFLW=0.0_r8
  ENDIF

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_MicpTranspFlxM_3D(idsalt,N,M6,M5,M4)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,M3,M2,M1))
  ENDDO

  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
    trcSalt_MicpTranspFlxM_3D(idsalt,N,M6,M5,M4) = VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,M3,M2,M1))*trcs_VLN_vr(ids_H1PO4,M3,M2,M1)
    trcSalt_MicpTranspFlxM_3D(ids,N,M6,M5,M4)    = VFLW*AZMAX1(trcSalt_solml2_vr(ids,M3,M2,M1))*trcs_VLN_vr(ids_H1PO4B,M3,M2,M1)
  ENDDO
  end subroutine SoluteMicporeLossXBoundaryM
!------------------------------------------------------------------------------------------

  subroutine AccumRofSaltFlxM(M,N,N1,N2,N4,N5,N4B,N5B,trcSalt_FloXSurRof_flxM_col)
!
!     Description:
!
  integer , intent(in) :: M,N,N1,N2,N4,N5,N4B,N5B
  real(r8), intent(inout) :: trcSalt_FloXSurRof_flxM_col(idsalt_beg:idsalt_end,JY,JX)
  integer :: NN,idsalt
!  
!     begin_execution
!

  D1202: DO NN=1,2
    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_FloXSurRof_flxM_col(idsalt,N2,N1)=trcSalt_FloXSurRof_flxM_col(idsalt,N2,N1)+trcSalt_FloXSurRof_flxM_2DH(idsalt,N,NN,N2,N1)
    ENDDO

    IF(IFLBM_2DH(M,N,NN,N5,N4).EQ.0)THEN   !only true when NN=iBehind
      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_FloXSurRof_flxM_col(idsalt,N2,N1)=trcSalt_FloXSurRof_flxM_col(idsalt,N2,N1)-trcSalt_FloXSurRof_flxM_2DH(idsalt,N,NN,N5,N4)
      ENDDO
    ENDIF

    IF(N4B.GT.0 .AND. N5B.GT.0 .AND. NN.EQ.iFront)THEN
      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_FloXSurRof_flxM_col(idsalt,N2,N1)=trcSalt_FloXSurRof_flxM_col(idsalt,N2,N1)-trcSalt_FloXSurRof_flxM_2DH(idsalt,N,NN,N5B,N4B)
      ENDDO
    ENDIF
  ENDDO D1202
  end subroutine AccumRofSaltFlxM
!------------------------------------------------------------------------------------------

  subroutine AccumSaltTranspInSnowM(M,NY,NX)
!
!     Accumulate accepted salt transfers through every resolved snow store.
!
  integer, intent(in) :: M,NY,NX
  integer :: LS,LS2,idsalt
!     begin_execution

  D1205: DO LS=1,JS
    IF(VLSnowHeatCapM_snvr(M,LS,NY,NX).GT.ZEROS(NY,NX))THEN
      LS2=MIN(JS,LS+1)
      ! IF LOWER LAYER IS IN THE SNOWPACK
      IF(LS.LT.JS.AND.VLSnowHeatCapM_snvr(M,LS2,NY,NX).GT.ZEROS(NY,NX))THEN
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_Aqua_flxM_snvr(idsalt,LS,NY,NX)=trcSalt_Aqua_flxM_snvr(idsalt,LS,NY,NX) &
            +trcSalt_AquaAdv_flxM_snvr(idsalt,LS,NY,NX)-trcSalt_AquaAdv_flxM_snvr(idsalt,LS2,NY,NX)
        ENDDO
        !     IF LOWER LAYER IS THE LITTER AND SOIL SURFACE        
      ELSE
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_Aqua_flxM_snvr(idsalt,LS,NY,NX)=trcSalt_Aqua_flxM_snvr(idsalt,LS,NY,NX)+trcSalt_AquaAdv_flxM_snvr(idsalt,LS,NY,NX) &
            -trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX)-trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)
        ENDDO

        DO idsalt=0,idsalt_nuts
          trcSalt_Aqua_flxM_snvr(idsalt_H0PO4+idsalt,LS,NY,NX)=trcSalt_Aqua_flxM_snvr(idsalt_H0PO4+idsalt,LS,NY,NX) &
            -trcSalt_AquaADV_Snow2Soil_flxM(idsalt_H0PO4B+idsalt,NY,NX)
        ENDDO
      ENDIF
    ELSEIF(LS.EQ.1)THEN
      ! First snowfall has no initial snow heat capacity. Accept its staged
      ! deposition once; donor-dependent discharge remains in the branch above.
      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_Aqua_flxM_snvr(idsalt,LS,NY,NX)=trcSalt_Aqua_flxM_snvr(idsalt,LS,NY,NX) &
          +trcSalt_AquaAdv_flxM_snvr(idsalt,LS,NY,NX)
      ENDDO
    ENDIF
  ENDDO D1205
  end subroutine AccumSaltTranspInSnowM

!------------------------------------------------------------------------------------------

  subroutine SnowSaltTranptM(M,NY,NX,band_snow_donor)
!
!     Route salt with drainage from resolved snow layers to snow and soil.
!
  implicit none
  integer, intent(in) :: M,NY,NX
  ! Present in the surface band pass; returns zero if there is no snow donor.
  integer, optional, intent(out) :: band_snow_donor

  character(len=*), parameter :: subname='SnowSaltTranptM'
  integer :: L,L2,idsalt,ids
  logical :: ICHKL
  real(r8) :: VFLWW   !fraction flowed into next snow layer
  real(r8) :: VFLWR,VFLWS,VFLWPO4,VFLWPOB
  
!     begin_execution
!
  call PrintInfo('beg '//subname)
  if(present(band_snow_donor))band_snow_donor=0
  ! Reset column discharge once; inactive layers must preserve bottom-layer fluxes.
  trcSalt_AquaADV_Snow2Litr_flxM(idsalt_beg:idsalt_end,NY,NX)=0._r8
  trcSalt_AquaADV_Snow2Soil_flxM(idsalt_beg:idsaltb_end,NY,NX)=0._r8
  ICHKL=.false.
  DO L=1,JS
    !
    IF(VLSnowHeatCapM_snvr(M,L,NY,NX).GT.ZEROS(NY,NX))THEN
      L2=MIN(JS,L+1)

      !Not bottom layer, next layer has significant snow
      IF(L.LT.JS .AND. VLSnowHeatCapM_snvr(M,L2,NY,NX).GT.ZEROS(NY,NX))THEN
        IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          VFLWW=AZMAX1(AMIN1(1.0_r8,WatFlowInSnowM_snvr(M,L2,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
        ELSE
          ! Intentional numerical redistribution for negligible liquid:
          ! transfer a small tracer fraction even when water flow is zero.
          ! This is not water-driven advection; its cumulative effect
          ! depends on the number of transport substeps.
          VFLWW=1.e-8_r8
        ENDIF

        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_AquaAdv_flxM_snvr(idsalt,L2,NY,NX) = trcSalt_ml2_snvr(idsalt,L,NY,NX)*VFLWW
        ENDDO

        !bottom snow layer or L2 has insignificant snow layer
      ELSE

        !Insignificant snow layer
        IF(L.LT.JS)THEN
          trcSalt_AquaAdv_flxM_snvr(idsalt_beg:idsalt_end,L2,NY,NX)=0.0_r8
        ENDIF
        !
        ! SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
        !
        IF(.not.ICHKL)THEN
          !flow into soil
          IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
            VFLWR=AZMAX1(AMIN1(1.0_r8,WatFlowSno2LitRM_col(M,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
            VFLWS=AZMAX1(AMIN1(1.0_r8,(WatFlowSno2MicPM_col(M,NY,NX)+WatFlowSno2MacPM_col(M,NY,NX))/VLWatSnow_snvr(L,NY,NX)))
          ELSE
            VFLWR=FracSurfByLitR_col(NY,NX)
            VFLWS=FracSurfBareSoil_col(NY,NX)
          ENDIF
          VFLWPO4=VFLWS*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
          VFLWPOB=VFLWS*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)

          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX) = trcSalt_ml2_snvr(idsalt,L,NY,NX)*VFLWR
          ENDDO

          DO idsalt=idsalt_beg,idsalt_KSO4
            trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX) = trcSalt_ml2_snvr(idsalt,L,NY,NX)*VFLWS
          ENDDO

          DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
            ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B  
            !soil
            trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX) = trcSalt_ml2_snvr(idsalt,L,NY,NX)*VFLWPO4
            !band
            trcSalt_AquaADV_Snow2Soil_flxM(ids,NY,NX) = trcSalt_ml2_snvr(idsalt,L,NY,NX)*VFLWPOB
          ENDDO
          if(present(band_snow_donor))band_snow_donor=L
          ICHKL=.true.
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  ! Snow-to-band transfer belongs only to the separate surface pass.
  if(.not.present(band_snow_donor))then
    trcSalt_AquaADV_Snow2Soil_flxM(idsalt_pband_beg:idsalt_pband_end,NY,NX)=0._r8
  endif
  call PrintInfo('end '//subname)
  end subroutine SnowSaltTranptM

!------------------------------------------------------------------------------------------

  subroutine Litr2SoilDifusTransptM(M,NY,NX,trcSalt_flx_diffusM)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8),intent(out) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)

  character(len=*), parameter :: subname='Litr2SoilDifusTransptM'
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcSaltDiffConductance(idsalt_beg:idsalt_KSO4)
  real(r8) :: trcSalt_conc_src(idsalt_beg:idsaltb_end)    !tracer solute concentration in source cell
  real(r8) :: trcSalt_conc_dest(idsalt_beg:idsaltb_end)   !tracer solutte concentration in destination cell
  real(r8) :: VOLWPB,VOLWPA,TORT0,TORT1
  real(r8) :: DLYR0,FLWRM1
  integer :: idsalt
!     begin_execution
!
  call PrintInfo('beg '//subname)

  FLWRM1=WatFLoLitr2SoilM_col(M,NY,NX)

  IF((VGeomLayer_vr(0,NY,NX).GT.ZEROS(NY,NX) .AND. VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
    .AND. (VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN

    VOLWPA=VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
    VOLWPB=VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)

    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,0,NY,NX)/VLWatMicPM_vr(M,0,NY,NX))
    ENDDO

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX))
    ENDDO

    IF(VOLWPA.GT.ZEROS2(NY,NX))THEN
      DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX)/VOLWPA)
      ENDDO
    ELSE
      trcSalt_conc_dest(idsalt_psoiL_beg:idsalt_psoil_end)=0.0_r8
    ENDIF
    
    IF(VOLWPB.GT.ZEROS2(NY,NX))THEN
      DO idsalt=idsalt_pband_beg,idsalt_pband_end
        trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX)/VOLWPB)
      ENDDO
    ELSE
      DO idsalt=0,idsalt_nuts
        trcSalt_conc_dest(idsalt_H0PO4B+idsalt)=trcSalt_conc_dest(idsalt_H0PO4+idsalt)
      ENDDO
    ENDIF
  !
  !     DIFFUSIVITIES IN RESIDUE AND SOIL SURFACE
  !
  !
    DLYR0 = AMAX1(ZERO2,DLYR_3D(3,0,NY,NX))
    TORT0 = TortMicPM_vr(M,0,NY,NX)*FracSurfByLitR_col(NY,NX)
    DLYR1 = AMAX1(ZERO2,DLYR_3D(3,NU_col(NY,NX),NY,NX))
    TORT1 = TortMicPM_vr(M,NU_col(NY,NX),NY,NX)
    TORTL = AMIN1(1.0_r8,(TORT0+TORT1)/(DLYR0+DLYR1))
    DISPN = (0.20_r8*(0.5_r8*DLYR1)**1.07_r8)*AMIN1(VFLWX,ABS(FLWRM1/AREA_3D(3,NU_col(NY,NX),NY,NX)))

    DIFPO=(SoluteDifusivitytscaledM_vr(NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)

    DIFAL=(AquaIonDifusivty2_vr(idsalt_Al,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance((/idsalt_Al,idsalt_AlOH,idsalt_AlOH2,idsalt_AlOH3,idsalt_AlOH4,idsalt_AlSO4/))=DIFAL

    DIFFE=(AquaIonDifusivty2_vr(idsalt_Fe,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance((/idsalt_Fe,idsalt_FeOH,idsalt_FeOH2,idsalt_FeOH3,idsalt_FeOH4,idsalt_FeSO4/))=DIFFE

    trcSaltDiffConductance(idsalt_Hp)=(AquaIonDifusivty2_vr(idsalt_Hp,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)

    DIFCA=(AquaIonDifusivty2_vr(idsalt_Ca,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance((/idsalt_Ca,idsalt_CaOH,idsalt_CaCO3,idsalt_CaHCO3,idsalt_CaSO4/))=DIFCA

    DIFMG=(AquaIonDifusivty2_vr(idsalt_Mg,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance((/idsalt_Mg,idsalt_MgOH2,idsalt_MgCO3,idsalt_MgHCO3,idsalt_MgSO4/))=DIFMG

    DIFNA=(AquaIonDifusivty2_vr(idsalt_Na,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance((/idsalt_Na,idsalt_NaCO3,idsalt_NaSO4/))=DIFNA

    DIFKA=(AquaIonDifusivty2_vr(idsalt_K,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance((/idsalt_K,idsalt_KSO4/))=DIFKA

    trcSaltDiffConductance(idsalt_OH)=(AquaIonDifusivty2_vr(idsalt_OH,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance(idsalt_SO4)=(AquaIonDifusivty2_vr(idsalt_SO4,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance(idsalt_Cl)=(AquaIonDifusivty2_vr(idsalt_Cl,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance(idsalt_CO3)=(AquaIonDifusivty2_vr(idsalt_CO3,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    trcSaltDiffConductance(idsalt_HCO3)=(AquaIonDifusivty2_vr(idsalt_HCO3,NU_col(NY,NX),NY,NX)*TORTL+DISPN)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    !
    !     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
    !     MICROPORES
    !
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_flx_diffusM(idsalt)=trcSaltDiffConductance(idsalt)*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_flx_diffusM(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt)) &
        *trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
    ENDDO

    DO idsalt=0,idsalt_nuts
      trcSalt_flx_diffusM(idsalt_H0PO4B+idsalt)=DIFPO*(trcSalt_conc_src(idsalt_H0PO4+idsalt)&
        -trcSalt_conc_dest(idsalt_H0PO4B+idsalt))*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)
    ENDDO
  else
      trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)=0.0_r8  
  endif
  call PrintInfo('end '//subname)
  end subroutine Litr2SoilDifusTransptM
!------------------------------------------------------------------------------------------

  subroutine LitrSoilAdvDifusTransptM(NY,NX,trcSalt_Adv2MicP_flx,trcSalt_flx_diffusM,trcSalt_Adv2MacP_flx,snow_mac_frac)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: snow_mac_frac ! fraction of soil-directed snow liquid entering macropores
  real(r8), intent(in) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSalt_Adv2MacP_flx(idsalt_beg:idsaltb_end)
  integer :: idsalt,ids

  !     begin_execution
  !
  !
  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_MicpTranspFlxM_3D(idsalt,3,0,NY,NX)=trcSalt_Precip2LitrM(idsalt,NY,NX)+trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX) &
      -trcSalt_Adv2MicP_flx(idsalt)-trcSalt_flx_diffusM(idsalt)
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MicpTranspFlxM_3D(idsalt,3,NU_col(NY,NX),NY,NX)=trcSalt_Precip2MicpM(idsalt,NY,NX) &
      +trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)*(1._r8-snow_mac_frac)+trcSalt_Adv2MicP_flx(idsalt)+trcSalt_flx_diffusM(idsalt)
  ENDDO

  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_MicpTranspFlxM_3D(idsalt,3,0,NY,NX)=trcSalt_MicpTranspFlxM_3D(idsalt,3,0,NY,NX)-trcSalt_Adv2MacP_flx(idsalt)
  ENDDO
  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
    ! Litter has no separate nutrient band. Both soil destinations draw
    ! from (or return to) the same unbanded litter phosphate species.
    trcSalt_MicpTranspFlxM_3D(idsalt,3,0,NY,NX)=trcSalt_MicpTranspFlxM_3D(idsalt,3,0,NY,NX) &
      -trcSalt_Adv2MicP_flx(ids)-trcSalt_flx_diffusM(ids)-trcSalt_Adv2MacP_flx(ids)
  ENDDO
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MacpTranspFlxM_3D(idsalt,3,NU_col(NY,NX),NY,NX)=trcSalt_MacpTranspFlxM_3D(idsalt,3,NU_col(NY,NX),NY,NX)+trcSalt_Adv2MacP_flx(idsalt) &
      +trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)*snow_mac_frac
  ENDDO
  end subroutine LitrSoilAdvDifusTransptM
!------------------------------------------------------------------------------------------

  subroutine AccumLitrSoilSaltTranspFluxM(NY,NX,trcSalt_Adv2MicP_flx,trcSalt_flx_diffusM)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  REAL(R8), intent(in) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  real(r8), intent(in) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)

  integer :: idsalt
  !
  !     begin_execution
  !     ACCUMULATE HOURLY FLUXES FOR USE IN REDIST.F
  !
  DO idsalt=idsalt_beg,idsalt_end
    trcSalt_TransptMicP_3D(idsalt,3,0,NY,NX)=trcSalt_TransptMicP_3D(idsalt,3,0,NY,NX)+trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX) &
      -trcSalt_Adv2MicP_flx(idsalt)-trcSalt_flx_diffusM(idsalt)
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_TransptMicP_3D(idsalt,3,NU_col(NY,NX),NY,NX)=trcSalt_TransptMicP_3D(idsalt,3,NU_col(NY,NX),NY,NX) &
      +trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)+trcSalt_Adv2MicP_flx(idsalt)+trcSalt_flx_diffusM(idsalt)
  ENDDO
  end subroutine AccumLitrSoilSaltTranspFluxM
!------------------------------------------------------------------------------------------

  subroutine MacPoreSaltAdvExM(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: trcSalt_Adv2MacP_flx(idsalt_beg:idsaltb_end)
  real(r8) :: VFLW
  integer :: idsalt
!     begin_execution

  IF(VLWatMacPM_vr(M,NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM_vr(M,NU_col(NY,NX),NY,NX)/VLWatMacPM_vr(M,NU_col(NY,NX),NY,NX)))
  ELSE
    VFLW=VFLWX
  ENDIF

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_Adv2MacP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,NU_col(NY,NX),NY,NX))
  ENDDO

  DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_Adv2MacP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,NU_col(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
  ENDDO

  DO idsalt=idsalt_pband_beg,idsalt_pband_end
    trcSalt_Adv2MacP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,NU_col(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flxM_vr(idsalt,NU_col(NY,NX),NY,NX)=trcSalt_Adv2MacP_flx(idsalt)
  ENDDO
  end subroutine MacPoreSaltAdvExM
!------------------------------------------------------------------------------------------

  subroutine MicPoreSaltAdvExM(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  integer :: idsalt
  real(r8) :: VFLW

!     begin_execution
  IF(VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM_vr(M,NU_col(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)))
  ELSE
    VFLW=-VFLWX
  ENDIF

  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX))
  ENDDO

  DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
    trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
  ENDDO

  DO idsalt=idsalt_pband_beg,idsalt_pband_end
    trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX)) &
      *trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)
  ENDDO

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flxM_vr(idsalt,NU_col(NY,NX),NY,NX)=trcSalt_Adv2MicP_flx(idsalt)
  ENDDO
  end subroutine MicPoreSaltAdvExM

!------------------------------------------------------------------------------------------

  subroutine MacMicPoreSaltDifusExM(M,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8)  :: trcSalt_Difus_Mac2Micp_flxM(idsalt_beg:idsaltb_end)
  real(r8) :: VOLWHS,VOLWT,VLWatMicPMNU
  integer :: idsalt
!     begin_execution
!

  IF(VLWatMacPM_vr(M,NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VLWatMicPMNU = VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)
    VOLWHS       = AMIN1(XFRS*VGeomLayer_vr(NU_col(NY,NX),NY,NX),VLWatMacPM_vr(M,NU_col(NY,NX),NY,NX))
    VOLWT        = VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)+VOLWHS

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,NU_col(NY,NX),NY,NX))*VLWatMicPMNU &
        -AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX))*VOLWHS)/VOLWT
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,NU_col(NY,NX),NY,NX))*VLWatMicPMNU &
        -AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,NU_col(NY,NX),NY,NX))*VLWatMicPMNU &
        -AZMAX1(trcSalt_solml2_vr(idsalt,NU_col(NY,NX),NY,NX))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)
    ENDDO
  !
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_Mac2MicPore_flxM_vr(idsalt,NU_col(NY,NX),NY,NX)=trcSalt_Mac2MicPore_flxM_vr(idsalt,NU_col(NY,NX),NY,NX)+trcSalt_Difus_Mac2Micp_flxM(idsalt)
    ENDDO
  else
   trcSalt_Difus_Mac2Micp_flxM(idsalt_beg:idsaltb_end)=0.0_r8
  endif
  end subroutine MacMicPoreSaltDifusExM

!------------------------------------------------------------------------------------------

  subroutine SaltFlowThruLitrRofM(M,N1,N2)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,N1,N2
  real(r8) :: VFLW
  INTEGER :: idsalt

  !     begin_execution
  !flow out of (N2,N1)
  IF(SurfRunoffPotentM_col(M,N2,N1).GT.ZEROS(N2,N1))THEN
    IF(VLWatMicPM_vr(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AMIN1(VFLWX,SurfRunoffPotentM_col(M,N2,N1)/VLWatMicPM_vr(M,0,N2,N1))
    ELSE
      VFLW=VFLWX
    ENDIF

    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_FloXSurRof_flxM(idsalt,N2,N1)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,0,N2,N1))
    ENDDO
  else
    trcSalt_FloXSurRof_flxM(idsalt_beg:idsalt_end,N2,N1)=0.0_r8
  endif  
  end subroutine SaltFlowThruLitrRofM

!------------------------------------------------------------------------------------------

  subroutine SaltSurfRofM(M,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='SaltSurfRofM'
  integer :: NY,NX,N2,N1,N,N3,N4,N5,N4B,N5B,NN
  real(r8) :: FQRM
  integer :: idsalt
  call PrintInfo('beg '//subname)

  trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,:,:,:,:)=0.0_r8
  !  
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      N2=NY; N1=NX
      IF(SurfRunoffPotentM_col(M,N2,N1).GT.ZEROS(N2,N1))THEN  
        DO N=1,2
          DO  NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              IF(NX.EQ.NHE .AND. NN.EQ.iFront .OR. &  !skip eastern boundary
                (NX.EQ.NHW .AND. NN.EQ.iBehind))THEN  !skip western boundary
                cycle
              ELSE
                N4  = NX+1; N5  = NY   !eastward
                N4B = NX-1;N5B = NY    !westward
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN
              IF(NY.EQ.NVS .AND. NN.EQ.iFront       &      !skip southern boundary
                .OR. (NY.EQ.NVN .AND. NN.EQ.iBehind))THEN  !skip northern boundary
                cycle
              ELSE
                N4  = NX; N5  = NY+1   !southward
                N4B = NX; N5B = NY-1   !northward
              ENDIF
            ENDIF
            !
            ! current grid has outflow 
            IF(SurfRunoffPotentM_col(M,N2,N1).GT.ZEROS(N2,N1))THEN
              IF(NN.EQ.iFront)THEN
                FQRM=QflxSurfRunoffM_2DH(M,N,iBehind,N5,N4)/SurfRunoffPotentM_col(M,N2,N1)
                DO idsalt=idsalt_beg,idsalt_end
                  trcSalt_FloXSurRof_flxM_2DH(idsalt,N,iBehind,N5,N4)=trcSalt_FloXSurRof_flxM(idsalt,N2,N1)*FQRM
                ENDDO
                !
              ELSEIF(NN.EQ.iBehind)THEN
                !legitimate inner grid
                IF(N4B.GT.0 .AND. N5B.GT.0)THEN
                  FQRM=QflxSurfRunoffM_2DH(M,N,iFront,N5B,N4B)/SurfRunoffPotentM_col(M,N2,N1)
                  DO idsalt=idsalt_beg,idsalt_end
                    trcSalt_FloXSurRof_flxM_2DH(idsalt,N,iFront,N5B,N4B)=trcSalt_FloXSurRof_flxM(idsalt,N2,N1)*FQRM
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  END subroutine SaltSurfRofM
!------------------------------------------------------------------------------------------

  subroutine SaltFlowThruSnowDriftM(M,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='SaltFlowThruSnowDriftM'
  integer :: N,NN,N4,N5,N4B,N5B,idsalt
  real(r8) :: FQRM,VFLW
  integer :: N1,N2,NY,NX

!     begin_execution

  call PrintInfo('beg '//subname)
  trcSalt_SnowDrift_flxM_2D(idsalt_beg:idsaltb_end,:,:,:,:)=0.0_r8  
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      N2=NY; N1=NX
      DO N=1,2
        DO  NN=1,2
          IF(N.EQ.iWestEastDirection)THEN
            IF(NX.EQ.NHE .AND. NN.EQ.iFront      &       !skip eastern boundary
              .OR. (NX.EQ.NHW .AND. NN.EQ.iBehind))THEN  !skip western boundary
              cycle
            ELSE
              N4  = NX+1  ;N5  = NY;!eastward
              N4B = NX-1  ;N5B = NY !westward
            ENDIF
          ELSEIF(N.EQ.iNorthSouthDirection)THEN
            IF(NY.EQ.NVS .AND. NN.EQ.iFront    &         !skip southern boundary
              .OR. (NY.EQ.NVN .AND. NN.EQ.iBehind))THEN  !skip northern boundary
              cycle
            ELSE
              N4  = NX;N5  = NY+1   !southward
              N4B = NX;N5B = NY-1   !northward
            ENDIF
          ENDIF
!
          IF(NN.EQ.iFront)THEN
            IF(DrySnoFlxByRedistM_2DH(M,N,iBehind,N5,N4).GT.0.0_r8)THEN
              IF(VcumSnoDWI_col(N2,N1).GT.ZEROS2(NY,NX))THEN
                VFLW=AZMAX1(AMIN1(VFLWX,DrySnoFlxByRedistM_2DH(M,N,iBehind,N5,N4)/VcumSnoDWI_col(N2,N1)))
              ELSE
                VFLW=VFLWX
              ENDIF

              ! Snow contains only unbanded salt species.
              DO idsalt=idsalt_beg,idsalt_end
                trcSalt_SnowDrift_flxM_2D(idsalt,N,iBehind,N5,N4) = VFLW*AZMAX1(trcSalt_ml2_snvr(idsalt,1,N2,N1))
              ENDDO
            ENDIF  
          ELSEIF(NN.EQ.iBehind .and. N4B.GT.0 .and. N5B.GT.0)then
            IF(DrySnoFlxByRedistM_2DH(M,N,iFront,N5B,N4B).GT.0.0_r8)THEN
              IF(VcumSnoDWI_col(N2,N1).GT.ZEROS2(NY,NX))THEN
                VFLW=AZMAX1(AMIN1(VFLWX,DrySnoFlxByRedistM_2DH(M,N,iFront,N5B,N4B)/VcumSnoDWI_col(N2,N1)))
              ELSE
                VFLW=VFLWX
              ENDIF
              DO idsalt=idsalt_beg,idsalt_end
                trcSalt_SnowDrift_flxM_2D(idsalt,N,iFront,N5B,N4B) = VFLW*AZMAX1(trcSalt_ml2_snvr(idsalt,1,N2,N1))
              ENDDO
            ENDIF  
          endif
        ENDDO
      ENDDO    
    enddo
  enddo
  call PrintInfo('end '//subname)
  end subroutine SaltFlowThruSnowDriftM
!------------------------------------------------------------------------------------------
  subroutine SoluteAdvMicroporeM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  character(len=*), parameter :: subname='SoluteAdvMicroporeM'
  real(r8) :: VFLW
  integer :: idsalt
  real(r8) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  
  call PrintInfo('beg '//subname)  

  !flow out of grid (N3,N2,N1)
  IF(WaterFlow2MicPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4).GT.0.0)THEN
    !
    IF(VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/VLWatMicPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF

    DO idsalt=idsalt_beg,idsalt_KSO4
        trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1))*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1))*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
    ENDDO    
    !flow out of grid (N6,N5,N4)
  ELSE
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MicPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      ! Intentional numerical redistribution of residual tracers in dry soil,
      ! including when the face water flow is zero.
      ! TODO: revisit the fallback magnitude and timestep dependence of VFLWX
      ! when revising dry-soil transport.
      VFLW=-VFLWX
    ENDIF
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ENDIF

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MicpTranspFlxM_3D(idsalt,N,SoilFaceSlot(N,N3,N6),N5,N4)=trcSalt_Adv2MicP_flx(idsalt)
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SoluteAdvMicroporeM
!------------------------------------------------------------------------------------------
  subroutine SoluteDifsMicroporeM(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2,N3
  integer, intent(in) :: N4,N5,N6
  real(r8) :: FaceArea,AreaOverDistance,Dispersivity
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  REAL(R8) :: VLWPA1,VLWPB1,VLWPA2,VLWPB2
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcSaltDiffConductance(idsalt_beg:idsalt_KSO4)
  real(r8) :: trcSalt_conc_src(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_conc_dest(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)
  integer :: idsalt

  !test both grids are legitimate
  IF(THETW1(N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1) &
    .AND.THETW1(N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4) &
    .AND.VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1) &
    .AND.VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
!

    VLWPA1=VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    VLWPB1=VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
    VLWPA2=VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    VLWPB2=VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
    ENDDO

    IF(VLWPA1.GT.ZEROS(N2,N1))THEN
      DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1)/VLWPA1)
      ENDDO
    ELSE
      trcSalt_conc_src(idsalt_psoil_beg:idsalt_psoil_end)=0.0_r8
    ENDIF

    IF(VLWPB1.GT.ZEROS(N2,N1))THEN
      DO idsalt=idsalt_pband_beg,idsalt_pband_end
        trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N3,N2,N1)/VLWPB1)
      ENDDO
    ELSE
      DO idsalt=0,idsalt_nuts
        trcSalt_conc_src(idsalt_H0PO4B+idsalt)=trcSalt_conc_src(idsalt_H0PO4+idsalt)
      ENDDO
    ENDIF

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
    ENDDO

    IF(VLWPA2.GT.ZEROS(N5,N4))THEN
      DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4)/VLWPA2)
      ENDDO
    ELSE
      trcSalt_conc_dest(idsalt_psoil_beg:idsalt_psoil_end)=0.0_r8
    ENDIF
    IF(VLWPB2.GT.ZEROS(N5,N4))THEN
      DO idsalt=idsalt_pband_beg,idsalt_pband_end
        trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4)/VLWPB2)
      ENDDO
    ELSE
      DO idsalt=0,idsalt_nuts
        trcSalt_conc_dest(idsalt_H0PO4B+idsalt)=trcSalt_conc_dest(idsalt_H0PO4+idsalt)
      ENDDO
    ENDIF
!
!     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
!

    DLYR1=AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2=AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL=(TortMicPM_vr(M,N3,N2,N1)*DLYR1+TortMicPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)

    call GetSoilFaceGeometry(N,N1,N2,N3,N4,N5,N6,FaceArea,AreaOverDistance,Dispersivity)
    DISPN=Dispersivity*AMIN1(VFLWX,ABS(WaterFlow2MicPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/FaceArea))
    DIFPO=(SoluteDifusivitytscaledM_vr(N6,N5,N4)*TORTL+DISPN)*AreaOverDistance

    DIFAL=(AquaIonDifusivty2_vr(idsalt_Al,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Al,idsalt_AlOH,idsalt_AlOH2,idsalt_AlOH3,idsalt_AlOH4,idsalt_AlSO4/))=DIFAL

    DIFFE=(AquaIonDifusivty2_vr(idsalt_Fe,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Fe,idsalt_FeOH,idsalt_FeOH2,idsalt_FeOH3,idsalt_FeOH4,idsalt_FeSO4/))=DIFFE

    trcSaltDiffConductance(idsalt_Hp)=(AquaIonDifusivty2_vr(idsalt_Hp,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance

    DIFCA=(AquaIonDifusivty2_vr(idsalt_Ca,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Ca,idsalt_CaOH,idsalt_CaCO3,idsalt_CaHCO3,idsalt_CaSO4/))=DIFCA

    DIFMG=(AquaIonDifusivty2_vr(idsalt_Mg,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Mg,idsalt_MgOH2,idsalt_MgCO3,idsalt_MgHCO3,idsalt_MgSO4/))=DIFMG

    DIFNA=(AquaIonDifusivty2_vr(idsalt_Na,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Na,idsalt_NaCO3,idsalt_NaSO4/))=DIFNA

    DIFKA=(AquaIonDifusivty2_vr(idsalt_K,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_K,idsalt_KSO4/))=DIFKA

    trcSaltDiffConductance(idsalt_OH)=(AquaIonDifusivty2_vr(idsalt_OH,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance(idsalt_SO4)=(AquaIonDifusivty2_vr(idsalt_SO4,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance(idsalt_Cl)=(AquaIonDifusivty2_vr(idsalt_Cl,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance(idsalt_CO3)=(AquaIonDifusivty2_vr(idsalt_CO3,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance(idsalt_HCO3)=(AquaIonDifusivty2_vr(idsalt_HCO3,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    !
    !     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
    !     MICROPORES
    !
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_flx_diffusM(idsalt)=trcSaltDiffConductance(idsalt)*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoiL_end
      trcSalt_flx_diffusM(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_flx_diffusM(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcSalt_flx_diffusM(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MicpTranspFlxM_3D(idsalt,N,SoilFaceSlot(N,N3,N6),N5,N4)=trcSalt_MicpTranspFlxM_3D(idsalt,N,SoilFaceSlot(N,N3,N6),N5,N4)+trcSalt_flx_diffusM(idsalt)
  ENDDO

  end subroutine SoluteDifsMicroporeM
!----------------------------------------------------------------------
  subroutine SoluteAdvMacroporeM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  character(len=*), parameter :: subname='SoluteAdvMacroporeM'
  integer :: idsalt
  real(r8) :: trcSalt_RFH(idsalt_beg:idsaltb_end)
  real(r8) :: VFLW


  call PrintInfo('beg '//subname)

  !flow out of (N3,N2,N1)
  IF(WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4).GT.0.0)THEN
    !
    IF(VLWatMacPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/VLWatMacPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    !
    !     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE IN VERTICAL FLUX
    !
    IF(N.EQ.3.AND.VLMacroPores_vr(N6,N5,N4).GT.VLWatMacPM_vr(M,N6,N5,N4))THEN
      DO idsalt=idsalt_beg,idsalt_KSO4
        trcSalt_RFH(idsalt)=VFLW*AZMAX1((trcSalt_soHml2_vr(idsalt,N3,N2,N1) &
          -AZMIN1(trcSalt_Mac2MicPore_flxM_vr(idsalt,N3,N2,N1))))
      ENDDO

      DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_RFH(idsalt)=VFLW*AZMAX1((trcSalt_soHml2_vr(idsalt,N3,N2,N1) &
          -AZMIN1(trcSalt_Mac2MicPore_flxM_vr(idsalt,N3,N2,N1))))*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
      ENDDO

      DO idsalt=idsalt_pband_beg,idsalt_pband_end
        trcSalt_RFH(idsalt)=VFLW*AZMAX1((trcSalt_soHml2_vr(idsalt,N3,N2,N1) &
          -AZMIN1(trcSalt_Mac2MicPore_flxM_vr(idsalt,N3,N2,N1))))*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
      ENDDO
      !
      !     OTHERWISE
      !
    ELSE
      DO idsalt=idsalt_beg,idsalt_KSO4
        trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N3,N2,N1))
      ENDDO

      DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
        trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N3,N2,N1)) &
          *trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
      ENDDO

      DO idsalt=idsalt_pband_beg,idsalt_pband_end
        trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N3,N2,N1)) &
          *trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
      ENDDO
    ENDIF
    !flow out of (N6,N5,N4)
  ELSEIF(WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4).LT.0.0)THEN
   !
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4)) &
        *trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_RFH(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4)) &
        *trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
    !     NO MACROPORE FLUX
    !
  ELSE
    trcSalt_RFH(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF
  !
  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MacpTranspFlxM_3D(idsalt,N,SoilFaceSlot(N,N3,N6),N5,N4)=trcSalt_RFH(idsalt)
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SoluteAdvMacroporeM

!----------------------------------------------------------------------
  subroutine SoluteDifsMacroporeM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  character(len=*), parameter :: subname='SoluteDifsMacroporeM'
  real(r8) :: FaceArea,AreaOverDistance,Dispersivity
  real(r8) :: DLYR1,DLYR2,TORTL,DISPN,DIFPO,DIFAL,DIFFE,DIFCA,DIFMG,DIFNA,DIFKA
  real(r8) :: trcSaltDiffConductance(idsalt_beg:idsalt_KSO4)
  real(r8) :: trcSalt_MacPore_Diffus(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_conc_src(idsalt_beg:idsaltb_end)
  real(r8) :: trcSalt_conc_dest(idsalt_beg:idsaltb_end)
  integer  :: idsalt

  call PrintInfo('beg '//subname)
  !test grid legitimacy
  IF(VLWatMacPM_vr(M,N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1)*VLMacroPores_vr(N3,N2,N1) &
    .AND.VLWatMacPM_vr(M,N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4)*VLMacroPores_vr(N6,N5,N4))THEN
    !
    !     MACROPORE CONCENTRATIONS IN CURRENT AND ADJACENT GRID CELLS
    !
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_conc_src(idsalt)=AZMAX1(trcSalt_soHml2_vr(idsalt,N3,N2,N1)/VLWatMacPM_vr(M,N3,N2,N1))
      trcSalt_conc_dest(idsalt)=AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4))
    ENDDO
    !
    !     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
    !

    DLYR1 = AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2 = AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL = (TortMacPM_vr(M,N3,N2,N1)*DLYR1+TortMacPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    call GetSoilFaceGeometry(N,N1,N2,N3,N4,N5,N6,FaceArea,AreaOverDistance,Dispersivity)
    DISPN = Dispersivity*AMIN1(VFLWX,ABS(WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/FaceArea))
    DIFPO = (SoluteDifusivitytscaledM_vr(N6,N5,N4)*TORTL+DISPN)*AreaOverDistance

    DIFAL=(AquaIonDifusivty2_vr(idsalt_Al,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Al,idsalt_AlOH,idsalt_AlOH2,idsalt_AlOH3,idsalt_AlOH4,idsalt_AlSO4/))=DIFAL

    DIFFE=(AquaIonDifusivty2_vr(idsalt_Fe,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Fe,idsalt_FeOH,idsalt_FeOH2,idsalt_FeOH3,idsalt_FeOH4,idsalt_FeSO4/))=DIFFE

    trcSaltDiffConductance(idsalt_Hp)=(AquaIonDifusivty2_vr(idsalt_Hp,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance

    DIFCA=(AquaIonDifusivty2_vr(idsalt_Ca,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Ca,idsalt_CaOH,idsalt_CaCO3,idsalt_CaHCO3,idsalt_CaSO4/))=DIFCA

    DIFMG=(AquaIonDifusivty2_vr(idsalt_Mg,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Mg,idsalt_MgOH2,idsalt_MgCO3,idsalt_MgHCO3,idsalt_MgSO4/))=DIFMG

    DIFNA=(AquaIonDifusivty2_vr(idsalt_Na,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_Na,idsalt_NaCO3,idsalt_NaSO4/))=DIFNA

    DIFKA=(AquaIonDifusivty2_vr(idsalt_K,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance((/idsalt_K,idsalt_KSO4/))=DIFKA

    trcSaltDiffConductance(idsalt_OH)=(AquaIonDifusivty2_vr(idsalt_OH,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance(idsalt_SO4)=(AquaIonDifusivty2_vr(idsalt_SO4,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance(idsalt_Cl)=(AquaIonDifusivty2_vr(idsalt_Cl,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance(idsalt_CO3)=(AquaIonDifusivty2_vr(idsalt_CO3,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    trcSaltDiffConductance(idsalt_HCO3)=(AquaIonDifusivty2_vr(idsalt_HCO3,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    !
    !     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
    !     MACROPORES
    !
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_MacPore_Diffus(idsalt)=trcSaltDiffConductance(idsalt)*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_MacPore_Diffus(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_MacPore_Diffus(idsalt)=DIFPO*(trcSalt_conc_src(idsalt)-trcSalt_conc_dest(idsalt))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcSalt_MacPore_Diffus(idsalt_beg:idsaltb_end)=0._r8
  ENDIF

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_MacpTranspFlxM_3D(idsalt,N,SoilFaceSlot(N,N3,N6),N5,N4)=trcSalt_MacpTranspFlxM_3D(idsalt,N,SoilFaceSlot(N,N3,N6),N5,N4)+trcSalt_MacPore_Diffus(idsalt)
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SoluteDifsMacroporeM

!----------------------------------------------------------------------
  subroutine SoluteAdvDifsMicMacporeM(M,N,N1,N2,N3,N4,N5,N6,THETW1)
  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  REAL(R8), INTENT(IN):: THETW1(JZ,JY,JX)
  integer :: idsalt


  !     SOLUTE TRANSPORT IN MICROPORES
  !
  call SoluteAdvMicroporeM(M,N,N1,N2,N3,N4,N5,N6)
  !
  call SoluteDifsMicroporeM(M,N,N1,N2,N3,N4,N5,N6,THETW1)

  !     SOLUTE TRANSPORT IN MACROPORES
  call SoluteAdvMacroporeM(M,N,N1,N2,N3,N4,N5,N6)
  !
  call SoluteDifsMacroporeM(M,N,N1,N2,N3,N4,N5,N6)

  end subroutine SoluteAdvDifsMicMacporeM
!----------------------------------------------------------------------
  subroutine SaltAdvDifuXMicMacporesM(M,N4,N5,N6)
  implicit none
  integer, intent(in) :: M,N4,N5,N6

  integer :: idsalt
  real(r8) :: VFLW
  real(r8) :: trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)
  real(r8) :: VOLWHS,VOLWT
  real(r8) :: trcSalt_Difus_Mac2Micp_flxM(idsalt_beg:idsaltb_end)

  !
  !     MACROPORE TO MICROPORE TRANSFER advection
  !
  IF(FWatExMacP2MicPM_vr(M,N6,N5,N4).GT.0.0)THEN
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
!
!     MICROPORE TO MACROPORE TRANSFER advection
!
  ELSEIF(FWatExMacP2MicPM_vr(M,N6,N5,N4).LT.0.0)THEN
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Adv2MicP_flx(idsalt)=VFLW*AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
    !
    !     NO MACROPORE TO MICROPORE TRANSFER
    !
  ELSE
    trcSalt_Adv2MicP_flx(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flxM_vr(idsalt,N6,N5,N4)=trcSalt_Adv2MicP_flx(idsalt)
  ENDDO

  !diffusion
  !
  IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
    VOLWHS = AMIN1(XFRS*VGeomLayer_vr(N6,N5,N4),VLWatMacPM_vr(M,N6,N5,N4))
    VOLWT  = VLWatMicPM_vr(M,N6,N5,N4)+VOLWHS

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*VOLWHS)/VOLWT
    ENDDO

    DO idsalt=idsalt_psoil_beg,idsalt_psoil_end
      trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    ENDDO

    DO idsalt=idsalt_pband_beg,idsalt_pband_end
      trcSalt_Difus_Mac2Micp_flxM(idsalt)=dts_HeatWatTP*(AZMAX1(trcSalt_soHml2_vr(idsalt,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcSalt_solml2_vr(idsalt,N6,N5,N4))*VOLWHS)/VOLWT*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    ENDDO
  ELSE
    trcSalt_Difus_Mac2Micp_flxM(idsalt_beg:idsaltb_end)=0.0_r8
  ENDIF

  DO idsalt=idsalt_beg,idsaltb_end
    trcSalt_Mac2MicPore_flxM_vr(idsalt,N6,N5,N4)=trcSalt_Mac2MicPore_flxM_vr(idsalt,N6,N5,N4)+trcSalt_Difus_Mac2Micp_flxM(idsalt)
  ENDDO

  ! Hourly exchange is accumulated after acceptance in AccumSaltFluxesM.

  end subroutine SaltAdvDifuXMicMacporesM
!----------------------------------------------------------------------
  subroutine SaltTranspXInterGridsM(M,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='SaltTranspXInterGridsM'
  integer :: L,N1,N2,N3,N4,N5,N6
  integer :: LL,N,NY,NX
  real(r8) :: VOLWHS,VOLWT
  real(r8) :: VFLW,VLWPA1,VLWPB1,VLWPA2,VLWPB2
  real(r8) :: THETW1(JZ,JY,JX)

!     begin_execution
!
  call PrintInfo('beg '//subname)
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      N1=NX;N2=NY
      DO  L=1,NL_col(NY,NX)
        !
        N3=L
        IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(NY,NX))call SaltAdvDifuXMicMacporesM(M,N1,N2,N3)

        DO  N=FlowDirIndicator_col(N2,N1),3

          IF(N.EQ.iWestEastDirection)THEN            
            IF(NX.EQ.NHE)THEN          !skip eatern boundary
              cycle
            ELSE
              N4 = NX+1;N5 = NY;N6 = L
            ENDIF
          ELSEIF(N.EQ.iNorthSouthDirection)THEN            
            IF(NY.EQ.NVS)THEN          !skip southern boundaries
              cycle
            ELSE
              N4 = NX;N5 = NY+1;N6 = L
            ENDIF
          ELSEIF(N.EQ.iVerticalDirection)THEN
            IF(L.EQ.NL_col(NY,NX))THEN    !skip bottom boundary
              cycle
            ELSE
              N4 = NX;N5 = NY;N6 = L+1
            ENDIF
          ENDIF

          DO LL=N6,NL_col(N5,N4)
            IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
              N6=LL
              exit
            ENDIF
          enddo
          ! Internal transport requires an active destination layer.
          if(LL>NL_col(N5,N4))cycle
          !
          !     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS FROM
          !     WATER CONTENTS AND WATER FLUXES 'FLQM' FROM 'WATSUB'
          !
          !
          IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(NY,NX))THEN
            IF(N3.GE.NUM_col(N2,N1) .AND. N6.GE.NUM_col(N5,N4) .AND. N3.LE.NL_col(N2,N1) .AND. N6.LE.NL_col(N5,N4))THEN
              THETW1(N3,N2,N1)=AZMAX1(VLWatMicPM_vr(M,N3,N2,N1)/VLSoilMicP_vr(N3,N2,N1))
              THETW1(N6,N5,N4)=AZMAX1(VLWatMicPM_vr(M,N6,N5,N4)/VLSoilMicP_vr(N6,N5,N4))
    
              call SoluteAdvDifsMicMacporeM(M,N,N1,N2,N3,N4,N5,N6,THETW1)        
            ELSEIF(N.NE.iVerticalDirection)THEN
              trcSalt_MicpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,SoilFaceSlot(N,N3,N6),N5,N4)=0.0_r8
              trcSalt_MacpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,SoilFaceSlot(N,N3,N6),N5,N4)=0.0_r8
            ENDIF
          ELSE

            trcSalt_MicpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,SoilFaceSlot(N,N3,N6),N5,N4)=0.0_r8
            trcSalt_MacpTranspFlxM_3D(idsalt_beg:idsaltb_end,N,SoilFaceSlot(N,N3,N6),N5,N4)=0.0_r8
          ENDIF
        enddo
      enddo
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SaltTranspXInterGridsM

!------------------------------------------------------------------------------------------
  subroutine SaltFall2Snowpack(I,NY,NX)
  implicit none
  integer, intent(in) :: I,NY,NX

  integer :: idsalt

  if(salt_model)then
    DO idsalt=idsalt_beg,idsalt_end
      trcSalt_AquaAdv_flxM_snvr(idsalt,1,NY,NX) = (Rain2SnowSurf_col(NY,NX)*trcsalt_rain_mole_conc_col(idsalt,NY,NX) &
        +Irrig2SnowSurf_col(NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX))*dts_HeatWatTP
    ENDDO
  endif

  end subroutine SaltFall2Snowpack
end module IngridTranspMod
