module ErosionBalMod
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use EcoSIMConfig,     only: nlbiomcp=>NumLiveMicrbCompts,NumMicbAFunGrupsPerCmplx,NumMicbHFunGrupsPerCmplx
  use EcoSiMParDataMod, only: micpar
  use SoilPropertyDataType
  use RootDataType
  use DebugToolMod
  USE EcoSIMCtrlDataType
  use EcoSIMCtrlMod
  use MicrobialDataType
  USE SOMDataType
  use SedimentDataType
  use GridDataType
  USE AqueChemDatatype
  use FlagDataType
  use FertilizerDataType
  use RedistDataMod
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: SinkSediments
  public :: ZeroErosionArray
  contains

!------------------------------------------------------------------------------------------

  subroutine SinkSediments(NY,NX)
  
  !Sedimentation
  implicit none
  integer, intent(in) :: NY,NX
  character(len=*), parameter :: subname = 'SinkSediments'
  integer :: L,LL,M,N,K,NGL,MID,idom,NE
  real(r8) :: FSINK,FSAN,FSIL,FCLA,FCEC,FAEC
  real(r8) :: FNX
  real(r8) :: FOMC,FOMN,FOMP
  real(r8) :: FORC,FORN,FORP,FOHC,FOHN,FOHP,FOHA
  real(r8) :: FOSC,FOSA,FOSN,FOSP
  integer  :: NTF,NTP,NTX
  real(r8) :: FPX,FSNX
! begin_execution
! SINK ALL SOLID C,N,P IN POND
!
! BKDS=bulk density
! DLYR=layer depth
! VLS=hourly sinking rate from hour1.f
! FSINK=hourly rate for sediment sinking
!
  call PrintInfo('beg '//subname)
  D9885: DO L=NL_col(NY,NX)-1,1,-1     !should this include litter layer? JYT, 2026 Mar 6th.
    IF(SoilBulkDensity_vr(L,NY,NX).LE.ZERO .AND. DLYR_3D(3,L,NY,NX).GT.ZERO)THEN
      !sinking from water to sediment layer
      D9880: DO LL=L+1,NL_col(NY,NX)
        IF(DLYR_3D(3,LL,NY,NX).GT.ZEROS(NY,NX))exit
      ENDDO D9880

      FSINK=AMIN1(1.0_r8,VLS_col(NY,NX)/DLYR_3D(3,L,NY,NX))
!
!     SOIL MINERALS
!
!     *ER=sediment flux from erosion.f
!     sediment code:XSED=total,XSAN=sand,XSIL=silt,XCLA=clay
!
      
      FSAN           = FSINK*SAND_vr(L,NY,NX)
      FSIL           = FSINK*SILT_vr(L,NY,NX)
      FCLA           = FSINK*CLAY_vr(L,NY,NX)
      SAND_vr(L,NY,NX)  = SAND_vr(L,NY,NX)-FSAN
      SILT_vr(L,NY,NX)  = SILT_vr(L,NY,NX)-FSIL
      CLAY_vr(L,NY,NX)  = CLAY_vr(L,NY,NX)-FCLA

      SAND_vr(LL,NY,NX) = SAND_vr(LL,NY,NX)+FSAN
      SILT_vr(LL,NY,NX) = SILT_vr(LL,NY,NX)+FSIL
      CLAY_vr(LL,NY,NX) = CLAY_vr(LL,NY,NX)+FCLA
      !
      !     FERTILIZER POOLS
      !
      !     *ER=sediment flux from erosion.f
      !     sediment code:NH4,NH3,NHU,NO3=NH4,NH3,urea,NO3
      !
      DO NTF=ifertn_beg,ifertn_end
        FNX                         = FSINK*FertN_mole_soil_vr(NTF,L,NY,NX)
        FertN_mole_soil_vr(NTF,L,NY,NX)  = FertN_mole_soil_vr(NTF,L,NY,NX)-FNX
        FertN_mole_soil_vr(NTF,LL,NY,NX) = FertN_mole_soil_vr(NTF,LL,NY,NX)+FNX
      ENDDO
!
!     EXCHANGEABLE CATIONS AND ANIONS
!
!     sediment code
!     :XN4=adsorbed NH4
!     :XHY,XAL,XFE,XCA,XMG,XNA,XKA,XHC,AL2,FE2
!      =adsorbed H,Al,Fe,Ca,Mg,Na,K,HCO3,AlOH2,FeOH2
!     :XOH0,XOH1,XOH2=adsorbed R-,R-OH,R-OH2
!     :XH1P,XH2P=adsorbed HPO4,H2PO4
!
! exclude banded nutrients
      DO NTX=idx_CEC,idx_cation_soil_end
        FSNX                        = FSINK*trcx_solml_vr(NTX,L,NY,NX)
        trcx_solml_vr(NTX,L,NY,NX)  = trcx_solml_vr(NTX,L,NY,NX)-FSNX
        trcx_solml_vr(NTX,LL,NY,NX) = trcx_solml_vr(NTX,LL,NY,NX)+FSNX
      ENDDO

      DO NTX=idx_AEC,idx_anion_soil_end
        FSNX                        = FSINK*trcx_solml_vr(NTX,L,NY,NX)
        trcx_solml_vr(NTX,L,NY,NX)  = trcx_solml_vr(NTX,L,NY,NX)-FSNX
        trcx_solml_vr(NTX,LL,NY,NX) = trcx_solml_vr(NTX,LL,NY,NX)+FSNX
      ENDDO

!
!     PRECIPITATES
!
!     sediment code
!     :PALO,PFEO=precip AlOH,FeOH
!     :PCAC,PCAS=precip CaCO3,CaSO4
!     :PALP,PFEP=precip AlPO4,FEPO4
!     :PCPM,PCPD,PCPH=precip CaH4P2O8,CaHPO4,apatite
!
! only for non-banded precipitates
      DO NTP=idsp_beg,idsp_beg_band-1
        FPX                           = FSINK*trcp_saltpml_vr(NTP,L,NY,NX)
        trcp_saltpml_vr(NTP,L,NY,NX)  = trcp_saltpml_vr(NTP,L,NY,NX)-FPX
        trcp_saltpml_vr(NTP,LL,NY,NX) = trcp_saltpml_vr(NTP,LL,NY,NX)+FPX
      ENDDO

!     MICROBIAL C,N,P
!
      D1970: DO K=1,micpar%NumOfLitrCmplxs

!         OMC,OMN,OMP=microbial C,N,P
!         ORC,ORN,ORP=microbial residue C,N,P
!         OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!         OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!
        D1960: DO N=1,NumMicbHFunGrupsPerCmplx
          DO NGL=JGniH(N),JGnfH(N)
            DO M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              DO NE=1,NumPlantChemElms
                FOMC                              = FSINK*mBiomeHeter_vr(NE,MID,K,L,NY,NX)
                mBiomeHeter_vr(NE,MID,K,LL,NY,NX) = mBiomeHeter_vr(NE,MID,K,LL,NY,NX)+FOMC
                mBiomeHeter_vr(NE,MID,K,L,NY,NX)  = mBiomeHeter_vr(NE,MID,K,L,NY,NX)-FOMC
              ENDDO
            enddo
          enddo
        ENDDO D1960
      ENDDO D1970

!         OMC,OMN,OMP=microbial C,N,P
!         ORC,ORN,ORP=microbial residue C,N,P
!         OHC,OHN,OHP,OHA=adsorbed C,N,P,acetate
!         OSC,OAA,OSN,OSP=SOC,colonized SOC,SON,SOP
!
      DO N=1,NumMicbAFunGrupsPerCmplx
        DO NGL=JGniA(N),JGnfA(N)
          DO M=1,nlbiomcp
            MID=micpar%get_micb_id(M,NGL)        
            DO NE=1,NumPlantChemElms              
              FOMC                            = FSINK*mBiomeAutor_vr(NE,MID,L,NY,NX)
              mBiomeAutor_vr(NE,MID,LL,NY,NX) = mBiomeAutor_vr(NE,MID,LL,NY,NX)+FOMC
              mBiomeAutor_vr(NE,MID,L,NY,NX)  = mBiomeAutor_vr(NE,MID,L,NY,NX)-FOMC
            ENDDO
          enddo
        enddo
      ENDDO

!
!     MICROBIAL RESIDUE C,N,P
!
      D1900: DO K=1,micpar%NumOfLitrCmplxs
        D1940: DO M=1,ndbiomcp
          DO NE=1,NumPlantChemElms                   
            FORC                           = FSINK*OMBioResdu_vr(NE,M,K,L,NY,NX)
            OMBioResdu_vr(NE,M,K,LL,NY,NX) = OMBioResdu_vr(NE,M,K,LL,NY,NX)+FORC
            OMBioResdu_vr(NE,M,K,L,NY,NX)  = OMBioResdu_vr(NE,M,K,L,NY,NX)-FORC
          ENDDO
        ENDDO D1940
!
!       ADSORBED C,N,P
!
        DO idom=idom_beg,idom_end
          FOHC                         = FSINK*SorbedOM_vr(idom,K,L,NY,NX)
          SorbedOM_vr(idom,K,LL,NY,NX) = SorbedOM_vr(idom,K,LL,NY,NX)+FOHC
          SorbedOM_vr(idom,K,L,NY,NX)  = SorbedOM_vr(idom,K,L,NY,NX)-FOHC
        ENDDO
!
!       SOC,N,P
!
        D1930: DO M=1,jsken
          FOSA=FSINK*SolidOMAct_vr(M,K,L,NY,NX)
          SolidOMAct_vr(M,K,LL,NY,NX)=SolidOMAct_vr(M,K,LL,NY,NX)+FOSA
          SolidOMAct_vr(M,K,L,NY,NX)=SolidOMAct_vr(M,K,L,NY,NX)-FOSA
          DO NE=1,NumPlantChemElms
            FOSC                        = FSINK*SolidOM_vr(NE,M,K,L,NY,NX)
            SolidOM_vr(NE,M,K,LL,NY,NX) = SolidOM_vr(NE,M,K,LL,NY,NX)+FOSC
            SolidOM_vr(NE,M,K,L,NY,NX)  = SolidOM_vr(NE,M,K,L,NY,NX)-FOSC
          ENDDO
        ENDDO D1930
      ENDDO D1900
    ENDIF
  ENDDO D9885

  call PrintInfo('end '//subname)
  end subroutine SinkSediments
!------------------------------------------------------------------------------------------

  subroutine ZeroErosionArray(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX

  IF(iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
    tErosionSedmLoss_col(NY,NX) = 0.0_r8
    TSandErosed_col(NY,NX)        = 0.0_r8
    TSiltErosed_col(NY,NX)        = 0.0_r8
    TCLAYErosed_col(NY,NX)        = 0.0_r8
    TNH4Erosed_molN_col(NY,NX)         = 0.0_r8
    TNH3Erosed_molN_col(NY,NX)         = 0.0_r8
    TNUreaErosed_molN_col(NY,NX)       = 0.0_r8
    TNO3Erosed_molN_col(NY,NX)         = 0.0_r8
    TNH4ErosBand_molN_col(NY,NX)     = 0.0_r8
    TNH3ErosBand_molN_col(NY,NX)     = 0.0_r8
    TNUreaErosBand_molN_col(NY,NX)   = 0.0_r8
    TNO3ErosBand_molN_col(NY,NX)     = 0.0_r8
    TPO4Erosed_molP_col(NY,NX)      = 0._r8
    TPO4ErosBand_molP_col(NY,NX)    = 0._r8

    trcx_TER_col(idx_beg:idx_end,NY,NX)   = 0.0_r8
    trcp_TER_col(idsp_beg:idsp_end,NY,NX) = 0.0_r8

    TOMEERhetr_col(1:NumPlantChemElms,1:NumLiveHeterBioms,1:jcplx,NY,NX) = 0.0_r8
    TOMEERauto_col(1:NumPlantChemElms,1:NumLiveAutoBioms,NY,NX)          = 0.0_r8

    TORMER_col(1:NumPlantChemElms,1:ndbiomcp,1:jcplx,NY,NX) = 0.0_r8
    TOHMER_col(idom_beg:idom_end,1:jcplx,NY,NX)             = 0.0_r8
    TOSMER_col(1:NumPlantChemElms,1:jsken,1:jcplx,NY,NX)    = 0.0_r8
    TOSAER_col(1:jsken,1:jcplx,NY,NX)                       = 0.0_r8
  ENDIF
  end subroutine ZeroErosionArray

end module ErosionBalMod
