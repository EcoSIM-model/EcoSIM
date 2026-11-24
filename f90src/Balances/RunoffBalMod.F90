module RunoffBalMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use minimathmod, only : AZMAX1
  use EcoSiMParDataMod    , only : micpar  
  use EcoSIMConfig, only : jcplx => jcplxc  
  use EcosimConst, only : patomw,natomw
  use DebugToolMod
  USE EcoSimSumDataType
  USE ChemTranspDataType
  use GridConsts
  use SnowDataType
  use FlagDataType
  use SoilBGCDataType
  USE SedimentDataType
  use MicrobialDataType
  USE SoilWaterDataType
  USE AqueChemDatatype
  USE EcoSIMCtrlDataType
  USE GridDataType
  use EcoSIMCtrlMod
  use ElmIDMod   
  use SoilHeatDatatype
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: XGridBoundSolutesRunoff

  contains

  subroutine XGridBoundSolutesRunoff(I,J,NY,NX,NHW,NHE,NVN,NVS)
  !
  !Diagnose solute fluxes across the landscape boundariers 
  !N-W-S-E, and bottom
  implicit none
  integer, intent(in) :: I,J, NY,NX,NHW,NHE,NVN,NVS
  integer :: N
  integer :: N1,N2     !source grid ()
  integer :: N4,N5,N6  !dest grid
  integer :: L,NN
  real(r8) :: XN
  real(r8) :: CXR,ZXR,PXR
  real(r8) :: ZGR

  CXR=0._r8;ZXR=0._r8;PXR=0._r8;ZGR=0._r8

  D9985: DO L=NU_col(NY,NX),NL_col(NY,NX)
!
!     LOCATE EXTERNAL BOUNDARIES
!
!     N2,N1=NY,NX of source grid cell
!     N5,N4=NY,NX of destination grid cell
!
!flow from west to east, north to south, up to down
    N1=NX;N2=NY
    D9980: DO N=FlowDirIndicator_col(NY,NX),3
      D9975: DO NN=1,2

        IF(N.EQ.iWestEastDirection)THEN !east-west direction          
          IF(NN.EQ.iFront)THEN
            IF(NX.EQ.NHE)THEN !eastern boundary              
              N4 = NX+1;N5 = NY;N6 = L
              XN = -1.0_r8   !going out of eastern boundary
            ELSE
              cycle
            ENDIF
          ELSEIF(NN.EQ.iBehind)THEN
            IF(NX.EQ.NHW)THEN  !western boundary              
              N4 = NX;N5 = NY;N6 = L
              XN = 1.0_r8       !coming in from western boundary
            ELSE
              cycle
            ENDIF
          ENDIF
        ELSEIF(N.EQ.iNorthSouthDirection)THEN  !north-south direction          
          IF(NN.EQ.iFront)THEN
            IF(NY.EQ.NVS)THEN  !south boundary              
              N4 = NX;N5 = NY+1;N6 = L
              XN = -1.0_r8      !going out of southern boundary
            ELSE
              cycle
            ENDIF
          ELSEIF(NN.EQ.iBehind)THEN
            IF(NY.EQ.NVN)THEN  !north boundary              
              N4 = NX;N5 = NY;N6 = L
              XN = 1.0_r8       !coming in from northern boundary
            ELSE
              cycle
            ENDIF
          ENDIF
        ELSEIF(N.EQ.iVerticalDirection)THEN !vertical direction          
          IF(NN.EQ.iFront)THEN  
            IF(L.EQ.NL_col(NY,NX))THEN  !at the bottom
              N4 = NX;N5 = NY;N6 = L+1
              XN = -1.0_r8       !going out from layer L into L+1
            ELSE
              cycle
            ENDIF
          ELSEIF(NN.EQ.iBehind)THEN
            cycle
          ENDIF
        ENDIF
!
        call XBoundarySurfRunoffs(I,J,L,N,NY,NX,N1,N2,N4,N5,NN,XN,CXR,ZXR,PXR,ZGR)
    !
        call XBoundarySubSurfRunoffs(I,J,N,NY,NX,N1,N2,N4,N5,N6,XN)

    !     WATER, HEAT, SOLUTES IN SNOW DRIFT
        IF(N.NE.3 .AND. L.EQ.NU_col(NY,NX)) &
          call WaterHeatSoluteBySnowDrift(N,NN,N4,N5,NY,NX,CXR,ZXR,PXR,ZGR,XN)
        
      ENDDO D9975
!
    ENDDO D9980
  ENDDO D9985

  end subroutine XGridBoundSolutesRunoff

!------------------------------------------------------------------------------------------

  subroutine XBoundarySurfRunoffs(I,J,L,N,NY,NX,N1,N2,N4,N5,NN,XN,CXR,ZXR,PXR,ZGR)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: L   !vertical layer
  integer, intent(in) :: N   !horizontal direction
  integer, intent(in) :: NY,NX  !geographic location
  integer, intent(in) :: N1,N2,N4,N5,NN
  real(r8), intent(in) :: XN  !-1 outgoing, 1 incoming
  real(r8), intent(out) :: CXR,ZXR,PXR,ZGR

  real(r8) :: SS1,SS2,SS3,SS4,SSR
  real(r8) :: OMRof(idom_beg:idom_end)
  real(r8) :: MXE(NumPlantChemElms),MOE(NumPlantChemElms)
  real(r8) :: ECHY,ECOH,ECAL,ECFE,ECCA,ECMG,ECNA,ECKA
  real(r8) :: ECCO,ECHC,ECSO,ECCL,ECNO
  real(r8) :: ECNDQ
  real(r8) :: SEF,SEX,SEP,SET
  real(r8) :: PPE
  real(r8) :: ZPE
  real(r8) :: WX,HGR,PSS
  real(r8) :: WQRN,HQRN
  real(r8) :: OXR,ER
  integer :: K,M,NO,NGL,NE,MID,idom

!     begin_execution
!     RUNOFF BOUNDARY FLUXES OF WATER AND HEAT
!
!     QR,QS,WatBySnowRedistrib,IceBySnowRedistrib=runoff from surface water, 
!     snowpack snow,water,ice from watsub.f
!     CRUN,Qrunoff_CumYr_col=cumulative water and snow runoff
!     HeatOut_lnds=cumulative heat loss through lateral and lower boundaries

! surface runoff
  IF(N.NE.iVerticalDirection .AND. L.EQ.NU_col(NY,NX))THEN
    !horizontal direction and surface layer
    WQRN                 = XN*XGridSurfRunoff_2DH(N,NN,N5,N4)    
!    QRunSurf_col(N2,N1) = QRunSurf_col(N2,N1)+WQRN
    
    IF(ABS(WQRN).GT.ZEROS(N5,N4))THEN
      CRUN                     = CRUN-WQRN
      Qrunoff_CumYr_col(NY,NX) = Qrunoff_CumYr_col(NY,NX)-WQRN
      HQRN                     = XN*HeatXGridBySurfRunoff_2DH(N,NN,N5,N4)
      HeatOut_lnds             = HeatOut_lnds-HQRN
!      HeatRunSurf_col(N2,N1)   = HeatRunSurf_col(N2,N1)+HQRN
!
!     RUNOFF BOUNDARY FLUXES OF C, N AND P
!
!     XN=direction indicator
!
      CXR = XN*(trcg_FloXSurRunoff_2D(idg_CO2,N,NN,N5,N4)+trcg_FloXSurRunoff_2D(idg_CH4,N,NN,N5,N4))
      ZXR = XN*(trcn_FloXSurRunoff_2D(ids_NH4,N,NN,N5,N4)+trcg_FloXSurRunoff_2D(idg_NH3,N,NN,N5,N4) &
        +trcn_FloXSurRunoff_2D(ids_NO3,N,NN,N5,N4)+trcn_FloXSurRunoff_2D(ids_NO2,N,NN,N5,N4))
      ZGR      = XN*(trcg_FloXSurRunoff_2D(idg_N2O,N,NN,N5,N4)+trcg_FloXSurRunoff_2D(idg_N2,N,NN,N5,N4))
      PXR      = XN*(trcn_FloXSurRunoff_2D(ids_H2PO4,N,NN,N5,N4)+trcn_FloXSurRunoff_2D(ids_H1PO4,N,NN,N5,N4))
      OMRof(:) = 0.0_r8
      D2575: DO K = 1, jcplx
        DO idom=idom_beg,idom_end
          OMRof(idom)=OMRof(idom)+XN*DOM_FloXSurRunoff_2DH(idom,K,N,NN,N5,N4)
        ENDDO
      ENDDO D2575
      TOMOU_lnds(ielmc)               = TOMOU_lnds(ielmc)-CXR-OMRof(ielmc)-OMRof(idom_acetate)
      TOMOU_lnds(ielmn)               = TOMOU_lnds(ielmn)-ZXR-OMRof(ielmn)-ZGR
      TOMOU_lnds(ielmp)               = TOMOU_lnds(ielmp)-PXR-OMRof(ielmp)
      HydroSufDOCFlx_col(NY,NX)       = -OMRof(ielmc)-OMRof(idom_acetate)
      HydroSufDICFlx_col(NY,NX)       = -CXR
      HydroSufDONFlx_col(NY,NX) = HydroSufDONFlx_col(NY,NX)-OMRof(ielmn)
      HydroSufDINFlx_CumYr_col(NY,NX) = HydroSufDINFlx_CumYr_col(NY,NX)-ZXR-ZGR
      HydroSufDOPFlx_col(NY,NX) = HydroSufDOPFlx_col(NY,NX)-OMRof(ielmp)
      HydroSufDIPFlx_CumYr_col(NY,NX) = HydroSufDIPFlx_CumYr_col(NY,NX)-PXR
      OXR                             = XN*trcg_FloXSurRunoff_2D(idg_O2,N,NN,N5,N4)
      OXYGOU                          = OXYGOU-OXR
      HGR                             = XN*trcg_FloXSurRunoff_2D(idg_H2,N,NN,N5,N4)
      H2GOU                           = H2GOU+HGR
!
!     RUNOFF BOUNDARY FLUXES OF SOLUTES
!
!     XN=direction indicator
!     TOMOU_lnds(ielmp),TIONOU=total P,salt loss through lateral and lower boundaries
!
      IF(salt_model)THEN
        PSS=XN*patomw*(trcSalt_FloXSurRunoff_2D(idsalt_H0PO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_CaPO4,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_FeHPO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_CaHPO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_MgHPO4,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_H3PO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_FeH2PO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_CaH4P2O8,N,NN,N5,N4))
        SS1=XN*(trcSalt_FloXSurRunoff_2D(idsalt_Al,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_Fe,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_Hp,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_Ca,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_Mg,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_Na,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_K,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_OH,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_SO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_Cl,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_CO3,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_H0PO4,N,NN,N5,N4))
        SS2=XN*2.0_r8*(trcSalt_FloXSurRunoff_2D(idsalt_HCO3,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_AlOH,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_AlSO4,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_FeOH,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_FeSO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_CaOH,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_CaCO3,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_CaSO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_MgOH2,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_MgCO3,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_MgSO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_NaCO3,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_NaSO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_KSO4,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_CaPO4,N,NN,N5,N4))
        SS3=XN*3.0_r8*(trcSalt_FloXSurRunoff_2D(idsalt_AlOH2,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_FeOH2,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_CaHCO3,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_MgHCO3,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_FeHPO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_CaHPO4,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_MgHPO4,N,NN,N5,N4))
        SS4=XN*4.0_r8*(trcSalt_FloXSurRunoff_2D(idsalt_AlOH3,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_FeOH3,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_H3PO4,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_FeH2PO4,N,NN,N5,N4) &
          +trcSalt_FloXSurRunoff_2D(idsalt_CaH4P2O8,N,NN,N5,N4)) &
          +XN*5.0_r8*(trcSalt_FloXSurRunoff_2D(idsalt_AlOH4,N,NN,N5,N4)+trcSalt_FloXSurRunoff_2D(idsalt_FeOH4,N,NN,N5,N4))
        PSS=PSS+XN*patomw*trcSalt_FloXSnow_2DH(idsalt_H0PO4,N,NN,N5,N4)
        TOMOU_lnds(ielmp)=TOMOU_lnds(ielmp)-PSS
        SSR=SS1+SS2+SS3+SS4
        TIONOU=TIONOU-SSR
        HydroIonFlx_CumYr_col(NY,NX)=HydroIonFlx_CumYr_col(NY,NX)-SSR
!     WRITE(20,3336)'SSR',I,J,N,N5,N4,SSR,SS1,SS2,SS3,SS4,TIONOU
!3336  FORMAT(A8,5I6,20F16.9)
!
!       SURFACE RUNOFF ELECTRICAL CONDUCTIVITY
!
!       QR=surface runoff from watsub.f
!       XQR*=solute in runoff from TranspSalt.f
!       ECNDQ=electrical conductivity of runoff
!       where are those coefficients from?
        WX=XGridSurfRunoff_2DH(N,NN,N5,N4)
        IF(ABS(WX).GT.ZEROS(N5,N4))THEN
          ECHY  = 0.337_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_Hp,N,NN,N5,N4)/WX)
          ECOH  = 0.192_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_OH,N,NN,N5,N4)/WX)
          ECAL  = 0.056_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_Al,N,NN,N5,N4)*3.0_r8/WX)
          ECFE  = 0.051_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_Fe,N,NN,N5,N4)*3.0_r8/WX)
          ECCA  = 0.060_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_Ca,N,NN,N5,N4)*2.0_r8/WX)
          ECMG  = 0.053_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_Mg,N,NN,N5,N4)*2.0_r8/WX)
          ECNA  = 0.050_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_Na,N,NN,N5,N4)/WX)
          ECKA  = 0.070_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_K,N,NN,N5,N4)/WX)
          ECCO  = 0.072_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_CO3,N,NN,N5,N4)*2.0_r8/WX)
          ECHC  = 0.044_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_HCO3,N,NN,N5,N4)/WX)
          ECSO  = 0.080_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_SO4,N,NN,N5,N4)*2.0_r8/WX)
          ECCL  = 0.076_r8*AZMAX1(trcSalt_FloXSurRunoff_2D(idsalt_Cl,N,NN,N5,N4)/WX)
          ECNO  = 0.071_r8*AZMAX1(trcn_FloXSurRunoff_2D(ids_NO3,N,NN,N5,N4)/(WX*natomw))
          ECNDQ = ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA+ECCO+ECHC+ECSO+ECCL+ECNO
!     WRITE(*,9991)'ECNDQ',iYearCurrent,I,J,N4,N5,N,NN,WX,ECNDQ
!9991  FORMAT(A8,7I4,2E12.4)
        ELSE
          ECNDQ=0.0_r8
        ENDIF
      ENDIF
!
!     RUNOFF BOUNDARY FLUXES OF SEDIMENT FROM EROSION
!
!     iErosionMode=erosion flag
!     *ER=sediment flux from erosion.f
!     sediment code:XSED=total,XSAN=sand,XSIL=silt,XCLA=clay
!     TSedmErossLoss_lnds,SedmErossLoss_CumYr_col=cumulative sediment loss through lateral and lower boundaries
!
      IF(N.NE.3.AND.iErosionMode.EQ.ieros_frzthaweros.OR.iErosionMode.EQ.ieros_frzthawsomeros)THEN
        IF(ABS(cumSed_Eros_2D(N,NN,N5,N4)).GT.ZEROS(N5,N4))THEN
          ER                             = XN*cumSed_Eros_2D(N,NN,N5,N4)
          TSedmErossLoss_lnds            = TSedmErossLoss_lnds-ER
          SedmErossLoss_CumYr_col(NY,NX) = SedmErossLoss_CumYr_col(NY,NX)-ER
!
!
!         RUNOFF BOUNDARY FLUXES OF ORGANIC MATTER FROM EROSION

!         MICROBIAL C IN RUNOFF SEDIMENT
!
          MXE(ielmc)=0.0_r8
          MXE(ielmn)=XN*natomw*(trcx_Eros_mol_2D(idx_NH4,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_NH4B,N,NN,N5,N4))
          ZPE=XN*natomw*(XNH4Soil_Eros_molN_2D(N,NN,N5,N4)+XNH3Soil_Eros_molN_2D(N,NN,N5,N4) &
            +XUreaSoil_Eros_molN_2D(N,NN,N5,N4)+XNO3Soil_Eros_molN_2D(N,NN,N5,N4)+XNH4Band_Eros_molN_2D(N,NN,N5,N4) &
            +XNH3Band_Eros_molN_2D(N,NN,N5,N4)+XUreaBand_Eros_molN_2D(N,NN,N5,N4)+XNO3Band_Eros_molN_2D(N,NN,N5,N4))

          MXE(ielmp)=XN*patomw*(trcx_Eros_mol_2D(idx_HPO4,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_H2PO4,N,NN,N5,N4) &
            +trcx_Eros_mol_2D(idx_HPO4B,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_H2PO4B,N,NN,N5,N4))

          PPE=XN*patomw*(1._r8*(trcp_Eros_mol_2D(idsp_AlPO4,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_FePO4,N,NN,N5,N4) &
            +trcp_Eros_mol_2D(idsp_CaHPO4,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_AlPO4B,N,NN,N5,N4) &
            +trcp_Eros_mol_2D(idsp_FePO4B,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_CaHPO4B,N,NN,N5,N4)) &
            +2.0_r8*(trcp_Eros_mol_2D(idsp_CaH4P2O8,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_CaH4P2O8B,N,NN,N5,N4)) &
            +3.0_r8*(trcp_Eros_mol_2D(idsp_Apatite,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_ApatiteBand,N,NN,N5,N4))&
            +XPO4Soil_Eros_molP_2D(N,NN,N5,N4)+XPO4Band_Eros_molP_2D(N,NN,N5,N4))

          MOE(:)=0.0_r8
          
          D3580: DO K=1,jcplx
            DO NO=1,NumMicbHFunGrupsPerCmplx
              DO M=1,nlbiomcp
                DO NGL=JGniH(NO),JGnfH(NO)
                  MID=micpar%get_micb_id(M,NGL)                
                  DO NE=1,NumPlantChemElms
                    MOE(NE)=MOE(NE)+XN*OMEERhetr_2D(NE,MID,K,N,NN,N5,N4)
                  ENDDO
                enddo
              enddo
            enddo
          ENDDO D3580
          DO NO=1,NumMicbAFunGrupsPerCmplx
            DO M=1,nlbiomcp
              DO NGL=JGniA(NO),JGnfA(NO)
                MID=micpar%get_micb_id(M,NGL)
                DO NE=1,NumPlantChemElms
                  MOE(NE)=MOE(NE)+XN*OMEERauto_2D(NE,MID,N,NN,N5,N4)
                ENDDO
              enddo
            enddo
          enddo

!
    !     MICROBIAL RESIDUE C IN RUNOFF SEDIMENT
!
          D3575: DO K=1,jcplx
            D3570: DO M=1,ndbiomcp
              DO NE=1,NumPlantChemElms
                MOE(NE)=MOE(NE)+XN*OMBioResdu_Eros_2D(NE,M,K,N,NN,N5,N4)
              ENDDO
            ENDDO D3570
!
        !   DOC, ADSORBED AND HUMUS C IN RUNOFF SEDIMENT
!
            DO NE=1,NumPlantChemElms            
              MOE(NE)=MOE(NE)+XN*SorbedOM_Eros_2D(NE,K,N,NN,N5,N4)
            ENDDO
            MOE(ielmc)=MOE(ielmc)+XN*SorbedOM_Eros_2D(idom_acetate,K,N,NN,N5,N4)
            D3565: DO M=1,jsken
              DO NE=1,NumPlantChemElms            
                MOE(NE)=MOE(NE)+XN*SolidOM_Eros_2D(NE,M,K,N,NN,N5,N4)
              ENDDO
            ENDDO D3565
          ENDDO D3575

          DO NE=1,NumPlantChemElms            
            TOMOU_lnds(NE)=TOMOU_lnds(NE)-MOE(NE)-MXE(NE)
          ENDDO
          TOMOU_lnds(ielmn) = TOMOU_lnds(ielmn)-ZPE
          TOMOU_lnds(ielmp) = TOMOU_lnds(ielmp)-PPE

          HydroSufDOCFlx_col(NY,NX)       = HydroSufDOCFlx_col(NY,NX)-MOE(ielmc)
          HydroSufDONFlx_col(NY,NX) = HydroSufDONFlx_col(NY,NX)-MOE(ielmn)
          HydroSufDOPFlx_col(NY,NX) = HydroSufDOPFlx_col(NY,NX)-MOE(ielmp)
          HydroSufDICFlx_col(NY,NX)       = HydroSufDICFlx_col(NY,NX)-MXE(ielmc)
          HydroSufDINFlx_CumYr_col(NY,NX) = HydroSufDINFlx_CumYr_col(NY,NX)-MXE(ielmn)-ZPE
          HydroSufDIPFlx_CumYr_col(NY,NX) = HydroSufDIPFlx_CumYr_col(NY,NX)-MXE(ielmp)-PPE
!
!         ADSORBED AND PRECIPITATED SALTS IN RUNOFF SEDIMENTS

!         *ER=sediment flux from erosion.f
!         sediment code
!         TIONOU,HydroIonFlx_CumYr_col=total salt loss through lateral and lower boundaries
!
          IF(salt_model)THEN
            SEF=XN*(XNH3Soil_Eros_molN_2D(N,NN,N5,N4)+XUreaSoil_Eros_molN_2D(N,NN,N5,N4)+XNO3Soil_Eros_molN_2D(N,NN,N5,N4) &
              +XNH3Band_Eros_molN_2D(N,NN,N5,N4)+XUreaBand_Eros_molN_2D(N,NN,N5,N4)+XNO3Band_Eros_molN_2D(N,NN,N5,N4)) &
              +2.0*(XNH4Soil_Eros_molN_2D(N,NN,N5,N4)+XNH4Band_Eros_molN_2D(N,NN,N5,N4))
            SEX=XN*(trcx_Eros_mol_2D(idx_Hp,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_Al,N,NN,N5,N4) &
              +trcx_Eros_mol_2D(idx_Fe,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_Ca,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_Mg,N,NN,N5,N4) &
              +trcx_Eros_mol_2D(idx_Na,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_K,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_COOH,N,NN,N5,N4) &
              +trcx_Eros_mol_2D(idx_OHe,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_OHeB,N,NN,N5,N4)) &
              +XN*2.0*(trcx_Eros_mol_2D(idx_NH4,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_NH4B,N,NN,N5,N4) &
              +trcx_Eros_mol_2D(idx_OH,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_OHB,N,NN,N5,N4)) &
              +XN*3.0*(trcx_Eros_mol_2D(idx_AlOH2,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_FeOH2,N,NN,N5,N4) &
              +trcx_Eros_mol_2D(idx_OHp,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_OHpB,N,NN,N5,N4) &
              +trcx_Eros_mol_2D(idx_HPO4,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_HPO4B,N,NN,N5,N4)) &
              +XN*4.0*(trcx_Eros_mol_2D(idx_H2PO4,N,NN,N5,N4)+trcx_Eros_mol_2D(idx_H2PO4B,N,NN,N5,N4))
            SEP=XN*2.0*(trcp_Eros_mol_2D(idsp_CaCO3,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_CaSO4,N,NN,N5,N4) &
              +trcp_Eros_mol_2D(idsp_AlPO4,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_FePO4,N,NN,N5,N4) &
              +trcp_Eros_mol_2D(idsp_AlPO4B,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_FePO4B,N,NN,N5,N4)) &
              +XN*3.0*(trcp_Eros_mol_2D(idsp_CaHPO4,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_CaHPO4B,N,NN,N5,N4)) &
              +XN*4.0*(trcp_Eros_mol_2D(idsp_AlOH3,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_FeOH3,N,NN,N5,N4)) &
              +XN*7.0*(trcp_Eros_mol_2D(idsp_CaH4P2O8,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_CaH4P2O8B,N,NN,N5,N4)) &
              +XN*9.0*(trcp_Eros_mol_2D(idsp_Apatite,N,NN,N5,N4)+trcp_Eros_mol_2D(idsp_ApatiteBand,N,NN,N5,N4))
              +XN*(XPO4Soil_Eros_molP_2D(N,NN,N5,N4)+XPO4Band_Eros_molP_2D(N,NN,N5,N4))

            SET=SEF+SEX+SEP
            TIONOU=TIONOU-SET
            HydroIonFlx_CumYr_col(NY,NX)=HydroIonFlx_CumYr_col(NY,NX)-SET
!     WRITE(*,3342)'SET',I,J,N4,N5,NN,N,SET,SEF,SEX,SEP,TIONOU
!3342  FORMAT(A8,6I4,12F18.6)
          ENDIF
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  end subroutine XBoundarySurfRunoffs
!------------------------------------------------------------------------------------------

  subroutine XBoundarySubSurfRunoffs(I,J,N,NY,NX,N1,N2,N4,N5,N6,XN)
  !
  !Description
  !Subsurface across boundary flow
  implicit none
  integer, intent(in) :: I,J,N,NY,NX
  integer, intent(in) :: N1,N2,N4,N5,N6
  real(r8), intent(in) :: XN
  real(r8) :: ECHY,ECOH,ECAL,ECFE,ECCA,ECMG,ECNA,ECKA
  real(r8) :: ECCO,ECHC,ECSO,ECCL,ECNO
  real(r8) :: ECNDX
  real(r8) :: WO,SO,HO
  real(r8) :: MOD(idom_beg:idom_end),HOD,OOD
  real(r8) :: MXD(NumPlantChemElms),ZGD
  real(r8) :: SSD,SHD,PHD,PQD
  real(r8) :: WX
  integer :: K,idom
!     begin_execution
!     SUBSURFACE BOUNDARY FLUXES OF WATER AND HEAT
!
!     FLW,WaterFlowMacP,HFLW=micropore,macropore,heat flux through lateral and lower boundaries from watsub.f
!     QH2OLoss_lnds,HeatOut_lnds=cumulative water, heat loss through lateral and lower boundaries
!     H2OLoss_CumYr_col,QDischarg2WTBL_col=cumulative,hourly water loss through lateral and lower boundaries
!

  IF(FlowDirIndicator_col(NY,NX).NE.iVerticalDirection .OR. N.EQ.iVerticalDirection)THEN
    HO           = XN*HeatFlow2Soil_3D(N,N6,N5,N4)
    HeatOut_lnds = HeatOut_lnds-HO
    WO           = XN*(WaterFlowSoiMicP_3D(N,N6,N5,N4)+WaterFlowSoiMacP_3D(N,N6,N5,N4))   !XN<0, going out grid

    IF(abs(WO)>0._r8)THEN

      QH2OLoss_lnds            = QH2OLoss_lnds-WO
      H2OLoss_CumYr_col(N2,N1) = H2OLoss_CumYr_col(N2,N1)-WO
!
!     SUBSURFACE BOUNDARY FLUXES OF N2O, N2, NH4, NH3, NO3, NO2 AND DON

!!
      MOD=0.0_r8

      D450: DO K=1,jcplx
        DO idom=idom_beg,idom_end
          MOD(idom)=MOD(idom)+XN*(DOM_MicpTransp_3D(idom,K,N,N6,N5,N4)+DOM_Macp_Transp_flx_3D(idom,K,N,N6,N5,N4))
        ENDDO
      ENDDO D450

      MXD(ielmc)=XN*(trcs_TransptMicP_3D(idg_CO2,N,N6,N5,N4)+trcs_TransptMacP_3D(idg_CO2,N,N6,N5,N4) &
        +Gas_AdvDif_Flx_3D(idg_CO2,N,N6,N5,N4)+trcs_TransptMicP_3D(idg_CH4,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(idg_CH4,N,N6,N5,N4)+Gas_AdvDif_Flx_3D(idg_CH4,N,N6,N5,N4))

      MXD(ielmn)=XN*(trcs_TransptMicP_3D(ids_NH4,N,N6,N5,N4)+trcs_TransptMicP_3D(idg_NH3,N,N6,N5,N4) &
        +trcs_TransptMicP_3D(ids_NO3,N,N6,N5,N4) &
        +trcs_TransptMicP_3D(ids_NH4B,N,N6,N5,N4)+trcs_TransptMicP_3D(idg_NH3B,N,N6,N5,N4)&
        +trcs_TransptMicP_3D(ids_NO3B,N,N6,N5,N4) &
        +trcs_TransptMicP_3D(ids_NO2,N,N6,N5,N4)+trcs_TransptMicP_3D(ids_NO2B,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(ids_NH4,N,N6,N5,N4)+trcs_TransptMacP_3D(idg_NH3,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(ids_NO3,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(ids_NH4B,N,N6,N5,N4)+trcs_TransptMacP_3D(idg_NH3B,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(ids_NO3B,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(ids_NO2,N,N6,N5,N4)+trcs_TransptMacP_3D(ids_NO2B,N,N6,N5,N4))

      ZGD=XN*(trcs_TransptMicP_3D(idg_N2,N,N6,N5,N4)+Gas_AdvDif_Flx_3D(idg_N2,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(idg_N2,N,N6,N5,N4) &
        +trcs_TransptMicP_3D(idg_N2O,N,N6,N5,N4)+Gas_AdvDif_Flx_3D(idg_N2O,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(idg_N2O,N,N6,N5,N4) &
        +Gas_AdvDif_Flx_3D(idg_NH3,N,N6,N5,N4))

      MXD(ielmp)=XN*(trcs_TransptMicP_3D(ids_H2PO4,N,N6,N5,N4)+trcs_TransptMicP_3D(ids_H2PO4B,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(ids_H2PO4,N,N6,N5,N4)+trcs_TransptMacP_3D(ids_H2PO4B,N,N6,N5,N4)&
        +trcs_TransptMicP_3D(ids_H1PO4,N,N6,N5,N4) &
        +trcs_TransptMicP_3D(ids_H1PO4B,N,N6,N5,N4)+trcs_TransptMacP_3D(ids_H1PO4,N,N6,N5,N4) &
        +trcs_TransptMacP_3D(ids_H1PO4B,N,N6,N5,N4))

      TOMOU_lnds(ielmc)          = TOMOU_lnds(ielmc)-MOD(ielmc)-MOD(idom_acetate)-MXD(ielmc)
      TOMOU_lnds(ielmn)          = TOMOU_lnds(ielmn)-MOD(ielmn)-MXD(ielmn)-ZGD
      TOMOU_lnds(ielmp)          = TOMOU_lnds(ielmp)-MOD(ielmp)-MXD(ielmp)
      HydroSubsDOCFlx_col(N2,N1) = -MOD(ielmc)-MOD(idom_acetate)
      HydroSubsDONFlx_col(N2,N1) = -MOD(ielmn)
      HydroSubsDOPFlx_col(N2,N1) = -MOD(ielmp)
      HydroSubsDICFlx_col(N2,N1) = -MXD(ielmc)
      HydroSubsDINFlx_col(N2,N1) = -MXD(ielmn)
      HydroSubsDIPFlx_col(N2,N1) = -MXD(ielmp)
!
!     SUBSURFACE BOUNDARY FLUXES OF O2
!
!     OXYGOU,H2GOU=cumulative O2,H2 loss through lateral and lower boundaries
!
      OOD    = XN*(trcs_TransptMicP_3D(idg_O2,N,N6,N5,N4)+trcs_TransptMacP_3D(idg_O2,N,N6,N5,N4)+Gas_AdvDif_Flx_3D(idg_O2,N,N6,N5,N4))
      OXYGOU = OXYGOU-OOD
      HOD    = XN*(trcs_TransptMicP_3D(idg_H2,N,N6,N5,N4)+trcs_TransptMacP_3D(idg_H2,N,N6,N5,N4)+Gas_AdvDif_Flx_3D(idg_H2,N,N6,N5,N4))
      H2GOU  = H2GOU-HOD
!
!     SUBSURFACE BOUNDARY FLUXES OF SOLUTES

!     TIONOU,HydroIonFlx_CumYr_col=total salt loss through lateral and lower boundaries
!
      IF(salt_model)THEN
        PQD=XN*patomw*(trcSalt_TransptMicP_3D(idsalt_H0PO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_H0PO4B,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_CaPO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaPO4B,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeHPO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_CaHPO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_MgHPO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeHPO4B,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_CaHPO4B,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_MgHPO4B,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_H3PO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeH2PO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaH4P2O8,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_H3PO4B,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeH2PO4B,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaH4P2O8B,N,N6,N5,N4))

        PHD=XN*patomw*(trcSalt_TransptMacP_3D(idsalt_H0PO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_H0PO4B,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_CaPO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaPO4B,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeHPO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_CaHPO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_MgHPO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeHPO4B,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_CaHPO4B,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_MgHPO4B,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_H3PO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeH2PO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaH4P2O8,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_H3PO4B,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeH2PO4B,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaH4P2O8B,N,N6,N5,N4))

        TOMOU_lnds(ielmp)=TOMOU_lnds(ielmp)-PQD-PHD
        SSD=XN*(trcSalt_TransptMicP_3D(idsalt_Al,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_Fe,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_Hp,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_Ca,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_Mg,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_Na,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_K,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_OH,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_SO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_Cl,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CO3,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_H0PO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_H0PO4B,N,N6,N5,N4)+2.0*(trcSalt_TransptMicP_3D(idsalt_HCO3,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_AlOH,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_AlSO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_FeOH,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeSO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_CaOH,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaCO3,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_CaSO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_MgOH2,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_MgCO3,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_MgSO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_NaCO3,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_NaSO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_KSO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaPO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_CaPO4B,N,N6,N5,N4)) &
          +3.0*(trcSalt_TransptMicP_3D(idsalt_AlOH2,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeOH2,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaHCO3,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_MgHCO3,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeHPO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaHPO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_MgHPO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeHPO4B,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaHPO4B,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_MgHPO4B,N,N6,N5,N4)) &
          +4.0*(trcSalt_TransptMicP_3D(idsalt_AlOH3,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_FeOH3,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_H3PO4,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeH2PO4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaH4P2O8,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_H3PO4B,N,N6,N5,N4) &
          +trcSalt_TransptMicP_3D(idsalt_FeH2PO4B,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_CaH4P2O8B,N,N6,N5,N4)) &
          +5.0*(trcSalt_TransptMicP_3D(idsalt_AlOH4,N,N6,N5,N4)+trcSalt_TransptMicP_3D(idsalt_FeOH4,N,N6,N5,N4)))

        SHD=XN*(trcSalt_TransptMacP_3D(idsalt_Al,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_Fe,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_Hp,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_Ca,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_Mg,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_Na,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_K,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_OH,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_SO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_Cl,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CO3,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_H0PO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_H0PO4B,N,N6,N5,N4) &
          +2.0*(trcSalt_TransptMacP_3D(idsalt_HCO3,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_AlOH,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_AlSO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_FeOH,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeSO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_CaOH,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaCO3,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_CaSO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_MgOH2,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_MgCO3,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_MgSO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_NaCO3,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_NaSO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_KSO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaPO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_CaPO4B,N,N6,N5,N4)) &
          +3.0*(trcSalt_TransptMacP_3D(idsalt_AlOH2,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeOH2,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaHCO3,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_MgHCO3,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeHPO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaHPO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_MgHPO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeHPO4B,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaHPO4B,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_MgHPO4B,N,N6,N5,N4)) &
          +4.0*(trcSalt_TransptMacP_3D(idsalt_AlOH3,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_FeOH3,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_H3PO4,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeH2PO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaH4P2O8,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_H3PO4B,N,N6,N5,N4) &
          +trcSalt_TransptMacP_3D(idsalt_FeH2PO4B,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CaH4P2O8B,N,N6,N5,N4)) &
          +5.0*(trcSalt_TransptMacP_3D(idsalt_AlOH4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_AlOH4,N,N6,N5,N4)))

        SO=SSD+SHD
        TIONOU=TIONOU-SO
        HydroIonFlx_CumYr_col(N2,N1)=HydroIonFlx_CumYr_col(N2,N1)-SO

!
!     SUBSURFACE FLUX ELECTRICAL CONDUCTIVITY
!
!     ECNDQ=electrical conductivity of water flux
!
        WX=WaterFlowSoiMicP_3D(N,N6,N5,N4)+WaterFlowSoiMacP_3D(N,N6,N5,N4)
        IF(ABS(WX).GT.ZEROS(N2,N1))THEN
          ECHY  = 0.337*AZMAX1((trcSalt_TransptMicP_3D(idsalt_Hp,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_Hp,N,N6,N5,N4))/WX)
          ECOH  = 0.192*AZMAX1((trcSalt_TransptMicP_3D(idsalt_OH,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_OH,N,N6,N5,N4))/WX)
          ECAL  = 0.056*AZMAX1((trcSalt_TransptMicP_3D(idsalt_Al,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_Ca,N,N6,N5,N4))*3.0/WX)
          ECFE  = 0.051*AZMAX1((trcSalt_TransptMicP_3D(idsalt_Fe,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_Fe,N,N6,N5,N4))*3.0/WX)
          ECCA  = 0.060*AZMAX1((trcSalt_TransptMicP_3D(idsalt_Ca,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_Ca,N,N6,N5,N4))*2.0/WX)
          ECMG  = 0.053*AZMAX1((trcSalt_TransptMicP_3D(idsalt_Mg,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_Mg,N,N6,N5,N4))*2.0/WX)
          ECNA  = 0.050*AZMAX1((trcSalt_TransptMicP_3D(idsalt_Na,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_Na,N,N6,N5,N4))/WX)
          ECKA  = 0.070*AZMAX1((trcSalt_TransptMicP_3D(idsalt_K,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_K,N,N6,N5,N4))/WX)
          ECCO  = 0.072*AZMAX1((trcSalt_TransptMicP_3D(idsalt_CO3,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_CO3,N,N6,N5,N4))*2.0/WX)
          ECHC  = 0.044*AZMAX1((trcSalt_TransptMicP_3D(idsalt_HCO3,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_HCO3,N,N6,N5,N4))/WX)
          ECSO  = 0.080*AZMAX1((trcSalt_TransptMicP_3D(idsalt_SO4,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_SO4,N,N6,N5,N4))*2.0/WX)
          ECCL  = 0.076*AZMAX1((trcSalt_TransptMicP_3D(idsalt_Cl,N,N6,N5,N4)+trcSalt_TransptMacP_3D(idsalt_Cl,N,N6,N5,N4))/WX)
          ECNO  = 0.071*AZMAX1((trcs_TransptMicP_3D(ids_NO3,N,N6,N5,N4)+trcs_TransptMacP_3D(ids_NO3,N,N6,N5,N4))/(WX*natomw))
          ECNDX = ECHY+ECOH+ECAL+ECFE+ECCA+ECMG+ECNA+ECKA+ECCO+ECHC+ECSO+ECCL+ECNO
        ELSE
          ECNDX=0.0_r8
        ENDIF
      ENDIF
!      SG=SG+trcs_TransptMicP_3D(idg_H2,N,N6,N5,N4)+Gas_AdvDif_Flx_3D(idg_H2,N,N6,N5,N4)
    ENDIF
  ENDIF
  end subroutine XBoundarySubSurfRunoffs

!------------------------------------------------------------------------------------------

  subroutine WaterHeatSoluteBySnowDrift(N,NN,N4,N5,NY,NX,CXR,ZXR,PXR,ZGR,XN)
  implicit none
  integer, intent(in) :: N,NY,NX,N4,N5
  REAL(R8), INTENT(IN) :: CXR,ZXR,PXR,ZGR
  real(r8), intent(in) :: XN
  integer, intent(in) :: NN
  character(len=*), parameter :: subname='WaterHeatSoluteBySnowDrift'
  real(r8) :: CXS,ZXS,ZGS,PXS
  real(r8) :: SS1,SS2,SS3,SS4,SSR
  real(r8) :: WQRS,HQRS,OXS,PSS
!     begin_execution
  call PrintInfo('beg '//subname)
  !check for it is at surface for lateral flow
  
  WQRS=XN*(DrySnoByRedistrib_2DH(N,NN,N5,N4)+WatSnoByRedist_2DH(N,NN,N5,N4)+IceSnoBySnowRedist_2DH(N,NN,N5,N4))
  IF(ABS(WQRS).GT.ZEROS(N5,N4))THEN
    CRUN                            = CRUN-WQRS
    Qrunoff_CumYr_col(NY,NX)        = Qrunoff_CumYr_col(NY,NX)-WQRS
    HQRS                            = XN*HeatSnoByRedist_2DH(N,NN,N5,N4)
    HeatOut_lnds                    = HeatOut_lnds-HQRS
    CXS                             = XN*(trcg_FloXSnow_2DH(idg_CO2,N,NN,N5,N4)+trcg_FloXSnow_2DH(idg_CH4,N,NN,N5,N4))
    ZXS                             = XN*(trcn_FloXSnow_2DH(ids_NH4,N,NN,N5,N4)+trcg_FloXSnow_2DH(idg_NH3,N,NN,N5,N4)+trcn_FloXSnow_2DH(ids_NO3,N,NN,N5,N4))
    ZGS                             = XN*(trcg_FloXSnow_2DH(idg_N2O,N,NN,N5,N4)+trcg_FloXSnow_2DH(idg_N2,N,NN,N5,N4))
    PXS                             = XN*(trcn_FloXSnow_2DH(ids_H2PO4,N,NN,N5,N4)+trcn_FloXSnow_2DH(ids_H1PO4,N,NN,N5,N4))
    TOMOU_lnds(ielmc)               = TOMOU_lnds(ielmc)-CXS
    TOMOU_lnds(ielmn)               = TOMOU_lnds(ielmn)-ZXS-ZGS
    TOMOU_lnds(ielmp)               = TOMOU_lnds(ielmp)-PXS
    HydroSufDICFlx_col(NY,NX)       = HydroSufDICFlx_col(NY,NX)-CXR
    HydroSufDINFlx_CumYr_col(NY,NX) = HydroSufDINFlx_CumYr_col(NY,NX)-ZXR-ZGR
    HydroSufDIPFlx_CumYr_col(NY,NX) = HydroSufDIPFlx_CumYr_col(NY,NX)-PXR
    OXS                             = XN*trcg_FloXSnow_2DH(idg_O2,N,NN,N5,N4)
    OXYGOU                          = OXYGOU-OXS
    IF(salt_model)THEN
      PSS=XN*patomw*(trcSalt_FloXSnow_2DH(idsalt_CaPO4,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_FeHPO4,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_CaHPO4,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_MgHPO4,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_H3PO4,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_FeH2PO4,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_CaH4P2O8,N,NN,N5,N4))
      TOMOU_lnds(ielmp)=TOMOU_lnds(ielmp)-PSS
      SS1=XN*(trcSalt_FloXSnow_2DH(idsalt_Al,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_Fe,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_Hp,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_Ca,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_Mg,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_Na,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_K,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_OH,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_SO4,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_Cl,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_CO3,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_H0PO4,N,NN,N5,N4))
      SS2=XN*2.0*(trcSalt_FloXSnow_2DH(idsalt_HCO3,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_AlOH,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_AlSO4,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_FeOH,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_FeSO4,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_CaOH,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_CaCO3,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_CaSO4,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_MgOH2,N,NN,N5,N4)&
        +trcSalt_FloXSnow_2DH(idsalt_MgCO3,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_MgSO4,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_NaCO3,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_NaSO4,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_KSO4,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_CaPO4,N,NN,N5,N4))
      SS3=XN*3.0*(trcSalt_FloXSnow_2DH(idsalt_AlOH2,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_FeOH2,N,NN,N5,N4)&
        +trcSalt_FloXSnow_2DH(idsalt_CaHCO3,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_MgHCO3,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_FeHPO4,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_CaHPO4,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_MgHPO4,N,NN,N5,N4))
      SS4=XN*4.0*(trcSalt_FloXSnow_2DH(idsalt_AlOH3,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_FeOH3,N,NN,N5,N4) &
        +trcSalt_FloXSnow_2DH(idsalt_H3PO4,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_FeH2PO4,N,NN,N5,N4)&
        +trcSalt_FloXSnow_2DH(idsalt_CaH4P2O8,N,NN,N5,N4)) &
        +XN*5.0*(trcSalt_FloXSnow_2DH(idsalt_AlOH4,N,NN,N5,N4)+trcSalt_FloXSnow_2DH(idsalt_FeOH4,N,NN,N5,N4))
      SSR=SS1+SS2+SS3+SS4
      TIONOU=TIONOU-SSR
      HydroIonFlx_CumYr_col(NY,NX)=HydroIonFlx_CumYr_col(NY,NX)-SSR
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine WaterHeatSoluteBySnowDrift

end module RunoffBalMod
