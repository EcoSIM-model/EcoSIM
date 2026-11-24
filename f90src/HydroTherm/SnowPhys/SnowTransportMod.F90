module SnowTransportMod
!
!code to do water and tracer transport in snowpack
!
  use data_kind_mod,  only: r8 => DAT_KIND_R8
  use data_const_mod, only: spval => DAT_CONST_SPVAL
  use abortutils,     only: endrun
  use EcoSimConst,    only: DENSICE, natomw, patomw
  use EcoSIMCtrlMod,  only: lverb
  use MiniMathMod,    only: fixEXConsumpFlux,AZERO,AZMAX1
  use TracerPropMod,  only: MolecularWeight
  use DebugToolMod
  use ClimForcDataType
  use IrrigationDataType
  use SurfSoilDataType
  use SnowPhysData
  use GridDataType
  use SnowDataType
  use SoilBGCDataType
  use AqueChemDatatype
  use SoilWaterDataType
  use SOMDataType
  use ChemTranspDataType
  use EcoSimSumDataType
implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: SaltPercolThruSnow
  public :: DiagSnowChemMass
  public :: SoluteTransportThruSnow

  contains
!------------------------------------------------------------------------------------------

  subroutine SaltPercolThruSnow(I,J,N1,N2)
  !
  !update solute due to snow flux
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: N1,N2

  character(len=*), parameter :: subname='SaltPercolThruSnow'
  integer :: LS, LS2
  integer :: idg,idn,idsalt,ids

!     begin_execution
   call PrintInfo('beg '//subname)

  ! NET WATER AND HEAT FLUXES THROUGH SNOWPACK  
  ! from top to bottom of snow layer 
  D1205: DO LS=1,JS
    !layer LS has significant snow
    IF(VLHeatCapSnow_snvr(LS,N2,N1).GT.VLHeatCapSnowMin_col(N2,N1))THEN
      !id of next snow layer
      LS2=MIN(JS,LS+1)
      
      ! Not surface layer, and is heat significant
      IF(LS.LT.JS .AND. VLHeatCapSnow_snvr(LS2,N2,N1).GT.VLHeatCapSnowMin_col(N2,N1))THEN
        
        ! NET SOLUTE FLUXES THROUGH SNOWPACK
        DO idg=idg_beg,idg_NH3
          trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)=trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)+trcg_AquaAdv_flx_snvr(idg,LS,N2,N1) &
            -trcg_AquaAdv_flx_snvr(idg,LS2,N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)=trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)+trcn_AquaAdv_flx_snvr(idn,LS,N2,N1) &
            -trcn_AquaAdv_flx_snvr(idn,LS2,N2,N1)
        ENDDO
        !
        !     NET SALT FLUXES THROUGH SNOWPACK
        !
        !
        IF(salt_model)THEN
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)=trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)+trcSalt_AquaAdv_flx_snvr(idsalt,LS,N2,N1)&
              -trcSalt_AquaAdv_flx_snvr(idsalt,LS2,N2,N1)
          ENDDO
        ENDIF
        !
        ! IF LOWER LAYER IS THE LITTER AND SOIL SURFACE
        !
      ELSE

        ! and NH3B
        DO idg=idg_beg,idg_NH3
          trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)=trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)+trcg_AquaAdv_flx_snvr(idg,LS,N2,N1) &
            -trcg_AquaADV_Snow2Litr_flx(idg,N2,N1)-trcg_AquaADV_Snow2Soil_flx(idg,N2,N1)
          trcg_snowMassloss_col(idg,N2,N1)  =trcg_snowMassloss_col(idg,N2,N1)+trcg_AquaADV_Snow2Litr_flx(idg,N2,N1)&
            +trcg_AquaADV_Snow2Soil_flx(idg,N2,N1)
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)=trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)+trcn_AquaAdv_flx_snvr(idn,LS,N2,N1) &
            -trcn_AquaADV_Snow2Litr_flx(idn,N2,N1)-trcn_AquaADV_Snow2Soil_flx(idn,N2,N1) 
          trcn_snowMassloss_col(idn,N2,N1)  =trcn_snowMassloss_col(idn,N2,N1)+ &
            trcn_AquaADV_Snow2Litr_flx(idn,N2,N1)+trcn_AquaADV_Snow2Soil_flx(idn,N2,N1) 
        ENDDO

        !add band flux
        trcg_AquaAdv_NetFlx_snvr(idg_NH3,LS,N2,N1)=trcg_AquaAdv_NetFlx_snvr(idg_NH3,LS,N2,N1)-trcg_AquaADV_Snow2Soil_flx(idg_NH3B,N2,N1)
        trcg_snowMassloss_col(idg_NH3,N2,N1)=trcg_snowMassloss_col(idg_NH3,N2,N1)+trcg_AquaADV_Snow2Soil_flx(idg_NH3B,N2,N1)
        DO ids=0,ids_nuts
          trcn_AquaAdv_NetFlx_snvr(ids_NH4+ids,LS,N2,N1)=trcn_AquaAdv_NetFlx_snvr(ids_NH4+ids,LS,N2,N1) &
            -trcn_AquaADV_Snow2Band_flx(ids_NH4B+ids,N2,N1)
          trcn_snowMassloss_col(ids_NH4+ids,N2,N1) = trcn_snowMassloss_col(ids_NH4+ids,N2,N1)  +  &
            trcn_AquaADV_Snow2Band_flx(ids_NH4B+ids,N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)=trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)+trcSalt_AquaAdv_flx_snvr(idsalt,LS,N2,N1) &
              -trcSalt_AquaADV_Snow2Litr_flx(idsalt,N2,N1)-trcSalt_AquaADV_Snow2Soil_flx(idsalt,N2,N1)
            trcSalt_snowMassloss_col(idsalt,N2,N1) = trcSalt_snowMassloss_col(idsalt,N2,N1) + &
              trcSalt_AquaADV_Snow2Litr_flx(idsalt,N2,N1)+trcSalt_AquaADV_Snow2Soil_flx(idsalt,N2,N1)
          ENDDO

          !add band flux
          DO idsalt=0,idsalt_nuts
            trcSalt_AquaAdv_NetFlx_snvr(idsalt_H0PO4+idsalt,LS,N2,N1)=trcSalt_AquaAdv_NetFlx_snvr(idsalt_H0PO4+idsalt,LS,N2,N1) &
              -trcSalt_AquaADV_Snow2Soil_flx(idsalt_H0PO4B+idsalt,N2,N1)

            trcSalt_snowMassloss_col(idsalt_H0PO4+idsalt,N2,N1) =trcSalt_snowMassloss_col(idsalt_H0PO4+idsalt,N2,N1) + &
              trcSalt_AquaADV_Snow2Soil_flx(idsalt_H0PO4B+idsalt,N2,N1)
          ENDDO
        ENDIF
      ENDIF
 
      ! layer LS does not have significant snow
      ! WATER,GAS,SOLUTE,SALT FLUXES INTO SNOWPACK SURFACE   
    ELSEIF(LS.EQ.1)THEN

        DO idg=idg_beg,idg_NH3
          trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)=trcg_AquaAdv_NetFlx_snvr(idg,LS,N2,N1)+trcg_AquaAdv_flx_snvr(idg,LS,N2,N1) &
            -trcg_AquaADV_Snow2Litr_flx(idg,N2,N1)-trcg_AquaADV_Snow2Soil_flx(idg,N2,N1)
          trcg_snowMassloss_col(idg,N2,N1)  =trcg_snowMassloss_col(idg,N2,N1)+trcg_AquaADV_Snow2Litr_flx(idg,N2,N1)&
            +trcg_AquaADV_Snow2Soil_flx(idg,N2,N1)          
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)=trcn_AquaAdv_NetFlx_snvr(idn,LS,N2,N1)+trcn_AquaAdv_flx_snvr(idn,LS,N2,N1) &
            -trcn_AquaADV_Snow2Litr_flx(idn,N2,N1)-trcn_AquaADV_Snow2Soil_flx(idn,N2,N1) 
          trcn_snowMassloss_col(idn,N2,N1)  =trcn_snowMassloss_col(idn,N2,N1)+ &
            trcn_AquaADV_Snow2Litr_flx(idn,N2,N1)+trcn_AquaADV_Snow2Soil_flx(idn,N2,N1) 
        ENDDO

        !add band flux
        trcg_AquaAdv_NetFlx_snvr(idg_NH3,LS,N2,N1)=trcg_AquaAdv_NetFlx_snvr(idg_NH3,LS,N2,N1)-trcg_AquaADV_Snow2Soil_flx(idg_NH3B,N2,N1)
        trcg_snowMassloss_col(idg_NH3,N2,N1)=trcg_snowMassloss_col(idg_NH3,N2,N1)+trcg_AquaADV_Snow2Soil_flx(idg_NH3B,N2,N1)
        DO ids=0,ids_nuts
          trcn_AquaAdv_NetFlx_snvr(ids_NH4+ids,LS,N2,N1)=trcn_AquaAdv_NetFlx_snvr(ids_NH4+ids,LS,N2,N1) &
            -trcn_AquaADV_Snow2Band_flx(ids_NH4B+ids,N2,N1)
          trcn_snowMassloss_col(ids_NH4+ids,N2,N1) = trcn_snowMassloss_col(ids_NH4+ids,N2,N1)  +  &
            trcn_AquaADV_Snow2Band_flx(ids_NH4B+ids,N2,N1)
        ENDDO

        IF(salt_model)THEN
          DO idsalt=idsalt_beg,idsalt_end
            trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)=trcSalt_AquaAdv_NetFlx_snvr(idsalt,LS,N2,N1)+trcSalt_AquaAdv_flx_snvr(idsalt,LS,N2,N1) &
              -trcSalt_AquaADV_Snow2Litr_flx(idsalt,N2,N1)-trcSalt_AquaADV_Snow2Soil_flx(idsalt,N2,N1)
            trcSalt_snowMassloss_col(idsalt,N2,N1) = trcSalt_snowMassloss_col(idsalt,N2,N1) + &
              trcSalt_AquaADV_Snow2Litr_flx(idsalt,N2,N1)+trcSalt_AquaADV_Snow2Soil_flx(idsalt,N2,N1)            
          ENDDO

          !add band flux
          DO idsalt=0,idsalt_nuts
            trcSalt_AquaAdv_NetFlx_snvr(idsalt_H0PO4+idsalt,LS,N2,N1)=trcSalt_AquaAdv_NetFlx_snvr(idsalt_H0PO4+idsalt,LS,N2,N1) &
              -trcSalt_AquaADV_Snow2Soil_flx(idsalt_H0PO4B+idsalt,N2,N1)

            trcSalt_snowMassloss_col(idsalt_H0PO4+idsalt,N2,N1) =trcSalt_snowMassloss_col(idsalt_H0PO4+idsalt,N2,N1) + &
              trcSalt_AquaADV_Snow2Soil_flx(idsalt_H0PO4B+idsalt,N2,N1)
          ENDDO
        ENDIF
      
    ENDIF
  ENDDO D1205
  call PrintInfo('end '//subname)
  end subroutine SaltPercolThruSnow

!------------------------------------------------------------------------------------------

  subroutine DiagSnowChemMass(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8) :: SSW,WS
  integer :: L,nsalts

  !
  !     UPDATE STATE VARIABLES WITH TOTAL FLUXES CALCULATED ABOVE
  !
  !     IF(J.EQ.24)THEN
  !
  !     SNOWPACK VARIABLES NEEDED FOR WATER, C, N, P, O, SOLUTE AND
  !     ENERGY BALANCES INCLUDING SUM OF ALL CURRENT STATE VARIABLES,
  !     CUMULATIVE SUMS OF ALL ADDITIONS AND REMOVALS
  !
  DO  L=1,JS
    WS                  = VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
    if(WS>0._r8)then
      TGasC_lnd           = TGasC_lnd+trcg_solsml_snvr(idg_CO2,L,NY,NX)+trcg_solsml_snvr(idg_CH4,L,NY,NX)
      DIC_mass_col(NY,NX) = DIC_mass_col(NY,NX)+trcg_solsml_snvr(idg_CO2,L,NY,NX)+trcg_solsml_snvr(idg_CH4,L,NY,NX)
      TSoilO2G_lnd        = TSoilO2G_lnd+trcg_solsml_snvr(idg_O2,L,NY,NX)
      TGasN_lnd           = TGasN_lnd+trcg_solsml_snvr(idg_N2,L,NY,NX)+trcg_solsml_snvr(idg_N2O,L,NY,NX)
      TDisolNH4_lnd       = TDisolNH4_lnd+trcn_solsml_snvr(ids_NH4,L,NY,NX)+trcg_solsml_snvr(idg_NH3,L,NY,NX)
      tNO3_lnd            = tNO3_lnd+trcn_solsml_snvr(ids_NO3,L,NY,NX)
      TDisolPi_lnd        = TDisolPi_lnd+trcn_solsml_snvr(ids_H1PO4,L,NY,NX)+trcn_solsml_snvr(ids_H2PO4,L,NY,NX)

      IF(salt_model)THEN
        SSW=0._r8
        do nsalts=idsalt_beg,idsalt_end
          SSW=SSW+trcSalt_ml_snvr(nsalts,L,NY,NX)*trcSaltIonNumber(nsalts)
        ENDDO  
        TION=TION+SSW
      ENDIF
    ENDIF
  ENDDO
  end subroutine DiagSnowChemMass
!------------------------------------------------------------------------------------------

  subroutine ChemicalBySnowRedistribution(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX

  integer :: idsalt,idg,ids
!     begin_execution
!     OVERLAND SNOW REDISTRIBUTION
!
  IF(abs(TDrysnoByRedist_col(NY,NX))>0._r8)THEN

    !loss of dissolved gases from surface snow
    DO idg=idg_beg,idg_NH3
      trcg_solsml_snvr(idg,1,NY,NX)=trcg_solsml_snvr(idg,1,NY,NX)+trcg_LossXSnowRedist_col(idg,NY,NX)
      GasHydroLoss_flx_col(idg,NY,NX)=GasHydroLoss_flx_col(idg,NY,NX)+trcg_LossXSnowRedist_col(idg,NY,NX)
    ENDDO

    DO ids=ids_nut_beg,ids_nuts_end
      trcn_solsml_snvr(ids,1,NY,NX)=trcn_solsml_snvr(ids,1,NY,NX)+trcn_LossXSnowRedist_col(ids,NY,NX)
    ENDDO

    IF(salt_model)THEN
      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_ml_snvr(idsalt,1,NY,NX)=trcSalt_ml_snvr(idsalt,1,NY,NX)+trcSalt_LossXSnowRedist_col(idsalt,NY,NX)
      ENDDO
    ENDIF
  ENDIF

  end subroutine ChemicalBySnowRedistribution

!------------------------------------------------------------------------------------------

  subroutine SoluteTransportThruSnow(I,J,L,NY,NX)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: L,NY,NX
  character(len=*), parameter :: subname='SoluteTransportThruSnow'
  integer :: idg,idn,idsalt
  !     begin_execution
  !     SNOWPACK SOLUTE CONTENT

  call PrintInfo('beg '//subname)
  DO idg=idg_beg,idg_NH3

    call fixEXConsumpFlux(trcg_solsml_snvr(idg,L,NY,NX),trcg_AquaAdv_NetFlx_snvr(idg,L,NY,NX),-1)
    !trcg_solsml_snvr(idg,L,NY,NX)=AZERO(trcg_solsml_snvr(idg,L,NY,NX))
    trcg_solsml_snvr(idg,L,NY,NX)=AZMAX1(trcg_solsml_snvr(idg,L,NY,NX))
  ENDDO

  if(trcg_solsml_snvr(idg_O2,L,NY,NX)<0._r8)then
    write(116,*)I*1000+J,'snow',L,trcg_solsml_snvr(idg_O2,L,NY,NX),SnowThickL_snvr(L,NY,NX),trcg_AquaAdv_NetFlx_snvr(idg_O2,L,NY,NX)  
    call endrun('negative tracer in '//trim(mod_filename)//' at line',__LINE__)
  endif

  DO idn =ids_nut_beg,ids_nuts_end
    call fixEXConsumpFlux(trcn_solsml_snvr(idn,L,NY,NX),trcn_AquaAdv_NetFlx_snvr(idn,L,NY,NX),-1)
    trcn_solsml_snvr(idn,L,NY,NX)=AZMAX1(trcn_solsml_snvr(idn,L,NY,NX))
    trcn_solsml_snvr(idn,L,NY,NX)=AZERO(trcn_solsml_snvr(idn,L,NY,NX))
  ENDDO
  !
  !
  IF(salt_model)THEN
    DO idsalt=idsalt_beg,idsalt_end
      call fixEXConsumpFlux(trcSalt_ml_snvr(idsalt,L,NY,NX),trcSalt_AquaAdv_NetFlx_snvr(idsalt,L,NY,NX),-1)
      trcSalt_ml_snvr(idsalt,L,NY,NX)=AZERO(trcSalt_ml_snvr(idsalt,L,NY,NX))
    ENDDO
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine SoluteTransportThruSnow

end module SnowTransportMod