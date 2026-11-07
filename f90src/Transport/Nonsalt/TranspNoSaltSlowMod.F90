module TranspNoSaltSlowMod
  ! Compute slow gas/solute transport and exchange among soil, litter and atmosphere.
  use SoilFaceMappingMod, only: SoilFaceSlot,SoilFaceIncomingRange,GetSoilFaceGeometry
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,       only: destroy, endrun, iulog
  use TracerPropMod,    only: MolecularWeight
  use EcoSIMCtrlMod,    only: iVerbLevel
  use EcoSiMParDataMod, only: micpar
  use NumericalAuxMod
  use PlantDataRateType
  use minimathmod    
  use IrrigationDataType
  use DebugToolMod
  use GridConsts
  use EcoSIMSolverPar
  use TranspNoSaltDataMod
  use AqueChemDatatype
  use GridDataType
  use SOMDataType
  use SoilBGCDataType
  use SurfSoilDataType
  use SurfLitterDataType
  use SnowDataType
  use EcoSIMCtrlDataType
  use SoilWaterDataType
  use LandSurfDataType
  use SoilPropertyDataType
  use ChemTranspDataType
  use ClimForcDataType
  USE EcoSimConst
  USE SoilHeatDataType
  use TracerIDMod
implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: TransptSlowNoSaltM
  public :: BubbleEffluxM

  contains
!------------------------------------------------------------------------------------------
  subroutine TransptSlowNoSaltM(I,J,M,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='TransptSlowNoSaltM'

  real(r8) :: pscal(ids_beg:ids_end), dpscal(ids_beg:ids_end)
  real(r8) :: pscal1_dom(idom_beg:idom_end),dpscal_dom(idom_beg:idom_end)
  real(r8) :: dpscal_max
  integer :: ids,idom,iterm
  integer :: snow_donor_col(JY,JX) ! selected snow discharge donor; zero means none
  
  call PrintInfo('beg '//subname)

  call EnterMassCheck(I,J,NHE,NHW,NVS,NVN)

  dpscal=1._r8;dpscal_dom=1._r8
  iterm=0; dpscal_max=1._r8

  call StageSlowTranspIterationM(I,J,M,NHE,NHW,NVS,NVN)

  ! Accept ID-changing surface transfers once, before independent soil transport.
  call ZeroFluxesM(I,J,M,NHE,NHW,NVS,NVN)
  call SnowSoluteVertDrainM(I,J,M,NHE,NHW,NVS,NVN,snow_donor_col,band_transfer_only=.true.)
  call SurfaceSoluteFluxM(I,J,M,NHE,NHW,NVS,NVN,snow_donor_col)

  DO while(dpscal_max>1.e-2_r8 .and. iterm<3)
    iterm=iterm+1
    pscal=1.000_r8
!    write(555,*)'============================'
    call ZeroFluxesM(I,J,M,NHE,NHW,NVS,NVN)  

    call SnowSoluteVertDrainM(I,J,M,NHE,NHW,NVS,NVN,snow_donor_col)

    call SurfaceSoluteFluxM(I,J,M,NHE,NHW,NVS,NVN)

    call TracerXBoundariesM(I,J,M,NHW,NHE,NVN,NVS,snow_donor_col)

    call TracerXInterGridsM(I,J,M,NHE,NHW,NVS,NVN)

    call GatherTranspFluxM(I,J,M,NHW,NHE,NVN,NVS)

    call SlowUpdateStateVarsM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal,pscal1_dom)

    call AccumSlowFluxesM(I,J,M,NHE,NHW,NVS,NVN,dpscal,dpscal_dom,pscal,pscal1_dom)

    dpscal_max=0._r8
    DO ids=ids_beg,ids_end      
      if(pscal(ids)>1.e-6_r8)then
        dpscal(ids)=dpscal(ids)*(1._r8-pscal(ids))
      else
        dpscal(ids)=0._r8
      endif
      dpscal_max=AMAX1(dpscal_max,dpscal(ids))
    enddo
    
    DO idom=idom_beg,idom_end
      dpscal_dom(idom)=dpscal_dom(idom)*(1._r8-pscal1_dom(idom))
      dpscal_max=AMAX1(dpscal_max,dpscal_dom(idom))
    ENDDO
  ENDDO
  call ApplySlowSourceSink(I,J,M,NHE,NHW,NVS,NVN,iterm)

  call ExitMassCheck(I,J,M,NHE,NHW,NVS,NVN,iterm)
  call PrintInfo('end '//subname)
  end subroutine TransptSlowNoSaltM
!------------------------------------------------------------------------------------------
  subroutine ApplySlowSourceSink(I,J,M,NHE,NHW,NVS,NVN,iterm)
  implicit none
  integer, intent(in) :: I,J,M,NHE,NHW,NVS,NVN,iterm
  integer :: NY,NX,idg,K,L,idom,ids
  real(r8) :: flux

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      idg=idg_NH3
      trcs_solml2_vr(idg,0,NY,NX)=trcs_solml2_vr(idg,0,NY,NX)+RBGCSrcSoluteM_vr(idg,0,NY,NX)
      call SubstrateDribbling(RBGCSinkSoluteM_vr(idg,0,NY,NX),trcs_solml_drib_vr(idg,0,NY,NX),trcs_solml2_vr(idg,0,NY,NX))

      flux=RBGCSrcSoluteM_vr(idg,0,NY,NX)-RBGCSinkSoluteM_vr(idg,0,NY,NX)
      trcs_NetProd_slow_flxM_col(idg,NY,NX) = trcs_NetProd_slow_flxM_col(idg,NY,NX)+flux
      trcs_netProd_lit_col(idg,NY,NX)       = trcs_netProd_lit_col(idg,NY,NX)+flux

      DO ids=ids_nut_beg,ids_nuts_end    
        trcs_solml2_vr(ids,0,NY,NX)=trcs_solml2_vr(ids,0,NY,NX)+RBGCSrcSoluteM_vr(ids,0,NY,NX)
        call SubstrateDribbling(RBGCSinkSoluteM_vr(ids,0,NY,NX),trcs_solml_drib_vr(ids,0,NY,NX),trcs_solml2_vr(ids,0,NY,NX))

        flux=RBGCSrcSoluteM_vr(ids,0,NY,NX)-RBGCSinkSoluteM_vr(ids,0,NY,NX)
        trcs_NetProd_slow_flxM_col(ids,NY,NX) = trcs_NetProd_slow_flxM_col(ids,NY,NX)+flux
        trcs_netProd_lit_col(ids,NY,NX)       = trcs_netProd_lit_col(ids,NY,NX)+flux
      ENDDO 

      DO K=1,micpar%NumOfLitrCmplxs
        DO idom=idom_beg,idom_end        
          call SubstrateDribbling(RBGCSink_DOM_micpM_vr(idom,K,0,NY,NX),DOM_MicP_drib_vr(idom,K,0,NY,NX),DOM_MicP2_vr(idom,K,0,NY,NX))
          DOM_NetProd_slow_flxM_col(idom,K,NY,NX)=DOM_NetProd_slow_flxM_col(idom,K,NY,NX)-RBGCSink_DOM_micpM_vr(idom,K,0,NY,NX)

        ENDDO
      ENDDO  

      DO L=NU_col(NY,NX),NL_col(NY,NX)    
        IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN        
          DO ids=idg_NH3,ids_end
            trcs_solml2_vr(ids,L,NY,NX)=trcs_solml2_vr(ids,L,NY,NX)+RBGCSrcSoluteM_vr(ids,L,NY,NX) 
            call SubstrateDribbling(RBGCSinkSoluteM_vr(ids,L,NY,NX),trcs_solml_drib_vr(ids,L,NY,NX),trcs_solml2_vr(ids,L,NY,NX)) 

            flux=RBGCSrcSoluteM_vr(ids,L,NY,NX)-RBGCSinkSoluteM_vr(ids,L,NY,NX)
            trcs_NetProd_slow_flxM_col(ids,NY,NX)  = trcs_NetProd_slow_flxM_col(ids,NY,NX)+flux
            trcs_NetProd_slow_soil_col(ids,NY,NX)  = trcs_NetProd_slow_soil_col(ids,NY,NX)+flux

            if(ids==idg_NH3)RGasNetProdSoil_col(ids,NY,NX)=RGasNetProdSoil_col(ids,NY,NX)+flux
          ENDDO
        ENDIF

        DO K=1,jcplx
          DO idom=idom_beg,idom_end        
            call SubstrateDribbling(RBGCSink_DOM_micpM_vr(idom,K,L,NY,NX),DOM_MicP_drib_vr(idom,K,L,NY,NX),DOM_MicP2_vr(idom,K,L,NY,NX))
            DOM_NetProd_slow_flxM_col(idom,K,NY,NX)=DOM_NetProd_slow_flxM_col(idom,K,NY,NX)-RBGCSink_DOM_micpM_vr(idom,K,L,NY,NX)
    
          ENDDO
        ENDDO  
      ENDDO

    ENDDO
  ENDDO  
  end subroutine ApplySlowSourceSink
!------------------------------------------------------------------------------------------
  subroutine EnterMassCheck(I,J,NHE,NHW,NVS,NVN)
  !
  !Description:
  !Prepare for mass conservation check 
  !
  implicit none
  integer, intent(in) :: I,J,NHE,NHW,NVS,NVN

  integer :: NY,NX,idg,L,K,idom,ids

!  write(555,*)'enter'
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      trcg_mass_begs(:,NY,NX)                 = 0._r8
      DOM_mass_begs(:,:,NY,NX)                = 0._r8
      DOM_Hydroloss_slow_flx_col(:,:,NY,NX)   = 0._r8
      TranspNetDOM_flx_col(:,:,NY,NX)         = 0._r8
      TranspNetDOM_flx2_col(:,:,NY,NX)        = 0._r8
      AtmGasDiff2Surf_slow_flx_col(:,NY,NX)   = 0._r8
      AtmGasDiff2Litr_slow_flx_col(:,NY,NX)   = 0._r8
      Gas_WetDeposit_slow_flx_col(:,NY,NX)    = 0._r8
      Gas_WetDepo2Snow_slow_flx_col(:,NY,NX)  = 0._r8
      Gas_Snowloss_slow_flx_col(:,NY,NX)      = 0._r8
      trcs_NetFlow2Litr_slow_flx_col(:,NY,NX) = 0._r8
      trcs_hydrloss_slow_flx_col(:,NY,NX)     = 0._r8
      trcs_NetProd_slow_flxM_col(:,NY,NX)     = 0._r8
      trcs_NetProd_slow_soil_col(:,NY,NX)     = 0._r8
      DOM_NetProd_slow_flxM_col(:,:,NY,NX)    = 0._r8
      trcg_mass_snow_begs(:,NY,NX)            = 0._r8
      trcg_mass_soil_begs(:,NY,NX)            = 0._r8
      trcs_netflow2soil_slow_flx_col(:,NY,NX) = 0._r8
      trcs_solml_dribM_beg_col(:,NY,NX)       = 0._r8
      trcs_drib_litr_slow_begf(:,NY,NX)       = 0._r8
      trcs_drib_soil_slow_begf(:,NY,NX)       = 0._r8
      GasHydroLoss_litr_slow_flx_col(:,NY,NX) = 0._r8
      trcs_drainage_slow_flx_col(:,NY,NX)     = 0._r8
      Gas_WetDepo2litr_slow_flx_col(:,NY,NX)  = 0._r8
      Gas_WeDepo2soil_slow_flx_col(:,NY,NX)   = 0._r8
      Gas_litr2Soil_slow_flx_col(:,NY,NX)     = 0._r8
      trcs_sno2litr_slow_flx_col(:,NY,NX)     = 0._r8
      dom_dribble_slow_beg_col(:,:,NY,NX)     = 0._r8
      do idg=idg_beg,idg_NH3
        trcs_drib_litr_slow_begf(idg,NY,NX)  = trcs_solml_drib_vr(idg,0,NY,NX)
      enddo

      do ids=ids_beg,ids_end
        trcs_solml_dribM_beg_col(ids,NY,NX) = trcs_solml_drib_vr(ids,0,NY,NX)
      enddo

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO ids=ids_beg,ids_end
          trcs_solml_dribM_beg_col(ids,NY,NX) = trcs_solml_dribM_beg_col(ids,NY,NX)+trcs_solml_drib_vr(ids,L,NY,NX)
          trcs_drib_soil_slow_begf(ids,NY,NX) = trcs_drib_soil_slow_begf(ids,NY,NX)+trcs_solml_drib_vr(ids,L,NY,NX)
        ENDDO  
      ENDDO

      DO idg=idg_beg,idg_NH3
        DO L=1,JS
          trcg_mass_snow_begs(idg,NY,NX)=trcg_mass_snow_begs(idg,NY,NX)+trcg_solsml2_snvr(idg,L,NY,NX)
        ENDDO

        trcg_mass_litr_begs(idg,NY,NX)     = trcs_solml2_vr(idg,0,NY,NX)

        DO L=NU_col(NY,NX),NL_col(NY,NX)
          trcg_mass_soil_begs(idg,NY,NX)     = trcg_mass_soil_begs(idg,NY,NX)+trcs_solml2_vr(idg,L,NY,NX)+trcs_soHml2_vr(idg,L,NY,NX)
        ENDDO 
      ENDDO

      idg=idg_NH3
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        trcg_mass_soil_begs(idg,NY,NX)=trcg_mass_soil_begs(idg,NY,NX)+trcs_solml2_vr(idg_NH3B,L,NY,NX)+trcs_soHml2_vr(idg_NH3B,L,NY,NX)
      ENDDO       

      DO idg=idg_beg,idg_NH3
        trcg_mass_begs(idg,NY,NX)=trcg_mass_snow_begs(idg,NY,NX)+trcg_mass_litr_begs(idg,NY,NX)+trcg_mass_soil_begs(idg,NY,NX)
      ENDDO

      DO K=1,micpar%NumOfLitrCmplxs
        DO idom=idom_beg,idom_end
          DOM_mass_begs(idom,K,NY,NX)=DOM_mass_begs(idom,K,NY,NX)+DOM_MicP2_vr(idom,K,0,NY,NX)            
          dom_dribble_slow_beg_col(idom,K,NY,NX)=dom_dribble_slow_beg_col(idom,K,NY,NX)+DOM_MicP_drib_vr(idom,K,0,NY,NX)
        ENDDO
      ENDDO

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO K=1,jcplx
          DO idom=idom_beg,idom_end
            DOM_mass_begs(idom,K,NY,NX)=DOM_mass_begs(idom,K,NY,NX)+DOM_MicP2_vr(idom,K,L,NY,NX)+DOM_MacP2_vr(idom,K,L,NY,NX)
            dom_dribble_slow_beg_col(idom,K,NY,NX)=dom_dribble_slow_beg_col(idom,K,NY,NX)+DOM_MicP_drib_vr(idom,K,L,NY,NX)            
          ENDDO
        ENDDO
      ENDDO
    ENDDO  
  ENDDO
  end subroutine EnterMassCheck
!------------------------------------------------------------------------------------------


  subroutine ExitMassCheck(I,J,M,NHE,NHW,NVS,NVN,iterm)
  !
  ! Description:
  ! Do slow transport mass conservation check.
  ! fluxes include:
  ! hydrological loss, snow drift, (sub)surface runoff, drainage : trcs_hydrloss_slow_flx_col(idg,NY,NX)
  ! net cheimcal production, and irrigation                      : trcs_NetProd_slow_flxM_col(idg,NY,NX)
  ! surface diffusion, dissolution into litter and soil          : AtmGasDiff2Surf_slow_flx_col(idg,NY,NX)
  ! wet deposition to snow, litter and topsoil                   : Gas_WetDeposit_slow_flx_col(idg,NY,NX)
  
  implicit none
  integer, intent(in) :: I,J,NHE,NHW,NVS,NVN,M,iterm

  integer :: NY,NX,idg,L,K,idom,ids
  real(r8) :: trcg_mass_now(idg_beg:idg_NH3)
  real(r8) :: trcg_mass_snow_now(idg_beg:idg_NH3)
  real(r8) :: trcg_mass_litr_now(idg_beg:idg_NH3)
  real(r8) :: trcg_mass_soil_now(idg_beg:idg_NH3)
  real(r8) :: dmass,err,err_rel
  real(r8) :: DOM_mass_now(idom_beg:idom_end,jcplx)
  real(r8) :: dom_dribble_col(idom_beg:idom_end,jcplx)
  real(r8) :: trcs_solml_drib_col(ids_beg:ids_end)
  real(r8) :: trcs_drib_soil_slow_col(ids_beg:ids_end)
  logical  :: err_test

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      dom_dribble_col         = 0._r8
      trcg_mass_snow_now      = 0._r8
      trcg_mass_soil_now      = 0._r8
      trcs_solml_drib_col     = 0._r8
      trcs_drib_soil_slow_col = 0._r8
      DO ids = ids_beg, ids_end
        trcs_solml_drib_col(ids) = trcs_solml_drib_vr(ids,0,NY,NX)
      ENDDO

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO ids=ids_beg,ids_end
          trcs_solml_drib_col(ids) = trcs_solml_drib_col(ids)+trcs_solml_drib_vr(ids,L,NY,NX)
          trcs_drib_soil_slow_col(ids)=trcs_drib_soil_slow_col(ids)+trcs_solml_drib_vr(ids,L,NY,NX)
        ENDDO
      ENDDO

      DO idg=idg_beg,idg_NH3
        DO L=1,JS
          trcg_mass_snow_now(idg)=trcg_mass_snow_now(idg)+trcg_solsml2_snvr(idg,L,NY,NX)
        ENDDO

        trcg_mass_litr_now(idg)  = trcs_solml2_vr(idg,0,NY,NX)

        DO L=NU_col(NY,NX),NL_col(NY,NX)
          trcg_mass_soil_now(idg)  = trcg_mass_soil_now(idg)+trcs_solml2_vr(idg,L,NY,NX)+trcs_soHml2_vr(idg,L,NY,NX)
        ENDDO 
      ENDDO
     
      idg=idg_NH3
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        trcg_mass_soil_now(idg)=trcg_mass_soil_now(idg)+trcs_solml2_vr(idg_NH3B,L,NY,NX)+trcs_soHml2_vr(idg_NH3B,L,NY,NX)
      ENDDO       

      DO idg=idg_beg,idg_NH3
        trcg_mass_now(idg)=trcg_mass_snow_now(idg)+trcg_mass_litr_now(idg)+trcg_mass_soil_now(idg)
      ENDDO
      DOM_mass_now(:,:)=0._r8
      DO K=1,micpar%NumOfLitrCmplxs
        DO idom=idom_beg,idom_end
          DOM_mass_now(idom,K)    = DOM_MicP2_vr(idom,K,0,NY,NX)
          dom_dribble_col(idom,K) = dom_dribble_col(idom,K)+DOM_MicP_drib_vr(idom,K,0,NY,NX)
        ENDDO        
      ENDDO
      
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO K=1,jcplx
          DO idom=idom_beg,idom_end
            DOM_mass_now(idom,K)    = DOM_mass_now(idom,K)+DOM_MicP2_vr(idom,K,L,NY,NX)+DOM_MacP2_vr(idom,K,L,NY,NX)
            dom_dribble_col(idom,K) = dom_dribble_col(idom,K)+DOM_MicP_drib_vr(idom,K,L,NY,NX)
          ENDDO
        ENDDO
      ENDDO
      do ids=ids_nuts_beg,ids_nuts_end
        trcnuts_hydrloss_flx_col(ids,NY,NX)=trcnuts_hydrloss_flx_col(ids,NY,NX)+trcs_hydrloss_slow_flx_col(ids,NY,NX)
      enddo
      DO idg=idg_beg,idg_NH3
        GasDiff2Surf_flx_col(idg,NY,NX)       = GasDiff2Surf_flx_col(idg,NY,NX)+AtmGasDiff2Surf_slow_flx_col(idg,NY,NX)
        GasHydroLoss_flx_col(idg,NY,NX)       = GasHydroLoss_flx_col(idg,NY,NX)+trcs_hydrloss_slow_flx_col(idg,NY,NX)
        Gas_WetDeposit_flx_col(idg,NY,NX)     = Gas_WetDeposit_flx_col(idg,NY,NX)+Gas_WetDeposit_slow_flx_col(idg,NY,NX)
        RGasNetProd_col(idg,NY,NX)            = RGasNetProd_col(idg,NY,NX)+trcs_NetProd_slow_flxM_col(idg,NY,NX)
        AtmGasDiff2Litr_flx_col(idg,NY,NX)    = AtmGasDiff2Litr_flx_col(idg,NY,NX)+AtmGasDiff2Litr_slow_flx_col(idg,NY,NX)
        Gas_litr2Soil_flx_col(idg,NY,NX)      = Gas_litr2Soil_flx_col(idg,NY,NX)+Gas_litr2Soil_slow_flx_col(idg,NY,NX)
        trcg_AquaADV_Snow2Litr_flx(idg,NY,NX) = trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)+trcs_sno2litr_slow_flx_col(idg,NY,NX)

        !loss through hydrological runoff
        GasHydroLoss_litr_flx_col(idg,NY,NX) = GasHydroLoss_litr_flx_col(idg,NY,NX)+GasHydroLoss_litr_slow_flx_col(idg,NY,NX)
        trcs_drainage_flx_col(idg,NY,NX)     = trcs_drainage_flx_col(idg,NY,NX)+trcs_drainage_slow_flx_col(idg,NY,NX)

        dmass = trcg_mass_begs(idg,NY,NX)-trcg_mass_now(idg)
        err   = dmass+trcs_hydrloss_slow_flx_col(idg,NY,NX)+trcs_NetProd_slow_flxM_col(idg,NY,NX) &
          +AtmGasDiff2Surf_slow_flx_col(idg,NY,NX)+Gas_WetDeposit_slow_flx_col(idg,NY,NX) &
          +trcs_solml_drib_col(idg)-trcs_solml_dribM_beg_col(idg,NY,NX)
          
        if(idg==idg_NH3)then
          !add flux from NH3B
          GasHydroLoss_flx_col(idg,NY,NX) = GasHydroLoss_flx_col(idg,NY,NX)+trcs_hydrloss_slow_flx_col(idg_NH3B,NY,NX)
          RGasNetProd_col(idg,NY,NX)      = RGasNetProd_col(idg,NY,NX)+trcs_NetProd_slow_flxM_col(idg_NH3B,NY,NX)

          err=err+trcs_hydrloss_slow_flx_col(idg_NH3B,NY,NX)+trcs_NetProd_slow_flxM_col(idg_NH3B,NY,NX)
          err=err+trcs_solml_drib_col(idg_NH3B)-trcs_solml_dribM_beg_col(idg_NH3B,NY,NX)

        endif

        errmass_slow(idg,NY,NX)=errmass_slow(idg,NY,NX)+err
        err_rel=safe_adb(err,AMAX1(trcg_mass_begs(idg,NY,NX),trcg_mass_now(idg)))
        err_test=AMIN1(abs(err),abs(err_rel))>1.e-5_r8 .or. abs(errmass_slow(idg,NY,NX))>1.e-5_r8 .OR. iVerbLevel==1 
        if(err_test)then          
            write(201,*)('-',L=1,50)
            write(201,*)(I*1000+J)*10+M,'iterm=',iterm,trcs_names(idg),'NY,NX=',NY,NX,'SLOW'
            write(201,*)'beg/end total mass',trcg_mass_begs(idg,NY,NX),trcg_mass_now(idg)
            write(201,*)'err,err_rel=',err,err_rel,'nsnol_col=',nsnol_col(NY,NX),'NU_col=',NU_col(NY,NX)
            write(201,*)'beg/end snow mass=',trcg_mass_snow_begs(idg,NY,NX),trcg_mass_snow_now(idg),trcg_mass_snow_begs(idg,NY,NX)-trcg_mass_snow_now(idg)
            write(201,*)'beg/end litr mass=',trcg_mass_litr_begs(idg,NY,NX),trcg_mass_litr_now(idg),trcg_mass_litr_begs(idg,NY,NX)-trcg_mass_litr_now(idg)
            write(201,*)'beg/end soil mass=',trcg_mass_soil_begs(idg,NY,NX),trcg_mass_soil_now(idg),trcg_mass_soil_begs(idg,NY,NX)-trcg_mass_soil_now(idg)
            write(201,*)'wetdepo    =',Gas_WetDeposit_slow_flx_col(idg,NY,NX)
            write(201,*)'diffus     =',AtmGasDiff2Surf_slow_flx_col(idg,NY,NX)
            write(201,*)'drib       =',trcs_solml_drib_col(idg)
            write(201,*)'netflx2litr=',trcs_NetFlow2Litr_slow_flx_col(idg,NY,NX)
            write(201,*)'dep2sno    =',Gas_WetDepo2Snow_slow_flx_col(idg,NY,NX)
            write(201,*)'snowloss   =',-Gas_Snowloss_slow_flx_col(idg,NY,NX)
            write(201,*)'snowdrift  =',trcg_SnowDrift_flx_col(idg,NY,NX)
            write(201,*)'snow mass  =',trcg_solsml2_snvr(idg,1:JS,NY,NX)
            if(idg==idg_NH3)then
              write(201,*)'netpro      =',trcs_NetProd_slow_flxM_col(idg,NY,NX)+trcs_NetProd_slow_flxM_col(idg_NH3B,NY,NX)
              write(201,*)'hydloss     =',trcs_hydrloss_slow_flx_col(idg,NY,NX)+trcs_hydrloss_slow_flx_col(idg_NH3B,NY,NX)
              write(201,*)'netflx2soil =',trcs_netflow2soil_slow_flx_col(idg_NH3,NY,NX)+trcs_netflow2soil_slow_flx_col(idg_NH3B,NY,NX)
            else
              write(201,*)'netpro      =',trcs_NetProd_slow_flxM_col(idg,NY,NX)
              write(201,*)'hydloss     =',trcs_hydrloss_slow_flx_col(idg,NY,NX)            
              write(201,*)'netflx2soil =',trcs_netflow2soil_slow_flx_col(idg,NY,NX)
            endif
            write(201,*)('-',L=1,50)
            write(201,*)'err snow       =',trcg_mass_snow_begs(idg,NY,NX)-trcg_mass_snow_now(idg)+Gas_WetDepo2Snow_slow_flx_col(idg,NY,NX) &
              -Gas_Snowloss_slow_flx_col(idg,NY,NX)
            write(201,*)('-',L=1,50)
            if(idg==idg_NH3)then
              write(201,*)'err soil      =',trcg_mass_soil_begs(idg,NY,NX)-trcg_mass_soil_now(idg)-trcs_drib_soil_slow_begf(idg,NY,NX) + &
                trcs_drib_soil_slow_col(idg)-trcs_drib_soil_slow_begf(idg_NH3B,NY,NX)+trcs_drib_soil_slow_col(idg_NH3B) + &
                trcs_netflow2soil_slow_flx_col(idg_NH3,NY,NX)+trcs_netflow2soil_slow_flx_col(idg_NH3B,NY,NX) + &
                trcs_NetProd_slow_soil_col(idg_NH3,NY,NX)+trcs_NetProd_slow_soil_col(idg_NH3B,NY,NX)
              write(201,*)'drib beg/end  =',trcs_drib_soil_slow_begf(idg,NY,NX)+trcs_drib_soil_slow_begf(idg_NH3B,NY,NX),&
                trcs_drib_soil_slow_col(idg)+trcs_drib_soil_slow_col(idg_NH3B)
              write(201,*)'netflo2soil   =',trcs_netflow2soil_slow_flx_col(idg,NY,NX)+trcs_netflow2soil_slow_flx_col(idg_NH3B,NY,NX)
              write(201,*)'netprod2soil NH3,NH3B =',trcs_NetProd_slow_soil_col(idg,NY,NX),trcs_NetProd_slow_soil_col(idg_NH3B,NY,NX)                 
              write(201,*)'wetdep2soil   =',Gas_WeDepo2soil_slow_flx_col(idg,NY,NX)              
            else
              write(201,*)'err soil      =',trcg_mass_soil_begs(idg,NY,NX)-trcg_mass_soil_now(idg)-trcs_drib_soil_slow_begf(idg,NY,NX) + &
                trcs_drib_soil_slow_col(idg) + trcs_netflow2soil_slow_flx_col(idg,NY,NX)+trcs_NetProd_slow_soil_col(idg,NY,NX) + &
                Gas_WeDepo2soil_slow_flx_col(idg,NY,NX)
              write(201,*)'drib beg/end  =',trcs_drib_soil_slow_begf(idg,NY,NX),trcs_drib_soil_slow_col(idg)
              write(201,*)'netflo2soil   =',trcs_netflow2soil_slow_flx_col(idg,NY,NX)
              write(201,*)'netprod2soil  =',trcs_NetProd_slow_soil_col(idg,NY,NX)  
              write(201,*)'wetdep2soil   =',Gas_WeDepo2soil_slow_flx_col(idg,NY,NX)
            endif
            write(201,*)('-',L=1,50)
            write(201,*)'err litr      =',trcg_mass_litr_begs(idg,NY,NX)-trcg_mass_litr_now(idg)-trcs_drib_litr_slow_begf(idg,NY,NX)+ &
              trcs_solml_drib_vr(idg,0,NY,NX)+AtmGasDiff2Litr_slow_flx_col(idg,NY,NX)+Gas_WetDepo2litr_slow_flx_col(idg,NY,NX)  &
              - Gas_litr2Soil_slow_flx_col(idg,NY,NX)+trcs_sno2litr_slow_flx_col(idg,NY,NX)
            write(201,*)'drib beg/end  =',trcs_drib_litr_slow_begf(idg,NY,NX),trcs_solml_drib_vr(idg,0,NY,NX)
            write(201,*)'atm dif2litr  =',AtmGasDiff2Litr_slow_flx_col(idg,NY,NX)
            write(201,*)'wetdep2litr   =',Gas_WetDepo2litr_slow_flx_col(idg,NY,NX)
            write(201,*)'litr2soil     =',Gas_litr2Soil_slow_flx_col(idg,NY,NX)
            write(201,*)'sno2litr      =',trcs_sno2litr_slow_flx_col(idg,NY,NX)  
          
          if(AMIN1(abs(err),abs(err_rel))>1.e-5_r8)call endrun(trim(mod_filename)//' at line',__LINE__)          
        endif
      ENDDO

      DO K=1,jcplx
        DO idom=idom_beg,idom_end
          dmass=DOM_mass_now(idom,K)-DOM_mass_begs(idom,K,NY,NX)
          err=dmass-DOM_Hydroloss_slow_flx_col(idom,K,NY,NX)-DOM_NetProd_slow_flxM_col(idom,K,NY,NX) &
            -dom_dribble_col(idom,K)+dom_dribble_slow_beg_col(idom,K,NY,NX)
          if(abs(err)>1.e-5)then
            write(201,*)('-',L=1,50)
            write(201,*)(I*1000+J)*10+M,'iterm=',iterm,'idom=',idom,'K=',K,NY,NX,NU_col(NY,NX),NUM_col(NY,NX)
            write(201,*)'dom mass beg/end 0',DOM_mass_begs(idom,K,NY,NX),DOM_mass_now(idom,K)
            write(201,*)'dom mass beg/end 3', DOM_mass3_col(idom,K,NY,NX),DOM_mass4_col(idom,K,NY,NX)
            write(201,*)'err           =',err
            write(201,*)'dmass         =',dmass
            write(201,*)'hydroloss.    =',DOM_Hydroloss_slow_flx_col(idom,K,NY,NX)
            write(201,*)'netflow       =',TranspNetDOM_flx_col(idom,K,NY,NX),TranspNetDOM_flx2_col(idom,K,NY,NX)
            write(201,*)'surface runoff=',DOM_SurfRunoff_flx_col(idom,K,NY,NX)
            write(201,*)'subsurf runoff=',DOM_transpFlx_2DH(idom,K,NY,NX)
            write(201,*)'domdrain      =',DOM_draing_col(idom,K,NY,NX)
            write(201,*)DOM_MicP2_vr(idom,K,0,NY,NX),(DOM_MicP2_vr(idom,K,L,NY,NX),L=NU_col(NY,NX),NL_col(NY,NX))
            write(201,*)(DOM_MacP2_vr(idom,K,L,NY,NX),L=NU_col(NY,NX),NL_col(NY,NX))
            write(201,*)'netpro        =',DOM_NetProd_slow_flxM_col(idom,K,NY,NX)
            write(201,*)'dribble,b/e   =',dom_dribble_slow_beg_col(idom,K,NY,NX),dom_dribble_col(idom,K)

           if(abs(err)>1.e-8_r8)call endrun(trim(mod_filename)//' at line',__LINE__)          
          endif
        ENDDO
      ENDDO  
    ENDDO
  ENDDO  

  end  subroutine ExitMassCheck
!------------------------------------------------------------------------------------------
  subroutine SlowUpdateStateVarsM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)
  implicit none
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS
  real(r8),intent(in) :: dpscal(ids_beg:ids_end)
  real(r8),intent(in) :: dpscal_dom(idom_beg:idom_end)
  real(r8),optional, intent(inout) :: pscal1(ids_beg:ids_end)
  real(r8),optional, intent(inout) :: pscal1_dom(idom_beg:idom_end)

  if(present(pscal1) .and. present(pscal1_dom))then

    pscal1     = 1.000_r8
    call UpdateSnowTracersM(I,J,M,NHW,NHE,NVN,NVS,dpscal, pscal1)

    pscal1_dom = 1.000_r8
    call UpdateSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)

    call UpdateSubSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)


  else
    DOM_mass3_col=0._r8
    DOM_mass4_col=0._r8

    call UpdateSnowTracersM(I,J,M,NHW,NHE,NVN,NVS,dpscal)

    call UpdateSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom)

    call UpdateSubSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom)
  endif
  end subroutine SlowUpdateStateVarsM
!------------------------------------------------------------------------------------------

  subroutine copyStateVars(I,J,M,NHE,NHW,NVS,NVN,iterm)
  implicit none
  integer, intent(in) :: I,J,M,NHE,NHW,NVS,NVN,iterm
  character(len=*), parameter :: subname='copyStateVars'
  integer :: NY,NX,idg,ids,L,K,idom
  
  call PrintInfo('beg '//subname)
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      !to litter
      DO  K=1,micpar%NumOfLitrCmplxs
        DO idom=idom_beg,idom_end
          DOM_MicP2_vr(idom,K,0,NY,NX)=AZERO(DOM_MicP_vr(idom,K,0,NY,NX))
        ENDDO
      ENDDO

      !exclude NH3B
      DO idg=ids_beg,ids_end
        trcs_solml2_vr(idg,0,NY,NX)=AZERO(trcs_solml_vr(idg,0,NY,NX))
      ENDDO
      
      DO L=NU_col(NY,NX),NL_col(NY,NX)    
        DO ids=ids_beg,ids_end
          trcs_solml2_vr(ids,L,NY,NX) = AZERO(trcs_solml_vr(ids,L,NY,NX))
          trcs_soHml2_vr(ids,L,NY,NX) = AZERO(trcs_soHml_vr(ids,L,NY,NX))
          if(trcs_solml2_vr(ids,L,NY,NX)<0._r8)then
            write(125,*)(I*1000+J)*100+M,iterm,trcs_names(ids),L,trcs_solml2_vr(ids,L,NY,NX),trcs_soHml2_vr(ids,L,NY,NX)
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
        ENDDO

        DO  K=1,jcplx
          DO idom=idom_beg,idom_end
            DOM_MicP2_vr(idom,K,L,NY,NX)=AZERO(DOM_MicP_vr(idom,K,L,NY,NX))
            DOM_MacP2_vr(idom,K,L,NY,NX)=AZERO(DOM_MacP_vr(idom,K,L,NY,NX))
          ENDDO
        ENDDO
      ENDDO

      DO L=1,JS
        DO idg=idg_beg,idg_NH3
          trcg_solsml2_snvr(idg,L,NY,NX)  = AZERO(trcg_solsml_snvr(idg,L,NY,NX))
        ENDDO

        DO ids=ids_nut_beg,ids_nuts_end
          trcn_solsml2_snvr(ids,L,NY,NX) = AZERO(trcn_solsml_snvr(ids,L,NY,NX))
        ENDDO        
      enddo
    ENDDO
  ENDDO    
  call PrintInfo('end '//subname)
  end subroutine copyStateVars
!------------------------------------------------------------------------------------------

  subroutine UpdateSubSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)
  implicit none
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS
  real(r8), intent(in):: dpscal(ids_beg:ids_end)
  real(r8), intent(in):: dpscal_dom(idom_beg:idom_end)
  real(r8),optional,intent(inout) :: pscal1(ids_beg:ids_end)
  real(r8),optional,intent(inout) :: pscal1_dom(idom_beg:idom_end)

  character(len=*), parameter :: subname = 'UpdateSubSurfTracerM'
  integer :: NY,NX,L,K,idg,ids,idom
  real(r8) :: pscal(ids_beg:ids_end)
  real(r8) :: pscal_dom(idom_beg:idom_end),flux
  logical :: lflux,lfupdate
  
  call PrintInfo('beg '//subname)
  lfupdate=.false.
  if(present(pscal1))then
    pscal=pscal1
  else
    lfupdate=.true.
    pscal=1.000_r8
  endif

  if(present(pscal1_dom))then
    pscal_dom=pscal1_dom
  else
    pscal_dom=1.000_r8
    lfupdate=.true.
  endif

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      DO L=NU_col(NY,NX)+1,NL_col(NY,NX)    
        IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN        

          DO idg=idg_beg,idg_end     
            flux=trcsol_Irrig_flxM_vr(idg,L,NY,NX)+trcs_Transp2Micp_flxM_vr(idg,L,NY,NX)+trcs_Mac2MicPore_flxM_vr(idg,L,NY,NX)
            flux=flux*dpscal(idg)    
            lflux=.true. .or. .not.isclose(flux,0._r8) .and. pscal(idg)>0._r8                 
            if(lflux .and. dpscal(idg)>tiny_p)then          
              call get_flux_scalar(trcs_solml2_vr(idg,L,NY,NX),flux,trcs_solml_vr(idg,L,NY,NX),pscal(idg))

              if(trcs_solml_vr(idg,L,NY,NX)<0._r8)then
                write(141,*)(I*1000+J)*100+M,trcs_names(idg),trcs_solml2_vr(idg,L,NY,NX),trcs_solml_vr(idg,L,NY,NX),pscal(idg)
                call endrun(trim(mod_filename)//' at line',__LINE__)
              endif

              if(lfupdate)then
                TranspNetSoil_slow_flxM_col(idg,NY,NX)=TranspNetSoil_slow_flxM_col(idg,NY,NX)+ &
                  trcsol_Irrig_flxM_vr(idg,L,NY,NX)+trcs_Transp2Micp_flxM_vr(idg,L,NY,NX) +&
                  trcs_Transp2Macp_flxM_vr(idg,L,NY,NX)
                trcs_solml2_vr(idg,L,NY,NX)=trcs_solml_vr(idg,L,NY,NX)
              endif  

              if(trcs_solml_vr(idg,L,NY,NX)<0._r8)then
                write(132,*)(I*1000+J)*100+M,trcs_names(idg),trcs_solml_vr(idg,L,NY,NX),pscal(idg)
                call endrun(trim(mod_filename)//' at line',__LINE__)
              endif  
            endif
          ENDDO

          DO ids=idg_NH3B+1,ids_end            
            flux=trcsol_Irrig_flxM_vr(ids,L,NY,NX)+trcs_Transp2Micp_flxM_vr(ids,L,NY,NX) + &
              trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)
            lflux=.true. .or. .not.isclose(flux,0._r8)  
            if(lflux .and. dpscal(ids)>tiny_p)then          
              flux=flux*dpscal(ids)   
              call get_flux_scalar(trcs_solml2_vr(ids,L,NY,NX),flux,trcs_solml_vr(ids,L,NY,NX),pscal(ids))            
              if(lfupdate)then
                trcs_solml2_vr(ids,L,NY,NX)=trcs_solml_vr(ids,L,NY,NX)
              endif
            endif
          ENDDO

          DO  K=1,jcplx
            DO idom=idom_beg,idom_end
              flux=DOM_Transp2Micp_flxM_vr(idom,K,L,NY,NX)+DOM_Mac2MicPore_flxM_vr(idom,K,L,NY,NX)
              lflux=.true. .or. .not.isclose(flux,0._r8)
              if(lflux .and. dpscal_dom(idom)>tiny_p)then          
                flux=flux*dpscal_dom(idom)            
                call get_flux_scalar(DOM_MicP2_vr(idom,K,L,NY,NX), flux, DOM_MicP_vr(idom,K,L,NY,NX),pscal_dom(idom))           

                flux=DOM_Transp2Macp_flxM_vr(idom,K,L,NY,NX)-DOM_Mac2MicPore_flxM_vr(idom,K,L,NY,NX)
                flux=flux*dpscal_dom(idom)
                call get_flux_scalar(DOM_MacP2_vr(idom,K,L,NY,NX),flux, DOM_MacP_vr(idom,K,L,NY,NX),pscal_dom(idom))

                if(lfupdate)then
                  DOM_mass3_col(idom,K,NY,NX)=DOM_mass3_col(idom,K,NY,NX)+DOM_MicP2_vr(idom,K,L,NY,NX)+DOM_MacP2_vr(idom,K,L,NY,NX)
                  DOM_MicP2_vr(idom,K,L,NY,NX)= DOM_MicP_vr(idom,K,L,NY,NX)
                  DOM_MacP2_vr(idom,K,L,NY,NX)= DOM_MacP_vr(idom,K,L,NY,NX)
                  DOM_mass4_col(idom,K,NY,NX) = DOM_mass4_col(idom,K,NY,NX)+DOM_MicP2_vr(idom,K,L,NY,NX)+DOM_MacP2_vr(idom,K,L,NY,NX)                  
                  TranspNetDOM_flxM_col(idom,K,NY,NX)=TranspNetDOM_flxM_col(idom,K,NY,NX)+ &
                    DOM_Transp2Micp_flxM_vr(idom,K,L,NY,NX)+DOM_Transp2Macp_flxM_vr(idom,K,L,NY,NX)  

                endif                  
              endif
            ENDDO
          ENDDO

          DO ids=ids_beg,ids_end            
            flux=trcs_Transp2Macp_flxM_vr(ids,L,NY,NX) - trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)
            lflux=.true. .or. .not.isclose(flux,0._r8)
            if(lflux .and. dpscal(ids)>tiny_p)then          
              flux=flux*dpscal(ids)   
              call get_flux_scalar(trcs_soHml2_vr(ids,L,NY,NX),flux, trcs_soHml_vr(ids,L,NY,NX),pscal(ids))            
              if(lfupdate)then
                trcs_soHml2_vr(ids,L,NY,NX)=trcs_soHml_vr(ids,L,NY,NX)
              endif  
            endif
          enddo
        ENDIF
      ENDDO  
    ENDDO
  ENDDO  
  if(present(pscal1))pscal1=pscal
  if(present(pscal1_dom))pscal1_dom=pscal_dom
  call PrintInfo('end '//subname)
  end subroutine UpdateSubSurfTracerM
!------------------------------------------------------------------------------------------
  subroutine AccumSlowFluxesM(I,J,M,NHE,NHW,NVS,NVN,dpscal,dpcsal_dom,pscal1,pscal1_dom)
  !
  !DESCRIPTION
  !Accumulate fluxes for mass conservation check
  !Lateral runoff, dissolution (aka dry gas deposition), 
  !wet deposition through rain and irrigation, 
  !snow fluxes to litter and soil
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  real(r8),intent(in) :: dpscal(ids_beg:ids_end)
  real(r8), intent(in):: pscal1(ids_beg:ids_end)
  real(r8), intent(in):: dpcsal_dom(idom_beg:idom_end)
  real(r8), intent(in):: pscal1_dom(idom_beg:idom_end)

  character(len=*), parameter :: subname='AccumSlowFluxesM'
  integer :: idom,K,ids,NY,NX,idg,idn,L
  real(r8):: ppscal(ids_beg:ids_end),pscal_max,FLUX
  real(r8):: ppscal_dom(idom_beg:idom_end)

  call PrintInfo('beg '//subname)
  pscal_max=0._r8
  DO ids=ids_beg,ids_end
    ppscal(ids) = dpscal(ids)*pscal1(ids)
    pscal_max=AMAX1(pscal1(ids),pscal_max)
  ENDDO

  do idom=idom_beg,idom_end    
    ppscal_dom(idom)=dpcsal_dom(idom)*pscal1_dom(idom)
    pscal_max=AMAX1(pscal1_dom(idom),pscal_max)
  enddo

  if(pscal_max<1.0001_r8)then  
    call SlowUpdateStateVarsM(I,J,M,NHW,NHE,NVN,NVS,ppscal,ppscal_dom)  
  endif

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      DO idg=idg_beg,idg_NH3      
        if(ppscal(idg)>tiny_p)then

          flux=ppscal(idg)*(RGasAtmDisol2LitrM_col(idg,NY,NX)+RGasAtmDisol2SoilM_col(idg,NY,NX))

          AtmGasDiff2Surf_slow_flx_col(idg,NY,NX) = AtmGasDiff2Surf_slow_flx_col(idg,NY,NX)+flux

          flux=ppscal(idg)*(trcg_Precip2LitrM_col(idg,NY,NX) + trcs_Precip2MicpM_col(idg,NY,NX) & 
            +trcg_AquaAdv_flxM_snvr(idg,1,NY,NX))

          Gas_WetDeposit_slow_flx_col(idg,NY,NX)= Gas_WetDeposit_slow_flx_col(idg,NY,NX)+flux
          Gas_WeDepo2soil_slow_flx_col(idg,NY,NX)=Gas_WeDepo2soil_slow_flx_col(idg,NY,NX)+ppscal(idg)*trcs_Precip2MicpM_col(idg,NY,NX)          
          Gas_WetDepo2litr_slow_flx_col(idg,NY,NX)= Gas_WetDepo2litr_slow_flx_col(idg,NY,NX)+ppscal(idg)*trcg_Precip2LitrM_col(idg,NY,NX)

          flux=ppscal(idg) *trcg_AquaAdv_flxM_snvr(idg,1,NY,NX)
          Gas_WetDepo2Snow_slow_flx_col(idg,NY,NX)=Gas_WetDepo2Snow_slow_flx_col(idg,NY,NX)+flux

          Gas_WetDepo2Snow_col(idg,NY,NX) = Gas_WetDepo2Snow_col(idg,NY,NX)+ flux

          flux=ppscal(idg)*(trcg_AquaADV_Snow2Litr_flxM(idg,NY,NX)+trcg_AquaADV_Snow2Soil_flxM(idg,NY,NX) &
            -trcg_SnowDrift_flxM(idg,NY,NX))

          trcs_NetFlow2Litr_slow_flx_col(idg,NY,NX)=trcs_NetFlow2Litr_slow_flx_col(idg,NY,NX)+ ppscal(idg) &
            *(trcs_MicpTranspFlxM_3D(idg,3,0,NY,NX)+RGasAtmDisol2LitrM_col(idg,NY,NX)+trcg_SurfRunoff_flxM(idg,NY,NX))

          Gas_Snowloss_flx_col(idg,NY,NX)      = Gas_Snowloss_flx_col(idg,NY,NX)+ flux

          Gas_Snowloss_slow_flx_col(idg,NY,NX) = Gas_Snowloss_slow_flx_col(idg,NY,NX)+ flux
            
          AtmGasDiff2Litr_slow_flx_col(idg,NY,NX)=AtmGasDiff2Litr_slow_flx_col(idg,NY,NX)+ppscal(idg)*RGasAtmDisol2LitrM_col(idg,NY,NX)

          Gas_WetDepo2Litr_col(idg,NY,NX)  =Gas_WetDepo2Litr_col(idg,NY,NX)+ppscal(idg) &
            *(trcg_Precip2LitrM_col(idg,NY,NX))

          Gas_litr2Soil_slow_flx_col(idg,NY,NX) =Gas_litr2Soil_slow_flx_col(idg,NY,NX)+ ppscal(idg) &
            *Gas_litr2Soil_flxM_col(idg,NY,NX)  

          !atmospheric disolution into soil
          flux                                      = ppscal(idg)*RGasAtmDisol2SoilM_col(idg,NY,NX)
          AtmGasDiff2Soil_flx_col(idg,NY,NX)        = AtmGasDiff2Soil_flx_col(idg,NY,NX)+flux
          trcs_netflow2soil_slow_flx_col(idg,NY,NX) = trcs_netflow2soil_slow_flx_col(idg,NY,NX)+flux

          Gas_WetDepo2Soil_col(idg,NY,NX)  =Gas_WetDepo2Soil_col(idg,NY,NX)+ppscal(idg)*trcs_Precip2MicpM_col(idg,NY,NX)

          !---------------------------------------------------------------------------
          L=NU_col(NY,NX)

          flux=ppscal(idg)*(RGasAtmDisol2SoilM_col(idg,NY,NX) &             
            +trcs_MicpTranspFlxM_3D(idg,3,L,NY,NX) +trcs_MacpTranspFlxM_3D(idg,3,L,NY,NX))

          trc_topsoil_flx_col(idg,NY,NX) =trc_topsoil_flx_col(idg,NY,NX)+flux

          TranspNetSoil_flx_col(idg,NY,NX)=TranspNetSoil_flx_col(idg,NY,NX)+ppscal(idg) &
            *(RGasAtmDisol2SoilM_col(idg,NY,NX)) 

          TranspNetSoil_flx2_col(idg,NY,NX)=TranspNetSoil_flx2_col(idg,NY,NX) + ppscal(idg) &
            * TranspNetSoil_slow_flxM_col(idg,NY,NX)

          DO L=NU_col(NY,NX),NL_col(NY,NX)
     
            TranspNetSoil_flx_col(idg,NY,NX)=TranspNetSoil_flx_col(idg,NY,NX)+ppscal(idg) &
              *(trcsol_Irrig_flxM_vr(idg,L,NY,NX)+trcs_Transp2Micp_flxM_vr(idg,L,NY,NX) &
              +trcs_Transp2Macp_flxM_vr(idg,L,NY,NX)) 

            transp_diff_slow_vr(idg,L,NY,NX) = transp_diff_slow_vr(idg,L,NY,NX) +  &
              (trcs_Transp2Micp_flxM_vr(idg,L,NY,NX)+trcs_Transp2Macp_flxM_vr(idg,L,NY,NX) &
              -trcs_MicpTranspFlxM_3D(idg,3,L,NY,NX) -trcs_MacpTranspFlxM_3D(idg,3,L,NY,NX)  &
              +trcs_MicpTranspFlxM_3D(idg,3,L+1,NY,NX) +trcs_MacpTranspFlxM_3D(idg,3,L+1,NY,NX))

            flux=ppscal(idg)*trcsol_Irrig_flxM_vr(idg,L,NY,NX)  
            trc_topsoil_flx_col(idg,NY,NX) =trc_topsoil_flx_col(idg,NY,NX)+flux

          ENDDO       
        endif
      ENDDO

      idg=idg_NH3B
      if(ppscal(idg)>tiny_p)then     
        flux                                     = trcn_AquaADV_Snow2Band_flxM(idg,NY,NX)*ppscal(idg)
        Gas_Snowloss_slow_flx_col(idg_NH3,NY,NX) = Gas_Snowloss_slow_flx_col(idg_NH3,NY,NX)+flux
        Gas_Snowloss_flx_col(idg_NH3,NY,NX)      = Gas_Snowloss_flx_col(idg_NH3,NY,NX)+flux

        flux                                          = ppscal(idg)*Gas_litr2Soil_flxM_col(idg,NY,NX)
        trcs_NetFlow2Litr_slow_flx_col(idg_NH3,NY,NX) = trcs_NetFlow2Litr_slow_flx_col(idg_NH3,NY,NX)-flux
        Gas_litr2Soil_slow_flx_col(idg_NH3,NY,NX)     = Gas_litr2Soil_slow_flx_col(idg_NH3,NY,NX)+ flux

        DO L=NU_col(NY,NX),NL_col(NY,NX)

          flux=ppscal(idg)*trcsol_Irrig_flxM_vr(idg,L,NY,NX)    
          trc_topsoil_flx_col(idg_NH3,NY,NX) =trc_topsoil_flx_col(idg_NH3,NY,NX)+flux

          TranspNetSoil_flx_col(idg_NH3,NY,NX)=TranspNetSoil_flx_col(idg_NH3,NY,NX)+ppscal(idg) &
            *(trcs_Transp2Micp_flxM_vr(idg,L,NY,NX)+trcs_Transp2Macp_flxM_vr(idg,L,NY,NX)) 
        ENDDO

        TranspNetSoil_flx2_col(idg_NH3,NY,NX)=TranspNetSoil_flx2_col(idg_NH3,NY,NX) + ppscal(idg) &
          * TranspNetSoil_slow_flxM_col(idg,NY,NX)          

        flux=ppscal(idg)*trcs_Precip2MicpM_col(idg,NY,NX)

        Gas_WetDeposit_slow_flx_col(idg_NH3,NY,NX)= Gas_WetDeposit_slow_flx_col(idg_NH3,NY,NX)+flux

        !atmospheric disolution into topsoil
        flux=ppscal(idg)*RGasAtmDisol2SoilM_col(idg,NY,NX)

        trcs_netflow2soil_slow_flx_col(idg,NY,NX)=trcs_netflow2soil_slow_flx_col(idg,NY,NX)+flux

        AtmGasDiff2Surf_slow_flx_col(idg_NH3,NY,NX) = AtmGasDiff2Surf_slow_flx_col(idg_NH3,NY,NX)+flux

        AtmGasDiff2Soil_flx_col(idg_NH3,NY,NX) = AtmGasDiff2Soil_flx_col(idg_NH3,NY,NX)+flux

        TranspNetSoil_flx_col(idg_NH3,NY,NX)=TranspNetSoil_flx_col(idg_NH3,NY,NX)+flux

        trc_topsoil_flx_col(idg_NH3,NY,NX) =trc_topsoil_flx_col(idg_NH3,NY,NX)+flux        

        L=NU_col(NY,NX)
        flux=ppscal(idg)*(trcs_MicpTranspFlxM_3D(idg,3,L,NY,NX) +trcs_MacpTranspFlxM_3D(idg,3,L,NY,NX))
        trc_topsoil_flx_col(idg_NH3,NY,NX) =trc_topsoil_flx_col(idg_NH3,NY,NX)+flux

      endif      

      RCH4PhysexchPrev_vr(0,NY,NX)= RCH4PhysexchPrev_vr(0,NY,NX) + ppscal(idg_CH4)*(RGasAtmDisol2LitrM_col(idg_CH4,NY,NX) &
        -trcg_Precip2LitrM_col(idg_CH4,NY,NX))

      RO2AquaSourcePrev_vr(0,NY,NX)   =RO2AquaSourcePrev_vr(0,NY,NX)+ppscal(idg_O2)*(RGasAtmDisol2LitrM_col(idg_O2,NY,NX) &
        -trcg_Precip2LitrM_col(idg_O2,NY,NX)) 

      RCH4PhysexchPrev_vr(NU_col(NY,NX),NY,NX)= RCH4PhysexchPrev_vr(NU_col(NY,NX),NY,NX)+ppscal(idg_CH4)*(RGasAtmDisol2SoilM_col(idg_CH4,NY,NX) &
        - trcs_Precip2MicpM_col(idg_CH4,NY,NX))

      RO2AquaSourcePrev_vr(NU_col(NY,NX),NY,NX)   =RO2AquaSourcePrev_vr(NU_col(NY,NX),NY,NX)+ppscal(idg_O2)*(RGasAtmDisol2SoilM_col(idg_O2,NY,NX) &
        - trcs_Precip2MicpM_col(idg_O2,NY,NX))
            
      do idg=idg_beg,idg_NH3
        if(ppscal(idg)>tiny_p)then      
          trcs_sno2litr_slow_flx_col(idg,NY,NX) = trcs_sno2litr_slow_flx_col(idg,NY,NX)+ppscal(idg)*trcg_AquaADV_Snow2Litr_flxM(idg,NY,NX)
        endif 
      ENDDO

      DO idg=idg_beg,idg_NH3
        if(ppscal(idg)>tiny_p)then
          trcg_AquaADV_Snow2Soil_flx(idg,NY,NX) = trcg_AquaADV_Snow2Soil_flx(idg,NY,NX)+trcg_AquaADV_Snow2Soil_flxM(idg,NY,NX)*ppscal(idg)
        endif
      enddo
      
      do ids=ids_nut_beg,ids_nuts_end
        if(ppscal(ids)>tiny_p)then
          trcn_AquaADV_Snow2Litr_flx(ids,NY,NX) = trcn_AquaADV_Snow2Litr_flx(ids,NY,NX)+trcn_AquaADV_Snow2Litr_flxM(ids,NY,NX)*ppscal(ids)
          trcn_AquaADV_Snow2Soil_flx(ids,NY,NX) = trcn_AquaADV_Snow2Soil_flx(ids,NY,NX)+trcn_AquaADV_Snow2Soil_flxM(ids,NY,NX)*ppscal(ids)
        endif
      enddo

      !from NH3B to ids_H2PO4B
      DO idn=ids_nutb_beg,ids_nutb_end
        if(ppscal(idn)>tiny_p)then
          trcn_AquaADV_Snow2Band_flx(idn,NY,NX) = trcn_AquaADV_Snow2Band_flx(idn,NY,NX)+trcn_AquaADV_Snow2Band_flxM(idn,NY,NX)*ppscal(idn)
        endif
      ENDDO

      L=NU_col(NY,NX)
      DO ids=ids_beg,ids_end        
        if(ppscal(ids)>tiny_p)then
          trcs_netflow2soil_slow_flx_col(ids,NY,NX)=trcs_netflow2soil_slow_flx_col(ids,NY,NX)+ppscal(ids) &
            *(trcs_MicpTranspFlxM_3D(ids,3,L,NY,NX) +trcs_MacpTranspFlxM_3D(ids,3,L,NY,NX))
        ENDIF
      ENDDO

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO ids=ids_beg,ids_end
          if(ppscal(ids)>tiny_p)then
            flux=ppscal(ids)*trcsol_Irrig_flxM_vr(ids,L,NY,NX)
            trcs_NetProd_slow_flxM_col(ids,NY,NX)      = trcs_NetProd_slow_flxM_col(ids,NY,NX)+flux
            trcs_NetProd_slow_soil_col(ids,NY,NX)      = trcs_NetProd_slow_soil_col(ids,NY,NX)+flux
          endif  
        ENDDO      
      ENDDO

      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          if(ppscal_dom(idom)>tiny_p)then

            TranspNetDOM_flx_col(idom,K,NY,NX)=TranspNetDOM_flx_col(idom,K,NY,NX)+TranspNetDOM_flxM_col(idom,K,NY,NX)*ppscal_dom(idom)

            FLUX=DOM_SurfRunoff_flxM(idom,K,NY,NX)*ppscal_dom(idom)
            DOM_Hydroloss_slow_flx_col(idom,K,NY,NX)  = DOM_Hydroloss_slow_flx_col(idom,K,NY,NX)+FLUX
            DOM_SurfRunoff_flx_col(idom,K,NY,NX)=DOM_SurfRunoff_flx_col(idom,K,NY,NX)+FLUX

            FLUX=DOM_transpFlxM_2DH(idom,K,NY,NX)*ppscal_dom(idom)
            DOM_Hydroloss_slow_flx_col(idom,K,NY,NX)  = DOM_Hydroloss_slow_flx_col(idom,K,NY,NX)+FLUX
            DOM_transpFlx_2DH(idom,K,NY,NX) = DOM_transpFlx_2DH(idom,K,NY,NX)+FLUX

            FLUX=ppscal_dom(idom)*(DOM_MicpTranspFlxM_3D(idom,K,3,NL_col(NY,NX)+1,NY,NX)+DOM_MacpTranspFlxM_3D(idom,K,3,NL_col(NY,NX)+1,NY,NX))
            DOM_draing_col(idom,K,NY,NX)             = DOM_draing_col(idom,K,NY,NX)+ FLUX
            DOM_Hydroloss_slow_flx_col(idom,K,NY,NX) = DOM_Hydroloss_slow_flx_col(idom,K,NY,NX)- FLUX

            TranspNetDOM_flx2_col(idom,K,NY,NX)=TranspNetDOM_flx2_col(idom,K,NY,NX)+ &
              (DOM_SurfRunoff_flxM(idom,K,NY,NX)+DOM_MicpTranspFlxM_3D(idom,K,3,0,NY,NX))*ppscal_dom(idom) 

            DO L=NU_col(NY,NX),NL_col(NY,NX)
              TranspNetDOM_flx2_col(idom,K,NY,NX)=TranspNetDOM_flx2_col(idom,K,NY,NX)+ &
                (DOM_Transp2Micp_flxM_vr(idom,K,L,NY,NX)+ DOM_Transp2Macp_flxM_vr(idom,K,L,NY,NX))*ppscal_dom(idom) 
            ENDDO  
            
          endif  
        ENDDO
      ENDDO

      DO ids=ids_beg,ids_end
        if(ppscal(ids)>tiny_p)then
          !lateral subsurface loss
          flux                                      = trcs_transpFlxM_2DH(ids,NY,NX)*ppscal(ids)
          trcs_netflow2soil_slow_flx_col(ids,NY,NX) = trcs_netflow2soil_slow_flx_col(ids,NY,NX)+flux
          trcs_SubsurTransp_flx_2DH(ids,NY,NX)      = trcs_SubsurTransp_flx_2DH(ids,NY,NX)+flux
          trcs_hydrloss_slow_flx_col(ids,NY,NX)     = trcs_hydrloss_slow_flx_col(ids,NY,NX)+flux

          !drainage loss at the bottom
          flux=ppscal(ids)*(trcs_MicpTranspFlxM_3D(ids,3,NL_col(NY,NX)+1,NY,NX)+trcs_MacpTranspFlxM_3D(ids,3,NL_col(NY,NX)+1,NY,NX))
          trcs_netflow2soil_slow_flx_col(ids,NY,NX) = trcs_netflow2soil_slow_flx_col(ids,NY,NX)-flux
          trcs_drainage_slow_flx_col(ids,NY,NX)     = trcs_drainage_slow_flx_col(ids,NY,NX) + flux
          trcs_hydrloss_slow_flx_col(ids,NY,NX)     = trcs_hydrloss_slow_flx_col(ids,NY,NX)-flux

        endif  
      ENDDO

      DO idg=idg_beg,idg_NH3
        if(ppscal(idg)>tiny_p)then
          trcg_SnowDrift_flx_col(idg,NY,NX)=trcg_SnowDrift_flx_col(idg,NY,NX)+trcg_SnowDrift_flxM(idg,NY,NX)*ppscal(idg)
        
          !surface loss thru runoff and snowdrift
          flux=ppscal(idg)*(trcg_SurfRunoff_flxM(idg,NY,NX)+trcg_SnowDrift_flxM(idg,NY,NX))
          GasHydroSurfLoss_flx_col(idg,NY,NX)   = GasHydroSurfLoss_flx_col(idg,NY,NX)+flux
          trcs_hydrloss_slow_flx_col(idg,NY,NX) = trcs_hydrloss_slow_flx_col(idg,NY,NX)+flux

          !lateral plus drainage 
          flux = ppscal(idg)*(trcs_transpFlxM_2DH(idg,NY,NX)  &
            -trcs_MicpTranspFlxM_3D(idg,3,NL_col(NY,NX)+1,NY,NX)-trcs_MacpTranspFlxM_3D(idg,3,NL_col(NY,NX)+1,NY,NX))

          GasHydroSubsLoss_flx_col(idg,NY,NX)=GasHydroSubsLoss_flx_col(idg,NY,NX)+flux

          GasHydroLoss_litr_slow_flx_col(idg,NY,NX)=GasHydroLoss_litr_slow_flx_col(idg,NY,NX)+ppscal(idg)  & 
            *trcg_SurfRunoff_flxM(idg,NY,NX)

        endif  
      ENDDO  

      !add NH3B to NH3
      if(ppscal(idg_NH3B)>tiny_p)then        
        flux=ppscal(idg_NH3B)*(trcs_transpFlxM_2DH(idg_NH3B,NY,NX)  &
          -trcs_MicpTranspFlxM_3D(idg_NH3B,3,NL_col(NY,NX)+1,NY,NX)-trcs_MacpTranspFlxM_3D(idg_NH3B,3,NL_col(NY,NX)+1,NY,NX))

        GasHydroSubsLoss_flx_col(idg_NH3,NY,NX) = GasHydroSubsLoss_flx_col(idg_NH3,NY,NX)+flux

      endif

      DO idn=ids_nut_beg,ids_nuts_end
        if(ppscal(idn)>tiny_p)then
          trcn_SnowDrift_flx_col(idn,NY,NX) =trcn_SnowDrift_flx_col(idn,NY,NX) + trcn_SnowDrift_flxM(idn,NY,NX)*ppscal(idn)
          trcn_SurfRunoff_flx_col(idn,NY,NX)=trcn_SurfRunoff_flx_col(idn,NY,NX)+ trcn_SurfRunoff_flxM(idn,NY,NX)*ppscal(idn)
        endif
      ENDDO  

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        RO2AquaSourcePrev_vr(L,NY,NX)   =RO2AquaSourcePrev_vr(L,NY,NX) +(trcs_Transp2Micp_flxM_vr(idg_O2,L,NY,NX) + &
          trcs_Mac2MicPore_flxM_vr(idg_O2,L,NY,NX)+trcsol_Irrig_flxM_vr(idg_O2,L,NY,NX))*ppscal(idg_O2)
        
        RCH4PhysexchPrev_vr(L,NY,NX) = RCH4PhysexchPrev_vr(L,NY,NX)+(trcs_Transp2Micp_flxM_vr(idg_CH4,L,NY,NX) + &
          trcs_Mac2MicPore_flxM_vr(idg_CH4,L,NY,NX)+trcsol_Irrig_flxM_vr(idg_CH4,L,NY,NX))*ppscal(idg_CH4)

        DO ids=ids_beg,ids_end
          if(ppscal(ids)>tiny_p)then
            trcs_irrig_flx_col(ids,NY,NX)=trcs_irrig_flx_col(ids,NY,NX)+ppscal(ids)*trcsol_Irrig_flxM_vr(ids,L,NY,NX)
          endif
        enddo          
      ENDDO
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine AccumSlowFluxesM
!------------------------------------------------------------------------------------------
  subroutine StageSlowTranspIterationM(I,J,M,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='StageSlowTranspIterationM'
  integer :: NY,NX,L, K

  call PrintInfo('beg '//subname)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        VLWatMicPMA_vr(L,NY,NX) = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
        VLWatMicPMB_vr(L,NY,NX) = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
        VLWatMicPXA_vr(L,NY,NX) = natomw*VLWatMicPMA_vr(L,NY,NX)
        VLWatMicPXB_vr(L,NY,NX) = natomw*VLWatMicPMB_vr(L,NY,NX)
        VLsoiAirPMA_vr(L,NY,NX) = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
        VLsoiAirPMB_vr(L,NY,NX) = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      ENDDO
    ENDDO  
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine StageSlowTranspIterationM
!------------------------------------------------------------------------------------------
  subroutine ZeroFluxesM(I,J,M,NHE,NHW,NVS,NVN)  

  implicit none
  integer, intent(in) :: I,J,M,NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='ZeroFluxesM'
  integer :: L,NY,NX,K

  call PrintInfo('beg '//subname)

  TranspNetDOM_flxM_col       = 0._r8
  TranspNetSoil_slow_flxM_col = 0._r8
  DOM_FloXSurRof_flxM_2DH     = 0.0_r8
  trcg_SurRof_flxM_2DH        = 0.0_r8
  trcn_SurRof_flxM_2DH        = 0.0_r8

  trcg_SnowDrift_flxM_2DH = 0.0_r8
  trcn_SnowDrift_flxM_2DH = 0.0_r8

  DOM_MacpTranspFlxM_3D   = 0._r8
  DOM_MicpTranspFlxM_3D   = 0._r8
  DOM_Mac2MicPore_flxM_vr = 0._r8

  trcs_MacpTranspFlxM_3D   = 0._r8
  trcs_MicpTranspFlxM_3D   = 0._r8
  trcs_Mac2MicPore_flxM_vr = 0._r8

  DOM_transpFlxM_2DH  = 0._r8
  trcs_transpFlxM_2DH = 0._r8

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      DO  K=1,jcplx
        DOM_SurfRunoff_flxM(idom_beg:idom_end,K,NY,NX)=0.0_r8
      ENDDO

      !initialize with zeros
      trcg_AquaADV_Snow2Litr_flxM(idg_beg:idg_NH3,NY,NX)           = 0.0_r8
      trcn_AquaADV_Snow2Litr_flxM(ids_nut_beg:ids_nuts_end,NY,NX)  = 0.0_r8
      trcg_AquaADV_Snow2Soil_flxM(idg_beg:idg_NH3,NY,NX)           = 0.0_r8
      trcn_AquaADV_Snow2Soil_flxM(ids_nut_beg:ids_nuts_end,NY,NX)  = 0.0_r8
      trcn_AquaADV_Snow2Band_flxM(ids_nutb_beg:ids_nutb_end,NY,NX) = 0.0_r8


      trcg_AquaAdv_flxM_snvr(idg_beg:idg_NH3,2:,NY,NX)            = 0.0_r8
      trcn_AquaAdv_flxM_snvr(ids_nut_beg:ids_nuts_end,2:JS,NY,NX) = 0.0_r8

      trcg_SurfRunoff_flxM(idg_beg:idg_NH3,NY,NX)          = 0.0_r8
      trcn_SurfRunoff_flxM(ids_nut_beg:ids_nuts_end,NY,NX) = 0.0_r8
      trcg_SnowDrift_flxM(idg_beg:idg_NH3,NY,NX)           = 0.0_r8
      trcn_SnowDrift_flxM(ids_nut_beg:ids_nuts_end,NY,NX)  = 0.0_r8

      !   INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
      DO  L=1,JS
        trcg_Aqua_flxM_snvr(idg_beg:idg_NH3,L,NY,NX)          = 0._r8
        trcn_Aqua_flxM_snvr(ids_nut_beg:ids_nuts_end,L,NY,NX) = 0._r8
      ENDDO

      DO L=NU_col(NY,NX),NL_col(NY,NX)    
        DO  K=1,jcplx
          DOM_Transp2Micp_flxM_vr(idom_beg:idom_end,K,L,NY,NX) = 0.0_r8
          DOM_Transp2Macp_flxM_vr(idom_beg:idom_end,K,L,NY,NX)= 0.0_r8
        ENDDO
        trcs_Transp2Micp_flxM_vr(ids_beg:ids_end,L,NY,NX) = 0.0_r8
        trcs_Transp2Macp_flxM_vr(ids_beg:ids_end,L,NY,NX) = 0._r8
      ENDDO  
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine ZeroFluxesM
!------------------------------------------------------------------------------------------
  subroutine TracerXInterGridsM(I,J,M,NHE,NHW,NVS,NVN)
  !
  !Description
  !Tracer exchange across the internal grids
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHE,NHW,NVS,NVN  
  character(len=*), parameter :: subname='TracerXInterGridsM'  
  integer :: NY,NX,L,N
  integer :: N3,N2,N1,LL
  integer :: N6,N5,N4
  LOGICAL :: LBOUND

  call PrintInfo('beg '//subname)
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      D125: DO L=1,NL_col(NY,NX)       
        N1=NX;N2=NY;N3=L    !source grid, avoid edges

        IF(N3.GE.NUM_col(N2,N1))call MicMacPoresSoluteExchM(M,N1,N2,N3)

        !identify the dest grid, which can be at the edges.                      
        D120: DO N=FlowDirIndicator_col(N2,N1),3
          IF(N.EQ.iWestEastDirection)THEN            
            IF(NX.EQ.NHE)CYCLE  !skip eastern boundary
            N4 = NX+1;N5 = NY;N6 = L
          ELSEIF(N.EQ.iNorthSouthDirection)THEN                        
            IF(NY.EQ.NVS)CYCLE  !skip southern boundary
            N4 = NX;N5 = NY+1;N6 = L
          ELSEIF(N.EQ.iVerticalDirection)THEN                        
            IF(L.EQ.NL_col(NY,NX))CYCLE  !skip bottom boundary
            N4 = NX;N5 = NY;N6 = L+1
          ENDIF
          
          DO LL=N6,NL_col(N5,N4)
            IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
              N6=LL
              exit
            ENDIF
          ENDDO
          if(LL>NL_col(N5,N4))cycle D120
          
          IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
            IF(N3.GE.NUM_col(N2,N1) .AND. N6.GE.NUM_col(N5,N4))THEN
              !
              call SoluteAdvDifusTranspM(I,J,M,N,N1,N2,N3,N4,N5,N6)
              !
            ENDIF                
          ENDIF
        ENDDO D120  
      ENDDO D125
    ENDDO
  ENDDO  
  call PrintInfo('end '//subname)    
  end subroutine TracerXInterGridsM


!------------------------------------------------------------------------------------------
  subroutine GatherTranspFluxM(I,J,M,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='GatherTranspFluxM'
  integer :: NY,NX,L,LL,N
  integer :: N6,N5,N4,N3,N2,N1

  call PrintInfo('beg '//subname)
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      DO L=NU_col(NY,NX),NL_col(NY,NX)

        N1=NX;N2=NY;N3=L          
        DO  N=FlowDirIndicator_col(NY,NX),3
          IF(N.EQ.iWestEastDirection)THEN
            !WEST-EAST
            N4 = NX+1; N5 = NY ;N6 = L
          ELSEIF(N.EQ.iNorthSouthDirection)THEN              
            N4  = NX;N5  = NY+1;N6  = L
          ELSEIF(N.EQ.iVerticalDirection)THEN !vertical
            N4 = NX;N5 = NY;N6 = L+1  !target
          ENDIF                      

          IF(FlowDirIndicator_col(N2,N1).NE.3 .OR. N.EQ.iVerticalDirection)THEN
            DO LL=N6,NL_col(N5,N4)
              IF(VLSoilPoreMicP_vr(LL,N5,N4).GT.ZEROS2(N5,N4))THEN
                N6=LL
                exit
              ENDIF
            ENDDO
            !grids (N3,N2,N1) and (N6,N5,N4) are not necessary at the boundary
            call NetTracerFlowXSoilPoresM(I,J,M,N,N1,N2,N3,N4,N5,N6)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO        
  call PrintInfo('end '//subname)
  end subroutine GatherTranspFluxM
!------------------------------------------------------------------------------------------

  subroutine SurfaceSoluteFluxM(I,J,M,NHE,NHW,NVS,NVN,band_snow_donor_col)
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  ! Present only in the once-per-substep ID-changing surface-transfer pass.
  integer, optional, intent(in) :: band_snow_donor_col(JY,JX)
  character(len=*), parameter :: subname='SurfaceSoluteFluxM'

  integer  :: NY,NX,K
  real(r8) :: FLWRM1,snow_soil_water,snow_mac_frac
  real(r8) :: trcs_cl_litr(ids_beg:ids_end)
  real(r8) :: trcs_cl_soil(ids_beg:ids_end)
  real(r8) :: solute_adv_Lit2Soil_flxM(ids_beg:ids_end)
  real(r8) :: trcs_Dif_Litr2Soil_flxM(ids_beg:ids_end)
  real(r8) :: DOM_Adv_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: DOM_Difus_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: solute_adv_Lit2MacP_flxM(ids_beg:ids_end)
  real(r8) :: DOM_Adv_Litr2MacP_flxM(idom_beg:idom_end,1:jcplx)


  call PrintInfo('beg '//subname)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      !-----------------------------------------------------------------
      !Vertical 1D exchange
      call AtmosLitterExchangeM(I,J,M,NY,NX,trcs_cl_litr)
      !
      call AtmosSoilExchangeM(I,J,M,NY,NX,trcs_cl_soil)
      !
      call Lit2SoilTracerAdvectM(I,J,M,NY,NX,FLWRM1,solute_adv_Lit2Soil_flxM,DOM_Adv_Litr2Soil_flxM,&
        solute_adv_Lit2MacP_flxM,DOM_Adv_Litr2MacP_flxM)
      !
      call LitterSoilTracerDiffusionM(I,J,M,NY,NX,FLWRM1,trcs_cl_litr,trcs_cl_soil,&
        trcs_Dif_Litr2Soil_flxM,DOM_Difus_Litr2Soil_flxM)

      ! Keep total snow discharge for donor accounting; split soil credits
      ! with the liquid-water destinations. Zero flow retains the old fallback.
      snow_soil_water=WatFlowSno2MicPM_col(M,NY,NX)+WatFlowSno2MacPM_col(M,NY,NX)
      snow_mac_frac=0._r8
      if(snow_soil_water>0._r8)snow_mac_frac=WatFlowSno2MacPM_col(M,NY,NX)/snow_soil_water
      !-----------------------------------------------------------------        
      !
      if(present(band_snow_donor_col))then
        call AcceptSurfaceBandTransfersM(NY,NX,band_snow_donor_col(NY,NX), &
          solute_adv_Lit2Soil_flxM+trcs_Dif_Litr2Soil_flxM,solute_adv_Lit2MacP_flxM,snow_mac_frac)
      else
        ! These transfers were already accepted; do not replay them in retries.
        solute_adv_Lit2Soil_flxM(ids_nutb_beg:ids_nutb_end)=0._r8
        trcs_Dif_Litr2Soil_flxM(ids_nutb_beg:ids_nutb_end)=0._r8
        solute_adv_Lit2MacP_flxM(ids_nutb_beg:ids_nutb_end)=0._r8
        call SurfLayerNetTracerFluxM(I,J,M,NY,NX,trcs_Dif_Litr2Soil_flxM,solute_adv_Lit2Soil_flxM,&
          DOM_Adv_Litr2Soil_flxM,DOM_Difus_Litr2Soil_flxM,solute_adv_Lit2MacP_flxM,DOM_Adv_Litr2MacP_flxM,snow_mac_frac)
      endif
    ENDDO
  ENDDO

  if(.not.present(band_snow_donor_col))then
    call TracerFlowThruSurRofM(I,J,M,NHE,NHW,NVS,NVN)
    call TracerFlowThruSnowRedist(I,J,M,NHE,NHW,NVS,NVN)
  endif

  call PrintInfo('end '//subname)
  end subroutine SurfaceSoluteFluxM

!------------------------------------------------------------------------------------------

  subroutine AcceptSurfaceBandTransfersM(NY,NX,LS,litter_to_micp,litter_to_macp,snow_mac_frac)
  ! Accept only transfers that change IDs between litter/snow and banded soil.
  ! All other transport uses the existing independent per-tracer retry solver.
  implicit none
  integer, intent(in) :: NY,NX,LS ! LS=0 means no bottom snow donor
  real(r8), intent(in) :: snow_mac_frac ! fraction of soil-directed snow liquid entering macropores
  real(r8), intent(in) :: litter_to_micp(ids_beg:ids_end) ! signed litter-to-band proposal
  real(r8), intent(in) :: litter_to_macp(ids_beg:ids_end) ! nonnegative litter-to-band proposal
  character(len=*), parameter :: subname='AcceptSurfaceBandTransfersM'
  integer :: ids,idsband,L
  real(r8) :: litter_mass,soil_mass,snow_mass,outflow,scale
  real(r8) :: micp_flux,macp_flux,snow_flux,soil_gain

  call PrintInfo('beg '//subname)
  L=NU_col(NY,NX)
  DO idsband=ids_nutb_beg,ids_nutb_end
    if(idsband==idg_NH3B)then
      ids=idg_NH3
    else
      ids=idsband-ids_NH4B+ids_NH4
    endif
    litter_mass=AZMAX1(trcs_solml2_vr(ids,0,NY,NX))
    soil_mass=AZMAX1(trcs_solml2_vr(idsband,L,NY,NX))
    ! Both litter destinations draw on the same donor snapshot. Incoming
    ! soil/snow transfers are not reused to fund outgoing transfers here.
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
      if(ids==idg_NH3)then
        snow_mass=AZMAX1(trcg_solsml2_snvr(ids,LS,NY,NX))
      else
        snow_mass=AZMAX1(trcn_solsml2_snvr(ids,LS,NY,NX))
      endif
      snow_flux=MIN(snow_mass,AZMAX1(trcn_AquaADV_Snow2Band_flxM(idsband,NY,NX)))
      if(ids==idg_NH3)then
        trcg_solsml2_snvr(ids,LS,NY,NX)=snow_mass-snow_flux
        trcg_solsml_snvr(ids,LS,NY,NX)=trcg_solsml2_snvr(ids,LS,NY,NX)
      else
        trcn_solsml2_snvr(ids,LS,NY,NX)=snow_mass-snow_flux
        trcn_solsml_snvr(ids,LS,NY,NX)=trcn_solsml2_snvr(ids,LS,NY,NX)
      endif
    endif

    trcs_solml2_vr(ids,0,NY,NX)=litter_mass-micp_flux-macp_flux
    trcs_solml2_vr(idsband,L,NY,NX)=soil_mass+micp_flux+snow_flux*(1._r8-snow_mac_frac)
    trcs_soHml2_vr(idsband,L,NY,NX)=trcs_soHml2_vr(idsband,L,NY,NX)+macp_flux+snow_flux*snow_mac_frac
    trcs_solml_vr(ids,0,NY,NX)=trcs_solml2_vr(ids,0,NY,NX)
    trcs_solml_vr(idsband,L,NY,NX)=trcs_solml2_vr(idsband,L,NY,NX)
    trcs_soHml_vr(idsband,L,NY,NX)=trcs_soHml2_vr(idsband,L,NY,NX)

    ! Record the accepted internal transfers once, outside ordinary retries.
    soil_gain=micp_flux+macp_flux+snow_flux
    trcs_netflow2soil_slow_flx_col(idsband,NY,NX)=trcs_netflow2soil_slow_flx_col(idsband,NY,NX)+soil_gain
    trcn_AquaADV_Snow2Band_flx(idsband,NY,NX)=trcn_AquaADV_Snow2Band_flx(idsband,NY,NX)+snow_flux
    if(ids==idg_NH3)then
      trcs_NetFlow2Litr_slow_flx_col(ids,NY,NX)=trcs_NetFlow2Litr_slow_flx_col(ids,NY,NX)-micp_flux
      Gas_litr2Soil_slow_flx_col(ids,NY,NX)=Gas_litr2Soil_slow_flx_col(ids,NY,NX)+micp_flux
      Gas_Snowloss_slow_flx_col(ids,NY,NX)=Gas_Snowloss_slow_flx_col(ids,NY,NX)+snow_flux
      Gas_Snowloss_flx_col(ids,NY,NX)=Gas_Snowloss_flx_col(ids,NY,NX)+snow_flux
      trc_topsoil_flx_col(ids,NY,NX)=trc_topsoil_flx_col(ids,NY,NX)+soil_gain
      TranspNetSoil_flx_col(ids,NY,NX)=TranspNetSoil_flx_col(ids,NY,NX)+soil_gain
      TranspNetSoil_flx2_col(ids,NY,NX)=TranspNetSoil_flx2_col(ids,NY,NX)+soil_gain
    endif
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine AcceptSurfaceBandTransfersM
!------------------------------------------------------------------------------------------

  subroutine AtmosLitterExchangeM(I,J,M,NY,NX,trcs_cl_litr)
  ! Exchange atmospheric gases with litter using the matching surface-substep coefficient.
  implicit none
  integer, intent(in) :: I,J,NY,NX,M
  real(r8),intent(out) :: trcs_cl_litr(ids_beg:ids_end) !tracer concentration

  character(len=*), parameter :: subname = 'AtmosLitterExchangeM'
  integer :: K,ids,idg,idom
  real(r8) :: trc_gasq(idg_beg:idg_NH3)  !equilibrium tracer concentation
  real(r8) :: DFGcc(idg_beg:idg_NH3)
  real(r8) :: DLYR0,TORT0

  call PrintInfo('beg '//subname)

  IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX) .AND. VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN

    DLYR0 = AMAX1(ZERO2,DLYR_3D(3,0,NY,NX)) !vertical layer thickness
    ! Use the producer's M-specific footprint, not the last surface substep.
    TORT0 = TortMicPM_vr(M,0,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)/(0.5_r8*DLYR0)*FracSurfByLitRM_col(M,NY,NX)

    DO idg=idg_beg,idg_NH3
      DFGcc(idg)=SoluteDifusivitytscaledM_vr(idg,0,NY,NX)*TORT0
    ENDDO

    !temporary solute concentration
    DO ids=ids_beg,ids_end
      trcs_cl_litr(ids)=(trcs_solml2_vr(ids,0,NY,NX)-0._r8*trcs_solml_drib_vr(ids,0,NY,NX))/VLWatMicPM_vr(M,0,NY,NX)
    ENDDO
    ! Both soil destinations draw from the same unbanded litter solution.
    trcs_cl_litr(idg_NH3B)=trcs_cl_litr(idg_NH3)
    trcs_cl_litr(ids_NH4B:ids_H2PO4B)=trcs_cl_litr(ids_NH4:ids_H2PO4)
!
!     GasXferCoeffLitr2Atm_col: litter-patch gas path, including litter area/time [m3]
!
    DO idg=idg_beg,idg_NH3
      ! At zero litter area both interface coefficients vanish. There is no
      ! atmospheric transfer; retain the litter concentration for consumers.
      if(DFGcc(idg)+GasXferCoeffLitr2Atm_col(idg,M,NY,NX).LE.0._r8)then
        RGasAtmDisol2LitrM_col(idg,NY,NX)=0._r8
        cycle
      endif
      !equilbrium concentration at the air-litter interface
      trc_gasq(idg)=(GasXferCoeffLitr2Atm_col(idg,M,NY,NX)*GasSolbility_vr(idg,0,NY,NX)*AtmGasCgperm3_col(idg,NY,NX) &
        +DFGcc(idg)*trcs_cl_litr(idg))/(DFGcc(idg)+GasXferCoeffLitr2Atm_col(idg,M,NY,NX))

      !Dissolution into litter water, based on two-point gradient flow
      RGasAtmDisol2LitrM_col(idg,NY,NX)=(trc_gasq(idg)-trcs_cl_litr(idg))*AMIN1(VLWatMicPM_vr(M,0,NY,NX),DFGcc(idg))

    ENDDO

  ELSE
    RGasAtmDisol2LitrM_col(idg_beg:idg_NH3,NY,NX)=0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine AtmosLitterExchangeM

!------------------------------------------------------------------------------------------
  subroutine AtmosSoilExchangeM(I,J,M,NY,NX,trcs_cl_soil)
  !
  !Description
  !do soil-atmosphere tracer exchange 
  implicit none
  integer,  intent(in) :: I,J  
  integer,  intent(in) :: M,NY,NX
  real(r8), intent(out) :: trcs_cl_soil(ids_beg:ids_end)

  character(len=*), parameter :: subname='AtmosSoilExchangeM'
  real(r8) :: DLYR1,TORT1,VLWatMicPOA,VLWatMicPOB,VLWatMicPPA,VLWatMicPPB
  real(r8) :: DiffusivitySolutEff
  real(r8) :: trc_gsolc,trc_gsolc2
  integer  :: K,idg,idom,SurfaceGasID
!
!     SURFACE EXCHANGE OF AQUEOUS CO2, CH4, O2, N2, NH3
!     THROUGH VOLATILIZATION-DISSOLUTION FROM AQUEOUS
!     DIFFUSIVITIES IN SURFACE SOIL LAYER
!
  call PrintInfo('beg '//subname)

  IF(VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
    VLWatMicPOA = VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)*trcs_VLN_vr(ids_NO3,NU_col(NY,NX),NY,NX)
    VLWatMicPOB = VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)*trcs_VLN_vr(ids_NO3B,NU_col(NY,NX),NY,NX)
    VLWatMicPPA = VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
    VLWatMicPPB = VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)
    DLYR1       = AMAX1(ZERO2,DLYR_3D(3,NU_col(NY,NX),NY,NX))
    TORT1       = TortMicPM_vr(M,NU_col(NY,NX),NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)/(0.5_r8*DLYR1)

   !excldue  NH3 and NH3B
    DO idg=idg_beg,idg_NH3-1
      trcs_cl_soil(idg)=(trcs_solml2_vr(idg,NU_col(NY,NX),NY,NX)-0._r8*trcs_solml_drib_vr(idg,NU_col(NY,NX),NY,NX))/VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)
    ENDDO

    !deal with inorganic nutrients, including NH3 and NH3B
    IF(VLWatMicPMA_vr(NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(idg_NH3) = AZMAX1(trcs_solml2_vr(idg_NH3,NU_col(NY,NX),NY,NX)/VLWatMicPMA_vr(NU_col(NY,NX),NY,NX))
      trcs_cl_soil(ids_NH4) = AZMAX1(trcs_solml2_vr(ids_NH4,NU_col(NY,NX),NY,NX)/VLWatMicPMA_vr(NU_col(NY,NX),NY,NX))
    ELSE
      trcs_cl_soil(idg_NH3) = 0.0_r8
      trcs_cl_soil(ids_NH4) = 0.0_r8
    ENDIF

    IF(VLWatMicPOA.GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(ids_NO3) = AZMAX1(trcs_solml2_vr(ids_NO3,NU_col(NY,NX),NY,NX)/VLWatMicPOA)
      trcs_cl_soil(ids_NO2) = AZMAX1(trcs_solml2_vr(ids_NO2,NU_col(NY,NX),NY,NX)/VLWatMicPOA)
    ELSE
      trcs_cl_soil(ids_NO3) = 0.0_r8
      trcs_cl_soil(ids_NO2) = 0.0_r8
    ENDIF

    IF(VLWatMicPPA.GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(ids_H1PO4) = AZMAX1(trcs_solml2_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)/VLWatMicPPA)
      trcs_cl_soil(ids_H2PO4) = AZMAX1(trcs_solml2_vr(ids_H2PO4,NU_col(NY,NX),NY,NX)/VLWatMicPPA)
    ELSE
      trcs_cl_soil(ids_H1PO4) = 0.0_r8
      trcs_cl_soil(ids_H2PO4) = 0.0_r8
    ENDIF
    IF(VLWatMicPMB_vr(NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(idg_NH3B) = AZMAX1(trcs_solml2_vr(idg_NH3B,NU_col(NY,NX),NY,NX)/VLWatMicPMB_vr(NU_col(NY,NX),NY,NX))
      trcs_cl_soil(ids_NH4B) = AZMAX1(trcs_solml2_vr(ids_NH4B,NU_col(NY,NX),NY,NX)/VLWatMicPMB_vr(NU_col(NY,NX),NY,NX))
    ELSE
      trcs_cl_soil(idg_NH3B) = trcs_cl_soil(idg_NH3)
      trcs_cl_soil(ids_NH4B) = trcs_cl_soil(ids_NH4)
    ENDIF

    IF(VLWatMicPOB.GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(ids_NO3B) = AZMAX1(trcs_solml2_vr(ids_NO3B,NU_col(NY,NX),NY,NX)/VLWatMicPOB)
      trcs_cl_soil(ids_NO2B)  = AZMAX1(trcs_solml2_vr(ids_NO2B,NU_col(NY,NX),NY,NX)/VLWatMicPOB)
    ELSE
      trcs_cl_soil(ids_NO3B) = trcs_cl_soil(ids_NO3)
      trcs_cl_soil(ids_NO2B) = trcs_cl_soil(ids_NO2)
    ENDIF

    IF(VLWatMicPPB.GT.ZEROS2(NY,NX))THEN
      trcs_cl_soil(ids_H1PO4B) = AZMAX1(trcs_solml2_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)/VLWatMicPPB)
      trcs_cl_soil(ids_H2PO4B) = AZMAX1(trcs_solml2_vr(ids_H2PO4B,NU_col(NY,NX),NY,NX)/VLWatMicPPB)
    ELSE
      trcs_cl_soil(ids_H1PO4B) = trcs_cl_soil(ids_H1PO4)
      trcs_cl_soil(ids_H2PO4B) = trcs_cl_soil(ids_H2PO4)
    ENDIF
    !
    !     SURFACE VOLATILIZATION-DISSOLUTION FROM DIFFERENCES
    !     BETWEEN ATMOSPHERIC AND top-layer SOIL EQUILIBRIUM
    !     CONCENTRATIONS
    ! include NH3B

    DO idg=idg_beg,idg_end
      ! Banded NH3 has the same gas transport properties as unbanded NH3.
      SurfaceGasID=idg
      if(idg.EQ.idg_NH3B)SurfaceGasID=idg_NH3
      DiffusivitySolutEff = SoluteDifusivitytscaledM_vr(idg,NU_col(NY,NX),NY,NX)*TORT1
      trc_gsolc           = (CondGasXSnowM_col(SurfaceGasID,M,NY,NX)*AtmGasCgperm3_col(idg,NY,NX)*GasSolbility_vr(idg,NU_col(NY,NX),NY,NX) &
        +DiffusivitySolutEff*trcs_cl_soil(idg))/(DiffusivitySolutEff+CondGasXSnowM_col(SurfaceGasID,M,NY,NX))

      RGasAtmDisol2SoilM_col(idg,NY,NX)=(trc_gsolc-trcs_cl_soil(idg)) &
        *AMIN1(VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)*trcs_VLN_vr(idg,NU_col(NY,NX),NY,NX),DiffusivitySolutEff)
    ENDDO
  ELSE
    RGasAtmDisol2SoilM_col(idg_beg:idg_end,NY,NX)      = 0.0_r8
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine AtmosSoilExchangeM

!------------------------------------------------------------------------------------------

  subroutine Lit2SoilTracerAdvectM(I,J,M,NY,NX,FLWRM1,solute_adv_Lit2Soil_flxM,DOM_Adv_Litr2Soil_flxM,&
    solute_adv_Lit2MacP_flxM,DOM_Adv_Litr2MacP_flxM)
  !
  ! Description:
  ! Litter and top soil tracer exchange by (upstream) advection
  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M, NY, NX
  real(r8),intent(out) :: solute_adv_Lit2Soil_flxM(ids_beg:ids_end)
  real(r8),intent(out) :: FLWRM1    !water flow from litter to soil at iteration M
  real(r8),intent(out) :: DOM_Adv_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)
  real(r8),intent(out) :: solute_adv_Lit2MacP_flxM(ids_beg:ids_end)
  real(r8),intent(out) :: DOM_Adv_Litr2MacP_flxM(idom_beg:idom_end,1:jcplx)

  character(len=*), parameter :: subname='Lit2SoilTracerAdvectM'
  REAL(R8) :: VFLW,FracMic,FracMac,GasFracMic
  integer :: K,idn,idom,idg,idn1,ids

  call PrintInfo('beg '//subname)

  FLWRM1=WatFLoLitr2SoilM_col(M,NY,NX)
  solute_adv_Lit2Soil_flxM=0._r8; DOM_Adv_Litr2Soil_flxM=0._r8
  solute_adv_Lit2MacP_flxM=0._r8; DOM_Adv_Litr2MacP_flxM=0._r8
  call LitterInfiltrationFractions(FLWRM1,WatFLoLitr2MacPM_col(M,NY,NX),SurfRunoffPotentM_col(M,NY,NX),&
    VLWatMicPM_vr(M,0,NY,NX),ZEROS2(NY,NX),VFLWX,FracMic,FracMac)
!
!
  !water flow out of litter to soil
  IF(FLWRM1.GT.0.0_r8)THEN
    VFLW=FracMic
    GasFracMic=0._r8
    IF(VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))GasFracMic=MIN(VFLWX,FLWRM1/VLWatMicPM_vr(M,0,NY,NX))

    DO  K=1,micpar%NumOfLitrCmplxs
      DO idom=idom_beg,idom_end
        DOM_Adv_Litr2Soil_flxM(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,0,NY,NX))
      ENDDO
    ENDDO

    !all volatiles
    DO idg=idg_beg,idg_end
      solute_adv_Lit2Soil_flxM(idg)=GasFracMic*AZMAX1(trcs_solml2_vr(idg,0,NY,NX))
    enddo
    ! Partition the shared litter NH3 solution between non-band and band soil.
    solute_adv_Lit2Soil_flxM(idg_NH3)=solute_adv_Lit2Soil_flxM(idg_NH3) &
      *trcs_VLN_vr(idg_NH3,NU_col(NY,NX),NY,NX)

    !from NH4 to H2PO4
    DO idn=ids_nut_beg,ids_nuts_end
      solute_adv_Lit2Soil_flxM(idn)=VFLW*AZMAX1(trcs_solml2_vr(idn,0,NY,NX))*trcs_VLN_vr(idn,NU_col(NY,NX),NY,NX)
    ENDDO

    !from NH4B to ids_H2PO4B
    DO idn=0,ids_nuts     
      idn1=ids_NH4B+idn
      solute_adv_Lit2Soil_flxM(idn1)=VFLW*AZMAX1(trcs_solml2_vr(idn+ids_NH4,0,NY,NX))*trcs_VLN_vr(idn1,NU_col(NY,NX),NY,NX)
    ENDDO

    solute_adv_Lit2Soil_flxM(idg_NH3B)=GasFracMic*AZMAX1(trcs_solml2_vr(idg_NH3,0,NY,NX))*trcs_VLN_vr(idg_NH3B,NU_col(NY,NX),NY,NX)

    ! flow from soil to litter
  ELSEIF(FLWRM1.LT.0._r8)THEN
    IF(VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FLWRM1/VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    DO K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_Adv_Litr2Soil_flxM(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,NU_col(NY,NX),NY,NX))
      ENDDO
    ENDDO

    DO ids=ids_beg,ids_end
      solute_adv_Lit2Soil_flxM(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,NU_col(NY,NX),NY,NX))
    ENDDO
  ENDIF
  ! Macropore infiltration draws from litter micropore solution. Keep gas
  ! coupling on its existing path; transport DOM and nonvolatile nutrients.
  associate(NumLitr => micpar%NumOfLitrCmplxs)
  DO K=1,NumLitr
    DO idom=idom_beg,idom_end
      DOM_Adv_Litr2MacP_flxM(idom,K)=FracMac*AZMAX1(DOM_MicP2_vr(idom,K,0,NY,NX))
    ENDDO
  ENDDO
  end associate
  DO idn=0,ids_nuts
    ids=ids_NH4+idn; idn1=ids_NH4B+idn
    solute_adv_Lit2MacP_flxM(ids)=FracMac*AZMAX1(trcs_solml2_vr(ids,0,NY,NX))*trcs_VLN_vr(ids,NU_col(NY,NX),NY,NX)
    solute_adv_Lit2MacP_flxM(idn1)=FracMac*AZMAX1(trcs_solml2_vr(ids,0,NY,NX))*trcs_VLN_vr(idn1,NU_col(NY,NX),NY,NX)
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine Lit2SoilTracerAdvectM
!------------------------------------------------------------------------------------------

  subroutine LitterSoilTracerDiffusionM(I,J,M,NY,NX,FLWRM1,trcs_cl_litr,trcs_cl_soil,&
    trcs_Dif_Litr2Soil_flxM,DOM_Difus_Litr2Soil_flxM)
  !
  !Description
  ! Diffusion flux between litter and soil  
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M, NY, NX
  real(r8), intent(in) :: FLWRM1                         !water flow from litter to topsoil
  real(r8), intent(in) :: trcs_cl_litr(ids_beg:ids_end)  !tracer concentation in litter
  real(r8), intent(in) :: trcs_cl_soil(ids_beg:ids_end)  !tracer concentation in topsoil
  real(r8), intent(out):: trcs_Dif_Litr2Soil_flxM(ids_beg:ids_end)   !tracer diffusion from litter to soil
  real(r8), intent(out):: DOM_Difus_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)  !tracer diffusion from litter to soil
  real(r8) :: CDOM_MicP_litr(idom_beg:idom_end,1:jcplx)  !DOM concentation in litter
  real(r8) :: CDOM_MicP_soil(idom_beg:idom_end,1:jcplx)  !DOM concentation in topsoil

  character(len=*), parameter :: subname='LitterSoilTracerDiffusionM'
  real(r8) :: TORT0,TORT1
  real(r8) :: DLYR0,DLYR1
  real(r8) :: DIFDOM(idom_beg:idom_end)

  real(r8) :: DIFDOM0
  real(r8) :: DIFDOM1
  real(r8) :: DISPN,DIF0,DIF1
  real(r8) :: SDifc(ids_beg:ids_end)
  integer  :: K,idn,idg,idom
!
!     VOLT,DLYR,AREA=soil surface volume, thickness, area
!     VLWatMicPM=micropore water-filled porosity from watsub.f
!
  call PrintInfo('beg '//subname)

  IF((VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX) .AND. VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX)) &
    .AND. (VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX)))THEN
!
    DO  K=1,micpar%NumOfLitrCmplxs
      do idom=idom_beg,idom_end
        CDOM_MicP_litr(idom,K)=AZMAX1(DOM_MicP2_vr(idom,K,0,NY,NX)/VLWatMicPM_vr(M,0,NY,NX))
      enddo
    ENDDO
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        CDOM_MicP_soil(idom,K)=AZMAX1(DOM_MicP2_vr(idom,K,NU_col(NY,NX),NY,NX)/VLWatMicPM_vr(M,NU_col(NY,NX),NY,NX))
      enddo
    ENDDO

!     DIFFUSIVITIES IN RESIDUE AND SOIL SURFACE

    DLYR0 = AMAX1(ZERO2,DLYR_3D(3,0,NY,NX))
    TORT0 = TortMicPM_vr(M,0,NY,NX)/DLYR0*FracSurfByLitR_col(NY,NX)
    DLYR1 = AMAX1(ZERO2,DLYR_3D(3,NU_col(NY,NX),NY,NX))
    TORT1 = TortMicPM_vr(M,NU_col(NY,NX),NY,NX)/DLYR1
    DISPN = (0.20_r8*(0.5_r8*DLYR1)**1.07_r8)*AMIN1(VFLWX,ABS(FLWRM1/AREA_3D(3,NU_col(NY,NX),NY,NX)))

    DO idom=idom_beg,idom_end
      DIFDOM0      = (DOM_diffusivitytscaledM_vr(idom,0,NY,NX)*TORT0+DISPN)
      DIFDOM1      = (DOM_diffusivitytscaledM_vr(idom,NU_col(NY,NX),NY,NX)*TORT1+DISPN)
      DIFDOM(idom) = DIFDOM0*DIFDOM1/(DIFDOM0+DIFDOM1)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    ENDDO

    DO idn=ids_beg,ids_end
      DIF0        = (SoluteDifusivitytscaledM_vr(idn,0,NY,NX)*TORT0+DISPN)
      DIF1        = (SoluteDifusivitytscaledM_vr(idn,NU_col(NY,NX),NY,NX)*TORT1+DISPN)
      SDifc(idn) = DIF0*DIF1/(DIF0+DIF1)*AREA_3D(3,NU_col(NY,NX),NY,NX)
    ENDDO
!
    DO  K=1,jcplx
      DO idom=idom_beg,idom_end
        DOM_Difus_Litr2Soil_flxM(idom,K)=DIFDOM(idom)*(CDOM_MicP_litr(idom,K)-CDOM_MicP_soil(idom,K))
      ENDDO
    ENDDO

    !exclude NH3B and NH3
    DO idg=idg_beg,idg_NH3-1
        trcs_Dif_Litr2Soil_flxM(idg)=SDifc(idg)*(trcs_cl_litr(idg)-trcs_cl_soil(idg))
    ENDDO

    !include NH3B and NH3
    DO idn=ids_nuts_beg,ids_nuts_end
      trcs_Dif_Litr2Soil_flxM(idn)=SDifc(idn)*(trcs_cl_litr(idn)-trcs_cl_soil(idn))*trcs_VLN_vr(idn,NU_col(NY,NX),NY,NX)
    ENDDO
  ELSE
    DO  K=1,jcplx
      DOM_Difus_Litr2Soil_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO
    trcs_Dif_Litr2Soil_flxM(ids_beg:ids_end)=0._r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine LitterSoilTracerDiffusionM
!------------------------------------------------------------------------------------------

  subroutine SurfLayerNetTracerFluxM(I,J,M,NY,NX,trcs_Dif_Litr2Soil_flxM,solute_adv_Lit2Soil_flxM,&
    DOM_Adv_Litr2Soil_flxM,DOM_Difus_Litr2Soil_flxM,solute_adv_Lit2MacP_flxM,DOM_Adv_Litr2MacP_flxM,snow_mac_frac)
  !
  !Description:
  !
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NY, NX
  real(r8), intent(in) :: snow_mac_frac ! fraction of soil-directed snow liquid entering macropores
  real(r8), intent(in) :: trcs_Dif_Litr2Soil_flxM(ids_beg:ids_end)              !diffusion flux
  real(r8), intent(in) :: solute_adv_Lit2Soil_flxM(ids_beg:ids_end)             !advection flux
  real(r8), intent(in) :: DOM_Adv_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)     !diffusion advection flux
  real(r8), intent(in) :: DOM_Difus_Litr2Soil_flxM(idom_beg:idom_end,1:jcplx)   !DOM diffusion flux
  real(r8), intent(in) :: solute_adv_Lit2MacP_flxM(ids_beg:ids_end)
  real(r8), intent(in) :: DOM_Adv_Litr2MacP_flxM(idom_beg:idom_end,1:jcplx)
  character(len=*), parameter :: subname ='SurfLayerNetTracerFluxM'   
  integer :: K,idg,idom,ids,ids1,ids2
  real(r8):: dom_flx,NutLitr2Soil

  call PrintInfo('beg '//subname)

  !litter layer, DOM from precipitation can be added here
  DO K=1,micpar%NumOfLitrCmplxs
    DO idom=idom_beg,idom_end
      dom_flx=DOM_Adv_Litr2Soil_flxM(idom,K)+DOM_Difus_Litr2Soil_flxM(idom,K)
      DOM_MicpTranspFlxM_3D(idom,K,3,0,NY,NX)         = -dom_flx
      DOM_MicpTranspFlxM_3D(idom,K,3,NU_col(NY,NX),NY,NX) = dom_flx
    ENDDO
  ENDDO

  DO idg=idg_beg,idg_NH3
    Gas_litr2Soil_flxM_col(idg,NY,NX)= solute_adv_Lit2Soil_flxM(idg)+trcs_Dif_Litr2Soil_flxM(idg) 

    trcs_MicpTranspFlxM_3D(idg,3,0,NY,NX)=trcg_Precip2LitrM_col(idg,NY,NX)+trcg_AquaADV_Snow2Litr_flxM(idg,NY,NX) &
      -Gas_litr2Soil_flxM_col(idg,NY,NX)

    trcs_MicpTranspFlxM_3D(idg,3,NU_col(NY,NX),NY,NX)=trcs_Precip2MicpM_col(idg,NY,NX)+trcg_AquaADV_Snow2Soil_flxM(idg,NY,NX)*(1._r8-snow_mac_frac) &
      +Gas_litr2Soil_flxM_col(idg,NY,NX)
  ENDDO
  !from NH4 to H2PO4
  do ids=ids_nut_beg,ids_nuts_end
    trcs_MicpTranspFlxM_3D(ids,3,0,NY,NX)=trcn_Precip2LitrM_col(ids,NY,NX)+trcn_AquaADV_Snow2Litr_flxM(ids,NY,NX) &
      -solute_adv_Lit2Soil_flxM(ids)-trcs_Dif_Litr2Soil_flxM(ids)

    trcs_MicpTranspFlxM_3D(ids,3,NU_col(NY,NX),NY,NX)=trcs_Precip2MicpM_col(ids,NY,NX)+trcn_AquaADV_Snow2Soil_flxM(ids,NY,NX)*(1._r8-snow_mac_frac) &
      +solute_adv_Lit2Soil_flxM(ids)+trcs_Dif_Litr2Soil_flxM(ids)
  enddo

  !take off the band 
  trcs_MicpTranspFlxM_3D(idg_NH3,3,0,NY,NX)   = trcs_MicpTranspFlxM_3D(idg_NH3,3,0,NY,NX) &
    -solute_adv_Lit2Soil_flxM(idg_NH3B)-trcs_Dif_Litr2Soil_flxM(idg_NH3B)  

  trcs_MicpTranspFlxM_3D(idg_NH3B,3,NU_col(NY,NX),NY,NX)   = trcs_MicpTranspFlxM_3D(idg_NH3B,3,NU_col(NY,NX),NY,NX) &
    +trcn_AquaADV_Snow2Band_flxM(idg_NH3B,NY,NX)*(1._r8-snow_mac_frac) &
    +solute_adv_Lit2Soil_flxM(idg_NH3B)+trcs_Dif_Litr2Soil_flxM(idg_NH3B)

  Gas_litr2Soil_flxM_col(idg_NH3B,NY,NX)= solute_adv_Lit2Soil_flxM(idg_NH3B)+trcs_Dif_Litr2Soil_flxM(idg_NH3B)  

  ![NH4,NO3,NO2,H1PO4,H2PO4]
  ![NH4B,NO3B,NO2B,H2PO4B,H2PO4B]
  do ids=0,ids_nuts
    ids1=ids+ids_NH4
    ids2=ids+ids_NH4B  
    NutLitr2Soil=solute_adv_Lit2Soil_flxM(ids2)+trcs_Dif_Litr2Soil_flxM(ids2)
    trcs_MicpTranspFlxM_3D(ids1,3,0,NY,NX)             = trcs_MicpTranspFlxM_3D(ids1,3,0,NY,NX)-NutLitr2Soil
    trcs_MicpTranspFlxM_3D(ids2,3,NU_col(NY,NX),NY,NX) = trcs_Precip2MicpM_col(ids2,NY,NX)+trcn_AquaADV_Snow2Band_flxM(ids2,NY,NX)*(1._r8-snow_mac_frac)+NutLitr2Soil
  enddo

  ! Paired litter debit / topsoil macropore credit. The soil flux gatherer
  ! consumes this top boundary just like the other macropore faces.
  associate(NumLitr => micpar%NumOfLitrCmplxs)
  DO K=1,NumLitr
    DO idom=idom_beg,idom_end
      DOM_MicpTranspFlxM_3D(idom,K,3,0,NY,NX)=DOM_MicpTranspFlxM_3D(idom,K,3,0,NY,NX)-DOM_Adv_Litr2MacP_flxM(idom,K)
      DOM_MacpTranspFlxM_3D(idom,K,3,NU_col(NY,NX),NY,NX)=DOM_MacpTranspFlxM_3D(idom,K,3,NU_col(NY,NX),NY,NX)+DOM_Adv_Litr2MacP_flxM(idom,K)
    ENDDO
  ENDDO
  end associate
  DO ids=0,ids_nuts
    ids1=ids_NH4+ids; ids2=ids_NH4B+ids
    trcs_MicpTranspFlxM_3D(ids1,3,0,NY,NX)=trcs_MicpTranspFlxM_3D(ids1,3,0,NY,NX) &
      -solute_adv_Lit2MacP_flxM(ids1)-solute_adv_Lit2MacP_flxM(ids2)
    trcs_MacpTranspFlxM_3D(ids1,3,NU_col(NY,NX),NY,NX)=trcs_MacpTranspFlxM_3D(ids1,3,NU_col(NY,NX),NY,NX)+solute_adv_Lit2MacP_flxM(ids1)
    trcs_MacpTranspFlxM_3D(ids2,3,NU_col(NY,NX),NY,NX)=trcs_MacpTranspFlxM_3D(ids2,3,NU_col(NY,NX),NY,NX)+solute_adv_Lit2MacP_flxM(ids2)
  ENDDO
  ! The total snow flux remains unchanged for the donor and hourly ledgers.
  DO idg=idg_beg,idg_NH3
    trcs_MacpTranspFlxM_3D(idg,3,NU_col(NY,NX),NY,NX)=trcs_MacpTranspFlxM_3D(idg,3,NU_col(NY,NX),NY,NX) &
      +trcg_AquaADV_Snow2Soil_flxM(idg,NY,NX)*snow_mac_frac
  ENDDO
  DO ids=ids_nut_beg,ids_nuts_end
    trcs_MacpTranspFlxM_3D(ids,3,NU_col(NY,NX),NY,NX)=trcs_MacpTranspFlxM_3D(ids,3,NU_col(NY,NX),NY,NX) &
      +trcn_AquaADV_Snow2Soil_flxM(ids,NY,NX)*snow_mac_frac
  ENDDO
  DO ids=ids_nutb_beg,ids_nutb_end
    trcs_MacpTranspFlxM_3D(ids,3,NU_col(NY,NX),NY,NX)=trcs_MacpTranspFlxM_3D(ids,3,NU_col(NY,NX),NY,NX) &
      +trcn_AquaADV_Snow2Band_flxM(ids,NY,NX)*snow_mac_frac
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SurfLayerNetTracerFluxM

!------------------------------------------------------------------------------------------

  subroutine TracerFlowThruSurRofM(I,J,M,NHE,NHW,NVS,NVN)
  !
  !Description:
  !Do inter-grid lateral transport at the surface
  !
  implicit none
  integer, intent(in) :: I,J, M, NHE, NHW, NVS, NVN
  character(len=*), parameter :: subname='TracerFlowThruSurRofM'
  
  real(r8) :: FQRM,VFLW
  integer :: K,N,NN,N1,N2,N4,N5,N4B,N5B,idg,idn,NY,NX
  integer :: idom
!
  call PrintInfo('beg '//subname)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      N1=NX;N2=NY

      !there is outflow from grid (N2,N1)
      IF(SurfRunoffPotentM_col(M,N2,N1).GT.ZEROS(N2,N1))THEN
        !Obtain potential tracer flux through surface runoff, include both inner grids and boundaries  
        IF(VLWatMicPM_vr(M,0,N2,N1).GT.ZEROS2(N2,N1))THEN
          VFLW=AMIN1(VFLWX,SurfRunoffPotentM_col(M,N2,N1)/VLWatMicPM_vr(M,0,N2,N1))
        ELSE
          VFLW=VFLWX
        ENDIF

        DO  K=1,micpar%NumOfLitrCmplxs
          DO idom=idom_beg,idom_end
            DOM_FloXSurRunoff_PotFlxM(idom,K,N2,N1)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,0,N2,N1))
          enddo
        ENDDO

        DO idg=idg_beg,idg_NH3
          trcg_FloXSurRunoff_PotFlxM(idg,N2,N1)=VFLW*AZMAX1(trcs_solml2_vr(idg,0,N2,N1))
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          trcn_FloXSurRunoff_PotFlxM(idn,N2,N1)=VFLW*AZMAX1(trcs_solml2_vr(idn,0,N2,N1))
        ENDDO
        !-----------------------------------------------------------------
        !     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
        !
        D4311: DO N=1,2
          D4306: DO NN=1,2

            IF(N.EQ.iWestEastDirection)THEN  
              IF((NX.EQ.NHE .AND. NN.EQ.iFront) &         !skip eastern boundary
                .OR. (NX.EQ.NHW .AND. NN.EQ.iBehind))THEN !skip western boundary
                cycle
              ELSE
                N4  = NX+1;N5  = NY
                N4B = NX-1;N5B = NY
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN 
              IF((NY.EQ.NVS .AND. NN.EQ.iFront) &          !skip southern boundary
                .OR. (NY.EQ.NVN .AND. NN.EQ.iBehind))THEN  !skip northern boundary
                cycle
              ELSE
                N4  = NX;N5  = NY+1
                N4B = NX;N5B = NY-1
              ENDIF
            ENDIF      
            !
            ! IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
            !
            IF(NN.EQ.iFront)THEN  !(N2,N1) is the front grid
              FQRM=QflxSurfRunoffM_2DH(M,N,iBehind,N5,N4)/SurfRunoffPotentM_col(M,N2,N1)
              !loss from grid (N2,N1) to (N5,N4), note (N2,N1) is behind (N5,N4)
              IF(FQRM.GT.ZEROS(N2,N1))then
                DO  K=1,jcplx
                  do idom=idom_beg,idom_end
                    DOM_FloXSurRof_flxM_2DH(idom,K,N,iBehind,N5,N4)=DOM_FloXSurRunoff_PotFlxM(idom,K,N2,N1)*FQRM
                  enddo
                ENDDO

                DO idg=idg_beg,idg_NH3
                  trcg_SurRof_flxM_2DH(idg,N,iBehind,N5,N4)=trcg_FloXSurRunoff_PotFlxM(idg,N2,N1)*FQRM
                ENDDO

                DO idn=ids_nut_beg,ids_nuts_end
                  trcn_SurRof_flxM_2DH(idn,N,iBehind,N5,N4)=trcn_FloXSurRunoff_PotFlxM(idn,N2,N1)*FQRM
                ENDDO
              ENDIF

              !  IF OVERLAND FLOW IS FROM CURRENT TO ADJACENT GRID CELL
            ELSEIF(NN.EQ.iBehind)THEN        
              !legitimate grid
              IF(N4B.GT.0 .AND. N5B.GT.0)THEN
                FQRM=QflxSurfRunoffM_2DH(M,N,iFront,N5B,N4B)/SurfRunoffPotentM_col(M,N2,N1)
                !loss from grid (N2,N1) to (N5B,N4B)
                IF(FQRM.GT.ZEROS(N2,N1))then
                  DO  K=1,jcplx
                    do idom=idom_beg,idom_end
                      DOM_FloXSurRof_flxM_2DH(idom,K,N,iFront,N5B,N4B)=DOM_FloXSurRunoff_PotFlxM(idom,K,N2,N1)*FQRM
                    enddo
                  ENDDO

                  DO idg=idg_beg,idg_NH3
                    trcg_SurRof_flxM_2DH(idg,N,iFront,N5B,N4B)=trcg_FloXSurRunoff_PotFlxM(idg,N2,N1)*FQRM
                  ENDDO

                  DO idn=ids_nut_beg,ids_nuts_end
                    trcn_SurRof_flxM_2DH(idn,N,iFront,N5B,N4B)=trcn_FloXSurRunoff_PotFlxM(idn,N2,N1)*FQRM
                  ENDDO
                ENDIF
              ENDIF  
            ENDIF        
          ENDDO  D4306
        ENDDO D4311
      ELSE
        DO K=1,jcplx
          DOM_FloXSurRunoff_PotFlxM(idom_beg:idom_end,K,N2,N1)=0.0_r8
        ENDDO
        trcg_FloXSurRunoff_PotFlxM(idg_beg:idg_NH3,N2,N1)          = 0.0_r8
        trcn_FloXSurRunoff_PotFlxM(ids_nut_beg:ids_nuts_end,N2,N1) = 0.0_r8
      ENDIF
    ENDDO
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine TracerFlowThruSurRofM
!------------------------------------------------------------------------------------------

  subroutine TracerFlowThruSnowRedist(I,J,M,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J,M,NHE,NHW,NVS,NVN
  integer :: NY,NX,N1,N2,N4,N5,N4B,N5B,N,NN
  integer :: idg, idn
  real(r8):: VFLW

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      N1=NX;N2=NY
      !----------------------------------
      !DO SNOW REDISTRIBUTION
      !
      D4311: DO N=1,2
        D4306: DO NN=1,2
          !
          IF(N.EQ.iWestEastDirection)THEN
            IF((NX.EQ.NHE .AND. NN.EQ.iFront) &         !skip eastern boundary
              .OR. (NX.EQ.NHW .AND. NN.EQ.iBehind))THEN !skip western boundary
              cycle
            ELSE
              N4  = NX+1;N5  = NY
              N4B = NX-1;N5B = NY
            ENDIF
          ELSEIF(N.EQ.iNorthSouthDirection)THEN            
            IF((NY.EQ.NVS .AND. NN.EQ.iFront) &          !skip southern boundary
              .OR. (NY.EQ.NVN .AND. NN.EQ.iBehind))THEN  !skip northern boundary
              cycle
            ELSE
              N4  = NX;N5  = NY+1
              N4B = NX;N5B = NY-1
            ENDIF
          ENDIF

          IF(NN.EQ.iFront)THEN
            !
            !   IF NO SNOW DRIFT THEN NO TRANSPORT
            !
            IF(DrySnoFlxByRedistM_2DH(M,N,iBehind,N5,N4).GT.ZEROS2(N2,N1))THEN
              IF(VcumSnoDWI_col(N2,N1).GT.ZEROS2(N2,N1))THEN
                VFLW=AZMAX1(AMIN1(VFLWX,DrySnoFlxByRedistM_2DH(M,N,iBehind,N5,N4)/VcumSnoDWI_col(N2,N1)))
              ELSE
                VFLW=VFLWX
              ENDIF

              DO idg=idg_beg,idg_NH3
                trcg_SnowDrift_flxM_2DH(idg,N,iBehind,N5,N4) = VFLW*AZMAX1(trcg_solsml2_snvr(idg,1,N2,N1))
              ENDDO

              DO idn=ids_nut_beg,ids_nuts_end
                trcn_SnowDrift_flxM_2DH(idn,N,iBehind,N5,N4) = VFLW*AZMAX1(trcn_solsml2_snvr(idn,1,N2,N1))
              ENDDO
            ENDIF  
            !check legitimate behind grid  
          ELSEIF(NN.EQ.iBehind .AND. N5B.GT.0 .AND. N4B.GT.0)THEN   

            IF(DrySnoFlxByRedistM_2DH(M,N,iFront,N5B,N4B).GT.ZEROS2(N2,N1))THEN
              IF(VcumSnoDWI_col(N2,N1).GT.ZEROS2(N2,N1))THEN
                VFLW=AZMAX1(AMIN1(VFLWX,DrySnoFlxByRedistM_2DH(M,N,iFront,N5B,N4B)/VcumSnoDWI_col(N2,N1)))
              ELSE
                VFLW=VFLWX
              ENDIF
              DO idg=idg_beg,idg_NH3
                trcg_SnowDrift_flxM_2DH(idg,N,iFront,N5B,N4B) = VFLW*AZMAX1(trcg_solsml2_snvr(idg,1,N2,N1))
              ENDDO

              DO idn=ids_nut_beg,ids_nuts_end
                trcn_SnowDrift_flxM_2DH(idn,N,iFront,N5B,N4B) = VFLW*AZMAX1(trcn_solsml2_snvr(idn,1,N2,N1))
              ENDDO          
            ENDIF
          ENDIF  
          
        ENDDO D4306
      ENDDO D4311
    ENDDO
  ENDDO  
  end subroutine TracerFlowThruSnowRedist
! ----------------------------------------------------------------------
  subroutine SnowSoluteVertDrainM(I,J,M,NHE,NHW,NVS,NVN,snow_donor_col,band_transfer_only)
  !
  !Description
  !Vertical drainage of solute in and out of snow.
  !out of snow fluxes are added to litter, soil and nutrient band, accordingly. 
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NHE,NHW,NVS,NVN
  ! Discharge donor within the fixed snow slots; zero means none.
  integer, intent(out) :: snow_donor_col(JY,JX)
  ! True only in the once-per-substep ID-changing surface transfer pass.
  logical, optional, intent(in) :: band_transfer_only
  logical :: band_pass
  integer :: NY, NX

  character(len=*), parameter :: subname='SnowSoluteVertDrainM'
  logical :: LCHKBOTML  !flag for dealing with the bottom snow layer
  integer :: L,L2
  real(r8) :: VFLWS,VFLWW,VFLWR
  real(r8) :: VFLWNOB,VFLWNO3,VFLWNHB
  real(r8) :: VFLWNH4,VFLWPO4,VFLWPOB
  integer  :: idg,idn
!
  call PrintInfo('beg '//subname)
  snow_donor_col=0
  band_pass=.false.
  if(present(band_transfer_only))band_pass=band_transfer_only

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      LCHKBOTML=.false.
      ! Hydrology can create snow in slot 1 before nsnol_col is updated.
      ! Scan fixed slots; substep heat capacity determines drainage donors.
      DO  L=1,JS
        !check having meaningful snow mass
        IF(VLSnowHeatCapM_snvr(M,L,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX)*1.e-4_r8)THEN
          !next layer index
          L2=MIN(JS,L+1)
          !L is an inner layer, because L2 is well defined
          IF(L.LT.JS .AND. VLSnowHeatCapM_snvr(M,L2,NY,NX).GT.VLHeatCapSnowMin_col(NY,NX))THEN
            !some flowout
            IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
              VFLWW=AZMAX1(AMIN1(1.0_r8,WatFlowInSnowM_snvr(M,L2,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
            ELSE          
              ! Intentional numerical redistribution for negligible liquid:
              ! transfer a small tracer fraction even when water flow is zero.
              ! This is not water-driven advection; its cumulative effect
              ! depends on the number of transport substeps.
              VFLWW=1.e-8_r8
            ENDIF
            !flux from layer L into layer L2
            DO idg=idg_beg,idg_NH3
              trcg_AquaAdv_flxM_snvr(idg,L2,NY,NX) = trcg_solsml2_snvr(idg,L,NY,NX)*VFLWW
            ENDDO

            DO idn=ids_nut_beg,ids_nuts_end
              trcn_AquaAdv_flxM_snvr(idn,L2,NY,NX) = trcn_solsml2_snvr(idn,L,NY,NX)*VFLWW
            ENDDO
            !L == JS: bottom layer or L<LS and next layer has no significant snow
          ELSE
            !nothing to add to L2
            IF(L.LT.JS)THEN
              trcg_AquaAdv_flxM_snvr(idg_beg:idg_NH3,L2,NY,NX)          = 0.0_r8
              trcn_AquaAdv_flxM_snvr(ids_nut_beg:ids_nuts_end,L2,NY,NX) = 0.0_r8
            ENDIF
            !
            !     SNOWPACK SOLUTE DISCHARGE TO SURFACE LITTER, SOIL SURFACE
            !
            IF(.not.LCHKBOTML)THEN
              IF(VLWatSnow_snvr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
                VFLWR=AZMAX1(AMIN1(1.0_r8,WatFlowSno2LitRM_col(M,NY,NX)/VLWatSnow_snvr(L,NY,NX)))
                VFLWS=AZMAX1(AMIN1(1.0_r8,(WatFlowSno2MicPM_col(M,NY,NX)+WatFlowSno2MacPM_col(M,NY,NX))/VLWatSnow_snvr(L,NY,NX)))
              ELSE
                VFLWR=FracSurfByLitR_col(NY,NX)
                VFLWS=FracSurfBareSoil_col(NY,NX)
              ENDIF
              VFLWNH4 = VFLWS*trcs_VLN_vr(ids_NH4,NU_col(NY,NX),NY,NX)
              VFLWNHB = VFLWS*trcs_VLN_vr(ids_NH4B,NU_col(NY,NX),NY,NX)
              VFLWNO3 = VFLWS*trcs_VLN_vr(ids_NO3,NU_col(NY,NX),NY,NX)
              VFLWNOB = VFLWS*trcs_VLN_vr(ids_NO3B,NU_col(NY,NX),NY,NX)
              VFLWPO4 = VFLWS*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)
              VFLWPOB = VFLWS*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)

              !snow flow to litter layer
              DO idg=idg_beg,idg_NH3
                trcg_AquaADV_Snow2Litr_flxM(idg,NY,NX)=trcg_solsml2_snvr(idg,L,NY,NX)*VFLWR
              ENDDO

              DO idn=ids_nut_beg,ids_nuts_end
                trcn_AquaADV_Snow2Litr_flxM(idn,NY,NX)=trcn_solsml2_snvr(idn,L,NY,NX)*VFLWR
              ENDDO
              !-----------------------
              !snow flow to topsoil 
              DO idg=idg_beg,idg_NH3-1
                trcg_AquaADV_Snow2Soil_flxM(idg,NY,NX)  = trcg_solsml2_snvr(idg,L,NY,NX)*VFLWS
              ENDDO

              trcg_AquaADV_Snow2Soil_flxM(idg_NH3,NY,NX)  = trcg_solsml2_snvr(idg_NH3,L,NY,NX)*VFLWNH4
              trcn_AquaADV_Snow2Band_flxM(idg_NH3B,NY,NX) = trcg_solsml2_snvr(idg_NH3,L,NY,NX)*VFLWNHB

              trcn_AquaADV_Snow2Soil_flxM(ids_NO2,NY,NX)    = 0._r8
              trcn_AquaADV_Snow2Soil_flxM(ids_NH4,NY,NX)    = trcn_solsml2_snvr(ids_NH4,L,NY,NX)*VFLWNH4
              trcn_AquaADV_Snow2Soil_flxM(ids_NO3,NY,NX)    = trcn_solsml2_snvr(ids_NO3,L,NY,NX)*VFLWNO3
              trcn_AquaADV_Snow2Soil_flxM(ids_H1PO4,NY,NX)  = trcn_solsml2_snvr(ids_H1PO4,L,NY,NX)*VFLWPO4
              trcn_AquaADV_Snow2Soil_flxM(ids_H2PO4,NY,NX)  = trcn_solsml2_snvr(ids_H2PO4,L,NY,NX)*VFLWPO4

              trcn_AquaADV_Snow2Band_flxM(ids_NO2B,NY,NX)   = 0._r8
              trcn_AquaADV_Snow2Band_flxM(ids_NH4B,NY,NX)   = trcn_solsml2_snvr(ids_NH4,L,NY,NX)*VFLWNHB
              trcn_AquaADV_Snow2Band_flxM(ids_NO3B,NY,NX)   = trcn_solsml2_snvr(ids_NO3,L,NY,NX)*VFLWNOB
              trcn_AquaADV_Snow2Band_flxM(ids_H1PO4B,NY,NX) = trcn_solsml2_snvr(ids_H1PO4,L,NY,NX)*VFLWPOB
              trcn_AquaADV_Snow2Band_flxM(ids_H2PO4B,NY,NX) = trcn_solsml2_snvr(ids_H2PO4,L,NY,NX)*VFLWPOB
              
              snow_donor_col(NY,NX)=L
              LCHKBOTML = .true.
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  ! Snow-to-band transfer is applied only in the separate surface pass.
  if(.not.band_pass)then
    trcn_AquaADV_Snow2Band_flxM(:,NVN:NVS,NHW:NHE)=0._r8
  endif
  call PrintInfo('end '//subname)
  end subroutine SnowSoluteVertDrainM

! ----------------------------------------------------------------------
  subroutine SoluteAdvDifusTranspM(I,J,M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: I,J,M,N
  integer, intent(in) :: N1,N2,N3  !source grid
  integer, intent(in) :: N4,N5,N6  !dest grid
  character(len=*), parameter :: subname='SoluteAdvDifusTranspM'
  real(r8)  :: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)  

  integer  :: K,ids,idom

  call PrintInfo('beg '//subname)

  !------------------------------------------------
  !     SOLUTE TRANSPORT IN MICROPORES
  !
  call MicroporeSoluteAdvectionM(I,J,M,N,N1,N2,N3,N4,N5,N6)
  !
  call MicroporeSoluteDiffusionM(I,J,M,N,N1,N2,N3,N4,N5,N6)
  !------------------------------------------------
  !     SOLUTE TRANSPORT IN MACROPORES
  !
  call MacroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6)
  !
  call MacroporeSoluteDispersionM(M,N,N1,N2,N3,N4,N5,N6)

  call PrintInfo('end '//subname)
  end subroutine SoluteAdvDifusTranspM

! ----------------------------------------------------------------------
  subroutine MicroporeSoluteDiffusionM(I,J,M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer , intent(in) :: I,J,M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)

  real(r8) :: FaceArea,AreaOverDistance,Dispersivity
  character(len=*), parameter :: subname='MicroporeSoluteDiffusionM'
  real(r8) :: THETW1_vr(JZ,JY,JX)  !soil saturation   
  real(r8) :: VLWatMicPOA,VLWatMicPOB,VLWatMicPPA,VLWatMicPPB
  real(r8) :: VLWatMicP2A,VLWatMicP2B,VLWatMicP3A,VLWatMicP3B,VLWatMicP4A,VLWatMicP4B
  real(r8) :: trcsolc1(ids_beg:ids_end)
  real(r8) :: trcsolc2(ids_beg:ids_end)
  real(r8) :: CDOM_MicP1(idom_beg:idom_end,1:jcplx)
  real(r8) :: CDOM_MicP2(idom_beg:idom_end,1:jcplx)
  real(r8) :: SDifc(ids_beg:ids_end),SDifFlx(ids_beg:ids_end)
  real(r8) :: DISPN,DOMDifc(idom_beg:idom_end)
  real(r8) :: DLYR1,DLYR2,TORTL
  integer  :: K,ids,idg,idom

  call PrintInfo('beg '//subname)
  THETW1_vr(N3,N2,N1)=AZMAX1(safe_adb(VLWatMicPM_vr(M,N3,N2,N1),VLSoilMicP_vr(N3,N2,N1)))
  THETW1_vr(N6,N5,N4)=AZMAX1(safe_adb(VLWatMicPM_vr(M,N6,N5,N4),VLSoilMicP_vr(N6,N5,N4)))

  IF(THETW1_vr(N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1) .AND. VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1) &         !grid (N3,N2,N1) OK
    .AND. THETW1_vr(N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4) .AND. VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN !grid (N6,N5,N4) OK

    !------------------------------------------------
    !     MICROPORE CONCENTRATIONS FROM WATER-FILLED POROSITY
    !     IN CURRENT AND ADJACENT GRID CELLS
    !
    D9810: DO K=1,jcplx
      do idom=idom_beg,idom_end
        !source
        CDOM_MicP1(idom,K)=AZMAX1(DOM_MicP2_vr(idom,K,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
        !dest
        CDOM_MicP2(idom,K)=AZMAX1(DOM_MicP2_vr(idom,K,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
      enddo
    ENDDO D9810

    DO idg=idg_beg,idg_NH3-1
      trcsolc1(idg)=AZMAX1(trcs_solml2_vr(idg,N3,N2,N1)/VLWatMicPM_vr(M,N3,N2,N1))
    ENDDO

    VLWatMicP4A = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4,N3,N2,N1)
    IF(VLWatMicP4A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4,N3,N2,N1)/VLWatMicP4A),tracerSolc_max(ids_NH4))
      trcsolc1(idg_NH3) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3,N3,N2,N1)/VLWatMicP4A),tracerSolc_max(idg_NH3))
    ELSE
      trcsolc1(ids_NH4)=0.0_r8
      trcsolc1(idg_NH3)=0.0_r8
    ENDIF

    VLWatMicP3A = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3,N3,N2,N1)
    IF(VLWatMicP3A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NO3,N3,N2,N1)/VLWatMicP3A),tracerSolc_max(ids_NO3))
      trcsolc1(ids_NO2) = AZMAX1(trcs_solml2_vr(ids_NO2,N3,N2,N1)/VLWatMicP3A)      
    ELSE
      trcsolc1(ids_NO3)=0.0_r8
      trcsolc1(ids_NO2)=0.0_r8
    ENDIF

    VLWatMicP2A = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    IF(VLWatMicP2A.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4) = AZMAX1(trcs_solml2_vr(ids_H1PO4,N3,N2,N1)/VLWatMicP2A)
      trcsolc1(ids_H2PO4) = AZMAX1(trcs_solml2_vr(ids_H2PO4,N3,N2,N1)/VLWatMicP2A)
    ELSE
      trcsolc1(ids_H1PO4)=0.0_r8
      trcsolc1(ids_H2PO4)=0.0_r8
    ENDIF

    VLWatMicP4B = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4B,N3,N2,N1)
    IF(VLWatMicP4B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NH4B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4B,N3,N2,N1)/VLWatMicP4B),tracerSolc_max(ids_NH4B))
      trcsolc1(idg_NH3B) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3B,N3,N2,N1)/VLWatMicP4B),tracerSolc_max(idg_NH3B))
    ELSE
      trcsolc1(ids_NH4B)=0.0_r8
      trcsolc1(idg_NH3B)=0.0_r8
    ENDIF

    VLWatMicP3B = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3B,N3,N2,N1)
    IF(VLWatMicP3B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_NO3B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NO3B,N3,N2,N1)/VLWatMicP3B),tracerSolc_max(ids_NO3B))
      trcsolc1(ids_NO2B) = AZMAX1(trcs_solml2_vr(ids_NO2B,N3,N2,N1)/VLWatMicP3B)
    ELSE
      trcsolc1(ids_NO3B) = trcsolc1(ids_NO3)
      trcsolc1(ids_NO2B) = trcsolc1(ids_NO2)
    ENDIF

    VLWatMicP2B = VLWatMicPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
    IF(VLWatMicP2B.GT.ZEROS2(N2,N1))THEN
      trcsolc1(ids_H1PO4B) = AZMAX1(trcs_solml2_vr(ids_H1PO4B,N3,N2,N1)/VLWatMicP2B)
      trcsolc1(ids_H2PO4B) = AZMAX1(trcs_solml2_vr(ids_H2PO4B,N3,N2,N1)/VLWatMicP2B)
    ELSE
      trcsolc1(ids_H1PO4B)=trcsolc1(ids_H1PO4)
      trcsolc1(ids_H2PO4B)=trcsolc1(ids_H2PO4)
    ENDIF

    DO idg=idg_beg,idg_NH3-1
      trcsolc2(idg)=AZMAX1(trcs_solml2_vr(idg,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4))
    ENDDO

    IF(VLWatMicPMA_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      !gN/m3 H2O, the maximum solubility is about 9.7gN/100g water at 25oC    
      trcsolc2(idg_NH3) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3,N6,N5,N4)/VLWatMicPMA_vr(N6,N5,N4)),tracerSolc_max(idg_NH3))
      trcsolc2(ids_NH4) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4,N6,N5,N4)/VLWatMicPMA_vr(N6,N5,N4)),tracerSolc_max(ids_NH4))
    ELSE
      trcsolc2(idg_NH3) = 0.0_r8
      trcsolc2(ids_NH4) = 0.0_r8
    ENDIF

    VLWatMicPOA = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3,N6,N5,N4)
    IF(VLWatMicPOA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3) = AZMAX1(trcs_solml2_vr(ids_NO3,N6,N5,N4)/VLWatMicPOA)
      trcsolc2(ids_NO2) = AZMAX1(trcs_solml2_vr(ids_NO2,N6,N5,N4)/VLWatMicPOA)
    ELSE
      trcsolc2(ids_NO3) = 0.0_r8
      trcsolc2(ids_NO2) = 0.0_r8
    ENDIF

    VLWatMicPPA = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    IF(VLWatMicPPA.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4) = AZMAX1(trcs_solml2_vr(ids_H1PO4,N6,N5,N4)/VLWatMicPPA)
      trcsolc2(ids_H2PO4) = AZMAX1(trcs_solml2_vr(ids_H2PO4,N6,N5,N4)/VLWatMicPPA)
    ELSE
      trcsolc2(ids_H1PO4) = 0.0_r8
      trcsolc2(ids_H2PO4) = 0.0_r8
    ENDIF

    IF(VLWatMicPMB_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      trcsolc2(idg_NH3B) = AMIN1(AZMAX1(trcs_solml2_vr(idg_NH3B,N6,N5,N4)/VLWatMicPMB_vr(N6,N5,N4)),tracerSolc_max(idg_NH3B))
      trcsolc2(ids_NH4B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NH4B,N6,N5,N4)/VLWatMicPMB_vr(N6,N5,N4)),tracerSolc_max(ids_NH4B))
    ELSE
      trcsolc2(idg_NH3B) = trcsolc2(idg_NH3)
      trcsolc2(ids_NH4B) = trcsolc2(ids_NH4)
    ENDIF

    VLWatMicPOB = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3B,N6,N5,N4)
    IF(VLWatMicPOB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_NO3B) = AMIN1(AZMAX1(trcs_solml2_vr(ids_NO3B,N6,N5,N4)/VLWatMicPOB),tracerSolc_max(ids_NO3B))
      trcsolc2(ids_NO2B) = AZMAX1(trcs_solml2_vr(ids_NO2B,N6,N5,N4)/VLWatMicPOB)
    ELSE
      trcsolc2(ids_NO3B) = trcsolc2(ids_NO3)
      trcsolc2(ids_NO2B) = trcsolc2(ids_NO2)
    ENDIF

    VLWatMicPPB = VLWatMicPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    IF(VLWatMicPPB.GT.ZEROS2(N5,N4))THEN
      trcsolc2(ids_H1PO4B) = AZMAX1(trcs_solml2_vr(ids_H1PO4B,N6,N5,N4)/VLWatMicPPB)
      trcsolc2(ids_H2PO4B) = AZMAX1(trcs_solml2_vr(ids_H2PO4B,N6,N5,N4)/VLWatMicPPB)
    ELSE
      trcsolc2(ids_H1PO4B) = trcsolc2(ids_H1PO4)
      trcsolc2(ids_H2PO4B) = trcsolc2(ids_H2PO4)
    ENDIF
    !------------------------------------------------
    !     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MICROPORES
    !
    DLYR1 = AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2 = AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL = (TortMicPM_vr(M,N3,N2,N1)*DLYR1+TortMicPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    call GetSoilFaceGeometry(N,N1,N2,N3,N4,N5,N6,FaceArea,AreaOverDistance,Dispersivity)
    DISPN = Dispersivity*AMIN1(VFLWX,ABS(WaterFlow2MicPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/FaceArea))

    DO idom=idom_beg,idom_end
      DOMDifc(idom)=(DOM_diffusivitytscaledM_vr(idom,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    ENDDO

    DO ids=ids_beg,ids_end
      SDifc(ids)=(SoluteDifusivitytscaledM_vr(ids,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    ENDDO
    !
    D9805: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Difus_Mac2Micp_flxM(idom,K)=DOMDifc(idom)*(CDOM_MicP1(idom,K)-CDOM_MicP2(idom,K))        
      enddo
    ENDDO D9805

    DO idg=idg_beg,idg_NH3-1
      SDifFlx(idg)=SDifc(idg)*(trcsolc1(idg)-trcsolc2(idg))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      SDifFlx(ids)=SDifc(ids)*(trcsolc1(ids)-trcsolc2(ids))*AMIN1(trcs_VLN_vr(ids,N3,N2,N1),trcs_VLN_vr(ids,N6,N5,N4))      
    ENDDO
  ELSE
    D9905: DO K=1,jcplx
      DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9905
    SDifFlx(ids_beg:ids_end)=0._r8
  ENDIF
  !diffusion flux into (N6,N5,N4)
  DO ids=ids_beg,ids_end
    trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)=trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)+SDifFlx(ids)

    if(abs(trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4))>1.e20)then
      write(*,*)I*1000+J,M
      write(*,*)ids,N,N6,N5,N4,trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4),SDifFlx(ids),trcs_names(ids)
      write(*,*)SDifc(ids),trcsolc1(ids),trcsolc2(ids),AMIN1(trcs_VLN_vr(ids,N3,N2,N1),trcs_VLN_vr(ids,N6,N5,N4))   
      write(*,*)trcs_solml2_vr(idg,N3,N2,N1),VLWatMicPM_vr(M,N3,N2,N1),trcs_solml2_vr(idg,N6,N5,N4),VLWatMicPM_vr(M,N6,N5,N4)
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDDO

  DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MicpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4) = DOM_MicpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4)+DOM_Difus_Mac2Micp_flxM(idom,K)
    ENDDO  
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine MicroporeSoluteDiffusionM

! ----------------------------------------------------------------------
  subroutine MacroporeSoluteAdvectionM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: DOM_Adv2MacP_flxM(idom_beg:idom_end,1:jcplx)
  character(len=*), parameter :: subname='MacroporeSoluteAdvectionM'
  real(r8) ::  trcs_Adv2MacP_flxM(ids_beg:ids_end)
  integer  :: K,idg,ids,idom
  real(r8) :: VFLW

  call PrintInfo('beg '//subname)
  !     WaterFlow2MacPM_3D=water flux through soil macropore from watsub.f
  !
  IF(WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4).GT.0.0_r8)THEN
    !
    !     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM CURRENT TO
    !     ADJACENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
    !     OF WATER FLUX AND MACROPORE SOLUTE CONCENTRATIONS IN CURRENT
    !     GRID CELL
    !
    IF(VLWatMacPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/VLWatMacPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    !
    !     ACCOUNT FOR MACROPORE-MICROPORE EXCHANGE
    !
    IF(N.EQ.iVerticalDirection .AND. VLMacroPores_vr(N6,N5,N4).GT.VLWatMacPM_vr(M,N6,N5,N4))THEN
      D9800: DO K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_Adv2MacP_flxM(idom,K)=VFLW*AZMAX1((DOM_MacP2_vr(idom,K,N3,N2,N1) &
            -AZMIN1(DOM_Mac2MicPore_flxM_vr(idom,K,N3,N2,N1))))
        enddo
      ENDDO D9800

      DO idg=idg_beg,idg_NH3-1
         trcs_Adv2MacP_flxM(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N3,N2,N1)-AZMIN1(trcs_Mac2MicPore_flxM_vr(idg,N3,N2,N1)))
      ENDDO

      DO ids=ids_nuts_beg,ids_nuts_end
         trcs_Adv2MacP_flxM(ids)=VFLW*AZMAX1((trcs_soHml2_vr(ids,N3,N2,N1) &
          -AZMIN1(trcs_Mac2MicPore_flxM_vr(ids,N3,N2,N1)*trcs_VLN_vr(ids,N3,N2,N1)))) &
          *trcs_VLN_vr(ids,N6,N5,N4)
      ENDDO
!
!     OTHERWISE
!
    ELSE
      D9850: DO K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_Adv2MacP_flxM(idom,K)=VFLW*AZMAX1(DOM_MacP2_vr(idom,K,N3,N2,N1))
        enddo
      ENDDO D9850
      !exclude NH3 and NH3B
      DO idg=idg_beg,idg_NH3-1
         trcs_Adv2MacP_flxM(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N3,N2,N1))
      ENDDO

      DO ids=ids_nuts_beg,ids_nuts_end
         trcs_Adv2MacP_flxM(ids)=VFLW*AZMAX1(trcs_soHml2_vr(ids,N3,N2,N1))*trcs_VLN_vr(ids,N6,N5,N4)
      ENDDO
    ENDIF
!
!     IF MACROPORE WATER FLUX FROM 'WATSUB' IS FROM ADJACENT TO
!     CURRENT GRID CELL THEN CONVECTIVE TRANSPORT IS THE PRODUCT
!     OF WATER FLUX AND MACROPORE SOLUTE CONCENTRATIONS IN ADJACENT
!     GRID CELL
!
  ELSEIF(WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4).LT.0.0_r8)THEN
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    D9665: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MacP_flxM(idom,K)=VFLW*AZMAX1(DOM_MacP2_vr(idom,K,N6,N5,N4))
      enddo
    ENDDO D9665

    DO idg=idg_beg,idg_end-2
       trcs_Adv2MacP_flxM(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
       trcs_Adv2MacP_flxM(ids)=VFLW*AZMAX1(trcs_soHml2_vr(ids,N6,N5,N4))*trcs_VLN_vr(ids,N6,N5,N4)
    ENDDO
  ELSE
!
!     NO MACROPORE FLUX
!
    D9795: DO K=1,jcplx
      DOM_Adv2MacP_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9795
     trcs_Adv2MacP_flxM(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_MacpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)= trcs_MacpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)+trcs_Adv2MacP_flxM(ids)
  ENDDO

  DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MacpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4) =DOM_MacpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4)+ DOM_Adv2MacP_flxM(idom,K)
    ENDDO
  ENDDO      
  call PrintInfo('end '//subname)
  end subroutine MacroporeSoluteAdvectionM

! ----------------------------------------------------------------------
  subroutine MacroporeSoluteDispersionM(M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6

  real(r8) :: FaceArea,AreaOverDistance,Dispersivity
  character(len=*), parameter :: subname='MacroporeSoluteDispersionM'
  real(r8) :: DOM_Difus_Macp_flxM(idom_beg:idom_end,1:jcplx)
  real(r8) :: VOLH2A,VOLH2B,VOLH3A,VOLH3B,VOLH4A,VOLH4B
  real(r8) :: VOLHMA,VOLHMB,VOLHOA,VOLHOB,VOLHPA,VOLHPB
  real(r8) :: trcs_coH1(ids_beg:ids_end)
  real(r8) :: trcs_coH2(ids_beg:ids_end)
  real(r8) :: SDifc(ids_beg:ids_end),TORTL
  real(r8) :: DOMDifc(idom_beg:idom_end)
  real(r8) :: DISPN,DLYR1,DLYR2
  real(r8) :: SDifHFlx(ids_beg:ids_end)
  integer  :: K,ids,idg,idom
  real(r8) :: CDOM_MacP1(idom_beg:idom_end,1:jcplx)
  real(r8) :: CDOM_MacP2_vr(idom_beg:idom_end,1:jcplx)

  call PrintInfo('beg '//subname)

  IF(VLWatMacPM_vr(M,N3,N2,N1).GT.SoilWatAirDry_vr(N3,N2,N1)*VLMacroPores_vr(N3,N2,N1) &
    .AND. VLWatMacPM_vr(M,N6,N5,N4).GT.SoilWatAirDry_vr(N6,N5,N4)*VLMacroPores_vr(N6,N5,N4))THEN
    !------------------------------------------------------------------------------------
    !     MACROPORE CONCENTRATIONS IN CURRENT AND ADJACENT GRID CELLS
    !
    D9790: DO K=1,jcplx
      do idom=idom_beg,idom_end
        CDOM_MacP1(idom,K)    = AZMAX1(DOM_MacP2_vr(idom,K,N3,N2,N1)/VLWatMacPM_vr(M,N3,N2,N1))
        CDOM_MacP2_vr(idom,K) = AZMAX1(DOM_MacP2_vr(idom,K,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4))
      enddo
    ENDDO D9790

    !exclude NH3 and NH3B
    DO idg=idg_beg,idg_NH3-1
      trcs_coH1(idg)=AZMAX1(trcs_soHml2_vr(idg,N3,N2,N1)/VLWatMacPM_vr(M,N3,N2,N1))
    ENDDO

    VOLH4A = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4,N3,N2,N1)
    IF(VOLH4A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4)=AZMAX1(trcs_soHml2_vr(ids_NH4,N3,N2,N1)/VOLH4A)
      trcs_coH1(idg_NH3)=AZMAX1(trcs_soHml2_vr(idg_NH3,N3,N2,N1)/VOLH4A)
    ELSE
      trcs_coH1(ids_NH4)=0.0_r8
      trcs_coH1(idg_NH3)=0.0_r8
    ENDIF

    VOLH3A = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3,N3,N2,N1)
    IF(VOLH3A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3)=AZMAX1(trcs_soHml2_vr(ids_NO3,N3,N2,N1)/VOLH3A)
      trcs_coH1(ids_NO2)=AZMAX1(trcs_soHml2_vr(ids_NO2,N3,N2,N1)/VOLH3A)
    ELSE
      trcs_coH1(ids_NO3)=0.0_r8
      trcs_coH1(ids_NO2)=0.0_r8
    ENDIF

    VOLH2A = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4,N3,N2,N1)
    IF(VOLH2A.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4)=AZMAX1(trcs_soHml2_vr(ids_H1PO4,N3,N2,N1)/VOLH2A)
      trcs_coH1(ids_H2PO4)=AZMAX1(trcs_soHml2_vr(ids_H2PO4,N3,N2,N1)/VOLH2A)
    ELSE
      trcs_coH1(ids_H1PO4)=0.0_r8
      trcs_coH1(ids_H2PO4)=0.0_r8
    ENDIF

    VOLH4B = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NH4B,N3,N2,N1)
    IF(VOLH4B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NH4B)=AZMAX1(trcs_soHml2_vr(ids_NH4B,N3,N2,N1)/VOLH4B)
      trcs_coH1(idg_NH3B)=AZMAX1(trcs_soHml2_vr(idg_NH3B,N3,N2,N1)/VOLH4B)
    ELSE
      trcs_coH1(ids_NH4B)=trcs_coH1(ids_NH4)
      trcs_coH1(idg_NH3B)=trcs_coH1(idg_NH3)
    ENDIF

    VOLH3B = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_NO3B,N3,N2,N1)
    IF(VOLH3B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_NO3B)=AZMAX1(trcs_soHml2_vr(ids_NO3B,N3,N2,N1)/VOLH3B)
      trcs_coH1(ids_NO2B)=AZMAX1(trcs_soHml2_vr(ids_NO2B,N3,N2,N1)/VOLH3B)
    ELSE
      trcs_coH1(ids_NO3B)=trcs_coH1(ids_NO3)
      trcs_coH1(ids_NO2B)=trcs_coH1(ids_NO2)
    ENDIF

    VOLH2B = VLWatMacPM_vr(M,N3,N2,N1)*trcs_VLN_vr(ids_H1PO4B,N3,N2,N1)
    IF(VOLH2B.GT.ZEROS2(N2,N1))THEN
      trcs_coH1(ids_H1PO4B)=AZMAX1(trcs_soHml2_vr(ids_H1PO4B,N3,N2,N1)/VOLH2B)
      trcs_coH1(ids_H2PO4B)=AZMAX1(trcs_soHml2_vr(ids_H2PO4B,N3,N2,N1)/VOLH2B)
    ELSE
      trcs_coH1(ids_H1PO4B)=trcs_coH1(ids_H1PO4)
      trcs_coH1(ids_H2PO4B)=trcs_coH1(ids_H2PO4)
    ENDIF

   !excldue NH3 and NH3B
    DO idg=idg_beg,idg_NH3-1
      trcs_coH2(idg)=AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4))
    ENDDO

    VOLHMA=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4,N6,N5,N4)
    IF(VOLHMA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4)=AZMAX1(trcs_soHml2_vr(ids_NH4,N6,N5,N4)/VOLHMA)
      trcs_coH2(idg_NH3)=AZMAX1(trcs_soHml2_vr(idg_NH3,N6,N5,N4)/VOLHMA)
    ELSE
      trcs_coH2(ids_NH4)=0.0_r8
      trcs_coH2(idg_NH3)=0.0_r8
    ENDIF

    VOLHOA = VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3,N6,N5,N4)
    IF(VOLHOA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3)=AZMAX1(trcs_soHml2_vr(ids_NO3,N6,N5,N4)/VOLHOA)
      trcs_coH2(ids_NO2)=AZMAX1(trcs_soHml2_vr(ids_NO2,N6,N5,N4)/VOLHOA)
    ELSE
      trcs_coH2(ids_NO3)=0.0_r8
      trcs_coH2(ids_NO2)=0.0_r8
    ENDIF

    VOLHPA = VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4,N6,N5,N4)
    IF(VOLHPA.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4)=AZMAX1(trcs_soHml2_vr(ids_H1PO4,N6,N5,N4)/VOLHPA)
      trcs_coH2(ids_H2PO4)=AZMAX1(trcs_soHml2_vr(ids_H2PO4,N6,N5,N4)/VOLHPA)
    ELSE
      trcs_coH2(ids_H1PO4)=0.0_r8
      trcs_coH2(ids_H2PO4)=0.0_r8
    ENDIF

    VOLHMB=VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NH4B,N6,N5,N4)
    IF(VOLHMB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NH4B)=AZMAX1(trcs_soHml2_vr(ids_NH4B,N6,N5,N4)/VOLHMB)
      trcs_coH2(idg_NH3B)=AZMAX1(trcs_soHml2_vr(idg_NH3B,N6,N5,N4)/VOLHMB)
    ELSE
      trcs_coH2(ids_NH4B)=trcs_coH2(ids_NH4)
      trcs_coH2(idg_NH3B)=trcs_coH2(idg_NH3)
    ENDIF

    VOLHOB = VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_NO3B,N6,N5,N4)
    IF(VOLHOB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_NO3B)=AZMAX1(trcs_soHml2_vr(ids_NO3B,N6,N5,N4)/VOLHOB)
      trcs_coH2(ids_NO2B)=AZMAX1(trcs_soHml2_vr(ids_NO2B,N6,N5,N4)/VOLHOB)
    ELSE
      trcs_coH2(ids_NO3B)=trcs_coH2(ids_NO3)
      trcs_coH2(ids_NO2B)=trcs_coH2(ids_NO2)
    ENDIF

    VOLHPB = VLWatMacPM_vr(M,N6,N5,N4)*trcs_VLN_vr(ids_H1PO4B,N6,N5,N4)
    IF(VOLHPB.GT.ZEROS2(N5,N4))THEN
      trcs_coH2(ids_H1PO4B)=AZMAX1(trcs_soHml2_vr(ids_H1PO4B,N6,N5,N4)/VOLHPB)
      trcs_coH2(ids_H2PO4B)=AZMAX1(trcs_soHml2_vr(ids_H2PO4B,N6,N5,N4)/VOLHPB)
    ELSE
      trcs_coH2(ids_H1PO4B)=trcs_coH2(ids_H1PO4)
      trcs_coH2(ids_H2PO4B)=trcs_coH2(ids_H2PO4)
    ENDIF

    !------------------------------------------------------------------------------------
    !
    !     DIFFUSIVITIES IN CURRENT AND ADJACENT GRID CELL MACROPORES
    !
    DLYR1 = AMAX1(ZERO2,DLYR_3D(N,N3,N2,N1))
    DLYR2 = AMAX1(ZERO2,DLYR_3D(N,N6,N5,N4))
    TORTL = (TortMacPM_vr(M,N3,N2,N1)*DLYR1+TortMacPM_vr(M,N6,N5,N4)*DLYR2)/(DLYR1+DLYR2)
    call GetSoilFaceGeometry(N,N1,N2,N3,N4,N5,N6,FaceArea,AreaOverDistance,Dispersivity)
    DISPN = Dispersivity*AMIN1(VFLWX,ABS(WaterFlow2MacPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/FaceArea))

    DO idom=idom_beg,idom_end
      DOMDifc(idom)=(DOM_diffusivitytscaledM_vr(idom,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    ENDDO

    DO ids=ids_beg,ids_end
      SDifc(ids)=(SoluteDifusivitytscaledM_vr(ids,N6,N5,N4)*TORTL+DISPN)*AreaOverDistance
    ENDDO
    !
    !     DIFFUSIVE FLUXES BETWEEN CURRENT AND ADJACENT GRID CELL
    !     MACROPORES
    !
    D9785: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Difus_Macp_flxM(idom,K)=DOMDifc(idom)*(CDOM_MacP1(idom,K)-CDOM_MacP2_vr(idom,K))
      enddo
    ENDDO D9785

    ! exclude NH3 and NH3B
    DO idg=idg_beg,idg_NH3-1
      SDifHFlx(idg)=SDifc(idg)*(trcs_coH1(idg)-trcs_coH2(idg))
    ENDDO

    DO ids=ids_nuts_beg,ids_end
      SDifHFlx(ids)=SDifc(ids)*(trcs_coH1(ids)-trcs_coH2(ids)) &
        *AMIN1(trcs_VLN_vr(ids,N3,N2,N1),trcs_VLN_vr(ids,N6,N5,N4))
    ENDDO
  ELSE
    D9780: DO K=1,jcplx
      DOM_Difus_Macp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9780
    SDifHFlx(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO ids = ids_beg,ids_end
    trcs_MacpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)=trcs_MacpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)+SDifHFlx(ids)
  ENDDO

  DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MacpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4) = DOM_MacpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4)+DOM_Difus_Macp_flxM(idom,K)
    enddo
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine MacroporeSoluteDispersionM

!----------------------------------------------------------------------

  subroutine MicroporeSoluteAdvectionM(I,J,M,N,N1,N2,N3,N4,N5,N6)
  implicit none
  integer, intent(in) :: I,J,M,N,N1,N2,N3,N4,N5,N6

  character(len=*), parameter :: subname='MicroporeSoluteAdvectionM'
  real(r8) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)
  real(r8) :: trcs_Adv2MicP_flx(ids_beg:ids_end)
  integer :: K,ids,idom
  real(r8) :: VFLW
   
  call PrintInfo('beg '//subname)
  !advection from (N3,N2,N1) to (N6,N5,N4)
  IF(WaterFlow2MicPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4).GT.0.0_r8)THEN
    IF(VLWatMicPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,SoilFaceSlot(N,N3,N6),N5,N4)/VLWatMicPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,N3,N2,N1))
      enddo
    ENDDO

    DO ids=ids_beg,ids_end
      trcs_Adv2MicP_flx(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,N3,N2,N1))
    ENDDO
    !advection from (N6,N5,N4) to (N3,N2,N1)
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
    D9815: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,N6,N5,N4))
      enddo
    ENDDO D9815
    DO ids=ids_beg,ids_end
      trcs_Adv2MicP_flx(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,N6,N5,N4))
    ENDDO
  ENDIF
  !add advection flux
  DO ids=ids_beg,ids_end
    trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)=trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)+trcs_Adv2MicP_flx(ids)

    if(abs(trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4))>1.e20)then
      write(*,*)I*1000+J,ids,N,N6,N5,N4,trcs_Adv2MicP_flx(ids)
      write(*,*)VFLW,trcs_solml2_vr(ids,N6,N5,N4),trcs_solml2_vr(ids,N3,N2,N1),trcs_names(ids)
      write(*,*)trcs_solml_vr(ids,N6,N5,N4),trcs_solml_vr(ids,N3,N2,N1)
      call endrun(trim(mod_filename)//' at line',__LINE__)
    endif
  ENDDO

  DO K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_MicpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4) = DOM_MicpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4)+DOM_Adv2MicP_flx(idom,K)
    enddo
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine MicroporeSoluteAdvectionM

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteExchM(M,N1,N2,N3)
  !
  !Description:
  !tracer exchange between micropore-macropore in each soil grid
  !
  implicit none
  integer, intent(in) :: M,N1,N2,N3
  character(len=*), parameter :: subname='MicMacPoresSoluteExchM'
  integer :: K,ids,idg,idom

  call PrintInfo('beg '//subname)

  do idg=idg_beg,idg_NH3-1
    trcg_VLWatMicP_vr(idg,N3,N2,N1)=VLWatMicPM_vr(M,N3,N2,N1)*GasSolbility_vr(idg,N3,N2,N1)
  enddo

  trcg_VLWatMicP_vr(idg_NH3,N3,N2,N1)=VLWatMicPMA_vr(N3,N2,N1)*GasSolbility_vr(idg_NH3,N3,N2,N1)
  trcg_VLWatMicP_vr(idg_NH3B,N3,N2,N1)=VLWatMicPMB_vr(N3,N2,N1)*GasSolbility_vr(idg_NH3,N3,N2,N1)

  call MicMacPoresSoluteAdvExchM(M,N1,N2,N3)
 
  call MicMacPoresSoluteDifExchM(M,N1,N2,N3)
  
  call PrintInfo('end '//subname)

  end subroutine MicMacPoresSoluteExchM

! ----------------------------------------------------------------------
  subroutine MicMacPoresSoluteAdvExchM(M,N4,N5,N6)
  !
  !Description
  !Do micropore-macropore tracer exchange by upstream advection
  implicit none
  integer, intent(in) :: M,N4,N5,N6

  character(len=*), parameter :: subname='MicMacPoresSoluteAdvExchM'
  real(r8) :: DOM_Adv2MicP_flx(idom_beg:idom_end,1:jcplx)
  real(r8) :: VFLW    !fraction in flow
  real(r8) :: trcs_Adv2MicP_flx(ids_beg:ids_end)
  integer :: K,idg,ids,idom

  call PrintInfo('beg '//subname)
!
!     MACROPORE TO MICROPORE TRANSFER
!
  IF(FWatExMacP2MicPM_vr(M,N6,N5,N4).GT.0.0_r8)THEN
    IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMAX1(AMIN1(VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMacPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=VFLWX
    ENDIF

    D9970: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MacP2_vr(idom,K,N6,N5,N4))
      enddo
    ENDDO D9970

    DO idg=idg_beg,idg_NH3-1
      trcs_Adv2MicP_flx(idg)=VFLW*AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4))
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_Adv2MicP_flx(ids)=VFLW*AZMAX1(trcs_soHml2_vr(ids,N6,N5,N4))*trcs_VLN_vr(ids,N6,N5,N4)
    ENDDO
!
!     MICROPORE TO MACROPORE TRANSFER (<0)
!
  ELSEIF(FWatExMacP2MicPM_vr(M,N6,N5,N4).LT.0.0_r8)THEN
    IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
      VFLW=AZMIN1(AMAX1(-VFLWX,FWatExMacP2MicPM_vr(M,N6,N5,N4)/VLWatMicPM_vr(M,N6,N5,N4)))
    ELSE
      VFLW=-VFLWX
    ENDIF

    D9965: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Adv2MicP_flx(idom,K)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,N6,N5,N4))
      enddo
    ENDDO D9965
    !
    DO idg=idg_beg,idg_NH3-1
      trcs_Adv2MicP_flx(idg)=VFLW*AZMAX1(trcs_solml2_vr(idg,N6,N5,N4))
    ENDDO

    !include NH3 and NH3B below
    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_Adv2MicP_flx(ids)=VFLW*AZMAX1(trcs_solml2_vr(ids,N6,N5,N4))*trcs_VLN_vr(ids,N6,N5,N4)
    ENDDO
    !
    !     NO MACROPORE TO MICROPORE TRANSFER
    !
  ELSE
    D9960: DO K=1,jcplx
      DOM_Adv2MicP_flx(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9960
    trcs_Adv2MicP_flx(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)=trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)+trcs_Adv2MicP_flx(ids)
  ENDDO
  !
  !     TOTAL CONVECTIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
  ! > 0, macropore to micropore
  ! < 0, micropore to macropore
  DO  K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)=DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)+DOM_Adv2MicP_flx(idom,K)
    enddo
  enddo
  call PrintInfo('end '//subname)
  end subroutine MicMacPoresSoluteAdvExchM

! ----------------------------------------------------------------------

  subroutine MicMacPoresSoluteDifExchM(M,N4,N5,N6)
  !
  !Description
  !Do macropore-micropore tracer exchange by diffusion
  !
  implicit none
  integer, intent(in) :: M,N4,N5,N6

  character(len=*), parameter :: subname='MicMacPoresSoluteDifExchM'
  real(r8) :: trcs_Difus_Mac2Micp_flxM(ids_beg:ids_end)
  real(r8) :: VLWatMacPS,VOLWT
  integer  :: K,ids,idg,idom
  real(r8) :: DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,1:jcplx)

  call PrintInfo('beg '//subname)
  IF(VLWatMacPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN    
    VLWatMacPS = AMIN1(XFRS*VGeomLayer_vr(N6,N5,N4),VLWatMacPM_vr(M,N6,N5,N4))
    VOLWT      = VLWatMicPM_vr(M,N6,N5,N4)+VLWatMacPS

    D9955: DO K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_Difus_Mac2Micp_flxM(idom,K)=dts_HeatWatTP*(AZMAX1(DOM_MacP2_vr(idom,K,N6,N5,N4)) &
          *VLWatMicPM_vr(M,N6,N5,N4)-AZMAX1(DOM_MicP2_vr(idom,K,N6,N5,N4))*VLWatMacPS)/VOLWT
      enddo
    ENDDO D9955

    DO idg=idg_beg,idg_NH3-1
      trcs_Difus_Mac2Micp_flxM(idg)=dts_HeatWatTP*(AZMAX1(trcs_soHml2_vr(idg,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcs_solml2_vr(idg,N6,N5,N4))*VLWatMacPS)/VOLWT
    ENDDO

    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_Difus_Mac2Micp_flxM(ids)=dts_HeatWatTP*(AZMAX1(trcs_soHml2_vr(ids,N6,N5,N4))*VLWatMicPM_vr(M,N6,N5,N4) &
        -AZMAX1(trcs_solml2_vr(ids,N6,N5,N4))*VLWatMacPS)/VOLWT &
        *trcs_VLN_vr(ids,N6,N5,N4)
    ENDDO
  ELSE
    D9975: DO K=1,jcplx
      DOM_Difus_Mac2Micp_flxM(idom_beg:idom_end,K)=0.0_r8
    ENDDO D9975
    trcs_Difus_Mac2Micp_flxM(ids_beg:ids_end)=0.0_r8
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)=trcs_Mac2MicPore_flxM_vr(ids,N6,N5,N4)+trcs_Difus_Mac2Micp_flxM(ids)
  ENDDO
  !
  ! TOTAL CONVECTIVE +DIFFUSIVE TRANSFER BETWEEN MACROPOES AND MICROPORES
  !
  DO  K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)=DOM_Mac2MicPore_flxM_vr(idom,K,N6,N5,N4)+DOM_Difus_Mac2Micp_flxM(idom,K)
    enddo
  enddo
  call PrintInfo('end '//subname)
  end subroutine MicMacPoresSoluteDifExchM

! ----------------------------------------------------------------------
  subroutine BubbleEffluxM(I,J,M,NHE,NHW,NVS,NVN)

  !
  !bubbling are diagnosed from the first unfrozen layer, upward, and deposited to the first layer with gas phase 
  implicit none
  integer, intent(in) :: I,J,M,NHE,NHW,NVS,NVN

  integer  :: N1,N2,N3
  logical  :: iFalseEbu                        !bubbling flag

  character(len=*), parameter :: subname='BubbleEffluxM'  
  real(r8) :: THETW1
  real(r8) :: BubbleEffluxStep(idg_beg:idg_end) !negative gas mass routed in this substep
  real(r8) :: GasMassSolubility(idg_beg:idg_end)
  real(r8) :: trcg_VOLG(idg_beg:idg_end)   !gas concentation corresponding to each volatile [mol d-2]
  integer  :: idg
  real(r8) :: VTATM,VTGAS,DVTGAS
  real(r8) :: FracEbu                      !fraction removed through ebullition
  real(r8) :: dPond
  integer  :: LG
  real(r8) :: ThetaIce=0.99_r8
  !
  call PrintInfo('beg '//subname)
  LG=0
  DO N1=NHW,NHE
    DO N2=NVN,NVS
      iFalseEbu = .false.
      BubbleEffluxStep=0._r8
      trcg_Ebu_flxM_vr(idg_beg:idg_end,NUM_col(N2,N1):NL_col(N2,N1),N2,N1)=0.0_r8      

      DO N3=NLF_col(N2,N1),NUM_col(N2,N1),-1  !sweep upward assuming bubbles go upward
        !
        THETW1=AZMAX1(safe_adb(VLWatMicPM_vr(M,N3,N2,N1),VLSoilMicP_vr(N3,N2,N1)))
        !
        IF(THETW1.GT.SoilWatAirDry_vr(N3,N2,N1)     & ! has significant water 
          .AND. .not.iFalseEbu)THEN                 ! ebullition is allowed

          dPond=AZMAX1(CumDepz2LayBottom_vr(N3,N2,N1)-WaterTableDepZi_col(N2,N1))*10._r8

          do idg=idg_beg,idg_NH3
            GasMassSolubility(idg) =MolecularWeight(idg)*GasSolbility_vr(idg,N3,N2,N1)  !conver into carbon g /mol
          enddo        
          GasMassSolubility(idg_NH3B)=GasMassSolubility(idg_NH3)
          !
          !     GASEOUS EQUIVALENT PARTIAL CONCENTRATIONS
          !
          DO idg=idg_beg,idg_end
            trcg_VOLG(idg)=AZMAX1(trcs_solml2_vr(idg,N3,N2,N1))/GasMassSolubility(idg)  !mol/d2 
          ENDDO 
          !
          !     GASEOUS EQUIVALENT ATMOSPHERIC CONCENTRATION
          !
          !     VTATM=molar gas concentration at atmospheric pressure
          !     VTGAS=total molar gas concentration, should the gaseous phase be included as well?     
          !    
          VTATM = AZMAX1((PBOT_col(N2,N1)+dPond)*VLWatMicPM_vr(M,N3,N2,N1)/(RGasC*TKS_vr(N3,N2,N1)))*1.E3_r8  !mol gas/ d2
          VTGAS = sum(trcg_VOLG(idg_beg:idg_end))

          !
          ! PROPORTIONAL REMOVAL OF EXCESS AQUEOUS GASES
          !
          IF(VTGAS.GT.VTATM)THEN
            DVTGAS  = 0.5_r8*(VTATM-VTGAS)         !<0, bubbling occurs
            FracEbu = AZMIN1(DVTGAS/VTGAS)
            DO idg=idg_beg,idg_end
              trcg_Ebu_flxM_vr(idg,N3,N2,N1)=FracEbu*AZMAX1(trcs_solml2_vr(idg,N3,N2,N1))  

              trcs_solml2_vr(idg,N3,N2,N1)=trcs_solml2_vr(idg,N3,N2,N1)+trcg_Ebu_flxM_vr(idg,N3,N2,N1)
              if(trcs_solml2_vr(idg,N3,N2,N1)<0._r8)then
                write(iulog,*)(I*1000+J)*100+M,trcs_names(idg),N3,trcs_solml2_vr(idg,N3,N2,N1),trcg_Ebu_flxM_vr(idg,N3,N2,N1)
                call endrun(trim(mod_filename)//' at line',__LINE__)
              endif

              if(N3==NUM_col(N2,N1) .or. dPond>0._r8)then !surface layer or below water table
                ! Accumulate current-substep bubbles pending upward routing.
                if(idg==idg_NH3B)then
                  BubbleEffluxStep(idg_NH3) = BubbleEffluxStep(idg_NH3)+trcg_Ebu_flxM_vr(idg,N3,N2,N1)
                else
                  BubbleEffluxStep(idg) = BubbleEffluxStep(idg)+trcg_Ebu_flxM_vr(idg,N3,N2,N1)
                endif
              else
                !add to the gas phase, which will then go through gaseous transport, triggering bubbling
                if(idg/=idg_NH3B)then
                  trcg_gasml2_vr(idg,N3,N2,N1)=trcg_gasml2_vr(idg,N3,N2,N1)-trcg_Ebu_flxM_vr(idg,N3,N2,N1)
                  if(AZERO(trcg_gasml2_vr(idg,N3,N2,N1))<0._r8)then
                    write(iulog,*)I*1000+J,N3,trcs_names(idg),trcg_gasml2_vr(idg,N3,N2,N1),trcg_gasml_vr(idg,N3,N2,N1),trcg_Ebu_flxM_vr(idg,N3,N2,N1)
                    call endrun(trim(mod_filename)//' at line',__LINE__)
                  endif
                else
                  trcg_gasml2_vr(idg_NH3,N3,N2,N1)=trcg_gasml2_vr(idg_NH3,N3,N2,N1)-trcg_Ebu_flxM_vr(idg,N3,N2,N1)    
                  if(AZERO(trcg_gasml2_vr(idg_NH3,N3,N2,N1))<0._r8)then
                    write(iulog,*)I*1000+J,M,N3,trcg_gasml2_vr(idg_NH3,N3,N2,N1),trcg_gasml_vr(idg_NH3,N3,N2,N1),trcg_Ebu_flxM_vr(idg,N3,N2,N1)
                    call endrun(trim(mod_filename)//' at line',__LINE__)
                  endif             
                endif                
              endif              
            ENDDO
            RO2AquaSourcePrev_vr(N3,N2,N1) = RO2AquaSourcePrev_vr(N3,N2,N1)+trcg_Ebu_flxM_vr(idg_O2,N3,N2,N1)
            RCH4PhysexchPrev_vr(N3,N2,N1)  = RCH4PhysexchPrev_vr(N3,N2,N1)+trcg_Ebu_flxM_vr(idg_CH4,N3,N2,N1)
          ENDIF
        ELSE
          ! Modeling approximation: a dry layer encountered in the upward scan
          ! suppresses ebullition in all layers above it. Its gas-filled pores
          ! can accommodate gas via the separate dissolution/diffusion pathways,
          ! potentially delaying efflux; this flag itself transfers no upper-layer gas.
          ! TODO: Revisit whether those pathways justify suppressing local bubble
          ! formation in a supersaturated wet layer above a dry layer.
          iFalseEbu=.true.
          ! Intercept only bubbles routed in this substep; earlier emissions
          ! have already left the column and remain in the hourly diagnostic.
          DO idg=idg_beg,idg_end
            !add to the gas phase, which will then go through gaseous transport, triggering bubbling
            if(BubbleEffluxStep(idg)<0._r8)then
              if(idg/=idg_NH3B)then
                trcg_gasml2_vr(idg,N3,N2,N1)=trcg_gasml2_vr(idg,N3,N2,N1)-BubbleEffluxStep(idg)
                if(AZERO(trcg_gasml2_vr(idg,N3,N2,N1))<0._r8)then
                  write(iulog,*)I*1000+J,N3,trcs_names(idg),trcg_gasml2_vr(idg,N3,N2,N1),trcg_gasml_vr(idg,N3,N2,N1),trcg_Ebu_flxM_vr(idg,N3,N2,N1)
                  call endrun(trim(mod_filename)//' at line',__LINE__)
                endif
              else
                trcg_gasml2_vr(idg_NH3,N3,N2,N1)=trcg_gasml2_vr(idg_NH3,N3,N2,N1)-BubbleEffluxStep(idg)
                if(AZERO(trcg_gasml2_vr(idg_NH3,N3,N2,N1))<0._r8)then
                  write(iulog,*)I*1000+J,M,N3,trcg_gasml2_vr(idg_NH3,N3,N2,N1),trcg_gasml_vr(idg_NH3,N3,N2,N1),trcg_Ebu_flxM_vr(idg,N3,N2,N1)
                  call endrun(trim(mod_filename)//' at line',__LINE__)
                endif             
              endif
              BubbleEffluxStep(idg)=0._r8
            endif 
          ENDDO      
        ENDIF
      ENDDO
      ! Only bubbles surviving the current upward scan contribute to emissions.
      ! Band NH3 was merged into the shared NH3 gas tracer above.
      trcg_ebu_flx_col(idg_beg:idg_NH3,N2,N1)=trcg_ebu_flx_col(idg_beg:idg_NH3,N2,N1) &
        +BubbleEffluxStep(idg_beg:idg_NH3)
    ENDDO
  ENDDO  
  
  call PrintInfo('end '//subname)
  end subroutine BubbleEffluxM
! ----------------------------------------------------------------------
  subroutine TracerThruRofSnowXYM(M,N,N1,N2,N4,N5,N4B,N5B)
  !
  !Tracer flux through surface runoff and 
  !snow drift
  implicit none
  integer, intent(in) :: M,N
  integer, intent(in) :: N1,N2,N4,N5,N4B,N5B

  character(len=*), parameter :: subname='TracerThruRofSnowXYM'
  integer :: NN,K,idg,ids,idom

  call PrintInfo('beg '//subname)
  DO NN=1,2
    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_SurfRunoff_flxM(idom,K,N2,N1)=DOM_SurfRunoff_flxM(idom,K,N2,N1)+DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,N2,N1)
      enddo
    enddo

    DO idg=idg_beg,idg_NH3
      trcg_SurfRunoff_flxM(idg,N2,N1)=trcg_SurfRunoff_flxM(idg,N2,N1)+trcg_SurRof_flxM_2DH(idg,N,NN,N2,N1)
    ENDDO
    
    DO ids=ids_nut_beg,ids_nuts_end
      trcn_SurfRunoff_flxM(ids,N2,N1)=trcn_SurfRunoff_flxM(ids,N2,N1)+trcn_SurRof_flxM_2DH(ids,N,NN,N2,N1)
    ENDDO
    !this includes boundary fluxes
    IF(IFLBM_2DH(M,N,NN,N5,N4).EQ.0)THEN
      D9551: DO  K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_SurfRunoff_flxM(idom,K,N2,N1)=DOM_SurfRunoff_flxM(idom,K,N2,N1)-DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,N5,N4)
        enddo
      enddo D9551

      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flxM(idg,N2,N1)=trcg_SurfRunoff_flxM(idg,N2,N1)-trcg_SurRof_flxM_2DH(idg,N,NN,N5,N4)
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flxM(ids,N2,N1)=trcn_SurfRunoff_flxM(ids,N2,N1)-trcn_SurRof_flxM_2DH(ids,N,NN,N5,N4)
      ENDDO
    ENDIF

    IF(N4B.GT.0 .AND. N5B.GT.0 .AND. NN.EQ.iFront)THEN
      DO  K=1,jcplx
        do idom=idom_beg,idom_end
          DOM_SurfRunoff_flxM(idom,K,N2,N1)=DOM_SurfRunoff_flxM(idom,K,N2,N1)-DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,N5B,N4B)
        enddo
      enddo
      DO idg=idg_beg,idg_NH3
        trcg_SurfRunoff_flxM(idg,N2,N1)=trcg_SurfRunoff_flxM(idg,N2,N1)-trcg_SurRof_flxM_2DH(idg,N,NN,N5B,N4B)
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_SurfRunoff_flxM(ids,N2,N1)=trcn_SurfRunoff_flxM(ids,N2,N1)-trcn_SurRof_flxM_2DH(ids,N,NN,N5B,N4B)
      ENDDO
    ENDIF
  
    !---------------
    !     NET OVERLAND SOLUTE FLUX IN SNOW
    !
    DO idg=idg_beg,idg_NH3
      trcg_SnowDrift_flxM(idg,N2,N1)=trcg_SnowDrift_flxM(idg,N2,N1)+trcg_SnowDrift_flxM_2DH(idg,N,NN,N2,N1)
    ENDDO

    do ids=ids_nut_beg,ids_nuts_end
      trcn_SnowDrift_flxM(ids,N2,N1)=trcn_SnowDrift_flxM(ids,N2,N1)+trcn_SnowDrift_flxM_2DH(ids,N,NN,N2,N1)
    enddo

    if(IFLBSM_2DH(M,N,NN,N5,N4).EQ.0)then
      DO idg=idg_beg,idg_NH3
        trcg_SnowDrift_flxM(idg,N2,N1)=trcg_SnowDrift_flxM(idg,N2,N1)-trcg_SnowDrift_flxM_2DH(idg,N,NN,N5,N4)
      ENDDO

      do ids=ids_nut_beg,ids_nuts_end
        trcn_SnowDrift_flxM(ids,N2,N1)=trcn_SnowDrift_flxM(ids,N2,N1)-trcn_SnowDrift_flxM_2DH(ids,N,NN,N5,N4)
      enddo
    endif

    IF(N4B.GT.0 .AND. N5B.GT.0 .AND. NN.EQ.iFront)THEN
      DO idg=idg_beg,idg_NH3
        trcg_SnowDrift_flxM(idg,N2,N1)=trcg_SnowDrift_flxM(idg,N2,N1)-trcg_SnowDrift_flxM_2DH(idg,N,NN,N5B,N4B)
      ENDDO

      do ids=ids_nut_beg,ids_nuts_end
        trcn_SnowDrift_flxM(ids,N2,N1)=trcn_SnowDrift_flxM(ids,N2,N1)-trcn_SnowDrift_flxM_2DH(ids,N,NN,N5B,N4B)
      enddo
    ENDIF
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine TracerThruRofSnowXYM

!------------------------------------------------------------------------------------------

  subroutine TracerFluxThruSnowZM(I,J,M,N1,N2,snow_donor)
  ! Assemble recorded transfers on all fixed snow slots, including new snow.
  ! Top-slot deposition is accepted even when its initial heat capacity is zero.
  ! The proposal routine selects the discharge donor; do not reclassify layers
  ! here using a different capacity threshold or substep snapshot.
  implicit none
  integer, intent(in) :: I,J,M,N1,N2
  integer, intent(in) :: snow_donor ! zero means no discharge proposal
  character(len=*), parameter :: subname='TracerFluxThruSnowZM'
  integer :: LS,idg,idn

  call PrintInfo('beg '//subname)

  DO LS=1,JS
    DO idg=idg_beg,idg_NH3
      trcg_Aqua_flxM_snvr(idg,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg,LS,N2,N1) &
        +trcg_AquaAdv_flxM_snvr(idg,LS,N2,N1)
      if(LS<JS)then
        trcg_Aqua_flxM_snvr(idg,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg,LS,N2,N1) &
          -trcg_AquaAdv_flxM_snvr(idg,LS+1,N2,N1)
      endif
    ENDDO
    DO idn=ids_nut_beg,ids_nuts_end
      trcn_Aqua_flxM_snvr(idn,LS,N2,N1)=trcn_Aqua_flxM_snvr(idn,LS,N2,N1) &
        +trcn_AquaAdv_flxM_snvr(idn,LS,N2,N1)
      if(LS<JS)then
        trcn_Aqua_flxM_snvr(idn,LS,N2,N1)=trcn_Aqua_flxM_snvr(idn,LS,N2,N1) &
          -trcn_AquaAdv_flxM_snvr(idn,LS+1,N2,N1)
      endif
    ENDDO

    if(LS==snow_donor)then
      DO idg=idg_beg,idg_NH3
        trcg_Aqua_flxM_snvr(idg,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg,LS,N2,N1) &
          -trcg_AquaADV_Snow2Litr_flxM(idg,N2,N1)-trcg_AquaADV_Snow2Soil_flxM(idg,N2,N1)
      ENDDO
      trcg_Aqua_flxM_snvr(idg_NH3,LS,N2,N1)=trcg_Aqua_flxM_snvr(idg_NH3,LS,N2,N1) &
        -trcn_AquaADV_Snow2Band_flxM(idg_NH3B,N2,N1)
      DO idn=0,ids_nuts
        trcn_Aqua_flxM_snvr(ids_NH4+idn,LS,N2,N1)=trcn_Aqua_flxM_snvr(ids_NH4+idn,LS,N2,N1) &
          -trcn_AquaADV_Snow2Litr_flxM(ids_NH4+idn,N2,N1) &
          -trcn_AquaADV_Snow2Soil_flxM(ids_NH4+idn,N2,N1) &
          -trcn_AquaADV_Snow2Band_flxM(ids_NH4B+idn,N2,N1)
      ENDDO
    endif
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine TracerFluxThruSnowZM
! ----------------------------------------------------------------------

  subroutine XBoundaryTracerTranspXYZM(M,N,NN,M1,M2,M3,M4,M5,M6)
  !
  !Description:
  !tracer advection across the boundaries in 3 different directions
  !
  implicit none
  integer, intent(in) :: M,N,NN
  integer, intent(in) :: M1,M2,M3
  integer, intent(in) :: M4,M5,M6    !note when NN==iBehind, (M3,M2,M1)==(M6,M5,M4)

  character(len=*), parameter :: subname='XBoundaryTracerTranspXYZM'
  integer :: K,idg,idn,ids,idom
  real(r8) :: VFLW

  call PrintInfo('beg '//subname)

  !------------------
  !Micropore
  !out of cell  (M3,M2,M1), horizontal flow + bottom, exclude top
  IF(NN.EQ.iFront .AND. WaterFlow2MicPM_3D(M,N,M6,M5,M4).GT.0.0_r8 &  
    .OR. (NN.EQ.iBehind .AND. WaterFlow2MicPM_3D(M,N,M6,M5,M4).LT.0.0_r8))THEN
    IF(VLWatMicPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MicPM_3D(M,N,M6,M5,M4)/VLWatMicPM_vr(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF

    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_MicpTranspFlxM_3D(idom,K,N,M6,M5,M4)=VFLW*AZMAX1(DOM_MicP2_vr(idom,K,M3,M2,M1))
      enddo
    enddo

    DO idn=ids_beg,ids_end
      trcs_MicpTranspFlxM_3D(idn,N,M6,M5,M4)=VFLW*AZMAX1(trcs_solml2_vr(idn,M3,M2,M1))*trcs_VLN_vr(idn,M3,M2,M1)

      if(abs(trcs_MicpTranspFlxM_3D(idn,N,M6,M5,M4))>1.e20)then
        write(*,*)VFLW,trcs_solml2_vr(idn,M3,M2,M1),trcs_VLN_vr(idn,M3,M2,M1)
        call endrun(trim(mod_filename)//' at line',__LINE__)                
      endif
    ENDDO
    !
    !  SOLUTE GAIN WITH SUBSURFACE MICROPORE WATER GAIN
    !coming into cell (N3,N2,N1) from the boundaries
  ELSE
    DO  K=1,jcplx
      DOM_MicpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
    enddo
    trcs_MicpTranspFlxM_3D(idg_beg:idg_NH3-1,N,M6,M5,M4)=0.0_r8

    !Assuming the nutrient tracer concentations follow that from irrigation
    !This seems not that reasonable when irrigation is not included as part of WaterFlow2MicPM_3D
    !double check  
    DO ids=ids_nuts_beg,ids_nuts_end
      trcs_MicpTranspFlxM_3D(ids,N,M6,M5,M4)=WaterFlow2MicPM_3D(M,N,M6,M5,M4)*trcn_irrig_vr(ids,M3,M2,M1)*trcs_VLN_vr(ids,M3,M2,M1)
      if(abs(trcs_MicpTranspFlxM_3D(ids,N,M6,M5,M4))>1.e20)then
        write(*,*)WaterFlow2MicPM_3D(M,N,M6,M5,M4), &
          trcn_irrig_vr(ids,M3,M2,M1),trcs_VLN_vr(ids,M3,M2,M1)
        call endrun(trim(mod_filename)//' at line',__LINE__)                
      endif        
    ENDDO
  ENDIF

  !------------------------------------------------------
  !     SOLUTE LOSS WITH SUBSURFACE MACROPORE WATER LOSS
  !
  IF(NN.EQ.iFront .AND. WaterFlow2MacPM_3D(M,N,M6,M5,M4).GT.0.0_r8 &
    .OR. (NN.EQ.iBehind .AND. WaterFlow2MacPM_3D(M,N,M6,M5,M4).LT.0.0_r8))THEN
    IF(VLWatMacPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=AMAX1(-VFLWX,AMIN1(VFLWX,WaterFlow2MacPM_3D(M,N,M6,M5,M4)/VLWatMacPM_vr(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF

    DO  K=1,jcplx
      do idom=idom_beg,idom_end
        DOM_MacpTranspFlxM_3D(idom,K,N,M6,M5,M4)=VFLW*AZMAX1(DOM_MacP2_vr(idom,K,M3,M2,M1))
      enddo
    enddo

    DO idn=ids_beg,ids_end
      trcs_MacpTranspFlxM_3D(idn,N,M6,M5,M4)=VFLW*AZMAX1(trcs_soHml2_vr(idn,M3,M2,M1))*trcs_VLN_vr(idn,M3,M2,M1)
    ENDDO
    !
    !     NO SOLUTE GAIN IN SUBSURFACE MACROPORES
    !
  ELSE
    DO  K=1,jcplx
      DOM_MacpTranspFlxM_3D(idom_beg:idom_end,K,N,M6,M5,M4)=0.0_r8
    enddo
    trcs_MacpTranspFlxM_3D(ids_beg:ids_end,N,M6,M5,M4)=0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine XBoundaryTracerTranspXYZM
!------------------------------------------------------------------------------------------

  subroutine XSurfBoundaryTracerRunoffXYM(M,N,NN,N1,N2,M4,M5,RCHQF)

  !Description
  !
  implicit none
  integer, intent(in) :: M,N,NN,N1,N2,M4,M5
  real(r8), intent(in) :: RCHQF
  character(len=*), parameter :: subname='XSurfBoundaryTracerRunoffXYM'
  real(r8) :: FQRM
  integer :: K,idg,idn,ids,idom
  logical :: LNOFlow

  call PrintInfo('beg '//subname)
  LNOFlow=.not.XGridRunoffFlag_2DH(NN,N,N2,N1) .OR. isclose(RCHQF,0.0_r8) .OR. SurfRunoffPotentM_col(M,N2,N1).LE.ZEROS(N2,N1)

  IF(.not. LNOFlow)THEN
    !
    ! SOLUTE LOSS FROM RUNOFF DEPENDING ON ASPECT
    ! AND BOUNDARY CONDITIONS SET IN SITE FILE
    ! use fluxes from subroutine XBoundSurfaceRunoffM in SurfPhysMod.F90
    ! Lose tracer to external environment
    IF((NN.EQ.iFront .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &          
      .OR. (NN.EQ.iBehind .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN  

      FQRM=QflxSurfRunoffM_2DH(M,N,NN,M5,M4)/SurfRunoffPotentM_col(M,N2,N1)
      DO  K=1,jcplx
        DO idom=idom_beg,idom_end
          DOM_FloXSurRof_flxM_2DH(idom,K,N,NN,M5,M4)=DOM_FloXSurRunoff_PotFlxM(idom,K,N2,N1)*FQRM
        ENDDO
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_SurRof_flxM_2DH(idg,N,NN,M5,M4)=trcg_FloXSurRunoff_PotFlxM(idg,N2,N1)*FQRM
      ENDDO

      DO ids=ids_nut_beg,ids_nuts_end
        trcn_SurRof_flxM_2DH(ids,N,NN,M5,M4)=trcn_FloXSurRunoff_PotFlxM(ids,N2,N1)*FQRM
      ENDDO
      !
      !     SOLUTE GAIN FROM RUNON DEPENDING ON ASPECT
      !     AND BOUNDARY CONDITIONS SET IN SITE FILE
      !
    ELSEIF((NN.EQ.iBehind .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).GT.ZEROS(N2,N1)) &
      .OR.(NN.EQ.iFront .AND. QflxSurfRunoffM_2DH(M,N,NN,M5,M4).LT.ZEROS(N2,N1)))THEN

      !For incoming flux, assume the DOM concentation is zero 
      DO  K=1,jcplx
        DOM_FloXSurRof_flxM_2DH(idom_beg:idom_end,K,N,NN,M5,M4)=0.0_r8
      enddo

      DO idg=idg_beg,idg_NH3
        trcg_SurRof_flxM_2DH(idg,N,NN,M5,M4) = QflxSurfRunoffM_2DH(M,N,NN,M5,M4)*trcg_rain_mole_conc_col(idg,M5,M4)*MolecularWeight(idg)
      ENDDO

      trcn_SurRof_flxM_2DH(ids_nut_beg:ids_nuts_end,N,NN,M5,M4)=0.0_r8
    ENDIF
  ENDIF
  !
  call PrintInfo('end '//subname)
  end subroutine XSurfBoundaryTracerRunoffXYM
!------------------------------------------------------------------------------------------

  subroutine TracerXBoundariesM(I,J,M,NHW,NHE,NVN,NVS,snow_donor_col)
  !
  !Description:
  !Do transport across the boundaries, upper surface, subsurface lateral,
  !and bottom
  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NHW, NHE, NVN, NVS

  integer, intent(in) :: snow_donor_col(JY,JX) ! discharge donors selected during proposal formation
  character(len=*), parameter :: subname='TracerXBoundariesM'
  integer :: NY,NX,L
  integer :: NN,N
  integer :: N1,N2,N3,N4,N5,N6,N4B,N5B  !inner grids
  integer :: M1,M2,M3,M4,M5,M6          !boundary exchange
  real(r8) :: RCHQF

  call PrintInfo('beg '//subname)
!
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      D9585: DO L=NU_col(NY,NX),NL_col(NY,NX)
        N1=NX;N2=NY;N3=L
        M1 = NX;M2 = NY;M3 = L             !source grid                            
        D9580: DO  N=FlowDirIndicator_col(NY,NX),3
          D9575: DO  NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              N4  = NX+1; N5  = NY          
              N4B = NX-1;N5B = NY;N6  = L
              IF(NN.EQ.iFront)THEN                   !eastward                
                IF(NX.EQ.NHE)THEN                    !eastern boundary
                  M4 = NX+1;M5 = NY;M6 = L           !dest
                  RCHQF  = RechargEastSurf_col(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN             !west                
                IF(NX.EQ.NHW)THEN                   !western boundary
                  M4 = NX;M5 = NY;M6 = L            !dest
                  RCHQF  = RechargWestSurf_col(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN              
              N4  = NX;N5  = NY+1  
              N4B = NX;N5B = NY-1; N6  = L
              IF(NN.EQ.iFront)THEN  
                IF(NY.EQ.NVS)THEN                     ! southern boundary
                  M4 = NX;M5 = NY+1;M6 = L            !target grid
                  RCHQF  = RechargSouthSurf_col(M2,M1)
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN              !north                
                IF(NY.EQ.NVN)THEN                    !northern boundary
                  M4 = NX;M5 = NY;M6 = L             !target
                  RCHQF  = RechargNorthSurf_col(M5,M4)
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN      !vertical
              N4 = NX;N5 = NY;N6 = L+1               !target
              IF(NN.EQ.iFront)THEN                
                IF(L.EQ.NL_col(NY,NX))THEN               !lower boundary
                  M4 = NX;M5 = NY;M6 = L+1           !target grid
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN              
                cycle
              ENDIF
            ENDIF            
            !
            call XBoundaryTracerFlowM(L,N,NN,M,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
            !
          ENDDO D9575
          !
          IF(L.EQ.NUM_col(N2,N1))THEN
            IF(N.NE.iVerticalDirection)THEN
              !horizontal flow
              call TracerThruRofSnowXYM(M,N,N1,N2,N4,N5,N4B,N5B)
              !
            ELSEIF(N.EQ.iVerticalDirection)THEN
              ! vertical direction
              call TracerFluxThruSnowZM(I,J,M,N1,N2,snow_donor_col(N2,N1))
            ENDIF
          ENDIF

        ENDDO D9580
      ENDDO D9585

    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine TracerXBoundariesM

! ----------------------------------------------------------------------
  subroutine XBoundaryTracerFlowM(L,N,NN,M,N1,N2,N3,M1,M2,M3,M4,M5,M6,RCHQF)
  implicit none
  !
  !Do cross boundary tracer transport, lateral+bottom, exclude top 
  !
  integer, intent(in) :: L,N, NN, M 
  integer, intent(in) :: N1, N2,N3
  integer, intent(in) :: M1, M2,M3,M4,M5,M6
  real(r8), intent(in):: RCHQF
  character(len=*), parameter :: subname='XBoundaryTracerFlowM'
  real(r8) :: FLGM,FQRM,VFLW
  integer :: K,idg

! begin_execution
  call PrintInfo('beg '//subname)

  !Surface/litter horizontal tracer flow
  IF(L.EQ.NUM_col(M2,M1) .AND. N.NE.iVerticalDirection)THEN        
    call XSurfBoundaryTracerRunoffXYM(M,N,NN,N1,N2,M4,M5,RCHQF)
  ENDIF
!
  IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    IF(FlowDirIndicator_col(M2,M1).NE.3 .OR. N.EQ.iVerticalDirection)THEN        
      call XBoundaryTracerTranspXYZM(M,N,NN,M1,M2,M3,M4,M5,M6)
    ENDIF
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine XBoundaryTracerFlowM

!------------------------------------------------------------------------------------------

  subroutine NetTracerFlowXSoilPoresM(I,J,M,N,N1,N2,N3,N4,N5,N6)
  !
  !Description
  !Vertical or lateral flux between grid grids
  !does not include litter layer
  implicit none

  integer, intent(in) :: I,J,N,M
  integer, intent(in) :: N1,N2,N3
  integer, intent(in) :: N4,N5,N6
  character(len=*), parameter :: subname='NetTracerFlowXSoilPoresM'
  integer :: K,LL,ids,idg,idom
  integer :: first_face,last_face

  call SoilFaceIncomingRange(N,N3,N2,N1,first_face,last_face)
  call PrintInfo('beg '//subname)
  
  DO  K=1,jcplx
    do idom=idom_beg,idom_end
      DOM_Transp2Micp_flxM_vr(idom,K,N3,N2,N1)=DOM_Transp2Micp_flxM_vr(idom,K,N3,N2,N1) &
        +sum(DOM_MicpTranspFlxM_3D(idom,K,N,first_face:last_face,N2,N1))-DOM_MicpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4)

      DOM_Transp2Macp_flxM_vr(idom,K,N3,N2,N1)=DOM_Transp2Macp_flxM_vr(idom,K,N3,N2,N1) &
        +sum(DOM_MacpTranspFlxM_3D(idom,K,N,first_face:last_face,N2,N1))-DOM_MacpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4)

      if(N.NE.iVerticalDirection)then    
        DOM_transpFlxM_2DH(idom,K,N2,N1)=DOM_transpFlxM_2DH(idom,K,N2,N1) &
          +sum(DOM_MicpTranspFlxM_3D(idom,K,N,first_face:last_face,N2,N1))-DOM_MicpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4) &
          +sum(DOM_MacpTranspFlxM_3D(idom,K,N,first_face:last_face,N2,N1))-DOM_MacpTranspFlxM_3D(idom,K,N,SoilFaceSlot(N,N3,N6),N5,N4)
      endif    
    enddo
  enddo
  !
  !incoming flux (N3,N2,N1), outgoing flux (N6,N5,N4), N is direction
  !
  DO ids=ids_beg,ids_end
    trcs_Transp2Micp_flxM_vr(ids,N3,N2,N1)=trcs_Transp2Micp_flxM_vr(ids,N3,N2,N1) &
      +sum(trcs_MicpTranspFlxM_3D(ids,N,first_face:last_face,N2,N1))-trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)

    trcs_Transp2Macp_flxM_vr(ids,N3,N2,N1)=trcs_Transp2Macp_flxM_vr(ids,N3,N2,N1) &
      +sum(trcs_MacpTranspFlxM_3D(ids,N,first_face:last_face,N2,N1))-trcs_MacpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)

    if(N.NE.iVerticalDirection)then
       trcs_transpFlxM_2DH(ids,N2,N1)=trcs_transpFlxM_2DH(ids,N2,N1) &
         +sum(trcs_MicpTranspFlxM_3D(ids,N,first_face:last_face,N2,N1))-trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4) &
         +sum(trcs_MacpTranspFlxM_3D(ids,N,first_face:last_face,N2,N1))-trcs_MacpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)
    endif

    if(abs(trcs_Transp2Micp_flxM_vr(ids,N3,N2,N1))>1.e20)then
      write(*,*)N,N3,N2,N1,N6,N5,N4,trcs_names(ids)
      write(*,*)sum(trcs_MicpTranspFlxM_3D(ids,N,first_face:last_face,N2,N1)),-trcs_MicpTranspFlxM_3D(ids,N,SoilFaceSlot(N,N3,N6),N5,N4)
      call endrun(trim(mod_filename)//' at line',__LINE__)                
    endif  
  ENDDO

  call PrintInfo('end '//subname)

  end subroutine NetTracerFlowXSoilPoresM

!------------------------------------------------------------------------------------------

  subroutine UpdateSurfTracerM(I,J,M,NHW,NHE,NVN,NVS,dpscal,dpscal_dom,pscal1,pscal1_dom)
  !
  !Description
  !
  implicit none
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS
  real(r8),intent(in) :: dpscal(ids_beg:ids_end)
  real(r8),intent(in) :: dpscal_dom(idom_beg:idom_end)
  real(r8), optional, intent(inout) :: pscal1(ids_beg:ids_end)
  real(r8), optional, intent(inout) :: pscal1_dom(idom_beg:idom_end)

  character(len=*), parameter :: subname='UpdateSurfTracerM'

  real(r8) :: pscal(ids_beg:ids_end)
  real(r8) :: pscal_dom(idom_beg:idom_end)
  integer :: K,idg,ids,idom,L,NY,NX
  real(r8) :: flux,flux0
  logical :: lflux,lfupdate

  call PrintInfo('beg '//subname)
  lfupdate=.false.
  if(present(pscal1))then
    pscal=pscal1
  else
    pscal=1.000_r8
    lfupdate=.true.
  endif  

  if(present(pscal1_dom))then
    pscal_dom=pscal1_dom
  else
    pscal_dom=1.000_r8
  endif  

  DO NX=NHW,NHE
    DO  NY=NVN,NVS  
      !-----------------------------------------------------------------------------------------------
      !to litter
      DO  K=1,micpar%NumOfLitrCmplxs
        DO idom=idom_beg,idom_end
          
          flux0=DOM_SurfRunoff_flxM(idom,K,NY,NX)+DOM_MicpTranspFlxM_3D(idom,K,3,0,NY,NX)          
          lflux=.true. .or. .not.isclose(flux0,0._r8)
          if(lflux .and. dpscal_dom(idom)>tiny_p)then  
            flux=flux0*dpscal_dom(idom)
            call get_flux_scalar(DOM_MicP2_vr(idom,K,0,NY,NX),flux,DOM_MicP_vr(idom,K,0,NY,NX),pscal_dom(idom))

            if(lfupdate)then
              DOM_mass3_col(idom,K,NY,NX)         = DOM_mass3_col(idom,K,NY,NX)+DOM_MicP2_vr(idom,K,0,NY,NX)
              DOM_MicP2_vr(idom,K,0,NY,NX)        = DOM_MicP_vr(idom,K,0,NY,NX)
              DOM_mass4_col(idom,K,NY,NX)         = DOM_mass4_col(idom,K,NY,NX)+DOM_MicP2_vr(idom,K,0,NY,NX)

              TranspNetDOM_flxM_col(idom,K,NY,NX) = TranspNetDOM_flxM_col(idom,K,NY,NX)+flux0                                       
            endif  
          endif
        ENDDO
      ENDDO

     !exclude NH3B
      DO idg=idg_beg,idg_NH3-1        
        flux = RGasAtmDisol2LitrM_col(idg,NY,NX)+trcg_SurfRunoff_flxM(idg,NY,NX)+trcs_MicpTranspFlxM_3D(idg,3,0,NY,NX)        
        lflux=.true. .or. .not.isclose(flux,0._r8)
        if(lflux .and. dpscal(idg)>tiny_p)then          
          flux=flux*dpscal(idg)
          call get_flux_scalar(trcs_solml2_vr(idg,0,NY,NX),flux,trcs_solml_vr(idg,0,NY,NX),pscal(idg))
          if(pscal(idg)<0._r8)then
            write(193,*)(I*1000+J)*100+M,trcs_names(idg),trcs_solml2_vr(idg,0,NY,NX),flux,trcs_solml_vr(idg,0,NY,NX),pscal(idg)
            call endrun(trim(mod_filename)//' at line',__LINE__)
          endif
          if(lfupdate)then
            trcs_solml2_vr(idg,0,NY,NX)=trcs_solml_vr(idg,0,NY,NX)
          endif  
        endif
      ENDDO

      idg=idg_NH3      
      flux = RGasAtmDisol2LitrM_col(idg,NY,NX)+trcg_SurfRunoff_flxM(idg,NY,NX)+trcs_MicpTranspFlxM_3D(idg,3,0,NY,NX)         
      lflux=.true. .or. .not.isclose(flux,0._r8)  
      if(lflux .and. dpscal(idg)>tiny_p)then          
        flux=flux*dpscal(idg)
        call get_flux_scalar(trcs_solml2_vr(idg,0,NY,NX),flux,trcs_solml_vr(idg,0,NY,NX),pscal(idg))
        if(lfupdate)then
          trcs_solml2_vr(idg,0,NY,NX)=trcs_solml_vr(idg,0,NY,NX)
        endif  
      endif

      DO ids=ids_nut_beg,ids_nuts_end        
        flux=trcs_MicpTranspFlxM_3D(ids,3,0,NY,NX)+trcn_SurfRunoff_flxM(ids,NY,NX)
        lflux=.true. .or. .not.isclose(flux,0._r8)
        if(lflux .and. dpscal(ids)>tiny_p)then          
          flux=flux*dpscal(ids)
          call get_flux_scalar(trcs_solml2_vr(ids,0,NY,NX),flux,trcs_solml_vr(ids,0,NY,NX),pscal(ids))
          if(lfupdate)then
            trcs_solml2_vr(ids,0,NY,NX)=trcs_solml_vr(ids,0,NY,NX)
          endif  
        endif
      ENDDO
        
      !----------------------------------------------------------------------------------------------------
      !topsoil
      L=NU_col(NY,NX)
      IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN            
        DO idg=idg_beg,idg_end
          flux0=RGasAtmDisol2SoilM_col(idg,NY,NX)+trcsol_Irrig_flxM_vr(idg,L,NY,NX) & 
            +trcs_Transp2Micp_flxM_vr(idg,L,NY,NX) + trcs_Mac2MicPore_flxM_vr(idg,L,NY,NX)
               
          lflux=.true. .or. .not.isclose(flux0,0._r8)              
          if(lflux .and. dpscal(idg)>tiny_p)then          
            flux=flux0*dpscal(idg)    
            call get_flux_scalar(trcs_solml2_vr(idg,L,NY,NX),flux,trcs_solml_vr(idg,L,NY,NX),pscal(idg))

            if(trcs_solml_vr(idg,L,NY,NX)<0._r8)then
              write(132,*)(I*1000+J)*100+M,trcs_names(idg),L,trcs_solml2_vr(idg,L,NY,NX),flux,trcs_solml_vr(idg,L,NY,NX),pscal(idg),dpscal(idg)
              call endrun(trim(mod_filename)//' at line',__LINE__)
            endif  
            if(lfupdate)then
              trcs_solml2_vr(idg,L,NY,NX)            = trcs_solml_vr(idg,L,NY,NX)
              TranspNetSoil_slow_flxM_col(idg,NY,NX) = TranspNetSoil_slow_flxM_col(idg,NY,NX)+ &
                RGasAtmDisol2SoilM_col(idg,NY,NX)+trcsol_Irrig_flxM_vr(idg,L,NY,NX) & 
                +trcs_Transp2Micp_flxM_vr(idg,L,NY,NX) + trcs_Transp2Macp_flxM_vr(idg,L,NY,NX)
            endif  
          endif
        ENDDO

        DO ids=ids_NH4B,ids_nuts_end          
          flux = trcsol_Irrig_flxM_vr(ids,L,NY,NX)+trcs_Transp2Micp_flxM_vr(ids,L,NY,NX) + trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)
            
          lflux=.true. .or. .not.isclose(flux,0._r8)  
          if(lflux .and. dpscal(ids)>tiny_p)then          
            flux=flux*dpscal(ids)     
            call get_flux_scalar(trcs_solml2_vr(ids,L,NY,NX),flux,trcs_solml_vr(ids,L,NY,NX),pscal(ids))
            if(lfupdate)then
              trcs_solml2_vr(ids,L,NY,NX)=trcs_solml_vr(ids,L,NY,NX)
            endif  
          endif
        ENDDO

        DO ids=ids_beg,ids_end          
          flux = trcs_Transp2Macp_flxM_vr(ids,L,NY,NX) - trcs_Mac2MicPore_flxM_vr(ids,L,NY,NX)
          lflux=.true. .or. .not.isclose(flux,0._r8)
          if(lflux .and. dpscal(ids)>tiny_p)then          
            flux=flux*dpscal(ids)     
            call get_flux_scalar(trcs_soHml2_vr(ids,L,NY,NX), flux, trcs_soHml_vr(ids,L,NY,NX), pscal(ids))
            if(lfupdate)then
              trcs_soHml2_vr(ids,L,NY,NX)=trcs_soHml_vr(ids,L,NY,NX)
            endif
          endif
        ENDDO

        DO  K=1,jcplx
          DO idom=idom_beg,idom_end            
            flux=DOM_Transp2Micp_flxM_vr(idom,K,L,NY,NX)+DOM_Mac2MicPore_flxM_vr(idom,K,L,NY,NX)
            lflux=.true. .or. .not.isclose(flux,0._r8)
            if(lflux .and. dpscal_dom(idom)>tiny_p)then          
              flux=flux*dpscal_dom(idom)
              call get_flux_scalar(DOM_MicP2_vr(idom,K,L,NY,NX),flux, DOM_MicP_vr(idom,K,L,NY,NX),pscal_dom(idom))

              flux=DOM_Transp2Macp_flxM_vr(idom,K,L,NY,NX)-DOM_Mac2MicPore_flxM_vr(idom,K,L,NY,NX)
              flux=flux*dpscal_dom(idom)
              call get_flux_scalar(DOM_MacP2_vr(idom,K,L,NY,NX), flux, DOM_MacP_vr(idom,K,L,NY,NX),pscal_dom(idom))

              if(lfupdate)then
                DOM_mass3_col(idom,K,NY,NX)=DOM_mass3_col(idom,K,NY,NX)+DOM_MicP2_vr(idom,K,L,NY,NX)+DOM_MacP2_vr(idom,K,L,NY,NX)
                DOM_MicP2_vr(idom,K,L,NY,NX)=DOM_MicP_vr(idom,K,L,NY,NX)
                DOM_MacP2_vr(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)
                DOM_mass4_col(idom,K,NY,NX)=DOM_mass4_col(idom,K,NY,NX)+DOM_MicP2_vr(idom,K,L,NY,NX)+DOM_MacP2_vr(idom,K,L,NY,NX)                
                TranspNetDOM_flxM_col(idom,K,NY,NX)=TranspNetDOM_flxM_col(idom,K,NY,NX)+ &
                  DOM_Transp2Micp_flxM_vr(idom,K,L,NY,NX)+DOM_Transp2Macp_flxM_vr(idom,K,L,NY,NX)

              endif
            endif
          ENDDO
        ENDDO
      ENDIF  
    ENDDO
  ENDDO
  if(present(pscal1))pscal1=pscal
  if(present(pscal1_dom))pscal1_dom=pscal_dom
  call PrintInfo('end '//subname)
  end subroutine UpdateSurfTracerM

!------------------------------------------------------------------------------------------

  subroutine UpdateSnowTracersM(I,J,M,NHW,NHE,NVN,NVS,dpscal,pscal1)
  implicit none
  integer, intent(in) :: I,J,M
  integer, intent(in) :: NHW,NHE,NVN,NVS
  real(r8),intent(in) :: dpscal(ids_beg:ids_end)
  real(r8),optional,intent(inout):: pscal1(ids_beg:ids_end)
  integer :: L,idg,idn,NY,NX
  real(r8) :: b,flux
  real(r8) :: pscal(ids_beg:ids_end)
  logical :: lflux,lfupdate
!
!
  if(present(pscal1))then
    pscal=pscal1
    lfupdate=.false.
  else
    pscal=1.000_r8
    lfupdate=.true.
!    write(555,*)'snowup'
  endif  

  DO NX=NHW,NHE
    DO  NY=NVN,NVS  

      DO idg=idg_beg,idg_NH3
        flux=trcg_SnowDrift_flxM(idg,NY,NX)+trcg_Aqua_flxM_snvr(idg,1,NY,NX)
        lflux=.true. .or. .not.isclose(flux,0._r8)
        if(lflux .and. dpscal(idg)>tiny_p)then          
          flux=flux*dpscal(idg)     
          call get_flux_scalar(trcg_solsml2_snvr(idg,1,NY,NX),flux,trcg_solsml_snvr(idg,1,NY,NX),pscal(idg))
          if(lfupdate)then
            trcg_solsml2_snvr(idg,1,NY,NX)=trcg_solsml_snvr(idg,1,NY,NX)
          endif  
        endif
      ENDDO

      DO idn=ids_nut_beg,ids_nuts_end
        flux=trcn_SnowDrift_flxM(idn,NY,NX)+trcn_Aqua_flxM_snvr(idn,1,NY,NX)      
        lflux=.true. .or. .not.isclose(flux,0._r8)
        if(lflux .and. dpscal(idn)>tiny_p)then          
          flux=flux*dpscal(idn)   
          call get_flux_scalar(trcn_solsml2_snvr(idn,1,NY,NX),flux,trcn_solsml_snvr(idn,1,NY,NX),pscal(idn))
          if(lfupdate)trcn_solsml2_snvr(idn,1,NY,NX)=trcn_solsml_snvr(idn,1,NY,NX)
        endif
      ENDDO

      DO  L=2,JS
        DO idg=idg_beg,idg_NH3
          flux=trcg_Aqua_flxM_snvr(idg,L,NY,NX)        
          lflux=.true. .or. .not.isclose(flux,0._r8)
          if(lflux .and. dpscal(idg)>tiny_p)then          
            flux=flux*dpscal(idg)   
            call get_flux_scalar(trcg_solsml2_snvr(idg,L,NY,NX),flux, trcg_solsml_snvr(idg,L,NY,NX),pscal(idg))
            if(pscal(idg)<0._r8 .or. trcg_solsml_snvr(idg,L,NY,NX) < 0._r8)then              
              if(lfupdate)then
                write(193,*)(I*1000+J)*100+M,'final',L,trcs_names(idg),trcg_solsml2_snvr(idg,L,NY,NX),flux, trcg_solsml_snvr(idg,L,NY,NX),pscal(idg),dpscal(idg)
                call endrun(trim(mod_filename)//' at line',__LINE__)
              endif  
            endif
            if(lfupdate)then
              trcg_solsml2_snvr(idg,L,NY,NX)=trcg_solsml_snvr(idg,L,NY,NX)
            endif
          endif  
        ENDDO

        DO idn=ids_nut_beg,ids_nuts_end
          flux=trcn_Aqua_flxM_snvr(idn,L,NY,NX)        
          lflux=.true. .or. .not.isclose(flux,0._r8)
          if(lflux .and. dpscal(idn)>tiny_p)then          
            flux=flux*dpscal(idn)             
            call get_flux_scalar(trcn_solsml2_snvr(idn,L,NY,NX),flux,trcn_solsml_snvr(idn,L,NY,NX),pscal(idn))
            if(lfupdate)trcn_solsml2_snvr(idn,L,NY,NX)=trcn_solsml_snvr(idn,L,NY,NX)
          endif  
        ENDDO
      ENDDO

    ENDDO
  ENDDO
  if(present(pscal1))pscal1=pscal

  end subroutine UpdateSnowTracersM

end module TranspNoSaltSlowMod
