module TranspNoSaltFastMod

  ! Description:
  ! Advance fast gas transport, gas-water exchange, and tracer source/sink accounting.

  use SoilFaceMappingMod, only: SoilFaceSlot,SoilFaceIncomingRange
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: destroy, endrun
  use EcoSIMCtrlMod, only: iVerbLevel
  use NumericalAuxMod
  use minimathmod
  use DebugToolMod
  use GridDataType  
  use TranspNoSaltDataMod
  use SoilBGCDataType
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SurfLitterDataType
  use AqueChemDatatype
  use SnowDataType
  use EcoSimConst
  use TracerIDMod
  use ChemTranspDataType
  use ClimForcDataType
  use EcoSIMSolverPar
  use LandSurfDataType
  use SurfSoilDataType
  use PlantDataRateType  
  use SoilPropertyDataType
implicit none
  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: TransptFastNoSaltMM
  contains


  subroutine TransptFastNoSaltMM(I,J,M,MM,NHE,NHW,NVS,NVN)

  !Description:
  !update trcs_solml2_vr and trcg_gasml2_vr
  !accounting for gaseous diffusion and dissolution
  !it does not consider litter-topsoil gas exchange, rather, their
  !changes are against the atmosphere, that is the gas concentation in 
  !litter is assumed to be equal to that in the atmosphere
  implicit none
  integer, intent(in) :: I,J,M,NHE,NHW,NVS,NVN,MM
  character(len=*), parameter :: subname='TransptFastNoSaltMM'
  real(r8) :: pscal(idg_beg:idg_end), dpscal(idg_beg:idg_end)
  real(r8) :: dpscal_max
  integer :: iterm,idg

  call PrintInfo('beg '//subname)
  call EnterMassCheck(I,J,NHE,NHW,NVS,NVN)
  dpscal=1._r8
  dpscal_max=1._r8
  iterm=0
  DO while(dpscal_max>1.e-2_r8 .and. iterm<3)
    iterm=iterm+1
    pscal=1.000_r8

    call ZeroTracerFluxMM(I,J,NHE,NHW,NVS,NVN)

    call LitterGasVolatilDissolMM(I,J,M,NHE,NHW,NVS,NVN)

    call SurfSoilFluxGasDifAdvMM(M,NHE,NHW,NVS,NVN)

    call TracerFlowXGridsMM(I,J,M,MM,NHE,NHW,NVS,NVN)

    call XBoundaryFluxMM(I,J,M,NHW,NHE,NVN,NVS)

    call GatherTranspFluxMM(I,J,M,NHW,NHE,NVN,NVS)

    call FastUpdateStateVarsMM(I,J,M,NHW, NHE, NVN, NVS,dpscal,pscal)

    call AccumFastFluxesMM(I,J,M,NHE,NHW,NVS,NVN,dpscal,pscal,MM)

    dpscal_max=0._r8
    DO idg=idg_beg,idg_end
      if(pscal(idg)>1.e-6_r8)then
        dpscal(idg)=dpscal(idg)*(1._r8-pscal(idg))
      else
        dpscal(idg)=0._r8
      endif
      dpscal_max=AMAX1(dpscal_max,dpscal(idg))
    enddo
    
  ENDDO

  call ApplyFastSourceSink(I,J,M,NHE,NHW,NVS,NVN,MM)

  call ExitMassCheck(I,J,M,NHE,NHW,NVS,NVN,MM)  
  call PrintInfo('end '//subname)
  end subroutine TransptFastNoSaltMM
!------------------------------------------------------------------------------------------  
  subroutine ApplyFastSourceSink(I,J,M,NHE,NHW,NVS,NVN,MM)  
  ! Description:
  ! Apply staged gas sources and sinks and accumulate their column mass budgets.
  implicit none
  integer, intent(in) :: I,J,NHE,NHW,NVS,NVN,M,MM

  integer :: NY,NX,L,idg
  real(r8) :: flux

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      DO idg=idg_beg,idg_NH3-1
        trcs_solml2_vr(idg,0,NY,NX)=trcs_solml2_vr(idg,0,NY,NX)+RBGCSrceGasMM_vr(idg,0,NY,NX)
        call SubstrateDribbling(RBGCSinkGasMM_vr(idg,0,NY,NX),trcs_solml_drib_vr(idg,0,NY,NX),trcs_solml2_vr(idg,0,NY,NX))
        flux=RBGCSrceGasMM_vr(idg,0,NY,NX)-RBGCSinkGasMM_vr(idg,0,NY,NX)
        trcg_NetPro_fast_col(idg,NY,NX)      = trcg_NetPro_fast_col(idg,NY,NX)+flux
        trcs_netProd_lit_fast_col(idg,NY,NX) = trcs_netProd_lit_fast_col(idg,NY,NX)+flux
      ENDDO

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO idg=idg_beg,idg_NH3-1
          trcs_solml2_vr(idg,L,NY,NX)=trcs_solml2_vr(idg,L,NY,NX)+RBGCSrceGasMM_vr(idg,L,NY,NX)

          call SubstrateDribbling(RBGCSinkGasMM_vr(idg,L,NY,NX),trcs_solml_drib_vr(idg,L,NY,NX),trcs_solml2_vr(idg,L,NY,NX))
          trcg_NetPro_fast_col(idg,NY,NX)    = trcg_NetPro_fast_col(idg,NY,NX)+RBGCSrceGasMM_vr(idg,L,NY,NX)-RBGCSinkGasMM_vr(idg,L,NY,NX)
          RGasNetProdSoil_col(idg,NY,NX) = RGasNetProdSoil_col(idg,NY,NX)+RBGCSrceGasMM_vr(idg,L,NY,NX)-RBGCSinkGasMM_vr(idg,L,NY,NX)

!          if(trcs_solml_drib_vr(idg,L,NY,NX)>1.e0_r8)then
!            write(*,*)(I*1000+J)*100+M,'idg=',trcs_names(idg),L,RBGCSinkGasMM_vr(idg,L,NY,NX),trcs_solml_drib_vr(idg,L,NY,NX),trcs_solml2_vr(idg,L,NY,NX)
!          endif
        ENDDO
        idg=idg_NH3
        ! NH3 fertilizer release enters the gas pool; other gas sources enter water.
        trcg_gasml2_vr(idg,L,NY,NX)=trcg_gasml2_vr(idg,L,NY,NX)+RBGCSrceGasMM_vr(idg,L,NY,NX)
        call SubstrateDribbling(RBGCSinkGasMM_vr(idg,L,NY,NX),trcs_solml_drib_vr(idg,L,NY,NX),trcg_gasml2_vr(idg,L,NY,NX))
        flux=RBGCSrceGasMM_vr(idg,L,NY,NX)-RBGCSinkGasMM_vr(idg,L,NY,NX)
        trcg_NetPro_fast_col(idg,NY,NX)=trcg_NetPro_fast_col(idg,NY,NX)+flux
        RGasNetProdSoil_col(idg,NY,NX)=RGasNetProdSoil_col(idg,NY,NX)+flux
      ENDDO  
    ENDDO
  ENDDO  

  end subroutine ApplyFastSourceSink
!------------------------------------------------------------------------------------------

  subroutine EnterMassCheck(I,J,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J,NHE,NHW,NVS,NVN

  integer :: NY,NX,idg,L

  DO NX=NHW,NHE
    DO  NY=NVN,NVS      
      trcg_NetPro_fast_col(:,NY,NX)       = 0._r8
      TranspNetSoil_fast_flx_col(:,NY,NX) = 0._r8
      RGas_Disol_FlxMM_vr(:,:,NY,NX)      = 0._r8
      trcg_mass3_fast_col(:,NY,NX)        = 0._r8
      GasDiff2Surf_fast_flx_col(:,NY,NX)  = 0._r8
      trcs_hydrloss_fast_flx_col(:,NY,NX) = 0._r8
      trcs_solml_dribM_beg_col(:,NY,NX)   = 0._r8
      trcs_drib_fast_beg_col(:,NY,NX)   = 0._r8
      trcs_drainage_fast_flx_col(:,NY,NX) = 0._r8
      trcs_netProd_lit_fast_col(:,NY,NX)  = 0._r8
      AtmGasDiff2Litr_fast_flx_col(:,NY,NX)=0._r8
      DO idg=idg_beg,idg_NH3
        trcg_mass_begf(idg,NY,NX)          = trcs_solml2_vr(idg,0,NY,NX)
        trcg_mass_litr_begf(idg,NY,NX)     = trcs_solml2_vr(idg,0,NY,NX)
        trcs_drib_fast_beg_col(idg,NY,NX)  = trcs_solml_drib_vr(idg,0,NY,NX)
        trcs_solml_dribM_beg_col(idg,NY,NX) = trcs_solml_drib_vr(idg,0,NY,NX)
        DO L=NU_col(NY,NX),NL_col(NY,NX)
          trcg_mass_begf(idg,NY,NX)          = trcg_mass_begf(idg,NY,NX)+trcg_gasml2_vr(idg,L,NY,NX)+trcs_solml2_vr(idg,L,NY,NX)
          trcs_solml_dribM_beg_col(idg,NY,NX) = trcs_solml_dribM_beg_col(idg,NY,NX)+trcs_solml_drib_vr(idg,L,NY,NX)        
        ENDDO 
      ENDDO

      idg=idg_NH3
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        trcg_mass_begf(idg,NY,NX)=trcg_mass_begf(idg,NY,NX)+trcs_solml2_vr(idg_NH3B,L,NY,NX)
      ENDDO       
    ENDDO  
  ENDDO
  end subroutine EnterMassCheck

!------------------------------------------------------------------------------------------

  subroutine ExitMassCheck(I,J,M,NHE,NHW,NVS,NVN,iterm)
  !
  !Description:
  ! Do mass conservation check for fast transport
  ! Fluxes include:
  ! net flux production: trcg_NetPro_fast_col(idg,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NHE,NHW,NVS,NVN,M,iterm

  integer :: NY,NX,idg,L
  real(r8) :: trcg_mass_now(idg_beg:idg_NH3)
  real(r8) :: trcg_mass_litr_mass_now(idg_beg:idg_NH3)
  real(r8) :: trcs_solml_drib_col(idg_beg:idg_NH3)
  real(r8) :: dmass,err

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      trcg_mass_now=0._r8
      DO idg=idg_beg,idg_NH3
        trcg_mass_now(idg)       = trcs_solml2_vr(idg,0,NY,NX)
        trcg_mass_litr_mass_now(idg)= trcs_solml2_vr(idg,0,NY,NX)
        trcs_solml_drib_col(idg) = trcs_solml_drib_vr(idg,0,NY,NX)
        DO L=NU_col(NY,NX),NL_col(NY,NX)
          trcg_mass_now(idg)       = trcg_mass_now(idg)+trcg_gasml2_vr(idg,L,NY,NX)+trcs_solml2_vr(idg,L,NY,NX)
          trcs_solml_drib_col(idg) = trcs_solml_drib_col(idg)+trcs_solml_drib_vr(idg,L,NY,NX)
        ENDDO 
      ENDDO

      idg=idg_NH3
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        trcg_mass_now(idg)=trcg_mass_now(idg)+trcs_solml2_vr(idg_NH3B,L,NY,NX)
      ENDDO

      DO idg=idg_beg,idg_NH3
        GasDiff2Surf_flx_col(idg,NY,NX)    = GasDiff2Surf_flx_col(idg,NY,NX)+GasDiff2Surf_fast_flx_col(idg,NY,NX)
        GasHydroLoss_flx_col(idg,NY,NX)    = GasHydroLoss_flx_col(idg,NY,NX)+trcs_hydrloss_fast_flx_col(idg,NY,NX)
        RGasNetProd_col(idg,NY,NX)         = RGasNetProd_col(idg,NY,NX)+trcg_NetPro_fast_col(idg,NY,NX)
        trcs_drainage_flx_col(idg,NY,NX)   = trcs_drainage_flx_col(idg,NY,NX)+trcs_drainage_fast_flx_col(idg,NY,NX)
        trcs_netProd_lit_col(idg,NY,NX)    = trcs_netProd_lit_col(idg,NY,NX)+trcs_netProd_lit_fast_col(idg,NY,NX)
        AtmGasDiff2Litr_flx_col(idg,NY,NX) = AtmGasDiff2Litr_flx_col(idg,NY,NX)+AtmGasDiff2Litr_fast_flx_col(idg,NY,NX)

        dmass = trcg_mass_now(idg)-trcg_mass_begf(idg,NY,NX)
        err   = dmass-trcg_NetPro_fast_col(idg,NY,NX)-GasDiff2Surf_fast_flx_col(idg,NY,NX)-trcs_hydrloss_fast_flx_col(idg,NY,NX) 
        err = err - trcs_solml_drib_col(idg)+trcs_solml_dribM_beg_col(idg,NY,NX)
        errmass_fast(idg,NY,NX)=errmass_fast(idg,NY,NX)+err
        if(abs(err)>1.e-5_r8 .OR. iVerbLevel==1 .or. trcs_solml_drib_col(idg)>1._r8)then
          if((iVerbLevel==1 .or. abs(err)>1.e-4_r8))then
            write(133,*)('-',L=1,50)
            write(133,*)(I*1000+J)*100+M,trcs_names(idg),'FAST'
            write(133,*)'init/final mass     =',trcg_mass_begf(idg,NY,NX),trcg_mass_now(idg),trcg_mass3_fast_col(idg,NY,NX)
            write(133,*)'dmass,              =',dmass,TranspNetSoil_fast_flx_col(idg,NY,NX)
            write(133,*)'dif                 =',GasDiff2Surf_fast_flx_col(idg,NY,NX)
            write(133,*)'hydrloss            =',trcs_hydrloss_fast_flx_col(idg,NY,NX) 
            write(133,*)'netpro              =',trcg_NetPro_fast_col(idg,NY,NX)
            write(133,*)'drib beg, end       =',trcs_solml_dribM_beg_col(idg,NY,NX),trcs_solml_drib_col(idg)
            write(133,*)'err                 =',err
            write(133,*)'drib vr.            =',trcs_solml_drib_vr(idg,0,NY,NX),trcs_solml_drib_vr(idg,NU_col(NY,NX):NL_col(NY,NX),NY,NX)
            write(133,*)'total sum'
            write(133,*)'GasDiff2Surf_flx_col=',GasDiff2Surf_flx_col(idg,NY,NX)
            write(133,*)'GasHydroLoss_flx    =',GasHydroLoss_flx_col(idg,NY,NX)            
            write(133,*)'RGasNetProd_col     =',RGasNetProd_col(idg,NY,NX)       
            write(133,*)'-----------------------------------'
            write(133,*) 'err litr           =',trcg_mass_litr_begf(idg,NY,NX)-trcg_mass_litr_mass_now(idg) &
              +trcs_netProd_lit_fast_col(idg,NY,NX)+ AtmGasDiff2Litr_fast_flx_col(idg,NY,NX) +trcs_solml_drib_vr(idg,0,NY,NX) &
              -trcs_drib_fast_beg_col(idg,NY,NX)                  
          endif
          if(abs(err)>1.e-4_r8)call endrun(trim(mod_filename)//' at line',__LINE__)          
        endif

      ENDDO

    ENDDO  
  ENDDO
  end subroutine ExitMassCheck
!------------------------------------------------------------------------------------------
  subroutine copyStateVars(NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: NHE,NHW,NVS,NVN
  integer :: NY,NX,idg,L

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      DO idg=idg_beg,idg_NH3      
       trcs_solml2_vr(idg,0,NY,NX) = AZERO(trcs_solml_vr(idg,0,NY,NX))
       DO L=NU_col(NY,NX),NL_col(NY,NX)
         trcg_gasml2_vr(idg,L,NY,NX)=AZERO(trcg_gasml_vr(idg,L,NY,NX))
       ENDDO
      enddo

      DO idg=idg_beg,idg_end 
        DO L=NU_col(NY,NX),NL_col(NY,NX)    
          trcs_solml2_vr(idg,L,NY,NX) = AZERO(trcs_solml_vr(idg,L,NY,NX))
        ENDDO
      enddo

    ENDDO
  ENDDO
  end subroutine copyStateVars

!------------------------------------------------------------------------------------------
  subroutine GatherTranspFluxMM(I,J,M,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I,J,M,NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='GatherTranspFluxMM'
  integer :: NY,NX,L,LL,N
  integer :: N6,N5,N4,N3,N2,N1,idg

  call PrintInfo('beg '//subname)
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      DO idg=idg_beg,idg_NH3
        RGasSinkScalar_vr(idg,0,NY,NX)=1._r8 !trcs_solml2_vr(idg,0,NY,NX)/(trcs_solml2_vr(idg,0,NY,NX) &
          !+trcs_solcoef_col(idg,NY,NX)*AMAX1(VLWatMicP_vr(0,NY,NX),1.e-5_r8))
      enddo

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO idg=idg_beg,idg_NH3
          RGasSinkScalar_vr(idg,L,NY,NX)=1._r8 !trcs_solml2_vr(idg,L,NY,NX)/(trcs_solml2_vr(idg,L,NY,NX) &
            !+trcs_solcoef_col(idg,NY,NX)*AMAX1(VLWatMicP_vr(L,NY,NX),1.E-5_r8))
        enddo

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
            call NetTracerFlowXSoilPoresMM(N,M,N1,N2,N3,N4,N5,N6)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO        
  call PrintInfo('end '//subname)
  end subroutine GatherTranspFluxMM

!------------------------------------------------------------------------------------------
  subroutine AccumFastFluxesMM(I,J,M,NHE,NHW,NVS,NVN,dpscal,pscal,iterm)
  implicit none
  integer, intent(in) :: I,J,M,iterm
  integer, intent(in) :: NHE,NHW,NVS,NVN
  real(r8), intent(in) :: dpscal(idg_beg:idg_end)
  real(r8), intent(in):: pscal(idg_beg:idg_end)
  character(len=*), parameter :: subname='AccumFastFluxesMM'
  real(r8) :: ppscal(idg_beg:idg_end)
  real(r8) :: pscal_max,flux
  integer :: idg,L,NY,NX

  call PrintInfo('beg '//subname)

  pscal_max=0._r8   
  DO idg=idg_beg,idg_end
    ppscal(idg) = dpscal(idg)*pscal(idg)
    pscal_max   = AMAX1(pscal(idg),pscal_max)
  ENDDO

  if(pscal_max<1.00001_r8)then
    call FastUpdateStateVarsMM(I,J,M,NHW, NHE, NVN, NVS,ppscal)
  endif

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      DO idg=idg_beg,idg_NH3   
        if(ppscal(idg)>tiny_p)then     
        
          GasDiff2Surf_fast_flx_col(idg,NY,NX) = GasDiff2Surf_fast_flx_col(idg,NY,NX)+(RGas_Disol_FlxMM_vr(idg,0,NY,NX) &
            +Gas_AdvDif_FlxMM_3D(idg,3,NU_col(NY,NX),NY,NX))*ppscal(idg)

          !litter layer
          flux                               = RGas_Disol_FlxMM_vr(idg,0,NY,NX)*ppscal(idg)
          AtmGasDiff2Litr_fast_flx_col(idg,NY,NX) = AtmGasDiff2Litr_fast_flx_col(idg,NY,NX)+flux

          flux                               = Gas_AdvDif_FlxMM_3D(idg,3,NU_col(NY,NX),NY,NX)*ppscal(idg)
          AtmGasDiff2Soil_flx_col(idg,NY,NX) = AtmGasDiff2Soil_flx_col(idg,NY,NX)+flux

          flux                                 = Gas_AdvDif_FlxMM_2DH(idg,NY,NX)*ppscal(idg)
          trcs_SubsurTransp_flx_2DH(idg,NY,NX) = trcs_SubsurTransp_flx_2DH(idg,NY,NX)+flux

          trcs_hydrloss_fast_flx_col(idg,NY,NX) = trcs_hydrloss_fast_flx_col(idg,NY,NX)+flux
          GasHydroSubsLoss_flx_col(idg,NY,NX)   = GasHydroSubsLoss_flx_col(idg,NY,NX)+flux
          RGasTranspFlxPrev_vr(idg,0,NY,NX)     = RGasTranspFlxPrev_vr(idg,0,NY,NX)+RGas_Disol_FlxMM_vr(idg,0,NY,NX)*ppscal(idg)

          !bottom
          flux=-Gas_AdvDif_FlxMM_3D(idg,3,NL_col(NY,NX)+1,NY,NX)*ppscal(idg)
          GasHydroSubsLoss_flx_col(idg,NY,NX)   = GasHydroSubsLoss_flx_col(idg,NY,NX)+flux
          trcs_hydrloss_fast_flx_col(idg,NY,NX) = trcs_hydrloss_fast_flx_col(idg,NY,NX)+flux
          trcs_drainage_fast_flx_col(idg,NY,NX) = trcs_drainage_fast_flx_col(idg,NY,NX) -flux

          trc_topsoil_flx_col(idg,NY,NX)        = trc_topsoil_flx_col(idg,NY,NX)+ppscal(idg)*Gas_AdvDif_FlxMM_3D(idg,3,NU_col(NY,NX),NY,NX)
          TranspNetSoil_flx2_col(idg,NY,NX)     = TranspNetSoil_flx2_col(idg,NY,NX) + ppscal(idg)*TranspNetSoil_fast_flxM_col(idg,NY,NX)
          TranspNetSoil_fast_flx_col(idg,NY,NX) = TranspNetSoil_fast_flx_col(idg,NY,NX)+ppscal(idg)*TranspNetSoil_fast_flxM_col(idg,NY,NX)

          !soil column
          DO L=NU_col(NY,NX),NL_col(NY,NX)
            RGasTranspFlxPrev_vr(idg,L,NY,NX) = RGasTranspFlxPrev_vr(idg,L,NY,NX)+ Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)*ppscal(idg)
            TranspNetSoil_flx_col(idg,NY,NX)  = TranspNetSoil_flx_col(idg,NY,NX)+ppscal(idg)*Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)
          ENDDO
        endif
      ENDDO
    ENDDO
  ENDDO  
  call PrintInfo('end '//subname)
  end subroutine AccumFastFluxesMM

!------------------------------------------------------------------------------------------
  subroutine ZeroTracerFluxMM(I,J,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J,NHE,NHW,NVS,NVN

  Gas_AdvDif_FlxMM_vr(idg_beg:idg_NH3,:,:,:)       = 0._r8
  Gas_AdvDif_FlxMM_3D(idg_beg:idg_NH3,:,:,:,:)     = 0.0_r8
  Gas_AdvDif_FlxMM_2DH(idg_beg:idg_NH3,:,:)        = 0._r8
  TranspNetSoil_fast_flxM_col(idg_beg:idg_NH3,:,:) = 0._r8
  end subroutine ZeroTracerFluxMM

!------------------------------------------------------------------------------------------

  subroutine LitterGasVolatilDissolMM(I,J,M,NHE,NHW,NVS,NVN)
  !
  !Description:
  !Gas dissolution flux into litter layer
  !This differs from the slow process in subroutine LitterAtmosExchangeM.
  !where the latter used aqueous diffusivity. 
  !These two can be viewed as non-eqiulbrium version of the problem tackled in
  !Tang and Riley (2013), Hydrol. Earth Syst. Sci, 17, 873–893, https://doi.org/10.5194/hess-17-873-2013.
  !
  implicit none

  integer, intent(in) :: I,J, M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  integer :: NY, NX
  real(r8) :: VOLGas
  real(r8) :: trc_gascl
  integer  :: idg
!
  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      IF(VGeomLayer_vr(0,NY,NX).GT.ZEROS2(NY,NX) .AND. VLsoiAirPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX) &
        .AND. VLWatMicPM_vr(M,0,NY,NX).GT.ZEROS2(NY,NX))THEN

        do idg=idg_beg,idg_NH3
          trcg_VLWatMicP_vr(idg,0,NY,NX)  = VLWatMicPM_vr(M,0,NY,NX)*GasSolbility_vr(idg,0,NY,NX)
          trc_gascl                       = trcg_gascl_vr(idg,0,NY,NX)*VLsoiAirPM_vr(M,0,NY,NX)
          VOLGas                          = trcg_VLWatMicP_vr(idg,0,NY,NX)+VLsoiAirPM_vr(M,0,NY,NX)
          RGas_Disol_FlxMM_vr(idg,0,NY,NX) = DiffusivitySolutEffM_vr(M,0,NY,NX) &
            *(AZMAX1(trc_gascl)*trcg_VLWatMicP_vr(idg,0,NY,NX) &
            - (trcs_solml2_vr(idg,0,NY,NX)-0._r8*trcs_solml_drib_vr(idg,0,NY,NX))*VLsoiAirPM_vr(M,0,NY,NX))/VOLGas
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  end subroutine LitterGasVolatilDissolMM

!------------------------------------------------------------------------------------------

  subroutine SurfSoilFluxGasDifAdvMM(M,NHE,NHW,NVS,NVN)
!
! DESCRIPTION:
! surface soil gaseous diffusion, advection, dissolution & volatilization
  implicit none

  integer, intent(in) :: M, NHE,NHW,NVS,NVN

  character(len=*), parameter :: subname='SurfSoilFluxGasDifAdvMM'
  real(r8) :: VFLW
  integer :: idg
  integer :: NY, NX

  call PrintInfo('beg '//subname)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      !Only deal with air-filled grids
      IF(FracAirFilledSoilPoreM_vr(M,NU_col(NY,NX),NY,NX).GT.AirFillPore_Min &
        .AND. SoilBulkDensity_vr(NU_col(NY,NX),NY,NX).GT.ZERO)THEN

        call Atm2TopSoilGasDifussionMM(M,NY,NX)

        call Atm2TopSoilGasAdvectionMM(M,NY,NX)

      ENDIF
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine SurfSoilFluxGasDifAdvMM

! ----------------------------------------------------------------------
  subroutine Atm2TopSoilGasAdvectionMM(M,NY,NX)
  !
  !Description
  !It assumes that air is not compressible, and the total soil pore volume 
  !is changing slowly (if it does change).
  implicit none
  integer, intent(in) :: M,NY,NX
  real(r8) :: VFLW
  real(r8) :: RGas_Adv_flxMM(idg_beg:idg_NH3)
  integer :: idg

  !topsoil gains water, so gases are squeezed out, and soil loses tracer
  IF(WaterFlow2SoilMM_3D(3,NU_col(NY,NX),NY,NX).GT.0.0_r8)THEN    
    IF(VLsoiAirPM_vr(M,NU_col(NY,NX),NY,NX).GT.ZEROS2(NY,NX))THEN
      VFLW=-AZMAX1(AMIN1(VFLWX,WaterFlow2SoilMM_3D(3,NU_col(NY,NX),NY,NX)/VLsoiAirPM_vr(M,NU_col(NY,NX),NY,NX)))
    ELSE
      VFLW=-VFLWX
    ENDIF
    
    DO idg=idg_beg,idg_NH3
      RGas_Adv_flxMM(idg)=VFLW*AZMAX1(trcg_gasml2_vr(idg,NU_col(NY,NX),NY,NX))
    ENDDO
    !topsoil gains tracers    
  ELSE
    DO idg=idg_beg,idg_NH3
      RGas_Adv_flxMM(idg)=-WaterFlow2SoilMM_3D(3,NU_col(NY,NX),NY,NX)*AtmGasCgperm3_col(idg,NY,NX)
    ENDDO
  ENDIF
  !
  !     TOTAL SOIL GAS FLUX + CONVECTIVE FLUX
  !
  DO idg=idg_beg,idg_NH3
    Gas_AdvDif_FlxMM_3D(idg,3,NU_col(NY,NX),NY,NX)=Gas_AdvDif_FlxMM_3D(idg,3,NU_col(NY,NX),NY,NX)+RGas_Adv_flxMM(idg)
  ENDDO
  end subroutine Atm2TopSoilGasAdvectionMM

! ----------------------------------------------------------------------
  subroutine GasTransportMM(M,N,N1,N2,N3,N4,N5,N6)

  !Description
  !do grid to grid gas transport (diffusion+advection)
  !the direction can be horizontal or vertical, DEPENDING
  !on the grid configuration.
  implicit none
  integer, intent(in) :: M         !iteration id
  integer, intent(in) :: N         !direction 
  integer, intent(in) :: N1,N2,N3  !source grid
  integer, intent(in) :: N4,N5,N6  !destination grid
  character(len=*), parameter :: subname='GasTransportMM'
  integer :: idg

  call PrintInfo('beg '//subname)

  IF(FracAirFilledSoilPoreM_vr(M,N3,N2,N1).GT.AirFillPore_Min      &    !grid supports gas
    .AND. FracAirFilledSoilPoreM_vr(M,N6,N5,N4).GT.AirFillPore_Min &    !grid supports gas  
    .AND. VLsoiAirPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1)               &    !source grid has significant air volume
    .AND. VLsoiAirPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN               !dest grid has significant air volume

    ! TOTAL SOIL GAS FLUX FROM DIFFUSIVE
    call GasDiffusionMM(M,N,N1,N2,N3,N4,N5,N6)

    ! TOTAL SOIL GAS FLUX FROM CONVECTIVE FLUX
    call UpstreamGasAdvectionMM(M,N,N1,N2,N3,N4,N5,N6)

  ENDIF
  call PrintInfo('end '//subname)
  end subroutine GasTransportMM

! ----------------------------------------------------------------------
  subroutine GasDissolutionMM(I,J,M,N,N4,N5,N6)
  !
  !Description
  !Compute gas dissolution 
  implicit none
  integer, intent(in) :: I,J,M,N
  integer, intent(in) :: N4,N5,N6
  character(len=*), parameter :: subname='GasDissolutionMM'
  integer :: idg
  real(r8) :: NH3GasMassZone
!
!
  call PrintInfo('beg '//subname)
  RGas_Disol_FlxMM_vr(idg_beg:idg_end,N6,N5,N4)=0._r8
  IF(VLWatMicPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN
    IF(FracAirFilledSoilPoreM_vr(M,N6,N5,N4).GT.AirFillPore_Min)THEN

      do idg=idg_beg,idg_NH3-1
        RGas_Disol_FlxMM_vr(idg,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
         (AZMAX1(trcg_gasml2_vr(idg,N6,N5,N4))*trcg_VLWatMicP_vr(idg,N6,N5,N4) &
          -(trcs_solml2_vr(idg,N6,N5,N4)-0._r8*trcs_solml_drib_vr(idg,N6,N5,N4))*VLsoiAirPM_vr(M,N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg,N6,N5,N4)+VLsoiAirPM_vr(M,N6,N5,N4))
      enddo    

      IF(VLsoiAirPMA_vr(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXA_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        ! Both dissolved zones exchange with one gas pool. Allocate its mass by
        ! zone air volume so the zone-specific equation uses the common gas concentration.
        NH3GasMassZone=AZMAX1(trcg_gasml2_vr(idg_NH3,N6,N5,N4)) &
          *VLsoiAirPMA_vr(N6,N5,N4)/VLsoiAirPM_vr(M,N6,N5,N4)
        RGas_Disol_FlxMM_vr(idg_NH3,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
         (NH3GasMassZone*trcg_VLWatMicP_vr(idg_NH3,N6,N5,N4) &
          -(trcs_solml2_vr(idg_NH3,N6,N5,N4)-0._r8*trcs_solml_drib_vr(idg_NH3,N6,N5,N4))*VLsoiAirPMA_vr(N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg_NH3,N6,N5,N4)+VLsoiAirPMA_vr(N6,N5,N4))
      ENDIF

      IF(VLsoiAirPMB_vr(N6,N5,N4).GT.ZEROS2(N5,N4).AND.VLWatMicPXB_vr(N6,N5,N4).GT.ZEROS2(N5,N4))THEN
        NH3GasMassZone=AZMAX1(trcg_gasml2_vr(idg_NH3,N6,N5,N4)) &
          *VLsoiAirPMB_vr(N6,N5,N4)/VLsoiAirPM_vr(M,N6,N5,N4)
        RGas_Disol_FlxMM_vr(idg_NH3B,N6,N5,N4)=DiffusivitySolutEffM_vr(M,N6,N5,N4)* &
          (NH3GasMassZone*trcg_VLWatMicP_vr(idg_NH3B,N6,N5,N4) &
          -(trcs_solml2_vr(idg_NH3B,N6,N5,N4)-0._r8*trcs_solml_drib_vr(idg_NH3B,N6,N5,N4))*VLsoiAirPMB_vr(N6,N5,N4)) &
          /(trcg_VLWatMicP_vr(idg_NH3B,N6,N5,N4)+VLsoiAirPMB_vr(N6,N5,N4))
      ENDIF

    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine GasDissolutionMM

! ----------------------------------------------------------------------
  subroutine UpstreamGasAdvectionMM(M,N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer, intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  character(len=*), parameter :: subname='UpstreamGasAdvectionMM'
  real(r8) :: RGasAdv
  real(r8) :: VFLW,FLQW
  integer :: idg
!
!     CONVECTIVE GAS TRANSFER DRIVEN BY SOIL WATER FLUXES
!     FROM 'WATSUB' AND GAS CONCENTRATIONS IN THE ADJACENT GRID CELLS
!     DEPENDING ON WATER FLUX DIRECTION
!
!     by assuming volume conservation, gases and water flow in opposite direction
!
  call PrintInfo('beg '//subname)
  FLQW=WaterFlow2SoilMM_3D(N,SoilFaceSlot(N,N3,N6),N5,N4)
  
  !water flow into grid (N6,N5,N4), so gas out of grid
  IF(FLQW.GT.0.0_r8)THEN
    !dest grid has air-filled pores
    IF(VLsoiAirPM_vr(M,N6,N5,N4).GT.ZEROS2(N5,N4))THEN   
      VFLW=-AZMAX1(AMIN1(VFLWX,FLQW/VLsoiAirPM_vr(M,N6,N5,N4)))   !negative flow 
    !dest grid is aturated  
    ELSE
      VFLW=-VFLWX
    ENDIF

    DO idg=idg_beg,idg_NH3
      RGasAdv                             = VFLW*AZMAX1(trcg_gasml2_vr(idg,N6,N5,N4))
      Gas_AdvDif_FlxMM_3D(idg,N,SoilFaceSlot(N,N3,N6),N5,N4) = Gas_AdvDif_FlxMM_3D(idg,N,SoilFaceSlot(N,N3,N6),N5,N4)+RGasAdv
    ENDDO
    !water flow out of source grid, gas into grid (N6,N5,N4)  
  ELSE
    IF(VLsoiAirPM_vr(M,N3,N2,N1).GT.ZEROS2(N2,N1))THEN
      VFLW=-AZMIN1(AMAX1(-VFLWX,FLQW/VLsoiAirPM_vr(M,N3,N2,N1)))
    ELSE
      VFLW=VFLWX
    ENDIF
    DO idg=idg_beg,idg_NH3
      RGasAdv                             = VFLW*AZMAX1(trcg_gasml2_vr(idg,N3,N2,N1))
      Gas_AdvDif_FlxMM_3D(idg,N,SoilFaceSlot(N,N3,N6),N5,N4) = Gas_AdvDif_FlxMM_3D(idg,N,SoilFaceSlot(N,N3,N6),N5,N4)+RGasAdv
    ENDDO
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine UpstreamGasAdvectionMM

!------------------------------------------------------------------------------------------

  subroutine XBoundaryFluxMM(I,J,M,NHW,NHE,NVN,NVS)
  !
  !Do transport across the boundaries
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M,NHW, NHE, NVN, NVS

  character(len=*), parameter :: subname='XBoundaryFluxMM'
  integer :: NY,NX,L
  integer :: NN,N
  integer :: M1,M2,M3,M4,M5,M6          !boundary exchange
  integer :: LL
  call PrintInfo('beg '//subname)
  !
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      D9585: DO L=NU_col(NY,NX),NL_col(NY,NX)
!
!     LOCATE ALL EXTERNAL BOUNDARIES AND SET BOUNDARY CONDITIONS
!     ENTERED IN 'READS'
!
        D9580: DO  N=FlowDirIndicator_col(NY,NX),3
          D9575: DO  NN=1,2
            IF(N.EQ.iWestEastDirection)THEN
              !WEST-EAST
              IF(NN.EQ.iFront)THEN   !eastward                
                IF(NX.EQ.NHE)THEN !eastern boundary
                  M4 = NX+1;M5 = NY;M6 = L
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN  !west                
                IF(NX.EQ.NHW)THEN  !western boundary
                  M4 = NX;M5 = NY;M6 = L
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iNorthSouthDirection)THEN              
              IF(NN.EQ.iFront)THEN  
                IF(NY.EQ.NVS)THEN   ! southern boundary
                  M4 = NX;M5 = NY+1;M6 = L   !target grid
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN  !north                
                IF(NY.EQ.NVN)THEN ! northern boundary
                  M4 = NX;M5 = NY;M6 = L !target
                ELSE
                  cycle
                ENDIF
              ENDIF
            ELSEIF(N.EQ.iVerticalDirection)THEN !vertical
              IF(NN.EQ.iFront)THEN                
                IF(L.EQ.NL_col(NY,NX))THEN      !lower boundary
                  M4 = NX;M5 = NY;M6 = L+1  !target grid
                ELSE
                  cycle
                ENDIF
              ELSEIF(NN.EQ.iBehind)THEN
                !nothing for the upper boundary
                cycle
              ENDIF
            ENDIF            
            M1 = NX;M2 = NY;M3 = L    !source grid    
            !        
            call XBoundaryTracerFlowMM(N,NN,M,M1,M2,M3,M4,M5,M6)
          ENDDO D9575
        ENDDO D9580
      ENDDO D9585

    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine XBoundaryFluxMM

! ----------------------------------------------------------------------
  subroutine XBoundaryTracerFlowMM(N,NN,M,M1,M2,M3,M4,M5,M6)
  implicit none

  integer, intent(in) :: N, NN, M
  integer, intent(in) :: M1, M2,M3  !source grid
  integer, intent(in) :: M4, M5,M6  !dest grid

  character(len=*), parameter :: subname='XBoundaryTracerFlowMM'
  real(r8) :: FLGM,FQRM,VFLW
  integer :: K,idg

! begin_execution
  call PrintInfo('beg '//subname)

!   dt_GasCyc=1/NPT, NPT is number of gas iterations per M

  FLGM=WaterFlow2SoilMM_3D(N,M6,M5,M4)
  
  !make sure out of grid (M3,M2,M1)
  IF(NN.EQ.iFront .AND. FLGM.LT.0.0_r8           &  !water coming in from front, so gas goes out 
    .OR. (NN.EQ.iBehind .AND. FLGM.GT.0.0_r8))THEN  !water coming in from behind, so gas goes out  
    IF(VLsoiAirPM_vr(M,M3,M2,M1).GT.ZEROS2(M2,M1))THEN
      VFLW=-AMAX1(-VFLWX,AMIN1(VFLWX,FLGM/VLsoiAirPM_vr(M,M3,M2,M1)))
    ELSE
      VFLW=0.0_r8
    ENDIF
    !
    ! add gas advection flux
    !
    DO idg=idg_beg,idg_NH3
      Gas_AdvDif_FlxMM_3D(idg,N,M6,M5,M4) = Gas_AdvDif_FlxMM_3D(idg,N,M6,M5,M4)+VFLW*AZMAX1(trcg_gasml2_vr(idg,M3,M2,M1))
    ENDDO
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine XBoundaryTracerFlowMM
!------------------------------------------------------------------------------------------

  subroutine NetTracerFlowXSoilPoresMM(N,M,N1,N2,N3,N4,N5,N6)
  !
  !Description
  !Vertical or lateral flux between grid grids
  implicit none

  integer, intent(in) :: N,M,N1,N2,N3,N4,N5,N6
  integer :: K,LL,idg,idom
  integer :: first_face,last_face
  real(r8) :: flux
!
!     NET GAS FLUX
!
  call SoilFaceIncomingRange(N,N3,N2,N1,first_face,last_face)
  IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
    DO idg=idg_beg,idg_NH3
      flux =sum(Gas_AdvDif_FlxMM_3D(idg,N,first_face:last_face,N2,N1))-Gas_AdvDif_FlxMM_3D(idg,N,SoilFaceSlot(N,N3,N6),N5,N4)
      Gas_AdvDif_FlxMM_vr(idg,N3,N2,N1)=Gas_AdvDif_FlxMM_vr(idg,N3,N2,N1)+flux

      if(N.NE.iVerticalDirection)then
        Gas_AdvDif_FlxMM_2DH(idg,N2,N1)=Gas_AdvDif_FlxMM_2DH(idg,N2,N1)+flux
      endif  
    ENDDO
  ENDIF
  end subroutine NetTracerFlowXSoilPoresMM

! ----------------------------------------------------------------------

  subroutine GasDiffusionMM(M,N,N1,N2,N3,N4,N5,N6)

  implicit none
  integer , intent(in) :: M,N,N1,N2,N3,N4,N5,N6
  real(r8) :: trc_gasc1(idg_beg:idg_end)
  real(r8) :: trc_gasc2(idg_beg:idg_end)
  character(len=*), parameter :: subname='GasDiffusionMM'
  real(r8) :: CNDC1   !conductance in source grid
  real(r8) :: CNDC2   !conductance in destination grid
  real(r8) :: DFLG2   !air-filled porosity effect on gaseous diffusivity in source
  real(r8) :: DFLGL   !air-filled porosity effect on gaseous diffusivity in destination
  real(r8) :: FaceArea
  integer :: idg
  real(r8):: GasDifuscoefMM_3D

  call PrintInfo('beg '//subname)

  ! The traversal fixes (N3,N2,N1) as the face owner, independently of the
  ! flux sign. Use its area on both sides, matching water/solute transport.
  ! This ownership convention does not calculate unequal layers' geometric overlap.
  FaceArea=AREA_3D(N,N3,N2,N1)
  DFLG2=2.0_r8*POROQ*AZMAX1(FracAirFilledSoilPoreM_vr(M,N3,N2,N1))**2/POROS_vr(N3,N2,N1)*FaceArea/DLYR_3D(N,N3,N2,N1)

  DFLGL=2.0_r8*POROQ*AZMAX1(FracAirFilledSoilPoreM_vr(M,N6,N5,N4))**2/POROS_vr(N6,N5,N4)*FaceArea/DLYR_3D(N,N6,N5,N4)
  !
  DO idg=idg_beg,idg_NH3
    !     GASOUS CONDUCTANCES
    CNDC1             = DFLG2*GasDifctScaledMM_vr(idg,N3,N2,N1)
    CNDC2             = DFLGL*GasDifctScaledMM_vr(idg,N6,N5,N4)
    GasDifuscoefMM_3D = (CNDC1*CNDC2)/(CNDC1+CNDC2)
    !
    !diffusion flux from (N3,N2,N1) into (N6,N5,N4)
    !
    trc_gasc1(idg)                      = AZMAX1(trcg_gasml2_vr(idg,N3,N2,N1)/VLsoiAirPM_vr(M,N3,N2,N1))
    trc_gasc2(idg)                      = AZMAX1(trcg_gasml2_vr(idg,N6,N5,N4)/VLsoiAirPM_vr(M,N6,N5,N4))
    Gas_AdvDif_FlxMM_3D(idg,N,SoilFaceSlot(N,N3,N6),N5,N4) = Gas_AdvDif_FlxMM_3D(idg,N,SoilFaceSlot(N,N3,N6),N5,N4)+GasDifuscoefMM_3D*(trc_gasc1(idg)-trc_gasc2(idg))
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine GasDiffusionMM
! ----------------------------------------------------------------------
  subroutine Atm2TopSoilGasDifussionMM(M,NY,NX)
  !
  !Description:
  !Gas diffusion between atmosphere and soil. (>0 into soil)
  implicit none
  integer, intent(in) :: M,NY,NX
  character(len=*), parameter :: subname='Atm2TopSoilGasDifussionMM'
  real(r8) :: DFLG2
  real(r8) :: trcg_cl2
  real(r8) :: DGQ_cef           !effective diffusivity for topsoil-atmosphere gas excchange
  real(r8) :: GasDifuscoefMM_3D
  integer  :: idg

!     GASEOUS DIFFUSIVITIES
!
  call PrintInfo('beg '//subname)
  DFLG2=AZMAX1(FracAirFilledSoilPoreM_vr(M,NU_col(NY,NX),NY,NX))*POROQ &
    *FracAirFilledSoilPoreM_vr(M,NU_col(NY,NX),NY,NX)/POROS_vr(NU_col(NY,NX),NY,NX) &
    *AREA_3D(3,NU_col(NY,NX),NY,NX)/AMAX1(ZERO2,DLYR_3D(3,NU_col(NY,NX),NY,NX))

  DO idg=idg_beg,idg_NH3
    GasDifuscoefMM_3D=DFLG2*GasDifctScaledMM_vr(idg,NU_col(NY,NX),NY,NX)
!
!     SURFACE GAS CONCENTRATIONS
!
    trcg_cl2=AZMAX1(trcg_gasml2_vr(idg,NU_col(NY,NX),NY,NX))/VLsoiAirPM_vr(M,NU_col(NY,NX),NY,NX)
!
!     EQUILIBRIUM CONCENTRATIONS AT SOIL SURFACE AT WHICH
!     GASEOUS DIFFUSION THROUGH SOIL SURFACE LAYER = GASEOUS
!     DIFFUSION THROUGH ATMOSPHERE BOUNDARY LAYER CALCULATED
!     FROM GASEOUS DIFFUSIVITY AND BOUNDARY LAYER CONDUCTANCE
!
    DGQ_cef=GasDifuscoefMM_3D*PARGas_CefMM(idg,NY,NX)/(GasDifuscoefMM_3D+PARGas_CefMM(idg,NY,NX))

    Gas_AdvDif_FlxMM_3D(idg,3,NU_col(NY,NX),NY,NX)=Gas_AdvDif_FlxMM_3D(idg,3,NU_col(NY,NX),NY,NX)+DGQ_cef*(AtmGasCgperm3_col(idg,NY,NX)-trcg_cl2)
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine Atm2TopSoilGasDifussionMM

!------------------------------------------------------------------------------------------
  subroutine FastUpdateStateVarsMM(I,J,M,NHW, NHE, NVN, NVS,dpscal,pscal)
  implicit none
  integer, intent(in) :: I,J,M  
  integer, intent(in) :: NHW, NHE, NVN, NVS
  real(r8),intent(in) :: dpscal(idg_beg:idg_end)
  real(r8),optional,intent(inout) :: pscal(idg_beg:idg_end)
  
  if(present(pscal))then
    pscal=1.000_r8
    call UpdateStateVarsMM(I,J,M,NHW, NHE, NVN, NVS,dpscal,pscal)
  else
    call UpdateStateVarsMM(I,J,M,NHW, NHE, NVN, NVS,dpscal)
  endif  
  end subroutine FastUpdateStateVarsMM
!------------------------------------------------------------------------------------------
  subroutine UpdateStateVarsMM(I,J,M,NHW, NHE, NVN, NVS,dpscal,pscal1)
  implicit none
  integer, intent(in) :: I,J,M  
  integer, intent(in) :: NHW, NHE, NVN, NVS
  real(r8),intent(in) :: dpscal(idg_beg:idg_end)
  real(r8), optional, intent(inout) :: pscal1(idg_beg:idg_end)

  character(len=*), parameter :: subname='UpdateStateVarsMM'
  integer :: NY,NX
  integer :: L,idg
  real(r8) :: pscal(idg_beg:idg_end)  
  real(r8) :: flux
  logical :: lflux,lfupdate

  call PrintInfo('beg '//subname)
  if(present(pscal1))then
    pscal=pscal1
    lfupdate=.false.
  else
    pscal=1.000_r8
    lfupdate=.true.
  endif  

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      !in litter layer
      DO idg=idg_beg,idg_NH3        
        flux = RGas_Disol_FlxMM_vr(idg,0,NY,NX)
        lflux=.true. .or. .not.isclose(flux,0._r8)
        if(lflux .and. dpscal(idg)>tiny_p)then  
          flux=flux*dpscal(idg)
          call get_flux_scalar(trcs_solml2_vr(idg,0,NY,NX),flux,trcs_solml_vr(idg,0,NY,NX),pscal(idg))    
          if(lfupdate)then
            trcg_mass3_fast_col(idg,NY,NX)         = trcg_mass3_fast_col(idg,NY,NX)+trcs_solml2_vr(idg,0,NY,NX)
            trcs_solml2_vr(idg,0,NY,NX)            = trcs_solml_vr(idg,0,NY,NX)
            TranspNetSoil_fast_flxM_col(idg,NY,NX) = TranspNetSoil_fast_flxM_col(idg,NY,NX)+RGas_Disol_FlxMM_vr(idg,0,NY,NX)
          endif  

        endif
      ENDDO
       
      !in the soil      
      !
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        !other gases are taken from aqueous phases
        DO idg=idg_beg,idg_NH3-1          
          flux = RGas_Disol_FlxMM_vr(idg,L,NY,NX)
          lflux=.true. .or. .not.isclose(flux,0._r8)

          if(lflux .and. dpscal(idg)>tiny_p)then  
            flux=flux*dpscal(idg)
            call get_flux_scalar(trcs_solml2_vr(idg,L,NY,NX),flux, trcs_solml_vr(idg,L,NY,NX),pscal(idg))
            
            if(lfupdate)then
              trcg_mass3_fast_col(idg,NY,NX) = trcg_mass3_fast_col(idg,NY,NX)+trcs_solml2_vr(idg,L,NY,NX)
              trcs_solml2_vr(idg,L,NY,NX)    = trcs_solml_vr(idg,L,NY,NX)
            endif  

            flux=-RGas_Disol_FlxMM_vr(idg,L,NY,NX)+Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)
            flux=flux*dpscal(idg)
            call get_flux_scalar(trcg_gasml2_vr(idg,L,NY,NX),flux, trcg_gasml_vr(idg,L,NY,NX),pscal(idg))

            if(lfupdate)then
              trcg_mass3_fast_col(idg,NY,NX)         = trcg_mass3_fast_col(idg,NY,NX)+trcg_gasml2_vr(idg,L,NY,NX)
              trcg_gasml2_vr(idg,L,NY,NX)            = trcg_gasml_vr(idg,L,NY,NX)
              TranspNetSoil_fast_flxM_col(idg,NY,NX) = TranspNetSoil_fast_flxM_col(idg,NY,NX) + Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)
            endif  
          endif
        ENDDO

        !nonband/band NH3
        DO idg=idg_NH3,idg_NH3B          
          flux=RGas_Disol_FlxMM_vr(idg,L,NY,NX)     
          lflux=.true. .or. .not.isclose(flux,0._r8)     
          if(lflux .and. dpscal(idg)>tiny_p)then  
            flux=flux*dpscal(idg)
            call get_flux_scalar(trcs_solml2_vr(idg,L,NY,NX),flux, trcs_solml_vr(idg,L,NY,NX),pscal(idg_NH3))
            if(lfupdate)then
              trcg_mass3_fast_col(idg_NH3,NY,NX) = trcg_mass3_fast_col(idg_NH3,NY,NX)+trcs_solml2_vr(idg,L,NY,NX)
              trcs_solml2_vr(idg,L,NY,NX)        = trcs_solml_vr(idg,L,NY,NX)
            endif  
          endif
        ENDDO

        !gaseous NH3 is consumed through geochemistry
        idg  = idg_NH3
        flux = -RGas_Disol_FlxMM_vr(idg,L,NY,NX)- RGas_Disol_FlxMM_vr(idg_NH3B,L,NY,NX)+Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)
           
        lflux=.true. .or. .not.isclose(flux,0._r8) 

        if(lflux .and. dpscal(idg)>tiny_p)then  
          flux=flux*dpscal(idg)
          call get_flux_scalar(trcg_gasml2_vr(idg,L,NY,NX),flux, trcg_gasml_vr(idg,L,NY,NX),pscal(idg))  
          if(lfupdate)then
            trcg_mass3_fast_col(idg,NY,NX)         = trcg_mass3_fast_col(idg,NY,NX)+trcg_gasml2_vr(idg,L,NY,NX)
            trcg_gasml2_vr(idg,L,NY,NX)            = trcg_gasml_vr(idg,L,NY,NX)
            TranspNetSoil_fast_flxM_col(idg,NY,NX) = TranspNetSoil_fast_flxM_col(idg,NY,NX)+Gas_AdvDif_FlxMM_vr(idg,L,NY,NX)
          endif
          pscal(idg_NH3B)=pscal(idg_NH3)  
        endif

      ENDDO
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  if(present(pscal1))pscal1=pscal  
  end subroutine UpdateStateVarsMM
!------------------------------------------------------------------------------------------

  subroutine TracerFlowXGridsMM(I,J,M,MM,NHE,NHW,NVS,NVN)
!
! DESCRIPTION
! exchanges tracers within (gaseous vs aqueous phase) and between
! grid cells.
  implicit none
  integer, intent(in) :: I,J,M,MM
  integer, intent(in) :: NHE,NHW,NVS,NVN

  character(len=*), parameter :: subname='TracerFlowXGridsMM'
  real(r8) :: VFLW
  real(r8) :: VOLH2A,VOLH2B
  real(r8) :: VLWatMacPS,VOLWT

  integer :: iFalseEbu,N,L,K,LL,NY,NX
  integer :: N1,N2,N3  !source grid index
  integer :: N4,N5,N6  !dest grid index

! begin_execution
  call PrintInfo('beg '//subname)

!     N3,N2,N1=L,NY,NX of source grid cell
!     N6,N5,N4=L,NY,NX of destination grid cell
!
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      iFalseEbu=ifalse
      D125: DO L=NU_col(NY,NX),NL_col(NY,NX)
        !source
        N1=NX;N2=NY;N3=L    
        IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))then
          call GasDissolutionMM(I,J,M,N,N1,N2,N3)
        endif  
        !
        !     LOCATE INTERNAL BOUNDARIES BETWEEN ADJACENT GRID CELLS
        !
        D120: DO N=FlowDirIndicator_col(N2,N1),3

          IF(N.EQ.iWestEastDirection)THEN
            IF(NX.EQ.NHE)THEN  !skip eastern boundary
              cycle
            ELSE
              N4 = NX+1;N5 = NY;N6 = L
            ENDIF
          ELSEIF(N.EQ.iNorthSouthDirection)THEN            
            IF(NY.EQ.NVS)THEN  !skip southern boundary
              cycle
            ELSE
              N4 = NX;N5 = NY+1;N6 = L
            ENDIF
          ELSEIF(N.EQ.iVerticalDirection)THEN            
            IF(L.EQ.NL_col(NY,NX))THEN  !skip bottom boundary
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
          ENDDO
          if(LL>NL_col(N5,N4))cycle D120
          !
          !     SOLUTE FLUXES BETWEEN ADJACENT GRID CELLS FROM
          !     WATER CONTEnsolutes AND WATER FLUXES 'WaterFlow2SoilMM_3D' FROM 'WATSUB'
          !
          IF(VLSoilPoreMicP_vr(N3,N2,N1).GT.ZEROS2(N2,N1))THEN
            IF(N3.GE.NU_col(N2,N1) .AND. N6.GE.NU_col(N5,N4) .AND. N3.LE.NL_col(N2,N1) .AND. N6.LE.NL_col(N5,N4))THEN               
              call GasTransportMM(M,N,N1,N2,N3,N4,N5,N6)
            ENDIF           
          ENDIF
        ENDDO D120
      ENDDO D125
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine TracerFlowXGridsMM

end module TranspNoSaltFastMod
