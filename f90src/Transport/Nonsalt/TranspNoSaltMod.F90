module TranspNoSaltMod
!!
! Description:
!
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,       only: destroy, endrun
  USE MiniMathMod,      ONLY: AZMAX1,  fixnegmass, flux_mass_limiter, AZERO, safe_adb
  use TracerPropMod,    only: MolecularWeight
  use EcoSiMParDataMod, only: micpar
  use EcoSIMCtrlMod,    only: iVerbLevel,ldo_transpt_bubbling
  use NumericalAuxMod
  use TranspNoSaltFastMod
  use InitNoSaltTransportMod
  use TranspNoSaltSlowMod
  use DebugToolMod
  use SOMDataType
  use ChemTranspDataType
  use GridConsts
  use EcoSIMSolverPar
  use SoilPhysDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use SoilBGCDataType
  use ClimForcDataType
  use FertilizerDataType
  use EcosimConst
  use SurfLitterDataType
  use SnowDataType
  use SurfSoilDataType
  use LandSurfDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use PlantDataRateType
  use GridDataType
  use TranspNoSaltDataMod
  use IrrigationDataType
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  real(r8), parameter :: tinyval=1.e-13_r8
  public :: TranspNoSalt
  public :: InitTranspNoSalt,DestructTranspNoSalt
  contains

  subroutine InitTranspNoSalt()
  implicit none

  call InitTransfrData

  end subroutine InitTranspNoSalt
!------------------------------------------------------------------------------------------
  subroutine DestructTranspNoSalt
  implicit none

  call DestructTransfrData

  end subroutine DestructTranspNoSalt
!------------------------------------------------------------------------------------------

  subroutine EnterMassCheck(NHW,NHE,NVN,NVS)
  !
  !Description:
  !Initialize mass balance check before staging no-salt transport
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX,idg,L,ids

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      trcs_mass_beg(:,NY,NX) = 0._r8
      errmass_fast(:,NY,NX)  = 0._r8
      errmass_slow(:,NY,NX)  = 0._r8
      trcs_solml_drib_beg_col(:,NY,NX)     = 0._r8
      trcs_solml_drib_soi_beg_col(:,NY,NX) = 0._r8
      trcs_netProd_lit_col(:,NY,NX)        = 0._r8

      do ids=ids_beg,ids_end
        trcs_solml_drib_beg_col(ids,NY,NX) = trcs_solml_drib_vr(ids,0,NY,NX)
      enddo
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO ids=ids_beg,ids_end
          trcs_solml_drib_soi_beg_col(ids,NY,NX) = trcs_solml_drib_soi_beg_col(ids,NY,NX)+trcs_solml_drib_vr(ids,L,NY,NX)
        ENDDO  
      ENDDO
      DO ids=ids_beg,ids_end
        trcs_solml_drib_beg_col(ids,NY,NX) = trcs_solml_drib_beg_col(ids,NY,NX)+ trcs_solml_drib_soi_beg_col(ids,NY,NX)
      ENDDO  

      DO idg=idg_beg,idg_NH3
        trcs_mass_litr_beg(idg,NY,NX)=trcs_solml_vr(idg,0,NY,NX)
        trcs_mass_snow_beg(idg,NY,NX)=0._r8
        trcs_mass_soil_beg(idg,NY,NX)=0._r8

        ! First-snow tracer mass exists before SnowpackLayering updates nsnol.
        DO L=1,JS
          trcs_mass_snow_beg(idg,NY,NX)=trcs_mass_snow_beg(idg,NY,NX)+trcg_solsml_snvr(idg,L,NY,NX)
        ENDDO
        
        DO L=NU_col(NY,NX),NL_col(NY,NX)          
          trcs_mass_soil_beg(idg,NY,NX) = trcs_mass_soil_beg(idg,NY,NX)+trcs_solml_vr(idg,L,NY,NX)+trcg_gasml_vr(idg,L,NY,NX)+ trcs_soHml_vr(idg,L,NY,NX)
        ENDDO      
      ENDDO

      idg=idg_NH3
      trcs_mass_litr_beg(idg,NY,NX) = trcs_mass_litr_beg(idg,NY,NX)+trcs_solml_vr(idg_NH3B,0,NY,NX)
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        trcs_mass_soil_beg(idg,NY,NX) = trcs_mass_soil_beg(idg,NY,NX)+trcs_solml_vr(idg_NH3B,L,NY,NX)+trcs_soHml_vr(idg_NH3B,L,NY,NX)
      ENDDO

      DO idg=idg_beg,idg_NH3
        trcs_mass_beg(idg,NY,NX)=trcs_mass_litr_beg(idg,NY,NX)+trcs_mass_snow_beg(idg,NY,NX)+trcs_mass_soil_beg(idg,NY,NX) 
      ENDDO  
    ENDDO
  ENDDO  
  end subroutine EnterMassCheck
!------------------------------------------------------------------------------------------

  subroutine ExitMassCheck(I,J,NHW,NHE,NVN,NVS,M)
  !
  !Description: 
  ! Overall mass conservation check for non-salt chemical transport
  ! fluxes inlcude:
  ! surface emission        : SurfGasEmiss_flx_col(idg,NY,NX), ebullition, surface diffusion and wet deposition
  ! hydrological loss       : GasHydroLoss_flx_col(idg,NY,NX) &
  ! biogeochemcal production: RGasNetProd_col(idg,NY,NX)
  ! uptake by plants        : trcs_Soil2plant_uptake_col(idg,NY,NX)

  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer, optional, intent(in) :: M

  integer :: NY,NX,idg,L,ids  
  real(r8) :: errmass,mass_litr,mass_snow,mass_soil
  real(r8) :: netflx2soil,delta_mass
  real(r8) :: trcs_mass_now(idg_beg:idg_NH3)
  real(r8) :: trcs_solml_drib_col(ids_beg:ids_end)
  real(r8) :: trcs_solml_drib_soil_col(ids_beg:ids_end)

  DO NX=NHW,NHE
    DO  NY=NVN,NVS

      trcs_mass_now(:)=0._r8

      trcs_solml_drib_col(:)   = 0._r8
      trcs_solml_drib_soil_col(:)=0._r8
      do ids=ids_beg,ids_end
        trcs_solml_drib_col(ids) = trcs_solml_drib_vr(ids,0,NY,NX)
      enddo
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO ids=ids_beg,ids_end
          trcs_solml_drib_soil_col(ids)=trcs_solml_drib_soil_col(ids)+trcs_solml_drib_vr(ids,L,NY,NX)
        ENDDO  
      ENDDO
      DO ids=ids_beg,ids_end
        trcs_solml_drib_col(ids) = trcs_solml_drib_col(ids)+trcs_solml_drib_soil_col(ids)
      enddo

      DO idg=idg_beg,idg_NH3
        SurfGasEmiss_flx_col(idg,NY,NX) =  trcg_ebu_flx_col(idg,NY,NX) &
          + GasDiff2Surf_flx_col(idg,NY,NX)+Gas_WetDeposit_flx_col(idg,NY,NX) 

        mass_litr = trcs_solml_vr(idg,0,NY,NX)
        mass_snow = 0._r8
        DO L=1,JS
          mass_snow=mass_snow+trcg_solsml_snvr(idg,L,NY,NX)
        ENDDO        

        mass_soil=0._r8
        DO L=NU_col(NY,NX),NL_col(NY,NX)
          mass_soil=mass_soil+trcg_gasml_vr(idg,L,NY,NX)+trcs_solml_vr(idg,L,NY,NX)+trcs_soHml_vr(idg,L,NY,NX)  
        ENDDO

        if(idg==idg_NH3)then
          !also include NH3B
          mass_litr=mass_litr+trcs_solml_vr(idg_NH3B,0,NY,NX)
          DO L=NU_col(NY,NX),NL_col(NY,NX)
            mass_soil=mass_soil+trcs_solml_vr(idg_NH3B,L,NY,NX)+trcs_soHml_vr(idg_NH3B,L,NY,NX)
          ENDDO
        endif

        trcs_mass_now(idg) = mass_snow+mass_litr+mass_soil
        delta_mass         = trcs_mass_now(idg)-trcs_mass_beg(idg,NY,NX)
        errmass            = -delta_mass+SurfGasEmiss_flx_col(idg,NY,NX)+GasHydroLoss_flx_col(idg,NY,NX) &
          +RGasNetProd_col(idg,NY,NX)+trcs_solml_drib_col(idg) 

        RGasNetProd_col(idg,NY,NX) = RGasNetProd_col(idg,NY,NX)+trcs_Soil2plant_uptake_col(idg,NY,NX)

        if(idg==idg_NH3)then
          RGasNetProd_col(idg,NY,NX) = RGasNetProd_col(idg,NY,NX)+trcs_Soil2plant_uptake_col(idg_NH3B,NY,NX)
          errmass=errmass+trcs_solml_drib_col(idg_NH3B) 
        endif
 
        if(abs(errmass)>1.e-4 .and. abs(safe_adb(errmass,delta_mass))>1.e-3_r8)then
          if(iVerbLevel==1 .or. abs(errmass)>1.e-4)then
            write(121,*)('-',L=1,50)
            if(present(M))then
              write(121,*)(I*1000+J)*10+M,trcs_names(idg),'total',M
            else
              write(121,*)(I*1000+J)*10,trcs_names(idg),'total'
            endif
            write(121,*)'NU   =',NU_col(NY,NX),trcs_solml_drib_vr(idg,1:NL_col(NY,NX),NY,NX)
            write(121,*)'errmass total, fast, slow=',errmass, errmass_fast(idg,NY,NX),errmass_slow(idg,NY,NX),'NY NX=',NY,NX
            write(121,*)'beg end delta mass       =',trcs_mass_beg(idg,NY,NX),trcs_mass_now(idg),-delta_mass
            write(121,*)'beg sno,litr,soil        =',trcs_mass_snow_beg(idg,NY,NX),trcs_mass_litr_beg(idg,NY,NX),trcs_mass_soil_beg(idg,NY,NX)
            write(121,*)'mass sno,litr,soil       =',mass_snow,mass_litr,mass_soil
            write(121,*)'Delta mass sno, litr,soil=',mass_snow-trcs_mass_snow_beg(idg,NY,NX),mass_litr-trcs_mass_litr_beg(idg,NY,NX), &
              mass_soil-trcs_mass_soil_beg(idg,NY,NX)
            write(121,*)'mass snow                =',trcg_solsml_snvr(idg,1:JS,NY,NX)
            write(121,*)'------------------------------------------------------------------------------'            
            write(121,*)'dif to litr              =',AtmGasDiff2Litr_flx_col(idg,NY,NX)
            write(121,*)'wet dep to litr.         =',Gas_WetDepo2Litr_col(idg,NY,NX)
            write(121,*)'gas litr2soil            =',Gas_litr2Soil_flx_col(idg,NY,NX)  
            write(121,*)'snow2litr                =',trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)
            write(121,*)'runoff loss              =',GasHydroLoss_litr_flx_col(idg,NY,NX) 
            write(121,*)'dribble litr             =',trcs_solml_drib_vr(idg,0,NY,NX)           
            write(121,*)'net prod litr            =',trcs_netProd_lit_col(idg,NY,NX) 
            write(121,*)'netflx2litr              =',trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)+AtmGasDiff2Litr_flx_col(idg,NY,NX) &
              +Gas_WetDepo2Litr_col(idg,NY,NX)-Gas_litr2Soil_flx_col(idg,NY,NX)+GasHydroLoss_litr_flx_col(idg,NY,NX)&
              +trcs_netProd_lit_col(idg,NY,NX) 
            write(121,*)'err in litr               =', trcg_AquaADV_Snow2Litr_flx(idg,NY,NX)+AtmGasDiff2Litr_flx_col(idg,NY,NX) &
              +Gas_WetDepo2Litr_col(idg,NY,NX)-Gas_litr2Soil_flx_col(idg,NY,NX)+GasHydroLoss_litr_flx_col(idg,NY,NX) &
              -mass_litr+trcs_mass_litr_beg(idg,NY,NX)+trcs_solml_drib_vr(idg,0,NY,NX)+trcs_netProd_lit_col(idg,NY,NX) 
            write(121,*)'------------------------------------------------------------------------------'
            write(121,*)'nsnol_col                 =',nsnol_col(NY,NX)
            write(121,*)'err in snow               =',Gas_WetDepo2Snow_col(idg,NY,NX)-Gas_Snowloss_flx_col(idg,NY,NX)-mass_snow+trcs_mass_snow_beg(idg,NY,NX)
            write(121,*)'snow depo loss drift      =',Gas_WetDepo2Snow_col(idg,NY,NX),Gas_Snowloss_flx_col(idg,NY,NX),trcg_SnowDrift_flx_col(idg,NY,NX)
            write(121,*)'ntransp2soil              =',TranspNetSoil_flx_col(idg,NY,NX),TranspNetSoil_flx2_col(idg,NY,NX),TranspNetSoil_flx_col(idg,NY,NX)-TranspNetSoil_flx2_col(idg,NY,NX)              
            write(121,*)'drain                     =',trcs_drainage_flx_col(idg,NY,NX)  
            write(121,*)'ebu                       =',trcg_ebu_flx_col(idg,NY,NX)
            write(121,*)'dif                       =',GasDiff2Surf_flx_col(idg,NY,NX)                                    
            write(121,*)'------------------------------------------------------------------------------'
            write(121,*)'surf all emis             =',SurfGasEmiss_flx_col(idg,NY,NX)          
            write(121,*)'hydloss                   =',GasHydroLoss_flx_col(idg,NY,NX)            
            write(121,*)'RGasNetProd               =',RGasNetProd_col(idg,NY,NX),RGasNetProdSoil_col(idg,NY,NX),&
              RGasNetProd_col(idg,NY,NX)-trcs_Soil2plant_uptake_col(idg,NY,NX)   
            write(121,*)'soil2roots                =',trcs_Soil2plant_uptake_col(idg,NY,NX)                          
            write(121,*)'drib                      =',trcs_solml_drib_col(idg)    
            write(121,*)'tflx.                     =',SurfGasEmiss_flx_col(idg,NY,NX)+GasHydroLoss_flx_col(idg,NY,NX) &
              +RGasNetProd_col(idg,NY,NX)-trcs_Soil2plant_uptake_col(idg,NY,NX)+trcs_solml_drib_col(idg),&
                SurfGasEmiss_flx_col(idg,NY,NX),GasHydroLoss_flx_col(idg,NY,NX), &
                RGasNetProd_col(idg,NY,NX),'drib',trcs_solml_drib_col(idg),trcs_Soil2plant_uptake_col(idg,NY,NX)
            write(121,*)'------------------------------------------------------------------------------'          
            write(121,*)'latloss                   =',trcs_SubsurTransp_flx_2DH(idg,NY,NX)
            write(121,*)'wetdep                    =',Gas_WetDeposit_flx_col(idg,NY,NX)          
            write(121,*)'sol                       =',trcs_solml_vr(idg,NU_col(NY,NX):NL_col(NY,NX),NY,NX)
            write(121,*)'gas                       =',trcg_gasml_vr(idg,NU_col(NY,NX):NL_col(NY,NX),NY,NX)          
            write(121,*)'Hml                       =',trcs_soHml_vr(idg,NU_col(NY,NX):NL_col(NY,NX),NY,NX)  
            if(idg==idg_NH3)then
              write(121,*)'sol_NH3B            =',trcs_solml_vr(idg_NH3B,NU_col(NY,NX):NL_col(NY,NX),NY,NX)
              write(121,*)'Hml_NH3B            =',trcs_soHml_vr(idg_NH3B,NU_col(NY,NX):NL_col(NY,NX),NY,NX)  
            endif
            write(121,*)'transp',NU_col(NY,NX),NL_col(NY,NX)
            write(121,*)transp_diff_slow_vr(idg,NU_col(NY,NX):NL_col(NY,NX),NY,NX)
            write(121,*)'--------------------------'          
            write(121,*)'err in soil               =',AtmGasDiff2Soil_flx_col(idg,NY,NX) &
              +Gas_WetDepo2Soil_col(idg,NY,nX)+Gas_litr2Soil_flx_col(idg,NY,NX)+trcg_AquaADV_Snow2Soil_flx(idg,NY,NX) &
              -trcs_drainage_flx_col(idg,NY,NX)+trcs_SubsurTransp_flx_2DH(idg,NY,NX)+trcg_ebu_flx_col(idg,NY,NX) &
              + trcs_irrig_flx_col(idg,NY,NX)+RGasNetProdSoil_col(idg,NY,NX) &
              - mass_soil+trcs_mass_soil_beg(idg,NY,NX) 
            write(121,*)'diff2soil                 =',AtmGasDiff2Soil_flx_col(idg,NY,NX)
            write(121,*)'wetdep2soil               =',Gas_WetDepo2Soil_col(idg,NY,nX)
            write(121,*)'litr2soil                 =',Gas_litr2Soil_flx_col(idg,NY,NX)
            write(121,*)'sno2soil                  =',trcg_AquaADV_Snow2Soil_flx(idg,NY,NX) 
            write(121,*)'soil2drainage             =',trcs_drainage_flx_col(idg,NY,NX)
            write(121,*)'lateral hydro loss        =',trcs_SubsurTransp_flx_2DH(idg,NY,NX)
            write(121,*)'ebu                       =',trcg_ebu_flx_col(idg,NY,NX)
            write(121,*)'irrig.                    =',trcs_irrig_flx_col(idg,NY,NX)
            write(121,*)'netpro(incl uptk2plt)     =',RGasNetProdSoil_col(idg,NY,NX)
            write(121,*)'drib                      =',trcs_solml_drib_soil_col(idg) 
            write(121,*)'netflx soil               =',AtmGasDiff2Soil_flx_col(idg,NY,NX) &
              +Gas_WetDepo2Soil_col(idg,NY,NX)+Gas_litr2Soil_flx_col(idg,NY,NX)+trcg_AquaADV_Snow2Soil_flx(idg,NY,NX) &
              -trcs_drainage_flx_col(idg,NY,NX)+trcs_SubsurTransp_flx_2DH(idg,NY,NX)+trcg_ebu_flx_col(idg,NY,NX) &
              + trcs_irrig_flx_col(idg,NY,NX)+RGasNetProdSoil_col(idg,NY,NX)
            write(121,*)'------------------------------------------------------------------------------------------------------'
            write(193,*)'tptgas',trcg_gasml_vr(idg,1:NL_col(NY,NX),NY,NX)
            write(193,*)'tptmic',trcs_solml_vr(idg,0:NL_col(NY,NX),NY,NX)
            write(193,*)'tptmac',trcs_soHml_vr(idg,1:NL_col(NY,NX),NY,NX)
            write(194,*)'begmas',trcs_mass_beg(idg,NY,NX)  
            write(194,*)'tptmas',sum(trcg_gasml_vr(idg,1:NL_col(NY,NX),NY,NX))+sum(trcs_solml_vr(idg,0:NL_col(NY,NX),NY,NX))+sum(trcs_soHml_vr(idg,1:NL_col(NY,NX),NY,NX))
          endif
          if(abs(errmass)>1.e-4 .and. abs(safe_adb(errmass,delta_mass))>1.e-3_r8) &
            call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
      ENDDO
    ENDDO
  ENDDO  
  end subroutine ExitMassCheck

!------------------------------------------------------------------------------------------

  SUBROUTINE TranspNoSalt(I,J,NHW,NHE,NVN,NVS)
!
!     THIS SUBROUTINE CALCULATES 3-DIMENSIONAL FLUXES OF ALL SOIL
!     NON-SALT SOLUTES AND GASES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='TranspNoSalt'
  integer :: MX,MM
  integer :: M,NY,NX

!     execution begins here
  call PrintInfo('beg '//subname)

  call EnterMassCheck(NHW,NHE,NVN,NVS)
!
  call InitTranspNoSaltModel(I,J,NHW,NHE,NVN,NVS)
!
! TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
! NPH=NPX,NPT=NPY
! NPH=no. of cycles per hour for water, heat and solute flux calculations
! NPG=NPG=NPH*NPT,number of cycles per hour for gas flux calculations, 
! dt_GasCyc=1.0_r8/NPT, number of cycles for gas flux calculations for 1 iteration in NPH,

  VFLWX=1._r8/NPH
  DO M=1,NPH
    
    !transport without reaction rates, update trcs_solml2_vr
    
    call TransptSlowNoSaltM(I,J,M,NHE,NHW,NVS,NVN)

    !DO gas dissolution-volatization, and solute uptake in each iteration MM        
    call PassSlow2FastIterationM(I,J,M,NHE,NHW,NVS,NVN)
    
    DO MM=1,NPT      
      call TransptFastNoSaltMM(I,J,M,MM,NHE,NHW,NVS,NVN)
    ENDDO
    
    !Do bubbling
    if(ldo_transpt_bubbling)call BubbleEffluxM(I,J,M,NHE,NHW,NVS,NVN)

  ENDDO

  call BackCopyStateVars(I,J,NHE,NHW,NVS,NVN)

  call SummaryFluxes(I,J,NHE,NHW,NVS,NVN)

  call ExitMassCheck(I,J,NHW,NHE,NVN,NVS)

  call PrintInfo('end '//subname)
  
  END subroutine TranspNoSalt

!------------------------------------------------------------------------------------------
  subroutine SummaryFluxes(I,J,NHE,NHW,NVS,NVN)
  implicit none
  integer, intent(in) :: I,J,NHE,NHW,NVS,NVN
  integer :: NY,NX,K
  integer :: n, nitrogen_ids(6), phosphorus_ids(4)
  real(r8) :: surface_n, surface_p

  nitrogen_ids=[ids_NH4,ids_NO3,ids_NO2,ids_NH4B,ids_NO3B,ids_NO2B]
  phosphorus_ids=[ids_H1PO4,ids_H2PO4,ids_H1PO4B,ids_H2PO4B]
  
  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      ! Accepted transport ledgers are signed gains; these inorganic
      ! diagnostics are signed exports (positive out, negative in).
      ! Preserve the legacy carbon grouping: CO2 + CH4, including lateral
      ! and bottom gas-phase transport in the subsurface total.
      HydroSufDICFlx_col(NY,NX)=HydroSufDICFlx_col(NY,NX) &
        -GasHydroSurfLoss_flx_col(idg_CO2,NY,NX)-GasHydroSurfLoss_flx_col(idg_CH4,NY,NX)
      HydroSubsDICFlx_col(NY,NX)=HydroSubsDICFlx_col(NY,NX) &
        -GasHydroSubsLoss_flx_col(idg_CO2,NY,NX)-GasHydroSubsLoss_flx_col(idg_CH4,NY,NX)

      ! DIN excludes N2/N2O. The NH3 subsurface ledger already includes
      ! NH3B and fast gas transport; do not add those contributions twice.
      surface_n=-GasHydroSurfLoss_flx_col(idg_NH3,NY,NX)
      HydroSubsDINFlx_col(NY,NX)=HydroSubsDINFlx_col(NY,NX)-GasHydroSubsLoss_flx_col(idg_NH3,NY,NX)
      do n=1,size(nitrogen_ids)
        HydroSubsDINFlx_col(NY,NX)=HydroSubsDINFlx_col(NY,NX)-trcnuts_hydrloss_flx_col(nitrogen_ids(n),NY,NX)
        if(n<=3)surface_n=surface_n-trcn_SurfRunoff_flx_col(nitrogen_ids(n),NY,NX) &
          -trcn_SnowDrift_flx_col(nitrogen_ids(n),NY,NX)
      enddo
      surface_p=0._r8
      do n=1,size(phosphorus_ids)
        HydroSubsDIPFlx_col(NY,NX)=HydroSubsDIPFlx_col(NY,NX)-trcnuts_hydrloss_flx_col(phosphorus_ids(n),NY,NX)
        if(n<=2)surface_p=surface_p-trcn_SurfRunoff_flx_col(phosphorus_ids(n),NY,NX) &
          -trcn_SnowDrift_flx_col(phosphorus_ids(n),NY,NX)
      enddo
      ! Add this hour only; these two fields retain their annual resets.
      HydroSufDINFlx_CumYr_col(NY,NX)=HydroSufDINFlx_CumYr_col(NY,NX)+surface_n
      HydroSufDIPFlx_CumYr_col(NY,NX)=HydroSufDIPFlx_CumYr_col(NY,NX)+surface_p
      DO K=1,jcplx
        HydroSufDONFlx_col(NY,NX) = HydroSufDONFlx_col(NY,NX)+DOM_SurfRunoff_flx_col(idom_don,K,NY,NX)
        HydroSufDOPFlx_col(NY,NX) = HydroSufDOPFlx_col(NY,NX)+DOM_SurfRunoff_flx_col(idom_dop,K,NY,NX)
        HydroSufDOCFlx_col(NY,NX) = HydroSufDOCFlx_col(NY,NX)+DOM_SurfRunoff_flx_col(idom_doc,K,NY,NX) &
          +DOM_SurfRunoff_flx_col(idom_acetate,K,NY,NX) 

        HydroSubsDONFlx_col(NY,NX)=HydroSubsDONFlx_col(NY,NX)+DOM_Hydroloss_slow_flx_col(idom_don,K,NY,NX)              
        HydroSubsDOPFlx_col(NY,NX)=HydroSubsDOPFlx_col(NY,NX)+DOM_Hydroloss_slow_flx_col(idom_dop,K,NY,NX)              
        HydroSubsDOCFlx_col(NY,NX)=HydroSubsDOCFlx_col(NY,NX)+DOM_Hydroloss_slow_flx_col(idom_doc,K,NY,NX) &
          +DOM_Hydroloss_slow_flx_col(idom_acetate,K,NY,NX)           

      ENDDO
      HydroSubsDOCFlx_col(NY,NX)=HydroSubsDOCFlx_col(NY,NX)-HydroSufDOCFlx_col(NY,NX)
      HydroSubsDONFlx_col(NY,NX)=HydroSubsDONFlx_col(NY,NX)-HydroSufDONFlx_col(NY,NX)
      HydroSubsDOPFlx_col(NY,NX)=HydroSubsDOPFlx_col(NY,NX)-HydroSufDOPFlx_col(NY,NX)
    ENDDO
  ENDDO
  end subroutine SummaryFluxes
!------------------------------------------------------------------------------------------
  subroutine PassSlow2FastIterationM(I,J,M,NHE,NHW,NVS,NVN)
  ! Stage slow-substep states and species coefficients for fast gas subcycles.
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: M
  integer, intent(in) :: NHE,NHW,NVS,NVN
  character(len=*), parameter :: subname='PassSlow2FastIterationM'
  integer :: NY,NX,L,N

  call PrintInfo('beg '//subname)
  WaterFlow2SoilMM_3D=(WaterFlow2MicPM_3D(M,:,:,:,:)+WaterFlow2MacPM_3D(M,:,:,:,:))*dt_GasCyc

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
      !
      ! GASEOUS BOUNDARY LAYER CONDUCTANCES
      !
      ! PARG=boundary layer conductance above soil surface from watsub.f
      !
      ! Species dependence was assembled in the molecular film/pore segments;
      ! aerodynamic transport stays common to gases. Apply only the subcycle time.
      PARGas_CefMM(:,NY,NX) = CondGasXSnowM_col(:,M,NY,NX)*dt_GasCyc

      VLWatMicPXA_vr(0,NY,NX)          = natomw*VLWatMicPM_vr(M,0,NY,NX)

      L=NU_col(NY,NX)
      tScalReductVLsoiAirPMM_vr(L,NY,NX) = ReductVLsoiAirPM_vr(M,L,NY,NX)*dt_GasCyc

      VLWatMicPMA_vr(L,NY,NX)          = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
      VLWatMicPMB_vr(L,NY,NX)          = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      VLWatMicPXA_vr(L,NY,NX)          = natomw*VLWatMicPMA_vr(L,NY,NX)
      VLWatMicPXB_vr(L,NY,NX)          = natomw*VLWatMicPMB_vr(L,NY,NX)
      VLsoiAirPMA_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
      VLsoiAirPMB_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)


      DO L=NU_col(NY,NX)+1,NL_col(NY,NX)
        tScalReductVLsoiAirPMM_vr(L,NY,NX) = ReductVLsoiAirPM_vr(M,L,NY,NX)*dt_GasCyc

        VLWatMicPMA_vr(L,NY,NX) = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
        VLWatMicPMB_vr(L,NY,NX) = VLWatMicPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)
        VLWatMicPXA_vr(L,NY,NX) = natomw*VLWatMicPMA_vr(L,NY,NX)
        VLWatMicPXB_vr(L,NY,NX) = natomw*VLWatMicPMB_vr(L,NY,NX)

        VLsoiAirPMA_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4,L,NY,NX)
        VLsoiAirPMB_vr(L,NY,NX)          = VLsoiAirPM_vr(M,L,NY,NX)*trcs_VLN_vr(ids_NH4B,L,NY,NX)

      ENDDO
    ENDDO
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine PassSlow2FastIterationM

end module TranspNoSaltMod
