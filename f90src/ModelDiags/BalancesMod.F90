module BalancesMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use CanopyDataType,    only: QVegET_col
  use GridDataType,      only: NU_col, NL_col
  use EcoSimConst,       only: DENSICE
  use abortutils,        only: endrun
  use EcoSIMCtrlMod,     only: fixWaterLevel,etimer
  use AqueChemDatatype,  only: trcg_mass_cumerr_col
  use PlantMgmtDataType, only: iDayPlanting_pft, iDayPlantHarvest_pft
  use EcoSIMCtrlMod,     only: iVerbLevel
  use NumericalAuxMod
  use RootDataType
  use SoilBGCDataType
  use GridDataType
  use SurfLitterDataType
  use CanopyDataType
  use SnowDataType
  use BalanceCheckDataType
  use SoilWaterDataType
  use ClimForcDataType
  use SurfSoilDataType
  use EcosimBGCFluxType
  use SoilHeatDataType
  use PlantDataRateType
  use EcoSimSumDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: BegCheckBalances
  public :: EndCheckBalances
  public :: SummarizeTracerMass
  public :: SummarizeSnowMass
  public :: SummarizeTracers
contains

  subroutine BegCheckBalances(I,J,NHW,NHE,NVN,NVS)  
  !
  !Description
  !Prepare for next mass balance check
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX,idg,ids,L
  
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      WaterErr_col(NY,NX)           = WatMass_col(NY,NX)
      TillageLitr2SoilH2O_col(NY,NX) = 0._r8
      HeatErr_col(NY,NX)            = HeatStore_col(NY,NX)
      SnowEngyBeg_col(NY,NX)        = SnowEngyEnd_col(NY,NX)
      CanopyWaterMassBeg_col(NY,NX) = CanopyWaterMassEnd_col(NY,NX)
      SnowMassBeg_col(NY,NX)        = SnowMassEnd_col(NY,NX)
      LitWatMassBeg_col(NY,NX)      = LitWatMassEnd_col(NY,NX)
      SoilWatMassBeg_col(NY,NX)     = SoilWatMassEnd_col(NY,NX)

      DO ids=ids_beg,ids_end
        trcs_solml_dribBeg_col(ids,NY,NX)=trcs_solml_drib_vr(ids,0,NY,NX)
      ENDDO

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO ids=ids_beg,ids_end
          trcs_solml_dribBeg_col(ids,NY,NX)=trcs_solml_dribBeg_col(ids,NY,NX)+trcs_solml_drib_vr(ids,L,NY,NX)
        ENDDO        
      ENDDO

      DO idg=idg_beg,idg_end
        trcg_TotalMass_beg_col(idg,NY,NX) = trcg_TotalMass_col(idg,NY,NX)
        trcg_soilMass_beg_col(idg,NY,NX) = trcg_soilMass_col(idg,NY,NX)
      enddo
      DO idg=idg_beg,idg_NH3      
        trcg_rootMass_beg_col(idg,NY,NX)  = trcg_rootMass_col(idg,NY,NX)  
        trcg_snowMass_beg_col(idg,NY,NX)  = trcg_snowMass_col(idg,NY,NX)
      enddo  
    ENDDO
  ENDDO  
  end subroutine BegCheckBalances
!------------------------------------------------------------------------------------------
  subroutine SummarizeStorage(I,J,NHW,NHE,NVN,NVS)  

  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX

  WatMassStore_lnd = 0._r8
  HeatStore_lnd    = 0._r8
  DO NX = NHW, NHE
    DO  NY=NVN,NVS

      call SumUpWaterStorage(NY,NX)

      call SumUpHeatStorage(NY,NX)

    ENDDO
  ENDDO
  end subroutine SummarizeStorage
!------------------------------------------------------------------------------------------

  subroutine SummarizeSnowMass(NY,NX)
  implicit none
  integer, intent(in) :: NY,NX
  integer :: L
  SnowMassEnd_col(NY,NX)=0._r8
  DO  L=1,JS
    if(AMAX1(VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX))<=0._r8)cycle
    SnowMassEnd_col(NY,NX) = SnowMassEnd_col(NY,NX)+VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
  ENDDO

  end subroutine SummarizeSnowMass
!------------------------------------------------------------------------------------------
  subroutine SumUpWaterStorage(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  !canopy
  CanopyWaterMassEnd_col(NY,NX) = WatHeldOnCanopy_col(NY,NX)+CanopyBiomWater_col(NY,NX)+SnowOnCanopy_col(NY,NX)
  !snow  
  SnowMassEnd_col(NY,NX)=0._r8
  DO  L=1,JS
    if(AMAX1(VLDrySnoWE_snvr(L,NY,NX),VLWatSnow_snvr(L,NY,NX),VLIceSnow_snvr(L,NY,NX))<=0._r8)cycle
    SnowMassEnd_col(NY,NX) = SnowMassEnd_col(NY,NX)+VLDrySnoWE_snvr(L,NY,NX)+VLWatSnow_snvr(L,NY,NX)+VLIceSnow_snvr(L,NY,NX)*DENSICE
  ENDDO
  !litter
  LitWatMassEnd_col(NY,NX) = VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)*DENSICE

  WatMass_col(NY,NX)       = CanopyWaterMassEnd_col(NY,NX)+SnowMassEnd_col(NY,NX)+LitWatMassEnd_col(NY,NX)
  SoilWatMassEnd_col(NY,NX)= 0._r8
  !add water in soil 
  DO L=NU_col(NY,NX),NL_col(NY,NX)
    SoilWatMassEnd_col(NY,NX) = SoilWatMassEnd_col(NY,NX)+ VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)+(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))*DENSICE
  ENDDO
!  write(444,*)'bal',('=x',L=1,50)
!  write(444,*)(VLWatMicP_vr(L,NY,NX),L=NU_col(NY,NX),NL_col(NY,NX))
!  write(444,*)(VLWatMacP_vr(L,NY,NX),L=NU_col(NY,NX),NL_col(NY,NX))
!  write(444,*)((VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))*DENSICE,L=NU_col(NY,NX),NL_col(NY,NX))

  WatMass_col(NY,NX)=WatMass_col(NY,NX)+SoilWatMassEnd_col(NY,NX)

  WatMassStore_lnd = WatMassStore_lnd+WatMass_col(NY,NX)

  end subroutine SumUpWaterStorage

!------------------------------------------------------------------------------------------
  subroutine SumUpHeatStorage(NY,NX)

  implicit none
  integer, intent(in) :: NY,NX
  integer :: L

  HeatStore_col(NY,NX) = 0._r8
  !add snow
  DO  L=1,JS
    if(VLHeatCapSnow_snvr(L,NY,NX)<=0._r8)cycle
    HeatStore_col(NY,NX) = HeatStore_col(NY,NX)+ VLHeatCapSnow_snvr(L,NY,NX)*TKSnow_snvr(L,NY,NX)    
  ENDDO
  SnowEngyEnd_col(NY,NX) = HeatStore_col(NY,NX)
  !add canopy + litter
  HeatStore_col(NY,NX)   = HeatStore_col(NY,NX)+CanopyHeatStor_col(NY,NX) + TKS_vr(0,NY,NX)*VHeatCapacity_vr(0,NY,NX)

  !add soil
  DO L=NU_col(NY,NX),NL_col(NY,NX)
    HeatStore_col(NY,NX) = HeatStore_col(NY,NX)+ VHeatCapacity_vr(L,NY,NX)*TKS_vr(L,NY,NX)
  ENDDO

  HeatStore_lnd = HeatStore_lnd+HeatStore_col(NY,NX)

  end subroutine SumUpHeatStorage
!------------------------------------------------------------------------------------------

  subroutine EndCheckBalances(I,J,NHW,NHE,NVN,NVS)  
  !
  !Description
  !Do mass balance check
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS
  integer :: NY,NX
  real(r8), parameter :: err_h2o=1.e-4_r8
  real(r8), parameter :: err_engy=1.e-6_r8
  real(r8) :: WaterErr_test
  real(r8) :: HeatErr_test
  real(r8) :: literH2Oerr_test
  real(r8) :: canopyH2Oerr_test
  real(r8) :: SnowMassErr_test
  real(r8) :: SoilWatErr_test
  real(r8) :: precipErr_test
  real(r8) :: prec2expSErr_test
  real(r8) :: prec2SnoErr_test
  real(r8) :: tracer_mass_err
  real(r8) :: tracer_rootmass_err
  real(r8) :: tracer_snowmass_err
  real(r8) :: dCO2err       !Inconsistency error due to asynchrnous update between soil and wihtin root CO2 [gC d-2]
  real(r8) :: dgaserr       !Inconsistency error due to asynchrnous update between soil and wihtin root gas [g d-2] 
  integer  :: ii,idg,ids,L
  real(r8) :: trcs_solml_drib_col(ids_beg:ids_end)

  call SummarizeTracerMass(I,J,NHW,NHE,NVN,NVS)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS

      if(fixWaterLevel)then
        SoilWatErr_test=SoilWatMassBeg_col(NY,NX)-SoilWatMassEnd_col(NY,NX)      
      else
        SoilWatErr_test=SoilWatMassBeg_col(NY,NX)-SoilWatMassEnd_col(NY,NX)+Qinflx2Soil_col(NY,NX) &
          -QDrain_col(NY,NX)-QDischarg2WTBL_col(NY,NX)+TPlantRootH2OUptake_col(NY,NX)+QLaterFlow2Cell_col(NY,NX) &
          +TillageLitr2SoilH2O_col(NY,NX)
      endif
      precipErr_test    = RainPrecThrufall_col(NY,NX)-Rain2LitR_col(NY,NX)-Rain2Soil_col(NY,NX)-RainPrec2Sno_col(NY,NX)
      prec2expSErr_test = Rain2ExposedSurf_col(NY,NX)-Rain2LitR_col(NY,NX)-Rain2Soil_col(NY,NX)
      prec2SnoErr_test  = Prec2Snow_col(NY,NX)-RainPrec2Sno_col(NY,NX)-SnoFalPrec_col(NY,NX)

      literH2Oerr_test = LitWatMassBeg_col(NY,NX)-LitWatMassEnd_col(NY,NX)+QRunSurf_col(NY,NX)+WatFLo2LitR_col(NY,NX) &
        -TillageLitr2SoilH2O_col(NY,NX)

      SnowMassErr_test = SnowMassBeg_col(NY,NX)-SnowMassEnd_col(NY,NX)+Prec2Snow_col(NY,NX)-QSnowH2Oloss_col(NY,NX)

      canopyH2Oerr_test=CanopyWaterMassBeg_col(NY,NX)-CanopyWaterMassEnd_col(NY,NX)+RainIntceptByCanopy_col(NY,NX) &
        +QVegET_col(NY,NX)-TPlantRootH2OUptake_col(NY,NX)-QCanopyWatLoss2Dist_col(NY,NX)+SnowIntceptByCanopy_col(NY,NX)

      if(fixWaterLevel)then
        WaterErr_test = WaterErr_col(NY,NX)-WatMass_col(NY,NX)
      else      
        WaterErr_test = WaterErr_col(NY,NX)-WatMass_col(NY,NX)+PrecAtm_col(NY,NX)+Irrigation_col(NY,NX)+QLaterFlow2Cell_col(NY,NX) &
          +RainLitr_col(NY,NX)+VapXAir2GSurf_col(NY,NX)+QVegET_col(NY,NX)+QRunSurf_col(NY,NX) &
          -QDrain_col(NY,NX)-QDischarg2WTBL_col(NY,NX)+TPlantRootH2OUptake_col(NY,NX)-QCanopyWatLoss2Dist_col(NY,NX)
      endif
      HeatErr_test = HeatErr_col(NY,NX)-HeatStore_col(NY,NX)+THeatRootRelease_col(NY,NX) &
        +HeatSource_col(NY,NX)+Eco_NetRad_col(NY,NX)+Eco_Heat_Latent_col(NY,NX)+Eco_Heat_Sens_col(NY,NX)&
        +PrecHeat_col(NY,NX)+THeatSoiThaw_col(NY,NX)+THeatSnowThaw_col(NY,NX)+HeatRunSurf_col(NY,NX) &
        -HeatDrain_col(NY,NX)-HeatDischar_col(NY,NX)-CanopyHeatLoss2Dist_col(NY,NX)

      if(abs(WaterErr_test)>err_h2o .and. etimer%get_nstep()>1)then
        if(iVerbLevel==1 .or. abs(SoilWatErr_test)>err_h2o)then
          write(110,*)('=',ii=1,50)
          write(110,*)I*1000+J/24.,'NY,NX ',NY,NX
          write(110,*)'init H2O         =',WaterErr_col(NY,NX)
          write(110,*)'final H2O        =',WatMass_col(NY,NX)
          write(110,*)'delta H2O        =',WatMass_col(NY,NX)-WaterErr_col(NY,NX)
          write(110,*)'prec, irrig      =',PrecAtm_col(NY,NX),Irrigation_col(NY,NX) 
          write(110,*)'litfall H2O      =',RainLitr_col(NY,NX)
          write(110,*)'surf evap        =',VapXAir2GSurf_col(NY,NX)    
          write(110,*)'RainIntceptByCan =',RainIntceptByCanopy_col(NY,NX)   
          write(110,*)'SnowIntceptByCan =',SnowIntceptByCanopy_col(NY,NX)
          write(110,*)'plant evap       =',QVegET_col(NY,NX)
          write(110,*)'run on           =',QRunSurf_col(NY,NX),WatFLo2LitR_col(NY,NX)
          write(110,*)('-',ii=1,50)
          write(110,*)'Qinflx2Soil_col  =',Qinflx2Soil_col(NY,NX)
          write(110,*)'drain            =',QDrain_col(NY,NX)
          write(110,*)'discharge        =',QDischarg2WTBL_col(NY,NX)
          write(110,*)'root uptake      =',TPlantRootH2OUptake_col(NY,NX)    
          write(110,*)'QWatIntLaterflow =',QLaterFlow2Cell_col(NY,NX)
          write(110,*)('-',ii=1,50)
          write(110,*)'SnowMassBeg,end  =',SnowMassBeg_col(NY,NX),SnowMassEnd_col(NY,NX)
          write(110,*)'CanopyWatbeg,end =',CanopyWaterMassBeg_col(NY,NX),CanopyWaterMassEnd_col(NY,NX)          
          write(110,*)'WatHeldOnCanopy  =',WatHeldOnCanopy_col(NY,NX)
          write(110,*)'SnowOnCanopy     =',SnowOnCanopy_col(NY,NX)
          write(110,*)'CanopyBiomWater_col    =',CanopyBiomWater_col(NY,NX)
          write(110,*)'SoilWatMassBegEnd=',SoilWatMassBeg_col(NY,NX),SoilWatMassEnd_col(NY,NX),&
            SoilWatMassBeg_col(NY,NX)-SoilWatMassEnd_col(NY,NX),NU_col(NY,NX)
          write(110,*)'snofall          =',Prec2Snow_col(NY,NX),SnoFalPrec_col(NY,NX)
          write(110,*)'nsnol_col        =',nsnol_col(NY,NX)
          write(110,*)'snowloss         =',QSnowH2Oloss_col(NY,NX)
          write(110,*)'Snowxfer         =',QSnoWatXfer2Soil_col(NY,NX)+QSnoIceXfer2Soil_col(NY,NX)*DENSICE
          write(110,*)('-',ii=1,50)
          write(110,*)'precipErr_test   =',precipErr_test,RainPrecThrufall_col(NY,NX),RainPrec2Sno_col(NY,NX)
          write(110,*)'prec2expSErr_test=',prec2expSErr_test,Rain2ExposedSurf_col(NY,NX),Rain2LitR_col(NY,NX),Rain2Soil_col(NY,NX)
          write(110,*)'prec2SnoErr_test =',prec2SnoErr_test        
          write(110,*)'canopyH2Oerr_test=',canopyH2Oerr_test
          write(110,*)'SnowMassErr_test =',SnowMassErr_test        
          write(110,*)'literH2Oerr_test =',literH2Oerr_test,LitWatMassBeg_col(NY,NX),LitWatMassEnd_col(NY,NX)
          write(110,*)'SoilWatErr_test  =',SoilWatErr_test
          write(110,*)'col WaterErr_test=',WaterErr_test
        endif
        if(abs(SoilWatErr_test)>err_h2o ) & ! .or. abs(SnowMassErr_test)>1.e-3_r8*SnowMassBeg_col(NY,NX)) &
          call endrun('H2O error test failure in '//trim(mod_filename)//' at line',__LINE__)
        
      endif

      DO ids=ids_beg,ids_end
        trcs_solml_drib_col(ids)=trcs_solml_drib_vr(ids,0,NY,NX)
      ENDDO

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO ids=ids_beg,ids_end
          trcs_solml_drib_col(ids)=trcs_solml_drib_col(ids)+trcs_solml_drib_vr(ids,L,NY,NX)
        ENDDO        
      ENDDO

      !Turn off the tracer mass conservation check momentarily, due to complication of 
      !grid change. It will be turned on when a safe strategy will be figured out later.
      
      Hydroloss_NH4_cumflx_col(NY,NX)=Hydroloss_NH4_cumflx_col(NY,NX)+ &
        trcnuts_hydrloss_flx_col(ids_NH4,NY,NX)+trcnuts_hydrloss_flx_col(ids_NH4B,NY,NX)
      Hydroloss_NO3_cumflx_col(NY,NX)=Hydroloss_NO3_cumflx_col(NY,NX)+ &
        trcnuts_hydrloss_flx_col(ids_NO3,NY,NX)+trcnuts_hydrloss_flx_col(ids_NO3B,NY,NX)

      DO idg=idg_beg,idg_NH3        
        
        tracer_mass_err = trcg_TotalMass_beg_col(idg,NY,NX) -trcg_TotalMass_col(idg,NY,NX) &
          +SurfGasEmiss_all_flx_col(idg,NY,NX)+GasHydroLoss_flx_col(idg,NY,NX) &
          +trcs_solml_drib_col(idg)+RGasNetProd_col(idg,NY,NX)-trcs_deadroot2soil_col(idg,NY,NX) 

        GasHydroLoss_cumflx_col(idg,NY,NX)=GasHydroLoss_cumflx_col(idg,NY,NX)+GasHydroLoss_flx_col(idg,NY,NX)

        tracer_rootmass_err = trcg_rootMass_beg_col(idg,NY,NX)-trcg_rootMass_col(idg,NY,NX)-trcs_deadroot2soil_col(idg,NY,NX) &
          +trcs_Soil2plant_uptake_col(idg,NY,NX)+trcg_air2root_flx_col(idg,NY,NX)+TRootGasLossDisturb_col(idg,NY,NX)

        if(idg==idg_NH3)then           
          tracer_mass_err = tracer_mass_err+trcg_TotalMass_beg_col(idg_NH3B,NY,NX)-trcg_TotalMass_col(idg_NH3B,NY,NX) 
                       
          tracer_mass_err=tracer_mass_err+trcs_solml_drib_col(idg_NH3B) 

          tracer_rootmass_err=tracer_rootmass_err+trcs_Soil2plant_uptake_col(idg_NH3B,NY,NX)
        endif

        if(idg.eq.idg_O2)then
          tracer_rootmass_err = tracer_rootmass_err-RootO2_TotSink_col(NY,NX)
        elseif(idg.eq.idg_CO2)then
          tracer_rootmass_err = tracer_rootmass_err+RootCO2Ar2Root_col(NY,NX)
        endif

        trcg_mass_cumerr_col(idg,NY,NX)=trcg_mass_cumerr_col(idg,NY,NX)+ tracer_mass_err         
        tracer_snowmass_err=trcg_snowMass_beg_col(idg,NY,NX)-trcg_snowMass_col(idg,NY,NX) + &
          Gas_WetDepo2Snow_col(idg,NY,NX)-Gas_Snowloss_flx_col(idg,NY,NX)          

        if(AMAX1(abs(tracer_mass_err),abs(tracer_rootmass_err))>1.e-5_r8 .and. etimer%get_nstep()>1)then
          if(iVerbLevel==1 .or. abs(tracer_mass_err)>1.e-1_r8)then
            write(111,*)('-',ii=1,50)
            write(111,*)'NU   =',NU_col(NY,NX)
            write(111,*)I*1000+J,'NY NX=',NY,NX,trcs_names(idg),iDayPlantHarvest_pft(1,NY,NX),iDayPlanting_pft(1,NY,NX),'Final'
            write(111,*)'beg end delta mass=',trcg_TotalMass_beg_col(idg,NY,NX),trcg_TotalMass_col(idg,NY,NX),&
              trcg_TotalMass_beg_col(idg,NY,NX)-trcg_TotalMass_col(idg,NY,NX)
            write(111,*)'mass_err         =',tracer_mass_err
            write(111,*)'gasdif,ebu       =',GasDiff2Surf_flx_col(idg,NY,NX),trcg_ebu_flx_col(idg,NY,NX)
            write(111,*)'phenoflx         =',TRootGasLossDisturb_col(idg,NY,NX)
            write(111,*)'wedepo           =',Gas_WetDeposit_flx_col(idg,NY,NX)          
            write(111,*)'----------------------'
            write(111,*)'surf all emis    =',SurfGasEmiss_all_flx_col(idg,NY,NX)          
            write(111,*)'surf emis        =',SurfGasEmiss_flx_col(idg,NY,NX), SurfGasEmiss_all_flx_col(idg,NY,NX) -SurfGasEmiss_flx_col(idg,NY,NX)                      
            write(111,*)'GasHydroloss     =',GasHydroLoss_flx_col(idg,NY,NX)
            if(idg.eq.idg_N2)then
              write(111,*)'RGasNetProd,rNfix,mup=',RGasNetProd_col(idg,NY,NX),RootN2Fix_col(NY,NX),trcs_RMicbUptake_col(idg,NY,NX)
            elseif(idg.eq.idg_CO2)then
              write(111,*)'RGasNetProd     =',RGasNetProd_col(idg,NY,NX),-RootCO2AutorPrev_col(NY,NX)  
            elseif(idg.eq.idg_NH3)then
              write(111,*)'RGasNetProd,micup =',RGasNetProd_col(idg,NY,NX)
            else
              write(111,*)'RGasNetProd,micup =',RGasNetProd_col(idg,NY,NX),trcs_RMicbUptake_col(idg,NY,NX)
            endif
            write(111,*)'drib                =',trcs_solml_drib_col(idg) 
            write(111,*)'------------------'
            write(111,*)'snow beg_end mass=',trcg_snowMass_beg_col(idg,NY,NX),trcg_snowMass_col(idg,NY,NX)
            write(111,*)'snowloss2SoilLitr=',trcg_snowMassloss_col(idg,NY,NX)
            write(111,*)'tracersnowfall   =',trcg_AquaAdv_flx_snvr(idg,1,NY,NX)
            write(111,*)'snowmass err     =',tracer_snowmass_err
            write(111,*)'------------------'
            write(111,*)'xroot beg_end_delta mass     =',trcg_TotalMass_beg_col(idg,NY,NX)-trcg_rootMass_beg_col(idg,NY,NX),&
              trcg_TotalMass_col(idg,NY,NX)-trcg_rootMass_col(idg,NY,NX),trcg_TotalMass_beg_col(idg,NY,NX)-trcg_rootMass_beg_col(idg,NY,NX) &
              -(trcg_TotalMass_col(idg,NY,NX)-trcg_rootMass_col(idg,NY,NX))  
!            write(194,*)'balmas',trcg_TotalMass_beg_col(idg,NY,NX)-trcg_rootMass_beg_col(idg,NY,NX),trcg_TotalMass_col(idg,NY,NX)-trcg_rootMass_col(idg,NY,NX)  
            write(111,*)'soil beg_end mass delta dmass=',trcg_soilMass_beg_col(idg,NY,NX), trcg_soilMass_col(idg,NY,NX),&
              trcg_soilMass_beg_col(idg,NY,NX)-trcg_soilMass_col(idg,NY,NX)
            write(111,*)'soil2atm         =',trcg_ebu_flx_col(idg,NY,NX)+GasDiff2Surf_flx_col(idg,NY,NX)+Gas_WetDeposit_flx_col(idg,NY,NX)  
            if(idg==idg_O2)then
              write(111,*)'netpro soil.   =',RGasNetProd_col(idg_O2,NY,NX),RootO2_TotSink_col(NY,NX)
              write(111,*)'tflx.          =',SurfGasEmiss_flx_col(idg,NY,NX)+GasHydroLoss_flx_col(idg,NY,NX) &
                +RGasNetProd_col(idg_O2,NY,NX)+RootO2_TotSink_col(NY,NX)+trcs_solml_drib_col(idg)-trcs_Soil2plant_uptake_col(idg,NY,NX),&
                SurfGasEmiss_flx_col(idg,NY,NX),GasHydroLoss_flx_col(idg,NY,NX), &
                RGasNetProd_col(idg_O2,NY,NX),'drib',trcs_solml_drib_col(idg),trcs_Soil2plant_uptake_col(idg,NY,NX),RootO2_TotSink_col(NY,NX)
            endif
            if(idg==idg_CO2)then
              write(111,*)'ar2soil          =',RootCO2Ar2Soil_col(NY,NX)
            endif    
            write(111,*)'col_mass loss      =',trcg_ebu_flx_col(idg,NY,NX)+GasDiff2Surf_flx_col(idg,NY,NX)+Gas_WetDeposit_flx_col(idg,NY,NX) &
                +trcs_Soil2plant_uptake_col(idg,NY,NX) 
            write(111,*)'=================='
            write(111,*)'root beg_end delta mass=',trcg_rootMass_beg_col(idg,NY,NX),trcg_rootMass_col(idg,NY,NX),&
              trcg_rootMass_beg_col(idg,NY,NX)-trcg_rootMass_col(idg,NY,NX)
            write(111,*)'root trcmass_err       =',tracer_rootmass_err
            write(111,*)'pltair2root            =',trcg_air2root_flx_col(idg,NY,NX)     

            if(idg==idg_O2)then
              write(111,*)'total O2Sink      =',RootO2_TotSink_col(NY,NX)
            endif  
            
            if(idg==idg_CO2)then
              dCO2err=RootCO2Emis2Root_col(NY,NX)-trcs_Soil2plant_uptake_col(idg_CO2,NY,NX)-RootCO2Ar2Root_col(NY,NX)
              write(111,*)'tplt2root,soi2root =',RootCO2Emis2Root_col(NY,NX),trcs_Soil2plant_uptake_col(idg_CO2,NY,NX)
              write(111,*)'ar2root, deltaErr  =',RootCO2Ar2Root_col(NY,NX),dCO2err
            else
              dgaserr=trcs_Soil2plant_uptake_col(idg,NY,NX)-trcs_Soil2plant_uptakep_col(idg,NY,NX)
              if(idg==idg_NH3)then
                dgaserr=dgaserr+trcs_Soil2plant_uptake_col(idg_NH3B,NY,NX)-trcs_Soil2plant_uptakep_col(idg_NH3B,NY,NX)
                write(111,*)'plt_soiluptake     =',trcs_Soil2plant_uptake_col(idg,NY,NX)+trcs_Soil2plant_uptake_col(idg_NH3B,NY,NX),&
                  trcs_Soil2plant_uptakep_col(idg,NY,NX)+trcs_Soil2plant_uptakep_col(idg_NH3B,NY,NX),dgaserr
              else
                write(111,*)'plt_soiluptake     =',trcs_Soil2plant_uptake_col(idg,NY,NX),trcs_Soil2plant_uptakep_col(idg,NY,NX),dgaserr
              endif
            endif
            write(111,*)'deadroot2soil    =',trcs_deadroot2soil_col(idg,NY,NX)
          endif

          if(abs(tracer_mass_err)>1.e-1_r8) &
            call endrun('tracer'//trcs_names(idg)//' error test failure in '//trim(mod_filename)//' at line',__LINE__)
        endif
      enddo      

      DO ids=ids_beg,ids_end
        trcs_solml_dribBeg_col(ids,NY,NX)=trcs_solml_drib_col(ids)
      ENDDO        
      cycle      
      if(abs(HeatErr_test)>err_engy)then
        write(110,*)('-',ii=1,50)
        write(110,*)I+J/24.,'NY,NX,Engy conservation error',NY,NX,HeatErr_test
        write(110,*)'Init Engy=',HeatErr_col(NY,NX)
        write(110,*)'Final Engy=',HeatStore_col(NY,NX)
        write(110,*)'Delta Engy=',HeatStore_col(NY,NX)-HeatErr_col(NY,NX)
        write(110,*)'Delta snow engy=',SnowEngyEnd_col(NY,NX)-SnowEngyBeg_col(NY,NX)
        write(110,*)'Heat source=',HeatSource_col(NY,NX)
        write(110,*)'Net Rad=',Eco_NetRad_col(NY,NX)
        write(110,*)'Latent Heat=',Eco_Heat_Latent_col(NY,NX)
        write(110,*)'Sensible Heat=',Eco_Heat_Sens_col(NY,NX)
        write(110,*)'Heat prec=',PrecHeat_col(NY,NX)
        write(110,*)'Heat phase change, soil, snow=',THeatSoiThaw_col(NY,NX),THeatSnowThaw_col(NY,NX)
        write(110,*)'Heat runoff=',HeatRunSurf_col(NY,NX)
        write(110,*)'Heat drain=',HeatDrain_col(NY,NX)
        write(110,*)'Heat disch=',HeatDischar_col(NY,NX)
        call endrun('Engy error test failure in '//trim(mod_filename)//' at line',__LINE__)
      endif

    ENDDO
  ENDDO
  end subroutine EndCheckBalances 


!------------------------------------------------------------------------------------------
  subroutine SummarizeTracers(I,J,NHW,NHE,NVN,NVS)
  !
  !Description
  !sum up mass of tracers
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NHW,NHE,NVN,NVS  
  integer :: NY,NX
  integer :: idg,L
  real(r8) :: trcg_snow(idg_beg:idg_NH3)
  real(r8) :: trcg_root(idg_beg:idg_NH3)
  real(r8) :: trcg_litr(idg_beg:idg_NH3)
  real(r8) :: trcg_soil(idg_beg:idg_end)

  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      trcg_TotalMass_col(idg_beg:idg_end,NY,NX)=0._r8
      trcg_snow=0._r8
      trcg_root=0._r8
      trcg_soil=0._r8
      trcg_litr=0._r8
      DO idg=idg_beg,idg_NH3
        trcg_litr(idg)=trcs_solml_vr(idg,0,NY,NX)
      enddo
      trcg_litr(idg_NH3)=trcg_litr(idg_NH3)+trcs_solml_vr(idg_NH3B,0,NY,NX)
      !note NL may change due to soil relayering, which can be tested using
      !by comparing the values of NL and NLI      

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO idg=idg_beg,idg_NH3
          trcg_soil(idg)=trcg_soil(idg) + trcg_gasml_vr(idg,L,NY,NX)
          trcg_root(idg)=trcg_root(idg) + trcg_root_vr(idg,L,NY,NX)
        ENDDO

        DO idg=idg_beg,idg_end  
          trcg_soil(idg)=trcg_soil(idg) + trcs_solml_vr(idg,L,NY,NX) + trcs_soHml_vr(idg,L,NY,NX)
        ENDDO
      ENDDO

      DO idg=idg_beg,idg_end  
        trcg_soilMass_col(idg,NY,NX)=trcg_soil(idg)
      ENDDO    

      !Because idg_NH3B does not exist in snow
      !add now      
      DO L=1,nsnol_col(NY,NX)          
        DO idg=idg_beg,idg_NH3
          trcg_snow(idg)=trcg_snow(idg)+trcg_solsml_snvr(idg,L,NY,NX)
        ENDDO
      ENDDO    
      DO idg=idg_beg,idg_NH3  
        trcg_snowMass_col(idg,NY,NX)=trcg_snow(idg)
      ENDDO    
      DO idg=idg_beg,idg_NH3
        trcg_rootMass_col(idg,NY,NX) = trcg_root(idg)
      ENDDO
      
      DO idg=idg_beg,idg_NH3
        trcg_TotalMass_col(idg,NY,NX)=trcg_snowMass_col(idg,NY,NX)+trcg_rootMass_col(idg,NY,NX)+trcg_soilMass_col(idg,NY,NX)+trcg_litr(idg)
      ENDDO
      trcg_TotalMass_col(idg_NH3B,NY,NX)=trcg_soilMass_col(idg_NH3B,NY,NX)
      
      if(.false.)then 
        idg=idg_CO2
        write(115,*)I*1000+J,NY,NX,'nsnol_col     trcname     snow   litr  root  soil  total',NU_col(NY,NX),NUI_col(NY,NX)
        write(115,*)nsnol_col(NY,NX),trcs_names(idg),trcg_snow(idg),trcg_litr(idg),trcg_root(idg),trcg_soilMass_col(idg,NY,NX),trcg_TotalMass_col(idg,NY,NX)
        DO L=NUI_col(NY,NX),NL_col(NY,NX)
          write(115,*)L,trcg_gasml_vr(idg,L,NY,NX),trcs_solml_vr(idg,L,NY,NX), trcs_soHml_vr(idg,L,NY,NX)
        ENDDO  
      endif

    ENDDO
  ENDDO
  end subroutine SummarizeTracers

!------------------------------------------------------------------------------------------

  subroutine SummarizeTracerMass(I,J,NHW,NHE,NVN,NVS)

  implicit none 
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS

  call SummarizeStorage(I,J,NHW,NHE,NVN,NVS)

  call SummarizeTracers(I,J,NHW,NHE,NVN,NVS)

  end subroutine SummarizeTracerMass

end module BalancesMod
