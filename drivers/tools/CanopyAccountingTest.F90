program CanopyAccountingTest
  ! Description:
  ! Exercise accepted live exchanges after complete population removal, including
  ! root water/heat, nutrients, NH3, canopy fluxes and retained-slot balance residuals.
  use data_kind_mod, only: r8 => DAT_KIND_R8,yearIJ_type
  use EcoSiMParDataMod, only: micpar,pltpar
  use InitAllocMod, only: InitPlantMorphSize
  use PlantAPIData
  use TracerIDMod
  use ElmIDMod, only: ifalse
  use EcoSimConst, only: cpw
  use EcoSIMCtrlMod, only: ldo_sp_mode
  use ExtractsMod, only: extracts
  use DebugToolMod, only: PrintInfo
  use abortutils, only: iulog
  implicit none
  type(yearIJ_type) :: yearIJ
  logical, allocatable :: LiveExchange(:),PassiveExchange(:)
  call micpar%Init()
  call InitPlantMorphSize()
  ! This isolated API test has no grid namelist to supply dimensions.
  pltpar%JP1=3;pltpar%JZ1=2;pltpar%NumCanopyLayers1=2
  pltpar%NumOfLeafAzimuthSectors=1;pltpar%NumOfSkyAzimuthSects1=1
  pltpar%NumLeafInclinationClasses1=1;pltpar%MaxNodesPerBranch1=2
  pltpar%jroots=2
  call InitTracerIDs(.false.)
  call InitPlantAPIData()
  allocate(LiveExchange(JP1));LiveExchange=.false.;LiveExchange(1)=.true.
  allocate(PassiveExchange(JP1));PassiveExchange=.false.
  yearIJ%I=1;yearIJ%J=12
  ldo_sp_mode=.false.
  plt_site%NP=1;plt_site%NP0=2;plt_site%NU=1;plt_site%NK=1
  plt_site%AREA3=1._r8;plt_site%PlantPopuLive_pft=0._r8
  plt_pheno%IsPlantActive_pft=ifalse
  plt_site%PlantElmBalCum_pft=0._r8
  plt_site%PlantElmBalCum_pft(:,1)=[1._r8,2._r8,3._r8]
  plt_site%PlantElmBalCum_pft(:,2)=[4._r8,5._r8,6._r8]
  plt_site%PlantElemntStoreLandscape=0._r8
  plt_morph%Myco_pft=1;plt_morph%MaxSoilLays4Root_pft=1
  plt_site%MaxNumRootLays=1
  plt_morph%CanopyLeafAreaZ_pft=0._r8;plt_morph%CanopyStemSurfAreaZ_pft=0._r8
  plt_morph%CanopyLeafArea_pft=0._r8;plt_morph%CanopyStemSurfArea_pft=0._r8
  plt_morph%totRootLenDens_vr=0._r8;plt_morph%RootLenDensPerPlant_pvr=0._r8
  plt_biom%CanopyLeafCLyr_pft=0._r8;plt_biom%StandDeadStrutElms_pft=0._r8
  plt_biom%StandDeadCompKElms_pft=0._r8
  plt_biom%StandingDeadStrutElms_col=0._r8
  plt_bgcr%LitrfallElms_pft=0._r8;plt_bgcr%LitrFallStrutElms_col=0._r8
  plt_bgcr%LitrfallElms_pvr=0._r8;plt_bgcr%LitrfalStrutElms_vr=0._r8
  plt_bgcr%CO2NetFix_pft=0._r8;plt_bgcr%CO2NetFix_pft(1)=.2_r8
  plt_bgcr%Canopy_NEE_col=0._r8
  plt_bgcr%NH3Dep2Can_pft=0._r8;plt_bgcr%NH3Dep2Can_pft(1)=.04_r8
  plt_bgcr%NH3Emis_CumYr_pft=0._r8
  plt_bgcr%RootGasLossDisturb_pft=0._r8;plt_rbgc%TRootGasLossDisturb_col=0._r8
  plt_rbgc%PlantRootSoilElmNetX_pft=0._r8
  plt_rbgc%PlantRootSoilElmNetX_pft(:,1)=[.1_r8,.2_r8,.3_r8]
  plt_rbgc%trcg_rootml_pvr=0._r8;plt_rbgc%trcs_rootml_pvr=0._r8
  plt_rbgc%trcg_air2root_flx_vr=0._r8;plt_rbgc%trcg_air2root_flx_pvr=0._r8
  plt_rbgc%RCO2Emis2Root_rpvr=0._r8;plt_rbgc%RootO2Uptk_pvr=0._r8
  plt_rbgc%RootUptkSoiSol_pvr=0._r8;plt_rbgc%RootUptkSoiSol_pvr(idg_CO2,1,1,1)=.02_r8
  plt_rbgc%RootNutUptake_pvr=0._r8;plt_rbgc%RootNutUptake_pvr(ids_NH4B,1,1,1)=.03_r8
  plt_rbgc%Soil2RootMycoExudE_pvr=0._r8
  plt_rbgc%RootO2Dmnd4Resp_pvr=0._r8
  plt_rbgc%RootNH4DmndBand_pvr=0._r8;plt_rbgc%RootNH4DmndSoil_pvr=0._r8
  plt_rbgc%RootNO3DmndBand_pvr=0._r8;plt_rbgc%RootNO3DmndSoil_pvr=0._r8
  plt_rbgc%RootH1PO4DmndBand_pvr=0._r8;plt_rbgc%RootH1PO4DmndSoil_pvr=0._r8
  plt_rbgc%RootH2PO4DmndBand_pvr=0._r8;plt_rbgc%RootH2PO4DmndSoil_pvr=0._r8
  plt_bgcr%RootCO2Emis2Root_vr=0._r8;plt_bgcr%RootCO2Emis2Root_pvr=0._r8
  plt_bgcr%tRootMycoExud2Soil_vr=0._r8
  plt_rad%Eco_NetRad_col=0._r8;plt_rad%RadNet2Canopy_pft=0._r8;plt_rad%RadNet2Canopy_pft(1)=.5_r8
  plt_rad%LWRadCanopy_pft=0._r8
  plt_ew%CanopyBiomWater_pft=0._r8;plt_ew%WatHeldOnCanopy_pft=0._r8;plt_ew%SnowOnCanopy_pft=0._r8
  plt_ew%CanopyBiomWater_col=0._r8;plt_ew%WatHeldOnCanopy_col=0._r8;plt_ew%SnowOnCanopy_col=0._r8
  plt_ew%Transpiration_pft=0._r8;plt_ew%Transpiration_pft(1)=-.001_r8
  plt_ew%VapXAir2Canopy_pft=0._r8;plt_ew%ETCanopy_CumYr_pft=0._r8
  plt_ew%QVegET_col=0._r8;plt_ew%VapXAir2Canopy_col=0._r8
  plt_ew%Eco_Heat_Sens_col=0._r8;plt_ew%HeatXAir2PCan_pft=0._r8;plt_ew%HeatXAir2PCan_pft(1)=-.1_r8
  plt_ew%Eco_Heat_Latent_col=0._r8;plt_ew%CanopyEvapTransLHeat_pft=0._r8
  plt_ew%CanopyEvapTransLHeat_pft(1)=-.05_r8
  plt_ew%Eco_Heat_GrndSurf_col=0._r8;plt_ew%HeatStorCanopy_pft=0._r8;plt_ew%HeatStorCanopy_pft(1)=-.35_r8
  plt_ew%CanopyHeatStor_col=0._r8;plt_ew%HeatFlx2Canopy_col=0._r8;plt_ew%LWRadCanG=0._r8
  plt_ew%TKC_pft=280._r8;plt_ew%VHeatCapCanopy_pft=0._r8;plt_ew%ENGYX_pft=0._r8
  plt_ew%VHeatCapCanopy_pft(1)=.2_r8;plt_ew%ENGYX_pft(1)=55._r8
  plt_ew%CanopyHeatLoss2Dist_col=0._r8
  plt_ew%RPlantRootH2OUptk_pvr=0._r8;plt_ew%RPlantRootH2OUptk_pvr(1,1,1)=-.001_r8
  plt_ew%TWaterPlantRoot2Soil_vr=0._r8;plt_ew%THeatLossRoot2Soil_vr=0._r8
  plt_ew%TKS_vr=288._r8;plt_ew%TKCanopy_pft=280._r8

  call extracts(yearIJ,LiveExchange,PassiveExchange)
  call check_close('removed dry thermal mass is exported once',plt_ew%CanopyHeatLoss2Dist_col,56._r8)
  call check_close('removed dry capacity is zero',plt_ew%VHeatCapCanopy_pft(1),0._r8)
  call check_close('accepted root water after complete removal',plt_ew%TWaterPlantRoot2Soil_vr(1),-.001_r8)
  call check_close('accepted root heat',plt_ew%THeatLossRoot2Soil_vr(1),-.001_r8*cpw*288._r8)
  call check_close('accepted root CO2',plt_rbgc%trcs_Soil2plant_uptake_vr(idg_CO2,1),.02_r8)
  call check_close('accepted root N',plt_rbgc%trcs_Soil2plant_uptake_vr(ids_NH4B,1),.03_r8)
  call check_close('accepted canopy NH3',plt_bgcr%NH3Emis_CumYr_pft(1),.04_r8)
  call check_close('accepted photosynthesis',plt_bgcr%Canopy_NEE_col,.2_r8)
  call check_close('accepted canopy radiation',plt_rad%Eco_NetRad_col,.5_r8)
  call check_close('accepted canopy sensible heat',plt_ew%Eco_Heat_Sens_col,-.1_r8)
  call check_close('accepted canopy ET',plt_ew%QVegET_col,-.001_r8)
  if(MAXVAL(ABS(plt_bgcr%LitrFallStrutElms_col+[.1_r8,.2_r8,.3_r8]))>1.e-12_r8)error stop 'lost root C/N/P exchange'
  if(MAXVAL(ABS(plt_site%PlantElemntStoreLandscape-[5._r8,7._r8,9._r8]))>1.e-12_r8)error stop 'lost retained residuals'

  ! On the next hour no new live exchange was produced. Old rates must not
  ! be aggregated again just because the trait slot remains active/retained.
  LiveExchange=.false.
  call extracts(yearIJ,LiveExchange,PassiveExchange)
  call check_close('no repeated dry energy export',plt_ew%CanopyHeatLoss2Dist_col,56._r8)
  call check_close('no repeated canopy ET without current live exchange',plt_ew%QVegET_col,-.001_r8)
  call check_close('no repeated NH3 without current live exchange',plt_bgcr%NH3Emis_CumYr_pft(1),.04_r8)
  ! A passive flux already entered the air stores before all dead C was removed.
  ! Eligibility comes from that accepted solve, rather than the post-growth stock.
  PassiveExchange(1)=.true.
  plt_ew%Transpiration_pft=0._r8
  plt_ew%VapXAir2Canopy_pft(1)=.002_r8
  call extracts(yearIJ,LiveExchange,PassiveExchange)
  call check_close('passive sensible heat after dead removal',plt_ew%Eco_Heat_Sens_col,-.2_r8)
  call check_close('passive water after dead removal',plt_ew%QVegET_col,.001_r8)
  call check_close('passive branch has no biological exchange',plt_bgcr%NH3Emis_CumYr_pft(1),.04_r8)

  ! Dead material arriving after uptake must not resurrect old passive rates.
  PassiveExchange=.false.
  plt_biom%StandDeadCompKElms_pft(1,1,1)=1._r8
  plt_biom%StandDeadStrutElms_pft(1,1)=1._r8
  call extracts(yearIJ,LiveExchange,PassiveExchange)
  call check_close('new dead stock does not reuse sensible heat',plt_ew%Eco_Heat_Sens_col,-.2_r8)
  call check_close('new dead stock does not reuse vapor',plt_ew%QVegET_col,.001_r8)

  ! An inactive retained slot outside NP still owns its intercepted inventory.
  plt_ew%WatHeldOnCanopy_col=0._r8;plt_ew%SnowOnCanopy_col=0._r8
  plt_ew%WatHeldOnCanopy_pft(2)=.003_r8;plt_ew%SnowOnCanopy_pft(2)=.004_r8
  call extracts(yearIJ,LiveExchange,PassiveExchange)
  call check_close('retained liquid without exchange',plt_ew%WatHeldOnCanopy_col,.003_r8)
  call check_close('retained snow without exchange',plt_ew%SnowOnCanopy_col,.004_r8)
  write(iulog,*)'PASS: live/passive decisions, no stale rates, retained water/snow inventories'

contains
!------------------------------------------------------------------------------------------
  subroutine check_close(label,actual,expected)
    ! Description:
    ! Compare an independently specified accepted exchange against its aggregated value.
    character(len=*), intent(in) :: label ! process being checked [-]
    real(r8), intent(in) :: actual,expected ! aggregated and expected values in the process units
    character(len=*), parameter :: subname='check_close'
    call PrintInfo('beg '//subname)
    if(ABS(actual-expected)>1.e-12_r8)then
      write(iulog,*)label,actual,expected
      error stop 'accepted exchange was lost or duplicated'
    endif
    call PrintInfo('end '//subname)
  end subroutine check_close
end program CanopyAccountingTest
