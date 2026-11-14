program RootUptakeSplitTest
  ! Description:
  ! Verify that hydraulic preparation preserves chemical bookkeeping and that
  ! the later chemical reset preserves accepted canopy/root water and heat.
  use data_kind_mod, only: r8 => DAT_KIND_R8,yearIJ_type
  use EcoSiMParDataMod, only: micpar,pltpar
  use InitAllocMod, only: InitPlantMorphSize
  use PlantAPIData
  use TracerIDMod, only: InitTracerIDs
  use UptakesMod, only: InitUptake,RootWaterUptake,RootChemicalUptake
  use EcoSIMCtrlMod, only: ldo_sp_mode
  use EcoSIMSolverPar, only: NPH
  use DebugToolMod, only: PrintInfo
  use abortutils, only: iulog
  implicit none
  type(yearIJ_type) :: yearIJ
  logical, allocatable :: Live(:),Passive(:)
  real(r8), allocatable :: DryMatter(:),Path(:,:,:),Radius(:,:,:),Share(:,:,:),Depth(:,:),Length(:,:,:)

  call micpar%Init()
  call InitPlantMorphSize()
  pltpar%JP1=2;pltpar%JZ1=2;pltpar%jroots=2
  pltpar%NumCanopyLayers1=2;pltpar%NumOfLeafAzimuthSectors=1;pltpar%NumOfSkyAzimuthSects1=1
  pltpar%NumLeafInclinationClasses1=1;pltpar%MaxNodesPerBranch1=2
  call InitTracerIDs(.false.)
  call InitPlantAPIData()
  call InitUptake()
  allocate(Live(JP1),Passive(JP1),DryMatter(JP1),Depth(JZ1,JP1))
  allocate(Path(jroots,JZ1,JP1),Radius(jroots,JZ1,JP1),Share(jroots,JZ1,JP1),Length(jroots,JZ1,JP1))
  yearIJ%year=1980;yearIJ%I=1;yearIJ%J=12
  NPH=1;ldo_sp_mode=.true.
  plt_site%NP=0;plt_site%NU=1;plt_site%NK=2;plt_site%ALT=0._r8;plt_site%ZERO=1.e-12_r8
  plt_site%VLWatMicPM_vr=.2_r8
  plt_ew%ElvAdjstedSoilH2OPSIMPa_vr=-.1_r8
  plt_soilchem%SoilBulkDensity_vr=1._r8;plt_soilchem%SoilWatAirDry_vr=.01_r8
  plt_soilchem%VLSoilMicP_vr=1._r8;plt_soilchem%VLMicP_vr=1._r8
  plt_soilchem%VLWatMicP_vr=.2_r8;plt_soilchem%VLiceMicP_vr=0._r8
  plt_rbgc%RootNutUptake_pvr=7._r8
  plt_rbgc%RootO2Dmnd4Resp_pvr=11._r8
  plt_rbgc%RootSoilCO2Debt_vr=13._r8
  plt_ew%RPlantRootH2OUptk_pvr=-.01_r8
  plt_ew%TWaterPlantRoot2Soil_vr=99._r8
  plt_ew%THeatLossRoot2Soil_vr=999._r8
  call RootWaterUptake(yearIJ,Live,Passive,DryMatter,Path,Radius,Share,Depth,Length)
  call check('water step preserves nutrient uptake',all(plt_rbgc%RootNutUptake_pvr==7._r8))
  call check('water step preserves previous O2 demand',all(plt_rbgc%RootO2Dmnd4Resp_pvr==11._r8))
  call check('water step preserves CO2 debt',all(plt_rbgc%RootSoilCO2Debt_vr==13._r8))
  call check('water step clears its own flux',all(plt_ew%RPlantRootH2OUptk_pvr==0._r8))
  call check('water step clears totals including unused layers',all(plt_ew%TWaterPlantRoot2Soil_vr==0._r8))
  call check('water step clears heat including unused layers',all(plt_ew%THeatLossRoot2Soil_vr==0._r8))
  call check('bare slots have no exchange',.not.any(Live.or.Passive))

  ! An inactive column still enters the chemical reset after hydraulic exchange.
  ! Accepted fluxes here represent a previous live/passive solve in the same hour.
  ldo_sp_mode=.false.;plt_site%NP=2
  plt_ew%RPlantRootH2OUptk_pvr=-.002_r8;plt_ew%RootH2OUptkStress_pvr=-.003_r8
  plt_ew%Transpiration_pft=-.004_r8;plt_ew%VapXAir2Canopy_pft=.001_r8
  plt_ew%HeatXAir2PCan_pft=.5_r8;plt_ew%HeatStorCanopy_pft=.6_r8
  plt_ew%CanopyEvapTransLHeat_pft=.7_r8;plt_rad%RadNet2Canopy_pft=.8_r8
  plt_rad%LWRadCanopy_pft=.9_r8
  plt_ew%TWaterPlantRoot2Soil_vr=-.002_r8
  plt_ew%THeatLossRoot2Soil_vr=-.003_r8
  call RootChemicalUptake(yearIJ,Live,DryMatter,Path,Radius,Share,Depth,Length)
  call check('chemical step preserves accepted root water',all(plt_ew%RPlantRootH2OUptk_pvr==-.002_r8))
  call check('chemical step preserves total root water',all(plt_ew%TWaterPlantRoot2Soil_vr==-.002_r8))
  call check('chemical step preserves total donor heat',all(plt_ew%THeatLossRoot2Soil_vr==-.003_r8))
  call check('chemical step preserves root stress diagnostic',all(plt_ew%RootH2OUptkStress_pvr==-.003_r8))
  call check('chemical step preserves transpiration',all(plt_ew%Transpiration_pft==-.004_r8))
  call check('chemical step preserves wet canopy exchange',all(plt_ew%VapXAir2Canopy_pft==.001_r8))
  call check('chemical step preserves sensible heat',all(plt_ew%HeatXAir2PCan_pft==.5_r8))
  call check('chemical step preserves stored heat',all(plt_ew%HeatStorCanopy_pft==.6_r8))
  call check('chemical step preserves latent heat',all(plt_ew%CanopyEvapTransLHeat_pft==.7_r8))
  call check('chemical step preserves net radiation',all(plt_rad%RadNet2Canopy_pft==.8_r8))
  call check('chemical step preserves longwave',all(plt_rad%LWRadCanopy_pft==.9_r8))
  call check('chemical step resets nutrient outputs',all(plt_rbgc%RootNutUptake_pvr==0._r8))
  call check('chemical step clears inactive O2 demand',all(plt_rbgc%RootO2Dmnd4Resp_pvr==0._r8))
  call check('chemical step resets transient debt',all(plt_rbgc%RootSoilCO2Debt_vr==0._r8))
contains
!------------------------------------------------------------------------------------------
  subroutine check(Label,Condition)
    ! Description:
    ! Reject a violated hydraulic/chemical state boundary.
    character(len=*), intent(in) :: Label ! diagnostic description [-]
    logical, intent(in) :: Condition ! true when the boundary contract holds [-]
    character(len=*), parameter :: subname='check'
    call PrintInfo('beg '//subname)
    if(.not.Condition)then
      write(iulog,*)Label
      error stop 1
    endif
    call PrintInfo('end '//subname)
  end subroutine check
end program RootUptakeSplitTest
