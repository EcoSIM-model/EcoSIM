module SnowDataType
  ! Store snow constituents, layer geometry, temperatures and exchange diagnostics.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod
  use EcoSIMCtrlMod, only : salt_model
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target, allocatable ::  VLSnowHeatCapM_snvr(:,:,:,:)            !volumetric heat capacity of snowpack, [MJ/K d-2]
  real(r8),target, allocatable ::  WatFlowInSnowM_snvr(:,:,:,:)            !snowpack water flux, [m3 d-2 h-1]
  real(r8),target, allocatable ::  DrySnoFlxByRedistM_2DH(:,:,:,:,:)       !runoff snow flux, [m3 d-2 t-1]
  REAL(R8),target, allocatable ::  SoilAlbedo_col(:,:)                     !snowpack albedo, [-]
  REAL(R8),target, allocatable ::  SnowAlbedo_col(:,:)                     !snow albedo, [-]
  real(r8),target, allocatable ::  NewSnowDens_col(:,:)                    !new snowpack density, [Mg m-3]
  real(r8),target, allocatable ::  TCSnow_snvr(:,:,:)                      !snow temperature, [oC]
  real(r8),target, allocatable ::  TKSnow_snvr(:,:,:)                      !snow temperature, [K]
  real(r8),target, allocatable ::  VLHeatCapSnow_snvr(:,:,:)               !snowpack heat capacity, [MJ m-3 K-1]
  real(r8),target, allocatable ::  VLDrySnoWE_snvr(:,:,:)                  !water equivalent dry snow in snowpack layer, [m3 d-2]
  real(r8),target, allocatable ::  VLWatSnow_snvr(:,:,:)                   !snow water volume in snowpack layer, [m3 d-2]
  real(r8),target, allocatable ::  VLIceSnow_snvr(:,:,:)                   !snow ice volume in snowpack layer, [m3 d-2]
  real(r8),target, allocatable ::  VLSnoDWIprev_snvr(:,:,:)                !snow volume in snowpack layer, [m3 d-2]
  real(r8),target, allocatable ::  SnoDens_snvr(:,:,:)                     !snowpack density, [Mg m-3]
  real(r8),target, allocatable ::  SnowThickL_snvr(:,:,:)                  !snowpack layer thickness, [m]
  real(r8),target, allocatable ::  WatXfer2SnoLay_snvr(:,:,:)              !hourly snow water transfer, [m3 d-2 h-1]
  real(r8),target, allocatable ::  SnoXfer2SnoLay_snvr(:,:,:)              !hourly snow transfer to each layer,       [m d-2 h-1]
  real(r8),target, allocatable ::  IceXfer2SnoLay_snvr(:,:,:)              !hourly snow ice transfer to each layer,   [m d-2 h-1]
  real(r8),target, allocatable ::  HeatXfer2SnoLay_snvr(:,:,:)             !hourly convective heat flux from water transfer, [MJ d-2 h-1]
  integer ,target, allocatable ::  nsnol_col(:,:)                          !number of snow layers in column,[-]
  real(r8),target, allocatable ::  cumSnowDepz_col(:,:,:)                  !cumulative depth to bottom of snowpack layer,[m]
  real(r8),target, allocatable ::  VLSnoDWIMax_snvr(:,:,:)                 !maximum snowpack volume allowed in each layer, [m3 d-2]
  real(r8),target, allocatable ::  SnowDepth_col(:,:)                      !snowpack depth, [m]
  real(r8),target, allocatable ::  VcumDrySnoWE_col(:,:)                   !snow volume in snowpack (water equivalent), [m3 d-2]
  real(r8),target, allocatable ::  VcumWatSnow_col(:,:)                    !water volume in snowpack, [m3 d-2]
  real(r8),target, allocatable ::  VcumIceSnow_col(:,:)                    !ice volume in snowpack, [m3 d-2]
  real(r8),target, allocatable ::  VcumSnoDWI_col(:,:)                     !dry snow volume, [m3 d-2]
  real(r8),target, allocatable ::  VcumSnowWE_col(:,:)                     !water equivalent snowpack, [m3 d-2]
  real(r8),target, allocatable ::  VLHeatCapSnowMin_col(:,:)               !legacy snow layer/surface-property threshold, [MJ K-1 per column]
  real(r8),target, allocatable ::  WatConvSno2MicP_snvr(:,:,:)             !water from snowpack to soil micropores,[m3 d-2 h-1]
  real(r8),target, allocatable ::  WatConvSno2MacP_snvr(:,:,:)             !water from snowpack to soil macropores,[m3 d-2 h-1]
  real(r8),target, allocatable ::  HeatConvSno2Soi_snvr(:,:,:)             !convective heat from snowpack to soil,[MJ d-2 h-1]
  real(r8),target, allocatable ::  WatConvSno2LitR_snvr(:,:,:)             !water flux from snowpack to litter, [m3 d-2 h-1]
  real(r8),target, allocatable ::  HeatConvSno2LitR_snvr(:,:,:)            !convective heat flux from snowpack to litter, [MJ d-2 h-1]
  real(r8),target, allocatable ::  DrySnoByRedistrib_2DH(:,:,:,:)          !snowpack runoff snow, [m3 d-2 h-1]
  real(r8),target, allocatable ::  WatSnoByRedist_2DH(:,:,:,:)             !snowpack runoff water, [m3 d-2 h-1]
  real(r8),target, allocatable ::  IceSnoBySnowRedist_2DH(:,:,:,:)         !snowpack runoff ice, [m3 d-2 h-1]
  real(r8),target, allocatable ::  HeatSnoByRedist_2DH(:,:,:,:)            !snowpack runoff heat, [MJ d-2 h-1]
  real(r8),target, allocatable ::  trcg_FloXSnow_2DH(:,:,:,:,:)            !snowpack runoff CO2 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  trcn_FloXSnow_2DH(:,:,:,:,:)            !snowpack runoff NH4 flux, [g d-2 h-1]
  real(r8),target, allocatable ::  THeatSnowThaw_col(:,:)                  !total heat associated with phase change in snow, [MJ/d2/h]
  real(r8),target, allocatable ::  trcg_solsml_snvr(:,:,:,:)               ! Disolved volatile tracers in snow, [g d-2]
  real(r8),target, allocatable ::  trcn_solsml_snvr(:,:,:,:)               ! Dissolved nutrient tracers in snow, [g d-2]
  real(r8),target, allocatable ::  trcSalt_ml_snvr(:,:,:,:)                ! snowpack salt dissolved tracers, [mol d-2]
  real(r8),target, allocatable ::  SnowEngyBeg_col(:,:)                    !total snow-held energy at the beginning of the time step, [MJ d-2]
  real(r8),target, allocatable ::  SnowEngyEnd_col(:,:)                    !total snow-held energy at the end of the time step, [MJ d-2]
  real(r8),target, allocatable ::  SnowMassBeg_col(:,:)                    !snow mass H2O eqv at the beginning of the time step, [m3 H2O d-2]
  real(r8),target, allocatable ::  SnowMassEnd_col(:,:)                    !snow mass H2O eqv at the end of the time step, [m3 H2O d-2]
  real(r8),target, allocatable ::  trcSalt_FloXSnow_2DH(:,:,:,:,:)         !total salt in snow drift, [mol d-2 h-1]
  real(r8),target, allocatable ::  Prec2Snow_col(:,:)                      !precipiation to snow, [m3 H2O d-2 h-1]
  real(r8),target, allocatable ::  PrecHeat2Snow_col(:,:)                  !precipitation heat to snow, [MJ d-2 h-1]
  real(r8),target, allocatable ::  QSnowH2Oloss_col(:,:)                   !snow water eqv loss to other storage, [m3 H2O d-2 h-1]
  real(r8),target, allocatable ::  QSnowHeatLoss_col(:,:)                  !total heatloss from snow, [MJ d-2 h-1]
  real(r8),target, allocatable ::  trcg_AquaADV_Snow2Litr_flx(:,:,:)       !aqeuous volatile tracer from snow to litter, [g d-2 h-1]
  real(r8),target, allocatable ::  trcn_AquaADV_Snow2Litr_flx(:,:,:)       !aqeuous nutrient tracer from snow to litter, [g d-2 h-1]
  real(r8),target, allocatable ::  trcg_AquaADV_Snow2Soil_flx(:,:,:)       !aqueous volatile tracer from snow to soil, [g d-2 h-1]
  real(r8),target, allocatable ::  trcn_AquaADV_Snow2Soil_flx(:,:,:)       !aqueous nutrient tracer from snow to soil, [g d-2 h-1]
  real(r8),target, allocatable ::  trcn_AquaADV_Snow2Band_flx(:,:,:)       !aqueous nutrient tracer from snow to band soil, [g d-2 h-1]
  real(r8),target, allocatable ::  trcSalt_AquaADV_Snow2Soil_flx(:,:,:)    !salt flux from snow to soil, [mol d-2 h-1]
  real(r8),target, allocatable ::  trcSalt_AquaADV_Snow2Litr_flx(:,:,:)    !salt flux from snow to litter, [mol d-2 h-1]
  real(r8),target, allocatable ::  trcg_snowMass_beg_col(:,:,:)            !total mass of valatile tracer in snow at previous time step, [g d-2]
  real(r8),target, allocatable ::  trcg_snowMass_col(:,:,:)                !total mass of valatile tracer in snow, [g d-2]
  real(r8),target, allocatable ::  trcg_snowMassloss_col(:,:,:)            !total volatile mass of tracer loss from snow, [g d-2 h-1]
  real(r8),target, allocatable ::  trcn_snowMassloss_col(:,:,:)            !total nutrient mass of tracer loss from snow, [g d-2 h-1]
  real(r8),target, allocatable ::  trcSalt_snowMassloss_col(:,:,:)         !total salt mass of tracer loss from snow [g d-2 h-1]
  real(r8),target,allocatable ::   trcg_AquaAdv_flx_snvr(:,:,:,:)          !aqueous volatile tracer flux in snow, [g/d2/h]
  real(r8),target,allocatable ::   trcn_AquaAdv_flx_snvr(:,:,:,:)          !aqueous nutrient tracer flux in snow, [g/d2/h]
  real(r8),target,allocatable ::   trcSalt_AquaAdv_flx_snvr(:,:,:,:)       !aqueous salt tracer flux through snow, [g/d2/h]  
  integer, target,allocatable ::   IFLBSM_2DH(:,:,:,:,:)                   !flag for snow redistribution in intermediate iterations [-]
  real(r8),target,allocatable ::   trcn_SnowDrift_flx_col(:,:,:)           !nutrient loss from snow due to drifting [g/d2/h]
  real(r8),target,allocatable ::   trcg_SnowDrift_flx_col(:,:,:)           !volatile loss from snow due to drifting [g/d2/h]
  real(r8),target,allocatable ::   trcSalt_SnowDrift_flx_col(:,:,:)        !salt loss through snow drift [mol/d2/h]
  real(r8),target,allocatable :: Gas_WetDepo2Snow_col(:,:,:)               !wet deposition to snow, [g d-2 h-1] 
  real(r8),target,allocatable :: Gas_Snowloss_flx_col(:,:,:)               !tracer loss from snow [g d-2 h-1]  

!----------------------------------------------------------------------

contains
  subroutine InitSnowData

  implicit none
  allocate(Gas_Snowloss_flx_col(idg_beg:idg_NH3,JY,JX)); Gas_Snowloss_flx_col=0._r8  
  allocate(trcg_AquaAdv_flx_snvr(idg_beg:idg_NH3,JS,JY,JX)); trcg_AquaAdv_flx_snvr=0._r8
  allocate(trcn_AquaAdv_flx_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX)); trcn_AquaAdv_flx_snvr=0._r8
  allocate(trcg_snowMassloss_col(idg_beg:idg_NH3,JY,JX)); trcg_snowMassloss_col=0._r8
  allocate(trcn_snowMassloss_col(ids_nut_beg:ids_nuts_end,JY,JX)); trcn_snowMassloss_col=0._r8
  allocate(IFLBSM_2DH(60,2,2,JY,JX));IFLBSM_2DH=-1
  allocate(trcn_SnowDrift_flx_col(ids_nut_beg:ids_nuts_end,JY,JX));trcn_SnowDrift_flx_col=0._r8
  allocate(trcg_SnowDrift_flx_col(idg_beg:idg_NH3,JY,JX)); trcg_SnowDrift_flx_col=0._r8
  allocate(Gas_WetDepo2Snow_col(idg_beg:idg_NH3,JY,JX)); Gas_WetDepo2Snow_col=0._r8
  allocate(trcg_AquaADV_Snow2Litr_flx(idg_beg:idg_NH3,JY,JX)) ;trcg_AquaADV_Snow2Litr_flx=0._r8 
  allocate(trcn_AquaADV_Snow2Litr_flx(ids_nut_beg:ids_nuts_end,JY,JX));trcn_AquaADV_Snow2Litr_flx=0._r8
  allocate(trcn_AquaADV_Snow2Soil_flx(ids_nut_beg:ids_nuts_end,JY,JX)); trcn_AquaADV_Snow2Soil_flx=0._r8
  allocate(trcn_AquaADV_Snow2Band_flx(ids_nutb_beg:ids_nutb_end,JY,JX)); trcn_AquaADV_Snow2Band_flx=0._r8
  allocate(trcg_AquaADV_Snow2Soil_flx(idg_beg:idg_end,JY,JX)); trcg_AquaADV_Snow2Soil_flx=0._r8
  allocate(trcg_snowMass_beg_col(idg_beg:idg_NH3,JY,JX)); trcg_snowMass_beg_col=0._r8
  allocate(trcg_snowMass_col(idg_beg:idg_NH3,JY,JX)); trcg_snowMass_col=0._r8
  allocate(QSnowHeatLoss_col(JY,JX)); QSnowHeatLoss_col=0._r8
  allocate(QSnowH2Oloss_col(JY,JX)); QSnowH2Oloss_col=0._r8
  allocate(PrecHeat2Snow_col(JY,JX)); PrecHeat2Snow_col=0._r8
  allocate(Prec2Snow_col(JY,JX)); Prec2Snow_col=0._r8
  allocate(SnowMassBeg_col(JY,JX)); SnowMassBeg_col=0._r8
  allocate(SnowMassEnd_col(JY,JX)); SnowMassEnd_col=0._r8
  allocate(SnowEngyEnd_col(JY,JX)); SnowEngyEnd_col=0._r8
  allocate(SnowEngyBeg_col(JY,JX)); SnowEngyBeg_col=0._r8
  allocate(THeatSnowThaw_col(JY,JX))   ; THeatSnowThaw_col=0._r8
  allocate(VLSnowHeatCapM_snvr(60,JS,JY,JX));VLSnowHeatCapM_snvr=0._r8
  allocate(WatFlowInSnowM_snvr(60,JS,JY,JX)); WatFlowInSnowM_snvr=0._r8
  allocate(DrySnoFlxByRedistM_2DH(60,2,2,JV,JH));    DrySnoFlxByRedistM_2DH=0._r8
  allocate(SoilAlbedo_col(JY,JX));        SoilAlbedo_col=0._r8
  allocate(SnowAlbedo_col(JY,JX));        SnowAlbedo_col=0._r8
  allocate(NewSnowDens_col(JY,JX));       NewSnowDens_col=0._r8
  allocate(TCSnow_snvr(JS,JY,JX));      TCSnow_snvr=0._r8
  allocate(TKSnow_snvr(JS,JY,JX));      TKSnow_snvr=0._r8
  allocate(VLHeatCapSnow_snvr(JS,JY,JX));    VLHeatCapSnow_snvr=0._r8
  allocate(VLDrySnoWE_snvr(JS,JY,JX));   VLDrySnoWE_snvr=0._r8
  allocate(VLWatSnow_snvr(JS,JY,JX));   VLWatSnow_snvr=0._r8
  allocate(VLIceSnow_snvr(JS,JY,JX));   VLIceSnow_snvr=0._r8
  allocate(VLSnoDWIprev_snvr(JS,JY,JX));    VLSnoDWIprev_snvr=0._r8
  allocate(SnoDens_snvr(JS,JY,JX));    SnoDens_snvr=0._r8
  allocate(SnowThickL_snvr(JS,JY,JX));    SnowThickL_snvr=0._r8
  allocate(WatXfer2SnoLay_snvr(JS,JY,JX));    WatXfer2SnoLay_snvr=0._r8
  allocate(SnoXfer2SnoLay_snvr(JS,JY,JX));    SnoXfer2SnoLay_snvr=0._r8
  allocate(IceXfer2SnoLay_snvr(JS,JY,JX));    IceXfer2SnoLay_snvr=0._r8
  allocate(HeatXfer2SnoLay_snvr(JS,JY,JX));   HeatXfer2SnoLay_snvr=0._r8
  allocate(nsnol_col(JY,JX)); nsnol_col=0
  allocate(cumSnowDepz_col(0:JS,JY,JX)); cumSnowDepz_col=0._r8
  allocate(VLSnoDWIMax_snvr(JS,JY,JX));    VLSnoDWIMax_snvr=0._r8
  allocate(SnowDepth_col(JY,JX));       SnowDepth_col=0._r8
  allocate(VcumSnowWE_col(JY,JX));   VcumSnowWE_col=0._r8
  allocate(VcumDrySnoWE_col(JY,JX));       VcumDrySnoWE_col=0._r8
  allocate(VcumWatSnow_col(JY,JX));       VcumWatSnow_col=0._r8
  allocate(VcumIceSnow_col(JY,JX));       VcumIceSnow_col=0._r8
  allocate(VcumSnoDWI_col(JY,JX));        VcumSnoDWI_col=0._r8
  allocate(VLHeatCapSnowMin_col(JY,JX));      VLHeatCapSnowMin_col=0._r8
  allocate(WatConvSno2MicP_snvr(JS,JY,JX));     WatConvSno2MicP_snvr=0._r8
  allocate(WatConvSno2MacP_snvr(JS,JY,JX));    WatConvSno2MacP_snvr=0._r8
  allocate(HeatConvSno2Soi_snvr(JS,JY,JX));    HeatConvSno2Soi_snvr=0._r8
  allocate(WatConvSno2LitR_snvr(JS,JY,JX));    WatConvSno2LitR_snvr=0._r8
  allocate(HeatConvSno2LitR_snvr(JS,JY,JX));   HeatConvSno2LitR_snvr=0._r8
  allocate(DrySnoByRedistrib_2DH(2,2,JV,JH));        DrySnoByRedistrib_2DH=0._r8
  allocate(WatSnoByRedist_2DH(2,2,JV,JH));        WatSnoByRedist_2DH=0._r8
  allocate(IceSnoBySnowRedist_2DH(2,2,JV,JH));        IceSnoBySnowRedist_2DH=0._r8
  allocate(HeatSnoByRedist_2DH(2,2,JV,JH));       HeatSnoByRedist_2DH=0._r8
  allocate(trcg_FloXSnow_2DH(idg_beg:idg_NH3,2,2,JV,JH));    trcg_FloXSnow_2DH=0._r8
  allocate(trcn_FloXSnow_2DH(ids_nut_beg:ids_nuts_end,2,2,JV,JH));    trcn_FloXSnow_2DH=0._r8

! exclude NH3B
  allocate(trcg_solsml_snvr(idg_beg:idg_NH3,JS,JY,JX));trcg_solsml_snvr=0._r8
  allocate(trcn_solsml_snvr(ids_nut_beg:ids_nuts_end,JS,JY,JX));trcn_solsml_snvr=0._r8
  if(salt_model)then
    allocate(trcSalt_ml_snvr(idsalt_beg:idsalt_end,JS,JY,JX)); trcSalt_ml_snvr=0._r8
    allocate(trcSalt_AquaADV_Snow2Litr_flx(idsalt_beg:idsalt_end,JY,JX));  trcSalt_AquaADV_Snow2Litr_flx=0._r8
    allocate(trcSalt_AquaADV_Snow2Soil_flx(idsalt_beg:idsaltb_end,JY,JX)); trcSalt_AquaADV_Snow2Soil_flx=0._r8
    allocate(trcSalt_FloXSnow_2DH(idsalt_beg:idsalt_end,2,2,JV,JH));     trcSalt_FloXSnow_2DH=0._r8
    allocate(trcSalt_snowMassloss_col(idsalt_beg:idsalt_end,JY,JX)); trcSalt_snowMassloss_col=0._r8
    allocate(trcSalt_AquaAdv_flx_snvr(idsalt_beg:idsalt_end,JS,JY,JX)); trcSalt_AquaAdv_flx_snvr=0._r8    
    allocate(trcSalt_SnowDrift_flx_col(idsalt_beg:idsalt_end,JY,JX));  trcSalt_SnowDrift_flx_col=0._r8
  endif


  end subroutine InitSnowData

!----------------------------------------------------------------------
  subroutine DestructSnowData
  use abortutils, only : destroy
  implicit none
  if(salt_model)then
    call destroy(trcSalt_ml_snvr)
    call destroy(trcSalt_FloXSnow_2DH)
    call destroy(trcSalt_AquaADV_Snow2Litr_flx)
    call destroy(trcSalt_AquaADV_Snow2Soil_flx)
    call destroy(trcSalt_snowMassloss_col)
    call destroy(trcSalt_AquaAdv_flx_snvr)
    call destroy(trcSalt_SnowDrift_flx_col)    
  endif

  call destroy(trcg_AquaAdv_flx_snvr)
  call destroy(trcn_AquaAdv_flx_snvr)
  call destroy(IFLBSM_2DH)
  call destroy(trcn_AquaADV_Snow2Soil_flx)
  call destroy(trcn_AquaADV_Snow2Band_flx)
  call destroy(trcg_AquaADV_Snow2Soil_flx)
  call destroy(trcg_snowMassloss_col)
  call destroy(trcn_snowMassloss_col)
  call destroy(trcg_AquaADV_Snow2Litr_flx)
  call destroy(trcn_AquaADV_Snow2Litr_flx)
  call destroy(QSnowHeatLoss_col)
  call destroy(QSnowH2Oloss_col)
  call destroy(PrecHeat2Snow_col)
  call destroy(Prec2Snow_col)
  call destroy(SnowMassEnd_col)
  call destroy(SnowMassBeg_col)
  call destroy(SnowEngyEnd_col)
  call destroy(SnowEngyBeg_col)
  call destroy(VLSnowHeatCapM_snvr)
  call destroy(WatFlowInSnowM_snvr)
  call destroy(DrySnoFlxByRedistM_2DH)
  call destroy(SoilAlbedo_col)
  call destroy(SnowAlbedo_col)
  call destroy(NewSnowDens_col)
  call destroy(trcg_solsml_snvr)
  call destroy(TCSnow_snvr)
  call destroy(TKSnow_snvr)
  call destroy(VLHeatCapSnow_snvr)
  call destroy(VLDrySnoWE_snvr)
  call destroy(VLWatSnow_snvr)
  call destroy(VLIceSnow_snvr)
  call destroy(VLSnoDWIprev_snvr)
  call destroy(SnoDens_snvr)
  call destroy(SnowThickL_snvr)
  call destroy(WatXfer2SnoLay_snvr)
  call destroy(SnoXfer2SnoLay_snvr)
  call destroy(IceXfer2SnoLay_snvr)
  call destroy(HeatXfer2SnoLay_snvr)
  call destroy(nsnol_col)
  call destroy(cumSnowDepz_col)
  call destroy(VLSnoDWIMax_snvr)
  call destroy(SnowDepth_col)
  call destroy(VcumSnowWE_col)
  call destroy(VcumDrySnoWE_col)
  call destroy(VcumWatSnow_col)
  call destroy(VcumIceSnow_col)
  call destroy(VcumSnoDWI_col)
  call destroy(VLHeatCapSnowMin_col)
  call destroy(WatConvSno2MicP_snvr)
  call destroy(WatConvSno2MacP_snvr)
  call destroy(HeatConvSno2Soi_snvr)
  call destroy(WatConvSno2LitR_snvr)
  call destroy(HeatConvSno2LitR_snvr)
  call destroy(DrySnoByRedistrib_2DH)
  call destroy(WatSnoByRedist_2DH)
  call destroy(IceSnoBySnowRedist_2DH)
  call destroy(HeatSnoByRedist_2DH)
  call destroy(trcg_FloXSnow_2DH)
  call destroy(trcn_FloXSnow_2DH)
  call destroy(trcn_solsml_snvr)
  call destroy(THeatSnowThaw_col)
  call destroy(trcg_snowMass_col)
  call destroy(trcg_snowMass_beg_col)
  call destroy(trcn_SnowDrift_flx_col)
  call destroy(trcg_SnowDrift_flx_col)
  call destroy(Gas_WetDepo2Snow_col)
  call destroy(Gas_Snowloss_flx_col)  
  end subroutine DestructSnowData

end module SnowDataType
