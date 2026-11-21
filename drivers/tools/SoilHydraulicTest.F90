program SoilHydraulicTest
  ! Description:
  ! Characterize the existing hydraulic scheme before introducing alternatives.
  ! Frozen pressure/conductivity samples come from pre-refactor GNU 13.4.0 output;
  ! independent bounds, monotonicity and domain checks exercise the public API.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: iulog
  use DebugToolMod, only: PrintInfo
  use EcoSIMCtrlMod, only: lverb,use_ecosys_drainage,soil_hydraulic_scheme
  use PhysPars, only: soil_bet_c,litter_bet_c
  use ElmIDMod, only: ielmc,isoi_set,isoi_unset,isoi_scnv,isoi_scnh
  use EcoSIMCtrlDataType
  use SoilPhysDataType
  use SoilWaterDataType
  use SoilPropertyDataType
  use HydroThermData
  use SoilHeatDataType
  use FlagDataType
  use SOMDataType
  use SurfLitterDataType
  use SoilHydroParaMod, only: SoilHydroProperty
  use LitterHydraulicMod, only: LitterHydroproperty,LitterHydraulicParameters, &
    GetLitterHydraulicParameters,UpdateLitterHydraulicLimits,LitterMatricPotential, &
    LitterWaterContent,LitterDryContent,LitterMatricPotentialDerivative,LitterMoistureClass,LitterResistance,LitterResistanceRetention
  use SoilPhysParaMod, only: ComputePsiMCM,ComputeSoilAirDryContent, &
    ComputeKsatReferenceContent,getMoistK
  use HydrologyDiagMod, only: ComputeWaterTableEquilContent
  use SoilHydraulicModelMod, only: SoilHydraulicModel,SoilHydraulicParameters
  use EcoSIMSoilHydraulicMod, only: EcoSIMSoilHydraulicModel
  use BETResistanceMod, only: BETResistanceParameters,MakeBETResistanceParameters,BETResistanceState, &
    BETResistanceHumidity,BETReferenceTemperature,SoilBETCutoffFCFraction,LitterBETCutoffFCFraction
  use PorousSurfaceResistanceMod, only: PorousSurfaceResistance
  use EcoSimConst, only: RGASC
  use MiniMathMod, only: vapsat
  use SoilHydraulicSchemeMod, only: InitializeSoilHydraulicScheme,FinalizeSoilHydraulicScheme, &
    GetSoilHydraulicModel,GetSoilHydraulicSchemeName,GetSoilHydraulicParameters
  implicit none
  integer, parameter :: classes(4)=[1,11,51,100]
  real(r8), parameter :: expected_pressure(4,4)=reshape([ &
    -0.50000000000000001e-3_r8, -0.14851842118345370e-2_r8, &
    -0.16086837886916330_r8, -2500.0000000000000_r8, &
    -0.50000000000000001e-3_r8, -0.73333162275748256e-2_r8, &
    -0.16086837886916330_r8, -2500.0000000000000_r8, &
    -0.50000000000000001e-3_r8, -0.99557317652258850e-2_r8, &
    -0.16086837886916330_r8, -2500.0000000000000_r8, &
    -0.50000000000000001e-3_r8, -0.27326365358774710e-2_r8, &
    -0.16086837886916330_r8, -2500.0000000000000_r8],[4,4])
  real(r8), parameter :: expected_vertical(4,4)=reshape([ &
    1.0000000000000000_r8, 0.81444774791937638e-1_r8, &
    0.30663774453083744e-5_r8, 0.19131935130694238e-17_r8, &
    1.0000000000000000_r8, 0.11533685387620114_r8, &
    0.38293720408472745e-4_r8, 0.23892459028120865e-16_r8, &
    1.0000000000000000_r8, 0.13075448922795219_r8, &
    0.55779091774690179e-4_r8, 0.34802041970247154e-16_r8, &
    1.0000000000000000_r8, 0.77327840747116947e-1_r8, &
    0.80430917995973638e-5_r8, 0.50182964525634386e-17_r8],[4,4])
  real(r8), parameter :: expected_entry(4)=[ &
    -0.15000000000000000e-2_r8, &
    -0.15000000000000000e-2_r8, &
    -0.15000000000000000e-2_r8, &
    -0.21358236677699745e-2_r8]
  real(r8), parameter :: expected_theta_sat(4)=[ &
    0.40050000000000002_r8, &
    0.44550000000000001_r8, &
    0.44550000000000001_r8, &
    0.41400000000000003_r8]
  real(r8), parameter :: litter_entry(2)=[ &
    -0.15000000000000000e-2_r8, &
    -0.11602176772093887e-2_r8]
  real(r8), parameter :: litter_theta_sat(2)=[ &
    0.66400000000000003_r8, &
    0.69600000000000006_r8]
  real(r8) :: theta,psi,content
  integer :: mode,k,j
  type(SoilHydraulicParameters) :: parameters
  class(SoilHydraulicModel), pointer :: model

  call initialize_state()
  call check_bet_resistance_state()
  call InitializeSoilHydraulicScheme()
  call require('default scheme identity',GetSoilHydraulicSchemeName()=='ecosim')
  model=>GetSoilHydraulicModel()
  select type(model)
  type is(EcoSIMSoilHydraulicModel)
    call require('default concrete scheme',.true.)
  class default
    call require('default concrete scheme',.false.)
  end select
  do mode=1,4
    CSoilOrgM_vr=0._r8
    if(mode==2)CSoilOrgM_vr(ielmc,1,1,1)=60000._r8
    if(mode==3)CSoilOrgM_vr(ielmc,1,1,1)=150000._r8
    use_ecosys_drainage=mode==4
    PSISoilAirEntry_vr=-.0015_r8;ThetaSat_vr=0._r8
    call SoilHydroProperty(1,1,1,1,1)
    do j=1,4
      k=classes(j)
      theta=.45_r8*(1._r8-real(k-1,r8)/100._r8)
      call ComputePsiMCM(1,1,1,theta,psi)
      call check_close('frozen pressure sample',psi,expected_pressure(j,mode))
      call check_close('frozen conductivity sample',HydroCond_3D(3,k,1,1,1),expected_vertical(j,mode))
    enddo
    call check_close('air-entry policy',PSISoilAirEntry_vr(1,1,1),expected_entry(mode))
    call check_close('ThetaSat policy',ThetaSat_vr(1,1,1),expected_theta_sat(mode))
    parameters=GetSoilHydraulicParameters(1,1,1,1.e-6_r8)
    if(mode==1)call check_resistance_cache(model,parameters)
    call check_retention_derivative(model,parameters)
    call check_soil_resistance(model,parameters)
    do j=1,11
      theta=.1_r8+.034_r8*(j-1)
      psi=model%matric_potential(parameters,theta)
      call check_close('forward/inverse retention',model%water_content(parameters,psi),theta)
    enddo
    call check_close('inverse dry bound',model%water_content(parameters,-25000._r8),SoilWatAirDry_vr(1,1,1))
    call check_close('inverse saturation bound',model%water_content(parameters,0._r8),.45_r8)
    call check_close('relative conductivity dispatch',model%relative_conductivity(parameters,.225_r8), &
      expected_vertical(3,mode))
    do k=1,100
      call check_close('horizontal anisotropy',HydroCond_3D(1,k,1,1,1),2._r8*HydroCond_3D(3,k,1,1,1))
      call require('finite positive conductivity',ieee_is_finite(HydroCond_3D(3,k,1,1,1)) &
        .and. HydroCond_3D(3,k,1,1,1)>0._r8)
      if(k>1)call require('wet-to-dry monotonicity',HydroCond_3D(3,k,1,1,1)<=HydroCond_3D(3,k-1,1,1,1))
    enddo
    theta=.30_r8
    call ComputePsiMCM(1,1,1,theta,psi)
    call check_close('FC anchor',psi,-.033_r8)
    theta=.15_r8
    call ComputePsiMCM(1,1,1,theta,psi)
    call check_close('WP anchor',psi,-1.5_r8)
    theta=1.e-6_r8
    call ComputePsiMCM(1,1,1,theta,psi)
    call check_close('dry pressure bound',psi,-2500._r8)
    theta=.45_r8-.5e-6_r8
    call ComputePsiMCM(1,1,1,theta,psi)
    call check_close('near-saturation pressure',psi,-.0005_r8)
    call require('near-saturation content adjustment',theta==.45_r8)
  enddo

  call ComputeSoilAirDryContent(1,1,1,content)
  call check_close('dry content characterization',content,.038992079977186643_r8)
  content=ComputeKsatReferenceContent(1,1,1)
  call check_close('Ksat reference content',content,.30_r8)
  call ComputeWaterTableEquilContent(1,1,1,-.033_r8,.44_r8,content)
  call check_close('legacy equilibrium approximation',content,.30_r8)
  call ComputeWaterTableEquilContent(1,1,1,-1.5_r8,.44_r8,content)
  call check_close('preserve extrapolation below FC',content,.2073503394434138_r8)
  call ComputeWaterTableEquilContent(1,1,1,0._r8,.44_r8,content)
  call check_close('saturated equilibrium capacity',content,.44_r8)
  call ComputeWaterTableEquilContent(1,1,1,-.033_r8,.20_r8,content)
  call check_close('equilibrium liquid capacity cap',content,.20_r8)
  call require('dry lookup clip',getMoistK(0._r8,.45_r8)==100)
  call require('lookup below half saturation',getMoistK(.499_r8*.45_r8,.45_r8)==52)
  call require('lookup at half saturation',getMoistK(.5_r8*.45_r8,.45_r8)==51)
  call require('saturated lookup',getMoistK(.45_r8,.45_r8)==1)

  ! Missing conductivity directions are estimated independently of supplied ones.
  ISOIL_vr(isoi_scnv,1,1,1)=isoi_unset
  call SoilHydroProperty(1,1,1,1,1)
  call check_close('estimated vertical Ksat',SatHydroCondVert_vr(1,1,1),.385_r8)
  call check_close('supplied horizontal Ksat retained',SatHydroCondHrzn_vr(1,1,1),2._r8)
  ISOIL_vr=isoi_set
  ISOIL_vr(isoi_scnh,1,1,1)=isoi_unset
  SatHydroCondVert_vr=1._r8
  call SoilHydroProperty(1,1,1,1,1)
  call check_close('supplied vertical Ksat retained',SatHydroCondVert_vr(1,1,1),1._r8)
  call check_close('estimated horizontal Ksat',SatHydroCondHrzn_vr(1,1,1),.385_r8)
  ISOIL_vr=isoi_set
  SatHydroCondHrzn_vr=2._r8

  ! Threshold scans without a crossing must retain the initialized/restored state.
  use_ecosys_drainage=.false.
  PSISoilAirEntry_vr=-3000._r8;ThetaSat_vr=.123_r8
  call SoilHydroProperty(1,1,1,1,1)
  call check_close('no pressure crossing retains ThetaSat',ThetaSat_vr(1,1,1),.123_r8)
  use_ecosys_drainage=.true.
  SatHydroCondVert_vr=0._r8
  call SoilHydroProperty(1,1,1,1,1)
  call check_close('zero Ksat retains ThetaSat',ThetaSat_vr(1,1,1),.123_r8)
  call check_close('zero Ksat retains air entry',PSISoilAirEntry_vr(1,1,1),-3000._r8)

  ! Litter uses its own porosity and r=1, even with the legacy drainage policy.
  fLitRPore4WatHold_col=.8_r8;POROS_vr(0,1,1)=.8_r8
  FieldCapacity_vr(0,1,1)=.4_r8;WiltPoint_vr(0,1,1)=.1_r8
  LOGPOROS_vr(0,1,1)=log(.8_r8);LOGFldCapacity_vr(0,1,1)=log(.4_r8)
  LOGWiltPoint_vr(0,1,1)=log(.1_r8)
  PSD_vr(0,1,1)=log(.8_r8)-log(.4_r8);FCD_vr(0,1,1)=log(.4_r8)-log(.1_r8)
  SRP_vr(0,1,1)=1._r8;SatHydroCondVert_vr(0,1,1)=.98_r8
  do mode=1,2
    use_ecosys_drainage=mode==2
    PSISoilAirEntry_vr(0,1,1)=-.0015_r8;ThetaSat_vr(0,1,1)=0._r8
    call LitterHydroproperty(1,1,use_ecosys_drainage)
    call check_close('litter air-entry policy',PSISoilAirEntry_vr(0,1,1),litter_entry(mode))
    call check_close('litter ThetaSat policy',ThetaSat_vr(0,1,1),litter_theta_sat(mode))
    call check_close('litter curvature retained',SRP_vr(0,1,1),1._r8)
    call require('no horizontal litter transport',all(HydroCond_3D(1:2,:,0,1,1)==0._r8))
    call check_close('saturated litter conductivity',HydroCond_3D(3,1,0,1,1),.98_r8)
  enddo
  call check_litter_independence()
  call FinalizeSoilHydraulicScheme()
  write(*,*) 'SoilHydraulicTest passed'
contains
!------------------------------------------------------------------------------------------
  subroutine check_resistance_cache(scheme,context)
    ! Description:
    ! Exercise initialization, adjustable C invalidation and deterministic native
    ! anchor reconstruction. Check the prescribed wet cutoff against native data.
    class(SoilHydraulicModel), intent(in) :: scheme ! selected native model [-]
    type(SoilHydraulicParameters), intent(in) :: context ! initialized one-layer view
    type(SoilHydraulicParameters) :: refreshed
    type(BETResistanceParameters) :: rebuilt
    real(r8) :: work,pressure,log_g,log_ratio,native_slope,extension_pressure,saved_fc(4)
    character(len=*), parameter :: subname='check_resistance_cache'

    call PrintInfo('beg '//subname)
    associate( &
      native => scheme, & ! input: selected native hydraulic implementation [-]
      material => context) ! input: initialized layer view and cached resistance material [-]
    call require('soil material cache initialized',material%resistance%ready)
    call check_close('soil cutoff is 0.6 field capacity',material%resistance%theta_wet,.6_r8*material%field_capacity)
    rebuilt=native%resistance_parameters(material)
    call require('soil cache rebuild is deterministic',all(rebuilt%coefficient==material%resistance%coefficient))
    work=.6_r8*material%field_capacity
    pressure=native%matric_potential(context,work)
    native_slope=native%matric_potential_derivative(context,work)
    call BETResistanceState(material%resistance,work,273.15_r8,extension_pressure,log_g,log_ratio)
    call require('soil cutoff pressure matches native at different T',ABS(extension_pressure/pressure-1._r8)<1.e-12_r8)
    call require('soil cutoff slope matches native at different T',ABS(EXP(log_g)/native_slope-1._r8)<1.e-12_r8)
    soil_bet_c=24._r8
    refreshed=GetSoilHydraulicParameters(1,1,1,material%saturation_tolerance)
    call require('soil C override refreshes cache',refreshed%resistance%adsorption_c==24._r8)
    call require('soil C override changes monolayer',refreshed%resistance%monolayer/=material%resistance%monolayer)
    soil_bet_c=18._r8
    refreshed=GetSoilHydraulicParameters(1,1,1,material%saturation_tolerance)
    call require('soil cache reset is deterministic',all(refreshed%resistance%coefficient==material%resistance%coefficient))
    saved_fc=[FieldCapacity_vr(1,1,1),LOGFldCapacity_vr(1,1,1),PSD_vr(1,1,1),FCD_vr(1,1,1)]
    FieldCapacity_vr(1,1,1)=.9_r8*material%field_capacity
    LOGFldCapacity_vr(1,1,1)=LOG(FieldCapacity_vr(1,1,1))
    PSD_vr(1,1,1)=LOGPOROS_vr(1,1,1)-LOGFldCapacity_vr(1,1,1)
    FCD_vr(1,1,1)=LOGFldCapacity_vr(1,1,1)-LOGWiltPoint_vr(1,1,1)
    refreshed=GetSoilHydraulicParameters(1,1,1,material%saturation_tolerance)
    call check_close('soil material change refreshes FC-scaled cutoff', &
      refreshed%resistance%theta_wet,.6_r8*FieldCapacity_vr(1,1,1))
    FieldCapacity_vr(1,1,1)=saved_fc(1);LOGFldCapacity_vr(1,1,1)=saved_fc(2)
    PSD_vr(1,1,1)=saved_fc(3);FCD_vr(1,1,1)=saved_fc(4)
    refreshed=GetSoilHydraulicParameters(1,1,1,material%saturation_tolerance)
    call require('restored FC gives identical cached coefficients', &
      all(refreshed%resistance%coefficient==material%resistance%coefficient))
    work=SoilBulkDensity_vr(1,1,1)
    SoilBulkDensity_vr(1,1,1)=0._r8
    refreshed=GetSoilHydraulicParameters(1,1,1,material%saturation_tolerance)
    call require('pond adapter does not initialize selected-soil resistance',.not.refreshed%resistance%ready)
    SoilBulkDensity_vr(1,1,1)=work
    end associate
    call PrintInfo('end '//subname)
  end subroutine check_resistance_cache
!------------------------------------------------------------------------------------------
  subroutine check_bet_resistance_state()
    ! Description:
    ! Verify both cubic endpoints, monotonicity across temperature, independent
    ! BET inversion/slopes, finite-difference cubic slopes and hyperdry stability.
    type(BETResistanceParameters) :: retention,flat
    real(r8), parameter :: temperatures(3)=[233.15_r8,293.15_r8,323.15_r8]
    real(r8) :: thermal,theta,potential,log_g,log_ratio,x,expected,factor,resistance
    real(r8) :: lo,hi,unused1,unused2,step,previous,temperature,slope_lo,slope_hi
    integer :: j,i
    character(len=*), parameter :: subname='check_bet_resistance_state'

    call PrintInfo('beg '//subname)
    call require('global soil and litter BET defaults',soil_bet_c==18._r8.and.litter_bet_c==13._r8)
    call require('requested material wet cutoffs',SoilBETCutoffFCFraction==.6_r8.and.LitterBETCutoffFCFraction==.5_r8)
    retention=MakeBETResistanceParameters(.1_r8,-2._r8,55._r8,.02_r8,soil_bet_c)
    do j=1,size(temperatures)
      temperature=temperatures(j)
      thermal=RGASC*temperature/18._r8
      call BETResistanceState(retention,.1_r8,temperature,potential,log_g,log_ratio)
      call require('cubic fixed native wet pressure',ABS(potential/(-2._r8)-1._r8)<1.e-12_r8)
      call require('cubic fixed native wet slope',ABS(EXP(log_g)/55._r8-1._r8)<1.e-12_r8)
      call BETResistanceState(retention,retention%theta_bet,temperature,potential,log_g,log_ratio)
      call require('temperature-updated BET joining pressure',ABS(potential/(thermal*LOG(.3_r8))-1._r8)<1.e-12_r8)
      call BETResistanceState(retention,retention%theta_bet*(1._r8-1.e-8_r8),temperature,lo,slope_lo,unused1)
      call BETResistanceState(retention,retention%theta_bet*(1._r8+1.e-8_r8),temperature,hi,slope_hi,unused2)
      call require('cubic BET pressure continuity',ABS(lo/hi-1._r8)<1.e-7_r8)
      call require('cubic BET slope continuity',ABS(slope_lo-slope_hi)<1.e-6_r8)
      previous=-huge(1._r8)
      do i=1,400
        theta=.1_r8*real(i,r8)/400._r8
        call BETResistanceState(retention,theta,temperature,potential,log_g,log_ratio)
        call require('temperature-dependent cubic BET monotonicity',potential>previous.and.ieee_is_finite(log_g))
        previous=potential
      enddo
      theta=.5_r8*(retention%theta_wet+retention%theta_bet)
      step=theta*1.e-5_r8
      call BETResistanceState(retention,theta,temperature,potential,log_g,log_ratio)
      call BETResistanceState(retention,theta-step,temperature,lo,unused1,unused2)
      call BETResistanceState(retention,theta+step,temperature,hi,unused1,unused2)
      call require('cubic analytic derivative against finite difference',ABS((hi-lo)/(2._r8*step)/EXP(log_g)-1._r8)<1.e-7_r8)
      x=.01_r8
      factor=(1._r8-x)*(1._r8+17._r8*x)
      theta=retention%monolayer*18._r8*x/factor
      call BETResistanceState(retention,theta,temperature,potential,log_g,log_ratio)
      call require('BET independent Kelvin inverse',ABS(potential/(thermal*LOG(x))-1._r8)<1.e-12_r8)
      call check_close('BET independent theta/RH',EXP(log_ratio),theta/x)
      expected=thermal*factor**2/(retention%monolayer*18._r8*x*(1._r8+17._r8*x*x))
      call check_close('BET independent analytic derivative',EXP(log_g),expected)
    enddo
    call BETResistanceState(retention,tiny(1._r8),293.15_r8,potential,log_g,log_ratio)
    call require('hyperdry logarithmic factors finite',all(ieee_is_finite([potential,log_g,log_ratio])))
    call check_close('BET finite dry theta/RH',EXP(log_ratio),retention%monolayer*18._r8)
    resistance=PorousSurfaceResistance(tiny(1._r8),.01_r8,1.e-300_r8,0._r8,.01_r8,1._r8,log_g+log_ratio)
    call require('hyperdry liquid resistance finite',ieee_is_finite(resistance).and.resistance>0._r8)
    flat=MakeBETResistanceParameters(.1_r8,-2500._r8,0._r8,.2_r8,18._r8)
    call BETResistanceState(flat,.1_r8,293.15_r8,potential,log_g,log_ratio)
    call require('capped native cutoff remains pressure continuous',ABS(potential/(-2500._r8)-1._r8)<1.e-12_r8)
    call require('capped native cutoff remains slope continuous',EXP(log_g)<1.e-7_r8)
    call PrintInfo('end '//subname)
  end subroutine check_bet_resistance_state
!------------------------------------------------------------------------------------------
  subroutine check_litter_resistance(litter)
    ! Description:
    ! Check transport against an analytic slope and independently formed bulk
    ! diffusivity, including patch geometry, wet films, ice and blocked limits.
    type(LitterHydraulicParameters), intent(in) :: litter ! litter retention and pore parameters
    real(r8), parameter :: vol=.01_r8,area=1._r8,k=1.e-7_r8,dgas=.077_r8,tau=.66_r8,rhov=1.e-5_r8
    real(r8) :: theta,psi,analytic,dv,expected,r,gas_only,icy,thin,liquid_only,source_vapor
    character(len=*), parameter :: subname='check_litter_resistance'

    call PrintInfo('beg '//subname)
    theta=.0334_r8
    psi=LitterMatricPotential(litter,theta)
    analytic=-psi*litter%log_dry_potential_ratio/(litter%log_dry_content_ratio*theta)
    dv=dgas*(litter%physical_porosity-theta)**2*tau/litter%physical_porosity
    expected=(vol/area)/(2._r8*(dv+theta*k*analytic/rhov))
    r=LitterResistance(litter,theta*vol,0._r8,vol,area,k,dgas,tau,rhov)
    call require('analytic parallel litter resistance',abs(r/expected-1._r8)<2.e-7_r8)
    call check_close('litter area-volume scaling',LitterResistance(litter,2._r8*theta*vol,0._r8, &
      2._r8*vol,2._r8*area,k,dgas,tau,rhov),r)
    call check_close('litter covered-area thickness',LitterResistance(litter,theta*vol,0._r8, &
      vol,.5_r8*area,k,dgas,tau,rhov),2._r8*r)
    thin=LitterResistance(litter,theta*vol*1.e-6_r8,0._r8,vol*1.e-6_r8,area,k,dgas,tau,rhov)
    call require('no artificial thin-litter thickness floor',abs(thin/r-1.e-6_r8)<1.e-17_r8)
    gas_only=LitterResistance(litter,0._r8,0._r8,vol,area,k,dgas,tau,rhov)
    call check_close('dry gas-only limit',gas_only,vol/(2._r8*area*dgas*litter%physical_porosity*tau))
    icy=LitterResistance(litter,0._r8,.2_r8*vol,vol,area,k,dgas,tau,rhov)
    expected=vol/(2._r8*area*dgas*(litter%physical_porosity-.2_r8)**2*tau/litter%physical_porosity)
    call check_close('ice reduces actual pore-air diffusion',icy,expected)
    call require('ice-filled pores increase porous resistance',icy>gas_only)
    liquid_only=LitterResistance(litter,theta*vol,0._r8,vol,area,k,0._r8,tau,rhov)
    expected=vol/(2._r8*area*theta*k*analytic/rhov)
    call require('liquid-only limit',abs(liquid_only/expected-1._r8)<2.e-7_r8)
    call require('parallel paths lower resistance',r<LitterResistance(litter,theta*vol,0._r8, &
      vol,area,0._r8,dgas,tau,rhov).and.r<liquid_only)
    call require('blocked pore limit',LitterResistance(litter,theta*vol,0._r8,vol,area, &
      0._r8,0._r8,tau,rhov)==huge(1._r8))
    call require('mobile liquid wet-film limit',LitterResistance(litter,1.01_r8*litter%holding_fraction*vol, &
      0._r8,vol,area,k,dgas,tau,rhov)==0._r8)
    call require('exposed ice limit',LitterResistance(litter,0._r8,litter%physical_porosity*vol, &
      vol,area,k,dgas,tau,rhov)==0._r8)
    call require('no footprint limit',LitterResistance(litter,theta*vol,0._r8,vol,0._r8,k,dgas,tau,rhov)==0._r8)
    call require('no porous volume limit',LitterResistance(litter,0._r8,0._r8,0._r8,area,k,dgas,tau,rhov)==0._r8)
    r=LitterResistance(litter,theta*vol,0._r8,vol,area,1.e-300_r8,dgas,tau,0._r8)
    call require('underflowed source vapor safe',ieee_is_finite(r).and.r>0._r8)
    theta=.001_r8
    psi=LitterMatricPotential(litter,theta)
    source_vapor=vapsat(293.15_r8)*EXP(18._r8*psi/(RGASC*293.15_r8))
    r=LitterResistance(litter,theta*vol,0._r8,vol,area,k,dgas,tau,source_vapor,temperature=293.15_r8,source_potential=psi)
    ! Independent 60-digit BET inversion with analytic native slope at 0.5*FC.
    expected=.000099656100869129922_r8
    call require('independent C13 cubic/BET litter resistance reference',ABS(r/expected-1._r8)<2.e-7_r8)
    call require('BET leaves native litter pressure unchanged',LitterMatricPotential(litter,theta)==psi)
    r=LitterResistance(litter,0._r8,0._r8,vol,area,k,dgas,tau,rhov,temperature=293.15_r8)
    call check_close('BET zero liquid leaves vapor path',r,gas_only)
    r=LitterResistance(litter,theta*vol,.2_r8*vol,vol,area,k,dgas,tau,rhov, &
      temperature=270._r8,source_potential=-2500._r8)
    expected=LitterResistance(litter,theta*vol,.2_r8*vol,vol,area,k,dgas,tau,rhov, &
      temperature=270._r8,source_potential=-1.5_r8)
    call check_close('BET retains ice-equilibrium vapor partition',r,expected)
    call PrintInfo('end '//subname)
  end subroutine check_litter_resistance
!------------------------------------------------------------------------------------------
  subroutine check_litter_independence()
    ! Description:
    ! Check physical retention anchors, inverse/slope consistency and invariance of
    ! the real litter cache/compatibility entry points under both soil selectors.
    type(LitterHydraulicParameters) :: litter,changed_geometry
    real(r8) :: reference_table(3,100,2),reference_entry(2),reference_sat(2)
    real(r8) :: reference_dry,reference_pressure,reference_slope,pressure,content,slope,analytic
    real(r8) :: held_water,held_ice,theta,reference_inverse,resistance,reference_resistance
    real(r8) :: wet_pressure,wet_log_slope,wet_log_ratio,native_wet_slope
    integer :: scheme,policy,k
    character(len=*), parameter :: subname='check_litter_independence'

    call PrintInfo('beg '//subname)
    POROS_vr(0,1,1)=.97_r8
    call UpdateLitterHydraulicLimits(1,1,.1336_r8)
    SatHydroCondVert_vr(0,1,1)=.98_r8
    litter=GetLitterHydraulicParameters(1,1)
    call require('litter resistance cache initialized',litter%resistance%ready)
    call check_close('litter cutoff is half field capacity',litter%resistance%theta_wet,.5_r8*litter%field_capacity)
    theta=litter%resistance%theta_wet
    call BETResistanceState(litter%resistance,theta,273.15_r8,wet_pressure,wet_log_slope,wet_log_ratio)
    native_wet_slope=LitterMatricPotentialDerivative(litter,theta)
    call require('litter cutoff pressure matches native at different T', &
      ABS(wet_pressure/LitterMatricPotential(litter,theta)-1._r8)<1.e-12_r8)
    call require('litter cutoff slope matches native at different T', &
      ABS(EXP(wet_log_slope)/native_wet_slope-1._r8)<1.e-12_r8)
    FieldCapacity_vr(0,1,1)=.9_r8*litter%field_capacity
    changed_geometry=GetLitterHydraulicParameters(1,1)
    call check_close('litter field capacity refreshes cutoff', &
      changed_geometry%resistance%theta_wet,.5_r8*FieldCapacity_vr(0,1,1))
    FieldCapacity_vr(0,1,1)=litter%field_capacity
    changed_geometry=GetLitterHydraulicParameters(1,1)
    call require('litter restored FC gives identical cached coefficients', &
      all(changed_geometry%resistance%coefficient==litter%resistance%coefficient))
    litter_bet_c=17._r8
    changed_geometry=GetLitterHydraulicParameters(1,1)
    call require('litter C override refreshes cache',changed_geometry%resistance%adsorption_c==17._r8)
    litter_bet_c=13._r8
    call UpdateLitterHydraulicLimits(1,1,.2_r8)
    changed_geometry=GetLitterHydraulicParameters(1,1)
    call check_close('litter capacity refreshes cutoff',changed_geometry%resistance%theta_wet,.05_r8)
    call require('litter capacity refreshes monolayer', &
      ABS(changed_geometry%resistance%monolayer/litter%resistance%monolayer-.2_r8/.1336_r8)<1.e-12_r8)
    call UpdateLitterHydraulicLimits(1,1,.1336_r8)
    changed_geometry=GetLitterHydraulicParameters(1,1)
    call require('litter restored cache deterministic',all(changed_geometry%resistance%coefficient==litter%resistance%coefficient))
    call check_litter_resistance(litter)
    call check_close('litter holding fraction',litter%holding_fraction,.1336_r8)
    call check_close('litter geometric porosity',litter%physical_porosity,.97_r8)
    call check_close('litter saturation anchor',LitterMatricPotential(litter,.1336_r8),-.0005_r8)
    call check_close('litter field anchor',LitterMatricPotential(litter,.0668_r8),-.033_r8)
    call check_close('litter wilt anchor',LitterMatricPotential(litter,.0167_r8),-1.5_r8)
    call check_close('litter zero water',LitterMatricPotential(litter,0._r8),-2500._r8)
    call check_close('litter dry inverse',LitterWaterContent(litter,-25000._r8),LitterDryContent(litter))
    call check_close('litter wet inverse',LitterWaterContent(litter,0._r8),.1336_r8)
    call require('empty litter class safe',LitterMoistureClass(0._r8,0._r8)==100)
    call require('dry litter class',LitterMoistureClass(0._r8,.1336_r8)==100)
    call require('wet litter class',LitterMoistureClass(.2_r8,.1336_r8)==1)
    changed_geometry=litter
    changed_geometry%physical_porosity=.4_r8
    do k=1,30
      theta=.1336_r8*real(k,r8)/30._r8
      pressure=LitterMatricPotential(litter,theta)
      call check_close('litter forward/inverse',LitterWaterContent(litter,pressure),theta)
      call require('litter curve independent of geometric porosity', &
        LitterMatricPotential(changed_geometry,theta)==pressure)
      call require('litter class preserves existing indexing', &
        LitterMoistureClass(theta,.1336_r8)==getMoistK(theta,.1336_r8))
    enddo
    theta=.0334_r8
    pressure=LitterMatricPotential(litter,theta)
    analytic=-pressure*litter%log_dry_potential_ratio/(litter%log_dry_content_ratio*theta)
    slope=LitterMatricPotentialDerivative(litter,theta)
    call require('litter dry-branch derivative',abs(slope/analytic-1._r8)<2.e-7_r8)
    theta=.1_r8
    pressure=LitterMatricPotential(litter,theta)
    analytic=-pressure*litter%log_wet_potential_ratio/(litter%log_wet_content_ratio*theta)
    slope=LitterMatricPotentialDerivative(litter,theta)
    call require('litter wet-branch derivative',abs(slope/analytic-1._r8)<2.e-7_r8)
    call require('litter dry plateau derivative',LitterMatricPotentialDerivative(litter,0._r8)==0._r8)
    call require('litter wet plateau derivative',LitterMatricPotentialDerivative(litter,.2_r8)==0._r8)
    ! Supply real ordinary-soil inputs so switching also exercises a soil refresh.
    allocate(PSIAtFldCapacity_col(1,1),PSIAtWiltPoint_col(1,1))
    PSIAtFldCapacity_col=-.033_r8;PSIAtWiltPoint_col=-1.5_r8
    allocate(VGThetaR_vr(1:1,1,1),VGAlpha_vr(1:1,1,1),VGN_vr(1:1,1,1),VGL_vr(1:1,1,1))
    VGThetaR_vr=.05_r8;VGAlpha_vr=100._r8;VGN_vr=1.6_r8;VGL_vr=.5_r8
    held_water=VLWatMicP_vr(0,1,1);held_ice=VLiceMicP_vr(0,1,1)
    do scheme=1,2
      if(scheme==1)then
        soil_hydraulic_scheme='ecosim'
      else
        soil_hydraulic_scheme='van_genuchten'
      endif
      call InitializeSoilHydraulicScheme()
      call require('soil selector actually switched',GetSoilHydraulicSchemeName()==soil_hydraulic_scheme)
      call SoilHydroProperty(1,1,1,1,1)
      do policy=1,2
        PSISoilAirEntry_vr(0,1,1)=-.0015_r8;ThetaSat_vr(0,1,1)=0._r8
        use_ecosys_drainage=policy==2
        call LitterHydroproperty(1,1,use_ecosys_drainage)
        ! The generic entry point must delegate before ordinary-soil estimators.
        call SoilHydroProperty(1,1,0,1,1)
        call check_close('generic litter holding fraction retained',fLitRPore4WatHold_col(1,1),.1336_r8)
        call check_close('generic litter curvature retained',SRP_vr(0,1,1),1._r8)
        if(scheme==1)then
          reference_table(:,:,policy)=HydroCond_3D(:,:,0,1,1)
          reference_entry(policy)=PSISoilAirEntry_vr(0,1,1)
          reference_sat(policy)=ThetaSat_vr(0,1,1)
        else
          call require('litter K invariant under soil scheme',all(HydroCond_3D(:,:,0,1,1)==reference_table(:,:,policy)))
          call require('litter air entry invariant',PSISoilAirEntry_vr(0,1,1)==reference_entry(policy))
          call require('litter saturation threshold invariant',ThetaSat_vr(0,1,1)==reference_sat(policy))
        endif
      enddo
      theta=.1_r8
      call ComputePsiMCM(0,1,1,theta,pressure)
      call ComputeSoilAirDryContent(0,1,1,content)
      litter=GetLitterHydraulicParameters(1,1)
      slope=LitterMatricPotentialDerivative(litter,.1_r8)
      resistance=LitterResistance(litter,.001_r8,0._r8,.01_r8,1._r8, &
        HydroCond_3D(3,LitterMoistureClass(.1_r8,litter%holding_fraction),0,1,1),.077_r8,.66_r8,1.e-5_r8)
      if(scheme==1)then
        reference_pressure=pressure;reference_dry=content;reference_slope=slope
        reference_inverse=LitterWaterContent(litter,-.05_r8)
        reference_resistance=resistance
      else
        call require('litter pressure dispatch invariant',pressure==reference_pressure)
        call require('litter dry content dispatch invariant',content==reference_dry)
        call require('litter slope invariant',slope==reference_slope)
        call require('litter resistance invariant',resistance==reference_resistance)
        call require('litter inverse invariant',LitterWaterContent(litter,-.05_r8)==reference_inverse)
      endif
      call require('litter liquid preserved',VLWatMicP_vr(0,1,1)==held_water)
      call require('litter ice preserved',VLiceMicP_vr(0,1,1)==held_ice)
    enddo
    call require('litter conductivity finite',all(ieee_is_finite(HydroCond_3D(:,:,0,1,1))))
    do k=2,100
      call require('litter conductivity monotonic',HydroCond_3D(3,k,0,1,1)<=HydroCond_3D(3,k-1,0,1,1))
    enddo
    call PrintInfo('end '//subname)
  end subroutine check_litter_independence
!------------------------------------------------------------------------------------------
  subroutine check_soil_resistance(scheme,context)
    ! Description:
    ! Verify the paper's parallel transport against an analytic retention slope,
    ! including the half-cell length, native units, dry-end slope extension,
    ! blocked paths and numeric limits.
    class(SoilHydraulicModel), intent(in) :: scheme      ! selected constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: context ! layer retention parameters
    real(r8), parameter :: dz=.02_r8,k=1.e-6_r8,dgas=.004_r8,bunsen=60000._r8
    real(r8) :: theta,work,pressure,slope,combined,rg,rw,expected,dry,boundary_resistance,source_bunsen
    character(len=*), parameter :: subname='check_soil_resistance'

    call PrintInfo('beg '//subname)
    theta=.20_r8;work=theta
    pressure=scheme%matric_potential(context,work)
    slope=-pressure*context%log_dry_potential_ratio/(context%log_dry_content_ratio*theta)
    expected=dz/(2._r8*(dgas+bunsen*theta*k*slope))
    combined=scheme%soil_resistance(context,theta,dz,k,dgas,bunsen)
    call require('paper soil resistance with analytic slope',ieee_is_finite(combined).and. &
      ABS(combined-expected)<2.e-7_r8*expected)
    call require('soil resistance preserves liquid content',theta==.20_r8)
    rg=scheme%soil_resistance(context,theta,dz,0._r8,dgas,bunsen)
    rw=scheme%soil_resistance(context,theta,dz,k,0._r8,bunsen)
    call check_close('vapor-only half-cell resistance',rg,dz/(2._r8*dgas))
    call check_close('parallel liquid and vapor paths',combined,1._r8/(1._r8/rg+1._r8/rw))
    call require('parallel resistance below either path',combined<MIN(rg,rw))
    call check_close('soil resistance thickness scaling', &
      scheme%soil_resistance(context,theta,2._r8*dz,k,dgas,bunsen),2._r8*combined)
    call check_close('liquid resistance conductivity scaling', &
      scheme%soil_resistance(context,theta,dz,2._r8*k,0._r8,bunsen),.5_r8*rw)
    call check_close('liquid resistance partition scaling', &
      scheme%soil_resistance(context,theta,dz,k,0._r8,2._r8*bunsen),.5_r8*rw)
    call check_close('zero B leaves vapor path',scheme%soil_resistance(context,theta,dz,k,dgas,0._r8),rg)
    call check_close('zero liquid leaves vapor path',scheme%soil_resistance(context,0._r8,dz,k,dgas,bunsen),rg)
    call require('zero thickness has no internal resistance', &
      scheme%soil_resistance(context,theta,0._r8,k,dgas,bunsen)==0._r8)
    call require('blocked soil transport', &
      scheme%soil_resistance(context,theta,dz,0._r8,0._r8,bunsen)==huge(1._r8))
    dry=scheme%dry_content(context)
    boundary_resistance=scheme%soil_resistance(context,dry,dz,k,0._r8,bunsen)
    call require('finite positive dry boundary resistance', &
      ieee_is_finite(boundary_resistance).and.boundary_resistance>0._r8.and.boundary_resistance<huge(1._r8))
    theta=.5_r8*dry;work=theta
    call check_close('extended slope retains actual liquid content', &
      scheme%soil_resistance(context,theta,dz,k,0._r8,bunsen),2._r8*boundary_resistance)
    expected=1._r8/(1._r8/rg+1._r8/(2._r8*boundary_resistance))
    call check_close('extended slope combines both transport paths', &
      scheme%soil_resistance(context,theta,dz,k,dgas,bunsen),expected)
    call check_close('dry extension preserves capped pressure',scheme%matric_potential(context,work),context%dry_potential)
    call require('dry extension preserves input content',theta==.5_r8*dry.and.work==theta)
    call require('dry retention derivative remains clipped',scheme%matric_potential_derivative(context,theta)==0._r8)
    theta=dry*(1._r8-1.e-8_r8)
    combined=scheme%soil_resistance(context,theta,dz,k,dgas,bunsen)
    theta=dry*(1._r8+1.e-8_r8)
    expected=scheme%soil_resistance(context,theta,dz,k,dgas,bunsen)
    call require('soil resistance continuous across dry cutoff',ABS(combined-expected)<1.e-6_r8*expected)
    theta=.5_r8*dry;work=theta
    pressure=scheme%matric_potential(context,work)
    source_bunsen=1._r8/(vapsat(293.15_r8)*EXP(18._r8*pressure/(RGASC*293.15_r8)))
    combined=scheme%soil_resistance(context,theta,dz,k,dgas,source_bunsen, &
      temperature=293.15_r8,source_potential=pressure)
    expected=.000027726926060740094_r8
    call require('independent C18 cubic/BET soil resistance reference',ABS(combined/expected-1._r8)<2.e-7_r8)
    call require('BET leaves native soil pressure unchanged',scheme%matric_potential(context,work)==pressure)
    combined=scheme%soil_resistance(context,theta,dz,k,dgas,bunsen, &
      temperature=270._r8,source_potential=-2500._r8,ice_volume=1.e-5_r8)
    expected=scheme%soil_resistance(context,theta,dz,k,dgas,bunsen, &
      temperature=270._r8,source_potential=-1.5_r8,ice_volume=1.e-5_r8)
    call check_close('BET soil retains ice-equilibrium partition',combined,expected)
    call require('large parallel conductance avoids product overflow', &
      scheme%soil_resistance(context,theta,dz,huge(1._r8)/4._r8,dgas,huge(1._r8)/4._r8)==0._r8)
    expected=scheme%soil_resistance(context,theta,huge(1._r8)/8._r8,0._r8,tiny(1._r8),bunsen)
    call require('large soil resistance remains finite',ieee_is_finite(expected).and.expected==huge(1._r8))
    call PrintInfo('end '//subname)
  end subroutine check_soil_resistance
!------------------------------------------------------------------------------------------
  subroutine check_retention_derivative(scheme,context)
    ! Description:
    ! Compare numerical slopes with analytic branch derivatives and exercise
    ! dry/wet endpoints, pressure clipping and the inout saturation adjustment.
    class(SoilHydraulicModel), intent(in) :: scheme          ! selected retention implementation [-]
    type(SoilHydraulicParameters), intent(in) :: context     ! current layer's retention parameters
    real(r8), parameter :: samples(6)=[.10_r8,.20_r8,.29_r8,.31_r8,.37_r8,.44_r8]
    type(SoilHydraulicParameters) :: snapped_context
    real(r8) :: sample,work,pressure,slope,expected,x,dry
    integer :: i
    character(len=*), parameter :: subname='check_retention_derivative'

    call PrintInfo('beg '//subname)
    do i=1,size(samples)
      sample=samples(i)
      slope=scheme%matric_potential_derivative(context,sample)
      call require('derivative leaves water content unchanged',sample==samples(i))
      work=sample
      pressure=scheme%matric_potential(context,work)
      if(sample<context%field_capacity)then
        expected=-pressure*context%log_dry_potential_ratio/(context%log_dry_content_ratio*sample)
      else
        x=LOG(context%porosity/sample)/context%log_wet_content_ratio
        expected=-pressure*context%curvature*context%log_wet_potential_ratio &
          *x**(context%curvature-1._r8)/(context%log_wet_content_ratio*sample)
      endif
      call require('analytic EcoSIM retention slope',ieee_is_finite(slope).and. &
        ABS(slope-expected)<2.e-7_r8*ABS(expected))
      call check_close('derivative preserves pressure evaluation',scheme%matric_potential(context,work),pressure)
    enddo
    dry=scheme%dry_content(context)
    call require('dry plateau slope',scheme%matric_potential_derivative(context,.5_r8*dry)==0._r8)
    call require('zero-water slope',scheme%matric_potential_derivative(context,0._r8)==0._r8)
    call require('oversaturated plateau slope', &
      scheme%matric_potential_derivative(context,1.1_r8*context%porosity)==0._r8)
    slope=scheme%matric_potential_derivative(context,dry)
    call require('finite one-sided dry endpoint',ieee_is_finite(slope).and.slope>0._r8)
    sample=context%porosity
    slope=scheme%matric_potential_derivative(context,sample)
    call require('finite one-sided saturated endpoint',ieee_is_finite(slope).and.slope>0._r8)
    call require('saturated derivative preserves content',sample==context%porosity)

    ! Both samples snap to porosity. A denominator made from the modified
    ! arguments would be zero even though the requested stencil has width.
    snapped_context=context
    snapped_context%saturation_tolerance=.01_r8
    sample=context%porosity-.005_r8
    slope=scheme%matric_potential_derivative(snapped_context,sample)
    call require('snapped plateau has finite zero slope',ieee_is_finite(slope).and.slope==0._r8)
    call require('snapped derivative preserves content',sample==context%porosity-.005_r8)
    call PrintInfo('end '//subname)
  end subroutine check_retention_derivative
!------------------------------------------------------------------------------------------
  subroutine initialize_state()
  ! Description:
  ! Allocate a minimal one-column state for public constitutive/property routines.
  character(len=*), parameter :: subname='initialize_state'
  call PrintInfo('beg '//subname)
  lverb=.false.;ZERO=1.e-16_r8;ZERO2=1.e-8_r8
  allocate(ZEROS(1,1),ZEROS2(1,1));ZEROS=ZERO;ZEROS2=ZERO2
  allocate(POROS_vr(0:1,1,1),FieldCapacity_vr(0:1,1,1),WiltPoint_vr(0:1,1,1))
  allocate(LOGPOROS_vr(0:1,1,1),LOGFldCapacity_vr(0:1,1,1),LOGWiltPoint_vr(0:1,1,1))
  allocate(PSD_vr(0:1,1,1),FCD_vr(0:1,1,1),SRP_vr(0:1,1,1),SoilWatAirDry_vr(0:1,1,1))
  allocate(LOGPSIFLD_col(1,1),LOGPSIAtSat(1,1),LOGPSIMXD_col(1,1),LOGPSIMND_col(1,1))
  allocate(PSI_at_Sautration_vr(0:1,1,1),PSISoilAirEntry_vr(0:1,1,1),ThetaSat_vr(0:1,1,1))
  allocate(SatHydroCondVert_vr(0:1,1,1),SatHydroCondHrzn_vr(0:1,1,1),HydroCond_3D(3,100,0:1,1,1))
  allocate(CSoilOrgM_vr(3,0:1,1,1),ISOIL_vr(4,0:1,1,1),SoilBulkDensity_vr(0:1,1,1))
  allocate(VLsoiAirP_vr(0:1,1,1),ThetaAir_vr(0:1,1,1),VGeomLayer_vr(0:1,1,1),VLSoilMicP_vr(0:1,1,1))
  allocate(VLMicP_vr(0:1,1,1),VLWatMicP_vr(0:1,1,1),VLiceMicP_vr(0:1,1,1),VLMacroPores_vr(0:1,1,1))
  allocate(VGeomLayert0_vr(0:1,1,1));VGeomLayert0_vr=1._r8
  allocate(VLWatMacP_vr(0:1,1,1),VLiceMacP_vr(0:1,1,1),TCS_vr(0:1,1,1))
  allocate(MacPoreRadius_vr(0:1,1,1),MacPoreNumbers_vr(0:1,1,1),PathLenMacPore_vr(0:1,1,1),HydroCondMacP_vr(0:1,1,1))
  POROS_vr=.45_r8;FieldCapacity_vr=.30_r8;WiltPoint_vr=.15_r8
  LOGPSIFLD_col=log(.033_r8);LOGPSIAtSat=log(.0005_r8)
  LOGPSIMXD_col=log(.033_r8/.0005_r8);LOGPSIMND_col=log(1.5_r8/.033_r8)
  PSI_at_Sautration_vr=-.0005_r8;SoilBulkDensity_vr=1.3_r8;ISOIL_vr=isoi_set
  SatHydroCondVert_vr=1._r8;SatHydroCondHrzn_vr=2._r8
  VGeomLayer_vr=1._r8;VLSoilMicP_vr=1._r8;VLMicP_vr=.45_r8;VLWatMicP_vr=.27_r8
  VLiceMicP_vr=0._r8;VLMacroPores_vr=0._r8;VLWatMacP_vr=0._r8;VLiceMacP_vr=0._r8;TCS_vr=10._r8
  allocate(fLitRPore4WatHold_col(1,1),VLSoilMicPMass_vr(0:1,1,1))
  VLSoilMicPMass_vr=1._r8
  CSoilOrgM_vr=0._r8
  call PrintInfo('end '//subname)
  end subroutine initialize_state
!------------------------------------------------------------------------------------------
  subroutine check_close(label,actual,expected)
  ! Description:
  ! Check frozen constitutive samples with a fixed roundoff allowance.
  character(len=*), intent(in) :: label ! assertion description [-]
  real(r8), intent(in) :: actual       ! evaluated scalar [units of the assertion]
  real(r8), intent(in) :: expected     ! reference scalar [same units as actual]
  character(len=*), parameter :: subname='check_close'
  call PrintInfo('beg '//subname)
  if(.not.ieee_is_finite(actual) .or. abs(actual-expected)>2.e-13_r8*max(abs(expected),tiny(1._r8)))then
    write(iulog,*) trim(label),': actual=',actual,' expected=',expected
    error stop 1
  endif
  call PrintInfo('end '//subname)
  end subroutine check_close
!------------------------------------------------------------------------------------------
  subroutine require(label,condition)
  ! Description:
  ! Fail a qualitative domain, clipping or ordering assertion.
  character(len=*), intent(in) :: label ! assertion description [-]
  logical, intent(in) :: condition      ! required condition [-]
  character(len=*), parameter :: subname='require'
  call PrintInfo('beg '//subname)
  if(.not.condition)then
    write(iulog,*) 'Failed: ',trim(label)
    error stop 1
  endif
  call PrintInfo('end '//subname)
  end subroutine require
end program SoilHydraulicTest
