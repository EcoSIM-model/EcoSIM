program SoilHydraulicTest
  ! Description:
  ! Characterize the existing hydraulic scheme before introducing alternatives.
  ! Frozen pressure/conductivity samples come from pre-refactor GNU 13.4.0 output;
  ! independent bounds, monotonicity and domain checks exercise the public API.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: iulog
  use DebugToolMod, only: PrintInfo
  use EcoSIMCtrlMod, only: lverb,use_ecosys_drainage
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
  use SoilHydroParaMod, only: SoilHydroProperty,LitterHydroproperty
  use SoilPhysParaMod, only: ComputePsiMCM,ComputeSoilAirDryContent, &
    ComputeKsatReferenceContent,getMoistK
  use HydrologyDiagMod, only: ComputeWaterTableEquilContent
  use SoilHydraulicModelMod, only: SoilHydraulicModel,SoilHydraulicParameters
  use EcoSIMSoilHydraulicMod, only: EcoSIMSoilHydraulicModel
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
  POROS0_col=.8_r8;POROS_vr(0,1,1)=.8_r8
  FieldCapacity_vr(0,1,1)=.4_r8;WiltPoint_vr(0,1,1)=.1_r8
  LOGPOROS_vr(0,1,1)=log(.8_r8);LOGFldCapacity_vr(0,1,1)=log(.4_r8)
  LOGWiltPoint_vr(0,1,1)=log(.1_r8)
  PSD_vr(0,1,1)=log(.8_r8)-log(.4_r8);FCD_vr(0,1,1)=log(.4_r8)-log(.1_r8)
  SRP_vr(0,1,1)=1._r8;SatHydroCondVert_vr(0,1,1)=.98_r8
  do mode=1,2
    use_ecosys_drainage=mode==2
    PSISoilAirEntry_vr(0,1,1)=-.0015_r8;ThetaSat_vr(0,1,1)=0._r8
    call LitterHydroproperty(1,1)
    call check_close('litter air-entry policy',PSISoilAirEntry_vr(0,1,1),litter_entry(mode))
    call check_close('litter ThetaSat policy',ThetaSat_vr(0,1,1),litter_theta_sat(mode))
    call check_close('litter curvature retained',SRP_vr(0,1,1),1._r8)
    call require('no horizontal litter transport',all(HydroCond_3D(1:2,:,0,1,1)==0._r8))
    call check_close('saturated litter conductivity',HydroCond_3D(3,1,0,1,1),.98_r8)
  enddo
  call FinalizeSoilHydraulicScheme()
  write(*,*) 'SoilHydraulicTest passed'
contains
!------------------------------------------------------------------------------------------
  subroutine check_soil_resistance(scheme,context)
    ! Description:
    ! Verify the paper's parallel transport against an analytic retention slope,
    ! including the half-cell length, native units, blocked paths and numeric limits.
    class(SoilHydraulicModel), intent(in) :: scheme      ! selected constitutive implementation [-]
    type(SoilHydraulicParameters), intent(in) :: context ! layer retention parameters
    real(r8), parameter :: dz=.02_r8,k=1.e-6_r8,dgas=.004_r8,bunsen=60000._r8
    real(r8) :: theta,work,pressure,slope,combined,rg,rw,expected,dry
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
    dry=.5_r8*scheme%dry_content(context)
    call check_close('clipped dry slope leaves vapor path',scheme%soil_resistance(context,dry,dz,k,dgas,bunsen),rg)
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
  allocate(POROS0_col(1,1),VLSoilMicPMass_vr(0:1,1,1))
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
