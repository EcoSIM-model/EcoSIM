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

  call initialize_state()
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
  call ComputeKsatReferenceContent(1,1,1,content)
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
  write(*,*) 'SoilHydraulicTest passed'
contains
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
