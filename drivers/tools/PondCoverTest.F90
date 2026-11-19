program PondCoverTest
  ! Description:
  ! Check that an ATS-owned pond reduces the exposed-soil cover fraction in the
  ! legacy (non-shared) ATS-coupled surface layer: no change without a pond or
  ! with the option off, continuity as the pond goes to zero, a monotonic
  ! decrease with pond volume, full cover at the surface-retention scale, and no
  ! effect where litter is not heat-wise significant.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: iulog
  use DebugToolMod, only: PrintInfo
  use GridConsts, only: JX,JY,JX0,JY0,JZ,JH,JV,JD
  use GridDataType, only: InitGridData,AREA_3D
  use ElmIDMod, only: ielmc
  use EcoSIMCtrlDataType, only: ats_shared_surface_layer
  use HydroThermData
  use SurfLitterDataType
  use SurfPhysData
  use SoilWaterDataType
  use SoilHeatDataType
  use SOMDataType, only: SoilOrgM_vr
  use SurfSoilDataType, only: FracSurfBareSoil_col,VLWatHeldCapSurf_col
  use SurfLitterPhysMod, only: PartLitSoilFractionM
  implicit none
  real(r8), parameter :: LitterC=60._r8           ! litter organic carbon [g C/cell]
  real(r8), parameter :: CellArea=1._r8           ! cell area [m2]
  real(r8), parameter :: RetentionScale=0.01_r8   ! surface water retention capacity [m3/cell]
  real(r8) :: BareNoPond,BarePrev,PondFraction
  integer :: i

  JX=1;JY=1;JX0=1;JY0=1;JZ=1;JH=1;JV=1;JD=2
  call InitHydroThermData()
  call InitSoilWater()
  call InitSoilHeatData()
  call InitSurfLitter(3)
  call InitSurfPhysData()
  call InitGridData()
  allocate(SoilOrgM_vr(3,0:1,1,1),FracSurfBareSoil_col(1,1),VLWatHeldCapSurf_col(1,1))
  ats_shared_surface_layer=.false.

  ! No pond: the pond-aware option must not change the cover, on or off.
  call reset_state()
  ats_pond_aware=.false.
  call PartLitSoilFractionM(1,1,1,1)
  BareNoPond=FracSurfBareSoil_col(1,1)
  call check_close('bare fraction with option off and no pond',BareNoPond, &
    exp(-0.8e-2_r8*LitterC/CellArea))
  ats_pond_aware=.true.
  call PartLitSoilFractionM(1,1,1,1)
  call check_close('option on, zero pond',FracSurfBareSoil_col(1,1),BareNoPond)

  ! Option off ignores any pond.
  ats_pond_aware=.false.
  XVLPondATS_col=0.5_r8*RetentionScale
  call PartLitSoilFractionM(1,1,1,1)
  call check_close('option off ignores pond',FracSurfBareSoil_col(1,1),BareNoPond)

  ! Option on: the covered fraction is removed from the exposed soil.
  ats_pond_aware=.true.
  do i=1,3
    PondFraction=0.25_r8*real(i,r8)
    XVLPondATS_col=PondFraction*RetentionScale
    call PartLitSoilFractionM(1,1,1,1)
    call check_close('bare fraction scales with uncovered area',FracSurfBareSoil_col(1,1), &
      BareNoPond*(1._r8-PondFraction))
    call check_close('litter and bare fractions sum to one', &
      FracSurfByLitR_col(1,1)+FracSurfBareSoil_col(1,1),1._r8)
  enddo

  ! Full cover at and beyond the retention scale.
  XVLPondATS_col=RetentionScale
  call PartLitSoilFractionM(1,1,1,1)
  call check_close('full cover at the retention scale',FracSurfBareSoil_col(1,1),0._r8)
  XVLPondATS_col=10._r8*RetentionScale
  call PartLitSoilFractionM(1,1,1,1)
  call check_close('cover saturates beyond the retention scale',FracSurfBareSoil_col(1,1),0._r8)
  call check_close('all area is litter pathway when fully covered',FracSurfByLitR_col(1,1),1._r8)

  ! Monotonic decrease and continuity with respect to the pond volume.
  BarePrev=BareNoPond
  do i=1,20
    XVLPondATS_col=RetentionScale*real(i,r8)/20._r8
    call PartLitSoilFractionM(1,1,1,1)
    if(FracSurfBareSoil_col(1,1)>BarePrev+1.e-14_r8)then
      write(iulog,*)'bare fraction increased with pond volume',i,FracSurfBareSoil_col(1,1),BarePrev
      error stop 'pond cover not monotonic'
    endif
    BarePrev=FracSurfBareSoil_col(1,1)
  enddo
  XVLPondATS_col=1.e-12_r8*RetentionScale
  call PartLitSoilFractionM(1,1,1,1)
  call check_close('continuous as the pond vanishes',FracSurfBareSoil_col(1,1),BareNoPond)

  ! No heat-wise significant litter: the pond has no litter pathway to join.
  call reset_state()
  VHeatCapacity1_vr(0,1,1)=0.5_r8*VHeatCapLitRMin_col(1,1)
  XVLPondATS_col=RetentionScale
  ats_pond_aware=.true.
  call PartLitSoilFractionM(1,1,1,1)
  call check_close('no litter: pond leaves the bare fraction alone',FracSurfBareSoil_col(1,1),1._r8)
  write(iulog,*)'PASS: ATS pond reduces the exposed-soil cover consistently'
contains
!------------------------------------------------------------------------------------------
  subroutine reset_state()
    ! Description:
    ! Set one column with significant litter, no mobile litter water and no pond.
    character(len=*), parameter :: subname='reset_state'
    call PrintInfo('beg '//subname)
    SoilOrgM_vr=0._r8;SoilOrgM_vr(ielmc,0,1,1)=LitterC
    AREA_3D=CellArea
    VHeatCapLitRMin_col=1.e-6_r8
    VHeatCapacity1_vr=0._r8;VHeatCapacity1_vr(0,1,1)=1.e-3_r8
    XVLMobileWaterLitR_col=0._r8;XVLPondATS_col=0._r8
    VLWatHeldCapSurf_col=RetentionScale
    call PrintInfo('end '//subname)
  end subroutine reset_state
!------------------------------------------------------------------------------------------
  subroutine check_close(Label,Actual,Expected)
    ! Description:
    ! Require a finite value equal to the reference within round-off.
    character(len=*), intent(in) :: Label ! invariant being checked [-]
    real(r8), intent(in) :: Actual ! computed cover fraction [-]
    real(r8), intent(in) :: Expected ! reference cover fraction [-]
    character(len=*), parameter :: subname='check_close'
    call PrintInfo('beg '//subname)
    if(.not.ieee_is_finite(Actual).or.abs(Actual-Expected)>1.e-12_r8)then
      write(iulog,*)Label,Actual,Expected
      error stop 'pond cover fraction is wrong'
    endif
    call PrintInfo('end '//subname)
  end subroutine check_close
end program PondCoverTest
