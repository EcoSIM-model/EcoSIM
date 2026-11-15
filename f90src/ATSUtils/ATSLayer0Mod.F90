module ATSLayer0Mod
  ! Description:
  ! Share native litter retention parameters with ATS without advancing water or energy.
  use iso_c_binding, only: c_double, c_int
  use data_kind_mod, only: r8=>DAT_KIND_R8
  use DebugToolMod, only: PrintInfo
  implicit none
  integer, parameter :: n_layer0_properties=14
contains
!------------------------------------------------------------------------------------------
  subroutine ecosim_layer0_material(column, properties, resistance) bind(C)
    ! Description:
    ! Export the native layer-zero retention curve and dry organic heat capacity.
    use SoilPhysDataType
    use SoilPropertyDataType
    use SurfLitterDataType
    use SharedDataMod, only: column_area,a_TCOND,a_dz
    use SurfLitterPhysMod, only: CalcLitRThermConductivity
    use HydroThermData, only: TKSoil1_vr,DLYRR_COL
    use SoilWaterDataType, only: PSI_at_Sautration_vr
    use SoilHeatDataType
    use GridDataType
    use EcoSimConst
    use SOMDataType
    use MiniFuncMod, only: gOC_to_m3_OM
    implicit none
    integer(c_int), value, intent(in) :: column ! zero-based local surface column
    real(c_double), intent(out) :: properties(n_layer0_properties) ! per-area native material parameters
    real(c_double), intent(out) :: resistance ! effective extra surface resistance [m2 K MW-1]
    integer :: ny
    real(r8) :: area,tc,dt,r0
    character(len=*), parameter :: subname='ecosim_layer0_material'
    call PrintInfo('beg '//subname)
    ny=column+1
    area=column_area(ny)
    properties=0._r8; resistance=0._r8
    if(area>0._r8)then
      properties(1)=VWatLitRHoldCapcity_col(ny,1)/area ! retained capacity [m H2O]
      properties(2)=cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,0,ny,1))/area ! dry capacity [MJ m-2 K-1]
      properties(3)=VLitR_col(ny,1)/area ! litter volume per area [m]
      if(DLYRR_COL(ny,1)>0._r8)then
        dt=abs(TKSoil1_vr(0,ny,1)-TKSoil1_vr(1,ny,1))*1.e-6_r8
        call CalcLitRThermConductivity(ny,1,dt,tc)
        if(tc>0._r8.and.a_TCOND(1,ny)>0._r8)then
          r0=1800._r8*DLYRR_COL(ny,1)/tc
          ! ATS applies covered area to this path; do not mix bare soil here.
          resistance=r0
        endif
      endif

      if(properties(1)>0._r8)then
      properties(4)=LOGFldCapacity_vr(0,ny,1)
      properties(5)=FCD_vr(0,ny,1)
      properties(6)=LOGPSIFLD_col(ny,1)
      properties(7)=LOGPSIMND_col(ny,1)
      properties(8)=LOGPOROS_vr(0,ny,1)
      properties(9)=PSD_vr(0,ny,1)
      properties(10)=SRP_vr(0,ny,1)
      properties(11)=LOGPSIAtSat(ny,1)
      properties(12)=LOGPSIMXD_col(ny,1)
      properties(13)=PSIHY ! dry-limit water potential [MPa]
      properties(14)=-exp(LOGPSIAtSat(ny,1)) ! saturation potential [MPa]
      endif
    endif
    call PrintInfo('end '//subname)
  end subroutine ecosim_layer0_material
!------------------------------------------------------------------------------------------
  subroutine InstallSharedLayer0State(ncolumns)
    ! Description:
    ! Install ATS conserved layer-zero water, ice and temperature as native mirrors.
    use SharedDataMod, only: column_area,a_WAT0,a_ICE0,a_TSURF,a_PRES0,a_POND0
    use SoilWaterDataType
    use SoilPropertyDataType
    use SoilHeatDataType
    use SoilPhysDataType
    use SurfPhysData, only: XVLMobileWaterLitR_col
    use SurfLitterPhysMod, only: PartLitSoilFractionM
    use SurfLitterDataType
    use GridDataType
    use EcoSimConst
    use SOMDataType
    use MiniFuncMod, only: gOC_to_m3_OM
    use HydroThermData, only: TKSoil1_vr,DLYRR_COL,VHeatCapacity1_vr,VLWatMicP1_vr,VLiceMicP1_vr, &
      FracSoiPAsWat_vr,FracSoiPAsIce_vr,FracAirFilledSoilPore_vr
    implicit none
    integer, intent(in) :: ncolumns ! number of local ATS surface columns
    integer :: ny
    character(len=*), parameter :: subname='InstallSharedLayer0State'
    call PrintInfo('beg '//subname)
    do ny=1,ncolumns
      AREA_3D(3,0,ny,1)=column_area(ny)
      VLWatMicP_vr(0,ny,1)=a_WAT0(ny)
      VLiceMicP_vr(0,ny,1)=a_ICE0(ny)
      VLWatMicP1_vr(0,ny,1)=a_WAT0(ny)
      VLiceMicP1_vr(0,ny,1)=a_ICE0(ny)
      TKS_vr(0,ny,1)=a_TSURF(ny)
      TKSoil1_vr(0,ny,1)=a_TSURF(ny)
      VHeatCapacity_vr(0,ny,1)=cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,0,ny,1)) &
        +cpw*a_WAT0(ny)+cpi*a_ICE0(ny)
      VHeatCapacity1_vr(0,ny,1)=VHeatCapacity_vr(0,ny,1)
      PSISoilMatricP_vr(0,ny,1)=min(0._r8,a_PRES0(ny)*1.e-6_r8)
      ! Native geometry is dry litter plus excess physical water/ice volume.
      ! Respect ATS's reported pond depth too (some inputs use a water-equivalent
      ! ice density); keep depth and volume consistent for thermal resistance.
      VGeomLayer_vr(0,ny,1)=VLitR_col(ny,1)+max(column_area(ny)*a_POND0(ny), &
        a_WAT0(ny)+a_ICE0(ny)-VWatLitRHoldCapcity_col(ny,1))
      DLYR_3D(3,0,ny,1)=VGeomLayer_vr(0,ny,1)/column_area(ny)
      VLMicP_vr(0,ny,1)=VWatLitRHoldCapcity_col(ny,1)
      DLYRR_COL(ny,1)=DLYR_3D(3,0,ny,1)
      if(VGeomLayer_vr(0,ny,1)>0._r8)then
        FracSoiPAsWat_vr(0,ny,1)=a_WAT0(ny)/VGeomLayer_vr(0,ny,1)
        FracSoiPAsIce_vr(0,ny,1)=a_ICE0(ny)/VGeomLayer_vr(0,ny,1)
        FracAirFilledSoilPore_vr(0,ny,1)=max(0._r8,1._r8 &
          -FracSoiPAsWat_vr(0,ny,1)-FracSoiPAsIce_vr(0,ny,1))
      else
        FracSoiPAsWat_vr(0,ny,1)=0._r8; FracSoiPAsIce_vr(0,ny,1)=0._r8
        FracAirFilledSoilPore_vr(0,ny,1)=1._r8
      endif
      ! Coverage/albedo preparation must also read today's ATS boundary,
      ! rather than the donor budget left by the previous EcoSIM advance.
      if(allocated(XVLMobileWaterLitR_col))then
      XVLMobileWaterLitR_col(ny,1)=max(0._r8,a_WAT0(ny)+a_ICE0(ny)-VWatLitRHoldCapcity_col(ny,1))
      call PartLitSoilFractionM(0,0,ny,1)
      endif
    enddo
    call PrintInfo('end '//subname)
  end subroutine InstallSharedLayer0State
!------------------------------------------------------------------------------------------
  subroutine ecosim_layer0_retention(properties, pressure, patm, depth, derivative, status) bind(C)
    ! Description:
    ! Invert the existing two-branch litter water-potential curve, with no state mutation.
    ! Depth is retained liquid-equivalent storage; mobile hydrostatic ponding is added by ATS.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(c_double), intent(in) :: properties(n_layer0_properties) ! native per-area material parameters
    real(c_double), value, intent(in) :: pressure,patm ! absolute surface and atmospheric pressure [Pa]
    real(c_double), intent(out) :: depth ! retained liquid-equivalent depth [m]
    real(c_double), intent(out) :: derivative ! d depth / d pressure [m Pa-1]
    integer(c_int), intent(out) :: status ! zero on success; one for inconsistent native parameters
    real(r8) :: psi,theta,dtheta,x,fcpsi
    character(len=*), parameter :: subname='ecosim_layer0_retention'
    call PrintInfo('beg '//subname)
    depth=0._r8; derivative=0._r8; status=0
    if(.not.all(ieee_is_finite(properties)).or.any(properties(1:3)<0._r8))then
      status=1
    else if(properties(1)>0._r8)then
      if(properties(3)<=0._r8.or.properties(5)<=0._r8.or.properties(7)<=0._r8.or. &
         properties(9)<=0._r8.or.properties(10)<=0._r8.or.properties(12)<=0._r8.or. &
         properties(13)>=properties(14).or.properties(14)>=0._r8)then
        status=1
      else
        psi=(pressure-patm)*1.e-6_r8
        if(psi>=properties(14))then
          depth=properties(1)
        else
          psi=max(psi,properties(13))
          fcpsi=-exp(properties(6))
          if(psi<fcpsi)then
            theta=exp(properties(4)-properties(5)*(log(-psi)-properties(6))/properties(7))
            dtheta=-theta*properties(5)/(properties(7)*psi)
          else
            x=max(0._r8,(log(-psi)-properties(11))/properties(12))
            theta=exp(properties(8)-properties(9)*x**(1._r8/properties(10)))
            dtheta=0._r8
            if(x>0._r8)dtheta=-theta*properties(9)/properties(10)* &
              x**(1._r8/properties(10)-1._r8)/(properties(12)*psi)
          endif
          depth=min(properties(1),theta*properties(3))
          if(depth<properties(1).and.(pressure-patm)*1.e-6_r8>properties(13)) &
            derivative=properties(3)*dtheta*1.e-6_r8
        endif
      endif
    endif
    call PrintInfo('end '//subname)
  end subroutine ecosim_layer0_retention
end module ATSLayer0Mod
