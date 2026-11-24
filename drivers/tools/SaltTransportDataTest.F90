program SaltTransportDataTest
  ! Exercise workspace bounds, deterministic initialization and complete cleanup.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils, only: iulog
  use DebugToolMod, only: PrintInfo
  use EcoSIMCtrlMod, only: lverb
  use GridConsts, only: JX,JY,JZ,JD,JV,JH,JS
  use TracerIDMod, only: InitTracerIDs,CleanUpTracerIDs,idsalt_beg,idsalt_end,idsaltb_end
  use TranspSaltDataMod
  implicit none
  integer :: cycle

  lverb=.false.
  call InitTracerIDs(.true.)
  do cycle=1,2
    JX=cycle; JY=3-cycle; JZ=cycle+2; JD=JZ; JV=JY+1; JH=JX+1
    call InitTranspSaltData()
    call check(allocated(AquaIonDifusivty2_vr),'allocated AquaIonDifusivty2_vr')
    call check(all(AquaIonDifusivty2_vr==0._r8),'zeroed AquaIonDifusivty2_vr')
    call check(allocated(trcSalt_FloXSurRof_flxM_2DH),'allocated trcSalt_FloXSurRof_flxM_2DH')
    call check(all(trcSalt_FloXSurRof_flxM_2DH==0._r8),'zeroed trcSalt_FloXSurRof_flxM_2DH')
    call check(allocated(trcSalt_SnowDrift_flxM_2D),'allocated trcSalt_SnowDrift_flxM_2D')
    call check(all(trcSalt_SnowDrift_flxM_2D==0._r8),'zeroed trcSalt_SnowDrift_flxM_2D')
    call check(allocated(trcSalt_FloXSurRof_flxM),'allocated trcSalt_FloXSurRof_flxM')
    call check(all(trcSalt_FloXSurRof_flxM==0._r8),'zeroed trcSalt_FloXSurRof_flxM')
    call check(allocated(trcSalt_AquaADV_Snow2Litr_flxM),'allocated trcSalt_AquaADV_Snow2Litr_flxM')
    call check(all(trcSalt_AquaADV_Snow2Litr_flxM==0._r8),'zeroed trcSalt_AquaADV_Snow2Litr_flxM')
    call check(allocated(trcSalt_AquaADV_Snow2Soil_flxM),'allocated trcSalt_AquaADV_Snow2Soil_flxM')
    call check(all(trcSalt_AquaADV_Snow2Soil_flxM==0._r8),'zeroed trcSalt_AquaADV_Snow2Soil_flxM')
    call check(allocated(trcSalt_solml2_vr),'allocated trcSalt_solml2_vr')
    call check(all(trcSalt_solml2_vr==0._r8),'zeroed trcSalt_solml2_vr')
    call check(allocated(trcSalt_RGeoChem_flxM_vr),'allocated trcSalt_RGeoChem_flxM_vr')
    call check(all(trcSalt_RGeoChem_flxM_vr==0._r8),'zeroed trcSalt_RGeoChem_flxM_vr')
    call check(allocated(trcSalt_AquaAdv_flxM_snvr),'allocated trcSalt_AquaAdv_flxM_snvr')
    call check(all(trcSalt_AquaAdv_flxM_snvr==0._r8),'zeroed trcSalt_AquaAdv_flxM_snvr')
    call check(allocated(trcSalt_Transp2Micp_flxM_vr),'allocated trcSalt_Transp2Micp_flxM_vr')
    call check(all(trcSalt_Transp2Micp_flxM_vr==0._r8),'zeroed trcSalt_Transp2Micp_flxM_vr')
    call check(allocated(trcSalt_Transp2Macp_flxM_vr),'allocated trcSalt_Transp2Macp_flxM_vr')
    call check(all(trcSalt_Transp2Macp_flxM_vr==0._r8),'zeroed trcSalt_Transp2Macp_flxM_vr')
    call check(allocated(trcSalt_Aqua_flxM_snvr),'allocated trcSalt_Aqua_flxM_snvr')
    call check(all(trcSalt_Aqua_flxM_snvr==0._r8),'zeroed trcSalt_Aqua_flxM_snvr')
    call check(allocated(SoluteDifusivitytscaledM_vr),'allocated SoluteDifusivitytscaledM_vr')
    call check(all(SoluteDifusivitytscaledM_vr==0._r8),'zeroed SoluteDifusivitytscaledM_vr')
    call check(allocated(trcSalt_MicpTranspFlxM_3D),'allocated trcSalt_MicpTranspFlxM_3D')
    call check(all(trcSalt_MicpTranspFlxM_3D==0._r8),'zeroed trcSalt_MicpTranspFlxM_3D')
    call check(allocated(trcSalt_MacpTranspFlxM_3D),'allocated trcSalt_MacpTranspFlxM_3D')
    call check(all(trcSalt_MacpTranspFlxM_3D==0._r8),'zeroed trcSalt_MacpTranspFlxM_3D')
    call check(allocated(trcSalt_soHml2_vr),'allocated trcSalt_soHml2_vr')
    call check(all(trcSalt_soHml2_vr==0._r8),'zeroed trcSalt_soHml2_vr')
    call check(allocated(trcSalt_Mac2MicPore_flxM_vr),'allocated trcSalt_Mac2MicPore_flxM_vr')
    call check(all(trcSalt_Mac2MicPore_flxM_vr==0._r8),'zeroed trcSalt_Mac2MicPore_flxM_vr')
    call check(allocated(trcSalt_ml2_snvr),'allocated trcSalt_ml2_snvr')
    call check(all(trcSalt_ml2_snvr==0._r8),'zeroed trcSalt_ml2_snvr')
    call check(allocated(trcSalt_Precip2LitrM),'allocated trcSalt_Precip2LitrM')
    call check(all(trcSalt_Precip2LitrM==0._r8),'zeroed trcSalt_Precip2LitrM')
    call check(allocated(trcSalt_Precip2MicpM),'allocated trcSalt_Precip2MicpM')
    call check(all(trcSalt_Precip2MicpM==0._r8),'zeroed trcSalt_Precip2MicpM')
    call check(allocated(trcSalt_FloXSurRof_flxM_col),'allocated trcSalt_FloXSurRof_flxM_col')
    call check(all(trcSalt_FloXSurRof_flxM_col==0._r8),'zeroed trcSalt_FloXSurRof_flxM_col')
    call check(allocated(trcSalt_SnowDrift_flxM_col),'allocated trcSalt_SnowDrift_flxM_col')
    call check(all(trcSalt_SnowDrift_flxM_col==0._r8),'zeroed trcSalt_SnowDrift_flxM_col')
    call check(allocated(trcSalt_Irrig_flxM_vr),'allocated trcSalt_Irrig_flxM_vr')
    call check(all(trcSalt_Irrig_flxM_vr==0._r8),'zeroed trcSalt_Irrig_flxM_vr')
    call check(lbound(trcSalt_solml2_vr,2)==0,'micropore workspace includes litter')
    call check(ubound(trcSalt_solml2_vr,2)==JZ,'micropore workspace includes deepest soil')
    call check(lbound(trcSalt_solml2_vr,1)==idsalt_beg,'first salt tracer')
    call check(ubound(trcSalt_solml2_vr,1)==idsaltb_end,'banded salt tracers')
    call check(all(shape(trcSalt_solml2_vr)==[idsaltb_end-idsalt_beg+1,JZ+1,JY,JX]),'mesh extents')
    call check(lbound(trcSalt_RGeoChem_flxM_vr,2)==1,'reaction flux begins in soil')
    call check(ubound(trcSalt_RGeoChem_flxM_vr,2)==JZ,'reaction flux reaches deepest soil')
    call check(size(trcSalt_AquaAdv_flxM_snvr,2)==JS,'snow flux layer extent')
    trcSalt_solml2_vr(:,0,:,:)=1._r8
    trcSalt_solml2_vr(idsaltb_end,JZ,JY,JX)=2._r8
    trcSalt_RGeoChem_flxM_vr(idsaltb_end,JZ,JY,JX)=3._r8
    trcSalt_AquaAdv_flxM_snvr(idsalt_end,JS,JY,JX)=4._r8
    call DestructTranspSaltData()
    call check(.not.allocated(AquaIonDifusivty2_vr),'released AquaIonDifusivty2_vr')
    call check(.not.allocated(trcSalt_FloXSurRof_flxM_2DH),'released trcSalt_FloXSurRof_flxM_2DH')
    call check(.not.allocated(trcSalt_SnowDrift_flxM_2D),'released trcSalt_SnowDrift_flxM_2D')
    call check(.not.allocated(trcSalt_FloXSurRof_flxM),'released trcSalt_FloXSurRof_flxM')
    call check(.not.allocated(trcSalt_AquaADV_Snow2Litr_flxM),'released trcSalt_AquaADV_Snow2Litr_flxM')
    call check(.not.allocated(trcSalt_AquaADV_Snow2Soil_flxM),'released trcSalt_AquaADV_Snow2Soil_flxM')
    call check(.not.allocated(trcSalt_solml2_vr),'released trcSalt_solml2_vr')
    call check(.not.allocated(trcSalt_RGeoChem_flxM_vr),'released trcSalt_RGeoChem_flxM_vr')
    call check(.not.allocated(trcSalt_AquaAdv_flxM_snvr),'released trcSalt_AquaAdv_flxM_snvr')
    call check(.not.allocated(trcSalt_Transp2Micp_flxM_vr),'released trcSalt_Transp2Micp_flxM_vr')
    call check(.not.allocated(trcSalt_Transp2Macp_flxM_vr),'released trcSalt_Transp2Macp_flxM_vr')
    call check(.not.allocated(trcSalt_Aqua_flxM_snvr),'released trcSalt_Aqua_flxM_snvr')
    call check(.not.allocated(SoluteDifusivitytscaledM_vr),'released SoluteDifusivitytscaledM_vr')
    call check(.not.allocated(trcSalt_MicpTranspFlxM_3D),'released trcSalt_MicpTranspFlxM_3D')
    call check(.not.allocated(trcSalt_MacpTranspFlxM_3D),'released trcSalt_MacpTranspFlxM_3D')
    call check(.not.allocated(trcSalt_soHml2_vr),'released trcSalt_soHml2_vr')
    call check(.not.allocated(trcSalt_Mac2MicPore_flxM_vr),'released trcSalt_Mac2MicPore_flxM_vr')
    call check(.not.allocated(trcSalt_ml2_snvr),'released trcSalt_ml2_snvr')
    call check(.not.allocated(trcSalt_Precip2LitrM),'released trcSalt_Precip2LitrM')
    call check(.not.allocated(trcSalt_Precip2MicpM),'released trcSalt_Precip2MicpM')
    call check(.not.allocated(trcSalt_FloXSurRof_flxM_col),'released trcSalt_FloXSurRof_flxM_col')
    call check(.not.allocated(trcSalt_SnowDrift_flxM_col),'released trcSalt_SnowDrift_flxM_col')
    call check(.not.allocated(trcSalt_Irrig_flxM_vr),'released trcSalt_Irrig_flxM_vr')
    ! Cleanup is safe when called again before reinitializing a different mesh.
    call DestructTranspSaltData()
  enddo
  call CleanUpTracerIDs()
  write(*,*)'PASS: salt workspace allocation, litter/band bounds, zeroing and cleanup on two meshes'
contains
!------------------------------------------------------------------------------------------
  subroutine check(Passed,Label)
  ! Fail if the salt workspace lifecycle contract is violated.
  logical, intent(in) :: Passed        !whether the lifecycle condition holds
  character(len=*), intent(in) :: Label !condition being checked
  character(len=*), parameter :: subname='check'
  call PrintInfo('beg '//subname)
  if(.not.Passed)then
    write(iulog,*)'FAIL: ',Label
    error stop 1
  endif
  call PrintInfo('end '//subname)
  end subroutine check
end program SaltTransportDataTest
