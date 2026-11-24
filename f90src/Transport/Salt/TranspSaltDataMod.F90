module TranspSaltDataMod
  ! Own temporary salt-transport state and flux arrays for a configured mesh.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use DebugToolMod, only: PrintInfo
  use GridConsts
  use SoilPropertyDataType
  use TracerIDMod
  use AqueChemDatatype
implicit none
  public
  CHARACTER(LEN=*),private,PARAMETER :: MOD_FILENAME=&
  __FILE__

  real(r8) :: VFLWX=0.5_r8
  real(r8),allocatable ::  AquaIonDifusivty2_vr(:,:,:,:)                      !aqueous diffusivity of ions [m2/hr]
  real(r8),allocatable ::  trcSalt_FloXSurRof_flxM_2DH(:,:,:,:,:)             !2D flow of salt through surface runoff at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_SnowDrift_flxM_2D(:,:,:,:,:)               !salt flux through snow drift at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_FloXSurRof_flxM(:,:,:)                     !salt through surface runoff at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_AquaADV_Snow2Litr_flxM(:,:,:)              !Salt flux from snow to litter [mol d-2]
  real(r8),allocatable ::  trcSalt_AquaADV_Snow2Soil_flxM(:,:,:)              !Salt flux from snow to soil [mol d-2]

  real(r8),allocatable ::  trcSalt_solml2_vr(:,:,:,:)                         !copy of salt tracers for transport [mol d-2]
  real(r8),allocatable ::  trcSalt_RGeoChem_flxM_vr(:,:,:,:)                  !salt reduction due to geochemistry [mol d-2]
  real(r8),allocatable ::  trcSalt_AquaAdv_flxM_snvr(:,:,:,:)                 !salt advection by water advection in snow at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_Transp2Micp_flxM_vr(:,:,:,:)               !total 3D micropore flux at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_Transp2Macp_flxM_vr(:,:,:,:)               !total 3D macropore flux at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_Aqua_flxM_snvr(:,:,:,:)                    !Aquatic salt flow into snow layer at iteration M [mol d-2]
  real(r8),allocatable ::  SoluteDifusivitytscaledM_vr(:,:,:)                 !time scaled solute diffusivity
  real(r8),allocatable ::  trcSalt_MicpTranspFlxM_3D(:,:,:,:,:)               !salt 3D mairopore flux at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_MacpTranspFlxM_3D(:,:,:,:,:)               !salt 3D macropore flux at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_soHml2_vr(:,:,:,:)                         !copy of salt in macropore [mol d-2]
  real(r8),allocatable ::  trcSalt_Mac2MicPore_flxM_vr(:,:,:,:)               !salt exchange from macropore and micropore at iteration M [mol d-2] 
  real(r8),allocatable ::  trcSalt_ml2_snvr(:,:,:,:)                          ! snowpack salt dissolved tracers [mol d-2]
  real(r8),allocatable ::  trcSalt_Precip2LitrM(:,:,:)                        !wet precipiation gas flux to litter at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_Precip2MicpM(:,:,:)                        !wet precipiation gas flux to topsoil at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_FloXSurRof_flxM_col(:,:,:)                 !salt flux through litter due to runoff [mol d-2]
  real(r8),allocatable ::  trcSalt_SnowDrift_flxM_col(:,:,:)                  !salt flux due to snow drifting at iteration M [mol d-2]
  real(r8),allocatable ::  trcSalt_Irrig_flxM_vr(:,:,:,:)                     !salt flux due to irrigation at iteration M [mol d-2]

!------------------------------------------------------------------------------------------

contains
  subroutine InitTranspSaltData
  ! Allocate zeroed transport workspaces after mesh and salt IDs are initialized.
  implicit none
  character(len=*), parameter :: subname='InitTranspSaltData'

  call PrintInfo('beg '//subname)

  ! Micropore state includes surface litter at L=0; reaction fluxes are soil-only.
  allocate(trcSalt_solml2_vr(idsalt_beg:idsaltb_end,0:JZ,JY,JX)); trcSalt_solml2_vr=0._r8
  allocate(trcSalt_RGeoChem_flxM_vr(idsalt_beg:idsaltb_end,JZ,JY,JX)); trcSalt_RGeoChem_flxM_vr=0._r8

  allocate(trcSalt_Irrig_flxM_vr(idsalt_beg:idsaltb_end,JZ,JY,JX));trcSalt_Irrig_flxM_vr=0._r8
  allocate(AquaIonDifusivty2_vr(idsalt_beg:idsalt_mend,JZ,JY,JX));   AquaIonDifusivty2_vr=0._r8
  allocate(trcSalt_FloXSurRof_flxM_col(idsalt_beg:idsalt_end,JY,JX));  trcSalt_FloXSurRof_flxM_col=0._r8
  allocate(trcSalt_SnowDrift_flxM_col(idsalt_beg:idsalt_end,JY,JX)); trcSalt_SnowDrift_flxM_col=0._r8
  allocate(trcSalt_AquaADV_Snow2Litr_flxM(idsalt_beg:idsalt_end,JY,JX));  trcSalt_AquaADV_Snow2Litr_flxM=0._r8
  allocate(trcSalt_AquaADV_Snow2Soil_flxM(idsalt_beg:idsaltb_end,JY,JX)); trcSalt_AquaADV_Snow2Soil_flxM=0._r8

  allocate(trcSalt_FloXSurRof_flxM_2DH(idsalt_beg:idsalt_end,2,2,JV,JH));   trcSalt_FloXSurRof_flxM_2DH=0._r8
  allocate(trcSalt_SnowDrift_flxM_2D(idsalt_beg:idsaltb_end,2,2,JV,JH));     trcSalt_SnowDrift_flxM_2D=0._r8
  allocate(trcSalt_FloXSurRof_flxM(idsalt_beg:idsalt_end,JY,JX));     trcSalt_FloXSurRof_flxM=0._r8
  allocate(trcSalt_AquaAdv_flxM_snvr(idsalt_beg:idsalt_end,JS,JY,JX));   trcSalt_AquaAdv_flxM_snvr=0._r8
  allocate(trcSalt_Aqua_flxM_snvr(idsalt_beg:idsalt_end,JS,JY,JX)); trcSalt_Aqua_flxM_snvr=0._r8

  allocate(SoluteDifusivitytscaledM_vr(JZ,JY,JX));   SoluteDifusivitytscaledM_vr=0._r8
  allocate(trcSalt_MicpTranspFlxM_3D(idsalt_beg:idsaltb_end,3,0:JD,JV,JH));trcSalt_MicpTranspFlxM_3D=0._r8
  allocate(trcSalt_MacpTranspFlxM_3D(idsalt_beg:idsaltb_end,3,JD,JV,JH)); trcSalt_MacpTranspFlxM_3D=0._r8
  allocate(trcSalt_soHml2_vr(idsalt_beg:idsaltb_end,JZ,JY,JX));    trcSalt_soHml2_vr=0._r8
  allocate(trcSalt_Mac2MicPore_flxM_vr(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_Mac2MicPore_flxM_vr=0._r8

  allocate(trcSalt_ml2_snvr(idsalt_beg:idsalt_end,JS,JY,JX)); trcSalt_ml2_snvr=0._r8
  allocate(trcSalt_Transp2Micp_flxM_vr(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_Transp2Micp_flxM_vr=0._r8
  allocate(trcSalt_Transp2Macp_flxM_vr(idsalt_beg:idsaltb_end,JZ,JY,JX));   trcSalt_Transp2Macp_flxM_vr=0._r8
  allocate(trcSalt_Precip2LitrM(idsalt_beg:idsalt_end,JY,JX));      trcSalt_Precip2LitrM=0._r8
  allocate(trcSalt_Precip2MicpM(idsalt_beg:idsaltb_end,JY,JX));     trcSalt_Precip2MicpM=0._r8
  call PrintInfo('end '//subname)
  end subroutine InitTranspSaltData

!------------------------------------------------------------------------------------------
  subroutine DestructTranspSaltData
  ! Release every workspace owned by InitTranspSaltData, permitting a new mesh.
  use abortutils, only : destroy
  implicit none
  character(len=*), parameter :: subname='DestructTranspSaltData'

  call PrintInfo('beg '//subname)

  call destroy(AquaIonDifusivty2_vr)
  call destroy(trcSalt_AquaAdv_flxM_snvr)
  call destroy(trcSalt_AquaADV_Snow2Litr_flxM)
  call destroy(trcSalt_AquaADV_Snow2Soil_flxM)

  call destroy(trcSalt_ml2_snvr)
  call destroy(trcSalt_MacpTranspFlxM_3D)
  call destroy(trcSalt_soHml2_vr)
  call destroy(trcSalt_Aqua_flxM_snvr)
  call destroy(trcSalt_MicpTranspFlxM_3D)
  call destroy(trcSalt_solml2_vr)
  call destroy(trcSalt_RGeoChem_flxM_vr)
  call destroy(trcSalt_SnowDrift_flxM_2D)
  call destroy(SoluteDifusivitytscaledM_vr)
  call destroy(trcSalt_FloXSurRof_flxM)
  call destroy(trcSalt_FloXSurRof_flxM_2DH)
  call destroy(trcSalt_Mac2MicPore_flxM_vr)
  call destroy(trcSalt_Transp2Micp_flxM_vr)
  call destroy(trcSalt_Transp2Macp_flxM_vr)
  call destroy(trcSalt_Precip2LitrM)
  call destroy(trcSalt_Precip2MicpM)

  call destroy(trcSalt_FloXSurRof_flxM_col)
  call destroy(trcSalt_SnowDrift_flxM_col)
  call destroy(trcSalt_Irrig_flxM_vr)
  call PrintInfo('end '//subname)
  end subroutine DestructTranspSaltData

end module TranspSaltDataMod
