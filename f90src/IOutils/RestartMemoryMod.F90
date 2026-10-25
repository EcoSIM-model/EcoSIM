module RestartMemoryMod
  ! Retain lagged process inputs and constitutive caches, in addition to pools.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL, ispval => DAT_CONST_ISPVAL
  use DebugToolMod, only : PrintInfo
  use abortutils, only : endrun, iulog
  use GridConsts, only : bounds
  use PlantMgmtDataType, only : NP_col
  use restUtilMod, only : restartvar, cpcol, cppft
  use ncdio_pio
  use SOMDataType, only : EPOC_vr
  use SOMDataType, only : EHUM_vr
  use FlagDataType, only : iResetSoilProf_col
  use SoilPhysDataType, only : LOGPOROS_vr
  use SoilPhysDataType, only : LOGFldCapacity_vr
  use SoilPhysDataType, only : PSD_vr
  use SoilPhysDataType, only : FCD_vr
  use SoilPhysDataType, only : SRP_vr
  use SoilHeatDataType, only : NumerSolidThermCond_vr
  use SoilHeatDataType, only : DenomSolidThermCond_vr
  use SoilWaterDataType, only : HydroCond_3D
  use SoilWaterDataType, only : HydroCondMacP_vr
  use SoilWaterDataType, only : ThetaSat_vr
  use PlantTraitDataType, only : BulkFactor4Snow_pft
  use PlantDataRateType, only : THeatLossRoot2Soil_vr
  use PlantDataRateType, only : RCH4EcoDmnd_vr
  use CanopyDataType, only : fSnowCanopy_col
  use CanopyDataType, only : fSnowCanopy_pft
  use CanopyDataType, only : TKCanopy_pft
  use SoilBGCDataType, only : trcg_TotalMass_col
  use RootDataType, only : RootMycoMassElm_vr
  use SoilPropertyDataType, only : PathLenMacPore_vr
  use SoilPropertyDataType, only : MacPoreRadius_vr
  use PlantTraitDataType, only : CanopyHeightDead_pft
  use PlantTraitDataType, only : StalkAveRadius_pft
  use PlantTraitDataType, only : BulkFactor4Snow_col
  use SoilPhysDataType, only : LOGWiltPoint_vr
  use PlantDataRateType, only : GrossRespC_CumYr_pft
  use SoilWaterDataType, only : WatMass_col
  use SoilHeatDataType, only : HeatStore_col
  use CanopyDataType, only : CanopyWaterMassEnd_col
  use SnowDataType, only : SnowMassEnd_col
  use SnowDataType, only : SnowEngyEnd_col
  use SurfLitterDataType, only : LitWatMassEnd_col
  use SoilWaterDataType, only : SoilWatMassEnd_col
  use SoilBGCDataType, only : trcg_soilMass_col
  use SoilBGCDataType, only : trcg_rootMass_col
  use SnowDataType, only : trcg_snowMass_col
  use ClimForcDataType, only : NWetDep_col
  use PlantDataRateType, only : RootUptk_Nmin_cumYr_pft
  use PlantDataRateType, only : RootUptk_Pmin_cumYr_pft
  use AqueChemDatatype, only : trcg_mass_cumerr_col
  use ClimForcDataType, only : TRAD_col
  use ClimForcDataType, only : HUDX_col
  use ClimForcDataType, only : HUDN_col
  use ClimForcDataType, only : TWIND_col
  use CanopyDataType, only : RubiscoActivity_brch
  implicit none
  private
  public :: restart_process_memory, has_restart_memory
  logical, save :: has_restart_memory=.false.
contains
!------------------------------------------------------------------------------------------
  subroutine restart_process_memory(ncid,flag)
  type(file_desc_t), intent(inout) :: ncid
  character(len=*), intent(in) :: flag
  character(len=*), parameter :: subname='restart_process_memory'
  integer :: varid,dimid
  logical :: found
  type(var_desc_t) :: vardesc
  real(r8), allocatable, target :: c1(:),c2(:,:),c3(:,:,:),c4(:,:,:,:)
  real(r8), pointer :: p1(:),p2(:,:),p3(:,:,:),p4(:,:,:,:)
  integer, allocatable, target :: i1(:)
  integer, pointer :: ip1(:)
  call PrintInfo("beg "//subname)
  associate(NHW=>bounds%NHW,NHE=>bounds%NHE,NVN=>bounds%NVN,NVS=>bounds%NVS, &
    begc=>bounds%begc,endc=>bounds%endc,begp=>bounds%begp,endp=>bounds%endp)
  if(flag=='read')then
    call ncd_inqvid(ncid,'continuation_memory_version',varid,vardesc,readvar=found)
    has_restart_memory=found
    if(.not.found)then
      write(iulog,*) 'Legacy checkpoint lacks process memory; exact continuation is unavailable.'
      call PrintInfo("end "//subname)
      return
    endif
  endif
  if(flag=='define')then
    call ncd_defdim(ncid,'memory_water_class',SIZE(HydroCond_3D,2),dimid)
    call ncd_defdim(ncid,'memory_gas_species',SIZE(trcg_TotalMass_col,1),dimid)
    call ncd_defdim(ncid,'memory_gas_unbanded',SIZE(trcg_rootMass_col,1),dimid)
  endif
  allocate(i1(begc:endc)); i1=1; ip1=>i1
  call restartvar(ncid,flag,varname='continuation_memory_version',dim1name='column', &
    long_name='Version of restart process-memory extension',units='1',interpinic_flag='skip', &
    data=ip1,missing_value=ispval,fill_value=ispval)
  if(flag=='read')then
    if(any(i1/=1))then
      write(iulog,*) 'Unsupported continuation-memory version: ',i1
      call endrun('Unsupported restart process-memory schema',__LINE__)
    endif
  endif
  deallocate(i1)

  ! EPOC_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(EPOC_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,EPOC_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_EPOC_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: EPOC_vr',units='-',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,EPOC_vr,c2)
  deallocate(c2)

  ! EHUM_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(EHUM_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,EHUM_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_EHUM_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: EHUM_vr',units='-',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,EHUM_vr,c2)
  deallocate(c2)

  ! iResetSoilProf_col: retain the value used by the next process call.
  allocate(i1(begc:endc)); i1=ispval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,iResetSoilProf_col,i1)
  ip1=>i1
  call restartvar(ncid,flag,varname='memory_iResetSoilProf_col', &
    dim1name='column', &
    long_name='Restart process memory: iResetSoilProf_col',units='-',interpinic_flag='skip', &
    data=ip1,missing_value=ispval,fill_value=ispval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,iResetSoilProf_col,i1)
  deallocate(i1)

  ! LOGPOROS_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(LOGPOROS_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LOGPOROS_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_LOGPOROS_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: LOGPOROS_vr',units='-',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,LOGPOROS_vr,c2)
  deallocate(c2)

  ! LOGFldCapacity_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(LOGFldCapacity_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LOGFldCapacity_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_LOGFldCapacity_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: LOGFldCapacity_vr',units='-',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,LOGFldCapacity_vr,c2)
  deallocate(c2)

  ! PSD_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(PSD_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PSD_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_PSD_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: PSD_vr',units='-',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,PSD_vr,c2)
  deallocate(c2)

  ! FCD_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(FCD_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,FCD_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_FCD_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: FCD_vr',units='-',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,FCD_vr,c2)
  deallocate(c2)

  ! SRP_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(SRP_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SRP_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_SRP_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: SRP_vr',units='-',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,SRP_vr,c2)
  deallocate(c2)

  ! NumerSolidThermCond_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(NumerSolidThermCond_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NumerSolidThermCond_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_NumerSolidThermCond_vr', &
    dim1name='column',dim2name='levsoi', &
    long_name='Restart process memory: NumerSolidThermCond_vr',units='MJ m h-1 K-1',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,NumerSolidThermCond_vr,c2)
  deallocate(c2)

  ! DenomSolidThermCond_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(DenomSolidThermCond_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,DenomSolidThermCond_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_DenomSolidThermCond_vr', &
    dim1name='column',dim2name='levsoi', &
    long_name='Restart process memory: DenomSolidThermCond_vr',units='MJ K-1',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,DenomSolidThermCond_vr,c2)
  deallocate(c2)

  ! HydroCond_3D: retain the value used by the next process call.
  allocate(c4(begc:endc,SIZE(HydroCond_3D,1),SIZE(HydroCond_3D,2),SIZE(HydroCond_3D,3))); c4=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroCond_3D,c4)
  p4=>c4
  call restartvar(ncid,flag,varname='memory_HydroCond_3D', &
    dim1name='column',dim2name='sdim',dim3name='memory_water_class',dim4name='levsoi1', &
    long_name='Restart process memory: HydroCond_3D',units='m MPa-1 h-1',interpinic_flag='skip', &
    data=p4,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,HydroCond_3D,c4)
  deallocate(c4)

  ! HydroCondMacP_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(HydroCondMacP_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HydroCondMacP_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_HydroCondMacP_vr', &
    dim1name='column',dim2name='levsoi', &
    long_name='Restart process memory: HydroCondMacP_vr',units='m MPa-1 h-1',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,HydroCondMacP_vr,c2)
  deallocate(c2)

  ! ThetaSat_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(ThetaSat_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,ThetaSat_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_ThetaSat_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: ThetaSat_vr',units='1',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,ThetaSat_vr,c2)
  deallocate(c2)

  ! BulkFactor4Snow_pft: retain the value used by the next process call.
  allocate(c1(begp:endp)); c1=spval
  if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,BulkFactor4Snow_pft,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_BulkFactor4Snow_pft', &
    dim1name='pft', &
    long_name='Restart process memory: BulkFactor4Snow_pft',units='m2 (kg SWE)-1',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,BulkFactor4Snow_pft,c1)
  deallocate(c1)

  ! THeatLossRoot2Soil_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(THeatLossRoot2Soil_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,THeatLossRoot2Soil_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_THeatLossRoot2Soil_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: THeatLossRoot2Soil_vr',units='MJ d-2',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,THeatLossRoot2Soil_vr,c2)
  deallocate(c2)

  ! RCH4EcoDmnd_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(RCH4EcoDmnd_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RCH4EcoDmnd_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_RCH4EcoDmnd_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: RCH4EcoDmnd_vr',units='gCd-2 h-1',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,RCH4EcoDmnd_vr,c2)
  deallocate(c2)

  ! fSnowCanopy_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,fSnowCanopy_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_fSnowCanopy_col', &
    dim1name='column', &
    long_name='Restart process memory: fSnowCanopy_col',units='-',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,fSnowCanopy_col,c1)
  deallocate(c1)

  ! fSnowCanopy_pft: retain the value used by the next process call.
  allocate(c1(begp:endp)); c1=spval
  if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,fSnowCanopy_pft,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_fSnowCanopy_pft', &
    dim1name='pft', &
    long_name='Restart process memory: fSnowCanopy_pft',units='-',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,fSnowCanopy_pft,c1)
  deallocate(c1)

  ! TKCanopy_pft: retain the value used by the next process call.
  allocate(c1(begp:endp)); c1=spval
  if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,TKCanopy_pft,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_TKCanopy_pft', &
    dim1name='pft', &
    long_name='Restart process memory: TKCanopy_pft',units='K',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,TKCanopy_pft,c1)
  deallocate(c1)

  ! trcg_TotalMass_col: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(trcg_TotalMass_col,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_TotalMass_col,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_trcg_TotalMass_col', &
    dim1name='column',dim2name='memory_gas_species', &
    long_name='Restart process memory: trcg_TotalMass_col',units='g d-2',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_TotalMass_col,c2)
  deallocate(c2)

  ! RootMycoMassElm_vr: retain the value used by the next process call.
  allocate(c4(begc:endc,SIZE(RootMycoMassElm_vr,1),SIZE(RootMycoMassElm_vr,2),SIZE(RootMycoMassElm_vr,3))); c4=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,RootMycoMassElm_vr,c4)
  p4=>c4
  call restartvar(ncid,flag,varname='memory_RootMycoMassElm_vr', &
    dim1name='column',dim2name='elmnts',dim3name='rootyps',dim4name='levsoi', &
    long_name='Restart process memory: RootMycoMassElm_vr',units='g d-2',interpinic_flag='skip', &
    data=p4,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,RootMycoMassElm_vr,c4)
  deallocate(c4)

  ! PathLenMacPore_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(PathLenMacPore_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,PathLenMacPore_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_PathLenMacPore_vr', &
    dim1name='column',dim2name='levsoi', &
    long_name='Restart process memory: PathLenMacPore_vr',units='m',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,PathLenMacPore_vr,c2)
  deallocate(c2)

  ! MacPoreRadius_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(MacPoreRadius_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,MacPoreRadius_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_MacPoreRadius_vr', &
    dim1name='column',dim2name='levsoi', &
    long_name='Restart process memory: MacPoreRadius_vr',units='m',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,MacPoreRadius_vr,c2)
  deallocate(c2)

  ! CanopyHeightDead_pft: retain the value used by the next process call.
  allocate(c1(begp:endp)); c1=spval
  if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,CanopyHeightDead_pft,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_CanopyHeightDead_pft', &
    dim1name='pft', &
    long_name='Restart process memory: CanopyHeightDead_pft',units='m',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,CanopyHeightDead_pft,c1)
  deallocate(c1)

  ! StalkAveRadius_pft: retain the value used by the next process call.
  allocate(c1(begp:endp)); c1=spval
  if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,StalkAveRadius_pft,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_StalkAveRadius_pft', &
    dim1name='pft', &
    long_name='Restart process memory: StalkAveRadius_pft',units='m',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,StalkAveRadius_pft,c1)
  deallocate(c1)

  ! BulkFactor4Snow_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,BulkFactor4Snow_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_BulkFactor4Snow_col', &
    dim1name='column', &
    long_name='Restart process memory: BulkFactor4Snow_col',units='m2 (kg SWE)-1',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,BulkFactor4Snow_col,c1)
  deallocate(c1)

  ! LOGWiltPoint_vr: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(LOGWiltPoint_vr,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LOGWiltPoint_vr,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_LOGWiltPoint_vr', &
    dim1name='column',dim2name='levsoi1', &
    long_name='Restart process memory: LOGWiltPoint_vr',units='-',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,LOGWiltPoint_vr,c2)
  deallocate(c2)

  ! GrossRespC_CumYr_pft: retain the value used by the next process call.
  allocate(c1(begp:endp)); c1=spval
  if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,GrossRespC_CumYr_pft,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_GrossRespC_CumYr_pft', &
    dim1name='pft', &
    long_name='Restart process memory: GrossRespC_CumYr_pft',units='g d-2 ',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,GrossRespC_CumYr_pft,c1)
  deallocate(c1)

  ! WatMass_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,WatMass_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_WatMass_col', &
    dim1name='column', &
    long_name='Restart process memory: WatMass_col',units='m3 d-2',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,WatMass_col,c1)
  deallocate(c1)

  ! HeatStore_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HeatStore_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_HeatStore_col', &
    dim1name='column', &
    long_name='Restart process memory: HeatStore_col',units='MJ d-2',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,HeatStore_col,c1)
  deallocate(c1)

  ! CanopyWaterMassEnd_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyWaterMassEnd_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_CanopyWaterMassEnd_col', &
    dim1name='column', &
    long_name='Restart process memory: CanopyWaterMassEnd_col',units='m3 d-2',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,CanopyWaterMassEnd_col,c1)
  deallocate(c1)

  ! SnowMassEnd_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnowMassEnd_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_SnowMassEnd_col', &
    dim1name='column', &
    long_name='Restart process memory: SnowMassEnd_col',units='m3 H2O d-2',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,SnowMassEnd_col,c1)
  deallocate(c1)

  ! SnowEngyEnd_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SnowEngyEnd_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_SnowEngyEnd_col', &
    dim1name='column', &
    long_name='Restart process memory: SnowEngyEnd_col',units='MJ d-2',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,SnowEngyEnd_col,c1)
  deallocate(c1)

  ! LitWatMassEnd_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,LitWatMassEnd_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_LitWatMassEnd_col', &
    dim1name='column', &
    long_name='Restart process memory: LitWatMassEnd_col',units='m3 H2O d-2',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,LitWatMassEnd_col,c1)
  deallocate(c1)

  ! SoilWatMassEnd_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,SoilWatMassEnd_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_SoilWatMassEnd_col', &
    dim1name='column', &
    long_name='Restart process memory: SoilWatMassEnd_col',units='1',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,SoilWatMassEnd_col,c1)
  deallocate(c1)

  ! trcg_soilMass_col: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(trcg_soilMass_col,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_soilMass_col,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_trcg_soilMass_col', &
    dim1name='column',dim2name='memory_gas_species', &
    long_name='Restart process memory: trcg_soilMass_col',units='g d-2',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_soilMass_col,c2)
  deallocate(c2)

  ! trcg_rootMass_col: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(trcg_rootMass_col,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_rootMass_col,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_trcg_rootMass_col', &
    dim1name='column',dim2name='memory_gas_unbanded', &
    long_name='Restart process memory: trcg_rootMass_col',units='g d-2',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_rootMass_col,c2)
  deallocate(c2)

  ! trcg_snowMass_col: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(trcg_snowMass_col,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_snowMass_col,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_trcg_snowMass_col', &
    dim1name='column',dim2name='memory_gas_unbanded', &
    long_name='Restart process memory: trcg_snowMass_col',units='g d-2',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_snowMass_col,c2)
  deallocate(c2)

  ! NWetDep_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,NWetDep_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_NWetDep_col', &
    dim1name='column', &
    long_name='Restart process memory: NWetDep_col',units='gN m-2 h-1',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,NWetDep_col,c1)
  deallocate(c1)

  ! RootUptk_Nmin_cumYr_pft: retain the value used by the next process call.
  allocate(c1(begp:endp)); c1=spval
  if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,RootUptk_Nmin_cumYr_pft,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_RootUptk_Nmin_cumYr_pft', &
    dim1name='pft', &
    long_name='Restart process memory: RootUptk_Nmin_cumYr_pft',units='gN d-2',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,RootUptk_Nmin_cumYr_pft,c1)
  deallocate(c1)

  ! RootUptk_Pmin_cumYr_pft: retain the value used by the next process call.
  allocate(c1(begp:endp)); c1=spval
  if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,RootUptk_Pmin_cumYr_pft,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_RootUptk_Pmin_cumYr_pft', &
    dim1name='pft', &
    long_name='Restart process memory: RootUptk_Pmin_cumYr_pft',units='gP d-2',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,RootUptk_Pmin_cumYr_pft,c1)
  deallocate(c1)

  ! trcg_mass_cumerr_col: retain the value used by the next process call.
  allocate(c2(begc:endc,SIZE(trcg_mass_cumerr_col,1))); c2=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_mass_cumerr_col,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_trcg_mass_cumerr_col', &
    dim1name='column',dim2name='memory_gas_unbanded', &
    long_name='Restart process memory: trcg_mass_cumerr_col',units='g/d2',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,trcg_mass_cumerr_col,c2)
  deallocate(c2)

  ! TRAD_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TRAD_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_TRAD_col', &
    dim1name='column', &
    long_name='Restart process memory: TRAD_col',units='MJ d-1',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,TRAD_col,c1)
  deallocate(c1)

  ! HUDX_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HUDX_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_HUDX_col', &
    dim1name='column', &
    long_name='Restart process memory: HUDX_col',units='kPa',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,HUDX_col,c1)
  deallocate(c1)

  ! HUDN_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,HUDN_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_HUDN_col', &
    dim1name='column', &
    long_name='Restart process memory: HUDN_col',units='kPa',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,HUDN_col,c1)
  deallocate(c1)

  ! TWIND_col: retain the value used by the next process call.
  allocate(c1(begc:endc)); c1=spval
  if(flag=='write')call cpcol(flag,NHW,NHE,NVN,NVS,TWIND_col,c1)
  p1=>c1
  call restartvar(ncid,flag,varname='memory_TWIND_col', &
    dim1name='column', &
    long_name='Restart process memory: TWIND_col',units='m d-1',interpinic_flag='skip', &
    data=p1,missing_value=spval,fill_value=spval)
  if(flag=='read')call cpcol(flag,NHW,NHE,NVN,NVS,TWIND_col,c1)
  deallocate(c1)

  ! RubiscoActivity_brch: retain the value used by the next process call.
  allocate(c2(begp:endp,SIZE(RubiscoActivity_brch,1))); c2=spval
  if(flag=='write')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,RubiscoActivity_brch,c2)
  p2=>c2
  call restartvar(ncid,flag,varname='memory_RubiscoActivity_brch', &
    dim1name='pft',dim2name='nbranches', &
    long_name='Restart process memory: RubiscoActivity_brch',units='-',interpinic_flag='skip', &
    data=p2,missing_value=spval,fill_value=spval)
  if(flag=='read')call cppft(flag,NHW,NHE,NVN,NVS,NP_col,RubiscoActivity_brch,c2)
  deallocate(c2)
  end associate
  call PrintInfo("end "//subname)
  end subroutine restart_process_memory
end module RestartMemoryMod
