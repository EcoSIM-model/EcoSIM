module ForcWriterMod
  ! Description:
  ! Export initialization data for the batch soil biogeochemistry driver.
  use abortutils, only: endrun
  use SoilPhysDataType
  use SOMDataType
  USE MicrobialDataType
  use SoilBGCDataType
  use GridDataType
  use GridConsts
  use AqueChemDatatype
  use SoilPropertyDataType
  use ClimForcDataType
  use SoilWaterDataType
  use ncdio_pio
  use EcoSIMConfig, only : jcplx => jcplxc
  use data_const_mod, only : spval  => DAT_CONST_SPVAL
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  type, public :: bgc_forc_config_type
    logical :: laddband
    integer :: year
    integer :: doy
    integer :: layer
    character(len=64) :: bgc_fname
  end type bgc_forc_config_type
  type(bgc_forc_config_type),public :: bgc_forc_conf

  logical, public :: do_bgcforc_write
  public :: WriteBBGCForc
  contains

!------------------------------------------------------------------------------------------
  subroutine WriteBBGCFORC(doy,year)
  ! Description:
  ! Write initial condition for the batch soil bgc model.
  ! Because the single layer bgc does not have band information.
  ! The band nutrients is set to zero, or add to the bulk soil
  implicit none
  integer,intent(in) :: doy
  integer,intent(in) :: year
  integer :: NY,NX,L
  integer :: k,kl,nl
  type(file_desc_t) :: ncf
  integer :: recordDimID
  real(r8),pointer :: data1d(:)
  real(r8), allocatable :: dat1d(:)

  ! The batch driver reconstructs the EcoSIM curve from legacy scalar fields.
  ! Reject VG export rather than create forcing with a different retention law.
  if(allocated(VGThetaR_vr))call endrun('Batch-soil forcing export supports only the ecosim hydraulic scheme')
  if(year==bgc_forc_conf%year .and. doy==bgc_forc_conf%doy)then
    NY=1;NX=1;L=bgc_forc_conf%Layer
    write(*,*)'write bbgc forc on day ',doy,'year',year
    call ncd_pio_createfile(ncf, trim(bgc_forc_conf%bgc_fname))
    call ncd_defdim(ncf,'jcplx',jcplx,recordDimID)
    call ncd_defdim(ncf,'jsken',jsken,recordDimID)
    call ncd_defdim(ncf,'ndbiomcp',2,recordDimID)
    call ncd_defdim(ncf,'nlbiomcp',3,recordDimID)
    call ncd_defdim(ncf,'NumLiveAutoBioms',NumLiveAutoBioms,recordDimID)
    call ncd_defdim(ncf,'NumLiveHeterBioms',NumLiveHeterBioms,recordDimID)
    call ncd_defdim(ncf,'NumHetetr1MicCmplx',NumHetetr1MicCmplx,recordDimID)
    call ncd_defdim(ncf,'NumMicrobAutoTrophCmplx',NumMicrobAutoTrophCmplx,recordDimID)
    call ncd_defdim(ncf,'element',NumPlantChemElms,recordDimID)
    call ncd_defdim(ncf,'ndoms',trc_confs%NDOMS,recordDimID)
    call ncd_defvar(ncf, 'pH', ncd_float, long_name='soil pH',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'VLSoilPoreMicP_vr', ncd_float, long_name='volume of soil layer',  &
            units='m3 d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ORGC', ncd_float, long_name='total soil organic C',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ElmAllocmatMicrblitr2POM', ncd_float, dim1name='ndbiomcp', &
            long_name='allocation coefficient to humus fractions',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'VLSoilMicP', ncd_float, long_name='micropore volume',  &
            units='m3 d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'BKVL', ncd_float, long_name='mass of soil layer',  &
            units='Mg d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'POROS', ncd_float, long_name='soil porosity',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'FC', ncd_float, long_name='water contents at field capacity',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'SRP', ncd_float, long_name='shape parameter for water desorption',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'WP', ncd_float, long_name='water contents at wilting point',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'EHUM',ncd_float, &
            long_name='partitioning coefficient between humus and microbial residue',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'DLYR3', ncd_float, long_name='thickness of soil layer',  &
            units='m', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CEC', ncd_float, long_name='soil cation exchange capacity',  &
            units='cmol kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'XCEC', ncd_float, long_name='soil cation exchange capacity',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'AEC', ncd_float, long_name='soil anion exchange capacity',  &
            units='cmol kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'XAEC', ncd_float, long_name='soil anion exchange capacity',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CFE', ncd_float, long_name='soil Fe content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCA', ncd_float, long_name='soil Ca content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CMG', ncd_float, long_name='soil Mg content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CNA', ncd_float, long_name='soil Na content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CKA', ncd_float, long_name='soil K content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CSO4', ncd_float, long_name='soil SO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCL', ncd_float, long_name='soil Cl content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CAL', ncd_float, long_name='soil Al content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZMG', ncd_float, long_name='soil aqueous Mg content micropore',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZNA', ncd_float, long_name='soil aqueous Na content micropore',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZKA', ncd_float, long_name='soil aqueous K content micropore',  &
            units='mol d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CALPO', ncd_float, long_name='soil AlPO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CFEPO', ncd_float, long_name='soil FePO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCAPD', ncd_float, long_name='soil CaHPO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'EPOC', ncd_float, &
            long_name='partitioning coefficient between POC and litter',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCAPH', ncd_float, long_name='soil apatite content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CALOH', ncd_float, long_name='soil AlOH3 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CFEOH', ncd_float, long_name='soil FeOH3 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCACO', ncd_float, long_name='soil CaCO3 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CCASO', ncd_float, long_name='soil CaSO4 content',  &
            units='mg kg-1', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'FracBulkSOMC_vr', ncd_float, dim1name='jcplx',&
            long_name='fraction of total organic C in complex',  &
            units='none', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OQM', ncd_float, dim1name='ndom',dim2name='jcplx',&
            long_name='dissolved organic matter micropore',  &
            units='g d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OHM', ncd_float, dim1name='ndoms',dim2name='jcplx',&
            long_name='adsorbed soil organic matter',  &
            units='g d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ATCS',ncd_float,long_name='Mean annual temperature',  &
            units='oC', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'CNOSC', ncd_float, dim1name='jsken',&
            dim2name='jcplx',long_name='N:C ratios of SOM kinetic components in each complex',  &
            units='gN/gC', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'CPOSC', ncd_float, dim1name='jsken',&
            dim2name='jcplx',long_name='P:C ratios of SOM kinetic components in each complex',  &
            units='gN/gC', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OSM', ncd_float, dim1name='element',dim2name='jsken',&
            dim3name='jcplx',long_name='humus soil C in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'OSA', ncd_float, dim1name='jsken',&
            dim2name='jcplx',long_name='colonized soil C in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ORM', ncd_float, dim1name='element',dim2name='ndbiomcp',&
            dim3name='jcplx',long_name='microbial residue C in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'mBiomeHeter', ncd_float,dim1name='element',&
            dim2name='NumLiveHeterBioms',dim3name='jcplx',&
            long_name='heterotrophic microbial biomass element in each complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'mBiomeAutor', ncd_float, dim1name='element',dim2name='NumLiveAutoBioms',&
            long_name='microbial biomass element in autotrophic complex',  &
            units='gC d-2', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'BKDS', ncd_float,long_name='soil bulk density',&
            units='Mg m-3', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'THETY', ncd_float,long_name='air-dry water content',&
            units='m3 m-3', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSIMX', ncd_float,long_name='log water potential at field capacity',&
            units='log(MPa)', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSIMS', ncd_float,long_name='log water potential at saturation',&
            units='log(MPa)', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSIMD', ncd_float, &
            long_name='log water potential at wilting point - log water potential at field capacity',&
            units='log(MPa)', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSISD', ncd_float, &
            long_name='log water potential at saturation - log water potential at field capacity',&
            units='log(MPa)', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'PSISE', ncd_float, long_name='water potential at saturation',&
            units='MPa', missing_value=spval, fill_value=spval)

    call ncd_defvar(ncf, 'ZNH4S', ncd_float,long_name='NH4 in micropore',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZNO3S', ncd_float,long_name='NO3(-) in micropore',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'ZNO2S', ncd_float,long_name='NO2(-) in micropore',  &
            units='gN d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'H2PO4', ncd_float,long_name='H2PO4 in micropore',  &
            units='gP d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf, 'H1PO4', ncd_float,long_name='H1PO4 in micropore',  &
            units='gP d-2', missing_value=spval, fill_value=spval)
    call ncd_defvar(ncf,'ATKA', ncd_float, long_name='mean annual air temperature',&
            units='K', missing_value=spval, fill_value=spval)
    call ncd_enddef(ncf)

    call ncd_putvar(ncf,'pH',PH_vr(L,NY,NX))
    call ncd_putvar(ncf,'VLSoilPoreMicP_vr',VLSoilPoreMicP_vr(L,NY,NX))
    call ncd_putvar(ncf,'ORGC',SoilOrgM_vr(ielmc,L,NY,NX))
    call ncd_putvar(ncf,'ElmAllocmatMicrblitr2POM',ElmAllocmatMicrblitr2POM_vr(:,L,NY,NX))
    call ncd_putvar(ncf,'VLSoilMicP',VLSoilMicP_vr(L,NY,NX))
    call ncd_putvar(ncf,'BKVL',VLSoilMicPMass_vr(L,NY,NX))
    call ncd_putvar(ncf,'POROS',POROS_vr(L,NY,NX))
    call ncd_putvar(ncf,'FC',FieldCapacity_vr(L,NY,NX))
    call ncd_putvar(ncf,'ATKA',TairKClimMean_col(NY,NX))
    call ncd_putvar(ncf,'WP',WiltPoint_vr(L,NY,NX))
    call ncd_putvar(ncf,'SRP',SRP_vr(L,NY,NX))
    call ncd_putvar(ncf,'EHUM',EHUM_vr(L,NY,NX))
    call ncd_putvar(ncf,'EPOC',EPOC_vr(L,NY,NX))
    call ncd_putvar(ncf,'THETY',SoilWatAirDry_vr(L,NY,NX))
    call ncd_putvar(ncf,'PSIMX',LOGPSIFLD_col(NY,NX))
    call ncd_putvar(ncf,'PSIMD',LOGPSIMND_col(NY,NX))
    call ncd_putvar(ncf,'PSIMS',LOGPSIAtSat(NY,NX))
    call ncd_putvar(ncf,'PSISD',LOGPSIMXD_col(NY,NX))
    call ncd_putvar(ncf,'PSISE',PSI_at_Sautration_vr(L,NY,NX))
    call ncd_putvar(ncf,'DLYR3',DLYR_3D(3,L,NY,NX))
    call ncd_putvar(ncf,'CEC',CEC_vr(L,NY,NX))
    call ncd_putvar(ncf,'XCEC',trcx_solml_vr(idx_CEC,L,NY,NX))
    call ncd_putvar(ncf,'AEC',AEC_vr(L,NY,NX))
    call ncd_putvar(ncf,'XAEC',trcx_solml_vr(idx_AEC,L,NY,NX))
    call ncd_putvar(ncf,'CFE',CFE_vr(L,NY,NX))
    call ncd_putvar(ncf,'CCA',CCA_vr(L,NY,NX))
    call ncd_putvar(ncf,'CMG',CMG_vr(L,NY,NX))
    call ncd_putvar(ncf,'CNA',CNA_vr(L,NY,NX))
    call ncd_putvar(ncf,'CKA',CKA_vr(L,NY,NX))
    call ncd_putvar(ncf,'CSO4',CSO4_vr(L,NY,NX))
    call ncd_putvar(ncf,'CCL',CCL_vr(L,NY,NX))
    call ncd_putvar(ncf,'CAL',CAL_vr(L,NY,NX))
    call ncd_putvar(ncf,'ZMG',trcSalt_solml_vr(idsalt_Mg,L,NY,NX))
    call ncd_putvar(ncf,'ZNA',trcSalt_solml_vr(idsalt_Na,L,NY,NX))
    call ncd_putvar(ncf,'ZKA',trcSalt_solml_vr(idsalt_K,L,NY,NX))
    call ncd_putvar(ncf,'CALPO',CALPO_vr(L,NY,NX))
    call ncd_putvar(ncf,'CFEPO',CFEPO_vr(L,NY,NX))
    call ncd_putvar(ncf,'CCAPD',CCAPD_vr(L,NY,NX))
    call ncd_putvar(ncf,'CCAPH',CCAPH_vr(L,NY,NX))
    call ncd_putvar(ncf,'CALOH',CALOH_vr(L,NY,NX))
    call ncd_putvar(ncf,'CFEOH',CFEOH_vr(L,NY,NX))
    call ncd_putvar(ncf,'CCACO',CCACO_vr(L,NY,NX))
    call ncd_putvar(ncf,'CCASO',CCASO_vr(L,NY,NX))
    call ncd_putvar(ncf,'BKDS',SoilBulkDensity_vr(L,NY,NX))

    call ncd_putvar(ncf,'FracBulkSOMC_vr',FracBulkSOMC_vr(:,L,NY,NX))
    call ncd_putvar(ncf,'OQM',DOM_MicP_vr(:,:,L,NY,NX))
    call ncd_putvar(ncf,'OHM',SorbedOM_vr(:,:,L,NY,NX))

    call ncd_putvar(ncf,'CNOSC',CNOSC_vr(:,:,L,NY,NX))
    call ncd_putvar(ncf,'CPOSC',CPOSC_vr(:,:,L,NY,NX))
    call ncd_putvar(ncf,'ATCS',ATCS_col(NY,NX))
    call ncd_putvar(ncf,'OSM',SolidOM_vr(:,:,:,L,NY,NX))
    call ncd_putvar(ncf,'OSA',SolidOMAct_vr(:,:,L,NY,NX))
    call ncd_putvar(ncf,'ORM',OMBioResdu_vr(:,:,:,L,NY,NX))
    call ncd_putvar(ncf,'mBiomeHeter',mBiomeHeter_vr(:,:,:,L,NY,NX))
    call ncd_putvar(ncf,'mBiomeAutor',mBiomeAutor_vr(:,:,L,NY,NX))

    if(bgc_forc_conf%laddband)then
      call ncd_putvar(ncf,'ZNH4S',trcs_solml_vr(ids_NH4,L,NY,NX)+trcs_solml_vr(ids_NH4B,L,NY,NX))
      call ncd_putvar(ncf,'ZNO3S',trcs_solml_vr(ids_NO3,L,NY,NX)+trcs_solml_vr(ids_NO3B,L,NY,NX))
      call ncd_putvar(ncf,'ZNO2S',trcs_solml_vr(ids_NO2,L,NY,NX)+trcs_solml_vr(ids_NO2B,L,NY,NX))
      call ncd_putvar(ncf,'H2PO4',trcs_solml_vr(ids_H2PO4,L,NY,NX)+trcs_solml_vr(ids_H2PO4B,L,NY,NX))
      call ncd_putvar(ncf,'H1PO4',trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H1PO4B,L,NY,NX))
    else
      call ncd_putvar(ncf,'ZNH4S',trcs_solml_vr(ids_NH4,L,NY,NX))
      call ncd_putvar(ncf,'ZNO3S',trcs_solml_vr(ids_NO3,L,NY,NX))
      call ncd_putvar(ncf,'ZNO2S',trcs_solml_vr(ids_NO2,L,NY,NX))
      call ncd_putvar(ncf,'H2PO4',trcs_solml_vr(ids_H2PO4,L,NY,NX))
      call ncd_putvar(ncf,'H1PO4',trcs_solml_vr(ids_H1PO4,L,NY,NX))
    endif
    call ncd_pio_closefile(ncf)
  endif
  end subroutine WriteBBGCForc

end module ForcWriterMod
