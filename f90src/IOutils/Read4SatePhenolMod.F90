module Read4SatePhenolMod
!
! Description:
! Read the PFT distribution and monthly structure prescribed in satellite phenology mode
! (ldo_sp_mode) from the NetCDF file sp_phenol_in.
!
! File layout (topographic units and PFT naming as in pft_mgmt_in):
!   dimensions: ntopou, maxpfts, nchar1 (=10), month (=12), year (unlimited)
!   int   NH1(ntopou), NV1(ntopou), NH2(ntopou), NV2(ntopou) : topo-unit bounds
!   int   NZ(ntopou)                        : tracked PFT slots (NP0_col)
!   int   year(year)                         : calendar year of each record
!   char  pft_type(year,ntopou,maxpfts,nchar1): PFT names; trailing blanks are fallow
!   double PlantPopuLive_pft(year,ntopou,maxpfts): population density, [plants m-2]
!   float LAI(year,ntopou,maxpfts,month)     : leaf area index, [m2 m-2]
!   float SAI(year,ntopou,maxpfts,month)     : stem area index, [m2 m-2]
!   float HTOP(year,ntopou,maxpfts,month)    : canopy top height, [m]
! A named slot with zero living population is standing dead only: LAI must
! be zero, SAI prescribes dead stem area, and HTOP is dead stem height.
! Its initial dead C/N/P comes from WTSTDI and the PFT litter composition.
! A file with a single year record is a climatology applied to every year;
! otherwise the record matching the forcing year is used.
!
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use netcdf,           only: nf90_get_var
  use abortutils,        only: endrun, iulog
  use DebugToolMod,      only: PrintInfo
  use ncdio_pio
  use GridDataType,      only: AREA_3D, NU_col
  use GridConsts,        only: JX, JY, JP
  use EcoSIMCtrlMod,     only: sp_phenol_in, first_topou, first_pft
  use EcoSIMHistMod,     only: DATAP, DATAPI
  use FlagDataType,      only: KoppenClimZone_col, IsPlantActive_pft
  use ElmIDMod,         only: iFalse
  use PlantTraitDataType, only: PlantPopuLive_pft, LeafStalkAreaAct_pft
  use PlantTraitTableMod, only: pftss_tab
  use PlantMgmtDataType, only: NP_col, NP0_col, LSG_pft, flag_active_pft
  use ClimForcDataType,  only: tlai_mon_pft, tsai_mon_pft, height_top_mon_pft, sp_popdens_pft
implicit none
  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: ReadSatePhenol
  contains

!------------------------------------------------------------------------------------------
  subroutine ReadSatePhenol(yeari,NHW,NHE,NVN,NVS)
  ! Description:
  ! Load PFT names, population densities, slot counts and monthly canopy structure for forcing year yeari.
  implicit none
  integer, intent(in) :: yeari    !forcing calendar year, [AD]
  integer, intent(in) :: NHW,NHE  !inclusive west/east column indices, [-]
  integer, intent(in) :: NVN,NVS  !inclusive north/south row indices, [-]

  character(len=*), parameter :: subname='ReadSatePhenol'
  type(file_desc_t) :: sp_nfid
  type(var_desc_t) :: vardesc
  integer :: ntopou,npfts,nmonths,nyears
  integer :: irec,iyear,year
  integer :: NTOPO,NH1,NV1,NH2,NV2,NX,NY,NZ,NS,NA,kk
  logical :: readvar
  character(len=2) :: climate
  character(len=10), allocatable :: pft_type(:)
  real(r8), allocatable :: popdens(:,:) !(maxpfts,ntopou), [plants m-2]
  real(r8), allocatable :: lai(:,:,:)    !(month,maxpfts,ntopou), [m2 m-2]
  real(r8), allocatable :: sai(:,:,:)    !(month,maxpfts,ntopou), [m2 m-2]
  real(r8), allocatable :: htop(:,:,:)   !(month,maxpfts,ntopou), [m]

  call PrintInfo('beg '//subname)
  if(len_trim(sp_phenol_in)==0)then
    call endrun('sp_phenol_in must be set when ldo_sp_mode is on in '//mod_filename,__LINE__)
  endif

  call ncd_pio_openfile(sp_nfid, sp_phenol_in, ncd_nowrite)

  ntopou  = get_dim_len(sp_nfid, 'ntopou')
  npfts   = get_dim_len(sp_nfid, 'maxpfts')
  nmonths = get_dim_len(sp_nfid, 'month')
  nyears  = get_dim_len(sp_nfid, 'year')
  if(npfts<1 .or. nyears<1 .or. get_dim_len(sp_nfid,'nchar1')/=10)then
    call endrun('sp_phenol_in requires PFT/year records and nchar1=10 in '//mod_filename,__LINE__)
  endif
  if(nmonths/=12)then
    write(iulog,*)subname//': month dimension is ',nmonths,' in ',trim(sp_phenol_in)
    call endrun('month dimension must be 12 in '//mod_filename,__LINE__)
  endif

  !select the record: a single record is a climatology used for every year
  if(nyears==1)then
    irec=1
  else
    irec=0
    DO iyear=1,nyears
      call ncd_getvar(sp_nfid,'year',iyear,year)
      if(year==yeari)then
        irec=iyear
        exit
      endif
    ENDDO
    if(irec==0)then
      write(iulog,*)subname//': year ',yeari,' not found in ',trim(sp_phenol_in)
      call endrun('missing forcing year in '//mod_filename,__LINE__)
    endif
  endif

  call check_var(sp_nfid,'PlantPopuLive_pft',vardesc,readvar)
  if(.not.readvar)call endrun('missing PlantPopuLive_pft in '//mod_filename,__LINE__)
  allocate(popdens(npfts,ntopou))
  call ncd_getvar(sp_nfid,'PlantPopuLive_pft',irec,popdens)
  allocate(lai(12,npfts,ntopou),sai(12,npfts,ntopou),htop(12,npfts,ntopou))
  call ncd_getvar(sp_nfid,'LAI',irec,lai)
  call ncd_getvar(sp_nfid,'SAI',irec,sai)
  call ncd_getvar(sp_nfid,'HTOP',irec,htop)
  allocate(pft_type(npfts))

  ! This file owns the standalone SP distribution, including bare columns.
  sp_popdens_pft(:,NVN:NVS,NHW:NHE)       = 0._r8
  PlantPopuLive_pft(:,NVN:NVS,NHW:NHE)    = 0._r8
  NP_col(NVN:NVS,NHW:NHE)                 = 0
  NP0_col(NVN:NVS,NHW:NHE)                = 0
  DATAP(:,NVN:NVS,NHW:NHE)                = 'NO'
  flag_active_pft(:,NVN:NVS,NHW:NHE)      = .false.
  tlai_mon_pft(:,:,NVN:NVS,NHW:NHE)       = 0._r8
  tsai_mon_pft(:,:,NVN:NVS,NHW:NHE)       = 0._r8
  height_top_mon_pft(:,:,NVN:NVS,NHW:NHE) = 0._r8

  if(first_topou)ntopou=1
  DO NTOPO=1,ntopou
    call ncd_getvar(sp_nfid,'NH1',NTOPO,NH1)
    call ncd_getvar(sp_nfid,'NV1',NTOPO,NV1)
    call ncd_getvar(sp_nfid,'NH2',NTOPO,NH2)
    call ncd_getvar(sp_nfid,'NV2',NTOPO,NV2)
    !nonpositive bounds end the populated topo units, as in pft_mgmt_in
    if(any((/NH1,NV1,NH2,NV2/)<=0))exit
    call ncd_getvar(sp_nfid,'NZ',NTOPO,NS)
    if(NS<0 .or. NS>JP .or. NS>npfts)then
      write(iulog,*)subname//': topo unit ',NTOPO,' NZ=',NS,' model/file capacity=',JP,npfts
      call endrun('invalid PFT slot count in '//mod_filename,__LINE__)
    endif
    if(NH1>NH2 .or. NH2>JX .or. NV1>NV2 .or. NV2>JY .or. &
       NH1<NHW .or. NH2>NHE .or. NV1<NVN .or. NV2>NVS)then
      write(iulog,*)subname//': topo unit ',NTOPO,' NH1,NH2,NV1,NV2=',NH1,NH2,NV1,NV2
      call endrun('invalid topo-unit bounds in '//mod_filename,__LINE__)
    endif

    call check_var(sp_nfid,'pft_type',vardesc,readvar)
    if(.not.readvar)call endrun('missing pft_type in '//mod_filename,__LINE__)
    call check_ret(nf90_get_var(sp_nfid%fh,vardesc%varid,pft_type, &
      start=(/1,1,NTOPO,irec/),count=(/10,npfts,1,1/)),subname//': reading pft_type')
    ! Normalize NetCDF NUL padding before testing names or adding a climate suffix.
    DO NZ=1,NS
      DO kk=1,len(pft_type(NZ))
        if(iachar(pft_type(NZ)(kk:kk))==0)pft_type(NZ)(kk:kk)=' '
      ENDDO
    ENDDO
    NA=0
    DO NZ=1,NS
      if(len_trim(pft_type(NZ))==0)exit
      NA=NZ
    ENDDO
    DO NZ=NA+1,NS
      if(len_trim(pft_type(NZ))>0)then
        call endrun('blank pft_type slots must follow named slots in '//mod_filename,__LINE__)
      endif
    ENDDO
    if(first_pft)then
      NS=MIN(NS,1)
      NA=MIN(NA,1)
    endif

    DO NX=NH1,NH2
      DO NY=NV1,NV2
        NP_col(NY,NX)               = NA
        NP0_col(NY,NX)              = NS
        LSG_pft(1:NS,NY,NX)         = 1
        flag_active_pft(1:NA,NY,NX) = .true.
        DO NZ=1,NA
          DATAP(NZ,NY,NX)=pft_type(NZ)
          if(KoppenClimZone_col(NY,NX)>0)then
            write(climate,'(I2)')KoppenClimZone_col(NY,NX)
            DATAP(NZ,NY,NX)=pft_type(NZ)(1:4)//climate
          endif
          if(any(lai(:,NZ,NTOPO)<0._r8) .or. any(sai(:,NZ,NTOPO)<0._r8) .or. any(htop(:,NZ,NTOPO)<0._r8) &
            .or. any(lai(:,NZ,NTOPO)>1.e10_r8) .or. any(sai(:,NZ,NTOPO)>1.e10_r8) .or. any(htop(:,NZ,NTOPO)>1.e10_r8) &
            .or. any(.not.ieee_is_finite(lai(:,NZ,NTOPO))) .or. any(.not.ieee_is_finite(sai(:,NZ,NTOPO))) &
            .or. any(.not.ieee_is_finite(htop(:,NZ,NTOPO))))then
            write(iulog,*)subname//': negative or missing value, topo unit ',NTOPO,' PFT ',NZ
            call endrun('bad LAI/SAI/HTOP in '//mod_filename,__LINE__)
          endif
          if(.not.ieee_is_finite(popdens(NZ,NTOPO)))then
            write(iulog,*)subname//': nonfinite population density, topo/PFT=',NTOPO,NZ
            call endrun('prescribed PlantPopuLive_pft must be finite and nonnegative in '//mod_filename,__LINE__)
          endif
          if(popdens(NZ,NTOPO)<0._r8 .or. popdens(NZ,NTOPO)>1.e10_r8)then
            write(iulog,*)subname//': invalid population density, topo/PFT/value=',NTOPO,NZ,popdens(NZ,NTOPO)
            call endrun('prescribed PlantPopuLive_pft must be finite and nonnegative in '//mod_filename,__LINE__)
          endif
          if(popdens(NZ,NTOPO)==0._r8)then
            if(any(lai(:,NZ,NTOPO)/=0._r8) .or. .not.any(sai(:,NZ,NTOPO)>0._r8) &
              .or. any(sai(:,NZ,NTOPO)>0._r8 .and. htop(:,NZ,NTOPO)<=0._r8))then
              write(iulog,*)subname//': invalid standing-dead geometry, topo/PFT=',NTOPO,NZ
              call endrun('zero-population SP slots require LAI=0 and positive dead SAI/HTOP',__LINE__)
            endif
          endif
          sp_popdens_pft(NZ,NY,NX)=popdens(NZ,NTOPO)
          PlantPopuLive_pft(NZ,NY,NX)=popdens(NZ,NTOPO)*AREA_3D(3,NU_col(NY,NX),NY,NX)
          tlai_mon_pft(1:12,NZ,NY,NX)       = lai(1:12,NZ,NTOPO)
          tsai_mon_pft(1:12,NZ,NY,NX)       = sai(1:12,NZ,NTOPO)
          height_top_mon_pft(1:12,NZ,NY,NX) = htop(1:12,NZ,NTOPO)
        ENDDO
        DO NZ=NA+1,NS
          ! Keep the last species' traits for a retained fallow/residue slot.
          if(DATAPI(NZ,NY,NX)>0)DATAP(NZ,NY,NX)=pftss_tab(DATAPI(NZ,NY,NX))
        ENDDO
        ! A blank prescribed slot has no living canopy, including when a
        ! previous year's canopy or a January checkpoint still contains one.
        PlantPopuLive_pft(NA+1:JP,NY,NX)=0._r8
        LeafStalkAreaAct_pft(NA+1:JP,NY,NX)=0._r8
        IsPlantActive_pft(NA+1:JP,NY,NX)=iFalse
      ENDDO
    ENDDO
  ENDDO

  call ncd_pio_closefile(sp_nfid)
  deallocate(lai,sai,htop,pft_type,popdens)
  call PrintInfo('end '//subname)
  end subroutine ReadSatePhenol

end module Read4SatePhenolMod
