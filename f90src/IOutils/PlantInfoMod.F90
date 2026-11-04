module PlantInfoMod
!
! DESCRIPTION
! code to read plant information
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use fileUtil,      only: open_safe, check_read, int2str, getavu, relavu, opnfil
  use minimathmod,   only: isLeap
  use abortutils,    only: endrun, iulog
  use DebugToolMod,  only: PrintInfo
  use StringToolsMod, only: count_delimited_items, string_compaction
  use PlantTraitTableMod
  use netcdf
  use ncdio_pio  
  use GridConsts
  use FlagDataType
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use CanopyDataType
  use PlantTraitDataType
  use PlantMathFuncMod
  use PlantMgmtDataType
  use EcosimConst
  use RootDataType
  use EcoSIMHistMod
  use CanopyRadDataType
  use GridDataType
  use EcoSIMCtrlMod
  use EcoSIMConfig

implicit none

  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  public :: ReadPlantInfo
  public :: ReadPlantTraitTable
  contains

  subroutine ReadPlantInfo(yearc,yeari,NHW,NHE,NVN,NVS)
!
! DESCRIPTION
!

  implicit none
  integer, intent(in) :: yearc, yeari
  integer, intent(in) :: NHW,NHE,NVN,NVS     !simulation domain specification

! RECOVER PLANT SPECIES DISTRIBUTION IN 'ROUTQ'
!
 
  CALL ReadPlantInfoNC(yeari,NHW,NHE,NVN,NVS)
!
!   READ INPUT DATA FOR PLANT SPECIES AND MANAGEMENT IN 'READQ'
!   AND SET UP OUTPUT AND CHECKPOINT FILES IN 'FOUTP'
!
  CALL READQ(yearc,yeari,NHW,NHE,NVN,NVS)

333   FORMAT(A8)

  end subroutine ReadPlantInfo
!------------------------------------------------------------------------------------------

  SUBROUTINE readq(yearc,yeari,NHW,NHE,NVN,NVS)
!!
! Description
! THIS SUBROUTINE READS INPUT DATA FROM PLANT SPECIES
!     AND MANAGEMENT FILES IDENTIFIED IN 'ROUTQ'
!
  implicit none
  integer, intent(in) :: yearc   !current model year
  integer, intent(in) :: yeari   !current forcing year 
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NX,NY,NZ,nu_plt
  character(len=*), parameter :: subname='READQ'
  character(len=128) :: fnm_loc
  logical :: pft_changed

! begin_execution
  call PrintInfo('beg '//subname)
  call ReadPlantManagementNC(yearc,yeari,NHW,NHE,NVN,NVS,pft_changed)

  if(disp_planttrait .and. pft_changed)then
    nu_plt=getavu()
    write(fnm_loc,'(A,I4)')trim(case_name)//'.plant_trait.desc.',etimer%get_curr_yearAD()
    call opnfil(fnm_loc,nu_plt,'f')    
  endif

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      !fallow slots NP_col+1..NP0_col reload the traits of their residue PFT
      D9985: DO NZ=1,NP0_col(NY,NX)
!
!       OPEN PFT(11), PFT MANAGEMENT(12) FILES FROM
!       FILE NAMES IN DATA ARRAYS LOADED IN MAIN.F
!       each column has its own management
!       PREFIX=path for files in current or higher level directory
!       DATAP=PFT file name
        call ReadPlantProperties(nu_plt,NZ,NY,NX,pft_changed)

      ENDDO D9985
    ENDDO D9990
  ENDDO D9995

  if(disp_planttrait .and. pft_changed)call relavu(nu_plt)
  call PrintInfo('end '//subname)
  
  RETURN
  END SUBROUTINE readq
!------------------------------------------------------------------------------------------
  subroutine readplantinginfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
  !
  !read planting information
  implicit none
  type(file_desc_t), intent(in) :: pftinfo_nfid
  integer, intent(in) :: ntopou
  integer, intent(in) :: iyear    !file record number
  integer, intent(in) :: yearc    !current model year
  integer, intent(in) :: NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='readplantinginfo'
  logical  :: readvar
  integer  :: NH1,NV1,NH2,NV2,NS
  integer :: NTOPO,NY,NX,NZ,lstr
  type(Var_desc_t) :: vardesc
  character(len=128) :: pft_pltinfo(JP),tstr,date_str
  integer :: IYR,IDY,ierr

  call PrintInfo('beg '//subname)
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO NZ=1,NP_col(NY,NX)
        PlantinDepz_pft(NZ,NY,NX)=1.e-6_r8
        PPmax_pft(NZ,NY,NX)=0._r8
      ENDDO
    ENDDO
  ENDDO

  DO NTOPO=1,ntopou
    call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
    call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NV1)
    call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH2)
    call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NV2)
    call ncd_getvar(pftinfo_nfid,'NZ',ntopo,NS)
    !
    ! Legacy plant files may reserve more ntopou slots than they populate.
    ! Nonpositive bounds (including NetCDF fill values) end the populated records;
    ! exit before validating NZ or reading planting/management data for unused slots.
    ! Keep this sentinel exit when changing validation. Populated units may span
    ! multiple cells, so their bounds are checked against JX/JY, not ntopou.
    !
    if(any((/NH1,NV1,NH2,NV2/)<=0))exit
    if(NS<0 .or. NS>JP)then
      write(iulog,*)'Invalid PFT count in '//subname,': topo unit=',NTOPO,' NZ=',NS,' capacity=',JP
      call endrun('PFT count must be between zero and JP in '//subname,__LINE__)
    endif
    !
    ! Topographic units may span several cells; validate against grid dimensions.
    if(NH1<1 .or. NH2>JX .or. NH1>NH2 .or. &
       NV1<1 .or. NV2>JY .or. NV1>NV2)then
      write(iulog,*)'Invalid bounds in '//subname,': topo unit=',NTOPO, &
        ' NH1,NH2,NV1,NV2=',NH1,NH2,NV1,NV2,' JX,JY=',JX,JY
      call endrun('Invalid topographic-unit bounds in '//subname,__LINE__)
    endif
    !planting information
    call check_var(pftinfo_nfid, 'pft_pltinfo', vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find pft_pltinfo in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_pltinfo, &
      start = (/1,1,NTOPO,iyear/),count = (/len(pft_pltinfo(1)),JP,1,1/)),&
      trim(mod_filename)//':: at line '//trim(int2str(__LINE__)))

    DO NX=NH1,NH2
      DO NY=NV1,NV2
        DO NZ=1,MIN(NS,NP_col(NY,NX))
          if(PlantPopuLive_pft(NZ,NY,NX).GT.1.e-2_r8)cycle          
          tstr=pft_pltinfo(NZ)
          if (string_compaction(tstr) .EQ. "") then
            cycle
          endif

          call CheckPlantRecordFields(tstr,3,4,NTOPO,NZ)
          date_str=''
          select case(count_delimited_items(tstr, ' '))
          case(3)
            lstr=3
            read(tstr,*,iostat=ierr)date_str,PPI_pft(NZ,NY,NX),PlantinDepz_pft(NZ,NY,NX)
          case(4)
            lstr=4
            read(tstr,*,iostat=ierr)date_str,PPI_pft(NZ,NY,NX),PlantinDepz_pft(NZ,NY,NX),PPmax_pft(NZ,NY,NX)
          case default
            lstr=0
            ierr=1
          end select
          if(ierr/=0)then
            write(iulog,*)'lstr=',lstr
            call endrun('Invalid planting record (expected DDMMYYYY population depth [maximum population]): '// &
              trim(tstr)//' in '//subname,__LINE__)
          endif
          PlantinDepz_pft(NZ,NY,NX)=AZMAX1(PlantinDepz_pft(NZ,NY,NX),1.e-6_r8)
          call ParsePlantEventDate(date_str,yearc,NTOPO,NZ,IDY,IYR)

          IF(IDY.GT.0.AND.IYR.GT.0)THEN
            IF(iselect_plantZ==-1 .OR. (iselect_plantZ>0 .and. NZ==iselect_plantZ))THEN
              iDayPlanting_pft(NZ,NY,NX)  = IDY
              IYR                         = yearc
              iYearPlanting_pft(NZ,NY,NX) = MIN(IYR,iYearCurrent)
              iPlantingDay_pft(NZ,NY,NX)  = iDayPlanting_pft(NZ,NY,NX)    !planting day
              iPlantingYear_pft(NZ,NY,NX) = iYearPlanting_pft(NZ,NY,NX)   !planting year
              PPatSeeding_pft(NZ,NY,NX)   = PPI_pft(NZ,NY,NX)             !population density              
            ENDIF           
          ENDIF          
        ENDDO
      ENDDO
    ENDDO
  ENDDO  
  CALL PrintInfo('end '//subname)
  end subroutine readplantinginfo

!------------------------------------------------------------------------------------------
  subroutine CheckPlantRecordFields(record,MinFields,MaxFields,NTOPO,NZ)
  ! Require explicit values before list-directed reads can retain old values.
  implicit none
  character(len=*), intent(inout) :: record
  integer, intent(in) :: MinFields,MaxFields,NTOPO,NZ
  character(len=*), parameter :: subname='CheckPlantRecordFields'
  character(len=len(record)) :: normalized
  integer :: i,first,nfields
  logical :: valid_record

  call PrintInfo('beg '//subname)
  normalized=string_compaction(record)
  do i=1,len(normalized)
    if(normalized(i:i)==achar(9))normalized(i:i)=' '
  enddo
  ! Slash terminates a read successfully; asterisk permits repeated/null values.
  valid_record=index(record,'/')==0.and.index(record,'*')==0
  first=1
  do i=1,len_trim(record)
    if(record(i:i)==',')then
      ! Check comma fields before normalization hides empty fields as whitespace.
      if(len_trim(normalized(first:i-1))==0)valid_record=.false.
      first=i+1
    endif
  enddo
  if(len_trim(normalized(first:))==0)valid_record=.false.
  nfields=count_delimited_items(normalized,' ')
  if(nfields<MinFields.or.nfields>MaxFields)valid_record=.false.
  if(.not.valid_record)then
    write(iulog,*)'Invalid plant record fields: topo unit=',NTOPO,' PFT=',NZ, &
      ' expected fields=',MinFields,MaxFields,' record=',trim(record)
    call endrun('Plant records require explicit values for every field in '//subname,__LINE__)
  endif
  record=normalized
  call PrintInfo('end '//subname)
  end subroutine CheckPlantRecordFields

!------------------------------------------------------------------------------------------
  subroutine ParsePlantEventDate(date_str,yearc,NTOPO,NZ,IDY,IYR)
  ! Validate DDMMYYYY before computing any calendar or daily-array index.
  implicit none
  character(len=*), intent(in) :: date_str
  integer, intent(in) :: yearc,NTOPO,NZ
  integer, intent(out) :: IDY,IYR
  character(len=*), parameter :: subname='ParsePlantEventDate'
  integer :: IDX,IMO,ierr
  integer :: MonthDays(12)
  logical :: valid_date

  call PrintInfo('beg '//subname)
  MonthDays=(/31,28,31,30,31,30,31,31,30,31,30,31/)
  ! Events are applied in yearc, including recurring management with year 0000.
  ! Preserve the input year for the caller's planting/harvest selection rules.
  if(isLeap(yearc))MonthDays(2)=29
  valid_date=len_trim(date_str)==8.and.verify(trim(date_str),'0123456789')==0
  if(valid_date)then
    read(date_str,'(I2,I2,I4)',iostat=ierr)IDX,IMO,IYR
    valid_date=ierr==0
  endif

  if(valid_date)valid_date=IMO>=1.and.IMO<=12
  ! Do not index MonthDays until the month has been validated.
  if(valid_date)valid_date=IDX>=1.and.IDX<=MonthDays(IMO)
  if(.not.valid_date)then
    write(iulog,*)'Invalid plant event date: topo unit=',NTOPO,' PFT=',NZ, &
      ' date=',trim(date_str),' simulation year=',yearc
    call endrun('Expected a valid DDMMYYYY date in '//subname,__LINE__)
  endif
  IDY=sum(MonthDays(1:IMO-1))+IDX
  call PrintInfo('end '//subname)
  end subroutine ParsePlantEventDate

!------------------------------------------------------------------------------------------
  subroutine InitPlantMgmnt(NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX,NZ,M

  ! initialize the disturbance arrays
  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO NZ=1,NP0_col(NY,NX)
        DO M=1,366
          iHarvstType_pft(NZ,M,NY,NX)                                   = -1
          jHarvstType_pft(NZ,M,NY,NX)                                   = 0     !do not terminate
          CanopyCutProxy_pft(NZ,M,NY,NX)                               = 1.0E+06_r8
          THIN_pft(NZ,M,NY,NX)                                          = -1.0_r8
          FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ,M,NY,NX)        = 1.0_r8
          FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ,M,NY,NX) = 1.0_r8
          FracBiomHarvsted(iHarvst_pft,iplthvst_stalk,NZ,M,NY,NX)       = 1.0_r8
          FracBiomHarvsted(iHarvst_pft,iplthvst_stdead,NZ,M,NY,NX)      = 1.0_r8
          FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ,M,NY,NX)        = 1.0_r8
          FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ,M,NY,NX) = 1.0_r8
          FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ,M,NY,NX)       = 1.0_r8
          FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ,M,NY,NX)      = 1.0_r8
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  end subroutine InitPlantMgmnt
!------------------------------------------------------------------------------------------
  subroutine ReadPlantMgmtinfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
  !
  !read plant management information
  implicit none
  type(file_desc_t), intent(in) :: pftinfo_nfid
  integer, intent(in) :: ntopou
  integer, intent(in) :: iyear    !file record number
  integer, intent(in) :: yearc    !current model year  
  integer, intent(in) :: NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='ReadPlantMgmtinfo'
  logical  :: readvar
  real(r8) :: ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23
  real(r8) :: ECUT24,HCUT,PCUT
  integer  :: pft_nmgnt(JP)
  integer  :: NH1,NV1,NH2,NV2,NS
  type(Var_desc_t) :: vardesc
  integer :: NTOPO,NY,NX,NZ
  character(len=128) :: pft_pltinfo(JP),tstr,date_str
  character(len=128) :: pft_mgmtinfo(24,JP)
  integer :: IYR,IDY,ICUT,IDYE,IDYG,JCUT,IDYS,ierr
  integer :: M,NN,N,nn1

  call PrintInfo('beg '//subname)
  call InitPlantMgmnt(NHW,NHE,NVN,NVS)
  
  DO NTOPO=1,ntopou
    call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
    call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NV1)
    call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH2)
    call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NV2)
    call ncd_getvar(pftinfo_nfid,'NZ',ntopo,NS)
    ! Legacy plant files may reserve more ntopou slots than they populate.
    ! Nonpositive bounds (including NetCDF fill values) end the populated records;
    ! exit before validating NZ or reading planting/management data for unused slots.
    ! Keep this sentinel exit when changing validation. Populated units may span
    ! multiple cells, so their bounds are checked against JX/JY, not ntopou.
    if(any((/NH1,NV1,NH2,NV2/)<=0))exit
    if(NS<0 .or. NS>JP)then
      write(iulog,*)'Invalid PFT count in '//subname,': topo unit=',NTOPO,' NZ=',NS,' capacity=',JP
      call endrun('PFT count must be between zero and JP in '//subname,__LINE__)
    endif

    call check_var(pftinfo_nfid, 'nmgnts', vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find pft_type in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_nmgnt, &
      start = (/1,NTOPO,iyear/),count = (/JP,1,1/)),trim(mod_filename)// &
      '::at line '//trim(int2str(__LINE__)))

    do NZ=1,NS
      if(pft_nmgnt(NZ)<0 .or. pft_nmgnt(NZ)>size(pft_mgmtinfo,1))then
        write(iulog,*)'Invalid event count in '//subname,': topo unit=',NTOPO, &
          ' PFT=',NZ,' nmgnts=',pft_nmgnt(NZ),' capacity=',size(pft_mgmtinfo,1)
        call endrun('Plant management event count exceeds the supported range in '//subname,__LINE__)
      endif
    enddo

    if(any(pft_nmgnt(1:NS)>0))then
      call check_var(pftinfo_nfid, 'pft_mgmt', vardesc, readvar)
      if(.not. readvar)then
        call endrun('fail to find pft_mgmt in '//trim(mod_filename), __LINE__)
      endif

      call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_mgmtinfo, &
        start = (/1,1,1,NTOPO,iyear/),count = (/len(pft_mgmtinfo(1,1)),24,JP,1,1/)),&
        trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
    endif
    ! Topographic units may span several cells; validate against grid dimensions.
    if(NH1<1 .or. NH2>JX .or. NH1>NH2 .or. &
       NV1<1 .or. NV2>JY .or. NV1>NV2)then
      write(iulog,*)'Invalid bounds in '//subname,': topo unit=',NTOPO, &
        ' NH1,NH2,NV1,NV2=',NH1,NH2,NV1,NV2,' JX,JY=',JX,JY
      call endrun('Invalid topographic-unit bounds in '//subname,__LINE__)
    endif

    DO NX=NH1,NH2
      DO NY=NV1,NV2
        DO NZ=1,MIN(NS,NP_col(NY,NX))          
          if(pft_nmgnt(NZ)>0)then
            NN=0
            DO nn1=1,pft_nmgnt(NZ)
              if(string_compaction(pft_mgmtinfo(NN1,NZ))=='')cycle
              tstr=trim(pft_mgmtinfo(NN1,NZ))
              call CheckPlantRecordFields(tstr,13,13,NTOPO,NZ)
              date_str=''
              read(tstr,*,iostat=ierr)date_str,ICUT,JCUT,HCUT,PCUT,ECUT11,ECUT12,ECUT13, &
                  ECUT14,ECUT21,ECUT22,ECUT23,ECUT24
              if(ierr/=0)then
                write(iulog,*)'Invalid plant management record: topo unit=',NTOPO,' PFT=',NZ, &
                  ' record=',trim(tstr)
                call endrun('Invalid management record in '//subname,__LINE__)
              endif
              call ParsePlantEventDate(date_str,yearc,NTOPO,NZ,IDY,IYR)

              IF(IDY.GT.0 .AND. JCUT.EQ.1)THEN
                iDayPlantHarvest_pft(NZ,NY,NX)=IDY
                IYR=yearc
                iYearPlantHarvest_pft(NZ,NY,NX)=MIN(IYR,iYearCurrent)
              ENDIF
              !     HVST=iHarvstType_pft=0-2:>0=cutting height,<0=fraction of LAI removed
              !          iHarvstType_pft=3:reduction of clumping factor, pruning
              !          iHarvstType_pft=4 or 6:animal or insect biomass(g LM m-2),
              !          iHarvstType_pft=5:fire
              !     THIN_pft=iHarvstType_pft=0-3,5: fraction of population removed,
              !          iHarvstType_pft=4 or 6:specific herbivory rate (g DM g-1 LM d-1)
              !
              !harvest is specified with two type numbers, first is large categroy, second is pft-specific operation
              iHarvstType_pft(NZ,IDY,NY,NX)    = ICUT
              jHarvstType_pft(NZ,IDY,NY,NX)    = JCUT
              CanopyCutProxy_pft(NZ,IDY,NY,NX) = HCUT
              THIN_pft(NZ,IDY,NY,NX)           = PCUT  
              
              !Note: pft-level-harvest minus ecosystem-level-harvest = litter to soil 
              !pft-level harvest
              FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ,IDY,NY,NX)        = ECUT11
              FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ,IDY,NY,NX) = ECUT12
              FracBiomHarvsted(iHarvst_pft,iplthvst_stalk,NZ,IDY,NY,NX)       = ECUT13
              FracBiomHarvsted(iHarvst_pft,iplthvst_stdead,NZ,IDY,NY,NX)      = ECUT14

              !ecosystem-level harvest
              !For a harvested pft component, a fraction by (1- ecosystem-level harvest fraction) becomes litter
              FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ,IDY,NY,NX)        = ECUT21
              FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ,IDY,NY,NX) = ECUT22
              FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ,IDY,NY,NX)       = ECUT23
              FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ,IDY,NY,NX)      = ECUT24

              IF(iHarvstType_pft(NZ,IDY,NY,NX).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ,IDY,NY,NX).EQ.iharvtyp_herbivo)THEN
                !animal or insect biomass
                NN=NN+1
                if(mod(nn,2)==0)then
                  IDYE=IDY
                  D580: DO IDYG=IDYS+1,IDYE-1
                    iHarvstType_pft(NZ,IDYG,NY,NX)    = ICUT
                    jHarvstType_pft(NZ,IDYG,NY,NX)    = JCUT
                    CanopyCutProxy_pft(NZ,IDYG,NY,NX) = HCUT !impact height by disturbance, mostly
                    THIN_pft(NZ,IDYG,NY,NX)           = PCUT
                    
                    FracBiomHarvsted(iHarvst_pft,iplthvst_leaf,NZ,IDYG,NY,NX)        = ECUT11
                    FracBiomHarvsted(iHarvst_pft,iplthvst_finenonleaf,NZ,IDYG,NY,NX) = ECUT12
                    FracBiomHarvsted(iHarvst_pft,iplthvst_stalk,NZ,IDYG,NY,NX)       = ECUT13
                    FracBiomHarvsted(iHarvst_pft,iplthvst_stdead,NZ,IDYG,NY,NX)      = ECUT14

                    FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ,IDYG,NY,NX)        = ECUT21
                    FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ,IDYG,NY,NX) = ECUT22
                    FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ,IDYG,NY,NX)       = ECUT23
                    FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ,IDYG,NY,NX)      = ECUT24
                  ENDDO D580
                endif
                IDYS=IDY
              ENDIF
            ENDDO
          endif
        ENDDO
      ENDDO
    ENDDO
  ENDDO  

  call PrintInfo('end '//subname)      
  end subroutine readplantmgmtinfo
!------------------------------------------------------------------------------------------

  subroutine ReadPlantManagementNC(yearc,yeari,NHW,NHE,NVN,NVS,pft_changed)
  !
  !
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: yeari   !forcing data year
  integer, intent(in) :: yearc   !current model year
  logical, optional, intent(out):: pft_changed

  character(len=*), parameter :: subname='ReadPlantManagementNC'
  type(file_desc_t) :: pftinfo_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar
  integer :: pft_dflag
  integer :: iyear     !data record year
  integer :: year      !year of the data
  integer :: ntopou,NY,NX,NZ
  logical :: pft_changed_loc

  call PrintInfo('beg '//subname)
  pft_changed_loc=.false.
  if(present(pft_changed))pft_changed=.false.
  if (len_trim(pft_mgmt_in)==0)return
  
  call ncd_pio_openfile(pftinfo_nfid, pft_mgmt_in, ncd_nowrite)

  call check_var(pftinfo_nfid, 'pft_dflag', vardesc, readvar)

  call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_dflag), &
    trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))

  if(pft_dflag==-1)then
    ! ReadPlantInfoNC has already disabled plants. No year, topographic-unit,
    ! planting, or management records are required for the no-plant flag.
    call ncd_pio_closefile(pftinfo_nfid)
    call PrintInfo('end '//subname)
    return
  endif

  !obtain the number of topographic units
  ntopou=get_dim_len(pftinfo_nfid, 'ntopou')    
  if(first_topou)ntopou=1
  
  if(pft_dflag==0)then
    if(lverb)write(iulog,*)'Constant pft data'
    if(IGO==0)then
      pft_changed_loc=.true.
      iyear=1    
      call readplantinginfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
      call InitPlantMgmnt(NHW,NHE,NVN,NVS)
    elseif(IGO>0)then
      iyear=1
      call readplantmgmtinfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
    endif
   
  elseif(pft_dflag==1)then
    pft_changed_loc=.true.
    if(lverb)write(iulog,*)'Transient pft data',yeari,yearc
    iyear=1
    DO while(.true.)
      call ncd_getvar(pftinfo_nfid,'year',iyear,year)
      if(year==yeari)exit   !when year is found matching the forcing data yeari.
      iyear=iyear+1
    ENDDO
    
    call readplantinginfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)    
    call readplantmgmtinfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
  else
    pft_changed_loc=.true.

    !transient pft data (flexible inputs)
    ! Intentionally use the simulation year for flexible planting/management,
    ! while ReadPlantInfoNC selects species using the forcing year. The two
    ! record years may differ when climate forcing is cycled.
    iyear=1
    DO while(.true.)
      call ncd_getvar(pftinfo_nfid,'year',iyear,year)
      if(year==yearc)exit  !when year is found matching the current model year.
      iyear=iyear+1
    ENDDO
    call readplantinginfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)
    call readplantmgmtinfo(pftinfo_nfid,ntopou,iyear,yearc,NHW,NHE,NVN,NVS)  
  endif
  call ncd_pio_closefile(pftinfo_nfid)

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      DO NZ=1,NP_col(NY,NX)
        iHarvestDay_pft(NZ,NY,NX)      = iDayPlantHarvest_pft(NZ,NY,NX)
        iHarvestYear_pft(NZ,NY,NX)     = iYearPlantHarvest_pft(NZ,NY,NX)
        iMaintPlantTrait_pft(NZ,NY,NX) = iFalse
      ENDDO
    ENDDO
  ENDDO
  if(present(pft_changed))pft_changed=pft_changed_loc
  call PrintInfo('end '//subname)
  end subroutine ReadPlantManagementNC
!------------------------------------------------------------------------------------------

  subroutine ReadPlantProperties(nu_plt,NZ,NY,NX,pft_changed)
  use GridConsts, only : MaxNumBranches
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  logical, intent(in) :: pft_changed
  character(len=*), parameter :: subname='ReadPlantProperties'
  integer :: N, NB
  real(r8) :: VRNXI,VRNLI

! begin_execution

!
! READ INPUTS FOR EACH PLANT SPECIES
!
  IF(DATAP(NZ,NY,NX).NE.'NO')THEN

    call SetPlantTraits(nu_plt,pft_changed,NZ,NY,NX,VRNLI,VRNXI)

    ! Drought leaf-out uses PSIMin4LeafOut(2:4) in PlantPhenolMod.
    IF ((iPlantPhenolType_pft(NZ,NY,NX).EQ.iphenotyp_drouhtdecidu .OR. &
         iPlantPhenolType_pft(NZ,NY,NX).EQ.4 .OR. &
         iPlantPhenolType_pft(NZ,NY,NX).EQ.5) .AND. &
        (iEmbryophyteType_pft(NZ,NY,NX).LT.2 .OR. &
         iEmbryophyteType_pft(NZ,NY,NX).GT.4)) THEN
      CALL ENDRUN('Unsupported botanical/phenology combination for '// &
        TRIM(DATAP(NZ,NY,NX))//': drought leaf-out requires IEBTYP=2:4; '// &
        TRIM(mod_filename), __LINE__)
    ENDIF

!
!   RE-CALCULATE PLANT INPUTS IN MODEL UNITS
!
!   LeafSWabsorptivity_pft,LeafPARabsorptivity_pft=leaf SW,PAR absorbtivity
!   2.5 converts g enzyme mass into g enzyme C
    VmaxSpecRubCarboxyRef_pft(NZ,NY,NX)  = 1.92_r8*VmaxSpecRubCarboxyRef_pft(NZ,NY,NX)
    VmaxRubOxyRef_pft(NZ,NY,NX)          = 1.92_r8*VmaxRubOxyRef_pft(NZ,NY,NX)
    VmaxPEPCarboxyRef_pft(NZ,NY,NX)      = 1.92_r8*VmaxPEPCarboxyRef_pft(NZ,NY,NX)
    SpecLeafChlAct_pft(NZ,NY,NX)         = 1.35_r8*SpecLeafChlAct_pft(NZ,NY,NX)
    
    SineBranchAngle_pft(NZ,NY,NX)        = SIN(BranchAngle_pft(NZ,NY,NX)*RadianPerDegree)
    SinePetolShethAngle_pft(NZ,NY,NX)    = SIN(PetolShethAngle_pft(NZ,NY,NX)*RadianPerDegree)
    MatureGroup_pft(NZ,NY,NX)            = GROUPX_pft(NZ,NY,NX)
    MatureGroup_pft(NZ,NY,NX)            = MatureGroup_pft(NZ,NY,NX)-ShootNodeNumAtPlanting_pft(NZ,NY,NX)
    
    IF(iPlantTurnoverPattern_pft(NZ,NY,NX).NE.0)THEN
      RefNodeInitRate_pft(NZ,NY,NX)        = RefNodeInitRate_pft(NZ,NY,NX)/MaxNodesPerBranch
      RateRefLeafAppearance_pft(NZ,NY,NX)  = RateRefLeafAppearance_pft(NZ,NY,NX)/MaxNodesPerBranch
      MatureGroup_pft(NZ,NY,NX)            = MatureGroup_pft(NZ,NY,NX)/MaxNodesPerBranch
      ShootNodeNumAtPlanting_pft(NZ,NY,NX) = ShootNodeNumAtPlanting_pft(NZ,NY,NX)/MaxNodesPerBranch
    ENDIF

    D5: DO NB=1,MaxNumBranches
      IF(iPlantPhenolType_pft(NZ,NY,NX).EQ.iphenotyp_evgreen .AND. iPlantPhenolPattern_pft(NZ,NY,NX).EQ.iplt_perennial)THEN
        !perennial evergreen
        HourReq4LeafOut_brch(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNLI+144.0_r8*PlantInitThermoAdaptZone_pft(NZ,NY,NX)*(NB-1))
        HourReq4LeafOff_brch(NB,NZ,NY,NX)=AMIN1(4380.0_r8,VRNXI+144.0_r8*PlantInitThermoAdaptZone_pft(NZ,NY,NX)*(NB-1))
      ELSE        
        HourReq4LeafOut_brch(NB,NZ,NY,NX)=VRNLI
        HourReq4LeafOff_brch(NB,NZ,NY,NX)=VRNXI
      ENDIF
    ENDDO D5
  ENDIF
! WRITE(*,1111)'CRITICAL DAYLENGTH',IGO,NZ,CriticPhotoPeriod_pft(NZ,NY,NX)
!1111    FORMAT(A20,2I8,E12.4)
  end subroutine ReadPlantProperties

!------------------------------------------------------------------------------------------
  subroutine ReadPlantTraitsNC(nu_plt,NZ,NY,NX,VRNLI,VRNXI)

  use ncdio_pio
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  real(r8), intent(out) :: VRNLI,VRNXI
  character(len=40) :: pft_lname
  character(len=64):: koppen_climl
  character(len=3) :: koppen_clims
  integer :: loc,N

  loc=get_pft_loc(KoppenClimZone_col(NY,NX),DATAP(NZ,NY,NX)(1:6),pft_lname,koppen_climl,koppen_clims)
  DATAPI(NZ,NY,NX)=loc
  call ncd_getvar(pft_nfid, 'ICTYP', loc, iPlantPhotosynsType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IGTYP', loc, iPlantRootProfile_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ISTYP', loc, iPlantPhenolPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IDTYP', loc, iPlantDevelopPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'INTYP', loc, iPlantNfixType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IWTYP', loc, iPlantPhenolType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IPTYP', loc, iPlantPhotoperiodType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IBTYP', loc, iPlantTurnoverPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'I2EXP', loc, iPlant2ndGrothPattern_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ISNTYP',loc, iPlantSnowIntercepType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'IRTYP', loc, iPlantGrainType_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'MY', loc, Myco_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ZTYPI', loc, PlantInitThermoAdaptZone_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'VCMX', loc, VmaxSpecRubCarboxyRef_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VOMX', loc, VmaxRubOxyRef_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VCMX4', loc, VmaxPEPCarboxyRef_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKCO2', loc, XKCO2_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKO2', loc, XKO2_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid,'ROOTMAGE',loc,RootMatureAge_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XKCO24', loc, Km4PEPCarboxy_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RUBP', loc, LeafRubisco2Protein_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PEPC', loc, LeafPEP2Protein_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ETMX', loc, SpecLeafChlAct_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CHL', loc, LeafProtein2Chl_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'fCHLMESO', loc, fMesophyllChlProtein_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'FCO2', loc, CanopyCi2CaRatio_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'ALBR', loc, RadSWLeafAlbedo_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ALBP', loc, RadPARLeafAlbedo_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'TAUR', loc, RadSWLeafTransmitance_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'TAUP', loc, RadPARLeafTransmitance_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'XRNI', loc, RefNodeInitRate_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XRLA', loc, RateRefLeafAppearance_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CTC', loc, TCChill4Seed_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'VRNLI', loc,VRNLI)
  call ncd_getvar(pft_nfid, 'VRNXI', loc,VRNXI)
  call ncd_getvar(pft_nfid, 'WDLF', loc,rLen2WidthLeaf_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PB', loc,NonstCMinConc2InitBranch_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'GROUPX', loc,GROUPX_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XTLI', loc,ShootNodeNumAtPlanting_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XDL', loc,CriticPhotoPeriod_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'XPPD', loc,PhotoPeriodSens_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'SLA1', loc,SLA1_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SSL1', loc,PetolShethLen2Mass_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SNL1', loc,NodeLenPergC_pft(NZ,NY,NX))


  call ncd_getvar(pft_nfid, 'CLASS', loc,LeafAngleClass_pft(1:NumLeafInclinationClasses,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CFI', loc,ClumpFactorInit_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ANGBR', loc,BranchAngle_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'ANGSH', loc,PetolShethAngle_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'STMX', loc,GrothStalkMaxSeedSites_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'SDMX', loc,MaxSeedNumPerSite_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GRMX', loc,SeedCMassMax_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GRDM', loc,SeedCMass_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GRW2L',loc,SeedWidth2LenRatio_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'GFILL', loc,GrainFillRate25C_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'WTSTDI', loc,StandingDeadInitC_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'RRAD1M', loc,Root1stMaxRadius_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RRAD2M', loc,Root2ndMaxRadius_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PORT', loc,RootPorosity_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'KLGMAX',loc,KLigMax_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'KLGMM',loc,KLigMM_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PR', loc,NonstCMinCon2InitRoot_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid,'RVSR',loc,RootVesselRadius_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSRR', loc,RootRadialResist_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSRA', loc,RootSingleVesselRstaxial_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'PTSHT', loc,ShootRootNonstElmConduts_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RTFQ', loc,FineRootBranchFreq_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'UPMXZH', loc,VmaxNH4Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPKMZH', loc,KmNH4Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPMNZH', loc,CMinNH4Root_pft(ipltroot,NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'UPMXZO', loc,VmaxNO3Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPKMZO', loc,KmNO3Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPMNZO', loc,CminNO3Root_pft(ipltroot,NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'UPMXPO', loc,VmaxPO4Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPKMPO', loc,KmPO4Root_pft(ipltroot,NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'UPMNPO', loc,CMinPO4Root_pft(ipltroot,NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'OSMO', loc,OrganOsmoPsi0pt_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RCS', loc,RCS_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'RSMX', loc,CuticleResist_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'DMLF', loc,LeafBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMSHE', loc,PetolShethBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMSTK', loc,StalkBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMRSV', loc,ReserveBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMHSK', loc,HuskBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMEAR', loc,EarBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMGR', loc,GrainBiomGrowthYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMRT', loc,RootBiomGrosYld_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'DMND', loc,NoduGrowthYield_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'CNWL', loc,rProteinC2LeafN_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPWL', loc,rProteinC2LeafP_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNWR', loc,rProteinC2RootN_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPWR', loc,rProteinC2RootP_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNLF', loc,rNCLeaf_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNSHE', loc,rNCSheath_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid,'CNRTLIG', loc,rNCLigRoot_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNSTK', loc,rNCStalk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNRSV', loc,rNCReserve_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNHSK', loc,rNCHusk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNEAR', loc,rNCEar_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNGR', loc,rNCGrain_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNRT', loc,rNCRoot_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CNND', loc,rNCNodule_pft(NZ,NY,NX))

  call ncd_getvar(pft_nfid, 'CPLF', loc,rPCLeaf_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPSHE', loc,rPCSheath_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPSTK', loc,rPCStalk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid,'CPRTLIG', loc,rPCLigRoot_pft(NZ,NY,NX))  
  call ncd_getvar(pft_nfid, 'CPRSV', loc,rPCReserve_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPHSK', loc,rPCHusk_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPEAR', loc,rPCEar_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPGR', loc,rPCGrain_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPRT', loc,rPCRootr_pft(NZ,NY,NX))
  call ncd_getvar(pft_nfid, 'CPND', loc,rPCNoduler_pft(NZ,NY,NX))

  if(disp_planttrait)then
    call pft_display(nu_plt,NZ,NY,NX,pft_lname,koppen_climl,koppen_clims)
    call photosyns_trait_disp(nu_plt,NZ,NY,NX)
    call plant_optic_trait_disp(nu_plt,NZ,NY,NX)
    call Phenology_trait_disp(nu_plt,NZ,NY,NX,VRNLI,VRNXI)
    call Shoot_morphology_trait_disp(nu_plt,NZ,NY,NX)
    call Root_trait_disp(nu_plt,NZ,NY,NX)
    call Root_nutrient_trait_disp(nu_plt,NZ,NY,NX)
    call plant_water_trait_disp(nu_plt,NZ,NY,NX)
    call plant_biomyield_trait_disp(nu_plt,NZ,NY,NX)
    call plant_biomstoich_trait_disp(nu_plt,NZ,NY,NX)
  endif

  end subroutine ReadPlantTraitsNC
!------------------------------------------------------------------------------------------
  subroutine SetPlantTraits(nu_plt,pft_changed,NZ,NY,NX,VRNLI,VRNXI)

  use ncdio_pio
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  logical, intent(in) :: pft_changed
  real(r8), intent(out) :: VRNLI,VRNXI

  character(len=*), parameter :: subname='SetPlantTraits'
  character(len=40) :: pft_lname
  character(len=64):: koppen_climl
  character(len=3) :: koppen_clims
  integer :: loc,N
  
  if(ats_cpl_mode)then
    loc=DATAPI(NZ,NY,NX)
  else
    loc=get_pft_loc(KoppenClimZone_col(NY,NX),DATAP(NZ,NY,NX)(1:6),pft_lname,koppen_climl,koppen_clims)
    !trait-table index of the PFT in this slot, kept in restart files so that a
    !fallow slot can reload the traits of the PFT whose residue it holds
    DATAPI(NZ,NY,NX)=loc
  endif
  StemSpecVolume_pft(NZ,NY,NX) = 1.e-6_r8/StemMassDensity_tab(loc)
  iEmbryophyteType_pft(NZ,NY,NX)     = iEmbryophyteType_pft_tab(loc)
  iPlantPhotosynsType_pft(NZ,NY,NX)  = iPlantPhotosynsType_pft_tab(loc)
  iPlantRootProfile_pft(NZ,NY,NX)     = iPlantRootProfile_tab(loc)
  iPlantPhenolPattern_pft(NZ,NY,NX)   = iPlantPhenolPattern_tab(loc)
  iPlantDevelopPattern_pft(NZ,NY,NX)  = iPlantDevelopPattern_tab(loc)
  iPlantNfixType_pft(NZ,NY,NX)        = iPlantNfixType_tab(loc)
  iPlantPhenolType_pft(NZ,NY,NX)      = iPlantPhenolType_tab(loc)
  iPlantPhotoperiodType_pft(NZ,NY,NX) = iPlantPhotoperiodType_tab(loc)
  iPlantTurnoverPattern_pft(NZ,NY,NX) = iPlantTurnoverPattern_tab(loc)
  iPlant2ndGrothPattern_pft(NZ,NY,NX) = iPlant2ndGrothPattern_tab(loc)
  iPlantGrainType_pft(NZ,NY,NX)       = iPlantGrainType_tab(loc)
  Myco_pft(NZ,NY,NX)                  = Myco_tab(loc)
  PlantInitThermoAdaptZone_pft(NZ,NY,NX)  = PlantInitThermoAdaptZone_tab(loc)
  xylemPhi_min_pft(NZ,NY,NX)          = xylemPhi_min_tab(loc)
  Radius95pctMature_pft(NZ,NY,NX)     = Radius95pctMature_tab(loc)
  xylemPhi_max_pft(NZ,NY,NX)          = xylemPhi_max_tab(loc)
  xylemPhi_mean_pft(NZ,NY,NX)         = xylemPhi_mean_tab(loc)

  VmaxSpecRubCarboxyRef_pft(NZ,NY,NX) = VmaxRubCarboxyRef_tab(loc)
  VmaxRubOxyRef_pft(NZ,NY,NX)         = VmaxRubOxyRef_tab(loc)
  VmaxPEPCarboxyRef_pft(NZ,NY,NX)     = VmaxPEPCarboxyRef_tab(loc)
  XKCO2_pft(NZ,NY,NX)                 = XKCO2_tab(loc)
  XKO2_pft(NZ,NY,NX)                  = XKO2_tab(loc)
  RootMatureAge_pft(NZ,NY,NX)         = RootMatureAge_tab(loc)
  Km4PEPCarboxy_pft(NZ,NY,NX)         = Km4PEPCarboxy_tab(loc)
  LeafRubisco2Protein_pft(NZ,NY,NX)   = LeafRubisco2Protein_tab(loc)
  LeafPEP2Protein_pft(NZ,NY,NX)       = LeafPEP2Protein_tab(loc)
  SpecLeafChlAct_pft(NZ,NY,NX)        = SpecChloryfilAct_tab(loc)
  iPlantSnowIntercepType_pft(NZ,NY,NX)= iPlantSnowIntercepType_tab(loc)
  LeafProtein2Chl_pft(NZ,NY,NX)       = LeafChl2Protein_tab(loc)
  fMesophyllChlProtein_pft(NZ,NY,NX)  = fChlMesophyll_tab(loc)
  CanopyCi2CaRatio_pft(NZ,NY,NX)      = CanopyCi2CaRatio_pft_tab(loc)

  RadSWLeafAlbedo_pft(NZ,NY,NX)        = RadSWLeafAlbedo_tab(loc)
  RadPARLeafAlbedo_pft(NZ,NY,NX)        = CanopyPARalbedo_tab(loc)
  RadSWLeafTransmitance_pft(NZ,NY,NX)  = RadSWLeafTransmis_tab(loc)
  RadPARLeafTransmitance_pft(NZ,NY,NX) = RadPARLeafTransmis_tab(loc)

  RefNodeInitRate_pft(NZ,NY,NX)          = RefNodeInitRate_tab(loc)
  RateRefLeafAppearance_pft(NZ,NY,NX)    = RefLeafAppearRate_tab(loc)
  TCChill4Seed_pft(NZ,NY,NX)             = TCChill4Seed_tab(loc)
  VRNLI                                  = VRNLI_tab(loc)
  VRNXI                                  = VRNXI_tab(loc)
  enh_cyto_pft(NZ,NY,NX)                 = enh_cyto_tab(Loc)
  StalkAxialResist_pft(NZ,NY,NX)         = StalkAxialResist_tab(loc)
  rLen2WidthLeaf_pft(NZ,NY,NX)           = rLen2WidthLeaf_tab(loc)
  NonstCMinConc2InitBranch_pft(NZ,NY,NX) = NonstCMinConc2InitBranch_tab(loc)

  GROUPX_pft(NZ,NY,NX)                 = GROUPX_tab(loc)
  ShootNodeNumAtPlanting_pft(NZ,NY,NX) = ShootNodeNumAtPlanting_tab(loc)
  CriticPhotoPeriod_pft(NZ,NY,NX)      = CriticPhotoPeriod_tab(loc)
  PhotoPeriodSens_pft(NZ,NY,NX)        = PhotoPeriodSens_tab(loc)

  SLA1_pft(NZ,NY,NX)         = SLA1_tab(loc)
  PetolShethLen2Mass_pft(NZ,NY,NX) = PetolShethLen2Mass_tab(loc)
  NodeLenPergC_pft(NZ,NY,NX) = NodeLenPergC_tab(loc)

  LeafAngleClass_pft(1:NumLeafInclinationClasses,NZ,NY,NX) = LeafAngleClass_tab(1:NumLeafInclinationClasses,loc)
  ClumpFactorInit_pft(NZ,NY,NX)                       = ClumpFactorInit_tab(loc)
  BranchAngle_pft(NZ,NY,NX)                           = BranchAngle_tab(loc)
  PetolShethAngle_pft(NZ,NY,NX)                       = PetolShethAngle_tab(loc)
  GrothStalkMaxSeedSites_pft(NZ,NY,NX)                = MaxPotentSeedNumber_tab(loc)
  MaxSeedNumPerSite_pft(NZ,NY,NX)                     = MaxSeedNumPerSite_tab(loc)
  SeedCMassMax_pft(NZ,NY,NX)                          = SeedCMassMax_tab(loc)
  SeedCMass_pft(NZ,NY,NX)                             = SeedCMass_tab(loc)
  SeedWidth2LenRatio_pft(NZ,NY,NX)                    = SeedWidth2LenRatio_tab(loc)
  GrainFillRate25C_pft(NZ,NY,NX)                      = GrainFillRate25C_tab(loc)
  StandingDeadInitC_pft(NZ,NY,NX)                     = StandingDeadInitC_tab(loc)

  rECLiveCRoot_pft(ielmn,NZ,NY,NX) = rNCLiveCRoot_tab(loc)
  rECLiveCRoot_pft(ielmp,NZ,NY,NX) = rPCLiveCRoot_tab(loc)
  rECDeadCRoot_pft(ielmn,NZ,NY,NX) = rNCDeadCRoot_tab(loc)
  rECDeadCRoot_pft(ielmp,NZ,NY,NX) = rPCDeadCRoot_tab(loc)

  Root1stMaxRadius_pft(1,NZ,NY,NX)       = Root1stMaxRadius_tab(loc)
  Root2ndMaxRadius_pft(1,NZ,NY,NX)       = Root2ndMaxRadius_tab(loc)
  RootPorosity_pft(1,NZ,NY,NX)           = RootPorosity_tab(loc)
  KLigMax_pft(NZ,NY,NX)                  = KLigMax_tab(loc)
  KLigMM_pft(NZ,NY,NX)                   = KLigMM_tab(loc)
  NonstCMinCon2InitRoot_pft(NZ,NY,NX)    = NonstCMinCon2InitRoot_tab(loc)
  RootRadialResist_pft(1,NZ,NY,NX)       = RootRadialResist_tab(loc)
  AlphaVesselAxialResist_pft(NZ,NY,NX)    = AlphaVesselAxialResist_tab(loc)
  RootVesselRadius_pft(NZ,NY,NX)         = RootVesselRadius_tab(loc)
  
  ShootRootNonstElmConduts_pft(NZ,NY,NX) = ShutRutNonstElmntConducts_tab(loc)
  FineRootBranchFreq_pft(NZ,NY,NX)       = FineRootBranchFreq_tab(loc)
  MediumRootBranchFreq_pft(NZ,NY,NX)     = MediumRootBranchFreq_tab(loc)
  VmaxNH4Root_pft(ipltroot,NZ,NY,NX)     = VmaxNH4Root_tab(loc)
  KmNH4Root_pft(ipltroot,NZ,NY,NX)       = KmNH4Root_tab(loc)
  CMinNH4Root_pft(ipltroot,NZ,NY,NX)     = CMinNH4Root_tab(loc)

  MorphogenBase_pft(NZ,NY,NX) = MorphogenBase_tab(loc)*0.01_r8
  VmaxNO3Root_pft(ipltroot,NZ,NY,NX) = VmaxNO3Root_tab(loc)
  KmNO3Root_pft(ipltroot,NZ,NY,NX)   = KmNO3Root_tab(loc)
  CminNO3Root_pft(ipltroot,NZ,NY,NX) = CminNO3Root_tab(loc)

  VmaxPO4Root_pft(ipltroot,NZ,NY,NX) = VmaxPO4Root_tab(Loc)
  KmPO4Root_pft(ipltroot,NZ,NY,NX)   = KmPO4Root_tab(loc)
  CMinPO4Root_pft(ipltroot,NZ,NY,NX) = CMinPO4Root_tab(loc)

  OrganOsmoPsi0pt_pft(NZ,NY,NX) = CanOsmoPsi0pt_tab(loc)
  RCS_pft(NZ,NY,NX)           = RCS_tab(loc)
  CuticleResist_pft(NZ,NY,NX) = CuticleResist_tab(loc)

  if(ats_cpl_mode)then
    irootType_pft(NZ,NY,NX)       = irootType_tab(loc)
    CanopyHeightLive_pft(NZ,NY,NX) = CanopyHeightLive_tab(loc)
  endif

  LeafBiomGrowthYld_pft(NZ,NY,NX)    = LeafBiomGrowthYld_tab(loc)
  PetolShethBiomGrowthYld_pft(NZ,NY,NX) = PetolShethBiomGrowthYld_tab(loc)
  StalkBiomGrowthYld_pft(NZ,NY,NX)   = StalkBiomGrowthYld_tab(loc)
  ReserveBiomGrowthYld_pft(NZ,NY,NX) = ReserveBiomGrowthYld_tab(loc)
  HuskBiomGrowthYld_pft(NZ,NY,NX)    = HuskBiomGrowthYld_tab(loc)
  EarBiomGrowthYld_pft(NZ,NY,NX)     = EarBiomGrowthYld_tab(loc)
  GrainBiomGrowthYld_pft(NZ,NY,NX)   = GrainBiomGrowthYld_tab(loc)
  RootBiomGrosYld_pft(NZ,NY,NX)      = RootBiomGrosYld_tab(loc)
  NoduGrowthYield_pft(NZ,NY,NX)      = NoduGrowthYield_tab(loc)  
  rProteinC2LeafN_pft(NZ,NY,NX) = rProteinC2LeafN_tab(loc)
  rProteinC2RootN_pft(NZ,NY,NX) = rProteinC2RootN_tab(loc)
  rNCLeaf_pft(NZ,NY,NX)    = rNCLeaf_tab(loc)
  rNCSheath_pft(NZ,NY,NX)  = rNCSheath_tab(loc)
  rNCStalk_pft(NZ,NY,NX)   = rNCStalk_tab(loc)
  rNCLigRoot_pft(NZ,NY,NX) = rNCLigRoot_tab(loc)
  rNCReserve_pft(NZ,NY,NX) = rNCReserve_tab(loc)
  rNCHusk_pft(NZ,NY,NX)    = rNCHusk_tab(loc)
  rNCEar_pft(NZ,NY,NX)     = rNCEar_tab(loc)
  rNCGrain_pft(NZ,NY,NX)   = rNCGrain_tab(loc)
  rNCRoot_pft(NZ,NY,NX)    = rNCRoot_tab(loc)
  rNCNodule_pft(NZ,NY,NX)  = rNCNodule_tab(loc)

  rPCLigRoot_pft(NZ,NY,NX) = rPCLigRoot_tab(loc)  
  rProteinC2LeafP_pft(NZ,NY,NX) = rProteinC2LeafP_tab(loc)
  rProteinC2RootP_pft(NZ,NY,NX) = rProteinC2RootP_tab(loc)
  rPCLeaf_pft(NZ,NY,NX)    = rPCLeaf_tab(loc)
  rPCSheath_pft(NZ,NY,NX)  = rPCSheath_tab(loc)
  rPCStalk_pft(NZ,NY,NX)   = rPCStalk_tab(loc)
  rPCReserve_pft(NZ,NY,NX) = rPCReserve_tab(loc)
  rPCHusk_pft(NZ,NY,NX)    = rPCHusk_tab(loc)
  rPCEar_pft(NZ,NY,NX)     = rPCEar_tab(loc)
  rPCGrain_pft(NZ,NY,NX)   = rPCGrain_tab(loc)
  rPCRootr_pft(NZ,NY,NX)   = rPCRootr_tab(loc)
  rPCNoduler_pft(NZ,NY,NX) = rPCNoduler_tab(loc)

  IF(CriticPhotoPeriod_pft(NZ,NY,NX).LT.0.0_r8)THEN
    CriticPhotoPeriod_pft(NZ,NY,NX)=DayLenthMax_col(NY,NX)
  ENDIF
  if(disp_planttrait.and.pft_changed)then
    call pft_display(nu_plt,NZ,NY,NX,pft_lname,koppen_climl,koppen_clims)
    call photosyns_trait_disp(nu_plt,NZ,NY,NX)
    call plant_optic_trait_disp(nu_plt,NZ,NY,NX)
    call Phenology_trait_disp(nu_plt,NZ,NY,NX,VRNLI,VRNXI)
    call Shoot_morphology_trait_disp(nu_plt,NZ,NY,NX)
    call Root_trait_disp(nu_plt,NZ,NY,NX)
    call Root_nutrient_trait_disp(nu_plt,NZ,NY,NX)
    call plant_water_trait_disp(nu_plt,NZ,NY,NX)
    call plant_biomyield_trait_disp(nu_plt,NZ,NY,NX)
    call plant_biomstoich_trait_disp(nu_plt,NZ,NY,NX)
  endif
  
  end subroutine SetPlantTraits

!------------------------------------------------------------------------------------------
  subroutine pft_display(nu_plt,NZ,NY,NX,pft_lname,koppen_climl,koppen_clims)
  use abortutils , only : endrun  
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  character(len=*),intent(in) :: pft_lname  
  character(len=*),intent(in) :: koppen_climl
  character(len=*),intent(in) :: koppen_clims
  integer :: j,id
  character(len=60) :: strval

!   iPlantPhotosynsType_pft=photosynthesis type:3=C3,4=C4
!   iPlantRootProfile_pft=root profile:0=shallow (eg bryophytes),1=intermediate(eg herbs),2=deep (eg trees)
!   iPlantPhenolPattern_pft=growth habit:0=annual,1=perennial
!   iPlantDevelopPattern_pft=growth habit:0=determinate,1=indeterminate
!   iPlantNfixType_pft=N2 fixation:1,2,3=rapid to slow root symbiosis (e.g.legumes),
!   4,5,6=rapid to slow canopy symbiosis (e.g. cyanobacteria)
!   iPlantPhenolType_pft=phenology type:0=evergreen,1=cold deciduous,2=drought deciduous,3=1+2
!   iPlantPhotoperiodType_pft=photoperiod type:0=day neutral,1=short day,2=long day
!   iPlantTurnoverPattern_pft=turnover:if iPlantRootProfile_pft=0 or 1:all above-ground:0,1=rapid(deciduous),2=very slow(evergreen),3=slow(semi-deciduous)
!                   :if iPlantRootProfile_pft=2:trees:1=rapid(deciduous),2=very slow(coniferous),3=slow(semi-deciduous)
!   iPlantGrainType_pft=storage organ:0=above ground,1=below ground
!   MY=mycorrhizal:1=no,2=yes
!   PlantInitThermoAdaptZone_pft=thermal adaptation zone:1=arctic,boreal,2=cool temperate,
!   3=warm temperate,4=subtropical,5=tropical
  write(nu_plt,*)('=',j=1,110)
  write(nu_plt,*)'PLANT traits for FUNCTIONAL TYPE (NZ,NY,NX)=',NZ,NY,NX,trim(DATAP(NZ,NY,NX)(1:6))
  id=0
  call writefixsl(nu_plt,'Plant name ',pft_lname,60)
  strval=koppen_clims//','//koppen_climl
  call writefixsl(nu_plt,'Koppen climate info',strval,60)
  write(nu_plt,*)('-',j=1,110)
  write(nu_plt,*)'PLANT CLASS INFORMATION'
  select CASE (iPlantPhotosynsType_pft(NZ,NY,NX))
  case (3)
    strval='C3'
  case (4)
    strval='C4'
  case default
    strval='Not defined'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'ICTYP: Photosynthesis pathway',strval,60,id)

  select case (iEmbryophyteType_pft(NZ,NY,NX))
  CASE (iembryotyp_Bryophytes)
    strval='Bryophytes'
  case (iembryotyp_Pteridophytes)
    strval='Pteridophytes'
  case (iembroytyp_Gymnosperms)
    strval='Gymnosperms'
  case (iembroytyp_Monocots)
    strval='Monocots'
  case (iembroytyp_Eudicots)
    strval='Eudicots'
  case default
     strval='Not defined'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'IEBTYP: Embrophyte type',strval,60,id)

  select case(iPlantRootProfile_pft(NZ,NY,NX))
  case (0)
    strval='Shallow root profile, like bryophytes (0)'
  case (1)
    strval='Intermediate deep fibrous root profile, like herbs (1)'
  case (2)
    strval='Deep penetrative root profile, like trees/drought tolerant grasses (2)'
  case (3)
    strval='Very deep tap root profile (3)'  
  case default
    strval='Not defined'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'IGTYP: Root profile pattern',strval,60,id)

  select case (iPlantPhenolPattern_pft(NZ,NY,NX))
  case (0)
    strval='Annual'
  case (1)
    strval='Perennial'
  case default
    strval='Not defined'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'ISTYP: Life cycle',strval,60,id)

  select case (iPlantDevelopPattern_pft(NZ,NY,NX))
  case (0)
    strval='Determinate'
  case (1)
    strval='Indeterminate'
  case default
    strval='Not defined'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'IDTYP: Growth pattern',strval,60,id)

  select case (iPlantNfixType_pft(NZ,NY,NX))
! 1,2, 3, e.g. legumes
  case (in2fixtyp_root_fast)
    strval='Rapid root N-fixation'
  case (in2fixtyp_root_medium)
    strval='Intermediate root N-fixation'
  case (in2fixtyp_root_slow)
    strval='Slow root N-fixation'
!4,5,6, e.g. cyanobacteria
  case (in2fixtyp_canopy_fast)
    strval='Rapid canopy N-fixation'
  case (in2fixtyp_canopy_medium)
    strval='Intermediate canopy N-fixation'
  case (in2fixtyp_canopy_slow)
    strval='Slow canopy N-fixation'
  case default
    strval='Not defined'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'INTYP: N-fixation symbiosis',strval,60,id)

  select case(iPlantPhenolType_pft(NZ,NY,NX))
  case (iphenotyp_evgreen)
    write(strval,'(A,I3)')'Evergreen ',iPlantPhenolType_pft(NZ,NY,NX)
  case (iphenotyp_coldecid)
    write(strval,'(A,I3)')'Cold deciduous ',iPlantPhenolType_pft(NZ,NY,NX)
  case (iphenotyp_drouhtdecidu)
    write(strval,'(A,I3)')'Drought deciduous ',iPlantPhenolType_pft(NZ,NY,NX)
  case (iphenotyp_coldroutdecid)
    write(strval,'(A,I3)')'Cold&drought-tolerant deciduous ',iPlantPhenolType_pft(NZ,NY,NX)
  case default
    strval='Not defined'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'IWTYP: Phenology type',strval,60,id)

  select case(iPlantPhotoperiodType_pft(NZ,NY,NX))
  case (iphotop_neutral)
    strval='Neutral day'
  case (iphotop_short)
    strval='Short day'
  case (iphotop_long)
    strval='Long day'
  case default
    strval='Not defined'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'IPTYP: Photoperiod',strval,60,id)

  if(is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX)))then
    select case(iPlantTurnoverPattern_pft(NZ,NY,NX))
    case (0, 1)
      write(strval,'(A,I2)')'Rapid, like deciduous trees ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case (2)
      write(strval,'(A,I2)')'Very slow, like evergreen trees ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case (3)
      write(strval,'(A,I2)')'Slow, like evergreen broadleaf trees',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case (4)
      write(strval,'(A,I2)')'Like semi-deciduous trees ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case (5)
      write(strval,'(A,I2)')'Like semi-evergreen trees ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case default
      write(strval,'(A,I2)')'Undefined tree-like pattern ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    end select
  else
    select case(iPlantTurnoverPattern_pft(NZ,NY,NX))
    case (0, 1)
      write(strval,'(A,I2)')'Rapid all aboveground biome, herbaceous ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    case (3)  
      write(strval,'(A,I2)')'Very slow, like evergreen grass ',iPlantTurnoverPattern_pft(NZ,NY,NX)      
    case default
      write(strval,'(A,I2)')'Undefined herbaceous pattern ',iPlantTurnoverPattern_pft(NZ,NY,NX)
    end select
  endif
  id=addone(id)
  call writefixsl(nu_plt,'IBTYP: Biome turnover pattern',strval,60,id)

  id=addone(id)
  if(iPlant2ndGrothPattern_pft(NZ,NY,NX).eq.0)THEN
    call writefixsl(nu_plt,'I2EXP: Secondary growth expression','No (0)',60,id)
  elseif(iPlant2ndGrothPattern_pft(NZ,NY,NX).eq.1)THEN
    call writefixsl(nu_plt,'I2EXP: Secondary growth expression','Yes (1)',60,id)
  ENDIF

  id=addone(id)
  select case(iPlantSnowIntercepType_pft(NZ,NY,NX))
  case (isnowcept_bryophyte)
    write(strval,'(A,I2)')'bryophyte ',iPlantSnowIntercepType_pft(NZ,NY,NX)
  case (isnowcept_grasses)
    write(strval,'(A,I2)')'grasse ',iPlantSnowIntercepType_pft(NZ,NY,NX)
  case (isnowcept_shrub)
    write(strval,'(A,I2)')'shrub ',iPlantSnowIntercepType_pft(NZ,NY,NX)  
  case (isnowcept_decidtree)
    write(strval,'(A,I2)')'deciduous tree ',iPlantSnowIntercepType_pft(NZ,NY,NX)  
  case (isnowcept_coniftree)
    write(strval,'(A,I2)')'confierous tree ',iPlantSnowIntercepType_pft(NZ,NY,NX)
  case default
    write(strval,'(A,I2)')'Undefined ',iPlantSnowIntercepType_pft(NZ,NY,NX)
  end select
  call writefixsl(nu_plt,'ISNTYP: snow interception pattern',strval,60,id)

  select case(iPlantGrainType_pft(NZ,NY,NX))
  case (igraintyp_abvgrnd)
    strval='Aboveground 0'
  case (igraintyp_blwgrnd)
    strval='Belowground 1'
  case default
    strval='Not defined'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'IRTYP: Storage organ',strval,60,id)

  select case(Myco_pft(NZ,NY,NX))
  case (1)
    strval='No'
  case (2)
    strval='YES'
  case default
    strval='Wrong option'
  end select
  id=addone(id)
  call writefixsl(nu_plt,'MY   : Mycorrhizal association MY',strval,60,id)
  
  if(PlantInitThermoAdaptZone_pft(NZ,NY,NX)<0 .or. PlantInitThermoAdaptZone_pft(NZ,NY,NX)>ithermozone_tropical)then
    write(strval,'(A,F6.2)')'Not defined ',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif(PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_arcboreal)then
    write(strval,'(A,F6.2)')'Arctic, boreal ',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif (PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_cooltempr)then
    write(strval,'(A,F6.2)')'Cool temperate ',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif (PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_warmtempr)then
    write(strval,'(A,F6.2)')'Warm temperate ',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif (PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_subtropic)then
    write(strval,'(A,F6.2)')'Subtropical ',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  elseif (PlantInitThermoAdaptZone_pft(NZ,NY,NX)<=ithermozone_tropical)then
    write(strval,'(A,F6.2)')'Tropical ',PlantInitThermoAdaptZone_pft(NZ,NY,NX)
  endif
  id=addone(id)
  call writefixsl(nu_plt,'ZTYPI: Growing season thermal adaptation zone',strval,60,id)
  end subroutine pft_display

!------------------------------------------------------------------------------------------
  subroutine photosyns_trait_disp(nu_plt,NZ,NY,NX)
  !
  !Description
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j, id
  character(len=81) :: line
  write(nu_plt,*)('-',j=1,110)
  write(nu_plt,*)'PHOTOSYNTHETIC PROPERTIES'
  id=0;
  if(iPlantPhotosynsType_pft(NZ,NY,NX).eq.ic3_photo)then
    id=addone(id)
    call writefixl(nu_plt,id,'VCMX','Saturated specific carboxylation rate by Rubisco at 25oC  [umol CO2 (gC rubisco)-1 s-1]',VmaxSpecRubCarboxyRef_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'VOMX','Saturated specific oxygenation rate  at 25oC by Rubisco  [umol O2 (gC rubisco)-1 s-1]',VmaxRubOxyRef_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'XKCO2','Reference Km for rubisco carboxylation of CO2  at 25oC  [uM] ',XKCO2_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'XKO2','Reference Km for rubisco oxygenation  at 25oC  [uM]',XKO2_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'RUBP','Fraction of total leaf protein as Rubisco enzyme [gC rubisco/(gC protein)]',LeafRubisco2Protein_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'ETMX','Specific chlorophyll activity  [umol e- (gC chl)-1 s-1]',SpecLeafChlAct_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'CHL','Fraction of total leaf protein as chlorophyll-bound for light-harvesting [-]',LeafProtein2Chl_pft(NZ,NY,NX),110)
  elseif(iPlantPhotosynsType_pft(NZ,NY,NX).eq.ic4_photo)then
    id=addone(id)
    call writefixl(nu_plt,id,'VCMX','Saturated specific carboxylation rate by Rubisco  at 25oC  [umol CO2 (gC rubisco)-1 s-1]',VmaxSpecRubCarboxyRef_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'VOMX','Saturated specific oxygenation rate by Rubisco  at 25oC  [umol O2 (gC rubisco)-1 s-1]',VmaxRubOxyRef_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'XKCO2','Reference Km for rubisco carboxylation of CO2  at 25oC  [uM]',XKCO2_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'XKO2','Reference Km for rubisco oxygenation at 25oC  [uM]',XKO2_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'RUBP','Fraction of total leaf protein as Rubisco enzyme  [gC rubisco/(gC protein)]',LeafRubisco2Protein_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'ETMX','Specific chlorophyll activity at 25oC [umol e- (gC chl)-1 s-1]',SpecLeafChlAct_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'VCMX4','Saturated specific PEP carboxylation rate  at 25oC  [umol CO2 (gC PEPC)-1 s-1]',VmaxPEPCarboxyRef_pft(NZ,NY,NX),110)  
    id=addone(id)
    call writefixl(nu_plt,id,'XKCO24','Reference KM for PEP carboxylation  at 25oC  [uM]',Km4PEPCarboxy_pft(NZ,NY,NX),110)    
    id=addone(id)
    call writefixl(nu_plt,id,'PEPC','Fraction of total leaf protein as PEPC enzyme [gC PEP/(gC protein)]',LeafPEP2Protein_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'CHL','Fraction of total leaf protein that as chlorophyll-bound for light-harvesting protein [-]',LeafProtein2Chl_pft(NZ,NY,NX),110)
    id=addone(id)
    call writefixl(nu_plt,id,'fCHLMESO','Fraction of chlorophyll-bound light-harvesting proteins in mesophyll cells [-]',fMesophyllChlProtein_pft(NZ,NY,NX),110)
  endif  
  id=addone(id)
  call writefixl(nu_plt,id,'FCO2','Intercellular-to-atmospheric CO2 concentration ratio  [-]',CanopyCi2CaRatio_pft(NZ,NY,NX),110)
  end subroutine photosyns_trait_disp

!------------------------------------------------------------------------------------------

  subroutine Phenology_trait_disp(nu_plt,NZ,NY,NX,VRNLI,VRNXI)

  implicit none
  integer, intent(in) :: NZ,NY,NX
  real(r8), intent(in) :: VRNLI,VRNXI
  integer, intent(in) :: nu_plt
  integer :: J,id

  write(nu_plt,*)('-',j=1,110)
  id=0
  write(nu_plt,*)'PHENOLOGICAL PROPERTIES'
  id=addone(id)
  call writefixl(nu_plt,id,'XRNI','Rate of node initiation at 25oC [h-1]',RefNodeInitRate_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'XRLA','Rate of leaf (dis)appearance at 25oC [h-1]',RateRefLeafAppearance_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CTC','Chilling temperature for CO2 fixation, seed loss [oC]',TCChill4Seed_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'VRNLI','Hours of favorable conditions (temp or water) required for spring leafout [h]',VRNLI,100)
  id=addone(id)
  call writefixl(nu_plt,id,'VRNXI','Hours of senescence conditions (temp or water) required for autumn leafoff [h]',VRNXI,100)
  id=addone(id)
  call writefixl(nu_plt,id,'PB','Nonstructural C concentration needed for shoot branching [gC nonst/gC structl]',NonstCMinConc2InitBranch_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'GROUPX','Maturity group, node number (on top of XTLI) required for floral initiation [-]',GROUPX_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'XTLI','Embryonic node number at planting [-]',ShootNodeNumAtPlanting_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'XDL','Critical photoperiod for leaf and flora development (<= maximum daylength;<0=maximum) [h]',CriticPhotoPeriod_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'XPPD','Critical difference between critical photoperiod and current day length [h]',PhotoPeriodSens_pft(NZ,NY,NX),100)

  end subroutine Phenology_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Shoot_morphology_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: N,j,id

  write(nu_plt,*)('-',j=1,110)  
  write(nu_plt,*)'MORPHOLOGICAL PROPERTIES'
  id=0;
  id=addone(id)
  call writefixl(nu_plt,id,'SLA1','Specific leaf area vs mass [m2 gC-1]',SLA1_pft(NZ,NY,NX),100)
  id=addone(id)  
  call writefixl(nu_plt,id,'SSL1','Specific Petiole/Sheath length vs mass [m gC-1]',PetolShethLen2Mass_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'SNL1','Specific stem internode length vs mass [m gC-1]',NodeLenPergC_pft(NZ,NY,NX),100)
  id=addone(id)
  call writeafixl(nu_plt,id,'CLASS','Fraction of leaf area in 0-22.5,22.5-45,45-67.5,67.5-90 degrees inclination classes (0 being flat) [-]',&
    LeafAngleClass_pft(1:NumLeafInclinationClasses,NZ,NY,NX),100)  
  id=addone(id)
  call writefixl(nu_plt,id,'STDENS','Stem dry mass density (including voids) [gC cm-3]',1.e-6_r8/StemSpecVolume_pft(NZ,NY,NX),100)    
  id=addone(id)
  call writefixl(nu_plt,id,'WDLF','Leaf length:width ratio [-]',rLen2WidthLeaf_pft(NZ,NY,NX),100)  
  id=addone(id)
  call writefixl(nu_plt,id,'CFI','Initial clumping factor seedling [-]',ClumpFactorInit_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'ANGBR','Stem angle from horizontal [degrees]',BranchAngle_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'ANGSH','PetolSheth angle from horizontal [degrees]',PetolShethAngle_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'STMX','Maximum potential seed sites for pre-anthesis stalk [sites (gC stalk)-1]',GrothStalkMaxSeedSites_pft(NZ,NY,NX),100)
  id=addone(id)  
  call writefixl(nu_plt,id,'SDMX','Maximum seed number per grain site [seeds site-1]',MaxSeedNumPerSite_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'GRMX','Maximum seed size per seed [gC seed-1]',SeedCMassMax_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'GRW2L','Seed width to length ratio at planting (assuming prolate spheroid)',SeedWidth2LenRatio_pft(NZ,NY,NX),100)   
  id=addone(id)
  call writefixl(nu_plt,id,'GRDM','Seed carbon mass at planting [gC seed-1]',SeedCMass_pft(NZ,NY,NX),100)    !could be greater than SeedCMassMax, accouting for seedling
  id=addone(id)
  call writefixl(nu_plt,id,'GFILL','Grain filling rate at 25 oC [gC seed-1 h-1]',GrainFillRate25C_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'WTSTDI','Mass of dead standing biomass at planting [gC m-2]',StandingDeadInitC_pft(NZ,NY,NX),100)
  if(is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX)))then
    id=addone(id)
    call writefixl(nu_plt,id,'CYTOS','Cytokinin sensitivity of corase root thickening [0-1]',enh_cyto_pft(NZ,NY,NX),100)    
  endif
  id=addone(id)
  call writefixl(nu_plt,id,'RSSA','Stalk axial resistance per m length for water uptake [MPa h m-4]',StalkAxialResist_pft(NZ,NY,NX),100)
  end subroutine Shoot_morphology_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j,id

  write(nu_plt,*)('-',j=1,110)
  write(nu_plt,*)'ROOT CHARACTERISTICS'
  id=0;id=addone(id)
  call writefixl(nu_plt,id,'RRAD1M','Maximum radius of young primary roots [m]',Root1stMaxRadius_pft(1,NZ,NY,NX),105)
  id=addone(id)
  call writefixl(nu_plt,id,'RRAD2M','Maximum radius of fine roots [m]',Root2ndMaxRadius_pft(1,NZ,NY,NX),105)
  id=addone(id)
  call writefixl(nu_plt,id,'PORT','Primary/fine root air porosity [m3 m-3]',RootPorosity_pft(1,NZ,NY,NX),105)

  if(is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX)))then
!    id=addone(id)
!    call writefixl(nu_plt,id,'KLGMAX','Maximum lignification rate of coarse roots [h-1]',KLigMax_pft(NZ,NY,NX),105)
!    id=addone(id)
!    call writefixl(nu_plt,id,'KLGMM','Half saturation parameter for coarse roots lignification [h-1]',KLigMM_pft(NZ,NY,NX),105)
    id=addone(id)
    call writefixl(nu_plt,id,'ROOTMAGE','Root age to trigger secondary growth [h]', RootMatureAge_pft(NZ,NY,NX),105)
    id=addone(id)
    call writefixl(nu_plt,id,'PhiMIN','The lumen fraction found in the youngest coarse root xylem of woody vascular [m2/m2]',xylemPhi_min_pft(NZ,NY,NX),105)
    id=addone(id)
    call writefixl(nu_plt,id,'PhiMAX','Asymptotic limit lumen fraction of the coarse root xyxlem area of woody vascular [m2/m2]', xylemPhi_max_pft(NZ,NY,NX),105)
    id=addone(id)
    call writefixl(nu_plt,id,'R95MAT','Critical radius where the woody root is considered 95% mature [m]', Radius95pctMature_pft(NZ,NY,NX),105)
    id=addone(id)    
    call writefixl(nu_plt,id,'MRTFQ','Medium root branching freq. on structrual roots [# m-1]', MediumRootBranchFreq_pft(NZ,NY,NX),105)    
  else
    id=addone(id)
    call writefixl(nu_plt,id,'PhiMean','The mean lumen area fraction found in the seminal roots of non-tree roots, [m2/m2]', xylemPhi_mean_pft(NZ,NY,NX),105)      
  endif
  id=addone(id)  
  call writefixl(nu_plt,id,'PR','Nonstructural C concentration needed for root'// &
    ' initiation (gC nonst/gC structl)',NonstCMinCon2InitRoot_pft(NZ,NY,NX),105)
  id=addone(id)  
  call writefixl(nu_plt,id,'RSRR','Radial resistcance per m2 root surface area for water uptake [MPa h m-1]',RootRadialResist_pft(1,NZ,NY,NX),105)
  id=addone(id)
  call writefixl(nu_plt,id,'ARSRA','Axial resistance correcting factor from lumen to actual vessel [-]',AlphaVesselAxialResist_pft(NZ,NY,NX),105)
  id=addone(id)
  call writefixl(nu_plt,id,'RVSR','Root vessel mean radius [m]',RootVesselRadius_pft(NZ,NY,NX),105)
  id=addone(id)
  call writefixl(nu_plt,id,'PTSHT','Rate scalar (<1) for equilibrating shoot-root '// &
    'nonstructural elemental concentrations [h-1]',ShootRootNonstElmConduts_pft(NZ,NY,NX),105)
  !as a rule of thumb, RTFQ often takes the value of .  
  id=addone(id)
  call writefixl(nu_plt,id,'RTFQ','Square root of (fine root branching frequency on 1st roots) X (root hair frequency on fine roots) [# m-1]',FineRootBranchFreq_pft(NZ,NY,NX),105)
  end subroutine Root_trait_disp

!------------------------------------------------------------------------------------------
  subroutine Root_nutrient_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j,id

  write(nu_plt,*)('-',j=1,110)
  write(nu_plt,*)'ROOT UPTAKE PARAMETERS'
  id=0;id=addone(id)
  call writefixl(nu_plt,id,'UPMXZH','Maximum rate for root uptake of NH4 [gN m-2 absorption area h-1]',VMaxNH4Root_pft(ipltroot,NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'UPKMZH','Half saturation Km for root uptake of NH4 [uM]',KmNH4Root_pft(ipltroot,NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'UPMNZH','Minimum concentration required for root uptake of NH4 [uM]',CMinNH4Root_pft(ipltroot,NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'UPMXZO','Maximum rate for root uptake of NO3 [gN m-2 absorption area h-1]',VMaxNO3Root_pft(ipltroot,NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'UPKMZO','Half saturation Km for root uptake of NO3 [uM]',KmNO3Root_pft(ipltroot,NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'UPMNZO','Minimum concentration for root uptake of NO3 [uM]',CMinNO3Root_pft(ipltroot,NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'UPMXPO','Maximum rate for root uptake of H2PO4 or H1PO4 [gP m-2 absorption area h-1]',VMaxPO4Root_pft(ipltroot,NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'UPKMPO','Half saturation Km for root uptake of H2PO4 or H1PO4 [uM]',KmPO4Root_pft(ipltroot,NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'UPMNPO','Minimum concentration required for root uptake of H2PO4 or H1PO4 [uM]',CMinPO4Root_pft(ipltroot,NZ,NY,NX),100)
  end subroutine Root_nutrient_trait_disp

!------------------------------------------------------------------------------------------
  subroutine plant_water_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j,id

  write(nu_plt,*)('-',j=1,110)    
  write(nu_plt,*)'WATER RELATIONS'
  id=0;id=addone(id)
  call writefixl(nu_plt,id,'OSMO','Organ osmotic potential at full turgor i.e. zero water potential [MPa]',OrganOsmoPsi0pt_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'RCS','e-folding turgor pressure for stomatal resistance [MPa]',RCS_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'RSMX','Cuticular resistance for H2O [s m-1]',CuticleResist_pft(NZ,NY,NX),100)
  end subroutine plant_water_trait_disp

!------------------------------------------------------------------------------------------
  subroutine plant_biomyield_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j,id

  write(nu_plt,*)('-',j=1,110)
  write(nu_plt,*)'ORGAN GROWTH YIELDS'
  id=0;id=addone(id)
  call writefixl(nu_plt,id,'DMLF','Leaf dry matter C production yield [gC leaf g-1 nonstrucal C]',LeafBiomGrowthYld_pft(NZ,NY,NX),101)
  id=addone(id)
  call writefixl(nu_plt,id,'DMSHE','PetolSheth dry matter C production yield [gC PetolSheth g-1 nonstrucal C]',PetolShethBiomGrowthYld_pft(NZ,NY,NX),101)
  id=addone(id)
  call writefixl(nu_plt,id,'DMSTK','Stalk dry matter C production yield [gC stalk g-1 nonstrucal C]',StalkBiomGrowthYld_pft(NZ,NY,NX),101)
  id=addone(id)
  call writefixl(nu_plt,id,'DMRSV','Stalk reserve C production yield [g reserve g-1 nonstrucal C]',ReserveBiomGrowthYld_pft(NZ,NY,NX),101)
  id=addone(id)
  call writefixl(nu_plt,id,'DMHSK','Husk dry matter C production yield [gC husk g-1 nonstrucal C]',HuskBiomGrowthYld_pft(NZ,NY,NX),101)
  id=addone(id)
  call writefixl(nu_plt,id,'DMEAR','Ear dry matter C production yiled [gC ear g-1 nonstrucal C]',EarBiomGrowthYld_pft(NZ,NY,NX),101)
  id=addone(id)
  call writefixl(nu_plt,id,'DMGR','Grain C production yiled [gC grain g-1 nonstrucal C]',GrainBiomGrowthYld_pft(NZ,NY,NX),101)
  id=addone(id)
  call writefixl(nu_plt,id,'DMRT','Root dry matter C production yield [gC root g-1 nonstrucal C]',RootBiomGrosYld_pft(NZ,NY,NX),101)
  
  if(iPlantNfixType_pft(NZ,NY,NX)/=0)then
    id=addone(id)
    call writefixl(nu_plt,id,'DMND','Canopy/root nodule bacteria dry matter C production yiled [gC g-1 nonstrucal C]',NoduGrowthYield_pft(NZ,NY,NX),101)
  endif
  end subroutine plant_biomyield_trait_disp

!------------------------------------------------------------------------------------------
  subroutine plant_biomstoich_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j,id

  write(nu_plt,*)('-',j=1,110)
  write(nu_plt,*)'ORGAN N AND P CONCENTRATIONS'
  id=0;id=addone(id)
  call writefixl(nu_plt,id,'CNWR','Ratio of root protein C to total root N [g protein C (g root N)-1]',rProteinC2RootN_pft(NZ,NY,NX),100)      
  id=addone(id)
  call writefixl(nu_plt,id,'CNWL','Ratio of leaf protein C to total leaf N [g protein C (g leaf N)-1]',rProteinC2LeafN_pft(NZ,NY,NX),100)    
  id=addone(id)
  call writefixl(nu_plt,id,'CNLF','Plant leaf maximum NC mass ratio [gN (gC)-1]',rNCLeaf_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CNSHE','Plant PetolSheth NC mass ratio [gN (gC)-1]',rNCSheath_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CNSTK','Plant stalk NC mass ratio [gN (gC)-1]',rNCStalk_pft(NZ,NY,NX),100)  
  if(is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX)))then  
!    id=addone(id)
!    call writefixl(nu_plt,id,'CNRTLIG','Plant lignified root NC mass ratio [gN (gC)-1]',rNCLigRoot_pft(NZ,NY,NX),100)
    id=addone(id)
    call writefixl(nu_plt,id,'CNLCRT','Live coarse root NC mass ratio [gN (gC)-1]',rECLiveCRoot_pft(ielmn,NZ,NY,NX),100)
    id=addone(id)
    call writefixl(nu_plt,id,'CNDCRT','Dead coarse root NC mass ratio [gN (gC)-1]',rECDeadCRoot_pft(ielmn,NZ,NY,NX),100)
  endif
  id=addone(id)
  call writefixl(nu_plt,id,'CNRSV','Plant stalk reserve NC mass ratio [gN (gC)-1]',rNCReserve_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CNHSK','Plant husk NC mass ratio [gN (gC)-1]',rNCHusk_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CNEAR','Plant ear NC mass ratio [gN (gC)-1]',rNCEar_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CNGR','Plant grain NC mass ratio [gN (gC)-1]',rNCGrain_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CNRT','Plant root NC mass ratio [gN (gC)-1]',rNCRoot_pft(NZ,NY,NX),100)
  if(iPlantNfixType_pft(NZ,NY,NX)/=0)then
    id=addone(id)  
    call writefixl(nu_plt,id,'CNND','Plant nodule NC mass ratio [gN (gC)-1]',rNCNodule_pft(NZ,NY,NX),100)
  endif
  id=addone(id)
  call writefixl(nu_plt,id,'CPWR','Ratio of total root protein C to total root P [g protein C (g root P)-1]',rProteinC2RootP_pft(NZ,NY,NX),100)      
  id=addone(id)
  call writefixl(nu_plt,id,'CPWL','Ratio of total leaf protein C to total P [g protein C (g leaf P)-1]',rProteinC2LeafP_pft(NZ,NY,NX),100)      
  id=addone(id)
  call writefixl(nu_plt,id,'CPLF','Plant leaf maximum PC mass ratio [gP (gC)-1]',rPCLeaf_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CPSHE','Plant PetolSheth PC mass ratio [gP (gC)-1]',rPCSheath_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CPSTK','Plant stalk PC mass ratio [gP (gC)-1]',rPCStalk_pft(NZ,NY,NX),100)
  if(is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX)))then  
!    id=addone(id)
!    call writefixl(nu_plt,id,'CPRTLIG','Plant lignified root PC mass ratio [gP (gC)-1]',rPCLigRoot_pft(NZ,NY,NX),100)
    id=addone(id)
    call writefixl(nu_plt,id,'CPLCRT','Live coarse root NC mass ratio [gN (gC)-1]',rECLiveCRoot_pft(ielmp,NZ,NY,NX),100)
    id=addone(id)
    call writefixl(nu_plt,id,'CPDCRT','Dead coarse root NC mass ratio [gN (gC)-1]',rECDeadCRoot_pft(ielmp,NZ,NY,NX),100)
  endif
  id=addone(id)
  call writefixl(nu_plt,id,'CPRSV','Plant stalk reserve PC mass ratio [gP (gC)-1]',rPCReserve_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CPHSK','Plant husk PC mass ratio [gP (gC)-1]',rPCHusk_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CPEAR','Plant ear PC mass ratio [gP (gC)-1]',rPCEar_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CPGR','Plant grain PC mass ratio [gP (gC)-1]',rPCGrain_pft(NZ,NY,NX),100)
  id=addone(id)
  call writefixl(nu_plt,id,'CPRT','Plant root PC mass ratio [gP (gC)-1]',rPCRootr_pft(NZ,NY,NX),100)
  if(iPlantNfixType_pft(NZ,NY,NX)/=0)then
    id=addone(id)  
    call writefixl(nu_plt,id,'CPND','Plant nodule PC mass ratio [gP (gC)-1]',rPCNoduler_pft(NZ,NY,NX),100)
  endif
  end subroutine plant_biomstoich_trait_disp

!------------------------------------------------------------------------------------------

  subroutine plant_optic_trait_disp(nu_plt,NZ,NY,NX)
  implicit none
  integer, intent(in) :: NZ,NY,NX
  integer, intent(in) :: nu_plt
  integer :: j,id

  write(nu_plt,*)('-',j=1,110)
  write(nu_plt,*)'OPTICAL PROPERTIES'
  id=0;id=addone(id)
  call writefixl(nu_plt,id,'ALBR','Leaf shortwave (0.3-2.5 um) radiation albedo [-]',RadSWLeafAlbedo_pft(NZ,NY,NX),80)
  id=addone(id)
  call writefixl(nu_plt,id,'ALBP','Leaf PAR (0.4-0.7 um) albedo [-]',RadPARLeafAlbedo_pft(NZ,NY,NX),80)
  id=addone(id)
  call writefixl(nu_plt,id,'TAUR','Leaf shortwave radiation (0.3-2.5 um) transmittance [-]',RadSWLeafTransmitance_pft(NZ,NY,NX),80)
  id=addone(id)
  call writefixl(nu_plt,id,'TAUP','Leaf PAR (0.4-0.7 um) transmittance [-]',RadPARLeafTransmitance_pft(NZ,NY,NX),80)
  end subroutine plant_optic_trait_disp

!------------------------------------------------------------------------------------------

  subroutine ReadPlantInfoNC(yeari,NHW,NHE,NVN,NVS)
  !
  !Description
  use netcdf
  use ncdio_pio
  use abortutils, only : endrun
  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer, intent(in) :: yeari    !forcing year information

  character(len=*), parameter :: subname='ReadPlantInfoNC'
  integer :: IDATE
  integer :: NPP(JY,JX)
  integer :: IYR,NX,NY,NZ,NN,NH1,NH2,NV1,NV2,NS
  integer :: NA        !number of PFT slots named for this year; the rest of the NS slots are fallow
  CHARACTER(len=2) :: CLIMATE
  integer :: ntopou,NTOPO
  type(file_desc_t) :: pftinfo_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar
  integer :: pft_dflag
  integer :: iyear,year
  integer :: ret,kk
  character(len=1) :: cc
  character(len=10) :: pft_gtype(JP)

  call PrintInfo('beg '//subname)
  IF(is_first_year .and. is_cold_run())THEN
    D9999: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        DO NZ=1,JP
          iDayPlanting_pft(NZ,NY,NX)      = -1E+06
          iYearPlanting_pft(NZ,NY,NX)     = -1E+06
          iDayPlantHarvest_pft(NZ,NY,NX)  = 1E+06
          iYearPlantHarvest_pft(NZ,NY,NX) = 1E+06
        enddo
      enddo
    ENDDO D9999

  ENDIF

  if (len_trim(pft_mgmt_in)==0)then
   !no plant management info to read
    D5995: DO NX=NHW,NHE
      DO  NY=NVN,NVS
        NP_col(NY,NX)=0
        DO NZ=1,NP0_col(NY,NX)
          DATAP(NZ,NY,NX)='NO'
        enddo
      enddo
    ENDDO D5995
  else

    call ncd_pio_openfile(pftinfo_nfid, pft_mgmt_in, ncd_nowrite)

    call check_var(pftinfo_nfid, 'pft_dflag', vardesc, readvar)

    if(.not. readvar)then
      call endrun('fail to find pft_dflag in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_dflag), &
      trim(mod_filename)//':: at line '//trim(int2str(__LINE__)))

    if (pft_dflag==-1)then
      !no pft data
      DO NX=NHW,NHE
        DO  NY=NVN,NVS
          NP_col(NY,NX)=0
          DO NZ=1,NP0_col(NY,NX)
            DATAP(NZ,NY,NX)='NO'
          enddo
        enddo
      ENDDO
    else
      if(pft_dflag==0)then
        ! constant pft data
        iyear=1
      else
        iyear=1
        ! Species intentionally follow the forcing year in all transient modes.
        ! Flexible planting/management separately follow the simulation year;
        ! do not make these year selections identical when forcing is cycled.
        DO while(.true.)
          call ncd_getvar(pftinfo_nfid,'year',iyear,year)
          if(year==yeari)exit
          iyear=iyear+1
        ENDDO
      endif
      !obtain the number of topographic units
      ntopou=get_dim_len(pftinfo_nfid, 'ntopou')
      if(first_topou)ntopou=1    
      DO NTOPO=1,ntopou
        call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
        call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NV1)
        call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH2)
        call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NV2)
        call ncd_getvar(pftinfo_nfid,'NZ',ntopo,NS)
        ! Legacy plant files may reserve more ntopou slots than they populate.
        ! Nonpositive bounds (including NetCDF fill values) end the populated records;
        ! exit before validating NZ or reading planting/management data for unused slots.
        ! Keep this sentinel exit when changing validation. Populated units may span
        ! multiple cells, so their bounds are checked against JX/JY, not ntopou.
        if(any((/NH1,NV1,NH2,NV2/)<=0))exit
        if(NS<0 .or. NS>JP)then
          write(iulog,*)'Invalid PFT count in '//subname,': topo unit=',NTOPO,' NZ=',NS,' capacity=',JP
          call endrun('PFT count must be between zero and JP in '//subname,__LINE__)
        endif

        NN=1
        ! Topographic units may span several cells; validate against grid dimensions.
        if(NH1<1 .or. NH2>JX .or. NH1>NH2 .or. &
           NV1<1 .or. NV2>JY .or. NV1>NV2)then
          write(iulog,*)'Invalid bounds in '//subname,': topo unit=',NTOPO, &
            ' NH1,NH2,NV1,NV2=',NH1,NH2,NV1,NV2,' JX,JY=',JX,JY
          call endrun('Invalid topographic-unit bounds in '//subname,__LINE__)
        endif
        D4995: DO NX=NH1,NH2
          D4990: DO NY=NV1,NV2
            NP0_col(NY,NX)=NS
            D4985: DO NZ=1,NS
              LSG_pft(NZ,NY,NX)=NN
            ENDDO D4985
          ENDDO D4990
        ENDDO D4995

        NA=0
        IF(NS.GT.0)THEN
          ! a vegetated topo unit
          call check_var(pftinfo_nfid, 'pft_type', vardesc, readvar)
          if(.not. readvar)then
            call endrun('fail to find pft_type in '//trim(mod_filename), __LINE__)
          endif
          call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, pft_gtype, &
            start = (/1,1,NTOPO,iyear/),count = (/len(pft_gtype(1)),JP,1,1/)), &
            trim(mod_filename)//':: at line '//trim(int2str(__LINE__)))

          ! A blank pft_type marks a fallow slot for this year. Named slots come first,
          ! so slots NA+1..NS are fallow and keep only the residue of their last crop.
          DO NZ=1,NS
            if(is_blank_pft(pft_gtype(NZ)))exit
            NA=NZ
          ENDDO
          DO NZ=NA+1,NS
            if(.not.is_blank_pft(pft_gtype(NZ)))then
              write(iulog,*)subname//': topo unit ',NTOPO,' year record ',iyear, &
                ' names PFT slot ',NZ,' after blank slot ',NA+1
              call endrun('fallow (blank) pft_type slots must follow the named slots in '//subname,__LINE__)
            endif
          ENDDO

          D4975: DO NX=NH1,NH2
            D4970: DO NY=NV1,NV2
              D4965: DO NZ=1,NA
                flag_active_pft(NZ,NY,NX)=.true.
                !modify pft name 
                IF(KoppenClimZone_col(NY,NX).GT.0)THEN
                  WRITE(CLIMATE,'(I2)')KoppenClimZone_col(NY,NX)
                  !the type of pft is specified by genra+Koppen climate zone
                  DATAX_pft(NZ)=trim(pft_gtype(NZ)(1:4)//CLIMATE)                  
                  !consider cases when koppen climate zone is zero  
                elseif(KoppenClimZone_col(NY,NX)==0)then
                  DATAX_pft(NZ)=''
                  DO kk=1,10
                    cc=pft_gtype(NZ)(kk:kk)
                    if(ichar(cc)==0)exit
                    DATAX_pft(NZ)(kk:kk)=cc
                  ENDDO
                ENDIF
              ENDDO D4965
              if(first_pft)then
                flag_active_pft(2:NS,NY,NX)=.false.
                NP0_col(NY,NX)=1
                NS=1
                NA=MIN(NA,1)
              endif  
            ENDDO D4970
          ENDDO D4975
        ENDIF

        D8995: DO NX=NH1,NH2
          D8990: DO NY=NV1,NV2
            call CheckPFTSlotChange(NTOPO,NA,NS,NY,NX)
            !NP_col slots grow; slots NP_col+1..NP0_col are fallow
            NP_col(NY,NX)=NA
            !PFT names for this year; their traits are read in readq
            D100: DO NZ=1,NA
              DATAP(NZ,NY,NX)=DATAX_pft(NZ)
            ENDDO D100
            !a fallow slot keeps the name of the PFT whose residue it holds
            D102: DO NZ=NA+1,NS
              if(DATAPI(NZ,NY,NX).GT.0)then
                DATAP(NZ,NY,NX)=pftss_tab(DATAPI(NZ,NY,NX))
              else
                DATAP(NZ,NY,NX)='NO'
              endif
            ENDDO D102

            D101: DO NZ=NS+1,JP
              DATAP(NZ,NY,NX)='NO'
            ENDDO D101
          ENDDO D8990
        ENDDO D8995
      ENDDO
    ENDIF 
    call ncd_pio_closefile(pftinfo_nfid)
  endif
  call PrintInfo('end '//subname)
  end subroutine ReadPlantInfoNC
!------------------------------------------------------------------------------------------

  logical function is_blank_pft(pft_name)
  !
  !Description
  !true when a pft_type entry is empty (blank or NetCDF null characters), i.e. a fallow slot
  implicit none
  character(len=*), intent(in) :: pft_name

  is_blank_pft=len_trim(pft_name)==0 .or. ichar(pft_name(1:1))==0
  end function is_blank_pft
!------------------------------------------------------------------------------------------

  subroutine CheckPFTSlotChange(NTOPO,NA,NS,NY,NX)
  !
  !Description
  !Check the PFT slots of a column before this year's PFT names are assigned.
  !Slots NA+1..NS turn fallow: their plants must already be dead, because fallow
  !slots get no growth, respiration or senescence. Slots 1..NA that change PFT while
  !holding standing dead must keep the same woody status, because the litter class
  !and decay rate of standing dead are taken from the current PFT traits.
  use ElmIDMod, only : ielmc
  implicit none
  integer, intent(in) :: NTOPO   !topographic unit
  integer, intent(in) :: NA      !number of named PFT slots this year
  integer, intent(in) :: NS      !number of PFT slots of the column
  integer, intent(in) :: NY,NX

  character(len=*), parameter :: subname='CheckPFTSlotChange'
  character(len=40) :: pft_lname
  character(len=64) :: koppen_climl
  character(len=3)  :: koppen_clims
  integer :: NZ,loc_old,loc_new
  logical :: woody_old,woody_new

  call PrintInfo('beg '//subname)
  DO NZ=NA+1,NS
    if(iPlantStateLive_pft(NZ,NY,NX).EQ.iTrue)then
      write(iulog,*)subname//': topo unit ',NTOPO,' column ',NY,NX,' PFT slot ',NZ, &
        ' (',trim(DATAP(NZ,NY,NX)),') is still alive at the start of its fallow year'
      call endrun('schedule a harvest that kills the plant before a fallow year in '//subname,__LINE__)
    endif
  ENDDO

  DO NZ=1,NA
    loc_old=DATAPI(NZ,NY,NX)
    if(loc_old.LE.0 .or. StandDeadStrutElms_pft(ielmc,NZ,NY,NX).LE.0._r8)cycle
    loc_new=get_pft_loc(KoppenClimZone_col(NY,NX),DATAX_pft(NZ)(1:6),pft_lname,koppen_climl,koppen_clims)
    if(loc_new.EQ.loc_old)cycle
    woody_old=iPlantTurnoverPattern_tab(loc_old).NE.0 .and. &
      is_plant_woody_vascular(iPlantRootProfile_tab(loc_old),iPlant2ndGrothPattern_tab(loc_old))
    woody_new=iPlantTurnoverPattern_tab(loc_new).NE.0 .and. &
      is_plant_woody_vascular(iPlantRootProfile_tab(loc_new),iPlant2ndGrothPattern_tab(loc_new))
    if(woody_old.NEQV.woody_new)then
      write(iulog,*)subname//': topo unit ',NTOPO,' column ',NY,NX,' PFT slot ',NZ,' changes from ', &
        trim(pftss_tab(loc_old)),' to ',trim(DATAX_pft(NZ)),' with standing dead C ', &
        StandDeadStrutElms_pft(ielmc,NZ,NY,NX)
      call endrun('woody and non-woody PFTs cannot share a slot holding standing dead in '//subname,__LINE__)
    endif
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine CheckPFTSlotChange

!------------------------------------------------------------------------------------------

  subroutine writefixl(nu_plt,id,varname, desc,value,width)
  implicit none
  integer,intent(in)  :: nu_plt
  integer, intent(in) :: id
  character(len=*), intent(in) :: varname
  character(len=*), intent(in) :: desc
  integer, intent(in) :: width
  real(r8),intent(in) :: value
  character(len=10) :: vname
  character(len=width) :: line
  character(len=256) :: lines
  line=desc
  line(width:width)=':'
  write(vname,*)varname
  write(lines,*)line,value
  write(nu_plt,'(I2,A,A10,A)')id,':',vname, '| '//trim(lines)
  end subroutine writefixl
!------------------------------------------------------------------------------------------

  subroutine writefixsl(nu_plt,desc,value,width,id)
  implicit none
  integer,intent(in)  :: nu_plt  
  integer, optional,intent(in) :: id
  character(len=*), intent(in) :: desc
  integer, intent(in) :: width
  character(len=*),intent(in) :: value

  character(len=width) :: line

  line=desc
  line(width:width)=': '
  if(present(id))then
    write(nu_plt,'(I2,A,A,A)')id,'| ',line,value   
  else
    write(nu_plt,'(A,A)')line,value   
  endif

  end subroutine writefixsl

!------------------------------------------------------------------------------------------

  subroutine writeafixl(nu_plt,id,varname,desc,values,width)
  implicit none
  integer,intent(in)  :: nu_plt
  integer, intent(in) :: id
  character(len=*), intent(in) :: varname
  character(len=*), intent(in) :: desc
  integer, intent(in) :: width
  real(r8),dimension(:),intent(in) :: values

  character(len=width) :: line
  character(len=256) :: lines
  character(len=10) :: vname

  line=desc
  line(width:width)=':'
  write(vname,*)varname
  write(lines,*)line,values
  write(nu_plt,'(I2,A,A10,A)')id,':',vname,'| '//trim(lines)
  end subroutine writeafixl
!------------------------------------------------------------------------------------------

  Subroutine ReadPlantTraitTable()

  implicit none
  character(len=*), parameter :: subname='ReadPlantTraitTable'

  call PrintInfo('beg '//subname)

  call ncd_pio_openfile(pft_nfid, pft_file_in, ncd_nowrite)
  call ncd_getvar(pft_nfid, 'pfts', pftss_tab)
  call ncd_getvar(pft_nfid,'pfts_long',pft_long_tab)
  call ncd_getvar(pft_nfid,'pfts_short',pft_short_tab)
  call ncd_getvar(pft_nfid,'koppen_clim_no',koppen_clim_ncode_tab)
  call ncd_getvar(pft_nfid,'koppen_clim_short',koppen_clim_short_tab)
  call ncd_getvar(pft_nfid,'koppen_clim_long',koppen_clim_long_tab)
  call ncd_getvar(pft_nfid, 'ICTYP', iPlantPhotosynsType_pft_tab)
  call ncd_getvar(pft_nfid,'IEBTYP',iEmbryophyteType_pft_tab)
  call ncd_getvar(pft_nfid, 'IGTYP', iPlantRootProfile_tab)  
  call ncd_getvar(pft_nfid, 'ISTYP', iPlantPhenolPattern_tab)
  call ncd_getvar(pft_nfid, 'IDTYP', iPlantDevelopPattern_tab)
  call ncd_getvar(pft_nfid, 'INTYP', iPlantNfixType_tab)
  call ncd_getvar(pft_nfid, 'IWTYP', iPlantPhenolType_tab)
  call ncd_getvar(pft_nfid, 'IPTYP', iPlantPhotoperiodType_tab)
  call ncd_getvar(pft_nfid, 'IBTYP', iPlantTurnoverPattern_tab)
  call ncd_getvar(pft_nfid, 'I2EXP', iPlant2ndGrothPattern_tab)
  call ncd_getvar(pft_nfid, 'IRTYP', iPlantGrainType_tab)
  call ncd_getvar(pft_nfid, 'MY',  Myco_tab)
  call ncd_getvar(pft_nfid,'STDENS',StemMassDensity_tab)
  call ncd_getvar(pft_nfid, 'ZTYPI', PlantInitThermoAdaptZone_tab)
  call ncd_getvar(pft_nfid, 'VCMX',  VmaxRubCarboxyRef_tab)
  call ncd_getvar(pft_nfid, 'VOMX',  VmaxRubOxyRef_tab)
  call ncd_getvar(pft_nfid, 'VCMX4', VmaxPEPCarboxyRef_tab)
  call ncd_getvar(pft_nfid, 'XKCO2', XKCO2_tab)
  call ncd_getvar(pft_nfid, 'XKO2',  XKO2_tab)
  call ncd_getvar(pft_nfid, 'ROOTMAGE', RootMatureAge_tab)
  call ncd_getvar(pft_nfid, 'XKCO24',Km4PEPCarboxy_tab)
  call ncd_getvar(pft_nfid, 'RUBP', LeafRubisco2Protein_tab)
  call ncd_getvar(pft_nfid, 'PEPC', LeafPEP2Protein_tab)
  call ncd_getvar(pft_nfid, 'ISNTYP',iPlantSnowIntercepType_tab)
  call ncd_getvar(pft_nfid, 'ETMX', SpecChloryfilAct_tab)
  call ncd_getvar(pft_nfid, 'CHL',  LeafChl2Protein_tab)
  call ncd_getvar(pft_nfid, 'fCHLMESO', fChlMesophyll_tab)
  call ncd_getvar(pft_nfid, 'FCO2', CanopyCi2CaRatio_pft_tab)
  call ncd_getvar(pft_nfid, 'ALBR', RadSWLeafAlbedo_tab)
  call ncd_getvar(pft_nfid, 'ALBP', CanopyPARalbedo_tab)
  call ncd_getvar(pft_nfid, 'TAUR', RadSWLeafTransmis_tab)
  call ncd_getvar(pft_nfid, 'TAUP', RadPARLeafTransmis_tab)
  call ncd_getvar(pft_nfid, 'XRNI', RefNodeInitRate_tab)
  call ncd_getvar(pft_nfid, 'XRLA', RefLeafAppearRate_tab)
  call ncd_getvar(pft_nfid, 'CTC',  TCChill4Seed_tab)
  call ncd_getvar(pft_nfid, 'VRNLI', VRNLI_tab)
  call ncd_getvar(pft_nfid, 'VRNXI', VRNXI_tab)
  call ncd_getvar(pft_nfid,'PhiMEAN',xylemPhi_mean_tab)
  call ncd_getvar(pft_nfid,'PhiMIN',xylemPhi_min_tab)
  call ncd_getvar(pft_nfid,'R95MAT',Radius95pctMature_tab)
  call ncd_getvar(pft_nfid,'PhiMAX',xylemPhi_max_tab)
  call ncd_getvar(pft_nfid, 'WDLF', rLen2WidthLeaf_tab)
  call ncd_getvar(pft_nfid, 'PB', NonstCMinConc2InitBranch_tab)
  call ncd_getvar(pft_nfid, 'GROUPX', GROUPX_tab)
  call ncd_getvar(pft_nfid, 'XTLI', ShootNodeNumAtPlanting_tab)
  call ncd_getvar(pft_nfid, 'XDL', CriticPhotoPeriod_tab)
  call ncd_getvar(pft_nfid, 'XPPD',PhotoPeriodSens_tab)
  call ncd_getvar(pft_nfid, 'SLA1', SLA1_tab)
  call ncd_getvar(pft_nfid, 'SSL1', PetolShethLen2Mass_tab)
  call ncd_getvar(pft_nfid, 'SNL1', NodeLenPergC_tab)
  call ncd_getvar(pft_nfid, 'CLASS', LeafAngleClass_tab)
  call ncd_getvar(pft_nfid, 'CFI', ClumpFactorInit_tab)
  call ncd_getvar(pft_nfid, 'ANGBR', BranchAngle_tab)
  call ncd_getvar(pft_nfid, 'ANGSH', PetolShethAngle_tab)
  call ncd_getvar(pft_nfid, 'STMX', MaxPotentSeedNumber_tab)
  call ncd_getvar(pft_nfid, 'SDMX', MaxSeedNumPerSite_tab)
  call ncd_getvar(pft_nfid, 'GRMX', SeedCMassMax_tab)
  call ncd_getvar(pft_nfid, 'GRDM', SeedCMass_tab)
  call ncd_getvar(pft_nfid, 'GRW2L',SeedWidth2LenRatio_tab)
  call ncd_getvar(pft_nfid, 'GFILL',GrainFillRate25C_tab)
  call ncd_getvar(pft_nfid, 'CYTOS',enh_cyto_tab)
  call ncd_getvar(pft_nfid,'RSSA',StalkAxialResist_tab)
  call ncd_getvar(pft_nfid, 'WTSTDI',StandingDeadInitC_tab)
  call ncd_getvar(pft_nfid, 'RRAD1M',Root1stMaxRadius_tab)
  call ncd_getvar(pft_nfid, 'RRAD2M',Root2ndMaxRadius_tab)
  call ncd_getvar(pft_nfid, 'PORT', RootPorosity_tab)
  call ncd_getvar(pft_nfid,'KLGMAX',KLigMax_tab)
  call ncd_getvar(pft_nfid,'KLGMM',KLigMM_tab)
  call ncd_getvar(pft_nfid,'CNLCRT',rNCLiveCRoot_tab)
  call ncd_getvar(pft_nfid,'CNDCRT',rNCDeadCRoot_tab)
  call ncd_getvar(pft_nfid,'CPLCRT',rPCLiveCRoot_tab)
  call ncd_getvar(pft_nfid,'CPDCRT',rPCDeadCRoot_tab)
  call ncd_getvar(pft_nfid, 'PR', NonstCMinCon2InitRoot_tab)
  call ncd_getvar(pft_nfid, 'RSRR', RootRadialResist_tab)
  call ncd_getvar(pft_nfid,'RVSR',RootVesselRadius_tab)
  call ncd_getvar(pft_nfid, 'ARSRA', AlphaVesselAxialResist_tab)
  call ncd_getvar(pft_nfid, 'PTSHT',ShutRutNonstElmntConducts_tab)
  call ncd_getvar(pft_nfid, 'RTFQ', FineRootBranchFreq_tab)
  call ncd_getvar(pft_nfid, 'MRTFQ',MediumRootBranchFreq_tab)
  call ncd_getvar(pft_nfid, 'UPMXZH',VmaxNH4Root_tab)
  call ncd_getvar(pft_nfid, 'UPKMZH',KmNH4Root_tab)
  call ncd_getvar(pft_nfid, 'UPMNZH',CMinNH4Root_tab)
  call ncd_getvar(pft_nfid, 'MOPHGEN',MorphogenBase_tab)
  call ncd_getvar(pft_nfid, 'UPMXZO', VmaxNO3Root_tab)
  call ncd_getvar(pft_nfid, 'UPKMZO', KmNO3Root_tab)
  call ncd_getvar(pft_nfid, 'UPMNZO', CminNO3Root_tab)
  call ncd_getvar(pft_nfid, 'UPMXPO', VmaxPO4Root_tab)
  call ncd_getvar(pft_nfid, 'UPKMPO', KmPO4Root_tab)
  call ncd_getvar(pft_nfid, 'UPMNPO', CMinPO4Root_tab)
  call ncd_getvar(pft_nfid, 'OSMO', CanOsmoPsi0pt_tab)
  call ncd_getvar(pft_nfid, 'RCS', RCS_tab)
  call ncd_getvar(pft_nfid, 'RSMX',CuticleResist_tab)
  call ncd_getvar(pft_nfid, 'DMLF', LeafBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMSHE',PetolShethBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMSTK',StalkBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMRSV',ReserveBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMHSK',HuskBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMEAR',EarBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMGR', GrainBiomGrowthYld_tab)
  call ncd_getvar(pft_nfid, 'DMRT', RootBiomGrosYld_tab)
  call ncd_getvar(pft_nfid, 'DMND', NoduGrowthYield_tab)
  call ncd_getvar(pft_nfid, 'CNWR', rProteinC2RootN_tab)
  call ncd_getvar(pft_nfid, 'CNWL', rProteinC2LeafN_tab)
  call ncd_getvar(pft_nfid, 'CNLF', rNCLeaf_tab)
  call ncd_getvar(pft_nfid, 'CNSHE', rNCSheath_tab)
  call ncd_getvar(pft_nfid, 'CNSTK', rNCStalk_tab)
  call ncd_getvar(pft_nfid,'CNRTLIG',rNCLigRoot_tab)
  call ncd_getvar(pft_nfid, 'CNRSV', rNCReserve_tab)
  call ncd_getvar(pft_nfid, 'CNHSK', rNCHusk_tab)
  call ncd_getvar(pft_nfid, 'CNEAR', rNCEar_tab)
  call ncd_getvar(pft_nfid, 'CNGR', rNCGrain_tab)
  call ncd_getvar(pft_nfid, 'CNRT', rNCRoot_tab)
  call ncd_getvar(pft_nfid, 'CNND', rNCNodule_tab)
  call ncd_getvar(pft_nfid,'CPRTLIG',rPCLigRoot_tab)
  call ncd_getvar(pft_nfid, 'CPWR', rProteinC2RootP_tab)  
  call ncd_getvar(pft_nfid, 'CPWL', rProteinC2LeafP_tab)    
  call ncd_getvar(pft_nfid, 'CPLF', rPCLeaf_tab)
  call ncd_getvar(pft_nfid, 'CPSHE',rPCSheath_tab)
  call ncd_getvar(pft_nfid, 'CPSTK', rPCStalk_tab)
  call ncd_getvar(pft_nfid, 'CPRSV', rPCReserve_tab)
  call ncd_getvar(pft_nfid, 'CPHSK', rPCHusk_tab)
  call ncd_getvar(pft_nfid, 'CPEAR', rPCEar_tab)
  call ncd_getvar(pft_nfid, 'CPGR', rPCGrain_tab)
  call ncd_getvar(pft_nfid, 'CPRT', rPCRootr_tab)
  call ncd_getvar(pft_nfid, 'CPND', rPCNoduler_tab)
  if(ats_cpl_mode)then
    call ncd_getvar(pft_nfid, 'IROOTJK', irootType_tab)
    call ncd_getvar(pft_nfid, 'CANHT', CanopyHeightLive_tab)
  endif
  call ncd_pio_closefile(pft_nfid)
  call PrintInfo('end '//subname)
  end Subroutine ReadPlantTraitTable

!------------------------------------------------------------------------------------------

  function get_pft_loc(koppen_def,pft_name,pft_lname,koppen_climl,koppen_clims)result(loc)
  !
  !!DESCRIPTION
  ! return the id of pft to be read
  use PlantTraitTableMod
  implicit none
  integer, intent(in) :: koppen_def            !koppen climate numerical code
  character(len=*), intent(in) :: pft_name     !pft short name
  character(len=40),intent(out):: pft_lname    !returning pft long name
  character(len=64),intent(out):: koppen_climl !koppen climate Description
  character(len=3), intent(out):: koppen_clims !koppen climate short code
  integer :: loc,loc1,len,k
  write(*,*)pft_name
  len=len_trim(pft_name)
  if(koppen_def==0)write(iulog,*)'Koppen climate not defined for the selected pft ',trim(pft_name)
  do 
    if(ichar(pft_name(len:len))==0 .or. pft_name(len:len)==' ')then
      len=len-1
    else
      exit
    endif
  enddo

  loc=1
  DO
    if(pftss_tab(loc)(1:len)==pft_name(1:len))exit
    loc=loc+1
  enddo

  loc1=1
  do 
    if(pft_name(1:4)==pft_short_tab(loc1))exit
    loc1=loc1+1
    if(loc1 > size(pft_short_tab))then
      exit      
    endif
  enddo
  if(loc1 <=  size(pft_short_tab))then
    pft_lname=pft_long_tab(loc1)
  else
    pft_lname=pft_name
  endif
  
  if(koppen_def==0)then
    koppen_climl='None'  
    koppen_clims='None'
    return
  endif
  loc1=1  
  do
    if(koppen_clim_ncode_tab(loc1)==pft_name(5:6))exit
    loc1=loc1+1
  enddo
  koppen_climl=koppen_clim_long_tab(loc1)
  koppen_clims=koppen_clim_short_tab(loc1)
  end function get_pft_loc
end module PlantInfoMod
