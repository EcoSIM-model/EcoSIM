module ReadManagementMod
  !!
  ! Description
  ! code to read Soil management info
  !
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: endrun
  use fileUtil,      only: open_safe, int2str
  use minimathmod,   only: isLeap
  use DebugToolMod,  only: PrintInfo
  use GridConsts
  use FlagDataType
  use FertilizerDataType
  use ClimForcDataType
  use EcoSIMCtrlMod, only : lverb, first_topou,pft_mgmt_in
  use SoilWaterDataType
  use LandSurfDataType
  use EcoSIMCtrlDataType
  use EcosimConst
  use EcoSIMHistMod
  use IrrigationDataType
  use GridDataType
  use EcoSIMConfig
  USE ncdio_pio
  use netcdf
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: ReadManagementFiles
  public :: ReadSoilManagementForYear
  public :: ReadFire
  contains

!------------------------------------------------------------------------------------------
  subroutine ReadTillageFile(soilmgmt_nfid,FileTillage,NH1,NH2,NV1,NV2,dimtill,EventYear)
  !
  !read soil tillage information. It also includes events due to fire and drainage
  implicit none
  character(len=*), intent(in) :: FileTillage
  integer, intent(in) :: NH1,NH2,NV1,NV2
  type(file_desc_t), intent(in) :: soilmgmt_nfid
  character(len=*), intent(in) :: dimtill
  integer, optional, intent(in) :: EventYear !configured fire year; otherwise use the record year
  integer :: NY,NX
  integer  :: IPLOW,IDY
  real(r8) :: DPLOW,DY
  integer  :: LPY,IDY1,IDY2,IDY3
  integer :: kk
  logical :: readvar
  type(Var_desc_t) :: vardesc
  character(len=24) :: tillf(367)
  integer :: ntill

  ntill=get_dim_len(soilmgmt_nfid,dimtill)
  if(ntill>367)then
    call endrun('Not enough memory size for array tillf in '//trim(mod_filename),__LINE__)
  endif

  call check_var(soilmgmt_nfid, FileTillage, vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find tillage file '//trim(FileTillage)//' in '//trim(mod_filename), __LINE__)
  endif
  
  call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, tillf(1:ntill)),&
      trim(mod_filename))

  !
  !     DY=date DDMMYYYY
  !     IPLOW,DPLOW=intensity,depth of disturbance
  !     iSoilDisturbType_col=soil disturbance type 1-20:tillage,21=litter removal,22=fire,23-24=drainage
  !     DepzCorp_col=intensity (fire) or depth (tillage,drainage) of disturbance
  !
  kk=1
  do while(len_trim(tillf(kk)).GT.0)

    read(tillf(kk),'(I2,I2,I4)')IDY1,IDY2,IDY3
    ! Fire records may contain placeholder years; use the configured event calendar.
    if(present(EventYear))IDY3=EventYear
    LPY=0
    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDY=IDY1
    else
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif
    read(tillf(kk),*)DY,IPLOW,DPLOW
    if(lverb)then
      print*,tillf(kk)
      print*,idy1,idy2,idy3,IPLOW,DPLOW
    endif
    
    D8995: DO NX=NH1,NH2
      D8990: DO NY=NV1,NV2
        iSoilDisturbType_col(IDY,NY,NX) = IPLOW
        DepzCorp_col(IDY,NY,NX)         = DPLOW
      ENDDO D8990
    ENDDO D8995
    if(kk.eq.ntill)exit
    kk=kk+1    
  enddo
  end subroutine ReadTillageFile
!------------------------------------------------------------------------------------------

  subroutine ReadIrrigationFile(soilmgmt_nfid,FileIrrig,NH1,NH2,NV1,NV2)
  use abortutils, only : iulog

  implicit none
  character(len=*), intent(in) :: FileIrrig
  integer, intent(in) :: NH1,NH2,NV1,NV2
  type(file_desc_t), intent(in) :: soilmgmt_nfid

  integer :: NY,NX,LPY,J,JEN
  integer :: IDY1,IDY2,IDY3,I,IDY
  integer :: IDYS,IHRS,iIrrigOptX,JST,IDYE,IHRE
  real(r8):: DY
  real(r8) :: DST,DEN,CIRRX,RR,FIRRX
  real(r8) :: DIRRX,PHQX,CKAQX,RRH,WDPTHI,CCLQX
  real(r8) :: CMGQX,CNAQX,CSOQX
  real(r8) :: CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX
  real(r8) :: IrrigProperties(13),ExistingProperties(13)
  integer :: kk
  logical :: readvar
  type(Var_desc_t) :: vardesc
  character(len=128) :: irrigf(24)
  integer :: nirrig

  nirrig=get_dim_len(soilmgmt_nfid,'nirri')
  if(nirrig==0)return

  if(nirrig>24)then
    call endrun("not enough size for array irrigf in "//trim(mod_filename), __LINE__)
  endif

  call check_var(soilmgmt_nfid, FileIrrig, vardesc, readvar)

  if(.not. readvar)then
    call endrun('fail to find irrigation file '//trim(FileIrrig)//' in '//trim(mod_filename), __LINE__)
  endif

  IF(FileIrrig(1:4).EQ.'auto')THEN
!
!       AUTOMATED IRRIGATION
!
!       DST,DEN=start,end dates,hours DDMMHHHH
!       TODO: Verify this input format; DDMMHHHH may be a documentation typo.
!       The final four digits are used both as a year in isLeap(IDY3) and
!       as a time in IHRS/IHRE. Check the intended input convention before
!       revising the parsing or leap-year calculation for either endpoint.
!       iIrrigOptX=flag for irrigation criterion,0=SWC,1=canopy water potl
!       FIRRX=depletion of SWC from CIRRX to WP(iIrrigOpt_col=0),or minimum canopy
!       water potential(iIrrigOpt_col=1), to trigger irrigation
!       CIRRX= fraction of FC to which irrigation will raise SWC
!       DIRRX= depth to which water depletion and rewatering is calculated
!       WDPTHI=depth at which irrigation is applied
!       PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,
!       CSOQX,CCLQX=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl
!       concentration in irrigation water
!

    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, irrigf(1)),&
      trim(mod_filename))

    READ(irrigf(1),*)DST,DEN,iIrrigOptX,FIRRX,CIRRX,DIRRX,WDPTHI &
      ,PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX &
      ,CSOQX,CCLQX
    READ(irrigf(1),'(I2,I2,I4)')IDY1,IDY2,IDY3

    IF(lverb)then
      print*,irrigf(1)
      print*,IDY1,IDY2,IDY3,DEN,iIrrigOptX,FIRRX,CIRRX,DIRRX,WDPTHI &
        ,PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX &
        ,CSOQX,CCLQX
    endif

    LPY=0
!   idy1: day, idy2:mon, idy3:hour
    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDYS=IDY1
    else
      IDYS=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif
    IHRS=IDY3
    LPY=0

    IDY1=INT(DEN/1.0E+06)
    IDY2=INT(DEN/1.0E+04-IDY1*1.0E+02)
    IDY3=INT(DEN-(IDY1*1.0E+06+IDY2*1.0E+04))
    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDYE=IDY1
    else
      IDYE=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif
    IHRE=IDY3
!
!   TRANSFER INPUTS TO MODEL ARRAYS
!
    D7965: DO NX=NH1,NH2
      D7960: DO NY=NV1,NV2
        Lirri_auto_col(NY,NX) = .true.
        iIrrigOpt_col(NY,NX) = iIrrigOptX
        IIRRA(1,NY,NX)       = IDYS
        IIRRA(2,NY,NX)       = IDYE
        IIRRA(3,NY,NX)       = INT(IHRS/100)
        IIRRA(4,NY,NX) = INT(IHRE/100)
        FIRRA_col(NY,NX)   = FIRRX
        CIRRA_col(NY,NX)   = CIRRX
        DIRRA(1,NY,NX) = DIRRX    !depth
        DIRRA(2,NY,NX) = WDPTHI   !width
        D220: DO I = 1, 366
          ! EcoSIM_ZM: initialize every day, not uninitialized scheduled IDY.
          PHQ(I,NY,NX)                                    = PHQX
          NH4_irrig_mole_conc(I,NY,NX)                    = CN4QX/14.0_r8
          NO3_irrig_mole_conc(I,NY,NX)                    = CNOQX/14.0_r8
          H2PO4_irrig_mole_conc(I,NY,NX)                  = CPOQX/31.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Al,I,NY,NX)  = CALQX/27.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Fe,I,NY,NX)  = CFEQX/55.8_r8
          trcsalt_irrig_mole_conc_col(idsalt_Ca,I,NY,NX)  = CCAQX/40.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Mg,I,NY,NX)  = CMGQX/24.3_r8
          trcsalt_irrig_mole_conc_col(idsalt_Na,I,NY,NX)  = CNAQX/23.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_K,I,NY,NX)   = CKAQX/39.1_r8
          trcsalt_irrig_mole_conc_col(idsalt_SO4,I,NY,NX) = CSOQX/32.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Cl,I,NY,NX)  = CCLQX/35.5_r8
        ENDDO D220
      ENDDO D7960
    ENDDO D7965
  ELSE
!
!    print*,'   SCHEDULED IRRIGATION'
!
    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, irrigf(1:nirrig)),&
      trim(mod_filename))

    do kk=1,nirrig
      if(len_trim(irrigf(kk))==0)exit
!
!     DY,RR,JST,JEN=date DDMMYYYY,amount (mm),start and end hours
!     PHQX,CN4QX,CNOQX,CPOQX,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,
!     CSOQX,CCLQX=pH,NH4,NO3,H2PO4,Al,Fe,Ca,Mg,Na,K,SO4,Cl
!     concentration in irrigation water
!
      READ(irrigf(kk),*)DY,RR,JST,JEN,WDPTHI,PHQX,CN4QX,CNOQX,CPOQX &
          ,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,CSOQX,CCLQX
      READ(irrigf(kk),'(I2,I2,I4)')IDY1,IDY2,IDY3

      IF(lverb)then
        write(*,*)irrigf(kk)
        write(*,*)IDY1,IDY2,IDY3,DY,RR,JST,JEN,WDPTHI,PHQX,CN4QX,CNOQX,CPOQX &
            ,CALQX,CFEQX,CCAQX,CMGQX,CNAQX,CKAQX,CSOQX,CCLQX
      endif
      LPY=0
      IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
      IF(IDY2.EQ.1)then
        IDY=IDY1
      else
        IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
      endif
      if(JST<1 .or. JEN>24 .or. JST>JEN .or. RR<0._r8)then
        write(iulog,*)'Invalid scheduled irrigation in ',trim(FileIrrig), &
          ': record=',kk,' day=',IDY,' start,end hours=',JST,JEN,' amount [mm]=',RR
        call endrun('Scheduled irrigation requires 1<=start<=end<=24 and a nonnegative amount',__LINE__)
      endif
      ! A zero-water record must not replace another event's daily properties.
      if(RR==0._r8)cycle
      ! Chemistry and application depth are daily fields; all pulses must agree.
      IrrigProperties=(/PHQX,CN4QX/14._r8,CNOQX/14._r8,CPOQX/31._r8, &
        CALQX/27._r8,CFEQX/55.8_r8,CCAQX/40._r8,CMGQX/24.3_r8, &
        CNAQX/23._r8,CKAQX/39.1_r8,CSOQX/32._r8,CCLQX/35.5_r8,WDPTHI/)
      RRH=RR/(JEN-(JST-1))
      D8965: DO NX=NH1,NH2
        D8960: DO NY=NV1,NV2
          if(any(RRIG(:,IDY,NY,NX)>0._r8))then
            ExistingProperties=(/PHQ(IDY,NY,NX),NH4_irrig_mole_conc(IDY,NY,NX), &
              NO3_irrig_mole_conc(IDY,NY,NX),H2PO4_irrig_mole_conc(IDY,NY,NX), &
              trcsalt_irrig_mole_conc_col(idsalt_Al,IDY,NY,NX), &
              trcsalt_irrig_mole_conc_col(idsalt_Fe,IDY,NY,NX), &
              trcsalt_irrig_mole_conc_col(idsalt_Ca,IDY,NY,NX), &
              trcsalt_irrig_mole_conc_col(idsalt_Mg,IDY,NY,NX), &
              trcsalt_irrig_mole_conc_col(idsalt_Na,IDY,NY,NX), &
              trcsalt_irrig_mole_conc_col(idsalt_K,IDY,NY,NX), &
              trcsalt_irrig_mole_conc_col(idsalt_SO4,IDY,NY,NX), &
              trcsalt_irrig_mole_conc_col(idsalt_Cl,IDY,NY,NX),WDPTH(IDY,NY,NX)/)
            if(any(IrrigProperties/=ExistingProperties))then
              write(iulog,*)'Conflicting scheduled irrigation in ',trim(FileIrrig), &
                ': record=',kk,' day=',IDY,' NY,NX=',NY,NX
              call endrun('Same-day irrigation requires identical chemistry and application depth',__LINE__)
            endif
          endif
          D2535: DO J=1,24
            IF(J.GE.JST.AND.J.LE.JEN)RRIG(J,IDY,NY,NX)=RRIG(J,IDY,NY,NX)+RRH/1000.0_r8
          ENDDO D2535
!
!           TRANSFER INPUTS TO MODEL ARRAYS
!
          PHQ(IDY,NY,NX)                                    = PHQX
          NH4_irrig_mole_conc(IDY,NY,NX)                    = CN4QX/14.0_r8
          NO3_irrig_mole_conc(IDY,NY,NX)                    = CNOQX/14.0_r8
          H2PO4_irrig_mole_conc(IDY,NY,NX)                  = CPOQX/31.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Al,IDY,NY,NX)  = CALQX/27.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Fe,IDY,NY,NX)  = CFEQX/55.8_r8
          trcsalt_irrig_mole_conc_col(idsalt_Ca,IDY,NY,NX)  = CCAQX/40.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Mg,IDY,NY,NX)  = CMGQX/24.3_r8
          trcsalt_irrig_mole_conc_col(idsalt_Na,IDY,NY,NX)  = CNAQX/23.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_K,IDY,NY,NX)   = CKAQX/39.1_r8
          trcsalt_irrig_mole_conc_col(idsalt_SO4,IDY,NY,NX) = CSOQX/32.0_r8
          trcsalt_irrig_mole_conc_col(idsalt_Cl,IDY,NY,NX)  = CCLQX/35.5_r8
          WDPTH(IDY,NY,NX)=WDPTHI
        ENDDO D8960
      ENDDO D8965
    enddo
  ENDIF
  end subroutine ReadIrrigationFile

!------------------------------------------------------------------------------------------

  subroutine ReadFertlizerFile(soilmgmt_nfid,FertFile,NH1,NH2,NV1,NV2)

  use abortutils, only : iulog
  implicit none
  character(len=*), intent(in) :: FertFile
  integer,intent(in) :: NH1,NH2,NV1,NV2
  type(file_desc_t), intent(in) :: soilmgmt_nfid

  integer :: NY,NX,LPY
  integer :: IDY1,IDY2,IDY3,IDY
  real(r8) :: AppDepth,DDMMYYYY
  real(r8) :: NH3Soil,UreaSoil,NO3Soil,NH4Band,NH3Band,UreaBand,NO3Band,NH4Soil
  real(r8) :: MonocalciumPhosphateSoil,MonocalciumPhosphateBand,hydroxyapatite,Gypsum,LimeStone
  real(r8) :: PlantResC,PlantResN,PlantResP,ManureC,ManureN,ManureP
  real(r8) :: BandWidth,PO4Soil,PO4Band
  integer  :: IR0,IR1,IR2
  integer  :: kk
  logical :: readvar
  type(Var_desc_t) :: vardesc
  character(len=128) :: fertf(12)
  integer :: nfert
  integer :: FirstRecord(366)

  nfert=get_dim_len(soilmgmt_nfid,'nfert')
  if(nfert==0)return
  if(nfert>12)then
    call endrun('Not enough size for array fertf in '//trim(mod_filename), __LINE__)
  endif
!
!     DY=date DDMMYYYY
!     *A,*B=broadGypsumt,banded fertilizer application
!     Z4,Z3,ZU,ZO=NH4,NH3,urea,NO3
!     PM*,PH*=Ca(H2PO4)2,apatite
!     LimeStone,Gypsum=CaCO3,CaSO4
!     *1,*2=litter,manure amendments
!     RSC,RSN,RSC=amendment C,N,P content
!     AppDepth=application depth
!     BandWidth=band row width
!     IRO,IR1,IR2=fertilizer,litter,manure type
!

  call check_var(soilmgmt_nfid, FertFile, vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find fertilizer '//trim(FertFile)//' in '//trim(mod_filename), __LINE__)
  endif

  call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, fertf(1:nfert)),&
      trim(mod_filename))

  FirstRecord=0
  do kk=1,nfert
    if(len_trim(fertf(kk))==0)exit
    !all nutrients are in g/m2
    READ(fertf(kk),*)DDMMYYYY,NH4Soil,NH3Soil,UreaSoil,NO3Soil,NH4Band,NH3Band,UreaBand,NO3Band &
      ,MonocalciumPhosphateSoil,MonocalciumPhosphateBand,hydroxyapatite,LimeStone,Gypsum &
      ,PlantResC,PlantResN,PlantResP,ManureC,ManureN,ManureP,AppDepth &
      ,BandWidth,PO4Soil,PO4Band,IR0,IR1,IR2

    LPY=0
    IDY1=INT(DDMMYYYY/1.0E+06_r8)
    !return for bad fertilization data
    if(IDY1==0)return
    IDY2=INT(DDMMYYYY/1.0E+04_r8-IDY1*1.0E+02_r8)
    IDY3=INT(DDMMYYYY-(IDY1*1.0E+06_r8+IDY2*1.0E+04_r8))       
     
    IF(LVERB)then
      print*,fertf(kk)
      PRINT*,IDY1,IDY2,IDY3,NH4Soil,NH3Soil,UreaSoil,NO3Soil,NH4Band,NH3Band,UreaBand,NO3Band &
        ,MonocalciumPhosphateSoil,MonocalciumPhosphateBand,hydroxyapatite,LimeStone,Gypsum &
        ,PlantResC,PlantResN,PlantResP,ManureC,ManureN,ManureP,AppDepth &
        ,BandWidth,PO4Soil,PO4Band,IR0,IR1,IR2
    endif

    IF(isLeap(IDY3).and.IDY2.GT.2)LPY=1
    IF(IDY2.EQ.1)then
      IDY=IDY1
    else
      IDY=30*(IDY2-1)+ICOR(IDY2-1)+IDY1+LPY
    endif

    ! One combined amendment record is applied at noon for each day.
    if(IDY<1 .or. IDY>size(FirstRecord))then
      write(iulog,*)'Invalid fertilizer day in ',trim(FertFile),': record=',kk,' day=',IDY
      call endrun('Fertilizer day is outside the daily schedule',__LINE__)
    endif
    if(FirstRecord(IDY)/=0)then
      write(iulog,*)'Duplicate fertilizer day in ',trim(FertFile),': day=',IDY, &
        ' records=',FirstRecord(IDY),kk
      call endrun('Combine same-day fertilizer and amendment inputs into one record',__LINE__)
    endif
    FirstRecord(IDY)=kk

    D8985: DO NX=NH1,NH2
      D8980: DO NY=NV1,NV2
!
!         ENTER AMENDMENTS INTO MODEL ARRAYS
!
!         NH4,NH3,UREA,NO3 BROADCAST (A) AND BANDED (B)
!
        FERT(ifert_N_nh4,IDY,NY,NX)       = NH4Soil        !NH4 broadcast
        FERT(ifert_N_nh3,IDY,NY,NX)       = NH3Soil        !NH3 broadcast
        FERT(ifert_N_urea,IDY,NY,NX)      = UreaSoil       !Urea broadcast
        FERT(ifert_N_no3,IDY,NY,NX)       = NO3Soil        !NO3 broadcast
        FERT(ifert_N_nh4_band,IDY,NY,NX)  = NH4Band        !NH4 band
        FERT(ifert_N_nh3_band,IDY,NY,NX)  = NH3Band        !NH3 band
        FERT(ifert_N_urea_band,IDY,NY,NX) = UreaBand       !Urea band
        FERT(ifert_N_no3_band,IDY,NY,NX)  = NO3Band        !NO3 band
!
!         MONOCALCIUM PHOSPHATE OR HYDROXYAPATITE BROADCAST (A)
!         AND BANDED (B)
!
        FERT(ifert_P_Ca_H2PO4_2_soil,IDY,NY,NX) = MonocalciumPhosphateSoil        !PM broadcast
        FERT(ifert_P_Ca_H2PO4_2_band,IDY,NY,NX) = MonocalciumPhosphateBand        !PM band
        FERT(ifert_P_apatite,IDY,NY,NX)         = hydroxyapatite                  !(Ca5(PO4)3OH)
!
!         LIME AND GYPSUM
!
        FERT(ifert_Ca_lime,IDY,NY,NX)   = LimeStone
        FERT(ifert_Ca_gypsum,IDY,NY,NX) = Gypsum
!
!         PLANT AND ANIMAL RESIDUE C, N AND P
!
        FERT(ifert_plant_resC,IDY,NY,NX)  = PlantResC
        FERT(ifert_plant_resN,IDY,NY,NX)  = PlantResN
        FERT(ifert_plant_resP,IDY,NY,NX)  = PlantResP
        FERT(ifert_plant_manuC,IDY,NY,NX) = ManureC
        ! Retain input provenance even if grazing later merges manure returns.
        ExternalManureC_col(IDY,NY,NX)=0._r8
        if(IR2<imanure_grazing)ExternalManureC_col(IDY,NY,NX)=ManureC
        FERT(ifert_plant_manuN,IDY,NY,NX) = ManureN
        FERT(ifert_plant_manuP,IDY,NY,NX) = ManureP
        FERT(ifert_PO4_soil,IDY,NY,NX)    = PO4Soil
        FERT(ifert_PO4_band,IDY,NY,NX)    = PO4Band
!
!         DEPTH AND WIDTH OF APPLICATION
!
        FDPTH(IDY,NY,NX) = AppDepth    !depth
        ROWI(IDY,NY,NX)  = BandWidth       !width
!
!         TYPE OF FERTILIZER,PLANT OR ANIMAL RESIDUE
!
        IYTYP(iamendtyp_fert,IDY,NY,NX)     = IR0
        IYTYP(iAmendtyp_plantRes,IDY,NY,NX) = IR1
        IYTYP(iAmendtyp_Manure,IDY,NY,NX)   = IR2
      ENDDO D8980
    ENDDO D8985
  enddo
  end subroutine ReadFertlizerFile

!------------------------------------------------------------------------------------------

  subroutine ReadSoilManagementForYear(simulation_year,forcing_year)
  ! Select only the schedule coordinate. Event doses, units and processes
  ! remain in the existing reader; missing years still stop explicitly.
  use EcoSIMCtrlMod, only : soil_mgmt_use_model_year
  implicit none
  integer, intent(in) :: simulation_year,forcing_year
  if(soil_mgmt_use_model_year)then
    call ReadManagementFiles(simulation_year)
  else
    call ReadManagementFiles(forcing_year)
  endif
  end subroutine ReadSoilManagementForYear

  subroutine ReadManagementFiles(yeari)
  !
  !DESCRIPTION
  !read in management information for soil: tillage
  !fertilizer and irrigation
  use EcoSIMCtrlMod, only : soil_mgmt_in
  implicit none
  integer, intent(in) :: yeari
  character(len=*), parameter :: subname='ReadManagementFiles'
  integer :: NH1,NV1,NH2,NV2
  type(file_desc_t) :: soilmgmt_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar
  integer :: ntopou
  integer :: NTOPO
  character(len=10) :: fertf
  character(len=10) :: tillf
  character(len=10) :: irrigf
  integer :: iyear,year,nyears
!
!   NH1,NV1,NH2,NV2=N,W and S,E corners of landscape unit
!   DATA1(8),DATA1(5),DATA1(6)=disturbance,fertilizer,irrigation files
!   PREFIX=path for files in current or higher level directory
  call PrintInfo('beg '//subname)
  ! Reset annual irrigation modes so absent selections cannot retain old settings.
  Lirri_auto_col=.false.
  call ncd_pio_openfile(soilmgmt_nfid, soil_mgmt_in, ncd_nowrite)
  nyears=get_dim_len(soilmgmt_nfid, 'year')
  if(nyears==0)then
    !no fertilization, tillage, or irrigation
    call ncd_pio_closefile(soilmgmt_nfid)
    return
  endif
  ! EcoSIM_ZM: stop at the declared year count rather than reading beyond
  ! the coordinate when a forcing year is absent from management input.
  DO iyear=1,nyears
    call ncd_getvar(soilmgmt_nfid,'year',iyear,year)
    if(year==yeari)exit
  ENDDO
  if(iyear>nyears)then
    call ncd_pio_closefile(soilmgmt_nfid)
    call endrun('management forcing year not present: '//trim(int2str(yeari)), __LINE__)
  endif

  ntopou=get_dim_len(soilmgmt_nfid, 'ntopou')
  if(ntopou==0)then
    call ncd_pio_closefile(soilmgmt_nfid)
    return
  endif
  if(first_topou)ntopou=1  
  DO NTOPO=1,ntopou
    call ncd_getvar(soilmgmt_nfid,'NH1',ntopo,NH1)
    call ncd_getvar(soilmgmt_nfid,'NV1',ntopo,NV1)
    call ncd_getvar(soilmgmt_nfid,'NH2',ntopo,NH2)
    call ncd_getvar(soilmgmt_nfid,'NV2',ntopo,NV2)

    if(any((/NH1,NH2,NV1,NV2/)<0))THEN
      call endrun('something wrong in NHX or NHX indices on file '//trim(soil_mgmt_in)&
        //' in '//trim(mod_filename), __LINE__)
    endif

    call check_var(soilmgmt_nfid, 'fertf', vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find fertf in '//trim(mod_filename), __LINE__)
    endif

    ! EcoSIM_ZM: route the schedule for this unit, not the final unit count.
    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, fertf, &
      start = (/1,ntopo,iyear/),count = (/len(fertf),1/)), &
      trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
    !!!!
    call check_var(soilmgmt_nfid, 'tillf', vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find tillf in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, tillf, &
      start = (/1,ntopo,iyear/),count = (/len(tillf),1/)), &
      trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
    !!!!
    call check_var(soilmgmt_nfid, 'irrigf', vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find irrigf in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, irrigf, &
      start = (/1,ntopo,iyear/),count = (/len(irrigf),1/)), &
      trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
!
!   READ TILLAGE INPUT FILE
!
    IF(trim(tillf).NE.'NO')THEN
      if(lverb)print*,'ReadTillageFile'
      call ReadTillageFile(soilmgmt_nfid,tillf,NH1,NH2,NV1,NV2,'ntill')
    ENDIF
!
!   READ FERTLIZER INPUT FILE
!
    IF(trim(fertf).NE.'NO')THEN
      if(lverb)print*,'ReadFertlizerFile'
      call ReadFertlizerFile(soilmgmt_nfid,fertf,NH1,NH2,NV1,NV2)
    ENDIF
!
!   READ IRRIGATION INPUT FILE
!
    Lirri_auto_col(NV1:NV2,NH1:NH2)=.false.
    IF(trim(irrigf).NE.'NO')THEN
      if(lverb)print*,'ReadIrrigationFile'
      call ReadIrrigationFile(soilmgmt_nfid,irrigf,NH1,NH2,NV1,NV2)
    ENDIF
  ENDDO
  call ncd_pio_closefile(soilmgmt_nfid)
  call PrintInfo('end '//subname)
  end subroutine ReadManagementFiles

!------------------------------------------------------------------------------------------

  subroutine ReadFire(yearc,fire_event_entry,NHW,NHE,NVN,NVS)
  use EcoSIMCtrlMod, only : soil_mgmt_in  
  use abortutils, only : iulog
  implicit none
  integer, intent(in) :: yearc
  character(len=21), intent(in) :: fire_event_entry
  integer, intent(in) :: NHW,NHE,NVN,NVS
  type(file_desc_t) :: soilmgmt_nfid
  character(len=10) :: firef
  integer :: ntopou,ntopou_plant,ntopo,pair_end
  logical :: readvar
  type(Var_desc_t) :: vardesc
  type(file_desc_t) :: pftinfo_nfid

  integer :: NH1,NH2,NV1,NV2

  call ncd_pio_openfile(soilmgmt_nfid, soil_mgmt_in, ncd_nowrite)
  call ncd_pio_openfile(pftinfo_nfid, pft_mgmt_in, ncd_nowrite)

  ntopou=get_dim_len(soilmgmt_nfid, 'ntopou')
  ntopou_plant=get_dim_len(pftinfo_nfid, 'ntopou')
  ! Both management files must describe the same number of topographic units.
  ! Check the full dimensions before first_topou restricts either fire loop.
  if(ntopou/=ntopou_plant)then
    write(iulog,*)'ReadFire: ntopou=',ntopou,' in soil file ',trim(soil_mgmt_in)
    write(iulog,*)'ReadFire: ntopou=',ntopou_plant,' in plant file ',trim(pft_mgmt_in)
    call endrun('ReadFire requires matching ntopou in soil and plant management files',__LINE__)
  endif

  if(ntopou==0)then
    call ncd_pio_closefile(soilmgmt_nfid)
    call ncd_pio_closefile(pftinfo_nfid)
    return
  endif
  if(first_topou)ntopou=1  
  !
  !fire disturbance to soil  
  pair_end = index(fire_event_entry, ',')  
  DO NTOPO=1,ntopou
    call ncd_getvar(soilmgmt_nfid,'NH1',ntopo,NH1)
    call ncd_getvar(soilmgmt_nfid,'NV1',ntopo,NV1)
    call ncd_getvar(soilmgmt_nfid,'NH2',ntopo,NH2)
    call ncd_getvar(soilmgmt_nfid,'NV2',ntopo,NV2)

    if(any((/NH1,NH2,NV1,NV2/)<0))THEN
      call endrun('something wrong in NHX or NHX indices on file '//trim(soil_mgmt_in)&
        //' in '//trim(mod_filename), __LINE__)
    endif

    call check_var(soilmgmt_nfid, trim(fire_event_entry(pair_end+1:)), vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find soil fire mgmt in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(soilmgmt_nfid%fh, vardesc%varid, firef, &
      start = (/1,ntopo/),count = (/len(firef),1/)), &
      trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
    call ReadTillageFile(soilmgmt_nfid,firef,NH1,NH2,NV1,NV2,'nfire',EventYear=yearc)
  ENDDO  
  call ncd_pio_closefile(soilmgmt_nfid)

  !fire disturbance to plants    
  DO NTOPO=1,ntopou
    call ncd_getvar(pftinfo_nfid,'NH1',ntopo,NH1)
    call ncd_getvar(pftinfo_nfid,'NV1',ntopo,NV1)
    call ncd_getvar(pftinfo_nfid,'NH2',ntopo,NH2)
    call ncd_getvar(pftinfo_nfid,'NV2',ntopo,NV2)

    if(any((/NH1,NH2,NV1,NV2/)<0))THEN
      call endrun('something wrong in NHX or NHX indices on file '//trim(pft_mgmt_in)&
        //' in '//trim(mod_filename), __LINE__)
    endif

    call check_var(pftinfo_nfid, trim(fire_event_entry(1:pair_end-1)), vardesc, readvar)
    if(.not. readvar)then
      call endrun('fail to find plant fire mgmt in '//trim(mod_filename), __LINE__)
    endif

    call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, firef, &
      start = (/1,ntopo/),count = (/len(firef),1/)), &
      trim(mod_filename)//'::at line '//trim(int2str(__LINE__)))
    
    call ReadPlantFireMgmt(yearc,pftinfo_nfid,firef,NH1,NH2,NV1,NV2)
  ENDDO

  call ncd_pio_closefile(pftinfo_nfid)

  end subroutine ReadFire

!------------------------------------------------------------------------------------------
  subroutine ReadPlantFireMgmt(yearc,pftinfo_nfid,fireFile,NH1,NH2,NV1,NV2)
  use PlantMgmtDataType 
  implicit none
  integer, intent(in) :: yearc
  type(file_desc_t), intent(in) :: pftinfo_nfid
  character(len=*), intent(in) :: fireFile
  integer, intent(in) :: NH1,NH2,NV1,NV2
  
  character(len=128) :: firefstr(5),tstr
  integer :: maxpfts
  logical :: readvar
  type(Var_desc_t) :: vardesc
  real(r8) :: DY,ECUT11,ECUT12,ECUT13,ECUT14,ECUT21,ECUT22,ECUT23
  real(r8) :: ECUT24,HCUT,PCUT
  integer :: LPY,IDX,IMO,IYR,IDY,ICUT,JCUT
  integer :: NY,NX,NZ

  maxpfts=get_dim_len(pftinfo_nfid,'maxpfts')

  call check_var(pftinfo_nfid, fireFile, vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find fire file '//trim(fireFile)//' in '//trim(mod_filename), __LINE__)
  endif
  
  call check_ret(nf90_get_var(pftinfo_nfid%fh, vardesc%varid, firefstr(1:maxpfts)),&
      trim(mod_filename))

  DO NX=NH1,NH2
    DO NY=NV1,NV2
      DO NZ=1,NP_col(NY,NX)
        if(len_trim(firefstr(NZ))==0)cycle
        tstr=trim(firefstr(NZ))

        !read day/month/year when management occurs
        read(tstr,'(I2,I2,I4)')IDX,IMO,IYR
!        write(*,*)TSTR
        READ(TSTR,*)DY,ICUT,JCUT,HCUT,PCUT,ECUT11,ECUT12,ECUT13,&
            ECUT14,ECUT21,ECUT22,ECUT23,ECUT24
        !"17050000,2,0,0.06,0.0,0.9,0.3,0.0,0.1,0.9,0.9,0.0,0.5               
        LPY=0
        ! FireEvents determines the execution year, regardless of the record year.
        if(isLeap(yearc) .and. IMO.GT.2)LPY=1
        !obtain the ordinal day
        IF(IMO.EQ.1)then
          IDY=IDX
        else
          IDY=30*(IMO-1)+ICOR(IMO-1)+IDX+LPY
        endif

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
        FracBiomHarvsted(iHarvst_col,iplthvst_leaf,NZ,IDY,NY,NX)        = ECUT21
        FracBiomHarvsted(iHarvst_col,iplthvst_finenonleaf,NZ,IDY,NY,NX) = ECUT22
        FracBiomHarvsted(iHarvst_col,iplthvst_stalk,NZ,IDY,NY,NX)       = ECUT23
        FracBiomHarvsted(iHarvst_col,iplthvst_stdead,NZ,IDY,NY,NX)      = ECUT24        
      ENDDO        
    ENDDO
  ENDDO

  end subroutine ReadPlantFireMgmt
end module ReadManagementMod
