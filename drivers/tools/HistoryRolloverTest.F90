program HistoryRolloverTest
  ! Description:
  ! Write synthetic hourly history through a leap-year boundary, then verify
  ! the files in a separate process so unflushed NetCDF buffers cannot hide loss.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use HistFileMod, only: hist_addfld1d,hist_htapes_build,hist_update_hbuf,hist_htapes_wrapup, &
    hist_nhtfrq,hist_mfilt,hist_empty_htapes
  use GridConsts, only: bounds,JX,JY,JZ,JP,MaxNumRootAxes,NumOfPlantMorphUnits
  use EcoSIMCtrlMod, only: etimer,lverb
  use EcoSIMConfig, only: case_name
  use DebugToolMod, only: PrintInfo
  use abortutils, only: iulog
  use netcdf
  implicit none
  character(len=32) :: Mode
  real(r8), pointer :: Sample(:)
  integer :: Step
  logical :: RestartWrites

  call get_command_argument(1,Mode)
  if(trim(Mode)=='verify')then
    call VerifyFile('history_rollover.ecosim.h0.1980-01-00-00000.nc',1,8760)
    call VerifyFile('history_rollover.ecosim.h0.1980-12-30-00000.nc',8761,24)
    call VerifyFile('history_rollover.ecosim.h0.1981-01-00-00000.nc',8785,24)
    write(iulog,*)'PASS: all 8808 hourly samples and time bounds survive rollover'
  else
    if(trim(Mode)/='plain'.and.trim(Mode)/='restart')error stop 'invalid fixture mode'
    RestartWrites=trim(Mode)=='restart'
    lverb=.false.
    JX=1;JY=1;JZ=1;MaxNumRootAxes=1;NumOfPlantMorphUnits=1
    bounds%ngrid=1;bounds%ntopou=1;bounds%ncols=1;bounds%npfts=JP
    bounds%NHW=1;bounds%NHE=1;bounds%NVN=1;bounds%NVS=1
    bounds%begg=1;bounds%endg=1;bounds%begt=1;bounds%endt=1
    bounds%begc=1;bounds%endc=1;bounds%begp=1;bounds%endp=JP
    case_name='history_rollover'
    hist_nhtfrq=1;hist_mfilt=8760;hist_empty_htapes=.false.
    allocate(Sample(1));Sample=0._r8
    call etimer%Init(year0=1980,nyears=2)
    call etimer%setClock(dtime=3600._r8,nelapstep=0)
    call hist_addfld1d(fname='SAMPLE',units='1',avgflag='I', &
      long_name='consecutive hourly sample number',ptr_col=Sample)
    call hist_htapes_build()
    do Step=1,8808
      Sample=real(Step,r8)
      call hist_update_hbuf(bounds)
      call etimer%update_time_stamp()
      ! Exercise reopening a partial tape, plus restart writing at year rollover.
      call hist_htapes_wrapup(RestartWrites.and.(Step==8761.or.Step==8770.or.Step==8784), &
        Step==8808,bounds,etimer%its_a_new_year())
    enddo
    deallocate(Sample)
  endif
contains
!------------------------------------------------------------------------------------------
  subroutine VerifyFile(Filename,FirstSample,ExpectedCount)
    ! Description:
    ! Read persisted samples and interval endpoints after the writer exits.
    character(len=*), intent(in) :: Filename ! history file to inspect [-]
    integer, intent(in) :: FirstSample ! first expected hourly sample number [-]
    integer, intent(in) :: ExpectedCount ! expected stored sample count [-]
    character(len=*), parameter :: subname='VerifyFile'
    integer :: FileID,DimID,VarID,Count,K
    real(r8), allocatable :: Values(:,:),TimeBounds(:,:)
    call PrintInfo('beg '//subname)
    call CheckNC(nf90_open(Filename,nf90_nowrite,FileID),'open '//Filename)
    call CheckNC(nf90_inq_dimid(FileID,'time',DimID),'find time dimension')
    call CheckNC(nf90_inquire_dimension(FileID,DimID,len=Count),'read sample count')
    if(Count/=ExpectedCount)then
      write(iulog,*)Filename,': expected/actual records',ExpectedCount,Count
      error stop 'missing history records at year rollover'
    endif
    allocate(Values(1,Count),TimeBounds(2,Count))
    call CheckNC(nf90_inq_varid(FileID,'SAMPLE',VarID),'find sample variable')
    call CheckNC(nf90_get_var(FileID,VarID,Values),'read samples')
    call CheckNC(nf90_inq_varid(FileID,'time_bounds',VarID),'find time bounds')
    call CheckNC(nf90_get_var(FileID,VarID,TimeBounds),'read time bounds')
    do K=1,Count
      if(Values(1,K)/=real(FirstSample+K-1,r8))error stop 'missing or duplicated sample'
      if(abs(TimeBounds(1,K)-real(FirstSample+K-2,r8)/24._r8)>1.e-10_r8.or. &
         abs(TimeBounds(2,K)-real(FirstSample+K-1,r8)/24._r8)>1.e-10_r8) &
        error stop 'gap or overlap in history time bounds'
    enddo
    call CheckNC(nf90_close(FileID),'close history file')
    deallocate(Values,TimeBounds)
    call PrintInfo('end '//subname)
  end subroutine VerifyFile
!------------------------------------------------------------------------------------------
  subroutine CheckNC(Status,Context)
    ! Description:
    ! Fail the regression test on a NetCDF I/O error.
    integer, intent(in) :: Status ! NetCDF return code [-]
    character(len=*), intent(in) :: Context ! operation being checked [-]
    character(len=*), parameter :: subname='CheckNC'
    call PrintInfo('beg '//subname)
    if(Status/=nf90_noerr)then
      write(iulog,*)Context,': ',trim(nf90_strerror(Status))
      error stop 'history rollover fixture I/O failed'
    endif
    call PrintInfo('end '//subname)
  end subroutine CheckNC
end program HistoryRolloverTest
