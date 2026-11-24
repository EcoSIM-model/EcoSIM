module timings
! Wall-clock profiling for the serial driver. Nested regions are inclusive.
! Module accumulators are not thread safe; call outside parallel regions.
  use iso_fortran_env, only : int64
  use data_kind_mod, only : r8 => DAT_KIND_R8
  implicit none
  private

  type, public :: timer_stamp_type
    integer(int64) :: count = 0_int64
    integer(int64) :: rate = 0_int64
    integer(int64) :: maximum = 0_int64
  end type timer_stamp_type

  type :: timer_entry_type
    character(len=:), allocatable :: name
    integer(int64) :: calls = 0_int64, step_calls = 0_int64
    real(r8) :: total = 0._r8, step_total = 0._r8
    real(r8) :: minimum = huge(0._r8), maximum = 0._r8
  end type timer_entry_type

  type(timer_entry_type), allocatable :: entries(:)
  type(timer_stamp_type) :: run_stamp
  integer :: used = 0, detail_unit, summary_unit
  integer(int64) :: step_number = 0_int64
  logical :: initialized = .false., detail_enabled = .false.
  public :: init_timer, start_timer, end_timer, end_timer_loop, finalize_timer
  public :: timer_elapsed, flush_timer_detail, timer_launch_label
contains

  function timer_launch_label() result(label)
! Local wall time as YYYYMMDD-HHMMSS, used to name the detail file.
    character(len=15) :: label
    character(len=8) :: date
    character(len=10) :: time
    call date_and_time(date=date, time=time)
    label = date//'-'//time(1:6)
  end function timer_launch_label

  subroutine init_timer(outdir, detail, run_start)
! run_start labels timing/time-<run_start>.txt; defaults to the current time.
    character(len=*), intent(in) :: outdir
    logical, optional, intent(in) :: detail
    character(len=*), optional, intent(in) :: run_start
    character(len=:), allocatable :: label
    integer :: status
    if (initialized) error stop 'timings: already initialized'
    if (allocated(entries)) deallocate(entries)
    allocate(entries(16))
    used = 0
    step_number = 0_int64
    detail_enabled = .false.
    if (present(detail)) detail_enabled = detail
    call execute_command_line('mkdir -p '//shell_quote(trim(outdir)//'/timing'), exitstat=status)
    if (status /= 0) error stop 'timings: cannot create output directory'
    open(newunit=summary_unit, file=trim(outdir)//'/timing/summary.csv', &
         status='replace', action='write', iostat=status)
    if (status /= 0) error stop 'timings: cannot open summary'
    if (detail_enabled) then
      if (present(run_start)) then
        label = trim(run_start)
      else
        label = timer_launch_label()
      endif
      open(newunit=detail_unit, file=trim(outdir)//'/timing/time-'//label//'.txt', &
           status='replace', action='write', iostat=status)
      if (status /= 0) error stop 'timings: cannot open timestep output'
      write(detail_unit,'(A)') 'step,timer,calls,seconds'
    endif
    initialized = .true.
    call start_timer(run_stamp)
  end subroutine init_timer

  subroutine start_timer(stamp)
    type(timer_stamp_type), intent(out) :: stamp
    call system_clock(stamp%count, stamp%rate, stamp%maximum)
    if (stamp%rate <= 0_int64 .or. stamp%count < 0_int64) &
      error stop 'timings: wall clock unavailable'
  end subroutine start_timer

  pure function timer_elapsed(stamp, stop_count) result(seconds)
! Also used with synthetic counts to test rollover without waiting for it.
    type(timer_stamp_type), intent(in) :: stamp
    integer(int64), intent(in) :: stop_count
    real(r8) :: seconds
    integer(int64) :: ticks
    if (stop_count >= stamp%count) then
      ticks = stop_count - stamp%count
    else
! Split the cycle to avoid integer overflow at maximum+1.
      ticks = stamp%maximum - stamp%count
      seconds = (real(ticks,r8) + real(stop_count,r8) + 1._r8) / real(stamp%rate,r8)
      return
    endif
    seconds = real(ticks,r8) / real(stamp%rate,r8)
  end function timer_elapsed

  subroutine end_timer(timer_name, stamp)
    character(len=*), intent(in) :: timer_name
    type(timer_stamp_type), intent(in) :: stamp
    type(timer_entry_type), allocatable :: grown(:)
    integer(int64) :: stop_count
    integer :: idx
    real(r8) :: seconds
    if (.not.initialized) error stop 'timings: not initialized'
    if (stamp%rate <= 0_int64) error stop 'timings: invalid stamp'
    call system_clock(stop_count)
    if (stop_count < 0_int64) error stop 'timings: wall clock unavailable'
    seconds = timer_elapsed(stamp, stop_count)
    do idx = 1, used
      if (entries(idx)%name == trim(timer_name)) exit
    enddo
    if (idx > used) then
      if (used == size(entries)) then
        allocate(grown(2*size(entries)))
        grown(1:used) = entries(1:used)
        call move_alloc(grown, entries)
      endif
      used = used + 1
      entries(idx)%name = trim(timer_name)
    endif
    entries(idx)%calls = entries(idx)%calls + 1_int64
    entries(idx)%total = entries(idx)%total + seconds
    entries(idx)%minimum = min(entries(idx)%minimum, seconds)
    entries(idx)%maximum = max(entries(idx)%maximum, seconds)
    if (detail_enabled) then
      entries(idx)%step_calls = entries(idx)%step_calls + 1_int64
      entries(idx)%step_total = entries(idx)%step_total + seconds
    endif
  end subroutine end_timer

  subroutine end_timer_loop()
    call flush_timer_detail(step_number)
    step_number = step_number + 1_int64
  end subroutine end_timer_loop

  subroutine flush_timer_detail(step)
! Setup records use -1; hourly step numbering starts at zero and spans years.
    integer(int64), optional, intent(in) :: step
    integer(int64) :: record_step
    integer :: idx
    if (.not.initialized) error stop 'timings: not initialized'
    record_step = -1_int64
    if (present(step)) record_step = step
    if (detail_enabled) then
      do idx = 1, used
        if (entries(idx)%step_calls == 0_int64) cycle
        write(detail_unit,'(I0,A,A,A,I0,A,ES24.16E3)') record_step, ',', &
          csv_name(entries(idx)%name), ',', entries(idx)%step_calls, ',', entries(idx)%step_total
        entries(idx)%step_calls = 0_int64
        entries(idx)%step_total = 0._r8
      enddo
    endif
  end subroutine flush_timer_detail

  subroutine finalize_timer()
    integer :: idx
    integer(int64) :: stop_count
    real(r8) :: run_seconds, percentage
    if (.not.initialized) return
    call system_clock(stop_count)
    if (stop_count < 0_int64) error stop 'timings: wall clock unavailable'
    run_seconds = timer_elapsed(run_stamp, stop_count)
    call flush_timer_detail()
    write(summary_unit,'(A)') '# Wall seconds; nested regions are inclusive and percentages may overlap.'
    write(summary_unit,'(A,ES24.16E3)') '# Whole run seconds: ', run_seconds
    write(summary_unit,'(A)') 'timer,calls,total_s,mean_s,min_s,max_s,run_percent'
    do idx = 1, used
      percentage = 0._r8
      if (run_seconds > 0._r8) percentage = 100._r8 * entries(idx)%total / run_seconds
      write(summary_unit,'(A,A,I0,5(A,ES24.16E3))') csv_name(entries(idx)%name), ',', &
        entries(idx)%calls, ',', entries(idx)%total, ',', entries(idx)%total/real(entries(idx)%calls,r8), &
        ',', entries(idx)%minimum, ',', entries(idx)%maximum, ',', percentage
    enddo
    close(summary_unit)
    if (detail_enabled) close(detail_unit)
    deallocate(entries)
    initialized = .false.
  end subroutine finalize_timer

  function csv_name(name) result(quoted)
    character(len=*), intent(in) :: name
    character(len=:), allocatable :: quoted
    integer :: idx
    quoted = '"'
    do idx = 1, len_trim(name)
      if (name(idx:idx) == '"') quoted = quoted//'"'
      quoted = quoted//name(idx:idx)
    enddo
    quoted = quoted//'"'
  end function csv_name

  function shell_quote(path) result(quoted)
    character(len=*), intent(in) :: path
    character(len=:), allocatable :: quoted
    integer :: idx
    quoted = "'"
    do idx = 1, len(path)
      if (path(idx:idx) == "'") then
        quoted = quoted//"'"//achar(34)//"'"//achar(34)//"'"
      else
        quoted = quoted//path(idx:idx)
      endif
    enddo
    quoted = quoted//"'"
  end function shell_quote
end module timings
