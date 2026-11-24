module FireMod

  use data_kind_mod, only: r8 => DAT_KIND_R8
  use EcoSIMCtrlMod, only: soil_mgmt_in
  use abortutils,    only: endrun
implicit none

  private
  save
  character(len=*), parameter :: mod_filename = &
  __FILE__

  logical, public :: use_fire =.false.
  integer, parameter :: max_fire_events=20
  integer :: fire_years(max_fire_events)
  integer :: num_fire_events=0
  integer :: fire_event_loc=1
  character(len=21) :: fire_files(max_fire_events)

  public :: config_fire
  public :: check_fire
  contains

  subroutine config_fire(FireEvents)
  !
  !configure fire set up
  !FireEvents should be of format: 'year1/file1;year2/file2'
  !where year1 is a four character string specifies year
  !and file1 specifies the fire specifics
  implicit none
  character(len=*), intent(in) :: FireEvents
  integer :: pair_start,pair_end,colon_pos,ios
  character(len=128) :: msg
  character(len=4) :: year_str
  character(len=30):: fire_event

  !initialization
  use_fire=.false.
  num_fire_events =0
  fire_event_loc=1
  if(len_trim(FireEvents)<5)return
  if(FireEvents(5:5)/='/')return
   
  pair_start = 1  
  ! Split by semicolons and process each pair
  do while (pair_start <= len_trim(FireEvents))
    if(num_fire_events>=max_fire_events)then
      call endrun('Too many fire events in '//trim(mod_filename),__LINE__)
    endif
    ! Find next semicolon or end of string
    pair_end = index(FireEvents(pair_start:), ';')
        
    if (pair_end == 0) then
      ! Last pair
      fire_event = FireEvents(pair_start:)
    else
      ! Extract pair
      pair_end = pair_start + pair_end - 2
      fire_event = FireEvents(pair_start:pair_end)
    end if
        
    ! Process this pair (year:file)
    colon_pos = index(fire_event, '/')
    if (colon_pos == 0 .or. colon_pos <= 4) then
      write(msg,*) 'Error: Invalid specification format for fire event: ', trim(fire_event)
      call endrun(trim(msg)//' in file '//trim(mod_filename),__LINE__)      
    end if
        
    ! Extract and validate year
    year_str = fire_event(colon_pos-4:colon_pos-1)
    read(year_str, '(I4)', iostat=ios) fire_years(num_fire_events + 1)
    if (ios /= 0) then
      write(msg,*) 'Error: Invalid year in pair: ', trim(fire_event)
      call endrun(trim(msg)//' in file '//trim(mod_filename),__LINE__)                  
    end if
    if(num_fire_events>0)then
      if(fire_years(num_fire_events+1)<=fire_years(num_fire_events))then
        call endrun('Fire event years must be strictly increasing in '//trim(mod_filename),__LINE__)
      endif
    endif
        
    ! Extract filename
    fire_files(num_fire_events + 1) = trim(fire_event(colon_pos+1:))
        
    ! Check if filename is empty
    if (len_trim(fire_files(num_fire_events + 1)) == 0) then
      write(msg,*) 'Error: Empty filename in pair: ', trim(fire_event)
      call endrun(trim(msg)//' in file '//trim(mod_filename),__LINE__)     
    end if
        
    num_fire_events = num_fire_events + 1
        
    ! Move to next pair
    if (pair_end == 0) exit

    pair_start = pair_end + 2  ! Skip semicolon

  end do  
  !the model uses fire
  use_fire=.true.
  if(use_fire)then
   if(soil_mgmt_in=='NO')then
     call endrun('fire mode is on, but soil_mgmt_in is not defined. Stop in file '//trim(mod_filename),__LINE__)     
   endif
  endif
  !set current fire event to 1
  fire_event_loc=1

  end subroutine config_fire
!------------------------------------------------------------------------------------------
  function check_fire(year_model,fire_event_entry)result(fire_yesno)
  implicit none
  integer, intent(in) :: year_model  !current model year
  character(len=21), intent(out) :: fire_event_entry  !retrieved fire_event
  logical :: fire_yesno

  fire_yesno=.false.
  fire_event_entry=''
  ! Skip events before the current run/restart year. Keep indexed accesses
  ! inside the bounds test; Fortran does not guarantee short-circuit .and.
  do while(fire_event_loc<=num_fire_events)
    if(fire_years(fire_event_loc)>=year_model)exit
    fire_event_loc=fire_event_loc+1
  enddo
  if(fire_event_loc>num_fire_events)return
  fire_yesno=(fire_years(fire_event_loc)==year_model)

  !move to the next fire event
  if(fire_yesno)then
    fire_event_entry=fire_files(fire_event_loc)  
    fire_event_loc=fire_event_loc+1
  endif  

  end function check_fire
end module FireMod
