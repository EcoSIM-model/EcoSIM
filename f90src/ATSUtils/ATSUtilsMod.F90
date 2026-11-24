module ATSUtilsMod
  use ChemTranspDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoilBGCDataType
  use GridConsts
  use GridDataType
  use SoilWaterDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: ComputeDatefromATS
  ! Add this interface block
  contains

!------------------------------------------------------------------------------------------
  subroutine ComputeDatefromATS(current_day, current_year, current_month, day_of_month, total_days_in_month)
    !This takes the current day (0-364) and current year
    !and returns the current month, day of that month and total days in that month
    implicit none
    integer, intent(in) :: current_day, current_year
    integer, intent(out) :: current_month, day_of_month, total_days_in_month

    integer, parameter :: days_in_month(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    integer :: cumulative_days, i, days_this_month
    logical :: is_leap

    ! Check if leap year
    !is_leap = (mod(current_year, 4) == 0 .and. mod(current_year, 100) /= 0) .or. (mod(current_year, 400) == 0)

    ! Find month and day
    cumulative_days = 0
    do i = 1, 12
      days_this_month = days_in_month(i)
      if (is_leap .and. i == 2) days_this_month = 29

      if (current_day < cumulative_days + days_this_month) then
        current_month = i
        day_of_month = current_day - cumulative_days + 1
        total_days_in_month = days_this_month
        exit
      endif
      cumulative_days = cumulative_days + days_this_month
    enddo

  end subroutine ComputeDatefromATS
end module ATSUtilsMod
