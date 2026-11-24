module PlantTestMod
 use ecosim_log_mod    , only : errMsg => shr_log_errMsg
 use bhistMod           , only : hist_freq_str_len
 use abortutils        , only : endrun
 use fileUtil          , only : stdout
 use ModelStatusType   , only : model_status_type
 use data_kind_mod     , only : r8 => DAT_KIND_R8
implicit none
  character(len=*), private,parameter :: mod_filename = &
  __FILE__

  public :: getvarllen
  public :: getvarlist
  public :: initmodel
  public :: runplant
  contains

! ----------------------------------------------------------------------
  function getvarllen()result(ans)
  implicit none

  integer :: ans

  ! 5 fluxes vars and 5 mass vars
  ans = 10
  end function getvarllen

! ----------------------------------------------------------------------

  subroutine getvarlist(nvars, varl, varlnml, unitl, vartypes)

  use bhistMod, only : hist_var_str_len,hist_unit_str_len, hist_var_lon_str_len
  use fileUtil, only :  var_flux_type, var_state_type
  implicit none
  integer, intent(in) :: nvars
  character(len=hist_var_str_len), intent(out) :: varl(:)     !variable name
  character(len=hist_var_lon_str_len), intent(out) :: varlnml(:)     !variable name
  character(len=hist_unit_str_len),intent(out) :: unitl(:)
  integer                         ,intent(out) :: vartypes(:)

  varl(1)='var1_con'; varlnml(1)='var1 concentration'; unitl(1)='mol m-3'
  varl(2)='var2_con'; varlnml(2)='var2 concentration'; unitl(2)='mol m-3'
  varl(3)='var3_con'; varlnml(3)='var3 concentration'; unitl(3)='mol m-3'
  varl(4)='var4_con'; varlnml(4)='var4 concentration'; unitl(4)='mol m-3'
  varl(5)='var5_con'; varlnml(5)='var5 concentration'; unitl(5)='mol m-3'
  varl(6)='flx1_con'; varlnml(6)='flx1'; unitl(6)='mol m-2 s-1'
  varl(7)='flx2_con'; varlnml(7)='flx1'; unitl(7)='mol m-2 s-1'
  varl(8)='flx3_con'; varlnml(8)='flx1'; unitl(8)='mol m-2 s-1'
  varl(9)='flx4_con'; varlnml(9)='flx1'; unitl(9)='mol m-2 s-1'
  varl(10)='flx5_con';varlnml(10)='flx1'; unitl(10)='mol m-2 s-1'

  end subroutine getvarlist


! ----------------------------------------------------------------------

  subroutine runplant(nvars,ystates0l, ystatesfl, err_status)

  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(in) :: ystates0l(nvars)
  real(r8), intent(out) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status

  integer :: jj

  call err_status%reset()

  do jj = 1, 5
    ystatesfl(jj)=ystates0l(jj)
    ystatesfl(jj+5)=ystates0l(jj+5)
  enddo


  end subroutine runplant


! ----------------------------------------------------------------------

  subroutine initmodel(nvars, ystatesfl, err_status)

  implicit none
  integer, intent(in) :: nvars
  real(r8), intent(inout) :: ystatesfl(nvars)
  type(model_status_type), intent(out) :: err_status

  call err_status%reset()
  ystatesfl(1:5)=0.01_r8
  ystatesfl(6:10)=1.e-5_r8

  end subroutine initmodel

! ----------------------------------------------------------------------

end module PlantTestMod
