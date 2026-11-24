module DebugToolMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: iulog
  use EcoSIMCtrlMod, only: lverb
  use GridDataType
  use EcoSimConst
  use SoilWaterDataType
implicit none
  private
  character(len=*), parameter :: mod_filename =&
   __FILE__  

  public :: PrintInfo
  public :: DebugPrint
  public :: GetSoilWaterMass
  interface DebugPrint
    module procedure DebugPrint_real_arr
    module procedure DebugPrint_real_sp
    module procedure DebugPrint_real_arrs
    module procedure DebugPrint_int
  end interface DebugPrint
  public :: do_debug
contains

  pure function do_debug()result(ans)
  implicit none
  logical :: ans
  
  ans= lverb
  end function do_debug
  !-----------------------------------------------------------------------
    subroutine PrintInfo(message)
    implicit none
    character(len=*), intent(in) :: message

    if(lverb)write(iulog,*)message
    end subroutine PrintInfo

  !-----------------------------------------------------------------------

   subroutine DebugPrint_real_arr(vnames,vals)
   implicit none
   character(len=32), intent(in) :: vnames(:)
   real(r8), intent(in) :: vals(:)
   integer :: nvars
   character(len=256) :: fmt
   nvars=size(vals)
   
   if(.not.lverb)return

   write(fmt,'(A,I2)')'(A,A,',nvars   
   fmt=trim(fmt)//'E14.6)'   
   write(*,fmt)vnames,'=',vals

   end subroutine DebugPrint_real_arr
  !-----------------------------------------------------------------------

   subroutine DebugPrint_real_arrs(vals)
   implicit none
   real(r8), intent(in) :: vals(:)
   integer :: nvars

   nvars=size(vals)
   
   if(.not.lverb)return
   write(*,*)vals

   end subroutine DebugPrint_real_arrs

  !-----------------------------------------------------------------------

   subroutine DebugPrint_real_sp(vname,val)
   implicit none
   character(len=*), intent(in) :: vname
   real(r8), intent(in) :: val
   
   if(.not.lverb)return
   write(iulog,*)vname,'=',val
   end subroutine DebugPrint_real_sp
  !-----------------------------------------------------------------------

   subroutine DebugPrint_int(vname,val)
   implicit none
   character(len=*), intent(in) :: vname
   integer, intent(in) :: val
   
   if(.not.lverb)return
   write(iulog,*)vname,'=',val
   end subroutine DebugPrint_int

  !-----------------------------------------------------------------------
  function GetSoilWaterMass(NY,NX)result(mass)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8) :: mass
  integer  :: L
  mass = 0._r8
  !add water in soil 
  DO L=NU_col(NY,NX),NL_col(NY,NX)
    mass = mass+ VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)+(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))*DENSICE
  ENDDO
  END function GetSoilWaterMass
end module DebugToolMod