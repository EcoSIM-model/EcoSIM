MODULE data_kind_mod
   implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
   !----------------------------------------------------------------------------
   ! precision/kind constants add data public
   !----------------------------------------------------------------------------
   public
   integer,parameter :: DAT_KIND_R8 = selected_real_kind(12) ! 8 byte real
   integer,parameter :: DAT_KIND_R4 = selected_real_kind( 6) ! 4 byte real
   integer,parameter :: DAT_KIND_RN = kind(1.0)              ! native real
   integer,parameter :: DAT_KIND_I8 = selected_int_kind (13) ! 8 byte integer
   integer,parameter :: DAT_KIND_I4 = selected_int_kind ( 6) ! 4 byte integer
   integer,parameter :: DAT_KIND_IN = kind(1)                ! native integer
   integer,parameter :: DAT_KIND_CS = 80                     ! short char
   integer,parameter :: DAT_KIND_CL = 256                    ! long char
   integer,parameter :: DAT_KIND_CX = 512                    ! extra-long char
   integer,parameter :: DAT_KIND_CXX= 4096                   ! extra-extra-long char

   type, public :: yearIJ_type
   integer :: year
   integer :: i
   integer :: J
   end type yearIJ_type

end module data_kind_mod 
