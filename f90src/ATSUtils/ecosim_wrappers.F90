!
!There needs to be a wrapper for the eocsim f90 driver
!as there are differences between how gfortran and intel compilers
!handle mangling conventions.
!
! Copied from the alquimia wrapper:
! **************************************************************************** !
!
! PFloTran Alquimia Inteface Wrappers
!
! Author: Benjamin Andre
!
! Different fortran compilers use different name mangling conventions
! for fortran modules:
!
!    gfortran : ___modulename_MOD_procedurename
!
!    intel : _modulename_mp_procedurename_
!
!    as a consequence we can't put the alquimia interface into a
!    module and call it directly from C/C++. Instead we use
!    some simple wrapper functions.
!
! Notes:
!
!  * Function call signatures are dictated by the alquimia API!
!
!  * alquimia data structures defined in AlquimiaContainers_module
!    (alquimia_containers.F90) are dictated by the alquimia API.
!
! **************************************************************************** !

subroutine EcoSIM_DataTest() bind(c)
    use, intrinsic :: iso_c_binding
    !Simple test for wrapper functionality
    implicit none

    write(*,*) "in the data test"

  end subroutine EcoSIM_DataTest

! **************************************************************************** !

subroutine EcoSIM_Setup(config, env, feedback, istate, sizes, num_iterations,&
                        num_columns, ncells_per_col_) bind(C)

  use, intrinsic :: iso_c_binding

  use EcoContainers_module
  use ATSCPLMod, only : ATS2EcoSIMData, Init_EcoSIM
  use ATSStateRegistryMod, only : PackInternalState

  implicit none

  ! function parameters
  type (EcoConfig), intent(in) :: config
  type (EcoEnvironment), intent(in) :: env
  type (EcoFeedback), intent(in) :: feedback
  type (EcoInternalState), intent(inout) :: istate
  type (EcoSizes), intent(in) :: sizes
  integer(c_int), VALUE :: num_columns
  integer(c_int), VALUE :: num_iterations
  integer(c_int), VALUE :: ncells_per_col_

  call ATS2EcoSIMData(env, feedback, sizes, config)

  call Init_EcoSIM(sizes)

  ! give ATS the initialized carried state; outputs are computed by the first advance
  call PackInternalState(istate, sizes, .false.)

end subroutine EcoSIM_Setup

! **************************************************************************** !

subroutine EcoSIM_Shutdown() bind(C)

  !For now this does nothing, but it should clear all
  !the data structures
  use, intrinsic :: iso_c_binding

  implicit none

end subroutine EcoSIM_Shutdown

! **************************************************************************** !

subroutine EcoSIM_Advance( &
     delta_t, &
     env, &
     feedback, &
     istate, &
     sizes, &
     num_iterations, &
     num_columns) bind(C)

  use, intrinsic :: iso_c_binding
  use EcoContainers_module
  use ATSCPLMod, only : Run_EcoSIM_one_step, ATS2EcoSIMData, EcoSIM2ATSData
  use ATSStateRegistryMod, only : PackInternalState, UnpackInternalState

  implicit none

  ! function parameters
  real (c_double), value, intent(in) :: delta_t
  type (EcoEnvironment), intent(in) :: env
  type (EcoFeedback), intent(in) :: feedback
  type (EcoInternalState), intent(inout) :: istate
  type (EcoSizes), intent(in) :: sizes
  integer(c_int), value, intent(in) :: num_iterations
  integer(c_int), value, intent(in) :: num_columns

  call ATS2EcoSIMData(env, feedback, sizes)

  ! restore carried state from ATS (identical unless ATS restarted)
  call UnpackInternalState(istate, sizes)

  call Run_EcoSIM_one_step(sizes)

  call EcoSIM2ATSData(feedback, sizes)

  call PackInternalState(istate, sizes, .true.)

end subroutine EcoSIM_Advance

! **************************************************************************** !
!
! Layout of the EcoSIM internal state (see ATSStateRegistryMod). Static: it
! depends only on sizes, so ATS can query it before EcoSIM_Setup.
!
! **************************************************************************** !

function EcoSIM_Internal_State_Layout_Version() &
     bind(C, name="ecosim_internal_state_layout_version") result(version)

  use, intrinsic :: iso_c_binding
  use ATSStateRegistryMod, only : kLayoutVersion

  implicit none
  integer(c_int) :: version

  version = kLayoutVersion

end function EcoSIM_Internal_State_Layout_Version

! **************************************************************************** !

function EcoSIM_Internal_State_Num_Entries(sizes) &
     bind(C, name="ecosim_internal_state_num_entries") result(num_entries)

  use, intrinsic :: iso_c_binding
  use EcoContainers_module, only : EcoSizes
  use ATSStateRegistryMod, only : NumStateEntries

  implicit none
  type (EcoSizes), intent(in) :: sizes
  integer(c_int) :: num_entries

  num_entries = NumStateEntries(sizes)

end function EcoSIM_Internal_State_Num_Entries

! **************************************************************************** !

subroutine EcoSIM_Internal_State_Entry(i, sizes, ats_name, ecosim_name, units, &
     description, ncomp, role) bind(C, name="ecosim_internal_state_entry")

  ! i is 0-based. The character buffers must hold kNameLen, kNameLen,
  ! kUnitsLen and kDescLen characters; they are returned null-terminated.
  use, intrinsic :: iso_c_binding
  use EcoContainers_module, only : EcoSizes
  use ATSStateRegistryMod, only : GetStateEntry, kNameLen, kUnitsLen, kDescLen

  implicit none
  integer(c_int), value, intent(in) :: i
  type (EcoSizes), intent(in) :: sizes
  character(kind=c_char), intent(out) :: ats_name(kNameLen), ecosim_name(kNameLen)
  character(kind=c_char), intent(out) :: units(kUnitsLen), description(kDescLen)
  integer(c_int), intent(out) :: ncomp, role

  character(len=kNameLen) :: f_ats_name, f_ecosim_name
  character(len=kUnitsLen) :: f_units
  character(len=kDescLen) :: f_description
  integer :: f_ncomp, f_role

  call GetStateEntry(i+1, sizes, f_ats_name, f_ecosim_name, f_units, f_description, &
       f_ncomp, f_role)
  call to_c_string(f_ats_name, ats_name)
  call to_c_string(f_ecosim_name, ecosim_name)
  call to_c_string(f_units, units)
  call to_c_string(f_description, description)
  ncomp = f_ncomp
  role = f_role

contains

  subroutine to_c_string(fstr, cstr)
    character(len=*), intent(in) :: fstr
    character(kind=c_char), intent(out) :: cstr(:)
    integer :: k, n
    n = min(len_trim(fstr), size(cstr)-1)
    do k = 1, n
      cstr(k) = fstr(k:k)
    enddo
    cstr(n+1:) = c_null_char
  end subroutine to_c_string

end subroutine EcoSIM_Internal_State_Entry

! **************************************************************************** !
!
! Sizes in bytes of the exchange containers as EcoSIM sees them, so ATS can
! check them against its C structs (EcoEngine::CheckContainerSizes). Order:
! EcoVectorDouble, EcoMatrixDouble, EcoTensorDouble, EcoSizes, EcoConfig,
! EcoEnvironment, EcoFeedback, EcoInternalState. On input n is the number of
! entries ATS can take; on output, the number EcoSIM knows.
!

subroutine EcoSIM_Container_Sizes(n, sizes) bind(C, name="ecosim_container_sizes")

  use, intrinsic :: iso_c_binding
  use EcoContainers_module

  implicit none
  integer(c_int), intent(inout) :: n
  integer(c_size_t), intent(out) :: sizes(*)

  integer, parameter :: num_types = 8
  integer(c_size_t) :: s(num_types)
  type (EcoVectorDouble) :: vector
  type (EcoMatrixDouble) :: matrix
  type (EcoTensorDouble) :: tensor
  type (EcoSizes) :: esizes
  type (EcoConfig) :: config
  type (EcoEnvironment) :: env
  type (EcoFeedback) :: feedback
  type (EcoInternalState) :: istate

  s = (/ c_sizeof(vector), c_sizeof(matrix), c_sizeof(tensor), c_sizeof(esizes), &
         c_sizeof(config), c_sizeof(env), c_sizeof(feedback), c_sizeof(istate) /)
  sizes(1:min(n, num_types)) = s(1:min(n, num_types))
  n = num_types

end subroutine EcoSIM_Container_Sizes
