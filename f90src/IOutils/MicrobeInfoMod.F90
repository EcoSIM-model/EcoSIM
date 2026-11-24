module MicrobeInfoMod

! define microbial parameters
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  use EcoSIMConfig
implicit none
  private
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  
contains

  subroutine WriteMicrobeTraits()
  implicit none
  
  end subroutine WriteMicrobeTraits

end module MicrobeInfoMod