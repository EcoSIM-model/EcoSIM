module MicrobeConfigMod
  use MicBGCPars
  use abortutils, only : endrun,iulog
  use fileUtil
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: ReadMicrobeNamelist
  contains

  subroutine ReadMicrobeNamelist(nml_buffer)
  implicit none
  character(len=ecosim_namelist_buffer_size), intent(in) :: nml_buffer
  character(len=*), parameter :: subname='ReadMicrobeNamelist'

  namelist /Microbes/NumGuild_Heter_Aerob_Bact,NumGuild_Heter_Aerob_Fung, &
    NumGuild_Heter_Facul_Dent,NumGuild_Heter_Aerob_N2Fixer,               &
    NumGuild_Heter_Anaer_N2Fixer,NumGuild_Heter_Anaer_Fermentor,          &
    NumGuild_Heter_AcetoMethanogen,NumGuild_Autor_H2genMethanogen,        &
    NumGuild_Autor_AmoniaOxidBact,NumGuild_Autor_NitritOxidBact,          &
    NumGuild_Autor_AerobMethOxid

  integer :: nml_error
  character(len=256) :: ioerror_msg

  !set all functional group having one guild by default

  NumGuild_Heter_Aerob_Bact      = 1
  NumGuild_Heter_Aerob_Fung      = 1
  NumGuild_Heter_Facul_Dent      = 1
  NumGuild_Heter_Aerob_N2Fixer   = 1
  NumGuild_Heter_Anaer_N2Fixer   = 1
  NumGuild_Heter_Anaer_Fermentor = 1
  NumGuild_Heter_AcetoMethanogen = 1
  NumGuild_Autor_H2genMethanogen = 1
  NumGuild_Autor_AmoniaOxidBact  = 1
  NumGuild_Autor_NitritOxidBact  = 1
  NumGuild_Autor_AerobMethOxid   = 1

  read(nml_buffer, nml=Microbes, iostat=nml_error, iomsg=ioerror_msg)
  if (nml_error /= 0) then
    write(iulog,'(a)')"ERROR reading ecosim namelist ",nml_error,ioerror_msg
    call endrun('stopped in '//trim(mod_filename), __LINE__)
  end if

  end subroutine ReadMicrobeNamelist

end module MicrobeConfigMod