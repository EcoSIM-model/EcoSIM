module EcoSIMHistMod

!!
! data types of plant characteristics
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  CHARACTER(len=16) :: DATAC(30,250,250)
  CHARACTER(len=16),target,allocatable :: DATAP(:,:,:)
  INTEGER, target, allocatable :: datapi(:,:,:)
  CHARACTER(len=16),target,allocatable :: DATAZ(:,:,:)
  CHARACTER(len=16),target,allocatable :: OUTFILS(:,:,:)
  CHARACTER(len=16),target,allocatable :: OUTFILP(:,:,:,:)
  CHARACTER(len=16) :: DATAX_pft(JP),DATAY(JP)
  CHARACTER(len=16) :: OUTS(10)
  CHARACTER(len=16) :: OUTP(10)
  CHARACTER(len=3) :: CHOICE(102,20)
  CHARACTER(len=8) :: CDATE
  CHARACTER(len=640):: PREFIX
  character(len=1280):: outdir
  integer :: IDATA(60),NOUTS(10),NOUTP(10)
  CHARACTER(len=16) :: DATA1(30)               !for model io related setup

  contains

  subroutine InitEcoSIMHistData

  implicit none

  allocate(DATAP(JP,JY,JX))
  allocate(datapi(jp,jy,jx)); datapi=0
  allocate(DATAZ(JP,JY,JX))
  allocate(OUTFILS(10,JY,JX))
  allocate(OUTFILP(10,JP,JY,JX))

  end subroutine InitEcoSIMHistData

!------------------------------------------------------------------------------------------

  subroutine DestructEcoSIMHistData

  use abortutils, only : destroy

  implicit none

  call destroy(datapi)
  call destroy(DATAP)
  call destroy(DATAZ)
  call destroy(OUTFILP)
  call destroy(OUTFILS)

  end subroutine DestructEcoSIMHistData

end module EcoSIMHistMod
