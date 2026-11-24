module GridMod

implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: SetMesh, SetMeshATS
  public :: get_col, get_pft
contains

!------------------------------------------------------------------------

  subroutine SetMesh(NHW,NVN,NHE,NVS)

  use EcoSIMConfig, only : column_mode
  use EcoSIMCtrlMod,only : grid_file_in,lverb,first_topou
  use ncdio_pio
  use netcdf
  use abortutils, only : endrun,iulog
  use GridConsts, only : JX,JY,JZ,JH,JV,JD,bounds,JP,JX0,JY0
!  set up the landscape rectangular mesh
!  beginning(NHW,NVN)
!  o--------------------------x
!  |                          |
!  |                          |
!  |                          |
!  |                          |
!  x--------------------------o
!                             end (NHE,NVS)

  implicit none
  integer, intent(out) :: NHW   !upper corner x index
  integer, intent(out) :: NVN   !upper corner y index
  integer, intent(out) :: NHE   !lower corner x index
  integer, intent(out) :: NVS   !lower corner y index
  integer :: nextra_grid
  INTEGER :: NZ,NY,NX,ic,ip
  type(file_desc_t) :: grid_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar

  call ncd_pio_openfile(grid_nfid, grid_file_in, ncd_nowrite)
  write(*,*)'open file ok',grid_file_in
  bounds%ngrid=get_dim_len(grid_nfid, 'ngrid')

  bounds%ntopou=get_dim_len(grid_nfid, 'ntopou')

  write(*,*)bounds%ngrid,' grids'
  call check_var(grid_nfid, 'NHW', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NHW in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NHW), 'in '//trim(mod_filename))

  call check_var(grid_nfid, 'NHE', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NHE in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NHE), 'in '//trim(mod_filename))

  call check_var(grid_nfid, 'NVN', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NVN in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NVN), 'in '//trim(mod_filename))

  call check_var(grid_nfid, 'NVS', vardesc, readvar)
  if(.not. readvar)then
    call endrun('fail to find NVS in '//trim(mod_filename), __LINE__)
  endif
  call check_ret(nf90_get_var(grid_nfid%fh, vardesc%varid, NVS), 'in '//trim(mod_filename))

  call ncd_pio_closefile(grid_nfid)
  if(lverb)write(*,*)'read grid finished'

  bounds%NHW =NHW
  bounds%NVN =NVN
  bounds%NHE =NHE
  bounds%NVS =NVS
  if(first_topou)then  
    NHW=1;NVN=1;NHE=1;NVS=1
    bounds%ntopou=1
    bounds%NHW =1
    bounds%NVN =1
    bounds%NHE =1
    bounds%NVS =1
  endif  
  JZ=20
  bounds%begg=1;bounds%endg=bounds%ngrid
  bounds%begt=1;bounds%endt=bounds%ntopou
  nextra_grid=1
  JX=(NHE-NHW)+1;JX0=JX
  JY=(NVS-NVN)+1;JY0=JY     

  bounds%ncols=JX*JY
  bounds%npfts=bounds%ncols*JP
  bounds%begc=1;bounds%endc=bounds%ncols
  bounds%begp=1;bounds%endp=bounds%npfts

  allocate(bounds%icol(JY,JX))
  allocate(bounds%ipft(JP,JY,JX))

  if(column_mode)nextra_grid=0
  JX=JX+nextra_grid
  JY=JY+nextra_grid
  
  ic=0;ip=0
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      ic=ic+1
      bounds%icol(NY,NX)=ic
      DO NZ=1,JP
        ip=ip+1
        bounds%ipft(NZ,NY,NX)=ip
      ENDDO
    ENDDO
  ENDDO
  
  JZ=20
  !read JZ from input data!
  JH=JX+nextra_grid
  JV=JY+nextra_grid
  JD=JZ+1

  end subroutine SetMesh

!------------------------------------------------------------------------

  !Version of the SetMesh function that does not attempt to load a
  !netcdf file as the mesh information is provided by ATS

  subroutine SetMeshATS(NHW,NVN,NHE,NVS)

  use EcoSIMConfig, only : column_mode
  use EcoSIMCtrlMod,only : grid_file_in
  use ncdio_pio
  use netcdf
  use abortutils, only : endrun,iulog
  use GridConsts, only : JX,JY,JZ,JH,JV,JD,bounds,JP,JX0,JY0
!  set up the landscape rectangular mesh
!  beginning(NHW,NVN)
!  o--------------------------x
!  |                          |
!  |                          |
!  |                          |
!  |                          |
!  x--------------------------o
!                             end (NHE,NVS)

  implicit none
  integer, intent(in) :: NHW   !upper corner x index
  integer, intent(in) :: NVN   !upper corner y index
  integer, intent(in) :: NHE   !lower corner x index
  integer, intent(in) :: NVS   !lower corner y index
  integer :: nextra_grid
  INTEGER :: NZ,NY,NX,ic,ip
  type(file_desc_t) :: grid_nfid
  type(Var_desc_t) :: vardesc
  logical :: readvar

  write(*,*) "(SetMeshATS) JX,JY,JZ: ", JX, JY, JZ

  bounds%NHW =NHW
  bounds%NVN =NVN
  bounds%NHE =NHE
  bounds%NVS =NVS

  bounds%begg=1;bounds%endg=bounds%ngrid
  bounds%begt=1;bounds%endt=bounds%ntopou
  
  nextra_grid=1
  JX=(NHE-NHW)+1;JX0=JX
  JY=(NVS-NVN)+1;JY0=JY 

  bounds%ncols=JX*JY
  bounds%npfts=bounds%ncols*JP
  bounds%begc=1;bounds%endc=bounds%ncols
  bounds%begp=1;bounds%endp=bounds%npfts

  if (associated(bounds%icol)) deallocate(bounds%icol)
  if (associated(bounds%ipft)) deallocate(bounds%ipft)
  allocate(bounds%icol(JY,JX))
  allocate(bounds%ipft(JP,JY,JX))

  if(column_mode)nextra_grid=0
  JX=JX+nextra_grid
  JY=JY+nextra_grid

  ic=0;ip=0
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      ic=ic+1
      bounds%icol(NY,NX)=ic
      DO NZ=1,JP
        ip=ip+1
        bounds%ipft(NZ,NY,NX)=ip
      ENDDO
    ENDDO
  ENDDO
  !read JZ from input data!

  JH=JX+nextra_grid
  JV=JY+nextra_grid
  JD=JZ+1

  write(*,*) "(SetMeshATS f) JX,JY,JZ: ", JX, JY, JZ

  end subroutine SetMeshATS

!------------------------------------------------------------------------

  integer function get_col(NY,NX)

  use GridConsts, only : JX0,JY0
  implicit none
  integer, intent(in) :: NY,NX

  get_col=(NX-1)*JY0+NY
  end function get_col
!------------------------------------------------------------------------
  integer function get_pft(NZ,NY,NX)

  use GridConsts, only : JX0,JY0,JP
  implicit none
  integer, intent(in) :: NZ,NY,NX

  get_pft=(NX-1)*(JP*JY0)+(NY-1)*JP+NZ

  end function get_pft

end module GridMod
