program EcoATSTest
  use ATSCPLMod
  use EcoContainers_module
  use SharedDataMod
  use GridDataType
  use RootDataType
  use EcoSIMCtrlMod
  
  !use ecosim_Time_Mod  , only : ecosim_time_type
  implicit none

  ! Declare variables
  type (EcoEnvironment) :: env
  type (EcoFeedback) :: feedback
  type (EcoSizes) :: sizes
  integer :: NY, NX, L, ii
  integer :: ncells_per_col_, ncol
  !real, dimension(6) :: rain_array
  real, dimension(6) :: rain_array = (/1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3, 1.0e-3/)

  !rain_array = 1.6e-5
  NX = 1
  NYS = 1
  ncells_per_col_ = 100
  ncol = 1

  ! Initialize sizes structure
  sizes%num_components = 1
  sizes%ncells_per_col_ = 100
  sizes%num_columns = 1
  a_bool = .false.
  pheno_bool = .true.
  !Set pft file manually
  pft_file_in = '../../../../input_data/ecosim_pftpar_20261021.nc'
  
  call Init_ATSEcoSIM_driver()

  call Init_EcoSIM(sizes)

  do ii = 1, size(rain_array)
    write(*,*) "For p_rain = ", rain_array(ii)

    do NY=1,NYS
      p_rain(NY) = rain_array(ii)
    enddo

    ! Call Run_EcoSIM_one_step to advance EcoSIM
    call Run_EcoSIM_one_step(sizes)
  end do

end program EcoATSTest

! ----------------------------------------------------------------------------

subroutine Init_ATSEcoSIM_driver()
  use EcoContainers_module
  use SharedDataMod
  use GridDataType
  use CanopyDataType
  implicit none

  ! Declare variables
  type (EcoEnvironment) :: env
  type (EcoFeedback) :: feedback
  type (EcoSizes) :: sizes
  integer :: NY, NX, L, NZ
  integer :: ncells_per_col_, ncol
  integer :: kmo,dofmon,ndaysmon
  real(r8) :: dist_step, dist_tot
  !type(ecosim_time_type) :: etimer

  NX = 1
  NYS = 1
  ncells_per_col_ = 100
  ncol = 1
  dist_step = 1.0
  dist_tot = 1.0 !IC or first layer is zero
  num_pfts = 5

  ! Initialize sizes structure
  sizes%num_components = 1
  sizes%ncells_per_col_ = 100
  sizes%num_columns = 1
  write(*,*) "setting time"
  !
  !call etimer%Init(year0=1804,nyears=12)
  !call etimer%setClock(dtime=3600._r8,nelapstep=0)
  !ndaysmon = etimer%get_curr_mon_days()  ! Would return 31 (days in March)
  !dofmon = etimer%get_curr_dom()         ! Would return 15 (15th day of month)
  !kmo = etimer%get_curr_mon()            ! Would return 3 (March)
  !or set by hand
  ndaysmon = 31
  dofmon = 15
  kmo = 3

  write(*,*) "Initalized timer with Month: ", kmo, " day: ", dofmon, " of ", ndaysmon

  ! needed in starts
  !pressure_at_field_capacity = 0.001
  !pressure_at_wilting_point = 0.001
  !heat_capacity = 7.5e-5

  pressure_at_field_capacity = -3.3e-2
  pressure_at_wilting_point = -1.5
  heat_capacity = 2.0e-2

  !need to allocate these because c_f_pointer sets them in the actual coupler
  !surface vars
  allocate(a_ASP(ncells_per_col_))
  allocate(tairc(1:ncells_per_col_))
  allocate(vpair(1:ncells_per_col_))
  allocate(uwind(1:ncells_per_col_))
  allocate(swrad(1:ncells_per_col_))
  allocate(sunrad(1:ncells_per_col_))
  allocate(p_rain(1:ncells_per_col_))
  allocate(p_snow(1:ncells_per_col_))
  allocate(a_SALB(1:ncells_per_col_))
  allocate(surf_e_source(1:ncells_per_col_))
  allocate(surf_w_source(1:ncells_per_col_))
  allocate(surf_snow_depth(1:ncells_per_col_))
  allocate(a_LAI(1:ncells_per_col_,ncol))
  allocate(a_SAI(1:ncells_per_col_,ncol))
  allocate(column_area(ncol))

  !subsurface vars
  allocate(a_TEMP(ncells_per_col_, ncol))
  allocate(a_CumDepz2LayBottom_vr(ncells_per_col_, ncol))
  allocate(a_AREA3(ncells_per_col_, ncol))
  allocate(a_BKDSI(ncells_per_col_, ncol))
  allocate(a_WC(ncells_per_col_, ncol))
  allocate(a_MATP(ncells_per_col_, ncol))
  allocate(a_PORO(ncells_per_col_, ncol))
  allocate(a_LDENS(ncells_per_col_, ncol))
  allocate(a_LSAT(ncells_per_col_, ncol))
  allocate(a_SSWS(ncells_per_col_,ncol))
  allocate(a_rdens(ncells_per_col_,ncol))
  
  !Pheno vars
  !allocate(LAI_col(ncells_per_col_,ncol))
  allocate(a_VEG(num_pfts))

  do NY=1,NYS
    do L=1,ncells_per_col_
      a_AREA3(L,NY) = 1.0
      !DH_col(NY,NX) = 0.316229
      !DV(NY,NX) = 0.316229
      a_CumDepz2LayBottom_vr(L,NY) = dist_tot
      dist_tot = dist_tot+dist_step
      a_BKDSI(L,NY) = 1.1*1.0e3_r8
      !a_CORGC(L,NY) = 0.0
      !a_CORGN(L,NY) = 0.0
      !a_CORGP(L,NY) = 0.0
      a_WC(L,NY) = 6.0e-4
      a_TEMP(L,NY) = 242.00
      a_MATP(L,NY) = -6.9
      a_PORO(L,NY) = 0.5
      a_LDENS(L,NY) = 1.0
      a_LSAT(L,NY) = 0.5
      a_rdens(L,NY) = 0.1
    enddo
    column_area(NY) = 1.0
    a_ASP(NY) = 0.0
    a_SALB(NY) = 0.5
    tairc(NY) = 242.13003959655759
    !double counting the conversions I think
    !vpair(NY) = 736.3/1.0e6_r8
    !uwind(NY) = 1.0*3600.0_r8
    !swrad(NY) = 400.0*0.0036_r8
    !sunrad(NY) =  219.78*0.0036_r8
    p_rain(NY) = 0.0
    p_snow(NY) = 0.0
    !p_rain(NY) = 3.e-8*1000.0_r8*3600.0_r8

    vpair(NY) = 3.9167352020740509E-002*1.0e3
    uwind(NY) = 1.1

    a_LAI(NY,1) = 0.2
    a_SAI(NY,1) = 0.05
    do NZ=1,num_pfts
      a_VEG(NZ) = 1
    enddo

  enddo

end subroutine Init_ATSEcoSIM_driver
