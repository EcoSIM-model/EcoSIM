module ATSCPLMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SharedDataMod
  use ATSEcoSIMInitMod
  use ATSEcoSIMAdvanceMod
  use EcoContainers_module
  use c_f_interface_module, only : c_f_string_ptr
  implicit none

  public
  character(len=*), private, parameter :: mod_filename=&
  __FILE__

contains
!------------------------------------------------------------------------------------------

  subroutine ATS2EcoSIMData(env, feedback, sizes, config)
  ! ATS -> EcoSIM. Called at setup (with config) and at every advance
  ! (without config, which is copied only once).
  implicit none
  type (EcoEnvironment), intent(in) :: env
  type (EcoFeedback), intent(in) :: feedback
  type (EcoSizes), intent(in) :: sizes
  type (EcoConfig), intent(in), optional :: config

  if (present(config)) call CopyConfigFromATS(config)
  call CopyEnvironmentFromATS(env, sizes)
  call CopyFeedbackFromATS(feedback, sizes)

  end subroutine ATS2EcoSIMData
!------------------------------------------------------------------------------------------

  subroutine EcoSIM2ATSData(feedback, sizes)
  ! EcoSIM -> ATS, after every advance
  implicit none
  type (EcoFeedback), intent(in) :: feedback
  type (EcoSizes), intent(in) :: sizes

  call CopyFeedbackToATS(feedback, sizes)

  end subroutine EcoSIM2ATSData
!------------------------------------------------------------------------------------------

  subroutine CopyConfigFromATS(config)
  ! run parameters and flags, once at setup; kept in SharedDataMod variables
  implicit none
  type (EcoConfig), intent(in) :: config
  character(len=512) :: f_string_buffer

  call c_f_string_ptr(config%pft_file, f_string_buffer)
  ecosim_pft_file_path = trim(f_string_buffer)

  heat_capacity = config%heat_capacity
  pressure_at_field_capacity = config%field_capacity
  pressure_at_wilting_point = config%wilting_point
  site_latitude = config%latitude
  site_solar_noon = config%solar_noon
  p_bool = config%p_bool
  a_bool = config%a_bool
  pheno_bool = config%pheno_bool
  wind_speed_ref_height = config%wind_speed_ref_height
  soil_surface_roughness = config%soil_surface_roughness

  end subroutine CopyConfigFromATS
!------------------------------------------------------------------------------------------

  subroutine CopyEnvironmentFromATS(env, sizes)
  ! soil state and properties, geometry, forcing, vegetation, clock;
  ! at setup and every advance
  implicit none
  type (EcoEnvironment), intent(in) :: env
  type (EcoSizes), intent(in) :: sizes

  real(r8), pointer :: data(:)
  real(r8), pointer :: data2D(:,:)
  real(r8), pointer :: data3D(:,:,:)
  integer :: size_col, num_cols, num_components
  type(c_ptr) :: data_ptr

  size_col = sizes%ncells_per_col_
  num_cols = env%shortwave_radiation%size
  num_components = env%mole_fraction%components
  num_pfts = sizes%num_pfts

  data_ptr = env%mole_fraction%data
  call c_f_pointer(data_ptr, data3D, [size_col, num_cols, num_components])
  !a_MFrac=data3D(:,:,:)

  data_ptr = env%temperature%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_TEMP=data2D(:,:)

  data_ptr = env%thermal_conductivity%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_TCOND=data2D(:,:)

  data_ptr = env%depth%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_CumDepz2LayBottom_vr = data2D(:,:)

  data_ptr = env%dz%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_dz = data2D(:,:)

  data_ptr = env%volume%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_Volume = data2D(:,:)

  data_ptr = env%water_content%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_WC = data2D(:,:)

  call c_f_pointer(env%column_area%data, data, (/num_cols/))
  column_area = data(:)

  data_ptr = env%bulk_density%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_BKDSI = data2D(:,:)

  data_ptr = env%liquid_density%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_LDENS = data2D(:,:)

  data_ptr = env%rock_density%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_RDENS = data2D(:,:)

  data_ptr = env%matric_pressure%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_MATP = data2D(:,:)

  data_ptr = env%porosity%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_PORO = data2D(:,:)

  data_ptr = env%liquid_saturation%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_LSAT = data2D(:,:)

  data_ptr = env%hydraulic_conductivity%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_HCOND = data2D(:,:)

  data_ptr = env%rooting_depth_fraction%data
  call c_f_pointer(data_ptr, data2D, [size_col, num_cols])
  a_FC = data2D(:,:)

  call c_f_pointer(env%shortwave_radiation%data, data, (/num_cols/))
  swrad = data(:)

  call c_f_pointer(env%longwave_radiation%data, data, (/num_cols/))
  sunrad = data(:)

  call c_f_pointer(env%air_temperature%data, data, (/num_cols/))
  tairc = data(:)

  call c_f_pointer(env%vapor_pressure_air%data, data, (/num_cols/))
  vpair = data(:)

  call c_f_pointer(env%wind_speed%data, data, (/num_cols/))
  uwind = data(:)

  call c_f_pointer(env%aspect%data, data, (/num_cols/))
  a_ASP = data(:)

  !call c_f_pointer(env%LAI%data, data, (/num_cols/))
  !a_LAI = data(:)

  call c_f_pointer(env%LAI%data, data2D, [size_col,num_cols])
  a_LAI = data2D(:,:)

  !call c_f_pointer(env%SAI%data, data, (/num_cols/))
  !a_SAI = data(:)

  call c_f_pointer(env%SAI%data, data2D, [size_col,num_cols])
  a_SAI = data2D(:,:)

  call c_f_pointer(env%vegetation_type%data, data, (/num_cols/))
  a_VEG = data(:)

  call c_f_pointer(env%snow_albedo%data, data, (/num_cols/))
  a_SALB = data(:)

  call c_f_pointer(env%surface_temperature%data, data, (/num_cols/))
  a_TSURF = data(:)

  ! PFT index of each column, in the first num_pfts entries of the column
  call c_f_pointer(env%plant_functional_type%data, data2D, [size_col,num_cols])
  a_PFT = data2D(:,:)

  atm_n2 = env%atm_n2
  atm_o2 = env%atm_o2
  atm_co2 = env%atm_co2
  atm_ch4 = env%atm_ch4
  atm_n2o = env%atm_n2o
  atm_h2 = env%atm_h2
  atm_nh3 = env%atm_nh3
  current_day = env%current_day
  current_year = env%current_year
  current_hour = env%current_hour

  ! p_bool was set from the config at setup
  if(p_bool)THEN
    call c_f_pointer(env%precipitation%data, data, (/num_cols/))
    p_total = data(:)
  else
    call c_f_pointer(env%precipitation%data, data, (/num_cols/))
    p_rain = data(:)

    call c_f_pointer(env%precipitation_snow%data, data, (/num_cols/))
    p_snow = data(:)
  endif

  end subroutine CopyEnvironmentFromATS
!------------------------------------------------------------------------------------------

  subroutine CopyFeedbackFromATS(feedback, sizes)
  ! snow depth is EcoSIM's snow state, kept by ATS between advances. The
  ! sources are recomputed by EcoSIM before use; they are copied so that the
  ! SharedDataMod arrays are allocated (allocation on assignment).
  implicit none
  type (EcoFeedback), intent(in) :: feedback
  type (EcoSizes), intent(in) :: sizes

  real(r8), pointer :: data(:)
  real(r8), pointer :: data2D(:,:)
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = feedback%snow_depth%size

  call c_f_pointer(feedback%subsurface_water_source%data, data2D, [size_col, num_cols])
  a_SSWS = data2D(:,:)

  call c_f_pointer(feedback%subsurface_energy_source%data, data2D, [size_col, num_cols])
  a_SSES = data2D(:,:)

  call c_f_pointer(feedback%surface_water_source%data, data, (/num_cols/))
  surf_w_source = data(:)

  call c_f_pointer(feedback%surface_energy_source%data, data, (/num_cols/))
  surf_e_source = data(:)

  call c_f_pointer(feedback%snow_depth%data, data, (/num_cols/))
  surf_snow_depth = data(:)

  call c_f_pointer(feedback%canopy_snow%data, data2D, [sizes%num_pfts, num_cols])
  a_CanSnow = data2D(:,:)

  end subroutine CopyFeedbackFromATS
!------------------------------------------------------------------------------------------

  subroutine CopyFeedbackToATS(feedback, sizes)
  ! sources and snow depth computed by this advance. The container is
  ! intent(in): only the arrays it points to are written.
  implicit none
  type (EcoFeedback), intent(in) :: feedback
  type (EcoSizes), intent(in) :: sizes

  real(r8), pointer :: data(:)
  real(r8), pointer :: data2D(:,:)
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = sizes%num_columns

  call c_f_pointer(feedback%subsurface_water_source%data, data2D, [size_col, num_cols])
  data2D(:,:) = a_SSWS

  call c_f_pointer(feedback%subsurface_energy_source%data, data2D, [size_col, num_cols])
  data2D(:,:) = a_SSES

  call c_f_pointer(feedback%surface_water_source%data, data, (/num_cols/))
  data(:) = surf_w_source

  call c_f_pointer(feedback%surface_energy_source%data, data, (/num_cols/))
  data(:) = surf_e_source

  call c_f_pointer(feedback%snow_depth%data, data, (/num_cols/))
  data(:) = surf_snow_depth

  call c_f_pointer(feedback%canopy_snow%data, data2D, [sizes%num_pfts, num_cols])
  data2D(:,:) = a_CanSnow

  end subroutine CopyFeedbackToATS

!------------------------------------------------------------------------------------------

  subroutine Run_EcoSIM_one_step(sizes)
  implicit none

  type (EcoSizes), intent(in) :: sizes

  !copy data from compuler to EcoSIM

  !run surface energy balance
  call SurfaceEBalance(sizes)

  !copy data back to coupler
  end subroutine Run_EcoSIM_one_step
!------------------------------------------------------------------------------------------

  subroutine Init_EcoSIM(sizes)
  !initialize ecosim
  implicit none

  type (EcoSizes), intent(in) :: sizes
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = sizes%num_columns

  call InitSharedData(size_col,num_cols)

  call Init_EcoSIM_Soil(num_cols)
  end subroutine Init_EcoSIM
!------------------------------------------------------------------------------------------

  subroutine SurfaceEBalance(sizes)
  implicit none

  integer :: K, vec_size
  type (EcoSizes), intent(in) :: sizes
  integer :: size_col, num_cols

  size_col = sizes%ncells_per_col_
  num_cols = sizes%num_columns

  !vec_size = size(surf_e_source)
  call RunEcoSIMSurfaceBalance(num_cols)

  end subroutine SurfaceEBalance

!------------------------------------------------------------------------------------------

  subroutine SetEcoSizes(sizes)

    use EcoContainers_module, only : EcoSizes

    implicit none

    type (EcoSizes), intent(out) :: sizes

    sizes%num_components = 1
    sizes%ncells_per_col_ = 100
    sizes%num_columns = 1

    write(*,*) "(SetEcoSizes f): N_cells, N_cols ", sizes%ncells_per_col_, sizes%num_columns
  end subroutine SetEcoSizes

!-----------------------------------------------------------------------------------------
end module ATSCPLMod
