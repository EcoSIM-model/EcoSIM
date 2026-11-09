
! This code is adapted from Alquimia with the credits below:
! Alquimia Copyright (c) 2013-2016, The Regents of the University of California,
! through Lawrence Berkeley National Laboratory (subject to receipt of any
! required approvals from the U.S. Dept. of Energy).  All rights reserved.
!
! Alquimia is available under a BSD license. See LICENSE.txt for more
! information.
!
! If you have questions about your rights to use or distribute this software,
! please contact Berkeley Lab's Technology Transfer and Intellectual Property
! Management at TTD@lbl.gov referring to Alquimia (LBNL Ref. 2013-119).
!
! NOTICE.  This software was developed under funding from the U.S. Department
! of Energy.  As such, the U.S. Government has been granted for itself and
! others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
! license in the Software to reproduce, prepare derivative works, and perform
! publicly and display publicly.  Beginning five (5) years after the date
! permission to assert copyright is obtained from the U.S. Department of Energy,
! and subject to any subsequent five (5) year renewals, the U.S. Government is
! granted for itself and others acting on its behalf a paid-up, nonexclusive,
! irrevocable, worldwide license in the Software to reproduce, prepare derivative
! works, distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Authors: Benjamin Andre <bandre@lbl.gov>
!        : Andrew Graus <agraus@lbl.gov>
!

! **************************************************************************** !
!
! ATS-EcoSIM Containers module
!
! Author: Andrew Graus
!
! WARNINGS:
!
!   * The data structures defined here are dictated by
!     the API! Do NOT change them unless you make
!     corresponding changes to the c containers (and doc).
!
!   * The order of the data members matters! If num_primary is the
!     first member and num_minerals is the third, then they must be in
!     those positions on both the c and fortran side of the interface!
!
!   * The names of the containers and their members are the same
!     for both c and fortran. The language interface doesn't require
!     this, but it makes it easier for the reader to understand what is
!     going on.
!
!   * The C side is data/EcoContainers.hh in the ATS EcoSIM PK, which
!     documents the role of each container. ecosim_container_sizes()
!     (ecosim_wrappers.F90) reports c_sizeof of these types so ATS can
!     check them at setup.
!
! **************************************************************************** !

module EcoContainers_module

  use, intrinsic :: iso_c_binding

  implicit none

  ! --------------------------------------------------------------------------
  ! primitive containers; matrices are cells x columns, column-major
  ! --------------------------------------------------------------------------
  type, public, bind(c) :: EcoVectorDouble
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type EcoVectorDouble

  type, public, bind(c) :: EcoVectorInt
     integer (c_int) :: size
     integer (c_int) :: capacity
     type (c_ptr) :: data
  end type EcoVectorInt

  type, public, bind(c) :: EcoMatrixDouble
     integer (c_int) :: cells
     integer (c_int) :: columns
     integer (c_int) :: capacity_cells
     integer (c_int) :: capacity_columns
     type (c_ptr) :: data
  end type EcoMatrixDouble

  type, public, bind(c) :: EcoMatrixInt
     integer (c_int) :: cells
     integer (c_int) :: columns
     integer (c_int) :: capacity_cells
     integer (c_int) :: capacity_columns
     type (c_ptr) :: data
  end type EcoMatrixInt

  type, public, bind(c) :: EcoTensorDouble
     integer (c_int) :: cells
     integer (c_int) :: columns
     integer (c_int) :: components
     integer (c_int) :: capacity_cells
     integer (c_int) :: capacity_columns
     integer (c_int) :: capacity_components
     type (c_ptr) :: data
  end type EcoTensorDouble

  type, public, bind(c) :: EcoTensorInt
     integer (c_int) :: cells
     integer (c_int) :: columns
     integer (c_int) :: components
     integer (c_int) :: capacity_cells
     integer (c_int) :: capacity_columns
     integer (c_int) :: capacity_components
     type (c_ptr) :: data
  end type EcoTensorInt

  type, public, bind(c) :: EcoSizes
     integer (c_int) :: ncells_per_col_
     integer (c_int) :: num_components
     integer (c_int) :: num_columns
     integer (c_int) :: num_pfts
  end type EcoSizes

  ! --------------------------------------------------------------------------
  ! ATS -> EcoSIM, once at setup: run parameters and flags
  ! --------------------------------------------------------------------------
  type, public, bind(c) :: EcoConfig
     real (c_double) :: heat_capacity
     real (c_double) :: field_capacity
     real (c_double) :: wilting_point
     real (c_double) :: latitude
     real (c_double) :: solar_noon
     logical (c_bool) :: p_bool
     logical (c_bool) :: a_bool
     logical (c_bool) :: pheno_bool
     logical (c_bool) :: microbe_bool
     real (c_double) :: wind_speed_ref_height
     real (c_double) :: soil_surface_roughness
     type (c_ptr) :: pft_file
  end type EcoConfig

  ! --------------------------------------------------------------------------
  ! ATS -> EcoSIM, at setup and every advance; never read back
  ! (members not filled by ATS are marked in data/EcoContainers.hh)
  ! --------------------------------------------------------------------------
  type, public, bind(c) :: EcoEnvironment
     ! per cell: ncells_per_col_ x num_columns
     type (EcoMatrixDouble) :: liquid_density
     type (EcoMatrixDouble) :: gas_density
     type (EcoMatrixDouble) :: ice_density
     type (EcoMatrixDouble) :: rock_density
     type (EcoMatrixDouble) :: porosity
     type (EcoMatrixDouble) :: water_content
     type (EcoMatrixDouble) :: matric_pressure
     type (EcoMatrixDouble) :: temperature
     type (EcoMatrixDouble) :: hydraulic_conductivity
     type (EcoMatrixDouble) :: bulk_density
     type (EcoMatrixDouble) :: liquid_saturation
     type (EcoMatrixDouble) :: gas_saturation
     type (EcoMatrixDouble) :: ice_saturation
     type (EcoMatrixDouble) :: relative_permeability
     type (EcoMatrixDouble) :: thermal_conductivity
     type (EcoMatrixDouble) :: volume
     type (EcoMatrixDouble) :: depth
     type (EcoMatrixDouble) :: dz
     type (EcoMatrixDouble) :: plant_wilting_factor
     type (EcoMatrixDouble) :: rooting_depth_fraction
     type (EcoMatrixDouble) :: plant_functional_type
     type (EcoMatrixDouble) :: LAI
     type (EcoMatrixDouble) :: SAI
     type (EcoTensorDouble) :: mole_fraction
     ! per column: num_columns
     type (EcoVectorDouble) :: column_area
     type (EcoVectorDouble) :: shortwave_radiation
     type (EcoVectorDouble) :: longwave_radiation
     type (EcoVectorDouble) :: air_temperature
     type (EcoVectorDouble) :: vapor_pressure_air
     type (EcoVectorDouble) :: wind_speed
     type (EcoVectorDouble) :: precipitation
     type (EcoVectorDouble) :: precipitation_snow
     type (EcoVectorDouble) :: elevation
     type (EcoVectorDouble) :: aspect
     type (EcoVectorDouble) :: slope
     type (EcoVectorDouble) :: vegetation_type
     type (EcoVectorDouble) :: snow_albedo
     type (EcoVectorDouble) :: surface_temperature
     ! atmosphere composition
     real (c_double) :: atm_n2
     real (c_double) :: atm_o2
     real (c_double) :: atm_co2
     real (c_double) :: atm_ch4
     real (c_double) :: atm_n2o
     real (c_double) :: atm_h2
     real (c_double) :: atm_nh3
     ! clock
     integer (c_int) :: current_day
     integer (c_int) :: current_year
     integer (c_int) :: current_hour
  end type EcoEnvironment

  ! --------------------------------------------------------------------------
  ! EcoSIM -> ATS, every advance (snow_depth and canopy_snow are also sent in)
  ! --------------------------------------------------------------------------
  type, public, bind(c) :: EcoFeedback
     type (EcoMatrixDouble) :: subsurface_water_source
     type (EcoMatrixDouble) :: subsurface_energy_source
     type (EcoVectorDouble) :: surface_water_source
     type (EcoVectorDouble) :: surface_energy_source
     type (EcoVectorDouble) :: snow_depth
     type (EcoMatrixDouble) :: canopy_snow
  end type EcoFeedback

  ! --------------------------------------------------------------------------
  ! EcoSIM-private carried state and outputs, stored by ATS (ATSStateRegistryMod)
  ! --------------------------------------------------------------------------
  type, public, bind(c) :: EcoInternalState
     integer (c_int) :: layout_version
     integer (c_int) :: num_entries
     integer (c_int) :: num_columns
     integer (c_int) :: values_per_column
     type (EcoMatrixDouble) :: values   ! values_per_column x num_columns
  end type EcoInternalState

end module EcoContainers_module
