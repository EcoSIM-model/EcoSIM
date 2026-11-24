module EcoSIMSolverPar

! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public

  real(r8) :: dts_wat        !time step for water transport calculation, [h]
  real(r8) :: dts_sno        !time step for snow process, [h]
  real(r8) :: XNPB
  real(r8) :: dt_watvap      !time step for water vapor fluxes, [h]
  real(r8) :: XNPD
  real(r8) :: dts_HeatWatTP   !time step for heat, water, solute transport , [h]
  real(r8) :: dt_GasCyc       !time step for gas iteration/cycling, [h]
  real(r8) :: dts_gas         !time step for gas flx update, [h]
  real(r8) :: XNPR
  real(r8) :: XNPS
  integer :: NPH             !number of model cycles per hour for heat and water fluxes
  integer :: NPT             !number of model cycles per hour for gas fluxes within each NPH
  integer :: NPG             !number of model cycles per hour for gas fluxes
  integer :: NPR             !number of cycles per time step used to calculate hat and water transfer in surface litter
  integer :: NPS             !number of cycles for solving snowpack heat and water fluxes
  real(r8) :: oscal_test=0._r8
end module EcoSIMSolverPar
