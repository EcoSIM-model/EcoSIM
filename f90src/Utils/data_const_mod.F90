!===============================================================================
! SVN $Id: DAT_const_mod.F90 6749 2007-10-04 20:58:20Z jwolfe $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_131231/shr/DAT_const_mod.F90 $
!===============================================================================

MODULE data_const_mod

   use data_kind_mod, only : DAT_KIND_IN
   use data_kind_mod, only : DAT_KIND_R8
   implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
   integer(DAT_KIND_IN),parameter,private :: R8 = DAT_KIND_R8 ! rename for local readability only

   !----------------------------------------------------------------------------
   ! physical constants (all data public)
   !----------------------------------------------------------------------------
   public
   real(R8),parameter :: DAT_CONST_YEARSECS=86400._R8*365._R8
   real(R8),parameter :: DAT_CONST_PI      = 3.14159265358979323846_R8  ! pi
   real(R8),parameter :: DAT_CONST_CDAY    = 86400.0_R8      ! sec in calendar day ~ sec
   real(R8),parameter :: DAT_CONST_SDAY    = 86164.0_R8      ! sec in siderial day ~ sec
   real(R8),parameter :: DAT_CONST_OMEGA   = 2.0_R8*DAT_CONST_PI/DAT_CONST_SDAY ! earth rot ~ rad/sec
   real(R8),parameter :: DAT_CONST_REARTH  = 6.37122e6_R8    ! radius of earth ~ m
   real(R8),parameter :: DAT_CONST_G       = 9.80616_R8      ! acceleration of gravity ~ m/s^2

   real(R8),parameter :: DAT_CONST_STEBOL  = 5.67e-8_R8      ! Stefan-Boltzmann constant ~ W/m^2/K^4
   real(R8),parameter :: DAT_CONST_BOLTZ   = 1.38065e-23_R8  ! Boltzmann's constant ~ J/K/molecule
   real(R8),parameter :: DAT_CONST_AVOGAD  = 6.02214e26_R8   ! Avogadro's number ~ molecules/kmole
   real(R8),parameter :: DAT_CONST_RGAS    = DAT_CONST_AVOGAD*DAT_CONST_BOLTZ       ! Universal gas constant ~ J/K/kmole
   real(R8),parameter :: DAT_CONST_MWDAIR  = 28.966_R8       ! molecular weight dry air ~ kg/kmole
   real(R8),parameter :: DAT_CONST_MWWV    = 18.016_R8       ! molecular weight water vapor
   real(R8),parameter :: DAT_CONST_RDAIR   = DAT_CONST_RGAS/DAT_CONST_MWDAIR        ! Dry air gas constant     ~ J/K/kg
   real(R8),parameter :: DAT_CONST_RWV     = DAT_CONST_RGAS/DAT_CONST_MWWV          ! Water vapor gas constant ~ J/K/kg
   real(R8),parameter :: DAT_CONST_ZVIR    = (DAT_CONST_RWV/DAT_CONST_RDAIR)-1.0_R8 ! RWV/RDAIR - 1.0
   real(R8),parameter :: DAT_CONST_KARMAN  = 0.4_R8          ! Von Karman constant
   real(R8),parameter :: DAT_CONST_PSTD    = 101325.0_R8     ! standard pressure ~ pascals
   real(R8),parameter :: DAT_CONST_PDB     = 0.0112372_R8    ! ratio of 13C/12C in Pee Dee Belemnite (C isotope standard)

   real(R8),parameter :: DAT_CONST_TKTRIP  = 273.16_R8       ! triple point of fresh water        ~ K
   real(R8),parameter :: DAT_CONST_TKFRZ   = 273.15_R8       ! freezing T of fresh water          ~ K
   real(R8),parameter :: DAT_CONST_TKFRZSW = DAT_CONST_TKFRZ - 1.8_R8 ! freezing T of salt water  ~ K

   real(R8),parameter :: DAT_CONST_RHODAIR = &               ! density of dry air at STP  ~ kg/m^3
                         DAT_CONST_PSTD/(DAT_CONST_RDAIR*DAT_CONST_TKFRZ)
   real(R8),parameter :: DAT_CONST_RHOFW   = 1.000e3_R8      ! density of fresh water     ~ kg/m^3
   real(R8),parameter :: DAT_CONST_RHOSW   = 1.026e3_R8      ! density of sea water       ~ kg/m^3
   real(R8),parameter :: DAT_CONST_RHOICE  = 0.917e3_R8      ! density of ice             ~ kg/m^3
   real(R8),parameter :: DAT_CONST_CPDAIR  = 1.00464e3_R8    ! specific heat of dry air   ~ J/kg/K
   real(R8),parameter :: DAT_CONST_CPWV    = 1.810e3_R8      ! specific heat of water vap ~ J/kg/K
   real(R8),parameter :: DAT_CONST_CPVIR   = (DAT_CONST_CPWV/DAT_CONST_CPDAIR)-1.0_R8 ! CPWV/CPDAIR - 1.0
   real(R8),parameter :: DAT_CONST_CPFW    = 4.188e3_R8      ! specific heat of fresh h2o ~ J/kg/K
   real(R8),parameter :: DAT_CONST_CPSW    = 3.996e3_R8      ! specific heat of sea h2o   ~ J/kg/K
   real(R8),parameter :: DAT_CONST_CPICE   = 2.11727e3_R8    ! specific heat of fresh ice ~ J/kg/K
   real(R8),parameter :: DAT_CONST_LATICE  = 3.337e5_R8      ! latent heat of fusion      ~ J/kg
   real(R8),parameter :: DAT_CONST_LATVAP  = 2.501e6_R8      ! latent heat of evaporation ~ J/kg
   real(R8),parameter :: DAT_CONST_LATSUB  = &               ! latent heat of sublimation ~ J/kg
                         DAT_CONST_LATICE + DAT_CONST_LATVAP
   real(R8),parameter :: DAT_CONST_OCN_REF_SAL = 34.7_R8     ! ocn ref salinity (psu)
   real(R8),parameter :: DAT_CONST_ICE_REF_SAL =  4.0_R8     ! ice ref salinity (psu)

   real(R8),parameter :: DAT_CONST_SPVAL   = 1.0e30_R8       ! special missing value
   integer ,parameter :: DAT_CONST_ISPVAL  = -999            ! special missing value

END MODULE data_const_mod
