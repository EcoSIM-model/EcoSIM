module EcoSiMParDataMod
 use PlantBGCPars, only : plant_bgc_par_type
 use MicBGCPars, only : MicParType
implicit none
  character(len=*),private, parameter :: mod_filename =&
   __FILE__
  type(plant_bgc_par_type), target, public :: pltpar
  type(MicParType), target, public :: micpar
end module EcoSiMParDataMod
