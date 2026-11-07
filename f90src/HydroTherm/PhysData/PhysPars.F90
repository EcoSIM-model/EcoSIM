module PhysPars
  use data_kind_mod , only : r8 => DAT_KIND_R8
  use data_const_mod, only : GravAcceleration=>DAT_CONST_G
  use EcoSimConst, only : ppmc  

implicit none

  character(len=*), private, parameter :: mod_filename=&
  __FILE__


! RAM=minimum boundary layer resistance [h m-1]
! MinSnowDepth=minimum snowpack depth for full cover (m)
! RZ=minimum resistance to evaporation of surface water [h m-1]
! RACX,RARX=minimum boundary layer resistances of canopy,litter (h m-1)

! Z1S,Z2SW,Z2SD,Z3SX=parameters for air-water gas transfers in soil
! Z1R,Z2RW,Z2RD,Z3RX=parameters for air-water gas transfers in litter
! Parameters for calculating convective effects on heat transfer
! in porous media (air and water)
! VISCW,VISCA=water,air viscosity (Mg m-1 s)
! CP[w,i,s]=heat capacity of water, ice, snow, kJ/kg/K
! CPo = heat capacity of organic matter, kJ/g/K
! FVOLAH=parameter for clay effect on macropore volume
! DTHETW=difference between saturation and effective saturation

  real(r8), parameter :: TRBA=0.000_r8 !threshold air-filled porosity for convective effects on heat transfer,	[m3 m-3]
  real(r8), parameter :: TRBW=0.375_r8 !threshold water-filled porosity for convective effects on heat transfer,	[m3 m-3]
  real(r8), parameter :: FCI=0.05_r8                      !field capacity of ice
  real(r8), parameter :: WPI=0.025_r8                     !wilting point of ice
  real(r8), parameter :: MinSnowDepth=0.075_r8
  real(r8), parameter :: RAM=1.39E-03_r8   
  real(r8), parameter :: RZ=0.0139_r8  
  real(r8), parameter :: EXPNW=2.07E-04_r8      !parameter used to calculate Nusselt number for water
  real(r8), parameter :: DIFFA=2.01E-05_r8
  real(r8), parameter :: DIFFW=1.45E-07_r8
  real(r8), parameter :: VISCW=ppmc             !water viscosity
  real(r8), parameter :: RYLXW=GravAcceleration*EXPNW/(VISCW*DIFFW)   
  real(r8), parameter :: VISCA=2.0E-08_r8  
  real(r8), parameter :: EXPNA=3.66E-03_r8     !parameter used to calculate Nusselt number for air
  real(r8), parameter :: RYLXA=GravAcceleration*EXPNA/(VISCA*DIFFA)   
  real(r8), parameter :: PRNTW=VISCW/DIFFW
  real(r8), parameter :: PRNTA=VISCA/DIFFA
  real(r8), parameter :: DNUSW=(1.0_r8+(0.492_r8/PRNTW)**0.5625_r8)**0.4444_r8
  real(r8), parameter :: DNUSA=(1.0_r8+(0.492_r8/PRNTA)**0.5625_r8)**0.4444_r8
  real(r8), parameter :: FVOLAH=0.0_r8
  real(r8), parameter :: DTHETW=ppmc
  real(r8), parameter :: ZEROL=1.E-3_r8        !minimum litter fraction
end module PhysPars
