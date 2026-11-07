module EcosimConst
  ! Physical constants and unit conversions shared across EcoSIM processes.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only: GravAcceleration=>DAT_CONST_G  
  implicit none
  character(len=*),private, parameter :: mod_filename = &
    __FILE__

  real(r8), parameter :: secsphour            = 3600._r8                   !seconds per hour,                                                           [s/hour]
  real(r8), parameter :: secspday             = 86400._r8                  !seconds per day,                                                            [s/day]
  real(r8), parameter :: secspyear            = 86400._r8*365._r8          !seconds in a normal earth year
  real(r8), parameter :: cpw                  = 4.186_r8                   !volumetric heat capacity for water,                                         [J/g/K]~[MJ/m3/K]
  real(r8), parameter :: cpi                  = 1.9274_r8                  !volumetric heat capacity for ice,                                           [MJ/m3/K]
  real(r8), parameter :: cpo                  = 2.496_r8                   !volumetric heat capacity for consolidated organic matter,  [MJ m-3 K-1]
  real(r8), parameter :: cpsand               = 2.128_r8                   !volumetric heat capacity for sand and rock,  [MJ m-3 K-1]
  real(r8), parameter :: cpmins               = 2.385_r8                   !volumetric heat capacity for silt and clay, [MJ m-3 K-1]
  real(r8), parameter :: cps                  = 2.095_r8                   !volumetric heat capacity for fresh snow,                                    [MJ/m3/K]
  real(r8), parameter :: TFice                = 273.15_r8                  !frozen temperature,                                                         [K]
  real(r8), parameter :: TC2K                 = 273.15_r8                  !temperature for converting celcius to Kelvin,                               [K]
  real(r8), parameter :: Tref                 = 273.15_r8                  !reference temperature for atmospheric variables,                            [K]
  real(r8), parameter :: VLHeatCapSnoMin      = 2.095E-04_r8               !legacy snow layer/surface-property threshold per unit area, [MJ m-2 K-1]
  real(r8), parameter :: VLHeatCapLitRMin     = 4.190E-05_r8                !minimum heat capacities for solving surface litter water and heat fluxes,  [MJ/K]
  real(r8), parameter :: VLHeatCapSoiMin      = 4.190E-03_r8                !minimum heat capacities for solving soil water and heat fluxes,            [MJ/K]
  real(r8), parameter :: PICON                = 3.14159265358979323846_r8   !Pi,                                                                        [-]
  real(r8), parameter :: PICON2h              = PICON*0.5_r8                !Pi/2,                                                                      [-]
  real(r8), parameter :: TwoPiCON             = PICON*2._r8                 !2Pi,                                                                       [-]
  real(r8), PARAMETER :: PSIPS                = -0.5E-03_r8                 !saturated water pressure,                                                  [Mpa]
  real(r8), parameter :: RadianPerDegree      = PICON/180._r8               !pi/180,                                                                    [rad/degree]
  real(r8), parameter :: LtHeatIceMelt        = 333.55_r8                   !latent heat of fusion release from water to ice,                           [kJ/kg]~[MJ/ton]
  real(r8), parameter :: EvapLHTC             = 2465.0_r8                   !latent heat of vaporization of water,                                      [kJ/kg]
  real(r8), parameter :: SublmHTC             = LtHeatIceMelt+EvapLHTC      !sublimation of ice,                                                        [kJ/kg]
  real(r8), parameter :: TCNDG                = 8.1E-03_r8                  !assumed thermal conductivity below lower soil boundary,                    [MJ m-1 h-1 K-1]
  real(r8), parameter :: RGASC                = 8.3143_r8                   !universal gas constant,                                                    [J (mole K)^-1]
  ! Leaf/surface boundary-layer resistance ratio: r_b(CO2)/r_b(H2O), [-].
  real(r8), parameter :: BndlResistCO2ToH2ORatio = 1.34_r8
  ! Chilton-Colburn analogy: r_b(H2O)/r_b(heat), [-]; r_heat = r_H2O/ratio.
  real(r8), parameter :: BndlResistH2OToHeatRatio = 0.92_r8
  real(r8), parameter :: OMCMassFrac          = 0.55_r8                     !C fraction in organic matter,                                              [gC /g OM] 
  real(r8), parameter :: DensitySolidOM       = 1.30_r8                     !Consolidated dry organic matter density, [g cm-3]~[Mg m-3]
  real(r8), parameter :: DensitySolidMineral  = 2.66_r8                     !mass density of mineral solids, [g cm-3]
  real(r8), parameter :: orgcden              = OMCMassFrac*1.E+06_r8       !volume density of organic carbon,   [cm3 C m-3 OM]
  real(r8), parameter :: hpresc               = 8.4334E+03_r8               !elapsing height for atmospheric pressure,                                  [m]
  real(r8), parameter :: POROQ                = 0.66_r8                     !soil porosity ^ 2/3,                                                       [-]
  real(r8), parameter :: FORGC                = 1.E+05_r8                   !minimum SOC for organic soil,                                              [gC Mg-1 soil]
  real(r8), parameter :: VolMaxSoilMoist4Fire = 1.0_r8                      !maximum soil water content for combustion,                                 [m3 m-3]
  real(r8), parameter :: FrcAsCH4byFire       = 0.01_r8                     !fraction of combusted C released as CH4,                                   [-]
  real(r8), parameter :: PSIHY                = -2500.0_r8                  !hygroscopic water potential, very dry (but not air dry),                   [MPa]
  real(r8), parameter :: OXKM                 = 0.080_r8                    !Km for heterotrophic O2 uptake,                                            [gO m-3]
  real(r8), parameter :: AirFillPore_Min      = 1.0E-03_r8                  !minimum air-filled porosity for gas transfer,                              [m3 m-3]
  real(r8), parameter :: THETPI               = 0.00_r8                     !air content of ice,                                                        [m3 m-3]
  real(r8), parameter :: DENSICE              = 0.92_r8-THETPI              !ice density,                                                               [g/cm3~ton/m3]
  real(r8), parameter :: ZW                   = 0.01_r8                     !snowpack surface roughness,                                                [m]
  real(r8), parameter :: Catomw               = 12._r8                      !C-12 molar mass,                                                           [g/mol]
  real(r8), parameter :: Natomw               = 14._r8                      !N-14 molar mass,                                                           [g/mol]
  real(r8), parameter :: Patomw               = 31._r8                      !P-31 molar mass,                                                           [g/mol]
  real(r8), parameter :: TWILGT               = 0.06976_r8                  !sine of solar inclination angle at twilight,                               [-]
  real(r8), parameter :: gC2MgOM             = 1.E-06_r8/OMCMassFrac        !multiplier to convert organic carbon (gC) to soil organic mass (Mg),       [Mg OM/gC]
  real(r8), parameter :: TKDif               = EvapLHTC/cpw                 !the maximum temperature enhancement due to vapor transport. [K]
  !terminate, [label for variable parsing ]
  real(r8), parameter :: ppmc=1.0E-06_r8                                    !part per million,[-]
  real(r8), parameter :: stefboltz_const=5.670374419e-8_r8*3600e-6          !stefan boltzman constant, [MJ/(hr m^2 K^4)]
  integer , parameter :: ICOR(12)=(/1,-1,0,0,1,1,2,3,3,4,4,5/)  
  real(r8), parameter :: mGravAccelerat=1.e-3_r8*GravAcceleration           !gravitational constant devided by 1000, [km/s2] 

end module EcosimConst
