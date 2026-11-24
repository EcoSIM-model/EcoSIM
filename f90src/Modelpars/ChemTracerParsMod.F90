module ChemTracerParsMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  implicit none
  public
  save
  character(len=*),private, parameter :: mod_filename = &
  __FILE__


  real(r8), parameter :: ARSG=6.8E-02_r8   !gaseous AR diffusivity	,[m2 h-1]
  real(r8), parameter :: ARSL=7.7E-06_r8   !aqueous AR diffusivity  ,[m2 h-1]
  real(r8), parameter :: CGSG=4.68E-02_r8  !gaseous CO2 diffusivity	,[m2 h-1]
  real(r8), parameter :: CLSG=4.25E-06_r8  !aqueous CO2 diffusivity,[m2 h-1]
  real(r8), parameter :: CHSG=7.80E-02_r8  !gaseous CH4 diffusivity,[m2 h-1]
  real(r8), parameter :: CQSG=7.08E-06_r8  !aqueous CH4 diffusivity,[m2 h-1]
  real(r8), parameter :: OGSG=6.43E-02_r8  !gaseous O2 diffusivity,[m2 h-1]
  real(r8), parameter :: OLSG=8.57E-06_r8  !aqueous O2 diffusivity,[m2 h-1]
  real(r8), parameter :: ZGSG=5.57E-02_r8  !gaseous N2 diffusivity,[m2 h-1]
  real(r8), parameter :: ZLSG=7.34E-06_r8  !aqueous N2 diffusivity,[m2 h-1]
  real(r8), parameter :: Z2SG=5.57E-02_r8  !gaseous N2O diffusivity,[m2 h-1]
  real(r8), parameter :: ZVSG=5.72E-06_r8  !aqueous N2O diffusivity,[m2 h-1]
  real(r8), parameter :: ZHSG=6.67E-02_r8  !gaseous NH3 diffusivity,[m2 h-1]
  real(r8), parameter :: ZNSG=4.00E-06_r8  !aqueous NH3 diffusivity,[m2 h-1]
  real(r8), parameter :: ZOSG=6.00E-06_r8  !aqueous NO3 diffusivity,[m2 h-1]
  real(r8), parameter :: POSG=3.00E-06_r8  !aqueous PO4 diffusivity,[m2 h-1]
  real(r8), parameter :: OCSG=1.0E-06_r8   !aqueous DOC diffusivity,[m2 h-1]
  real(r8), parameter :: ONSG=1.0E-06_r8   !aqueous DON diffusivity,[m2 h-1]
  real(r8), parameter :: OPSG=1.0E-06_r8   !aqueous DOP diffusivity,[m2 h-1]
  real(r8), parameter :: OASG=3.64E-06_r8  !aqueous acetate diffusivity,[m2 h-1]
  real(r8), parameter :: WGSG=7.70E-02_r8 !water vapor diffusivity,[m2 h-1]
  real(r8), parameter :: ALSG=5.0E-06_r8  !aqueous Al diffusivity,[m2 h-1]
  real(r8), parameter :: FESG=5.0E-06_r8  !aqueous Fe diffusivity,[m2 h-1]
  real(r8), parameter :: HYSG=5.0E-06_r8  !aqueous H diffusivity,[m2 h-1]
  real(r8), parameter :: CASG=5.0E-06_r8  !aqueous Ca diffusivity,[m2 h-1]
  real(r8), parameter :: GMSG=5.0E-06_r8  !aqueous Mg diffusivity,[m2 h-1]
  real(r8), parameter :: ANSG=5.0E-06_r8  !aqueous Na diffusivity,[m2 h-1]
  real(r8), parameter :: AKSG=5.0E-06_r8  !aqueous K diffusivity,[m2 h-1]
  real(r8), parameter :: OHSG=5.0E-06_r8  !aqueous OH diffusivity,[m2 h-1]
  real(r8), parameter :: C3SG=5.0E-06_r8  !aqueous CO3 diffusivity,[m2 h-1]
  real(r8), parameter :: HCSG=5.0E-06_r8  !aqueous HCO3 diffusivity,[m2 h-1]
  real(r8), parameter :: SOSG=5.0E-06_r8  !aqueous SO4 diffusivity,[m2 h-1]
  real(r8), parameter :: CLSX=5.0E-06_r8  !aqueous Cl diffusivity,[m2 h-1]
  real(r8), parameter :: HGSG=5.57E-02_r8 !gaseous H2 diffusivity,[m2 h-1]
  real(r8), parameter :: HLSG=7.34E-06_r8 !aqueous H2 diffusivity,[m2 h-1]

  real(r8), parameter :: SARX =3.421E-02_r8  !Ar solubility coefficient at 25oC, [g solute /g gas]
  real(r8), parameter :: SCO2X=7.391E-01_r8  !CO2 solubility coeficient at 25oC, [g solute /g gas]
  real(r8), parameter :: SCH4X=3.156E-02_r8  !CH4 solubility coeficient at 25oC, [g solute /g gas]
  real(r8), parameter :: SOXYX=2.925E-02_r8  !O2 solubility coeficient at 25oC, [g solute /g gas]
  real(r8), parameter :: SN2GX=1.510E-02_r8  !N2 solubility coeficient at 25oC, [g solute /g gas]
  real(r8), parameter :: SN2OX=5.241E-01_r8  !N2O solubility coeficient  at 25oC, [g solute /g gas]
  real(r8), parameter :: SNH3X=2.852E+02_r8  !NH3 solubility coeficient at 25oC, [g solute /g gas]
  real(r8), parameter :: SH2GX=3.156E-02_r8  !H2 solubility coeficient at 25oC, [g solute /g gas]

  !terminate [label for variable parsing]
end module ChemTracerParsMod
