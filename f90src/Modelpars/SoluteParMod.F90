module SoluteParMod
!     EQUILIBRIUM CONSTANTS
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  implicit none
  public
  character(len=*),private, parameter :: mod_filename = &
  __FILE__

  real(r8), parameter :: RUreaInhibtorConst(0:2)=real((/10.0E-02,1.0E-02,0.5E-02/),r8) !rate constants for decline in urea hydrolysis inhibition, [h^-1]
  real(r8), PARAMETER :: DPH2O=6.5E-09_r8                !equilibrium constant for H2O=H(+)+OH(-), [mol^2 m^-6] 
  real(r8), PARAMETER :: SPALO=6.5E-22_r8                !equilibrium constant for AlOH3(s)=Al(3+)+3OH(-),[mol^3 m^-9]
  real(r8), PARAMETER :: SPFEO=6.5E-27_r8                !equilibrium constant for FeOH3(s)=Fe(3+)+3OH(-),[mol^3 m^-9]
  real(r8), PARAMETER :: SPCAC=3.8E-03_r8                !equilibrium constant for CaCO3(s)=Ca(2+)+CO3(2-),[mol m^-3]
  real(r8), PARAMETER :: SPCAS=1.4E+01_r8                !equilibrium constant for CaSO4(s)=Ca(2+)+SO4(2-),[mol m^-3]
  real(r8), PARAMETER :: SPALP=1.0E-15_r8                !equilibrium constant for AlPO4(s)=Al(3+)+PO4(3-),[mol m^-3]
  real(r8), PARAMETER :: SPFEP=1.0E-20_r8                !equilibrium constant for FePO4(s)=Fe(3+)+PO4(3-),[mol m^-3]
  real(r8), PARAMETER :: SPCAM=7.0E+07_r8                !equilibrium constant for Ca(H2PO4)2(s)=Ca(2+)+2H2PO4(-),[mol^3 m^-9]
  real(r8), PARAMETER :: SPCAD=1.0E-01_r8                !equilibrium constant for CaHPO4(s)=Ca(2+)+HPO4(2-),[mol^3 m^-9]
  real(r8), PARAMETER :: SPCAH=2.3E-31_r8                !equilibrium constant for hydroxyapatite,Ca5(PO4)3OH(s)=5Ca(2+)+3PO4(3-)+OH(-),[mol^8 m^-24]
  real(r8), PARAMETER :: SXOH2=4.5E-05_r8                !equilibrium constant for X-OH2(+)=X-OH+H(+),[mol m^-3]
  real(r8), PARAMETER :: SXOH1=1.1E-06_r8                !equilibrium constant for X-OH1=X-O(-)+H(+),[mol m^-3]
  real(r8), PARAMETER :: SXH2P=2.0E+07_r8                !equilibrium constant for X-H2PO4+OH(-)=X-OH+H2PO4(-), dimensionless
  real(r8), PARAMETER :: SXH1P=2.0E+07_r8                !equilibrium constant for X-H1PO4(-)=X-OH+HPO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPCO2=4.2E-04_r8                !equilibrium constant for CO2 + H2O = HCO3(-) + H(+), [mol m^-3]
  real(r8), PARAMETER :: DPHCO=5.6E-08_r8                !equilibrium constant for HCO3(-)  = CO3(2-) + H(+), [mol m^-3]
  real(r8), PARAMETER :: DPN4=5.7E-07_r8                 !equilibrium constant for NH4(+) = NH3 + H(+), [mol m^-3]
  real(r8), PARAMETER :: DPAL1=4.6E-07_r8                !equilibrium constant for AlOH(2+)=Al(3+)+OH(-),[mol m^-3]
  real(r8), PARAMETER :: DPAL2=7.3E-07_r8                !equilibrium constant for Al(OH)2(+)=AlOH(2+)+OH(-), [mol m^-3]
  real(r8), PARAMETER :: DPAL3=1.8E-05_r8                !equilibrium constant for Al(OH)3(aq)=Al(OH)2(+)+OH(-),[mol m^-3]
  real(r8), PARAMETER :: DPAL4=1.2E-05_r8                !equilibrium constant for Al(OH)4(-)=Al(OH)3(aq)+OH(-),[mol m^-3]
  real(r8), PARAMETER :: DPALS=0.16_r8                   !equilibrium constant for Al(SO4)(+)=Al(3+)+SO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPFE1=4.6E-07_r8                !equilibrium constant for Fe(OH)(2+)=Fe(3+)+OH(-),[mol m^-3]
  real(r8), PARAMETER :: DPFE2=7.3E-07_r8                !equilibrium constant for Fe(OH)2(+)=FeOH(2+)+OH(-),[mol m^-3]
  real(r8), PARAMETER :: DPFE3=1.8E-05_r8                !equilibrium constant for Fe(OH)3(aq)=Fe(OH)2(+)+OH(-),[mol m^-3]
  real(r8), PARAMETER :: DPFE4=1.2E-05_r8                !equilibrium constant for Fe(OH)4(-)=Fe(OH)3(aq)+OH(-),[mol m^-^3]
  real(r8), PARAMETER :: DPFES=7.1E-02_r8                !equilibrium constant for Fe(SO4)(+)=Fe(3+)+SO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPCAO=12.5_r8                   !equilibrium constant for CaOH(+)=Ca(2+)+OH(-),[mol m^-3]
  real(r8), PARAMETER :: DPCAC=4.2E-02_r8                !equilibrium constant for CaCO3(aq)=Ca(2+)+CO3(2-),[mol m^-3]
  real(r8), PARAMETER :: DPCAH=13.5_r8                   !equilibrium constant for CaHCO3(+)=Ca(2+)+HCO3(-),[mol m^-3]
  real(r8), PARAMETER :: DPCAS=1.2_r8                    !equilibrium constant for CaSO4=Ca(2+)+SO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPMGO=0.7_r8                    !equilibrium constant for MgOH(+)=Mg(2+)+OH(-),[mol m^-3]
  real(r8), PARAMETER :: DPMGC=0.3_r8                    !equilibrium constant for MgCO3=Mg(2+)+CO3(2-),[mol m^-3]
  real(r8), PARAMETER :: DPMGH=67.0_r8                   !equilibrium constant for MgHCO3(-)=Mg(2+)+HCO3(-),[mol m^-3]
  real(r8), PARAMETER :: DPMGS=2.1_r8                    !equilibrium constant for MgSO4=Mg(2+)+SO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPNAC=0.45_r8                   !equilibrium constant for NaCO3(-)=Na(+)+CO3(2-),[mol m^-3]
  real(r8), PARAMETER :: DPNAS=3.3E+02_r8                !equilibrium constant for NaSO4(-)=Na(+)+SO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPKAS=5.0E+01_r8                !equilibrium constant for KSO4(-)=K(+)+SO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPH1P=4.5E-10_r8                !equilibrium constant for HPO4(2-)=H(+)+PO4(3-),[mol m^-3]
  real(r8), PARAMETER :: DPH2P=6.3E-05_r8                !equilibrium constant for H2PO4(-)=H(+)+HPO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPH3P=7.1_r8                    !equilibrium constant for H3PO4=H(+)+H2PO4(-),[mol m^-3]
  real(r8), PARAMETER :: DPF1P=4.5E-02_r8                !equilibrium constant for FeHPO4(+)=Fe(3+)+HPO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPF2P=3.7E-03_r8                !equilibrium constant for FeH2PO4(2+)=Fe(3+)+H2PO4(-),[mol m^-3]
  real(r8), PARAMETER :: DPC0P=3.5E-04_r8                !equilibrium constant for CaPO4(-)=Ca(2+)+PO4(3-),[mol m^-3]
  real(r8), PARAMETER :: DPC1P=1.82_r8                   !equilibrium constant for CaHPO4=Ca(2+)+HPO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPC2P=40.0_r8                   !equilibrium constant for CaH2PO4(+)=Ca(2+)+H2PO4(-),[mol m^-3]
  real(r8), PARAMETER :: DPM1P=1.23_r8                   !equilibrium constant for MgHPO4=Mg(2+)+HPO4(2-),[mol m^-3]
  real(r8), PARAMETER :: DPCOH=1.0E-02_r8                !equilibrium constant for X-COOH=X-COO(-)+H(+),[mol m^-3]
  real(r8), PARAMETER :: DPALO=6.3E+04_r8                !equilibrium constant for X-COO-Al(OH)2+H(+)=Al(OH)2(+)+X-COOH,[mol m^-3]
  real(r8), PARAMETER :: DPFEO=6.3E+04_r8                !equilibrium constant for X-COO-Fe(OH)2+H(+)=Fe(OH)2(+)+X-COOH,[mol m^-3]

  real(r8), PARAMETER :: DPCO3=DPCO2*DPHCO               !equilibrium constant for H2CO3= CO3(2-)+2H(+), [mol^2 m^-6] 
  real(r8), PARAMETER :: SHALO=SPALO/DPH2O**3            !equilibrium constant for Al(OH)3(s)+3H(+) =Al(3+)+3H2O, [m^6 mol^-2]
  real(r8), PARAMETER :: SYAL1=SPALO/DPAL1               !equilibrium constant for Al(OH)3(s) = Al(OH)(2+)+2OH(-), [mol^2 m^-6]
  real(r8), PARAMETER :: SHAL1=SYAL1/DPH2O**2            !equilibrium constant for Al(OH)3(s) + 2H(+)=Al(3+)+OH(-)+H2O, [m^3 mol^-1]
  real(r8), PARAMETER :: SYAL2=SYAL1/DPAL2               !equilibrium constant for Al(OH)3(s) = Al(OH)2(+)+OH(-), [mol m^-3]
  real(r8), PARAMETER :: SHAL2=SYAL2/DPH2O               !equilibrium constant for Al(OH)3(s)+H(+) = Al(OH)2(+) + H2O, [m^3 mol^-1]
  real(r8), PARAMETER :: SPAL3=SYAL2/DPAL3               !equilibrium constant for Al(OH)3(s) = Al(OH)3(aq),[-]
  real(r8), PARAMETER :: SYAL4=SPAL3/DPAL4               !equilibrium constant for Al(OH)3(s)+OH(-) = Al(OH)4(-), [m^3 mol^-1]
  real(r8), PARAMETER :: SHAL4=SYAL4*DPH2O               !equilibrium constant for Al(OH)3(s)+H2O = Al(OH)4(-)+H(+), [mol m^-3]
  real(r8), PARAMETER :: SHFEO=SPFEO/DPH2O**3            !equilibrium constant for Fe(OH)3(s)+3H(+) = Fe(3+)+3H2O, [m^6 mol^-2]
  real(r8), PARAMETER :: SYFE1=SPFEO/DPFE1               !equilibrium constant for Fe(OH)3(s) = Fe(OH)(2+)+2OH(-), [mol^2 m^-6]
  real(r8), PARAMETER :: SHFE1=SYFE1/DPH2O**2            !equilibrium constant for Fe(OH)3(s) + 2H(+)=Fe(3+)+OH(-)+H2O, [m^3 mol^-1]
  real(r8), PARAMETER :: SYFE2=SYFE1/DPFE2               !equilibrium constant for Fe(OH)3(s) = Fe(OH)2(+)+OH(-), [mol m^-3]
  real(r8), PARAMETER :: SHFE2=SYFE2/DPH2O               !equilibrium constant for Fe(OH)3(s)+H(+) = Fe(OH)2(+) + H2O, [m^3 mol^-1]
  real(r8), PARAMETER :: SPFE3=SYFE2/DPFE3               !equilibrium constant for Fe(OH)3(s) = Fe(OH)3(aq),[-]
  real(r8), PARAMETER :: SYFE4=SPFE3/DPFE4               !equilibrium constant for Fe(OH)3(s)+OH(-) = Fe(OH)4(-), [m^3 mol^-1]
  real(r8), PARAMETER :: SHFE4=SYFE4*DPH2O               !equilibrium constant for Fe(OH)3(s)+H2O = Fe(OH)4(-)+H(+), [mol m^-3]
  real(r8), PARAMETER :: SHCAC1=SPCAC/DPHCO              !equilibrium constant for CaCO3(s)+H(+) = Ca(2+) + HCO3(-), [-]
  real(r8), PARAMETER :: SYCAC1=SHCAC1*DPH2O             !equilibrium constant for CaCO3(s)+H2O = Ca(2+)+OH(-)+HCO3(-), [mol^2 m^-6]
  real(r8), PARAMETER :: SHCAC2=SHCAC1/DPCO2             !equilibrium constant for CaCO3(s)+2H(+)=Ca(2+)+CO2(aq)+H2O, [m^3 mol^-1]
  real(r8), PARAMETER :: SYCAC2=SHCAC2*DPH2O**2          !equilibrium constant for CaCO3(s)+H2O = Ca(2+)+2OH(-)+CO2(aq),[mol^3 m^-9]
  real(r8), PARAMETER :: SHA0P1=SPALP/DPH1P              !equilibrium constant for AlPO4(s)+H(+) = Al(3+)+HPO4(2-), [-]
  real(r8), PARAMETER :: SYA0P1=SHA0P1*DPH2O             !equilibrium constant for AlPO4(s)+H2O =Al(3+)+HPO4(2-)+OH(-), [mol^2 m^-6]
  real(r8), PARAMETER :: SPA1P1=SYA0P1/DPAL1             !equilibrium constant for AlPO4(s)+H2O =Al(OH)(2+)+HPO4(2-), [mol m^-3]
  real(r8), PARAMETER :: SYA2P1=SPA1P1/DPAL2             !equilibrium constant for AlPO4(s)+H2O+OH(-)=HPO4(2-)+Al(OH)2(+), [-]
  real(r8), PARAMETER :: SHA2P1=SYA2P1*DPH2O             !equilibrium constant for AlPO4(s)+2H2O = HPO4(2-)+Al(OH)2(+)+H(+),[mol^2 m^-6]
  real(r8), PARAMETER :: SYA3P1=SYA2P1/DPAL3             !equilibrium constant for AlPO4(s)+H2O+2OH(-)=Al(OH)3(aq)+HPO4(2-),[m^3 mol^-1]
  real(r8), PARAMETER :: SHA3P1=SYA3P1*DPH2O**2          !equilibrium constant for AlPO4(s)+3H2O=Al(OH)3(aq)+HPO4(2-)+2H(+),[mol^3 m^-9]
  real(r8), PARAMETER :: SYA4P1=SYA3P1/DPAL4             !equilibrium constant for AlPO4(s)+H2O+3OH(-)=Al(OH)4(-)+HPO4(2-),[m^6 mol^-2]
  real(r8), PARAMETER :: SHA4P1=SYA4P1*DPH2O**3          !equilibrium constant for AlPO4(s)+4H2O=Al(OH)4(-)+HPO4(2-)+3H(+),[mol^4 m^-12]
  real(r8), PARAMETER :: SHA0P2=SHA0P1/DPH2P             !equilibrium constant for AlPO4(s)+2H(+)=Al(3+)+H2PO4(-), [m^3 mol^-1]
  real(r8), PARAMETER :: SYA0P2=SHA0P2*DPH2O**2          !equilibrium constant for AlPO4(s)+2H2O=Al(3+)+H2PO4(-)+2OH(-),[mol^3 m^-9]
  real(r8), PARAMETER :: SYA1P2=SYA0P2/DPAL1             !equilibrium constant for AlPO4(s)+2H2O=Al(OH)(2+)+H2PO4(-)+OH(-),[mol^2 m^-6]
  real(r8), PARAMETER :: SHA1P2=SYA1P2/DPH2O             !equilibrium constant for AlPO4(s)+H2O+H(+)=Al(OH)(2+)+H2PO4(-),[-]
  real(r8), PARAMETER :: SPA2P2=SYA1P2/DPAL2             !equilibrium constant for AlPO4(s)+2H2O=Al(OH)2(+)+H2PO4(-),[mol m^-3]
  real(r8), PARAMETER :: SYA3P2=SPA2P2/DPAL3             !equilibrium constant for AlPO4(s)+2H2O+OH(-)=Al(OH)3(aq)+H2PO4(-),[-]
  real(r8), PARAMETER :: SHA3P2=SYA3P2*DPH2O             !equilibrium constant for AlPO4(s)+3H2O=Al(OH)3(aq)+H2PO4(-)+H(+),[mol^2 m^-6]
  real(r8), PARAMETER :: SYA4P2=SYA3P2/DPAL4             !equilibrium constant for AlPO4(s)+2H2O+2OH(-)=Al(OH)4(-)+H2PO4(-),[m^3 mol^-1]
  real(r8), PARAMETER :: SHA4P2=SYA4P2*DPH2O**2          !equilibrium constant for AlPO4(s)+4H2O=Al(OH)4(-)+H2PO4(-)+2H(+),[mol^3 m^-9]

  real(r8), PARAMETER :: SHF0P1=SPFEP/DPH1P              !equilibrium constant for FePO4(s)+H(+)= Fe(3+)+HPO4(2-), [-]
  real(r8), PARAMETER :: SYF0P1=SHF0P1*DPH2O             !equilibrium constant for FePO4(s)+H2O = Fe(3+)+HPO4(2-)+OH(-), [mol^2 m^-6]
  real(r8), PARAMETER :: SPF1P1=SYF0P1/DPFE1             !equilibrium constant for FePO4(s)+H2O = Fe(OH)(2+)+HPO4(2-), [mol m^-3]
  real(r8), PARAMETER :: SYF2P1=SPF1P1/DPFE2             !equilibrium constant for FePO4(s)+H2O+OH(-)=HPO4(2-)+Fe(OH)2(+), [-]
  real(r8), PARAMETER :: SHF2P1=SYF2P1*DPH2O             !equilibrium constant for FePO4(s)+2H2O = HPO4(2-)+Fe(OH)2(+)+H(+),[mol^2 m^-6]
  real(r8), PARAMETER :: SYF3P1=SYF2P1/DPFE3             !equilibrium constant for FePO4(s)+H2O+2OH(-)=Fe(OH)3(aq)+HPO4(2-),[m^3 mol^-1]
  real(r8), PARAMETER :: SHF3P1=SYF3P1*DPH2O**2          !equilibrium constant for FePO4(s)+3H2O=Fe(OH)3(aq)+HPO4(2-)+2H(+),[mol^3 m^-9]
  real(r8), PARAMETER :: SYF4P1=SYF3P1/DPFE4             !equilibrium constant for FePO4(s)+H2O+3OH(-)=Fe(OH)4(-)+HPO4(2-),[m^6 mol^-2]
  real(r8), PARAMETER :: SHF4P1=SYF4P1*DPH2O**3          !equilibrium constant for FePO4(s)+4H2O=Fe(OH)4(-)+HPO4(2-)+3H(+),[mol^4 m^-12]
  real(r8), PARAMETER :: SHF0P2=SHF0P1/DPH2P             !equilibrium constant for FePO4(s)+2H(+)=Fe(3+)+H2PO4(-), [m^3 mol^-1]
  real(r8), PARAMETER :: SYF0P2=SHF0P2*DPH2O**2          !equilibrium constant for FePO4(s)+2H2O =Fe(3+)+H2PO4(-)+2OH(-),[mol^3 m^-9]
  real(r8), PARAMETER :: SYF1P2=SYF0P2/DPFE1             !equilibrium constant for FePO4(s)+2H2O=Fe(OH)(2+)+H2PO4(-)+OH(-),[mol^2 m^-6]
  real(r8), PARAMETER :: SHF1P2=SYF1P2/DPH2O             !equilibrium constant for FePO4(s)+H2O+H(+)=Fe(OH)(2+)+H2PO4(-),[-]
  real(r8), PARAMETER :: SPF2P2=SYF1P2/DPFE2             !equilibrium constant for FePO4(s)+2H2O=Fe(OH)2(+)+H2PO4(-),[mol m^-3]
  real(r8), PARAMETER :: SYF3P2=SPF2P2/DPFE3             !equilibrium constant for FePO4(s)+2H2O+OH(-)=Fe(OH)3(aq)+H2PO4(-),[-]
  real(r8), PARAMETER :: SHF3P2=SYF3P2*DPH2O             !equilibrium constant for FePO4(s)+3H2O=Fe(OH)3(aq)+H2PO4(-)+H(+),[mol^2 m^-6]
  real(r8), PARAMETER :: SYF4P2=SYF3P2/DPFE4             !equilibrium constant for FePO4(s)+2H2O+2OH(-)=Fe(OH)4(-)+H2PO4(-),[m^3 mol^-1]
  real(r8), PARAMETER :: SHF4P2=SYF4P2*DPH2O**2          !equilibrium constant for FePO4(s)+4H2O=Fe(OH)4(-)+H2PO4(-)+2H(+),[mol^3 m^-9]
  real(r8), PARAMETER :: SHCAD2=SPCAD/DPH2P              !equilibrium constant for CaHPO4(s)+H(+)=Ca(2+)+H2PO4(-),[-]
  real(r8), PARAMETER :: SYCAD2=SHCAD2*DPH2O             !equilibrium constant for CaHPO4(s)+H2O=Ca(2+)+H2PO4(-)+OH(-),[mol^2 m^-6] 
  real(r8), PARAMETER :: SHCAH1=SPCAH/(DPH2O*DPH1P**3)   !equilibrium constant for Ca5(PO4)3OH(s)+4H(+)=3HPO4(2-)+5Ca(2+)+H2O,[mol^3 m^-9]
  real(r8), PARAMETER :: SYCAH1=SHCAH1*DPH2O**4          !equilibrium constant for Ca5(PO4)3OH(s)+3H2O=5Ca(2+)+3HPO4(2-)+4OH(-),[mol^11 m^-33]
  real(r8), PARAMETER :: SHCAH2=SHCAH1/DPH2P**3          !equilibrium constant for Ca5(PO4)3OH(s)+7H(+)+3OH(-)=3HPO4(2-)+5Ca(2+)+4H2O,[m^9 mol^-3]
  real(r8), PARAMETER :: SYCAH2=SHCAH2*DPH2O**7          !equilibrium constant for Ca5(PO4)3OH(s)+3H2O=3HPO4(2-)+5Ca(2+)+4OH(-),[mol^11 mol^-33]

  real(r8), PARAMETER :: DUKM=1.0_r8                      !Km for Urea hydrolysis, [g m^-3]
  real(r8), PARAMETER :: DUKI=2.5_r8                      !Ki for urea hydrolysis, [g m^-3]
  real(r8), PARAMETER :: A0=1.0_r8                        !activity for solids, [-]
  real(r8), PARAMETER :: COOH=2.5E-04_r8                  !Cation exchange capacity of SOC, [eqv (gC)^-1]
  real(r8), PARAMETER :: CCAMX=10.0_r8                    !maximum Ca concentration, [mol m^-3]
  real(r8), PARAMETER :: RFertNH4SpecReleaz=1.0E-00_r8    !specific rate constants for NH4 release after fertilizer application, [h^-1]
  real(r8), PARAMETER :: RFertNH3SpecReleaz=1.0E-00_r8    !specific rate constants for NH3 release after fertilizer application, [h^-1]
  real(r8), PARAMETER :: RFertUreaSpecHydrol=1.0E-01_r8   !specific rate constants for Urea release after fertilizer application, [h^-1]
  real(r8), PARAMETER :: RFertNO3SpecReleaz=1.0E-00_r8    !specific rate constants for NO3 release after fertilizer application, [h^-1]
  real(r8), PARAMETER :: SPPO4=5.0E-02_r8                 !specific rate constants for H2PO4 release after fertilizer application, [h^-1]

  !terminate [label for variable parsing]
!
!     MRXN=number of cycles for solving reaction equilibria
!     TPD,TADA,TADC,TSL=adsorpn,solute equilibria
!
  integer, PARAMETER :: MRXN=1
  real(r8), parameter :: TPD=2.5E-03_r8   !1/h,rate constants for dissoln
  real(r8), PARAMETER :: TPDX=TPD/MRXN    !1/h
  real(r8), PARAMETER :: TADA=5.0E-02_r8
  real(r8), PARAMETER :: TADAX=TADA/MRXN
  real(r8), PARAMETER :: TADC=5.0E-02_r8
  real(r8), PARAMETER :: TADCX=TADC/MRXN
  real(r8), PARAMETER :: TADC0=TADC*1.0E-02_r8
  real(r8), PARAMETER :: TSL=5.0E-01_r8
  real(r8), PARAMETER :: TSLX=TSL/MRXN
!

end module SoluteParMod
