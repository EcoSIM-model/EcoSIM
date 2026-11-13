module NitroPars
!
! DESCRIPTION:
! code defining parameters for nitro
! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSIMCtrlMod, only : micpar_file_in
  use netcdf
  use ncdio_pio  
  use EcoSimConst
! --- Jing Tao --- needed for the SPOSC_in carrier below, which is dimensioned (jskenc,jcplxc)
!          = (4 kinetic components, 5 organo-microbial complexes), matching micpar%SPOSC.
  use EcoSIMConfig, only : jskenc, jcplxc
  implicit none
  public
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  save
! SUBSTRATE DECOMPOSITION BY MICROBIAL POPULATIONS

  real(r8) :: ORAD         !microbial radius, [m]
  real(r8) :: BIOS         !microbial number density, [n m-3]
  real(r8) :: BIOA         !microbial surface area, [m2 m-3]
  real(r8) :: DCKI         !inhibition of decomposition by microbial concentration, [g C m-3]
  real(r8) :: RCCX         !maximum remobilization of microbial N, [-]
  real(r8) :: RCCQ         !maximum P recycling fractions , [-]
  real(r8) :: RCCZ         !minimum fractions for bacteria C recycling,[-]
  real(r8) :: RCCY         !maximum remobilization of microbial P, [-]
  real(r8) :: FPRIM        !fraction of nonstructural transferred with priming, [-] 
  real(r8) :: FPRIMM       !fraction of microbial C,N,P transferred with priming, [-]
  real(r8) :: OMGR         !rate constant for transferring nonstructural to structural microbial C,[h-1]
  real(r8) :: OQKI         !DOC product inhibition constant for decomposition, [g C m-3]
  real(r8) :: H2KI         !H2 normalizing mass concentation for Gibbs free energy calculation, [gH m-3]
  real(r8) :: OAKI         !acetate normalizing mass concentation for Gibbs free energy calculation, [gC m-3]
  real(r8) :: COMKI        !Km to slow microbial decomposition with low microbial C, [g micr C g-1 subs C]
  real(r8) :: COMKM        !Km to slow microbial maintenance respiration with low microbial C, [g micr C g-1 subs C]
  real(r8) :: CKC          !controls C remobilization of microbial C, [g C g-1 C]
  real(r8) :: FOSCZ0       !rate for mixing surface litter,  [h-1]
  real(r8) :: FOSCZL       !rate for mixing subsurface litter, [h-1]
  real(r8) :: FMN          !minimum ratio of total biological demand for any substrate by any microbial population,[-]
  real(r8) :: DCKM0        !Km for SOC decomposition, [g C g-1 soil]
  real(r8) :: DCKML        !Km for SOC decomposition, [g C g-1 soil] 
  real(r8) :: VMXO         !specific oxidation rates for all bacteria, [g C g-1C h-1]
  real(r8) :: VMXF         !specific oxidation rates for all fungi, [g C g-1C h-1] 
  real(r8) :: VMXCH4gAcet         !specific acetate oxidation rates for acetoclastic methanogens,[g C g-1C h-1] 
  real(r8) :: VMXNH3Oxi         !specific oxidation rates for ammonia oxidizers, [g  g-1C h-1] 
  real(r8) :: VMXNO2Oxi           !specific oxidation rates for nitrite oxidizers, [g  g-1C h-1] 
  real(r8) :: VMXCH4OxiAero     !specific oxidation rates for methanotrophs, [g   g-1C h-1] 
  real(r8) :: VMXCH4gH2         !specific oxidation rates for hydrogenotrophic methanogens,[g g-1C h-1] 
  real(r8) :: OQKM         !Km for DOC uptake by heterotrophs bacteria and fungi, [g C m-3]
  real(r8) :: OQKA         !Km for acetate uptake by heterotrophic fermenters, [g C m-3]
  real(r8) :: OQKAM        !Km for acetate uptake  by acetotrophic methanogens,[g C m-3]
  real(r8) :: CCKM         !Km for CO2 uptake,  [g C m-3]
  real(r8) :: CCK4         !Km for CH4 uptake, [g C m-3]
  real(r8) :: ZHKM         !Km for NH4 uptake by nitrifiers, [gN m-3]
  real(r8) :: ZNKM         !Km for NO2 uptake by nitrifiers, [gN m-3]
  real(r8) :: Z3KM         !Km for NO3 uptake by denitrifiers,[gN m-3]
  real(r8) :: Z2KM         !Km for NO2 uptake by denitrifiers, [gN m-3]
  real(r8) :: Z1KM         !Km for N2O uptake by denitrifiers, [gN m-3]
  real(r8) :: Z4MX         !maximum uptake rate for NH4 uptake kinetics by all microbial functional groups,[g N m-2 h-1]
  real(r8) :: Z4KU         !Km for NH4 uptake kinetics by all microbial functional groups,[g N m-3]
  real(r8) :: Z4MN         !minimum concentration for NH4 uptake kinetics by all microbial functional groups,[g N m-3]
  real(r8) :: ZOMX         !maximum uptake rate for NO3 uptake kinetics by all microbial functional groups, [g N m-2 h-1] 
  real(r8) :: ZOKU         !Km for NO3 uptake kinetics by all microbial functional groups, [g N m-3]
  real(r8) :: ZOMN         !minimum concentration for NO3 uptake kinetics by all microbial functional groups, [g N m-3]
  real(r8) :: HPMX         !maximum rate for H2PO4 uptake kinetics by all microbial functional groups, [g P m-2 h-1]
  real(r8) :: HPKU         !Km for H2PO4 uptake kinetics by all microbial functional groups,[g P m-3]
  real(r8) :: HPMN         !Minimum concentration for H2PO4 uptake kinetics by all microbial functional groups,[g P m-3]
  real(r8) :: ZFKM         !Km for N2 uptake by diazotrophs, [g N m-3]
  real(r8) :: H2KM         !Km for H2 uptake by hydrogenotrophic methanogens, [g H m-3]
  real(r8) :: EAMO10       !efficiency of CO2 conversion to biomass by AMO NC10,[-]
  real(r8) :: EAMO2D       !efficiency of CH4 conversion to biomass by AMO ANME2d,[-]
  real(r8) :: ECNH         !efficiency of CO2 conversion to biomass (CH2O) by ammonia oxidizers,[-]
  real(r8) :: ECNO         !efficiency of CO2 conversion to biomass (CH2O) by nitrite oxidizers,[-]
  real(r8) :: ECHO         !efficiency of CO2 conversion to biomass by methane oxidizers,[-]
  real(r8) :: VMX2AMONC10  !maximum specific uptake rate of NO2(-) by AMO NC-10, [h-1]
  real(r8) :: VMX3AMO2D    !maximum specific uptake rate of NO2(-) by AMO ANME-2D, [h-1]

  real(r8) :: RNFNI        !parameter for nitrification inhibition,[-]
  real(r8) :: ZHKI         !inhibition of nitrification inhibition by NH3, [g N m-3]
  real(r8) :: VMKI         !product inhibn for NOx reduction by denitrifiers, [g N m-3]
  real(r8) :: VHKI         !product inhibn for NH3 oxidation by nitrifiers, [g N m-3]
  real(r8) :: OXKA         !Km for O2 uptake by nitrifiers, [g O m-3]
  real(r8) :: EOMC         !energy requirements for microbial growth of aerobic bacteria, [kJ g-1 C]
  real(r8) :: EOMD         !energy requirements for microbial growth of denitrifiers, [kJ g-1 C]
  real(r8) :: EOMG         !energy requirements for microbial growth of fungi, [kJ g-1 C]
  real(r8) :: EOMF         !energy requirements for microbial growth of fermenters, [kJ g-1 C]
  real(r8) :: EOMH         !energy requirements for microbial biomass growth of methanogens, [kJ g-1 C]
  real(r8) :: EOMN         !energy requirements for microbial biomass growth of free-living diazotrophs, [kJ g-1 C]
  
  real(r8) :: GO2X         !free energy yields of redox reactions for DOC->CO2, [kJ g-1 C]
  real(r8) :: GH4X         !free energy yields of redox reactions for CH4/O2->CO2, [kJ g-1 C]
  real(r8) :: GCHX         !free energy yields of redox reactions for DOC-acetate, [kJ g-1 C]
  real(r8) :: GO2A         !free energy yields of O2-driven redox reactions for acetate-CO2, [kJ g-1 C]
  real(r8) :: GC4X         !free energy yields of redox reactions for acetate-CH4, [kJ g-1 C]
  real(r8) :: GCOX         !free energy yields of redox reactions for CO2-CH4, [kJ g-1 C]
  real(r8) :: GC2HX        !free energy yields of redox reactions for CH2O-CH4, [kJ g-1 C]
  real(r8) :: GNOX         !free energy yields of redox reactions for NO3-NO2, NO2-N2O,N2O-N2, [kJ gN-1]
  real(r8) :: GN3CX        !free energy yields of redox reaction for NO3-CH4, [kJ gC-1]
  real(r8) :: GN2X         !free energy yields of redox reactions for N2-NH3, [kJ g-1 N]
  real(r8) :: EN2D         !growth respiraiton efficiency for AMO ANME-2D, [-]

  real(r8) :: EN2X            !growth respiration efficiency for aerobic N2 fixation, [-]
  real(r8) :: EN2Y            !growth respiration efficiency for anaerobic N2 fixation, [-]
  real(r8) :: EO2X            !growth respiration efficiency for aerobic bacteria (DOC), [-], 1gC Energy/(1gC Energy+1gC biom)
  real(r8) :: EH4X            !growth respiration efficiency for aerobic Ch4 oxidizer, [-]
  real(r8) :: EO2G            !growth respiration efficiency for fungi, [-]
  real(r8) :: EO2D            !growth respiration efficiency for denitrifiers (aerobic), [-], 1gC Energy/(1gC Energy+1gC biom)
  real(r8) :: ENFX            !growth respiration efficiency for diazotrophs, [-], 1gC Energy/(1gC Energy+1gC biom)
  real(r8) :: ENOX            !growth respiration efficiency for denitrifiers (anaerobic), [-]
  real(r8) :: EO2A            !growth respiration efficiency for aerobic bacteria (acetate), [-]
  real(r8) :: TSORP           !sorption rate constant for OHC, [h-1]
  real(r8) :: HSORP           !sorption volume coefficient for OHC, [-]
  real(r8) :: SPOHC           !specific decomposition rate constant for adsorbed SOC, [g subs. C g-1 micr. C]
  real(r8) :: SPOHA           !specific decomposition rate constant for adsorbed acetate, [g subs. C g-1 micr. C]
  real(r8) :: RMOM            !specific maintenance respiration, [g C g-1 N h-1]
  real(r8) :: SPORC(2)        !specific decomposition rate constant microbial residue,  [g C g-1 N h-1]
  real(r8) :: SPOMC(2)        !specific decomposition rate constant microbial biomass,  [g C g-1 N h-1]
  real(r8) :: EN2F(8)         !N fixation yield from C oxidation, [g N g-1 C]
  real(r8) :: EFIRE(2,21:22)  !partition coefficient for N loss as NH3 and P loss as PO4 during combustion, [g gC-1]
  real(r8) :: CyanoInocC      !inoculation cyanobacteria biomass [gC m-2]
  real(r8) :: CyanoInoc_depth !inocoluation depth of cyanobacteria, [m]
! --- Jing Tao --- NitroPars-side carrier for the SOLID-SOM specific decomposition rate, promoted from a
!          hardcoded literal in MicBGCPars::SetPars to an optional MicrobePars.nc input.
!          Units h-1. Shape (jskenc,jcplxc) = (4 kinetic components: protein/carbohydrate/
!          cellulose/lignin) x (5 complexes: woody-litr/fine-litr/manure/POM/humus).
!          Named SPOSC_in, NOT SPOSC, because MicBGCPars::SetPars already aliases
!          SPOSC => this%SPOSC in an associate block and the names would clash.
!          Holds the BASE (pre-boost) rates: SetPars applies the x1.5 litter boost AFTER
!          copying this in, exactly as it did to the hardcoded values.
!          Why promote it: SPOSC is the first-order control on the solid-SOM hydrolysis flux
!          (MicBGCFGMod.F90:1442, Fs linear in SPOSC up to a 0.5 cap), while its sibling rate
!          constants SPORC/SPOHC/SPOMC above are ALREADY .nc inputs -- an undocumented split
!          that forced calibration to compensate through pool size instead (6x CORGC in R2).
  real(r8) :: SPOSC_in(jskenc,jcplxc)
  contains

  subroutine initNitroPars
  implicit none

  ORAD = 1.e-6_r8     !
  BIOS = 1.e-6_r8/(4._r8/3._r8*PICON*ORAD**3)
  BIOA = BIOS*4.0_r8*PICON*ORAD**2
  FMN  = 5.0E-03_r8
  H2KI = 1.0_r8
  OAKI = 12._r8
  CyanoInoc_depth = 1.e-3_r8   ![1mm]
  CyanoInocC = 1.e-2_r8        ![gC m-2]
! --- Jing Tao --- SPOSC_in default, set UNCONDITIONALLY -- deliberately BEFORE the .not.ReadPars() branch
!  below, not inside it. The other parameters only need a default when the whole read FAILS, but
!  SPOSC is a NEW variable that no existing MicrobePars.nc carries: a file can be read SUCCESSFULLY
!  and still not contain SPOSC. Setting it inside the else-branch would leave it uninitialised in
!  exactly the common case (an old .nc reads fine, SPOSC absent), which is a garbage decomposition
!  rate rather than a fallback. Setting it here is what makes the promotion default-off.
!
!  These 20 values are the BASE (pre-boost) rates copied verbatim from the literal they replace in
!  MicBGCPars::SetPars; the x1.5 litter-complex boost still happens there, after the copy. Column-
!  major, so each row of 4 below is one complex K over the 4 kinetic components M:
!    K=1 woody-litr, K=2 fine-litr, K=3 manure, K=4 POM, K=5 humus
!    M=1 protein,    M=2 carbohydrate, M=3 cellulose, M=4 lignin        [h-1]
  SPOSC_in = reshape((/7.5_r8,7.5_r8,1.5_r8,0.5_r8 &
    ,7.5_r8,7.5_r8,1.5_r8,0.5_r8     &
    ,7.5_r8,7.5_r8,1.5_r8,0.5_r8     &
    ,0.05_r8,0.0_r8,0.0_r8,0.0_r8 &
    ,0.05_r8,0.0167_r8,0.0_r8,0.0_r8/),shape(SPOSC_in))

  if(.not.ReadPars())then
    DCKI          = 2.5_r8
    RCCX          = 0.833_r8
    RCCQ          = 0.833_r8
    RCCZ          = 0.167_r8
    RCCY          = 0.833_r8
    FPRIM         = 5.0E-02_r8
    FPRIMM        = 1.e-6_r8
    OMGR          = 2.5E-01_r8  !
    OQKI          = 1.2E+03_r8

    OAKI          = 12.0_r8
    COMKI         = 1.0E-03_r8
    COMKM         = 1.0E-04_r8
    CKC           = 1.0E-03_r8
    FOSCZ0        = 2.0E-02_r8
    FOSCZL        = 2.0E-06_r8
    DCKM0         = 1.0E+03_r8
    DCKML         = 1.0E+03_r8
    VMX2AMONC10   = 0.01_R8
    VMX3AMO2D     = 0.01_R8 
    VMXO          = 0.125_r8
    VMXF          = 0.125_r8
    VMXCH4gAcet   = 0.125_r8
    VMXCH4gH2     = 0.125_r8
    VMXNH3Oxi     = 0.375_r8
    VMXNO2Oxi     = 0.25_r8
    VMXCH4OxiAero = 0.375_r8
    OQKM          = 1.2E+01_r8
    OQKA          = 1.2E+01_r8
    OQKAM         = 1.2E+01_r8
    CCKM          = 0.15_r8
    CCK4          = 1.2E-03_r8
    ZHKM          = 0.014_r8
    ZNKM          = 1.4_r8
    Z3KM          = 1.4_r8
    Z2KM          = 1.4_r8
    Z1KM          = 0.014_r8
    Z4MX          = 5.0E-03_r8
    Z4KU          = 0.40_r8
    Z4MN          = 0.0125_r8
    ZOMX          = 5.0E-03_r8
    ZOKU          = 0.35_r8
    ZOMN          = 0.03_r8
    HPMX          = 1.0E-03_r8
    HPKU          = 0.075_r8
    HPMN          = 0.002_r8
    ZFKM          = 0.14_r8
    H2KM          = 0.01_r8
    EAMO10        = 0.1_r8
    EAMO2d        = 0.1_r8
    ECNH          = 0.30_r8
    ECNO          = 0.10_r8
    ECHO          = 0.75_r8
    RNFNI         = 2.0E-04_r8
    ZHKI          = 7.0E+03_r8
    VMKI          = 0.25_r8
    VHKI          = 15.0_r8
    OXKA          = 0.16_r8
    EOMC          = 25.0_r8
    EOMD          = 37.5_r8
    EOMG          = 37.5_r8
    EOMF          = 75.0_r8
    EOMH          = 25.0_r8
    EOMN          = 75.0_r8
    GO2X          = 37.5_r8
    GH4X          = 66.5_r8
    GCHX          = 4.50_r8
    GC4X          = 3.00_r8
    GCOX          = 11.00_r8
    GC2HX         = 13.5_r8
    GNOX          = 10.0_r8
    GN3CX         = 41.08_r8
    GN2X          = 187.5_r8
    TSORP         = 0.5_r8
    HSORP         = 1.0_r8
    SPOHC         = 0.25_r8
    SPOHA         = 0.25_r8
    RMOM          = 0.010_r8
    SPORC         = (/7.5_r8,1.5_r8/)                                 !hydrolysis rate of microbial residue
    SPOMC         = (/1.0E-02_r8,0.1E-02_r8/)                         !basal mortality rates
  endif

  GO2A       = GO2X-GCHX
  EN2X       = GO2X/GN2X
  EN2Y       = GCHX/GN2X
  EO2X       = 1._r8/(1._r8+GO2X/EOMC)   !1gC produces GO2X energy, fixes GO2X/EOMC gC, so the fraction EO2X is used for producing energy through respiraiton
  EH4X       = 1._r8/(1._r8+GH4X/EOMC)   !aerobic methane oxidizer respiraiton efficiency
  EO2G       = 1._r8/(1._r8+GO2X/EOMG)   !
  EO2D       = 1._r8/(1._r8+GO2X/EOMD)
  ENFX       = 1._r8/(1._r8+GO2X/EOMN)
  ENOX       = 1._r8/(1._r8+GNOX/EOMC)
  EO2A       = 1._r8/(1._r8+GO2A/EOMC)
  EN2D       = 1._r8/(1._r8+GN3CX/EOMH)

  EN2F  = (/0.0_r8,0.0_r8,0.0_r8,0.0_r8,0.0_r8,EN2X,EN2Y,EN2X/)
  EFIRE = reshape((/1.0_r8,1.0_r8,0.917_r8,0.167_r8/),shape(EFIRE))

  end subroutine initNitroPars

!------------------------------------------------------------------------------------------
  function ReadPars()result(ans)

  implicit none
  logical :: ans
  type(file_desc_t) :: ncfid
! --- Jing Tao --- descriptor + presence flag for the guarded optional SPOSC read below.
  type(Var_desc_t)  :: vardesc_SPOSC
  logical           :: readv_SPOSC

  ans =.false.

  if (len_trim(micpar_file_in)==0)return
  
  call ncd_pio_openfile(ncfid, micpar_file_in, ncd_nowrite)

  call ncd_getvar(ncfid,'DCKI',DCKI)
  call ncd_getvar(ncfid,'RCCX',RCCX)
  call ncd_getvar(ncfid,'RCCQ',RCCQ)
  call ncd_getvar(ncfid,'RCCZ',RCCZ)
  call ncd_getvar(ncfid,'RCCY',RCCY)
  call ncd_getvar(ncfid,'FPRIM',FPRIM)
  call ncd_getvar(ncfid,'FPRIMM',FPRIMM)
  call ncd_getvar(ncfid,'OMGR',OMGR)
  call ncd_getvar(ncfid,'OQKI',OQKI)
  call ncd_getvar(ncfid,'COMKI',COMKI)
  call ncd_getvar(ncfid,'COMKM',COMKM)
  call ncd_getvar(ncfid,'CKC',CKC)
  call ncd_getvar(ncfid,'FOSCZ0',FOSCZ0)
  call ncd_getvar(ncfid,'FOSCZL',FOSCZL)
  call ncd_getvar(ncfid,'DCKM0',DCKM0)
  call ncd_getvar(ncfid,'DCKML',DCKML)
  call ncd_getvar(ncfid,'VMXO',VMXO)
  call ncd_getvar(ncfid,'VMXF',VMXF)
  call ncd_getvar(ncfid,'VMXCH4gAcet',VMXCH4gAcet)
  call ncd_getvar(ncfid,'VMXNO2Oxi',VMXNO2Oxi)
  call ncd_getvar(ncfid,'VMXNH3Oxi',VMXNH3Oxi)
  call ncd_getvar(ncfid,'VMXCH4OxiAero',VMXCH4OxiAero)
  call ncd_getvar(ncfid,'VMXCH4gH2',VMXCH4gH2)
  call ncd_getvar(ncfid,'OQKM',OQKM)
  call ncd_getvar(ncfid,'OQKA',OQKA)
  call ncd_getvar(ncfid,'OQKAM',OQKAM)
  call ncd_getvar(ncfid,'CCKM',CCKM)
  call ncd_getvar(ncfid,'CCK4',CCK4)
  call ncd_getvar(ncfid,'ZHKM',ZHKM)
  call ncd_getvar(ncfid,'ZNKM',ZNKM)
  call ncd_getvar(ncfid,'Z3KM',Z3KM)
  call ncd_getvar(ncfid,'Z2KM',Z2KM)
  call ncd_getvar(ncfid,'Z1KM',Z1KM)
  call ncd_getvar(ncfid,'Z4MX',Z4MX)
  call ncd_getvar(ncfid,'Z4KU',Z4KU)
  call ncd_getvar(ncfid,'Z4MN',Z4MN)
  call ncd_getvar(ncfid,'ZOMX',ZOMX)
  call ncd_getvar(ncfid,'ZOKU',ZOKU)
  call ncd_getvar(ncfid,'ZOMN',ZOMN)
  call ncd_getvar(ncfid,'HPMX',HPMX)
  call ncd_getvar(ncfid,'HPKU',HPKU)
  call ncd_getvar(ncfid,'HPMN',HPMN)
  call ncd_getvar(ncfid,'ZFKM',ZFKM)
  call ncd_getvar(ncfid,'H2KM',H2KM)
  call ncd_getvar(ncfid,'ECNH',ECNH)
  call ncd_getvar(ncfid,'ECNO',ECNO)        
  call ncd_getvar(ncfid,'ECHO',ECHO)
  call ncd_getvar(ncfid,'RNFNI',RNFNI)
  call ncd_getvar(ncfid,'ZHKI',ZHKI)
  call ncd_getvar(ncfid,'VMKI',VMKI)
  call ncd_getvar(ncfid,'VHKI',VHKI)
  call ncd_getvar(ncfid,'OXKA',OXKA)
  call ncd_getvar(ncfid,'EOMC',EOMC)
  call ncd_getvar(ncfid,'EOMD',EOMD)
  call ncd_getvar(ncfid,'EOMG',EOMG)
  call ncd_getvar(ncfid,'EOMF',EOMF)
  call ncd_getvar(ncfid,'EOMH',EOMH)
  call ncd_getvar(ncfid,'EOMN',EOMN)
  call ncd_getvar(ncfid,'GO2X',GO2X)
  call ncd_getvar(ncfid,'GH4X',GH4X)
  call ncd_getvar(ncfid,'GCHX',GCHX)
  call ncd_getvar(ncfid,'GC4X',GC4X)
  call ncd_getvar(ncfid,'GCOX',GCOX)          
  call ncd_getvar(ncfid,'GNOX',GNOX)
  call ncd_getvar(ncfid,'GN2X',GN2X)
  call ncd_getvar(ncfid,'TSORP',TSORP)
  call ncd_getvar(ncfid,'HSORP',HSORP)
  call ncd_getvar(ncfid,'SPOHC',SPOHC)
  call ncd_getvar(ncfid,'SPOHA',SPOHA)
  call ncd_getvar(ncfid,'RMOM',RMOM)
  call ncd_getvar(ncfid,'SPORC',SPORC)
  call ncd_getvar(ncfid,'SPOMC',SPOMC)

! --- Jing Tao --- SPOSC is read only IF PRESENT -- this must not be a bare ncd_getvar like the lines above.
!  ncd_getvar -> check_ret -> endrun ABORTS the run when a variable is missing (ncdio_pio.F90:418-422),
!  and every MicrobePars.nc written before this change lacks SPOSC. A bare read would therefore turn
!  every existing input file into a hard crash. check_var only sets a flag, so this leaves the
!  initNitroPars default in place when the variable is absent and the run is bit-identical to before.
!  print_err=.false. keeps a normal old-format file from logging a spurious error.
  call check_var(ncfid,'SPOSC',vardesc_SPOSC,readv_SPOSC,print_err=.false.)
  if(readv_SPOSC) call ncd_getvar(ncfid,'SPOSC',SPOSC_in)
  call ncd_getvar(ncfid,'VMX2AMONC10',VMX2AMONC10)
  call ncd_getvar(ncfid,'VMX3AMO2D',VMX3AMO2D)

  call ncd_pio_closefile(ncfid)
  ans=.true.
  end function ReadPars
end module NitroPars
