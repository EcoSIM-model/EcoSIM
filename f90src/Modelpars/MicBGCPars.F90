module MicBGCPars

! define microbial parameters
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils, only : destroy
  use ElmIDMod
  use EcoSIMConfig
implicit none
  private
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  
  integer, public :: NumGuild_Heter_Aerob_Bact      = 1
  integer, public :: NumGuild_Heter_Aerob_Fung      = 1
  integer, public :: NumGuild_Heter_Facul_Dent      = 1
  integer, public :: NumGuild_Heter_Aerob_N2Fixer   = 1
  integer, public :: NumGuild_Heter_Anaer_N2Fixer   = 1
  integer, public :: NumGuild_Heter_Anaer_Fermentor = 1
  integer, public :: NumGuild_Heter_AcetoMethanogen = 1
  integer, public :: NumGuild_Autor_H2genMethanogen = 1
  integer, public :: NumGuild_Autor_AmoniaOxidBact  = 1
  integer, public :: NumGuild_Autor_NitritOxidBact  = 1
  integer, public :: NumGuild_Autor_AerobMethOxid   = 1
  integer, public :: NumGuild_Autor_ANMO_ANME2d     = 1
  integer, public :: NumGuild_Autor_ANMO_ANMENC10   = 1
  integer, public :: NumGuild_Mixor_CYNO_Bacter     = 1 !while cyanobacteria is a mixotroph, it is now grouped into the autotrohpic complex

  type, public :: MicParType
  real(r8), pointer :: ORCI(:,:)                   !allocation of initial residue to kinetic components, [-]
  real(r8), pointer :: FL(:)                       !allocation to microbial kinetic fractions, [-]
  real(r8), pointer :: rNCOMC(:,:,:)               !maximum/minimum mass based heterotrophic microbial N:C, [gN gC-1]
  real(r8), pointer :: rPCOMC(:,:,:)               !maximum/minimum mass based heterotrophic microbial P:C,  [gP gC-1]
  real(r8), pointer :: rNCOMCAutor(:,:)            !maximum/minimum mass based autotrophic microbial N:C, [gN gC-1]
  real(r8), pointer :: rPCOMCAutor(:,:)            !maximum/minimum mass based autotrophic microbial P:C, [gP gC-1]
  real(r8), pointer :: rNCOMC_ave(:,:,:)           !group average maximum/minimum mass based microbial N:C, [gN gC-1]
  real(r8), pointer :: rPCOMC_ave(:,:,:)           !group average maximum/minimum mass based microbial P:C, [gP gC-1]
!  real(r8), pointer :: rNCOMCAutor_ave(:,:)        !group average maximum/minimum mass based microbial N:C, [gN gC-1]
!  real(r8), pointer :: rPCOMCAutora_ave(:,:)       !group average maximum/minimum mass based microbial P:C, [gP gC-1]
  real(r8), pointer :: DOSA(:)                     !rate constant for litter colonization by heterotrophs, [h-1]
  real(r8), pointer :: SPOSC(:,:)                  !specific decomposition rate constant, [h-1]
  real(r8), pointer :: CNOFC(:,:)                  !Fractions of initial litter to allocate N to kinetic components,[-]
  real(r8), pointer :: CPOFC(:,:)                  !Fractions of iniital litter to allocate P to kinetic components, [-]
  real(r8), pointer :: CNRH(:)                     !Default N:C ratios in SOC complexes,[gN gC-1]
  real(r8), pointer :: CPRH(:)                     !Default P:C ratios in SOC complexes, [gN gC-1]
  real(r8), pointer :: OMCF(:)                     !Initial fractional composition of heterotrophic microbial biomass, [gC gC-1]
  real(r8), pointer :: OMCA(:)                     !Initial fractional composition of autotrophic microbial biomass, [gC gC-1]
  integer  :: FG_guilds_heter(NumMicbHFunGrupsPerCmplx)  !# of guilds
  integer  :: FG_guilds_autor(NumMicbAFunGrupsPerCmplx)  !# of guilds

  !terminate  [label for variable parsing]
  integer :: jcplx   !# of microbe-substrate complexes
  integer :: jsken   !# of kinetic components of the substrates
  integer :: NumMicbHFunGrupsPerCmplx             !# of heterotrophic functional groups
  integer :: NumMicbAFunGrupsPerCmplx             !# of autotrohpic functional groups
  integer :: k_woody_comp
  integer :: k_fine_comp
  integer :: k_manure
  integer :: k_POM
  integer :: k_humus
  integer :: mid_AutoAmmoniaOxidBacter
  integer :: mid_AutoNitriteOxidBacter
  integer :: mid_AutoAeroCH4OxiBacter
  integer :: mid_AutoH2GenoCH4GenArchea
  integer :: mid_HeterAerobBacter
  integer :: mid_Facult_DenitBacter
  integer :: mid_Aerob_Fungi
  integer :: mid_fermentor
  integer :: mid_HeterAcetoCH4GenArchea
  integer :: mid_HeterAerobN2Fixer
  integer :: mid_HeterAnaerobN2Fixer
  integer :: mid_AutoAMOANME2D  !anaerobic methane oxidizer ANME-2d
  integer :: mid_AutoAMONC10    !anaerobic methane oxidizer NC10
  integer :: mid_HeterMixtCynoBacter     !cyano bacteria
  integer :: ndbiomcp   !number of necrobiomass components
  integer :: nlbiomcp   !number of living biomass components

  real(r8), pointer :: OMCI(:,:)                    !fraction  of initial SOC/litterfall in microbial biomass
  real(r8), pointer :: OHCK(:)                      !fractions of initial SOC in adsorbed C
  real(r8), pointer :: OMCK(:)                      !fractions of initial SOC as living microibal biomass
  real(r8), pointer :: ORCK(:)                      !fractions of initial SOC as litter
  real(r8), pointer :: OQCK(:)                      !fractions of initial SOC as DOC
  logical,  pointer :: is_activeMicrbFungrpAutor(:) !logical switch for autotrophic group
  logical,  pointer :: is_activeMicrbFungrpHeter(:) !logical switch for heterotrophic group
  logical,  pointer :: is_aerobic_hetr(:)           !logical flag for aerobic heterotrophs
  logical,  pointer :: is_anaerobic_hetr(:)         !logical flag for anaerobic heterotrophs
  logical,  pointer :: is_aerobic_autor(:)          !logical flag for aerobic autotrophs
  logical,  pointer :: is_litter(:)
  logical,  pointer :: is_finelitter(:)
  logical,  pointer :: is_CO2_autotroph(:)
  character(len=16) :: kiname(jskenc)
  character(len=16) :: cplxname(1:jcplxc)
  character(len=16) :: hmicname(NumMicbHFunGrupsPerCmplx)
  character(len=16) :: amicname(NumMicbAFunGrupsPerCmplx)
  character(len=16) :: micresb(NumDeadMicrbCompts)      !residual biomass name
  character(len=16) :: micbiom(1:NumLiveMicrbCompts)        !microbial biomass pool name
  integer, pointer :: JGniH(:)   !hetetroph guid indices for organic-microbial complex
  integer, pointer :: JGnfH(:)   !hetetroph guid indices for organic-microbial complex
  integer, pointer :: JGniA(:)   !autotroph guid indices for autotrophic-microbial complex
  integer, pointer :: JGnfA(:)   !autotroph guid indices for autotrophic-microbial complex

  integer :: NumMicrobAutoTrophCmplx        !total number of microbial guilds in the autotrophic complex
  integer :: NumHetetr1MicCmplx         !total number of microbial guilds in one organic-microbial complex
  integer :: NumOfLitrCmplxs                !number of litter organo-microbial complexes, plant litter + manure
  integer :: NumOfPlantLitrCmplxs           !number of plant litter complexs, woody + fine litter
  integer :: NumLiveHeterBioms              !total number of biomass components in one heterotrophic OM complexes
  integer :: NumLiveAutoBioms               !total number of biomass components in the autotroph complex
  integer :: iprotein
  integer :: icarbhyro
  integer :: icellulos
  integer :: ilignin
  contains
    procedure, public  :: Init
    procedure, public  :: SetPars
    procedure, public  :: get_micb_id
    procedure, private :: InitAllocate
    procedure, public  :: destroy      =>DestructMicBGCPar
  end type MicParType


contains

  subroutine Init(this)
  implicit none
  class(MicParType) :: this

  !organic matter is grouped into five complexes, including woody(1),
  ! non-woody(2), manure(3), litter, POC(4) and humus(5) (g Mg-1)

  this%ndbiomcp = NumDeadMicrbCompts  !number of necrobiomass components
  this%nlbiomcp = NumLiveMicrbCompts  !number of living biomass components
  
  this%jcplx                   = jcplxc         !number of microbe-substrate complexes
  this%jsken                   = jskenc         !number of kinetic components of the substrates

  this%NumMicbAFunGrupsPerCmplx = NumMicbAFunGrupsPerCmplx
  this%NumMicbHFunGrupsPerCmplx = NumMicbHFunGrupsPerCmplx

  !woody, non_woody litter and manure are defined as litter
  allocate(this%is_litter(1:this%jcplx));this%is_litter(:)=.false.
  allocate(this%is_finelitter(1:this%jcplx));this%is_finelitter(:)=.false.

  this%k_woody_comp                = 1;                   this%is_litter(this%k_woody_comp) = .true.
  this%k_fine_comp                 = this%k_woody_comp+1; this%is_litter(this%k_fine_comp) = .true.
  this%NumOfPlantLitrCmplxs        = this%k_fine_comp
  this%k_manure                    = this%k_fine_comp+1;   this%is_litter(this%k_manure)   = .true.
  this%is_litter(this%k_fine_comp) = .true.
  this%is_litter(this%k_manure)    = .true.
  this%is_finelitter(this%k_fine_comp)=.true.
  this%iprotein  = 1
  this%icarbhyro = 2
  this%icellulos = 3
  this%ilignin   = 4

  this%NumOfLitrCmplxs = this%k_manure
  this%k_POM           = this%k_manure+1
  this%k_humus         = this%k_POM+1
  
  this%kiname(1)   = 'protein'
  this%kiname(2)   = 'carbhydro'
  this%kiname(3)   = 'cellulose'
  this%kiname(4)   = 'lignin'
  this%cplxname(1) = 'woodylitr'
  this%cplxname(2) = 'folialitr'
  this%cplxname(3) = 'manure'
  this%cplxname(4) = 'pom'
  this%cplxname(5) = 'humus'
  this%hmicname(1) = 'aerohetrob'
  this%hmicname(2) = 'faculdenit'
  this%hmicname(3) = 'aerofungi'
  this%hmicname(4) = 'aneroferm'
  this%hmicname(5) = 'acetmethg'
  this%hmicname(6) = 'aeron2fix'
  this%hmicname(7) = 'aneron2fix'
  this%hmicname(8) = 'mixocyanob'

  this%amicname(1) = 'amoniaoxib'
  this%amicname(2) = 'nititeoxib'
  this%amicname(3) = 'aeromethtp'
  this%amicname(4) = 'anme2d'
  this%amicname(5) = 'hydromethg'  
  this%amicname(6) = 'anmenc10'
  
  this%micresb(1)  = 'kinetic'
  this%micresb(2)  = 'recalcitrant'
  this%micbiom(1)  = 'kinetic'
  this%micbiom(2)  = 'recalcitrant'
  this%micbiom(3)  = 'reserve'

  !set up functional group ids
  ! five om-complexes
  this%mid_HeterAerobBacter     = 1
  this%mid_Facult_DenitBacter     = 2
  this%mid_Aerob_Fungi            = 3
  this%mid_Fermentor              = 4
  this%mid_HeterAcetoCH4GenArchea = 5
  this%mid_HeterAerobN2Fixer      = 6
  this%mid_HeterAnaerobN2Fixer    = 7
  this%mid_HeterMixtCynoBacter    = 8

  !the autotrophic complex
  this%mid_AutoAmmoniaOxidBacter  = 1
  this%mid_AutoNitriteOxidBacter  = 2
  this%mid_AutoAeroCH4OxiBacter   = 3
  this%mid_AutoAMOANME2D          = 4  
  this%mid_AutoH2GenoCH4GenArchea = 5
  this%mid_AutoAMONC10            = 6

  this%FG_guilds_heter = 0
  this%FG_guilds_autor = 0
  this%FG_guilds_heter(this%mid_HeterAerobBacter)       = NumGuild_Heter_Aerob_Bact
  this%FG_guilds_heter(this%mid_Aerob_Fungi)            = NumGuild_Heter_Aerob_Fung
  this%FG_guilds_heter(this%mid_Facult_DenitBacter)     = NumGuild_Heter_Facul_Dent
  this%FG_guilds_heter(this%mid_HeterAerobN2Fixer)      = NumGuild_Heter_Aerob_N2Fixer
  this%FG_guilds_heter(this%mid_HeterAnaerobN2Fixer)    = NumGuild_Heter_Anaer_N2Fixer
  this%FG_guilds_heter(this%mid_fermentor)              = NumGuild_Heter_Anaer_Fermentor
  this%FG_guilds_heter(this%mid_HeterAcetoCH4GenArchea) = NumGuild_Heter_AcetoMethanogen
  this%FG_guilds_heter(this%mid_HeterMixtCynoBacter)    = NumGuild_Mixor_CYNO_Bacter

  this%FG_guilds_autor(this%mid_AutoH2GenoCH4GenArchea) = NumGuild_Autor_H2genMethanogen
  this%FG_guilds_autor(this%mid_AutoAmmoniaOxidBacter)  = NumGuild_Autor_AmoniaOxidBact
  this%FG_guilds_autor(this%mid_AutoNitriteOxidBacter)  = NumGuild_Autor_NitritOxidBact
  this%FG_guilds_autor(this%mid_AutoAeroCH4OxiBacter)   = NumGuild_Autor_AerobMethOxid
  this%FG_guilds_autor(this%mid_AutoAMOANME2D)          = NumGuild_Autor_ANMO_ANME2d
  this%FG_guilds_autor(this%mid_AutoAMONC10)            = NumGuild_Autor_ANMO_ANMENC10

  call this%Initallocate()

  this%is_aerobic_hetr(this%mid_HeterAerobBacter) = .true.
  this%is_aerobic_hetr(this%mid_Facult_DenitBacter) = .true.
  this%is_aerobic_hetr(this%mid_Aerob_Fungi)        = .true.
  this%is_aerobic_hetr(this%mid_HeterAerobN2Fixer)  = .true.
  this%is_aerobic_hetr(this%mid_HeterMixtCynoBacter) = .true.
  this%is_anaerobic_hetr(this%mid_fermentor)       = .true.
  this%is_anaerobic_hetr(this%mid_HeterAnaerobN2Fixer) = .true.

  this%is_aerobic_autor(this%mid_AutoAmmoniaOxidBacter) =.true.
  this%is_aerobic_autor(this%mid_AutoNitriteOxidBacter) =.true.
  this%is_aerobic_autor(this%mid_AutoAeroCH4OxiBacter)=.true.
  
  this%is_activeMicrbFungrpAutor(this%mid_AutoAmmoniaOxidBacter)  = .true.
  this%is_activeMicrbFungrpAutor(this%mid_AutoNitriteOxidBacter)  = .true.
  this%is_activeMicrbFungrpAutor(this%mid_AutoAeroCH4OxiBacter)   = .true.
  this%is_activeMicrbFungrpAutor(this%mid_AutoH2GenoCH4GenArchea) = .true.
  this%is_activeMicrbFungrpAutor(this%mid_AutoAMOANME2D)          = .true.
  this%is_activeMicrbFungrpAutor(this%mid_AutoAMONC10)            = .true.

  this%is_activeMicrbFungrpHeter(this%mid_HeterAerobBacter)       = .true.
  this%is_activeMicrbFungrpHeter(this%mid_Facult_DenitBacter)     = .true.
  this%is_activeMicrbFungrpHeter(this%mid_Aerob_Fungi)            = .true.
  this%is_activeMicrbFungrpHeter(this%mid_fermentor)              = .true.
  this%is_activeMicrbFungrpHeter(this%mid_HeterAcetoCH4GenArchea) = .true.
  this%is_activeMicrbFungrpHeter(this%mid_HeterAerobN2Fixer)      = .true.
  this%is_activeMicrbFungrpHeter(this%mid_HeterAnaerobN2Fixer)    = .true.
  this%is_activeMicrbFungrpHeter(this%mid_HeterMixtCynoBacter)    = .true.

  this%is_CO2_autotroph(this%mid_AutoAmmoniaOxidBacter)  = .true.
  this%is_CO2_autotroph(this%mid_AutoNitriteOxidBacter)  = .true.
  this%is_CO2_autotroph(this%mid_AutoH2GenoCH4GenArchea) = .true.
  this%is_CO2_autotroph(this%mid_AutoAMONC10)            = .true.
  
  end subroutine Init
!------------------------------------------------------------------------------------------

  subroutine SetPars(this)
! --- Jing Tao --- imported at SUBROUTINE scope, not module scope, so MicBGCPars does not acquire a
!  module-wide dependency on NitroPars for one array. NitroPars does not use MicBGCPars, so there
!  is no cycle. Boot order makes this safe: initNitroPars (which sets/reads SPOSC_in) runs from
!  InitModules -> initNitro -> MicBGCFGMod.F90:61, BEFORE readi -> MicPar%SetPars assigns it here.
  use NitroPars, only : SPOSC_in

  implicit none
  class(MicParType) :: this

  real(r8) :: OMCI1(NumLiveMicrbCompts,1:this%jcplx)  !allocation of living biomass to different components
  integer :: K,M,NGL,N
  associate(                         &
    OHCK        => this%OHCK,        &
    OMCK        => this%OMCK,        &
    ORCK        => this%ORCK,        &
    OQCK        => this%OQCK,        &
    ORCI        => this%ORCI,        &
    OMCI        => this%OMCI,        &
    CNOFC       => this%CNOFC,       &
    CPOFC       => this%CPOFC,       &
    rNCOMC      => this%rNCOMC,      &
    rPCOMC      => this%rPCOMC,      &
    rNCOMCAutor => this%rNCOMCAutor, &
    rPCOMCAutor => this%rPCOMCAutor, &
    FL          => this%FL,          &
    DOSA        => this%DOSA,        &
    SPOSC       => this%SPOSC,       &
    CNRH        => this%CNRH,        &
    CPRH        => this%CPRH,        &
    OMCF        => this%OMCF,        &
    OMCA        => this%OMCA         &
  )
  OHCK=real((/0.05,0.05,0.05,0.05,0.05/),r8)      
  OMCK=real((/0.01,0.01,0.01,0.01,0.01/),r8)      
  ORCK=real((/0.25,0.25,0.25,0.25,0.25/),r8)      
  OQCK=real((/0.005,0.005,0.005,0.005,0.005/),r8) 

  OMCI1=reshape(real(   &
    (/0.010,0.050,0.005,&
      0.050,0.050,0.005,&
      0.050,0.050,0.005,&
      0.010,0.050,0.005,&
      0.010,0.050,0.005/),r8),shape(OMCI1))

  ORCI=reshape(real((/0.01,0.05,0.01,0.05,0.01,0.05 &
     ,0.001,0.005,0.001,0.005/),r8),shape(ORCI))


  CNRH = (/3.33E-02_r8,3.33E-02_r8,3.33E-02_r8,5.00E-02_r8,12.50E-02_r8/)
  CPRH = (/3.33E-03_r8,3.33E-03_r8,3.33E-03_r8,5.00E-03_r8,12.50E-03_r8/)
  OMCF = (/0.20_r8,0.20_r8,0.30_r8,0.20_r8,0.050_r8,0.025_r8,0.025_r8,0.0_r8/)
  OMCA = (/0.6_r8,0.2_r8,0.05_r8,0.025_r8,0.1_r8,0.025_r8/)*0.1_r8

  OMCI(1:NumLiveMicrbCompts,:)=OMCI1

  ! CNOFC,CPOFC=fractions to allocate N,P to kinetic components
  ! rNCOMC,rPCOMC=maximum N:C and P:C ratios in microbial biomass

  CNOFC(1:jskenc,this%k_woody_comp) = real((/0.0050,0.0050,0.0050,0.0200/),r8)  !woody
  CPOFC(1:jskenc,this%k_woody_comp) = real((/0.0005,0.0005,0.0005,0.0020/),r8)  !woody
  CNOFC(1:jskenc,this%k_fine_comp)  = real((/0.0200,0.0200,0.0200,0.0200/),r8)  !non-woody
  CPOFC(1:jskenc,this%k_fine_comp)  = real((/0.0020,0.0020,0.0020,0.0020/),r8)  !non-woody
  CNOFC(1:jskenc,this%k_manure)     = real((/0.0200,0.0200,0.0200,0.0200/),r8)   !manure
  CPOFC(1:jskenc,this%k_manure)     = real((/0.0020,0.0020,0.0020,0.0020/),r8)   !manure  

  !microbial traits
  DOSA=(/0.25E-03_r8,0.25_r8,0.25_r8,0.25_r8,0.25_r8/)

! --- Jing Tao --- SPOSC (solid-SOM specific decomposition rate, h-1) is now sourced from NitroPars::SPOSC_in
!  instead of the literal that used to sit here. SPOSC_in defaults to those exact same 20 values in
!  initNitroPars and is overwritten ONLY if the MicrobePars.nc actually carries an SPOSC variable, so
!  with any pre-existing input file this line reproduces the previous behaviour bit-for-bit.
!  The x1.5 litter boost below is deliberately KEPT and still applied AFTER the copy: the .nc (like the
!  old literal) holds BASE, pre-boost rates. Moving the boost into the .nc, or dropping it, would either
!  double-boost the litter complexes or silently slow them by a third.
!  Rationale for making this tunable at all: it is the first-order control on the hydrolysis flux
!  (MicBGCFGMod.F90:1442, this file's sole consumer), and its siblings SPORC/SPOHC/SPOMC were already
!  .nc inputs.
  SPOSC = SPOSC_in

  SPOSC(:,1:this%NumOfLitrCmplxs)=SPOSC(:,1:this%NumOfLitrCmplxs)*1.5_r8

  FL(1:2)=real((/0.55,0.45/),r8)
  !set stoichiometry of heterotrophs
  D95: DO K=1,this%jcplx
    DO  N=1,this%NumMicbHFunGrupsPerCmplx
      !the set up of rNCOMC and rPCOMC may be 
      !set up based on external input
      IF(N.EQ.this%mid_Aerob_Fungi)THEN
        !Fungi      
        DO NGL=this%JGniH(n),this%JGnfH(n)
          rNCOMC(iLbiom_kinetic,NGL,K) = 0.15_r8           !NC ratio of kinetic biomass
          rNCOMC(iLbiom_struct,NGL,K) = 0.09_r8            !NC ratio of structural biomass
          rPCOMC(iLbiom_kinetic,NGL,K) = 0.015_r8          !PC ratio of kinetic biomass
          rPCOMC(iLbiom_struct,NGL,K) = 0.009_r8           !PC ratio of structural biomass
        ENDDO
        this%rNCOMC_ave(iLbiom_kinetic,N,K)=0.15_r8           
        this%rNCOMC_ave(iLbiom_struct,N,K)=0.09_r8            
        this%rPCOMC_ave(iLbiom_kinetic,N,K)=0.015_r8
        this%rPCOMC_ave(iLbiom_struct,N,K)=0.009_r8

        !bacteria  
      ELSE
        do NGL=this%JGniH(n),this%JGnfH(n)
          rNCOMC(iLbiom_kinetic,NGL,K)=0.225_r8
          rNCOMC(iLbiom_struct,NGL,K)=0.135_r8
          rPCOMC(iLbiom_kinetic,NGL,K)=0.0225_r8
          rPCOMC(iLbiom_struct,NGL,K)=0.0135_r8
        enddo
        this%rNCOMC_ave(iLbiom_kinetic,N,K)=0.225_r8
        this%rNCOMC_ave(iLbiom_struct,N,K)=0.135_r8
        this%rPCOMC_ave(iLbiom_kinetic,N,K)=0.0225_r8
        this%rPCOMC_ave(iLbiom_struct,N,K)=0.0135_r8
      ENDIF

      !reserve biomass
      do NGL=this%JGniH(n),this%JGnfH(n)
        rNCOMC(iLbiom_reserve,NGL,K)=DOT_PRODUCT(FL,rNCOMC(1:2,NGL,K))     !NC ratio
        rPCOMC(iLbiom_reserve,NGL,K)=DOT_PRODUCT(FL,rPCOMC(1:2,NGL,K))     !PC ratio
      enddo
      this%rNCOMC_ave(iLbiom_reserve,N,K)=DOT_PRODUCT(FL,this%rNCOMC_ave(1:2,N,K))
      this%rPCOMC_ave(iLbiom_reserve,N,K)=DOT_PRODUCT(FL,this%rPCOMC_ave(1:2,N,K))
    enddo
  ENDDO D95

  !set stoichiometry of autotrophs
  !the following may also be set based on external input
  DO  N=1,this%NumMicbAFunGrupsPerCmplx
    do NGL=this%JGniA(n),this%JGnfA(n)
      rNCOMCAutor(iLbiom_kinetic,NGL) = 0.225_r8
      rNCOMCAutor(iLbiom_struct,NGL)  = 0.135_r8
      rPCOMCAutor(iLbiom_kinetic,NGL) = 0.0225_r8
      rPCOMCAutor(iLbiom_struct,NGL)  = 0.0135_r8
    enddo
 !   this%rNCOMCAutor_ave(iLbiom_kinetic,N)=0.225_r8
 !   this%rNCOMCAutor_ave(iLbiom_struct,N)=0.135_r8

!    this%rPCOMCAutora_ave(iLbiom_kinetic,N)=0.0225_r8
!    this%rPCOMCAutora_ave(iLbiom_struct,N)=0.0135_r8

    do NGL=this%JGniA(n),this%JGnfA(n)
      rNCOMCAutor(iLbiom_reserve,NGL)=DOT_PRODUCT(FL,rNCOMCAutor(1:2,NGL))
      rPCOMCAutor(iLbiom_reserve,NGL)=DOT_PRODUCT(FL,rPCOMCAutor(1:2,NGL))
    enddo
!    this%rNCOMCAutor_ave(iLbiom_reserve,N)=DOT_PRODUCT(FL,this%rNCOMCAutor_ave(1:2,N))
!    this%rPCOMCAutora_ave(iLbiom_reserve,N)=DOT_PRODUCT(FL,this%rPCOMCAutora_ave(1:2,N))
  enddo

  end associate
  end subroutine SetPars
!------------------------------------------------------------------------------------------

  subroutine Initallocate(this)
  implicit none
  class(MicParType) :: this

  integer :: NumMicbAFunGrupsPerCmplx
  integer :: NumMicbHFunGrupsPerCmplx
  integer :: jcplx
  integer :: jsken
  integer :: n, kh,ka
  
  NumMicbAFunGrupsPerCmplx  =this%NumMicbAFunGrupsPerCmplx
  NumMicbHFunGrupsPerCmplx  =this%NumMicbHFunGrupsPerCmplx
  jcplx =this%jcplx
  jsken =this%jsken
  allocate(this%JGniH(NumMicbHFunGrupsPerCmplx))
  allocate(this%JGnfH(NumMicbHFunGrupsPerCmplx))

  allocate(this%JGniA(NumMicbAFunGrupsPerCmplx))
  allocate(this%JGnfA(NumMicbAFunGrupsPerCmplx))

  kh=1
  ka=1
  this%NumMicrobAutoTrophCmplx=0
  this%NumHetetr1MicCmplx=0
  !replace the functional group specification with external input later
  do N=1,NumMicbHFunGrupsPerCmplx
    this%JGniH(n)                = kh
    kh                           = kh+this%FG_guilds_heter(N)
    this%JGnfH(n)                = kh-1
    if(this%JGnfH(n)>=this%JGniH(n))then    
      this%NumHetetr1MicCmplx      = this%NumHetetr1MicCmplx+this%JGnfH(n)-this%JGniH(n)+1
    endif
  enddo

  do N=1,NumMicbAFunGrupsPerCmplx
    this%JGniA(n)                = ka
    ka                           = ka+this%FG_guilds_autor(N)
    this%JGnfA(n)                = ka-1
    if(this%JGnfA(n)>=this%JGniA(n))then
      this%NumMicrobAutoTrophCmplx = this%NumMicrobAutoTrophCmplx+this%JGnfA(n)-this%JGniA(n)+1
    endif
  enddo

  this%NumLiveHeterBioms=this%nlbiomcp*this%NumHetetr1MicCmplx
  this%NumLiveAutoBioms=this%nlbiomcp*this%NumMicrobAutoTrophCmplx

  allocate(this%DOSA(1:jcplx))
  allocate(this%SPOSC(jsken,1:jcplx))
  allocate(this%OHCK(1:jcplx))
  allocate(this%OMCK(1:jcplx))
  allocate(this%ORCK(1:jcplx))
  allocate(this%OQCK(1:jcplx))
  allocate(this%ORCI(NumDeadMicrbCompts,1:jcplx))
  allocate(this%OMCI(NumLiveMicrbCompts,1:jcplx))
  allocate(this%rNCOMC(NumLiveMicrbCompts,this%NumHetetr1MicCmplx,1:jcplx))
  allocate(this%rPCOMC(NumLiveMicrbCompts,this%NumHetetr1MicCmplx,1:jcplx))
  allocate(this%rNCOMCAutor(NumLiveMicrbCompts,this%NumMicrobAutoTrophCmplx))
  allocate(this%rPCOMCAutor(NumLiveMicrbCompts,this%NumMicrobAutoTrophCmplx))
  allocate(this%rNCOMC_ave(NumLiveMicrbCompts,NumMicbHFunGrupsPerCmplx,1:jcplx))
  allocate(this%rPCOMC_ave(NumLiveMicrbCompts,NumMicbHFunGrupsPerCmplx,1:jcplx))
!  allocate(this%rNCOMCAutor_ave(NumLiveMicrbCompts,NumMicbAFunGrupsPerCmplx))
!  allocate(this%rPCOMCAutora_ave(NumLiveMicrbCompts,NumMicbAFunGrupsPerCmplx))

  allocate(this%CNOFC(jsken,1:this%NumOfLitrCmplxs))
  allocate(this%CPOFC(jsken,1:this%NumOfLitrCmplxs))
  allocate(this%CNRH(1:jcplx))
  allocate(this%CPRH(1:jcplx))
  allocate(this%OMCF(NumMicbHFunGrupsPerCmplx))
  allocate(this%OMCA(NumMicbAFunGrupsPerCmplx))
  allocate(this%FL(2))
  allocate(this%is_activeMicrbFungrpAutor(NumMicbAFunGrupsPerCmplx)); this%is_activeMicrbFungrpAutor=.false.
  allocate(this%is_activeMicrbFungrpHeter(NumMicbHFunGrupsPerCmplx)); this%is_activeMicrbFungrpHeter=.false.
  allocate(this%is_CO2_autotroph(NumMicbAFunGrupsPerCmplx)); this%is_CO2_autotroph=.false.
  allocate(this%is_aerobic_hetr(NumMicbHFunGrupsPerCmplx)); this%is_aerobic_hetr=.false.
  allocate(this%is_anaerobic_hetr(NumMicbHFunGrupsPerCmplx));this%is_anaerobic_hetr=.false.
  allocate(this%is_aerobic_autor(NumMicbAFunGrupsPerCmplx));this%is_aerobic_autor=.false.
  end subroutine InitAllocate
!------------------------------------------------------------------------------------------

  subroutine DestructMicBGCPar(this)
  implicit none
  class(MicParType) :: this

  call destroy(this%OMCI)
  call destroy(this%OHCK)
  call destroy(this%OMCK)
  call destroy(this%ORCK)
  call destroy(this%OQCK)
  call destroy(this%ORCI)
  call destroy(this%OMCI)
  call destroy(this%rNCOMC)
  call destroy(this%rPCOMC)
  call destroy(this%rNCOMCAutor)
  call destroy(this%rPCOMCAutor)
  call destroy(this%CNOFC)
  call destroy(this%CPOFC)
  call destroy(this%FL)
  call destroy(this%DOSA)
  call destroy(this%SPOSC)
  call destroy(this%CNOFC)
  call destroy(this%CPOFC)
  call destroy(this%CNRH)
  call destroy(this%CPRH)
  call destroy(this%OMCF)
  call destroy(this%OMCA)
  call destroy(this%is_activeMicrbFungrpAutor)
  call destroy(this%is_CO2_autotroph)
  call destroy(this%is_aerobic_hetr)
  call destroy(this%is_anaerobic_hetr)
  end subroutine DestructMicBGCPar

!------------------------------------------------------------------------------------------
  pure function get_micb_id(this,M,NGL)result(id)
  !return the id of biomass component M for
  !hetetroph guild NGL
  !1,2,3: labile, recalcitrant, reserve
  implicit none
  class(MicParType), intent(in) :: this
  integer, intent(in) :: M   !biomass component
  integer, intent(in) :: NGL !guild id
  integer :: id

  id=this%nlbiomcp*(NGL-1)+M

  end function get_micb_id

!------------------------------------------------------------------------------------------

  pure function is_group_defined(this,igroup,isauto)result(isdef)
  implicit none
  class(MicParType), intent(in) :: this
  integer, intent(in) :: igroup
  logical, intent(in) :: isauto  
  logical :: isdef

  if(isauto)then
    isdef=igroup == this%mid_AutoAmmoniaOxidBacter     .or. &
       igroup == this%mid_AutoNitriteOxidBacter        .or. & 
       igroup == this%mid_AutoAeroCH4OxiBacter         .or. &
       igroup == this%mid_AutoH2GenoCH4GenArchea
  else
    isdef=igroup == this%mid_HeterAerobBacter  .or. &
       igroup == this%mid_Facult_DenitBacter     .or. &
       igroup == this%mid_Aerob_Fungi            .or. &
       igroup == this%mid_fermentor              .or. &
       igroup == this%mid_HeterAcetoCH4GenArchea    .or. &
       igroup == this%mid_HeterAerobN2Fixer          .or. &
       igroup == this%mid_HeterAnaerobN2Fixer
  endif
  end function is_group_defined
end module MicBGCPars
