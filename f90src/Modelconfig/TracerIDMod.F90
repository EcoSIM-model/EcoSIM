module TracerIDMod
! !
! code to initialize tracers
  use MiniMathMod, only : addone
  use abortutils, only : destroy
  use data_kind_mod,    only: r8 => DAT_KIND_R8  
  use ElmIDMod
implicit none
  save
  CHARACTER(LEN=*), private, PARAMETER :: MOD_FILENAME=&
  __FILE__

  integer, parameter :: idg_CO2=1
  integer, parameter :: idg_CH4=2
  integer, parameter :: idg_O2 =3
  integer, parameter :: idg_N2 =4
  integer, parameter :: idg_N2O=5
  integer, parameter :: idg_H2 =6
  integer, parameter :: idg_AR =7
  integer, parameter :: idg_NH3=8
  integer, parameter :: idg_beg = idg_CO2  
  character(len=10), allocatable :: trcs_names(:)

  integer :: idg_NH3B
  integer :: ids_NH4,ids_NH4B
  integer :: ids_NO3,ids_NO3B
  integer :: ids_NO2,ids_NO2B
  integer :: ids_H1PO4,ids_H1PO4B
  integer :: ids_H2PO4,ids_H2PO4B

! be careful about the dom tracers and chemical elements associated
! with plant and soil organic matter
  integer :: idom_DOC,idom_DON,idom_DOP,idom_acetate
  integer :: idom_beg,idom_end

  integer :: idg_end
  integer :: ids_beg,ids_end
  integer :: ids_nut_beg,ids_nuts_beg,ids_nuts_end
  integer :: ids_nutb_beg,ids_nutb_end
! salt tracers
  integer :: idsalt_beg,idsalt_end,idsalt_mend

  integer :: idsalt_Al      ! Al
  integer :: idsalt_Fe      ! Fe
  integer :: idsalt_Hp      ! H+
  integer :: idsalt_Ca      ! Ca
  integer :: idsalt_Mg      ! Mg
  integer :: idsalt_Na      ! Na
  integer :: idsalt_K       ! K
  integer :: idsalt_OH      ! OH
  integer :: idsalt_SO4     ! SO4
  integer :: idsalt_Cl      ! Cl
  integer :: idsalt_CO3     ! CO3,
  integer :: idsalt_HCO3    ! HCO3
  integer :: idsalt_AlOH    ! AlOH
  integer :: idsalt_AlOH2   ! AlOH2
  integer :: idsalt_AlOH3   ! AlOH3
  integer :: idsalt_AlOH4   ! AlOH4
  integer :: idsalt_AlSO4   ! AlSO4
  integer :: idsalt_FeOH    ! FeOH
  integer :: idsalt_FeOH2   ! FeOH2
  integer :: idsalt_FeOH3   ! FeOH3
  integer :: idsalt_FeOH4   ! FeOH4
  integer :: idsalt_FeSO4   ! FeSO4
  integer :: idsalt_CaOH   ! CaOH
  integer :: idsalt_CaCO3   ! CaCO3
  integer :: idsalt_CaHCO3  ! CaHCO3
  integer :: idsalt_CaSO4   ! CaSO4
  integer :: idsalt_MgOH2   ! MgOH2
  integer :: idsalt_MgCO3   ! MgCO3
  integer :: idsalt_MgHCO3  ! MgHCO3
  integer :: idsalt_MgSO4   ! MgSO4
  integer :: idsalt_NaCO3   ! NaCO3
  integer :: idsalt_NaSO4   ! NaSO4
  integer :: idsalt_KSO4    ! KSO4
  integer :: idsalt_H0PO4   ! PO4
  integer :: idsalt_H3PO4   ! H3PO4
  integer :: idsalt_FeHPO4  ! FeHPO4
  integer :: idsalt_FeH2PO4 ! FeH2PO4
  integer :: idsalt_CaPO4   ! CaPO4
  integer :: idsalt_CaHPO4  ! CaHPO4
  integer :: idsalt_CaH4P2O8 ! Legacy name: dissolved CaH2PO4(+), one P atom; retain ID for compatibility.
  integer :: idsalt_MgHPO4 ! MgHPO4
!band salt
  integer :: idsalt_H0PO4B   ! PO4
  integer :: idsalt_H3PO4B   ! H3PO4
  integer :: idsalt_FeHPO4B  ! FeHPO4
  integer :: idsalt_FeH2PO4B ! FeH2PO4
  integer :: idsalt_CaPO4B   ! CaPO4
  integer :: idsalt_CaHPO4B  ! CaHPO4
  integer :: idsalt_CaH4P2O8B ! Legacy name: band dissolved CaH2PO4(+), one P atom.
  integer :: idsalt_MgHPO4B ! MgHPO4
  integer :: idsaltb_beg    !band begin
  integer :: idsaltb_end    !band end
  integer :: idsalt_psoil_beg,idsalt_psoil_end
  integer :: idsalt_pband_beg,idsalt_pband_end

! precipitated tracers
  integer :: idsp_beg,idsp_end
  integer :: idsp_beg_band
  integer :: idsp_CaHPO4 !
  integer :: idsp_Apatite     !hydroxyapatite
  integer :: idsp_AlOH3  !
  integer :: idsp_FeOH3  !
  integer :: idsp_CaCO3  !
  integer :: idsp_CaSO4
  integer :: idsp_AlPO4
  integer :: idsp_FePO4
  integer :: idsp_CaH4P2O8 ! Solid Ca(H2PO4)2, two P atoms per formula unit.
  integer :: idsp_CaHPO4B
  integer :: idsp_ApatiteBand
  integer :: idsp_AlPO4B
  integer :: idsp_FePO4B
  integer :: idsp_CaH4P2O8B ! Band solid Ca(H2PO4)2, two P atoms per formula unit.
  integer :: idsp_psoi_beg,idsp_psoi_end
  integer :: idsp_p_beg, idsp_p_end
  
! exchangeable tracers
  integer :: idx_CEC    ! XCEC,  cation exchange capacity, [mol d-2]
  integer :: idx_NH4    ! XN4, exchangeable NH4 non-band, [mol d-2]
  integer :: idx_Hp     ! XHY, exchangeable H , [mol d-2]
  integer :: idx_Al     ! XAL, exchangeable Al, [mol d-2]
  integer :: idx_Fe     ! XFE, exchangeable Fe, [mol d-2]
  integer :: idx_Ca     ! XCA, exchangeable Ca, [mol d-2]
  integer :: idx_Mg     ! XMG, exchangeable Mg , [mol d-2]
  integer :: idx_Na     ! XNA, exchangeable Na, [mol d-2]
  integer :: idx_K      ! XKA, exchangeable K, [mol d-2]
  integer :: idx_COOH   ! XHC, exchangeable COOH , [mol d-2]
  integer :: idx_AlOH2  ! XALO2, exchangeable AlOH2 , [mol d-2]
  integer :: idx_FeOH2  ! XFEO2, exchangeable Fe(OH)2, [mol d-2]
  integer :: idx_NH4B   ! XNB, exchangeable NH4 band, [mol d-2]
  integer :: idx_cation_end
  integer :: idx_cation_soil_end

  integer :: idx_AEC    ! XAEC, anion exchange capacity, [mol d-2]
  integer :: idx_OHe    ! XOH0, exchangeable OH- non-band, [mol d-2]
  integer :: idx_OH     ! XOH1, exchangeable OH  non-band, [mol d-2]
  integer :: idx_OHp    ! XOH2, exchangeable OH2  non-band, [mol d-2]
  integer :: idx_HPO4   ! XH1P, exchangeable HPO4  non-band, [mol d-2]
  integer :: idx_H2PO4  ! XH2P, exchangeable H2PO4  non-band, [mol d-2]
  integer :: idx_OHeB   ! XOH0B, exchangeable OH- band, [mol d-2]
  integer :: idx_OHB    ! XOH1B, exchangeable OH  band, [mol d-2]
  integer :: idx_OHpB   ! XOH2B, exchangeable OH2  band, [mol d-2]
  integer :: idx_HPO4B  ! XH1PB, exchangeable HPO4  band, [mol d-2]
  integer :: idx_H2PO4B ! XH2PB, exchangeable H2PO4  band, [mol d-2]
  integer :: idx_beg, idx_end
  integer :: idx_anion_soil_end
  integer :: ifertn_beg,ifertn_end
  integer :: ifertnb_beg,ifertnb_end
  integer :: ids_nuts
  integer :: idsalt_nuts

  real(r8), allocatable :: tracerSolc_max(:)  !maximum tracer concentration, solublity, g/m3 H2O
  type, public :: trc_def_type
   integer :: NGasTracers    !number of gas tracers
   integer :: NSolutTracers    !number of solute tracers
   integer :: NSaltTracers   !number of salt tracers
   integer :: NPrecipTracers    !number of precipitate tracers
   integer :: nxtracers    !number of exchangeable tracers
   integer :: NnutrientTracers    !number of nutrient tracers
   integer :: NFertNitro
   integer :: NFertNitrob
   integer :: NDOMS
  end type trc_def_type

  type(trc_def_type), public :: trc_confs
  public :: CleanUpTracerIDs
  contains

  subroutine InitTracerIDs(lsalt_model)
  implicit none
  logical, intent(in) :: lsalt_model


  ifertn_beg=ifert_N_nh4;ifertn_end=ifert_N_no3
  ifertnb_beg=ifert_N_nh4_band;ifertnb_end=ifert_N_no3_band

  ids_beg=idg_beg
! for better array manipulation of land-atmosphere exchange,
! banded NH3 is considered as (potential) gas too. However, there
! is no banded gas NH3 concentration.

  idg_NH3B = idg_NH3+1
  idg_end  = idg_NH3B;

  ids_nuts_beg = idg_NH3;  !the first nutrient tracer, including band
  ids_end      = idg_end        !initalize the solute counter
  ids_NH4B     = addone(ids_end);
  ids_NO3B     = addone(ids_end);
  ids_NO2B     = addone(ids_end);
  ids_H1PO4B   = addone(ids_end);
  ids_H2PO4B   = addone(ids_end);
  ids_nuts     = ids_H2PO4B-ids_NH4B

  ids_nutb_beg=idg_NH3B;ids_nutb_end=ids_H2PO4B

  ids_NH4   = addone(ids_end);
  ids_NO3   = addone(ids_end);
  ids_NO2   = addone(ids_end);
  ids_H1PO4 = addone(ids_end);
  ids_H2PO4 = addone(ids_end);

  ids_nut_beg  = ids_NH4;  !the first non-band non-gaseous nutrient tracer
  ids_nuts_end = ids_H2PO4;!the last non-band nutrient tracer

  allocate(tracerSolc_max(ids_beg:ids_end));tracerSolc_max=1.e3_r8; !set a large enough number by default

  tracerSolc_max(idg_NH3)=320._r8*14._r8/(1.e-3_r8*17._r8)   !based on NH3 solublity at 25C 
  tracerSolc_max(ids_NH4)=370._r8*14._r8/(1.e-3_r8*53.49_r8) !based on NH4Cl solubility at 25oC
  tracerSolc_max(ids_NO3)=912_r8*14._r8/(1.e-3_r8*84.995_r8)    !based on NaNO3 solubility at 25oC 

  tracerSolc_max(idg_NH3B)=tracerSolc_max(idg_NH3)
  tracerSolc_max(ids_NH4B)=tracerSolc_max(ids_NH4)
  tracerSolc_max(ids_NO3B)=tracerSolc_max(ids_NO3)

  allocate(trcs_names(ids_beg:ids_end))
  trcs_names(idg_CO2)    = 'CO2';trcs_names(idg_CH4)     = 'CH4'
  trcs_names(idg_O2)     = 'O2'; trcs_names(idg_N2)      = 'N2'
  trcs_names(idg_N2O)    = 'N2O';trcs_names(idg_H2)      = 'H2'
  trcs_names(idg_Ar)     = 'Ar'
  trcs_names(idg_NH3)    = 'NH3';trcs_names(idg_NH3B)    = 'NH3B'
  trcs_names(ids_NH4B)   = 'NH4B';trcs_names(ids_NO3B)   = 'NO3B'
  trcs_names(ids_NO2B)   = 'NO2B';trcs_names(ids_H1PO4B) = 'H1PO4B'
  trcs_names(ids_H2PO4B) = 'H2PO4B';trcs_names(ids_NH4)  = 'NH4'
  trcs_names(ids_NO3)    = 'NO3';trcs_names(ids_NO2)     = 'NO2'
  trcs_names(ids_H1PO4)  = 'H1PO4';trcs_names(ids_H2PO4) = 'H2PO4'

  
  idom_beg=1;idom_end=0
  idom_doc     = addone(idom_end)
  idom_don     = addone(idom_end)
  idom_dop     = addone(idom_end)
  idom_acetate = addone(idom_end)

  idsalt_beg=1;idsalt_end=0  
  !tracers involved in irrigation
  idsalt_Al     = addone(idsalt_end)
  idsalt_Fe     = addone(idsalt_end)
  idsalt_Ca     = addone(idsalt_end)
  idsalt_Mg     = addone(idsalt_end)
  idsalt_Na     = addone(idsalt_end)
  idsalt_K      = addone(idsalt_end)
  idsalt_SO4    = addone(idsalt_end)
  idsalt_Cl     = addone(idsalt_end)
  idsaltb_end   = idsalt_end
  if(lsalt_model)then    
    idsalt_Hp     = addone(idsalt_end)  
    idsalt_OH     = addone(idsalt_end)  
    idsalt_CO3    = addone(idsalt_end)
    idsalt_HCO3   = addone(idsalt_end)
    idsalt_mend   = idsalt_end
    idsalt_AlOH   = addone(idsalt_end)
    idsalt_AlOH2  = addone(idsalt_end)
    idsalt_AlOH3  = addone(idsalt_end)
    idsalt_AlOH4  = addone(idsalt_end)
    idsalt_AlSO4  = addone(idsalt_end)
    idsalt_FeOH   = addone(idsalt_end)
    idsalt_FeOH2  = addone(idsalt_end)
    idsalt_FeOH3  = addone(idsalt_end)
    idsalt_FeOH4  = addone(idsalt_end)
    idsalt_FeSO4  = addone(idsalt_end)
    idsalt_CaOH   = addone(idsalt_end)
    idsalt_CaCO3  = addone(idsalt_end)
    idsalt_CaHCO3 = addone(idsalt_end)
    idsalt_CaSO4  = addone(idsalt_end)
    idsalt_MgOH2  = addone(idsalt_end)
    idsalt_MgCO3  = addone(idsalt_end)
    idsalt_MgHCO3 = addone(idsalt_end)
    idsalt_MgSO4  = addone(idsalt_end)
    idsalt_NaCO3  = addone(idsalt_end)
    idsalt_NaSO4  = addone(idsalt_end)
    idsalt_KSO4   = addone(idsalt_end)

    idsalt_H0PO4    = addone(idsalt_end)
    idsalt_H3PO4    = addone(idsalt_end)
    idsalt_FeHPO4   = addone(idsalt_end)
    idsalt_FeH2PO4  = addone(idsalt_end)
    idsalt_CaPO4    = addone(idsalt_end)
    idsalt_CaHPO4   = addone(idsalt_end)
    idsalt_CaH4P2O8 = addone(idsalt_end)
    idsalt_MgHPO4   = addone(idsalt_end)

    idsalt_end  = idsalt_MgHPO4
    idsaltb_beg = idsalt_end+1
    idsaltb_end = idsalt_end
    idsalt_psoil_beg=idsalt_H0PO4;idsalt_psoil_end=idsalt_MgHPO4

    idsalt_nuts=idsalt_MgHPO4-idsalt_H0PO4

    idsalt_H0PO4B    = addone(idsaltb_end)   ! PO4
    idsalt_H3PO4B    = addone(idsaltb_end)   ! H3PO4
    idsalt_FeHPO4B   = addone(idsaltb_end)   ! FeHPO4
    idsalt_FeH2PO4B  = addone(idsaltb_end)   ! FeH2PO4
    idsalt_CaPO4B    = addone(idsaltb_end)   ! CaPO4
    idsalt_CaHPO4B   = addone(idsaltb_end)   ! CaHPO4
    idsalt_CaH4P2O8B = addone(idsaltb_end)   ! CaH4P2O8
    idsalt_MgHPO4B   = addone(idsaltb_end)   ! MgHPO4
    idsalt_pband_beg=idsalt_H0PO4B;idsalt_pband_end=idsalt_MgHPO4B

    idsaltb_end=idsalt_MgHPO4B

  endif
!  double check
  idsp_beg=1
  idsp_end=0;

  idsp_AlOH3=addone(idsp_end)     !Al(3+)+3OH(-)
  idsp_FeOH3=addone(idsp_end)     !Fe(3+)+3OH(-)
  idsp_CaCO3=addone(idsp_end)     !Ca(2+)+CO3(2-)
  idsp_CaSO4=addone(idsp_end)     !Ca(2+)+SO4(2-)

  !precipitate
  idsp_Apatite   =addone(idsp_end)     !hydroxyapatite
  idsp_AlPO4=addone(idsp_end)     !Al(3+)+PO4(3-)
  idsp_FePO4=addone(idsp_end)     !Fe(3+)+PO4(3-)
  idsp_CaHPO4=addone(idsp_end)    !Ca(2+)+H(+)+PO4(3-)
  idsp_CaH4P2O8=addone(idsp_end)  !Ca(2+)+4H(+)+2PO4(3-)
  idsp_psoi_beg=idsp_Apatite;idsp_psoi_end=idsp_CaH4P2O8

  idsp_ApatiteBand  =addone(idsp_end)
  idsp_AlPO4B=addone(idsp_end)
  idsp_FePO4B=addone(idsp_end)
  idsp_CaHPO4B=addone(idsp_end)
  idsp_CaH4P2O8B=addone(idsp_end)
  idsp_end=idsp_CaH4P2O8B
  idsp_beg_band=idsp_ApatiteBand
  idsp_p_beg=idsp_Apatite;idsp_p_end=idsp_CaH4P2O8B

  idx_beg=1
  idx_end=0

  !cations
  idx_CEC   = addone(idx_end)    ! XCEC, cation exchange capacity,  [mol d-2]
  idx_NH4   = addone(idx_end)    ! XN4,  exchangeable NH4 non-band, [mol d-2], X-NH4
  idx_Hp    = addone(idx_end)    ! XHY,  exchangeable H,            [mol d-2], X-H
  idx_Al    = addone(idx_end)    ! XAL,  exchangeable Al,           [mol d-2], X-Al
  idx_Fe    = addone(idx_end)    ! XFE,  exchangeable Fe,           [mol d-2], X-Fe
  idx_Ca    = addone(idx_end)    ! XCA,  exchangeable Ca,           [mol d-2], X-Ca
  idx_Mg    = addone(idx_end)    ! XMG,  exchangeable Mg,           [mol d-2], X-Mg
  idx_Na    = addone(idx_end)    ! XNA,  exchangeable Na,           [mol d-2], X-Na
  idx_K     = addone(idx_end)    ! XKA,  exchangeable K,            [mol d-2], X-K
  idx_COOH  = addone(idx_end)  ! XHC,    exchangeable COOH,         [mol d-2], X-COOH
  idx_AlOH2 = addone(idx_end) ! XALO2,   exchangeable AlOH2,        [mol d-2], X-Al(OH)2
  idx_FeOH2 = addone(idx_end) ! XFEO2,   exchangeable Fe(OH)2,      [mol d-2], X-Fe(OH)2

  idx_NH4B            = addone(idx_end)   ! XNB, exchangeable NH4 band, [mol d-2]
  idx_cation_end      = idx_NH4B
  idx_cation_soil_end = idx_FeOH2

  !anions
  idx_AEC   = addone(idx_end)   ! XAEC, anion exchange capacity,      [mol d-2]
  idx_OHe   = addone(idx_end)   ! XOH0, exchangeable OH- non-band,    [mol d-2], X-O(-)
  idx_OH    = addone(idx_end)   ! XOH1, exchangeable OH  non-band,    [mol d-2], X-OH
  idx_OHp   = addone(idx_end)   ! XOH2, exchangeable OH2  non-band,   [mol d-2], X-OH2(+)
  idx_HPO4  = addone(idx_end)  ! XH1P,  exchangeable HPO4  non-band,  [mol d-2], X-HPO4(-)
  idx_H2PO4 = addone(idx_end) ! XH2P,   exchangeable H2PO4  non-band, [mol d-2], X-H2PO4

  idx_OHeB           = addone(idx_end) ! XOH0B,  exchangeable OH- band,    [mol d-2],
  idx_OHB            = addone(idx_end)  ! XOH1B, exchangeable OH  band,    [mol d-2]
  idx_OHpB           = addone(idx_end) ! XOH2B,  exchangeable OH2  band,   [mol d-2]
  idx_HPO4B          = addone(idx_end) ! XH1PB,  exchangeable HPO4  band,  [mol d-2]
  idx_H2PO4B         = addone(idx_end)! XH2PB,   exchangeable H2PO4  band, [mol d-2]
  idx_end            = idx_H2PO4B
  idx_anion_soil_end = idx_H2PO4

  trc_confs%NGasTracers              = idg_end-idg_beg
  trc_confs%NSolutTracers            = ids_end-ids_beg+1
  if(lsalt_model)trc_confs%NSaltTracers = idsaltb_end-idsalt_beg+1
  trc_confs%NPrecipTracers           = idsp_end-idsp_beg+1
  trc_confs%nxtracers                = idx_end-idx_beg+1
  trc_confs%NnutrientTracers         = ids_nuts_end-ids_nut_beg+1
  trc_confs%NFertNitro               = ifertn_end-ifertn_beg+1
  trc_confs%NFertNitrob              = ifertnb_end-ifertnb_beg+1
  trc_confs%NDOMS                    = idom_end-idom_beg+1

  return
  write(104,*)'ids_beg=',ids_beg,'ids_end=',ids_end
  write(104,*)'idg_CO2   =',idg_CO2
  write(104,*)'idg_CH4   =',idg_CH4
  write(104,*)'idg_O2    =',idg_O2
  write(104,*)'idg_N2    =',idg_N2
  write(104,*)'idg_N2O   =',idg_N2O
  write(104,*)'idg_H2    =',idg_H2
  write(104,*)'idg_NH3   =',idg_NH3

  write(104,*)'idg_NH3B  =',idg_NH3B
  write(104,*)'ids_NH4B  =',ids_NH4B
  write(104,*)'ids_NO3B  =',ids_NO3B
  write(104,*)'ids_NO2B  =',ids_NO2B
  write(104,*)'ids_H1PO4B=',ids_H1PO4B
  write(104,*)'ids_H2PO4B=',ids_H2PO4B

  write(104,*)'ids_NH4   =',ids_NH4
  write(104,*)'ids_NO3   =',ids_NO3
  write(104,*)'ids_NO2   =',ids_NO2
  write(104,*)'ids_H1PO4 =',ids_H1PO4
  write(104,*)'ids_H2PO4 =',ids_H2PO4

  write(104,*)'ids_nuts_beg=',ids_nuts_beg,'ids_nuts_beg=',ids_nuts_end

  stop
  end subroutine InitTracerIDs

  subroutine CleanUpTracerIDs()
  implicit none

  call destroy(trcs_names)
  end subroutine CleanUpTracerIDs
end module TracerIDMod

