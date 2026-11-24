module ElmIDMod
!
! DESCRIPTION:
!  Chemical element ids
implicit none
  character(len=*),private, parameter :: mod_filename =&
   __FILE__
  integer, parameter :: ielmc=1    !carbon element
  integer, parameter :: ielmn=2    !nitrogen element
  integer, parameter :: ielmp=3    !phosphorus element
  integer, parameter :: NumPlantChemElms=3   !totally three elements
  ! IDS for microbial biomass components
  integer, parameter :: iLbiom_kinetic = 1 !live kinetic biomass 
  integer, parameter :: iLbiom_struct  = 2 !live structural biomas
  integer, parameter :: iLbiom_reserve = 3 !live reserve biomass
  
  integer, parameter :: iDbiom_labile  = 1
  integer, parameter :: iDbiom_recalc  = 2
! erosion model options
  integer, parameter :: ieros_noaction       = -1
  integer, parameter :: ieros_frzthawelv     = 0
  integer, parameter :: ieros_frzthaweros    = 1
  integer, parameter :: ieros_frzthawsom     = 2
  integer, parameter :: ieros_frzthawsomeros = 3
! irrigation
  integer, parameter :: iIrrig_swc=0  !by soil water content
  integer, parameter :: iIrrig_cwp=1  !by canopy water potential
! water flux direction
  integer, parameter :: iWestEastDirection   = 1   !east-west direction
  integer, parameter :: iNorthSouthDirection = 2   !north-south direction
  integer, parameter :: iVerticalDirection   = 3   !vertical direction
  integer, parameter :: iFront             = 1
  integer, parameter :: iBehind              = 2
! soil properties
  integer, parameter :: isoi_fc    = 1   !field capacity
  integer, parameter :: isoi_wp    = 2   !wilting point
  integer, parameter :: isoi_scnv  = 3  !vertical hydraulic conductivity
  integer, parameter :: isoi_scnh  = 4  !horizontal hydraulic conductivity
  integer, parameter :: isoi_set   = 0
  integer, parameter :: isoi_unset = 1
! plant harvest
  integer, parameter :: iplthvst_leaf       =1 !leaf
  integer, parameter :: iplthvst_finenonleaf=2 !fine non-leaf
  integer, parameter :: iplthvst_stalk      =3 !woody
  integer, parameter :: iplthvst_stdead     =4 !standing dead
! photosynthesis
  integer, parameter :: ic4_photo=4
  integer, parameter :: ic3_photo=3
! fertilizer
  integer, parameter :: iamendtyp_fert   = 1
  integer, parameter :: iAmendtyp_plantRes  = 2
  integer, parameter :: iAmendtyp_Manure = 3
  integer, parameter :: ifert_on         = 1
  integer, parameter :: ifert_off               = 0
  integer, parameter :: ifert_N_nh4               = 1
  integer, parameter :: ifert_N_nh3               = 2
  integer, parameter :: ifert_N_urea              = 3
  integer, parameter :: ifert_N_no3               = 4
  integer, parameter :: ifert_N_nh4_band          = 5
  integer, parameter :: ifert_N_nh3_band          = 6
  integer, parameter :: ifert_N_urea_band         = 7
  integer, parameter :: ifert_N_no3_band          = 8
  integer, parameter :: ifert_P_Ca_H2PO4_2_soil   = 9
  integer, parameter :: ifert_P_Ca_H2PO4_2_band = 10
  integer, parameter :: ifert_P_apatite         = 11
  integer, parameter :: ifert_Ca_lime           = 12
  integer, parameter :: ifert_Ca_gypsum         = 13
  integer, parameter :: ifert_plant_resC        = 14
  integer, parameter :: ifert_plant_resN        = 15
  integer, parameter :: ifert_plant_resP        = 16
  integer, parameter :: ifert_plant_manuC       = 17
  integer, parameter :: ifert_plant_manuN       = 18
  integer, parameter :: ifert_plant_manuP       = 19
  integer, parameter :: ifert_PO4_soil          = 20
  integer, parameter :: ifert_PO4_band          = 21

  integer, parameter :: imanure_ruminant    = 1
  integer, parameter :: imanure_nonruminant = 2
  integer, parameter :: imanure_grazing     = 3

  integer, parameter :: iPlantRes_maize      = 1
  integer, parameter :: iPlantRes_wheat      = 2
  integer, parameter :: iPlantRes_soybean    = 3
  integer, parameter :: iPlantRes_oldStraw   = 4
  integer, parameter :: iPlantRes_Straw      = 5
  integer, parameter :: iPlantRes_compost    = 6
  integer, parameter :: iPlantRes_GreeManure = 7
  integer, parameter :: iPlantRes_simple     = 10
! root order
  integer, parameter :: iroot_1st =1
  integer, parameter :: iroot_2nd =2
! root profile
  integer, parameter :: irootp_shallow=0
! mycorrhizae
  integer, parameter :: imycorr_none=1
  integer, parameter :: imycorr_arbu=2   !arbuscular

  integer, parameter :: ipltroot=1
  integer, parameter :: imycorrhz=2

  integer, parameter :: iliving_branch=0

  integer, parameter :: itrue=1
  integer, parameter :: ifalse=0

  integer, parameter :: jharvtyp_noaction  = 0
  integer, parameter :: jharvtyp_terminate = 1
  integer, parameter :: jharvtyp_tmareseed = 2

  integer, parameter :: iharvtyp_none    = 0
  integer, parameter :: iharvtyp_grain   = 1
  integer, parameter :: iharvtyp_allabvg = 2
  integer, parameter :: iharvtyp_pruning = 3
  integer, parameter :: iharvtyp_grazing = 4
  integer, parameter :: iharvtyp_fire    = 5
  integer, parameter :: iharvtyp_herbivo = 6
  
  integer, parameter :: iplt_annual    = 0
  integer, parameter :: iplt_perennial = 1

  integer, parameter :: iplt_bryophyte=0
  integer, parameter :: iplt_grasslike=1
  integer, parameter :: iplt_treelike =2
  
  integer, parameter :: igraintyp_abvgrnd=0
  integer, parameter :: igraintyp_blwgrnd=1

  integer, parameter :: ideterminate          = 0
  integer, parameter :: inondeterminate       = 1
  integer, parameter :: ipltcal_Planting      = 0
  integer, parameter :: ipltcal_Emerge        = 1  !germination, ->0.9
  integer, parameter :: ipltcal_InitFloral    = 2  !main shoot or parent shoot leaf production, ->1.9
  integer, parameter :: ipltcal_Jointing      = 3  !tiller production, ->2.9
  integer, parameter :: ipltcal_Elongation    = 4  !stem elongation and booting, ->4.9
  integer, parameter :: ipltcal_heading       = 5  !heading, ->5.9
  integer, parameter :: ipltcal_Anthesis      = 6  !anthesis, ->6.9
  integer, parameter :: ipltcal_BeginSeedFill = 7  !grain milk stage, -7.9
  integer, parameter :: ipltcal_SetSeedNumber = 8  !grain doug s1, -8.5
  integer, parameter :: ipltcal_SetSeedMass   = 9  !grain dough s2, ->9.0
  integer, parameter :: ipltcal_EndSeedFill   = 10 !ripening, 9.9

  integer, parameter :: iphotop_neutral=0
  integer, parameter :: iphotop_short  =1
  integer, parameter :: iphotop_long   =2

  integer, parameter :: iembryotyp_Bryophytes   =0   !small, non-vascular land plants, like moss
  integer, parameter :: iembryotyp_Pteridophytes=1   !a vascular plant with xylem and phloem that reproduces by means of spores
  integer, parameter :: iembroytyp_Gymnosperms  =2   !example: conifers, cycads, ginkgos, and gnetophytes, woody, perennial seed-producing plants, typically lacking the protective outer covering which surrounds the seeds in flowering plants
  integer, parameter :: iembroytyp_Monocots     =3   !example: corn, grasses, wheat,onion, flowering plants whose seeds contain only one embryonic leaf, or cotyledon, they are angiosperms
  integer, parameter :: iembroytyp_Eudicots     =4   !example: soybean, rose, flowering plants that have two seed leaves upon germination, they are angiosperms
  
  integer, parameter :: isnowcept_bryophyte = 0
  integer, parameter :: isnowcept_grasses   = 1
  integer, parameter :: isnowcept_shrub     = 2
  integer, parameter :: isnowcept_decidtree = 3
  integer, parameter :: isnowcept_coniftree = 4

  integer, parameter :: iphenotyp_evgreen       = 0
  integer, parameter :: iphenotyp_coldecid      = 1
  integer, parameter :: iphenotyp_drouhtdecidu  = 2
  integer, parameter :: iphenotyp_coldroutdecid = 3

  integer, parameter :: in2fixtyp_none          = 0
  integer, parameter :: in2fixtyp_root_fast     = 1
  integer, parameter :: in2fixtyp_root_medium   = 2
  integer, parameter :: in2fixtyp_root_slow     = 3
  integer, parameter :: in2fixtyp_canopy_fast   = 4
  integer, parameter :: in2fixtyp_canopy_medium = 5
  integer, parameter :: in2fixtyp_canopy_slow   = 6

  integer, parameter :: ithermozone_arcboreal = 1
  integer, parameter :: ithermozone_cooltempr = 2
  integer, parameter :: ithermozone_warmtempr = 3
  integer, parameter :: ithermozone_subtropic = 4
  integer, parameter :: ithermozone_tropical  = 5

  integer, parameter :: itill_rmlitr = 21
  integer, parameter :: itill_fire   = 22
  integer, parameter :: iHarvst_pft  = 1    !indicator for pft specific harvest
  integer, parameter :: iHarvst_col  = 2    !indicator for col mean harvest
  public :: StriHarvtype,StrjHarvtype
contains
  
  function StrjHarvtype(jhavtype)result(ans)
  implicit none
  integer, intent(in) :: jhavtype
  character(len=48) :: ans

  select case (jhavtype)
  case (1)
    ans='terminate plant'
  case (2)
    ans= 'terminate&seed'
  case default
    ans='no touch'
  end select  
  end function StrjHarvtype
!------------------------------------------------------------------------------------------
  function StriHarvtype(iharvtyp)result(ans)
  implicit none
  integer, intent(in) :: iharvtyp
  character(len=48) :: ans

  select case (iharvtyp)
  case (1)
    ans='harvest grain'
  case (2)
    ans='harvest shoot'
  case(3)
    ans='shoot pruning'
  case (4)
    ans='shoot grazing'
  case (5)
    ans='fire burning'
  case (6)
    ans='animal grazing'
  case default
    ans= 'no touch'
  end select  
  end function StriHarvtype
end module ElmIDMod
