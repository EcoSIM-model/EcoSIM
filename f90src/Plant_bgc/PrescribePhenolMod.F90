module PrescribePhenolMod
  !
  ! !Description:
  ! Prescribe canopy structure and per-PFT root profiles for standalone SP and ATS.
  use abortutils, only: endrun, iulog
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use DebugToolMod, only: PrintInfo
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use EcoSimConst,   only: PICON2h
  use minimathmod,   only: AZMAX1
  use EcoSIMCtrlMod, only: etimer, ats_cpl_mode
  use EcoSIMCtrlDataType, only : ZEROS, current_month, day_of_month, total_days_in_month
  use EcoSIMHistMod, only : DATAPI
  use ElmIDMod
  use GridDataType
  use RootDataType
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantMorphologyAPIData, only : plt_morph
  use PlantBGCPars, only: Root2ndTipLen4uptk
  use ClimForcDataType
  use PlantTraitDataType
  use PlantMgmtDataType
  use CanopyDataType
  use CanopyRadDataType
  use PlantMathFuncMod, only: get_FDM, is_plant_woody_vascular
  use EcoSimConst,   only: ppmc
implicit none

  !pft information from CLM4.5
  !pft                                                 ra             rb           ztop(m)      zbot(m)
  ! 1. Needleleaf evergreen tree/NET Temperate         7.0            2.0          17           8.5
  ! 2. Needleleaf evergreen tree/NET Boreal            7.0            2.0          17           8.5
  ! 3. Needleleaf deciduous tree/NDT Boreal            7.0            2.0          14           7
  ! 4. Broadleave evergreen tree/BET tropical          7.0            1.0          35           1
  ! 5. Broadleave evergreen tree/BET temperate         7.0            1.0          35           1
  ! 6. Broadleave deciduous tree/BDT tropical          6.0            2.0          18           10
  ! 7. Broadleave deciduous tree/BDT temperate         6.0            2.0          20           11.5
  ! 8. Broadleave deciduous tree/BDT Boreal            6.0            2.0          20           11.5
  ! 9. Broadleave evergreen shrub/BES temperate        7.0            1.5          0.5          0.1
  !10. Broadleave deciduous shrub/BDS temperate        7.0            1.5          0.5          0.1
  !11. Broadleave deciduous shrub/BDS Boreal           7.0            1.5          0.5          0.1
  !12. C3 Arctic grass                                11.0            2.0          0.5          0.01
  !13. C3 Grass                                       11.0            2.0          0.5          0.01
  !14. C4 grass                                       11.0            2.0          0.5          0.01
  !15. C3 Crop Rainfed                                6.0             3.0          0.5          0.01
  !C3 Crop irrigated                              6.0             3.0          0.5          0.01
  !Corn Rainfed                                   6.0             3.0          0.5          0.01
  !Corn irrigated                                 6.0             3.0          0.5          0.01
  !Rainfed temperate cereals                      6.0             3.0          0.5          0.01
  !irrigated temperate cereals                    6.0             3.0          0.5          0.01
  !Rainfed winter cereals                         6.0             3.0          0.5          0.01
  !irrigated winter cereals                       6.0             3.0          0.5          0.01
  !Rainfed soybean                                6.0             3.0          0.5          0.01
  !irrigated soybean                              6.0             3.0          0.5          0.01
  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  public :: GetRootProfile
  public :: SetCanopyProfile
  public :: PrescribePhenologyInterp
  contains

  subroutine GetRootProfile(ra,rb,NL,ZH,RootProfl)
  implicit none
  real(r8), intent(in) :: ra
  real(r8), intent(in) :: rb
  integer,  intent(in) :: NL
  real(r8), intent(in) :: Zh(0:NL)
  real(r8), intent(out) :: RootProfl(1:NL)
  real(r8) :: rz(0:NL-1)
  integer :: L

  DO L=0,NL
    rz(L)=exp(-ra*zh(L))+exp(-rb*zh(L))
  enddo
  DO L=1,NL-1
    rootProfl(L)=0.5_r8*(rz(L-1)-rz(L))
  ENDDO
  rootProfl(NL)=0.5_r8*rz(NL-1)

  end subroutine GetRootProfile

!------------------------------------------------------------------------

  subroutine get_LAI_profile(zc,NZ)
  !Description:
  !Distribute leaf area in through the canopy
  !
  implicit none
  real(r8), intent(in) :: ZC  !plant height
  integer, intent(in) :: NZ   !plant index

  end subroutine get_LAI_profile

!------------------------------------------------------------------------

  subroutine get_StemArea_profile(NZ)
  !
  !Description:
  !Distribute steam area in through the canopy
  !
  implicit none
  integer, intent(in) :: NZ   !plant index


  end subroutine get_StemArea_profile

!------------------------------------------------------------------------

  function TreeStem_diameter_taperEq(h,p_stump_rate,p_stemp_rate,top_ht_rate,h_dbh,dbh,h_cb,h_top)result(d)

  ! Estimate diameter at height h using segment-specific taper rates as in Larsen (2017)
  implicit none

  real(r8), intent(in) :: h             ! - height at which to estimate diameter
  real(r8), intent(in) :: p_stump_rate  ! - stump taper rate [unitless]
  real(r8), intent(in) :: p_stemp_rate  ! - stem taper rate  [uitless]
  real(r8), intent(in) :: top_ht_rate   ! - crown base height/crown length [unitless]
  real(r8), intent(in) :: h_dbh         ! - breast height [m]
  real(r8), intent(in) :: dbh           ! - diameter at breast height [m]
  real(r8), intent(in) :: h_cb          ! - height to crown base [m]
  real(r8), intent(in) :: h_top         ! - total tree height (real)
  !Output
  real(r8) ::   d                       ! - estimated diameter at height h [m]

  if(h<=h_dbh)then
    !below breast height
    d=dbh+p_stump_rate*(h-h_dbh)
  elseif(h<=h_cb)then
    !below crown base
    d=dbh+p_stemp_rate*(h-h_dbh)
  else
    !above crown base
    d=(h_top-h)*top_ht_rate
  endif

  end function TreeStem_diameter_taperEq
!------------------------------------------------------------------------

  subroutine SetCanopyProfile(I,J,LeafAreaZsecLive_lpft,StemAreaZsecLive_lpft)
  !
  !Description:
  !distribute leaf area and stem area within the crown defined by top and bottom heights.
  !For simplicity, unifrom distribution is adopted
  !
  implicit none
  integer , intent(in) :: I,J
  real(r8),intent(out) :: LeafAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)     !leaf area in different zenith sectors in different canopy layers
  real(r8),intent(out) :: StemAreaZsecLive_lpft(NumLeafInclinationClasses1,NumCanopyLayers1,JP1)     !stem area in different zenith sectors in different canopy layers
  real(r8) :: dangle
  integer  :: N,L,NZ,NB

  associate(                                                        &
    CanopyStemSurfAreaZ_pft  => plt_morph%CanopyStemSurfAreaZ_pft  ,& !input  :plant canopy layer stem area, [m2 d-2]
    CanopyLeafAreaZ_pft      => plt_morph%CanopyLeafAreaZ_pft      ,& !inoput :total leaf area, [m2 d-2]
    LeafStalkAreaAct_pft     => plt_morph%LeafStalkAreaAct_pft     ,& !inoput :plant leaf+stem/stalk area, [m2 d-2]
    LeafStalkAreaAll_col     => plt_morph%LeafStalkAreaAll_col     ,& !inoput :stalk area of combined, each PFT canopy,[m^2 d-2]
    StemArea_col             => plt_morph%StemArea_col             ,& !input  :grid canopy stem area, [m2 d-2]
    LeafAngleClass_pft       => plt_morph%LeafAngleClass_pft       ,& !input  :fractionction of leaves in different angle classes, [-]
    CanopyLeafArea_col       => plt_morph%CanopyLeafArea_col       ,& !input  :grid canopy leaf area, [m2 d-2]
    LeafAreaZsec_brch        => plt_morph%LeafAreaZsec_brch        ,& !input  :leaf surface area, [m2 d-2]
    NP                       => plt_site%NP                         & !input  :current number of plant species,[-]
  )


  StemAreaZsecLive_lpft  = 0._r8
  LeafAreaZsecLive_lpft  = 0._r8
  dangle            = PICON2h/real(NumLeafInclinationClasses1,r8)         !the angle section width
  LeafStalkAreaAll_col = StemArea_col+CanopyLeafArea_col
  DO NZ=1,NP
    LeafStalkAreaAct_pft(NZ)=0._r8
    NB=1
    DO L=1,NumCanopyLayers1
      DO N=1,NumLeafInclinationClasses1
        LeafAreaZsecLive_lpft(N,L,NZ)=LeafAngleClass_pft(N,NZ)*CanopyLeafAreaZ_pft(L,NZ)/real(NumOfLeafAzimuthSectors1,r8)
      ENDDO
      StemAreaZsecLive_lpft(:,L,NZ)=0._r8
      StemAreaZsecLive_lpft(NumLeafInclinationClasses1,L,NZ)=CanopyStemSurfAreaZ_pft(L,NZ)/real(NumLeafInclinationClasses1,kind=r8)
      LeafStalkAreaAct_pft(NZ)=LeafStalkAreaAct_pft(NZ)+CanopyLeafAreaZ_pft(L,NZ)+CanopyStemSurfAreaZ_pft(L,NZ)

      !Assuming uniform azimuth desitribution for leaves
      DO N=1,NumLeafInclinationClasses1
        LeafAreaZsec_brch(N,L,L,NB,NZ)=LeafAreaZsecLive_lpft(N,L,NZ)
      ENDDO
    ENDDO

  ENDDO

  end associate
  end subroutine SetCanopyProfile

!------------------------------------------------------------------------------------------
  subroutine PrescribePhenologyInterp(I, NHW, NHE, NVN, NVS)
  ! Interpolate prescribed canopy structure and assign each named PFT's root profile.
  implicit none
  integer, intent(in) :: I
  Integer, intent(in) ::  NHW, NHE, NVN, NVS
  real(r8) :: t
  integer :: it(2)
  integer :: months(2)
  integer :: kmo,dofmon,ndaysmon,NY,NX,NZ,L
  real(r8) :: timwt(2)
  real(r8) :: PopulationDensity !prescribed population per ground area, [plants m-2]
  real(r8) :: RootCPerPlant    !prescribed live fine root C per plant, [gC]
  real(r8) :: NumPrimRootAxes  !primary root axes of the PFT, [d-2]
  real(r8) :: RootCBelow       !root C per plant at or below the current layer, [gC]
  real(r8) :: ZL1(0:NumCanopyLayers)
  real(r8) :: AreaInterval,AreaL
  real(r8) :: ARX  !interval canopy area: leaf+stem
  real(r8) :: DZL  !canopy interval height
  REAL(R8) :: PerPlantFineRootC_vr(1:JZ)
  REAL(R8) :: PerPlantRootLen_vr(1:JZ)

  if(ats_cpl_mode)then
    !ATS supplies daily LAI, SAI and height; keep uniform leaf inclination
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        DO NZ=1,NP_col(NY,NX)
          LeafAngleClass_pft(:,NZ,NY,NX) = 1._r8/real(NumLeafInclinationClasses,kind=r8)
        ENDDO
      ENDDO
    ENDDO
    dofmon = day_of_month
    kmo = current_month
    ndaysmon = total_days_in_month
  else
    !monthly LAI, SAI and height are read from sp_phenol_in (Read4SatePhenolMod);
    !leaf inclination is the PFT trait from pft_file_in
    ndaysmon = etimer%get_curr_mon_days()
    dofmon   = etimer%get_curr_dom()
    kmo      = etimer%get_curr_mon()
  end if

  t = (dofmon-0.5_r8) / ndaysmon
  it(1) = t + 0.5_r8
  it(2) = it(1) + 1
  months(1) = kmo + it(1) - 1
  months(2) = kmo + it(2) - 1
  if (months(1) <  1) months(1) = 12
  if (months(2) > 12) months(2) = 1

  timwt(1) = (it(1)+0.5_r8) - t
  timwt(2) = 1._r8-timwt(1)

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      ! SP prescribes secondary root C without simulating structural axes.
      ! Clear retained/empty slots, other elements and mycorrhizal groups too.
      RootMyco2ndStrutElms_pvr(:,:,:,1:NP0_col(NY,NX),NY,NX) = 0._r8
      Root2ndXNumL_rpvr(:,:,1:NP0_col(NY,NX),NY,NX)          = 0._r8
      Root2ndEffLen4uptk_rpvr(:,:,1:NP0_col(NY,NX),NY,NX)    = Root2ndTipLen4uptk
      CanopyLeafArea_col(NY,NX)                              = 0._r8
      StemArea_col(NY,NX)                                    = 0._r8
      !CanopyHeight_col(NY,NX)                               = 0.0_r8

      DO NZ=1,NP_col(NY,NX)
        if(.not.ats_cpl_mode)then
          ! Native daily canopy areas are cell totals. ATS already converts its inputs.
          tlai_day_pft(NZ,NY,NX) = &
            (timwt(1)*tlai_mon_pft(months(1),NZ,NY,NX)+timwt(2)*tlai_mon_pft(months(2),NZ,NY,NX)) &
            *AREA_3D(3,NU_col(NY,NX),NY,NX)
          tsai_day_pft(NZ,NY,NX) = &
            (timwt(1)*tsai_mon_pft(months(1),NZ,NY,NX)+timwt(2)*tsai_mon_pft(months(2),NZ,NY,NX)) &
            *AREA_3D(3,NU_col(NY,NX),NY,NX)
          CanopyHeightLive_pft(NZ,NY,NX) = timwt(1)*height_top_mon_pft(months(1),NZ,NY,NX)+timwt(2)*height_top_mon_pft(months(2),NZ,NY,NX)
        endif
        
        !empty ATS PFT slot: no plants, so it stays inactive in the uptake model
        if(ats_cpl_mode .and. DATAPI(NZ,NY,NX).EQ.0)then
          PlantPopuLive_pft(NZ,NY,NX)=0._r8
        elseif(.not.ats_cpl_mode)then
          PlantPopuLive_pft(NZ,NY,NX)=sp_popdens_pft(NZ,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
        endif

        CanopyLeafArea_col(NY,NX)  = CanopyLeafArea_col(NY,NX)+tlai_day_pft(NZ,NY,NX)
        StemArea_col(NY,NX)        = StemArea_col(NY,NX)+tsai_day_pft(NZ,NY,NX)
        CanopyHeight_col(NY,NX)    = AMAX1(CanopyHeight_col(NY,NX),CanopyHeightLive_pft(NZ,NY,NX))
        CanopyLeafAreaZ_pft(1:NumCanopyLayers,NZ,NY,NX)=tlai_day_pft(NZ,NY,NX)/real(NumCanopyLayers,kind=r8)
        CanopyStemSurfAreaZ_pft(1:NumCanopyLayers,NZ,NY,NX)=tsai_day_pft(NZ,NY,NX)/real(NumCanopyLayers,kind=r8)

        !per-PFT leaf and stem areas, and the canopy carbon they imply. The carbon
        !sets the plant water storage capacity in ROOTUPTAKES (CanopyMassC4H2OStorage)
        !  leaf + sheath C = leaf area / specific leaf area SLA1 [m2 gC-1]
        !  sapwood C = stem area x sapwood thickness (2.5 mm, as assumed in
        !             UptakesMod) / specific stalk volume [m3 gC-1]

        CanopyLeafArea_pft(NZ,NY,NX)     = tlai_day_pft(NZ,NY,NX)
        CanopyStemSurfArea_pft(NZ,NY,NX) = tsai_day_pft(NZ,NY,NX)
        if(SLA1_pft(NZ,NY,NX).GT.0._r8)then
          CanopyLeafSheathC_pft(NZ,NY,NX) = CanopyLeafArea_pft(NZ,NY,NX)/SLA1_pft(NZ,NY,NX)
        else
          CanopyLeafSheathC_pft(NZ,NY,NX) = 0._r8
        endif
        if(StemSpecVolume_pft(NZ,NY,NX).GT.0._r8)then
          CanopySapwoodC_pft(NZ,NY,NX) = CanopyStemSurfArea_pft(NZ,NY,NX)*2.5E-03_r8/StemSpecVolume_pft(NZ,NY,NX)
        else
          CanopySapwoodC_pft(NZ,NY,NX) = 0._r8
        endif
        !plant tissue water not yet set: start it in equilibrium with the leaf
        !carbon and the canopy water potential, as StartqMod does for the full model
        if(CanopyBiomWater_pft(NZ,NY,NX).LE.0._r8 .and. CanopyLeafSheathC_pft(NZ,NY,NX).GT.0._r8)then
          CanopyBiomWater_pft(NZ,NY,NX) = ppmc*CanopyLeafSheathC_pft(NZ,NY,NX)/get_FDM(PSICanopy_pft(NZ,NY,NX))
        endif
      ENDDO

      !set vertical desitribution of LAI and
      CanopyLeafAareZ_col(1:NumCanopyLayers,NY,NX)=CanopyLeafArea_col(NY,NX)/real(NumCanopyLayers,kind=r8)
      CanopyStemAareZ_col(1:NumCanopyLayers,NY,NX)=StemArea_col(NY,NX)/real(NumCanopyLayers,kind=r8)

      !divide canopy height
      CanopyHeightZ_col(NumCanopyLayers,NY,NX) = CanopyHeight_col(NY,NX)+0.01_r8
      ZL1(NumCanopyLayers)                     = CanopyHeightZ_col(NumCanopyLayers,NY,NX)
      ZL1(0)                                   = 0.0_r8

      !for simplicity, right now unifrom division is used
      !divide total are into NumCanopyLayers1, from top to bottom
      AreaInterval = (CanopyLeafArea_col(NY,NX)+StemArea_col(NY,NX))/NumCanopyLayers
      DZL          = CanopyHeightZ_col(NumCanopyLayers,NY,NX)/real(NumCanopyLayers,kind=r8)

      IF(AreaInterval.GT.ZEROS(NY,NX))THEN
         DO L=NumCanopyLayers,1,-1
           CanopyHeightZ_col(L-1,NY,NX)=CanopyHeightZ_col(L,NY,NX)-DZL
         ENDDO
      ENDIF

      DO NZ=1,NP_col(NY,NX)
        !empty ATS PFT slot: zero its roots rather than inheriting the
        !profile left in PerPlantRoot*_vr by the previous slot
        if(ats_cpl_mode .and. DATAPI(NZ,NY,NX).EQ.0)then
          DO L=NU_col(NY,NX),NL_col(NY,NX)
            RootTotLenPerPlant_pvr(ipltroot,L,NZ,NY,NX)    = 0._r8
            RootLenDensPerPlant_pvr(ipltroot,L,NZ,NY,NX)   = 0._r8
            PopuRootMycoC_pvr(ipltroot,L,NZ,NY,NX)         = 0._r8
            RootAbsorbLenPerPlant_pvr(ipltroot,L,NZ,NY,NX) = 0._r8
          ENDDO
          cycle
        endif
        ! Validate the PFT's prescribed axis length before using it as a divisor.
        if(.not.ieee_is_finite(MeanSecondaryAxisLength_pft(NZ,NY,NX)))then
          write(iulog,*)'Invalid MFRUTL for PFT/cell:',NZ,NY,NX,MeanSecondaryAxisLength_pft(NZ,NY,NX)
          call endrun('prescribed PFT MFRUTL must be finite and at least Root2ndTipLen4uptk',__LINE__)
        endif

        if(MeanSecondaryAxisLength_pft(NZ,NY,NX)<Root2ndTipLen4uptk)then
          write(iulog,*)'MFRUTL and minimum length:',MeanSecondaryAxisLength_pft(NZ,NY,NX),Root2ndTipLen4uptk
          call endrun('prescribed PFT MFRUTL must be finite and at least Root2ndTipLen4uptk',__LINE__)
        endif

        ! Each named PFT supplies its own valid root type; no column fallback.
        if(ats_cpl_mode)then
          call SetRootProfileZ(irootType_pft(NZ,NY,NX),NL_col(NY,NX), &
            CumDepz2LayBottom_vr(1:NL_col(NY,NX),NY,NX),PerPlantFineRootC_vr, &
            PerPlantRootLen_vr,PopulationDensity)
          PlantPopuLive_pft(NZ,NY,NX)=PopulationDensity
        else
          call SetRootProfileZ(irootType_pft(NZ,NY,NX),NL_col(NY,NX), &
            CumDepz2LayBottom_vr(1:NL_col(NY,NX),NY,NX),PerPlantFineRootC_vr, &
            PerPlantRootLen_vr,PopulationDensity,sp_popdens_pft(NZ,NY,NX))
        endif

        !number of primary root axes, by the pipe-model allometry of the full plant
        !model (GrosubsMod): max(1, root C per plant**b) x population. Each layer
        !gets the axes that reach it: the total scaled by the share of root C at
        !or below the layer.
        RootCPerPlant = sum(PerPlantFineRootC_vr(NU_col(NY,NX):NL_col(NY,NX)))

        IF(iPlantPhenolPattern_pft(NZ,NY,NX).EQ.iplt_annual)THEN
          NumPrimRootAxes = AMAX1(1.0_r8,RootCPerPlant**0.833_r8)*PlantPopuLive_pft(NZ,NY,NX)
        elseif(is_plant_woody_vascular(iPlantRootProfile_pft(NZ,NY,NX),iPlant2ndGrothPattern_pft(NZ,NY,NX)))then
          NumPrimRootAxes = AMAX1(1.0_r8,RootCPerPlant**0.7143_r8)*PlantPopuLive_pft(NZ,NY,NX)
        else
          NumPrimRootAxes = AMAX1(1.0_r8,RootCPerPlant**1.053_r8)*PlantPopuLive_pft(NZ,NY,NX)
        endif

        RootCBelow = RootCPerPlant
        DO L=NU_col(NY,NX),NL_col(NY,NX)
          RootTotLenPerPlant_pvr(ipltroot,L,NZ,NY,NX)    = PerPlantRootLen_vr(L)
          RootLenDensPerPlant_pvr(ipltroot,L,NZ,NY,NX)   = PerPlantRootLen_vr(L)/DLYR_3D(3,L,NY,NX)
          PopuRootMycoC_pvr(ipltroot,L,NZ,NY,NX)         = PerPlantFineRootC_vr(L)*PopulationDensity

          ! Map all prescribed live fine-root C to secondary structural C.
          ! gC plant-1 x plants m-2 x m2 cell-1 = gC cell-1.
          RootMyco2ndStrutElms_pvr(ielmc,ipltroot,L,NZ,NY,NX) = &
            PerPlantFineRootC_vr(L)*PopulationDensity*AREA_3D(3,NU_col(NY,NX),NY,NX)
          RootAbsorbLenPerPlant_pvr(ipltroot,L,NZ,NY,NX) = PerPlantRootLen_vr(L)

          !primary axes reaching layer L (was a filler of 2)
          if(RootCPerPlant.GT.0._r8)then
            Root1stXNumL_pvr(L,NZ,NY,NX) = AMAX1(1.0_r8,NumPrimRootAxes*RootCBelow/RootCPerPlant)
          else
            Root1stXNumL_pvr(L,NZ,NY,NX) = 2._r8
          endif
          RootCBelow = AZMAX1(RootCBelow-PerPlantFineRootC_vr(L))
          ! Prescribed length [m plant-1] x population [plants m-2] x cell area
          ! [m2 cell-1] / MFRUTL [m axis-1] gives secondary axes per cell.
          Root2ndXNumL_rpvr(ipltroot,L,NZ,NY,NX) = &
            PerPlantRootLen_vr(L)*PopulationDensity*AREA_3D(3,NU_col(NY,NX),NY,NX) &
            /MeanSecondaryAxisLength_pft(NZ,NY,NX)
          Root2ndEffLen4uptk_rpvr(ipltroot,L,NZ,NY,NX)=MeanSecondaryAxisLength_pft(NZ,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  end subroutine PrescribePhenologyInterp
!------------------------------------------------------------------------------------------
  subroutine SetRootProfileZ(irootType,NL,tmp_cdepthz,tmp_rootc,tmp_rootl,PlantPop,PrescribedPop)
  !
  !Description:
  ! Calculate per-individual live fine-root C and length by layer for a Jackson biome type.
  !Reference: Jackson et al. (1997), A global budget for fine root biomass, surface area, and nutrient contents, PNAS.
  !|Biome | Total fine root biomass (kg m$^{-2}$) | Live fine root biomass (kg m$^{-2}$) | Live fine root length (km m$^{-2}$) | Live fine root area index (m2 m$^{-2}$) |$\beta$|
  !|----------|----------|----------|----------|----------|----------|
  !|1. Boreal forest|0.60 | 0.23 | 2.6 | 4.6 | 0.943|
  !|2. Desert bushes | 0.27|0.13|4.0 | 5.5 | 0.970|
  !|3. Sclerophyllous shrubs and tress |0.52|0.28|8.4|11.6|0.950|
  !|4. Temperate confierous forest|0.82|0.50|6.1|11.0|0.980|
  !|5. Temperature deciduous forest|0.78|0.44|5.4|9.8|0.967|
  !|6. Temperate grassland|1.51|0.95|112|79.1|0.943|
  !|7. Tropical deciduous forest|0.57|0.28|3.5|6.3|0.982|
  !|8. Tropical evergreen forest|0.57|0.33|4.1|7.4|0.972|
  !|9. Tropical grassland/savanna|0.99|0.51|60.4|42.5|0.972|
  !|10. Tundra|0.96|0.34|7.4|5.2|0.909|

  implicit none
  integer, intent(in) :: iRootType        !Jackson biome index, 1:10, [-]
  integer, intent(in) :: NL               !number of soil layers to populate, [-]
  real(r8),intent(in) :: tmp_cdepthz(:)   !layer-bottom depths, size >= NL, [m]
  real(r8),intent(out) :: tmp_rootc(:)    !live fine-root C per layer and individual, size >= NL, [gC plant-1]
  real(r8),intent(out) :: tmp_rootl(:)    !live fine-root length per layer and individual, size >= NL, [m plant-1]
  real(r8),intent(out) :: PlantPop        !population density used for the per-plant conversion, [plants m-2]
  real(r8),intent(in),optional :: PrescribedPop !standalone SP population density from file, [plants m-2]
  character(len=*), parameter :: subname='SetRootProfileZ'
  real(r8), parameter :: beta(10)=real((/0.943,0.970,0.950,0.980,0.967,0.943,0.982,0.972,0.972,0.909/),kind=r8)
  ! Live fine-root C [gC m-2]: Jackson et al. Table 1 live biomass [kg dry mass m-2]
  ! x 0.488 gC/g dry mass x 1000 g/kg. Dead roots are excluded from prescribed plant C.
  real(r8), parameter :: livefrootC(10) = &
    (/0.23_r8,0.13_r8,0.28_r8,0.50_r8,0.44_r8, &
      0.95_r8,0.28_r8,0.33_r8,0.51_r8,0.34_r8/)*0.488_r8*1.e3_r8
  real(r8), parameter :: frootLen(10)=(/2.6,4.0,8.4,6.1,5.4,112.,3.5,4.1,60.4,7.4/)*1.e3_r8  !live fine root length
  real(r8), parameter :: PltPopDef(10)=(/0.6,1.,1.0,0.6,0.6,40.,0.6,0.6,40.,40./) !default plant population [1/m2]
  integer :: L
  real(r8) :: CumRootFrac_vr(JZ)  !Cumfraction of root in soil layers
  real(r8) :: RootFrac_vr(JZ)     !fraction of root in soil layers

  call PrintInfo('beg '//subname)
  if(irootType<1 .or. irootType>size(beta))then
    write(iulog,*)subname//': invalid IROOTJK=',irootType
    call endrun('prescribed PFT IROOTJK must be between 1 and 10 in '//subname,__LINE__)
  endif
  if(NL<1 .or. NL>size(tmp_cdepthz) .or. NL>size(tmp_rootc) .or. NL>size(tmp_rootl))then
    write(iulog,*)subname//': NL and depth/C/length sizes=',NL,size(tmp_cdepthz),size(tmp_rootc),size(tmp_rootl)
    call endrun('invalid prescribed root-profile dimensions in '//subname,__LINE__)
  endif
  tmp_rootc=0._r8
  tmp_rootl=0._r8
  PlantPop = PltPopDef(irootType) !ATS retains its existing biome defaults.
  if(present(PrescribedPop))PlantPop=PrescribedPop
  if(.not.ieee_is_finite(PlantPop))then
    call endrun('prescribed population density must be finite and positive in '//subname,__LINE__)
  endif
  if(PlantPop<=0._r8)then
    call endrun('prescribed population density must be finite and positive in '//subname,__LINE__)
  endif

  CumRootFrac_vr(1)=1._r8-beta(irootType)**(tmp_cdepthz(1)*100._r8)
  RootFrac_vr(1)=CumRootFrac_vr(1)
  tmp_rootl(1)=frootLen(irootType)*RootFrac_vr(1)/PlantPop
  tmp_rootc(1)=RootFrac_vr(1)*livefrootC(irootType)/PlantPop

  DO L=2,NL
    CumRootFrac_vr(L)=1._r8-beta(irootType)**(tmp_cdepthz(L)*100._r8)
    RootFrac_vr(L) = CumRootFrac_vr(L)-CumRootFrac_vr(L-1)
    tmp_rootl(L)=frootLen(irootType)*RootFrac_vr(L)/PlantPop
    tmp_rootc(L)=RootFrac_vr(L)*livefrootC(irootType)/PlantPop
  enddo

  call PrintInfo('end '//subname)
  end subroutine SetRootProfileZ
end module PrescribePhenolMod
