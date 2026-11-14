module PlantMathFuncMod
  ! Description:
  ! Provide shared plant process functions, including canopy thermal and radiative geometry.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use abortutils,    only: endrun, iulog
  use PlantAPICommonData
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantEnergyWaterAPIData, only : plt_ew
  use DebugToolMod
  use EcoSimConst
  use MiniMathMod
  use MiniFuncMod
  use ElmIDMod
  use PlantBGCPars, only: CanopySnowCapacityPerArea,CanopySnowAreaCoeff
implicit none
  character(len=*), parameter, private :: mod_filename=&
  __FILE__

  type, public  :: PlantSoluteUptakeConfig_type
    real(r8) :: SolAdvFlx  
    real(r8) :: SolDifusFlx    
    real(r8) :: UptakePerPlantRateMax
    real(r8) :: O2Stress        
    real(r8) :: PlantPopulation 
    real(r8) :: CAvailStress    
    real(r8) :: SoluteMassMax   
    real(r8) :: SoluteConc      
    real(r8) :: SoluteKM      
    real(r8) :: SoluteConcMin       
  end type PlantSoluteUptakeConfig_type
contains

!------------------------------------------------------------------------------------------
  subroutine CanopySnowOptics(SnowType,PhysicalArea,SnowVolume,SnowFraction,EffectiveArea)
    ! Description:
    ! Derive snow cover and radiative area from retained SWE and current exposed
    ! live/dead area. Buried material retains its snow but has no radiative area.
    integer, intent(in) :: SnowType ! isnowcept_* interception class [0:4]
    real(r8), intent(in) :: PhysicalArea ! exposed leaf plus live/dead stem area [m2/cell]
    real(r8), intent(in) :: SnowVolume ! retained intercepted snow SWE [m3/cell]
    real(r8), intent(out) :: SnowFraction ! snow-covered fraction of exposed material [-]
    real(r8), intent(out) :: EffectiveArea ! physical plus snow-added radiative area [m2/cell]
    character(len=*), parameter :: subname='CanopySnowOptics'

    call PrintInfo('beg '//subname)
    if(SnowType<LBOUND(CanopySnowAreaCoeff,1).or.SnowType>UBOUND(CanopySnowAreaCoeff,1) &
      .or.any(.not.ieee_is_finite([PhysicalArea,SnowVolume])).or.MIN(PhysicalArea,SnowVolume)<0._r8)then
      write(iulog,*)subname,': invalid class, exposed area or snow SWE: ',SnowType,PhysicalArea,SnowVolume
      call endrun('Invalid canopy snow geometry',__LINE__)
    endif
    SnowFraction=0._r8
    EffectiveArea=0._r8
    if(PhysicalArea>0._r8)then
      if(SnowVolume>0._r8)then
        SnowFraction=1._r8
        if(SnowVolume<CanopySnowCapacityPerArea(SnowType)*PhysicalArea) &
          SnowFraction=(SnowVolume/PhysicalArea/CanopySnowCapacityPerArea(SnowType))**0.15_r8
      endif
      EffectiveArea=PhysicalArea+CanopySnowAreaCoeff(SnowType)*SnowVolume
    endif
    call PrintInfo('end '//subname)
  end subroutine CanopySnowOptics

!------------------------------------------------------------------------------------------
  subroutine CanopyRadiationFractions(LiveArea,DeadArea,SnowVolume,SnowType,GroundArea, &
    CanopyFraction,LiveFraction,GroundFraction)
    ! Description:
    ! Partition the column longwave interception using the same PFT snow-added
    ! areas as shortwave optics. All inputs refer to the current exposed canopy.
    real(r8), intent(in) :: LiveArea(:),DeadArea(:) ! exposed physical live/dead area [m2/cell, PFT]
    real(r8), intent(in) :: SnowVolume(:) ! retained canopy snow SWE [m3/cell, PFT]
    integer, intent(in) :: SnowType(:) ! isnowcept_* interception class [PFT]
    real(r8), intent(in) :: GroundArea ! column ground area [m2/cell]
    real(r8), intent(out) :: CanopyFraction(:) ! total interception fraction [PFT]
    real(r8), intent(out) :: LiveFraction(:) ! live share of total interception [PFT]
    real(r8), intent(out) :: GroundFraction ! transmitted fraction [-]
    real(r8) :: EffectiveArea(SIZE(LiveArea)),SnowFraction(SIZE(LiveArea))
    real(r8) :: TotalArea,MeanSnowFraction,Intercepted,PhysicalArea
    integer :: NZ
    character(len=*), parameter :: subname='CanopyRadiationFractions'

    call PrintInfo('beg '//subname)
    if(.not.ieee_is_finite(GroundArea).or.GroundArea<=0._r8 &
      .or.SIZE(DeadArea)/=SIZE(LiveArea).or.SIZE(SnowVolume)/=SIZE(LiveArea) &
      .or.SIZE(SnowType)/=SIZE(LiveArea).or.SIZE(CanopyFraction)/=SIZE(LiveArea) &
      .or.SIZE(LiveFraction)/=SIZE(LiveArea))then
      write(iulog,*)subname,': invalid ground area or PFT dimensions: ',GroundArea
      call endrun('Invalid canopy radiation partition',__LINE__)
    endif
    if(any(.not.ieee_is_finite(LiveArea)).or.any(.not.ieee_is_finite(DeadArea)) &
      .or.any(LiveArea<0._r8).or.any(DeadArea<0._r8))then
      write(iulog,*)subname,': invalid exposed live or dead areas'
      call endrun('Invalid physical canopy area',__LINE__)
    endif
    CanopyFraction=0._r8
    LiveFraction=0._r8
    GroundFraction=1._r8
    do NZ=1,SIZE(LiveArea)
      call CanopySnowOptics(SnowType(NZ),LiveArea(NZ)+DeadArea(NZ),SnowVolume(NZ), &
        SnowFraction(NZ),EffectiveArea(NZ))
    enddo
    TotalArea=SUM(EffectiveArea)
    if(TotalArea>0._r8)then
      MeanSnowFraction=SUM(SnowFraction*EffectiveArea)/TotalArea
      Intercepted=1._r8-EXP(-(0.65_r8+0.15_r8*MeanSnowFraction)*TotalArea/GroundArea)
      CanopyFraction=Intercepted*EffectiveArea/TotalArea
      GroundFraction=1._r8-Intercepted
      do NZ=1,SIZE(LiveArea)
        PhysicalArea=LiveArea(NZ)+DeadArea(NZ)
        if(PhysicalArea>0._r8)LiveFraction(NZ)=CanopyFraction(NZ)*LiveArea(NZ)/PhysicalArea
      enddo
    endif
    call PrintInfo('end '//subname)
  end subroutine CanopyRadiationFractions

!------------------------------------------------------------------------------------------
  function CanopyStemThermalCarbon(Carbon,SurfaceArea,SpecificVolume) result(ActiveCarbon)
    ! Description:
    ! Limit live or dead stem C to the same thermally represented outer shell.
    ! Total material area retains this capacity when snow hides the exchange surface.
    real(r8), intent(in) :: Carbon ! available sapwood or standing-dead C [gC/cell]
    real(r8), intent(in) :: SurfaceArea ! total stem surface area, including buried material [m2/cell]
    real(r8), intent(in) :: SpecificVolume ! stem volume per carbon mass [m3/gC]
    real(r8) :: ActiveCarbon ! C represented in the canopy thermal store [gC/cell]
    character(len=*), parameter :: subname='CanopyStemThermalCarbon'
    real(r8), parameter :: CanopyStemShellDepth=2.5e-3_r8 ! shell thickness [m]
    call PrintInfo('beg '//subname)
    if(.not.ieee_is_finite(Carbon))then
      write(iulog,*)subname,': nonfinite stem carbon: ',Carbon
      call endrun('Invalid canopy stem carbon',__LINE__)
    endif
    ActiveCarbon=0._r8
    if(Carbon>0._r8)then
      if(.not.all(ieee_is_finite([Carbon,SurfaceArea,SpecificVolume])) &
        .or. SurfaceArea<0._r8 .or. SpecificVolume<=0._r8)then
        write(iulog,*)subname,': invalid stem C, area, specific volume: ',Carbon,SurfaceArea,SpecificVolume
        call endrun('Invalid canopy stem thermal geometry',__LINE__)
      endif
      ActiveCarbon=MIN(Carbon,CanopyStemShellDepth*SurfaceArea/SpecificVolume)
    endif
    call PrintInfo('end '//subname)
  end function CanopyStemThermalCarbon

!------------------------------------------------------------------------------------------

  function get_FDM(PSIOrgan,FDMP0)result(FDMP)
  !
  !compute the Ratio of leaf+sheath dry mass to symplasmic water (g g–1)
  !as a function of absolute value of leaf water potential (MPa)
  
  implicit none
  real(r8), intent(in) :: PSIOrgan !canopy water potential, MPa
  real(r8), optional, intent(out) ::  FDMP0
  
  real(r8) :: APSILT    !abosolute value of psi
  real(r8) :: FDMP      !=dry matter/water

  APSILT = ABS(PSIOrgan)
  FDMP   = 0.16_r8+0.10_r8*APSILT/(0.05_r8*APSILT+2.0_r8)
  if(present(FDMP0))FDMP0=0.16_r8

  end function get_FDM

!--------------------------------------------------------------------------------
  subroutine update_osmo_turg_pressure(PSIOrgan,CCPOLT,OSMO,TKP,PSIOsmo,PSITurg,FDMP1)
  !
  !DESCRIPTION
  !update the osmotic and turgor pressure of a plant organ
  implicit none
  real(r8), intent(in) :: PSIOrgan    !plant orgran pressure, [MPa]
  real(r8), intent(in) :: CCPOLT      !total organ non-structrual elemental concentration, g/g
  real(r8), intent(in) :: OSMO        !Organ osmotic potential when water potential = 0 MPa
  real(r8), intent(in) :: TKP         !organ temperature, Kelvin
  real(r8), intent(out) :: PSIOsmo    !osmotic pressure of the organ, MPa
  real(r8), intent(out) :: PSITurg    !turgor pressure of the organ, MPa
  real(r8), optional, intent(out) :: FDMP1

  real(r8) :: OSWT   !average molecular weight of CCPOLT 
  real(r8) :: FDMP   !Ratio of leaf+sheath dry mass to symplasmic water (g/g)
  real(r8) :: FDMP0  !FDMP at zero canopy water potential. 
  
  FDMP    = get_FDM(PSIOrgan,FDMP0)
  OSWT    = 36.0_r8+840.0_r8*AZMAX1(CCPOLT)
  PSIOsmo = FDMP*(OSMO/FDMP0-RGASC*TKP*CCPOLT/OSWT)
  PSITurg = AZMAX1(PSIOrgan-PSIOsmo)

  if(present(fdmp1))FDMP1=FDMP
  end subroutine update_osmo_turg_pressure
!--------------------------------------------------------------------------------

  function get_zero_turg_ccpolt(PSIOrg,OSMO,TKP)result(CCPOLT)
  implicit none
  real(r8), intent(in) :: PSIOrg      !plant orgran pressure, [MPa]
  real(r8), intent(in) :: OSMO        !Organ osmotic potential [MPa], when water potential = 0 MPa
  real(r8), intent(in) :: TKP         !organ temperature, [Kelvin]

  real(r8) :: CCPOLT      !total organ non-structrual elemental concentration, g/g
  !--- local variables ---
  real(r8) :: FDMP,FDMP0,ratio

  FDMP  = get_FDM(PSIOrg,FDMP0)

  ratio = RGASC*TKP/(OSMO/FDMP0-PSIOrg/FDMP)

  CCPOLT = 36._r8/(ratio-840._r8)
  end function get_zero_turg_ccpolt
!--------------------------------------------------------------------------------

  subroutine calc_seed_geometry(SeedCMass,rwidth2lenSeed,SeedVolumeMean,SeedLengthMean,SeedArea)
  !
  !DESCRIPTION
  !assuming the seed is spherical, compute its volume, diameter(=length), and surface area
  !     SeedVolumeMean,SeedLengthMean,SeedAreaMean=seed volume(m3),length(m),AREA_3D(m2)
  !     SeedCMass=seed C mass (g) from PFT file
  !

  implicit none
  real(r8), intent(in)  :: SeedCMass        !carbon mass per seed, gC/seed
  real(r8), intent(in)  :: rwidth2lenSeed   !seed width to length ratio
  real(r8), intent(out) :: SeedVolumeMean,SeedLengthMean,SeedArea
  real(r8) :: seedHalfLength
  real(r8), parameter :: pp=1.6075_r8

  SeedVolumeMean = SeedCMass*5.0E-06_r8
  seedHalfLength=(0.75_r8*SeedVolumeMean/PICON)**0.33_r8*rwidth2lenSeed**(-0.667_r8) !assume prolate
  SeedLengthMean = 2.0_r8*seedHalfLength  
  SeedArea       = 4.0_r8*PICON*seedHalfLength**2*((2._r8*rwidth2lenSeed**pp+rwidth2lenSeed**(2._r8*pp))/3._r8)**(1._r8/pp)
  
  end subroutine calc_seed_geometry
!--------------------------------------------------------------------------------
  pure function calc_root_grow_tempf(TKSO)result(fT_root)
  !
  !DESCRIPTION
  !compute the temperature dependence for plant root growth
  implicit none
  real(r8), intent(in) :: TKSO   !apparent temperature felt by the root
  real(r8) :: fT_root
  real(r8) :: RTK,STK,ACTV

  RTK=RGASC*TKSO
  STK=710.0_r8*TKSO
  ACTV=1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  FT_ROOT=EXP(25.229_r8-62500._r8/RTK)/ACTV

  end function calc_root_grow_tempf
!--------------------------------------------------------------------------------
  pure function calc_canopy_grow_tempf(TKGO)result(fT_canp)
  !
  !DESCRIPTION
  !compute the temperature dependence for plant canopy growth
  implicit none
  real(r8), intent(in) :: TKGO   !apparent temperature felt by the canopy
  real(r8) :: fT_canp
  real(r8) :: RTK,STK,ACTV

  RTK     = RGASC*TKGO
  STK     = 710.0_r8*TKGO
  ACTV    = 1+EXP((197500._r8-STK)/RTK)+EXP((STK-222500._r8)/RTK)
  FT_canp = EXP(25.229_r8-62500._r8/RTK)/ACTV

  end function calc_canopy_grow_tempf

!--------------------------------------------------------------------------------
  pure function calc_leave_grow_tempf(TKCO)result(TFNP)
  !DESCRIPTION
  !compute the temperature dependence for leave growth
  implicit none
  real(r8), intent(in) :: TKCO   !apparent temperature felt by the leave
  real(r8) :: TFNP
  real(r8) :: RTK,STK,ACTV
  
  RTK  = RGASC*TKCO
  STK  = 710.0_r8*TKCO
  ACTV = 1+EXP((197500_r8-STK)/RTK)+EXP((STK-218500._r8)/RTK)
  TFNP = EXP(24.269_r8-60000._r8/RTK)/ACTV
  end function calc_leave_grow_tempf

!--------------------------------------------------------------------------------
  pure function calc_plant_maint_tempf(TKCM)result(TFN5)
  implicit none
  real(r8), intent(in) :: TKCM

  real(r8) :: TFN5
  real(r8) :: RTK,STK,ACTVM

  RTK   = RGASC*TKCM
  STK   = 710.0_r8*TKCM
  ACTVM = 1._r8+EXP((195000._r8-STK)/RTK)+EXP((STK-232500._r8)/RTK)
  TFN5  = EXP(25.214_r8-62500._r8/RTK)/ACTVM
  END function calc_plant_maint_tempf
!--------------------------------------------------------------------------------

  subroutine SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig, PltUptake_Ol, &
    PltUptake_Sl, PltUptake_OSl, PltUptake_OSCl,ldebug)
  !
  !DESCRIPTION
  ! Solve uptake from soil-solution supply and root uptake capacity.
  ! The caller already includes C availability and nutrient feedback in UptakePerPlantRateMax.
  ! Both quadratic solutions retain these constraints and soil transport/concentration limits.
  ! SoluteMassMax is a separate population-level available-mass cap.
  ! Output rates are per grid cell [g N or g P cell-1 h-1].
  !
  !Q^2−(v+X-Y+DK)Q+(X−Y)v=0
  !Q is uptake rate per plant
  !v is uptake capacity per plant, with or without the O2Stress multiplier
  !K is the Michaelis-Menten half-saturation parameter
  !X=(q+D)C, with C as micropore solute concentration
  !Y=D*Cm, with Cm being the minimum concentration for uptake

  implicit none
  type(PlantSoluteUptakeConfig_type), intent(in) :: PlantSoluteUptakeConfig
  real(r8), intent(out) :: PltUptake_Ol     !O2- and C-limited uptake before the available-mass cap; soil supply still limits uptake
  real(r8), intent(out) :: PltUptake_Sl     !O2-unlimited uptake with caller C/nutrient constraints and the available-mass cap
  real(r8), intent(out) :: PltUptake_OSl    !approximate C-unlimited diagnostic, obtained by rescaling actual uptake; may exceed the mass cap
  real(r8), intent(out) :: PltUptake_OSCl   !actual uptake with O2, C/nutrient feedback, soil supply, and available-mass constraints
  logical, optional, intent(in) :: ldebug
  real(r8) :: UptakePerPlantRateMax_Ol   !caller C/nutrient-constrained capacity multiplied by O2Stress, [g element plant-1 h-1]
  real(r8) :: X, Y, B, C, BP, CP, delta
  real(r8) :: Uptake, Uptake_Ol
  logical :: lldebug
  ! Inputs apply to one root/mycorrhizal population, soil layer, and band/nonband zone.
  ! Solute mass is expressed as g N or g P. The available-mass cap is used for the hourly uptake step.
  associate(                                                                     &
    SolAdvFlx               => PlantSoluteUptakeConfig%SolAdvFlx               , & !input :advective water flow per plant (q), [m3 plant-1 h-1]
    SolDifusFlx             => PlantSoluteUptakeConfig%SolDifusFlx             , & !input :diffusive transport conductance per plant (D), [m3 plant-1 h-1]
    UptakePerPlantRateMax   => PlantSoluteUptakeConfig%UptakePerPlantRateMax   , & !input :uptake capacity before O2 limitation; caller includes C/nutrient feedback, [g element plant-1 h-1]
    O2Stress                => PlantSoluteUptakeConfig%O2Stress                , & !input :remaining uptake-capacity fraction under O2 limitation, 0=no capacity, 1=unlimited, [-]
    PlantPopulation         => PlantSoluteUptakeConfig%PlantPopulation         , & !input :live plant count converting per-plant uptake to grid-cell uptake, [plants cell-1]
    CAvailStress            => PlantSoluteUptakeConfig%CAvailStress            , & !input :carbon-availability factor (FCUP); positive on entry, 1=no C limitation, [-]
    SoluteMassMax           => PlantSoluteUptakeConfig%SoluteMassMax           , & !input :available solute allocated to this population after competition and minimum-concentration reserve, [g element cell-1]
    SoluteConc              => PlantSoluteUptakeConfig%SoluteConc              , & !input :bulk soil-solution concentration in the selected zone (C), [g element m-3]
    SoluteKM                => PlantSoluteUptakeConfig%SoluteKM                , & !input :Michaelis-Menten uptake half-saturation parameter (K), [g element m-3]
    SoluteConcMin           => PlantSoluteUptakeConfig%SoluteConcMin             & !input :minimum soil-solution concentration for uptake (Cm), [g element m-3]
  )
  lldebug=.false.
  if(present(ldebug))lldebug=ldebug

  UptakePerPlantRateMax_Ol=UptakePerPlantRateMax*O2Stress

  X=(SolDifusFlx+SolAdvFlx)*SoluteConc
  Y=SolDifusFlx*SoluteConcMin

  ! Solve with O2 limitation and caller C/nutrient constraints, before the available-mass cap.
  ! For u^2+B*u+C=0, C=0 gives delta=B*B and the selected root u=0 because B<=0.
  B     = -AZMAX1(UptakePerPlantRateMax_Ol+X-Y+SolDifusFlx*SoluteKM)
  C     = AZMAX1(X-Y)*UptakePerPlantRateMax_Ol
  delta = B*B-4.0_r8*C

  if(delta<0._r8)then
    Uptake_Ol=0._r8
  else
    Uptake_Ol=AZMAX1(-B-SQRT(delta))/2.0_r8
  endif

  ! Solve without O2 limitation; retain caller C/nutrient constraints and soil-solution supply.
  BP    = -AZMAX1(UptakePerPlantRateMax+X-Y+SolDifusFlx*SoluteKM)
  CP    = AZMAX1(X-Y)*UptakePerPlantRateMax
  delta = BP*BP-4.0_r8*CP
  if(delta<0._r8)then
    Uptake=0._r8
  else
    Uptake=AZMAX1(-BP-SQRT(delta))/2.0_r8
  endif

  ! Convert the O2- and C-limited solution to population uptake before the available-mass cap.
  PltUptake_Ol=AZMAX1(Uptake_Ol*PlantPopulation)

  ! Actual uptake: apply the available-mass cap to the O2- and C-limited population rate.
  PltUptake_OSCl=AMIN1(SoluteMassMax,PltUptake_Ol)

  ! O2-unlimited diagnostic: retain C/nutrient constraints and apply the available-mass cap.
  PltUptake_Sl=AMIN1(SoluteMassMax,Uptake*PlantPopulation)

  ! Approximate C-unlimited diagnostic; CAvailStress must be positive.
  ! Division does not exactly undo MIN(FCUP,FZUP/FPUP) in the caller or the nonlinear solver.
  ! This diagnostic is not recapped after rescaling and can exceed SoluteMassMax.
  PltUptake_OSl=PltUptake_OSCl/CAvailStress

  end associate
  end subroutine SoluteUptakeByPlantRoots
!--------------------------------------------------------------------------------

  pure function is_plant_woody_vascular(iPlantRootProfile_pft,iPlant2ndGrothPattern_pft)result(ans)
!
! currently, there are three plant growth types defined as
! iplt_bryophyte=0
! iplt_grasslike=1
!  iplt_treelike=2

  implicit none
  integer, intent(in) :: iPlantRootProfile_pft       !root profile type
  integer, intent(in) :: iPlant2ndGrothPattern_pft   !toggle for expressing secondary growth
  logical :: ans

  ans=iPlantRootProfile_pft > 1 .and. iPlant2ndGrothPattern_pft == 1
  end function is_plant_woody_vascular

!--------------------------------------------------------------------------------

  pure function is_root_bryophyte(iPlantRootProfile_pft)result(ans)
!
! currently, there are three plant growth types defined as
! iplt_bryophyte=0
! iplt_grasslike=1
!  iplt_treelike=2
! only bryophyte is considered as shallow roots
  implicit none
  integer, intent(in) :: iPlantRootProfile_pft
  logical :: ans

  ans=(iPlantRootProfile_pft == iplt_bryophyte)
  end function is_root_bryophyte

!--------------------------------------------------------------------------------
  pure function is_root_N2fix(iPlantNfixType_pft)result(yesno)
  implicit none
  integer, intent(in) :: iPlantNfixType_pft

  logical :: yesno

  yesno=iPlantNfixType_pft.GE.in2fixtyp_root_fast.AND.iPlantNfixType_pft.LE.in2fixtyp_root_slow

  end function is_root_N2fix
!--------------------------------------------------------------------------------
  pure function is_canopy_N2fix(iPlantNfixType_pft)result(yesno)
  implicit none
  integer, intent(in) :: iPlantNfixType_pft

  logical :: yesno

  yesno=iPlantNfixType_pft.GE.in2fixtyp_canopy_fast.AND.iPlantNfixType_pft.LE.in2fixtyp_canopy_slow

  end function is_canopy_N2fix
!--------------------------------------------------------------------------------
  pure function is_plant_N2fix(iPlantNfixType_pft)result(yesno)
  implicit none
  integer, intent(in) :: iPlantNfixType_pft

  logical :: yesno

  yesno=iPlantNfixType_pft.NE.iN2fixtyp_none

  end function is_plant_N2fix
!--------------------------------------------------------------------------------

  pure function fRespWatSens(WFN,iPlantRootProfile)result(ans)

  implicit none
  real(r8), intent(in) :: WFN               !turgor based leaf/root elongation
  integer, intent(in) :: iPlantRootProfile
  real(r8) :: ans

  IF(is_root_bryophyte(iPlantRootProfile))THEN
    ans=WFN**0.10_r8
  ELSE
    ans=WFN**0.25_r8
  ENDIF

  end function fRespWatSens

!----------------------------------------------------------------------------------------------------

  subroutine ExchFluxLimiter(fromState,toState,XFRE)
  implicit none
  real(r8),intent(in) :: fromState
  real(r8),intent(in) :: toState
  real(r8), intent(inout) :: XFRE

  IF(XFRE>0._r8)then
    XFRE=AMIN1(fromState*0.9999_r8,XFRE)
  ELSE  
    XFRE=AMAX1(-toState*0.9999_r8,XFRE)
  ENDIF
  end subroutine ExchFluxLimiter
!----------------------------------------------------------------------------------------------------

  subroutine advect_remap_mass_loss(n, dt, xr, c, Areas, ur, c_new, xL, lost_mass)
    !--------------------------------------------------------------------
    ! Conservative 1D advect-remap with proportional mass loss at right boundary
    ! Workspace arrays are local automatic arrays.
    !
    ! Inputs:
    !   n   - number of cells
    !   xr  - right-edge locations of each cell (size n), strictly increasing
    !   c   - cell-average concentration in each cell (size n)
    !   areas - cross-section area of each cell (size n)
    !   ur  - velocities at right-edge of each cell (size n), nonnegative
    !   dt  - time step (positive)
    !   xL  - left boundary (defaults to 0.0 when omitted)
    !
    ! Outputs:
    !   c_new    - updated cell-average concentration (size n)
    !   lost_mass- total mass lost this step
    !
    ! Local workspace:
    !   xE, xE_star : real(r8), size n+1
    !   dx, m : real(r8), size n
    !   M_star, M_on_fixed : real(r8), size n+1
    !--------------------------------------------------------------------
  implicit none

  integer, intent(in) :: n
  real(r8), intent(in) :: areas(n),xr(n), c(n), ur(n), dt 
  real(r8), optional, intent(in):: xL
  real(r8), intent(out) :: c_new(n)
  real(r8), optional, intent(out) :: lost_mass
  character(len=*), parameter :: subname='advect_remap_mass_loss'
  ! local workspace arrays
  real(r8)  :: xE(n+1), xE_star(n+1)
  real(r8)  :: dx(n), m(n)
  real(r8)  :: M_star(n+1), M_on_fixed(n+1)
  
  ! local scalars
  integer :: i
  real(r8) :: xR_most
  real(r8) :: total_initial
  real(r8) :: dt_res,dt_loc
  logical :: lhalf
  real(r8), parameter :: tiny = 1.0e-14_r8

  call PrintInfo('beg '//subname)
  ! Basic checks (lightweight)
  if (n <= 0) then
      if(present(lost_mass))lost_mass = 0._r8
      return
  end if
  if (dt <= 0._r8) then
      call endrun('Error: dt must be positive.  in '//trim(mod_filename)//' at line',__LINE__)                    
  end if

  do i = 2, n
      if (xr(i) <= xr(i-1)) then
        call endrun('Error: xr must be strictly increasing.  in '//trim(mod_filename)//' at line',__LINE__)                              
      end if
  end do
  do i = 1, n
      if (ur(i) < 0._r8) then
        call endrun('Error: ur must be nonnegative.  in '//trim(mod_filename)//' at line',__LINE__)                    
      end if
  end do

  ! rightmost fixed boundary
  xR_most = xr(n)

  ! fixed edges: xE(1) = xL, xE(2:n+1) = xr(1:n)
  if(present(xL))then
    xE(1) = xL
  else
    xE(1) = 0._r8
  endif
  do i = 1, n
      xE(i+1) = xr(i)
  end do

  ! compute dx and masses
  total_initial = 0._r8
  do i = 1, n
      dx(i) = xE(i+1) - xE(i)
      if (dx(i) <= 0._r8) then
        write(iulog,*)xE(i+1), xE(i), i
        call endrun('Error: non-positive cell width encountered.  in '//trim(mod_filename)//' at line',__LINE__)                    
      end if
  end do
  m = c * dx * Areas * 1.e6_r8
  total_initial = total_initial + sum(m)
  
  dt_res=dt; dt_loc=dt  
  DO
    ! compute moved edges (left edge fixed velocity = 0)

    do
      lhalf=.false.    
      xE_star(1) = xE(1) + 0._r8 * dt_loc      
      do i = 1, n
        xE_star(i+1) = xE(i+1) + ur(i) * dt_loc
        if (xE_star(i+1) + tiny < xE_star(i) .and. xE_star(i+1)<xR_most) then
          dt_loc=dt_loc*0.5_r8
          lhalf=.true.
          exit
        endif
      end do
      if(.not.lhalf)exit 
    enddo

    ! monotonicity check (advected edges should be increasing)
    do i = 2, n+1
      if (xE_star(i) + tiny < xE_star(i-1) .and. xE_star(i)<xR_most) then
        write(iulog,*)'ur',ur
        write(iulog,*)'xR_most',xR_most,dt_loc
        write(iulog,*)'xE_star',xE_star(2:n+1) 
        call endrun('Error: advected edges are not monotone increasing. Reduce dt.  in '//trim(mod_filename)//' at line',__LINE__)                    
      end if
    end do

    ! Remap the full mass on the displaced mesh. Sampling at the fixed
    ! domain edges accounts for boundary outflow exactly once; clipping
    ! masses before interpolation would apply the overlap fraction twice.
    M_star(1) = 0._r8
    do i = 1, n
       M_star(i+1) = M_star(i) + m(i)
    end do

    call interp_linear_clamped(n+1, xE_star, M_star, n+1, xE, M_on_fixed, 0._r8, M_star(n+1))

    ! new cell masses and concentrations'
    
    do i = 1, n
       m(i) = M_on_fixed(i+1) - M_on_fixed(i)    ! reuse m for new masses       
       c_new(i) = m(i)/ (dx(i)*areas(i))
    end do
    dt_res=dt_res-dt_loc
    if(dt_res<dt*1.e-2_r8)exit
    dt_loc=dt_res
  enddo  
  ! Include outflow from every substep, using the existing internal mass scale.
  if(present(lost_mass))lost_mass = total_initial - sum(m)
  c_new=c_new*1.e-6_r8
  call PrintInfo('end '//subname)
  end subroutine advect_remap_mass_loss

  !--------------------------------------------------------------------
  ! Simple piecewise-linear interpolation with clamped extrapolation:
  !   Inputs:
  !     nx  - length of x array
  !     x   - array of x nodes, length nx (must be nondecreasing)
  !     y   - array of y nodes, length nx
  !     nq  - number of query points (size of xq)
  !     xq  - query points (length nq)
  !   Outputs:
  !     yq  - interpolated values at xq (length nq)
  !   Extrapolation:
  !     xq <= x(1) -> y_left
  !     xq >= x(nx) -> y_right
  !--------------------------------------------------------------------
  subroutine interp_linear_clamped(nx, x, y, nq, xq, yq, y_left, y_right)
    implicit none
    integer, intent(in) :: nx, nq
    real(r8), intent(in) :: x(nx), y(nx), xq(nq)
    real(r8), intent(out) :: yq(nq)
    real(r8), intent(in) :: y_left, y_right

    integer :: iq, k
    real(r8) :: t

    do iq = 1, nq
       if (xq(iq) <= x(1)) then
          yq(iq) = y_left
       else if (xq(iq) >= x(nx)) then
          yq(iq) = y_right
       else
          ! find k s.t. x(k) <= xq < x(k+1)
          k = 1
          do while (k < nx - 1 .and. xq(iq) >= x(k+1))
             k = k + 1
          end do

          if (x(k+1) == x(k)) then
             yq(iq) = y(k)
          else
             t = safe_adb(xq(iq) - x(k), x(k+1) - x(k))             
             yq(iq) = (1.0_r8 - t) * AZERO(y(k)) + t * AZERO(y(k+1))
          end if
       end if
    end do
  end subroutine interp_linear_clamped

!--------------------------------------------------------------------
  pure function CalcStomataResist4H2O(NZ)result(Stomata_Resist)
  implicit none
  integer, intent(in) :: NZ
  real(r8) :: Stomata_Stress
  real(r8) :: Stomata_Resist
  
  associate(                                                               &
    RCS_pft                     => plt_photo%RCS_pft                      ,& !input  :e-folding turgor pressure for stomatal resistance, [MPa]
    PSICanopyTurg_pft           => plt_ew%PSICanopyTurg_pft               ,& !input  :plant canopy turgor water potential, [MPa]  
    CanopyMinStomaResistH2O_pft => plt_photo%CanopyMinStomaResistH2O_pft  ,& !input  :canopy minimum stomatal resistance, [s m-1]
    H2OCuticleResist_pft        => plt_photo%H2OCuticleResist_pft          & !input  :maximum stomatal resistance to vapor, [s h-1]
  )
  !greater value of RCS_pft(NZ), more sensitive to turgor change.

  Stomata_Stress = EXP(-PSICanopyTurg_pft(NZ)/RCS_pft(NZ))
  Stomata_Resist = CanopyMinStomaResistH2O_pft(NZ)+(H2OCuticleResist_pft(NZ)-CanopyMinStomaResistH2O_pft(NZ))*Stomata_Stress
  end associate
  end function CalcStomataResist4H2O

!--------------------------------------------------------------------
  pure function is_drought_deciduos(iphenotype)result(ans)
  implicit none
  integer, intent(in) :: iphenotype
  logical :: ans

  ans = iphenotype .GT. iphenotyp_coldecid

  end function is_drought_deciduos
!--------------------------------------------------------------------
  pure function is_cold_deciduos(iphenotype)result(ans)
  implicit none
  integer, intent(in) :: iphenotype
  logical :: ans
  
  ans = iphenotype .EQ.iphenotyp_coldecid .OR. iphenotype.EQ.iphenotyp_coldroutdecid
  end function is_cold_deciduos

!--------------------------------------------------------------------

  SUBROUTINE solve_root_diffusion_step(N, dt, c_init, lumen_areas, layer_thicknesses, &
                                        d_effective, c_next)
  ! ======================================================================
  ! Solves vertical diffusion of cytokinin inside root plumbing for one timestep.
  ! Uses the Implicit Backward Euler method with an integrated Thomas Algorithm.
  ! Zero-flux boundary conditions at both ends of the root profile.
  ! ======================================================================
  INTEGER, INTENT(IN) :: N                                 ! Number of vertical layers
  REAL(r8), INTENT(IN) :: dt                                   ! Timestep size (hour)  
  REAL(r8), DIMENSION(N), INTENT(IN) :: c_init                 ! Initial concentration (mg/m3)
  REAL(r8), DIMENSION(N), INTENT(IN) :: lumen_areas            ! Cumulative lumen area (m2)
  REAL(r8), DIMENSION(N), INTENT(IN) :: layer_thicknesses      ! Thickness of each layer (m)
  REAL(r8), dimension(N),INTENT(IN) :: d_effective                          ! Diffusivity constant (m2/hour)
  REAL(r8), DIMENSION(N), INTENT(OUT) :: c_next                ! Output updated concentration (mg/m3)
  character(len=*), parameter :: subname='solve_root_diffusion_step'
  ! Local Variables for Tridiagonal Matrix: A(i)*C(i-1) + B(i)*C(i) + C(i)*C(i+1) = D(i)
  REAL(r8), DIMENSION(N) :: a_diag, b_diag, c_diag, d_rhs
  REAL(r8), DIMENSION(N) :: lumen_volumes
  REAL(r8) :: area_interface, gamma_interface
  INTEGER :: i

  ! Local variables for the Thomas Algorithm solver
  REAL(r8), DIMENSION(N) :: c_prime, d_prime
  REAL(r8) :: m

  call PrintInfo('beg '//subname)
  if(N==0)return

  ! 1. Calculate layer volumes for mass tracking (Volume = Area * Thickness)
  lumen_volumes = lumen_areas * layer_thicknesses

  ! 2. Initialize storage terms and off-diagonal coefficients.
  a_diag = 0._r8
  b_diag = lumen_volumes
  c_diag = 0._r8
  d_rhs  = c_init * lumen_volumes

  ! 3. Assemble each interface once. The two half-layer diffusion
  ! resistances act in series through a shared interface area.
  ! Use the same conductance in both rows to conserve cytokinin mass.
  DO i = 1, N-1
    area_interface = Harmonicmean_safe(lumen_areas(i), lumen_areas(i+1))
    gamma_interface = dt * area_interface * Harmonicmean_safe( &
      d_effective(i)/layer_thicknesses(i), d_effective(i+1)/layer_thicknesses(i+1))

    b_diag(i)   = b_diag(i)   + gamma_interface
    b_diag(i+1) = b_diag(i+1) + gamma_interface
    c_diag(i)   = c_diag(i)   - gamma_interface
    a_diag(i+1) = a_diag(i+1) - gamma_interface
  END DO
  ! No exterior interfaces: both boundary fluxes are zero.

  ! 4. Execute Thomas Algorithm (Forward Elimination Phase)
  c_prime(1) = c_diag(1) / b_diag(1)
  d_prime(1) = d_rhs(1)  / b_diag(1)

  DO i = 2, N
    m = b_diag(i) - a_diag(i) * c_prime(i-1)
    IF (i < N) c_prime(i) = c_diag(i) / m
    d_prime(i) = (d_rhs(i) - a_diag(i) * d_prime(i-1)) / m
  END DO

  ! 5. Execute Thomas Algorithm (Back Substitution Phase)
  c_next(N) = d_prime(N)
  DO i = N-1, 1, -1
    c_next(i) = d_prime(i) - c_prime(i) * c_next(i+1)
  END DO
  call PrintInfo('end '//subname)
  END SUBROUTINE solve_root_diffusion_step

 !--------------------------------------------------------------------
  pure function smoothstep(v_on, v_full, v)result(ans)
  implicit none

  real(r8), intent(in) :: v_on
  real(r8), intent(in) :: v_full
  real(r8), intent(in) :: v
  real(r8) :: x
  real(r8) :: ans

  if (v_full <= v_on) then
    if (v >= v_full) then
      ans = 1._r8
    else
      ans = 0._r8
    endif
    return
  endif

  x = (v - v_on) / (v_full - v_on)
  x = max(0._r8, min(1._r8, x))
  ans = x * x * (3._r8 - 2._r8 * x)

  end function smoothstep
end module PlantMathFuncMod
