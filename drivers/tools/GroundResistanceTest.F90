program GroundResistanceTest
  ! Check surface gas networks, species scaling, patch areas and substep storage.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use abortutils, only: iulog
  use DebugToolMod, only: PrintInfo
  use SurfPhysMod, only: UpdateSurfaceGasConductancesM
  use HydroThermData, only: RawH2OSoil2CanAir_col,RawH2OSnow2CanAir_col,RawIsoTSoil2CanAir_col
  use SurfPhysData, only: RawH2OLitr2CanAir_col,ResistGrndGas2CanAir_col, &
    FracAsExposedSoil_col,FracEffAsLitR_col
  use LandSurfDataType, only: RawTAtm2CanopySinkZ_col
  use GridDataType, only: AREA_3D,NUM_col
  use EcoSIMSolverPar, only: dts_HeatWatTP
  use SurfSoilDataType, only: CondGasXSnowM_col,FracSurfSnowCovered_col
  use SoilPhysDataType, only: CondGasXSurf_col
  use SurfLitterDataType, only: GasXferCoeffLitr2Atm_col,FracSurfByLitR_col,FracSurfByLitRM_col
  use TracerIDMod, only: idg_beg,idg_NH3,idg_CO2,idg_CH4
  use ChemTracerParsMod, only: WGSG,CGSG,CHSG,OGSG,ZGSG,Z2SG,HGSG,ARSG,ZHSG
  use EcoSimConst, only: BndlResistCO2ToH2ORatio
  use MiniFuncMod, only: TEFGASDIF
  implicit none
  real(r8) :: PoreRatio(idg_beg:idg_NH3),FilmRatio(idg_beg:idg_NH3)
  real(r8) :: SoilBefore(idg_beg:idg_NH3),LitterBefore(idg_beg:idg_NH3)
  real(r8) :: GroundG,LitterG,FullExchange,PartialExchange,FullLitterG,TFAC
  integer :: idg,it

  allocate(RawH2OSoil2CanAir_col(1,1),RawH2OLitr2CanAir_col(1,1),RawH2OSnow2CanAir_col(1,1))
  allocate(RawIsoTSoil2CanAir_col(1,1),ResistGrndGas2CanAir_col(1,1),RawTAtm2CanopySinkZ_col(1,1))
  allocate(FracAsExposedSoil_col(1,1),FracEffAsLitR_col(1,1))
  allocate(FracSurfByLitR_col(1,1),FracSurfByLitRM_col(2,1,1),FracSurfSnowCovered_col(1,1))
  allocate(AREA_3D(3,0:1,1,1),NUM_col(1,1),CondGasXSurf_col(1,1))
  allocate(CondGasXSnowM_col(idg_beg:idg_NH3,2,1,1),GasXferCoeffLitr2Atm_col(idg_beg:idg_NH3,2,1,1))
  PoreRatio=[WGSG/CGSG,WGSG/CHSG,WGSG/OGSG,WGSG/ZGSG,WGSG/Z2SG,WGSG/HGSG,WGSG/ARSG,WGSG/ZHSG]
  FilmRatio=PoreRatio**(2._r8/3._r8)
  FilmRatio(idg_CO2)=BndlResistCO2ToH2ORatio

  ! Bare snow-free aerodynamic transport must be identical for every gas.
  call reset_surface()
  call UpdateSurfaceGasConductancesM(1,1,1,0.5_r8,0._r8)
  call check_close('bare H2O reference',ResistGrndGas2CanAir_col(1,1),0.004_r8)
  do idg=idg_beg,idg_NH3
    call check_close('common bare aerodynamic path',CondGasXSnowM_col(idg,1,1,1),0.25_r8/0.034_r8)
    call check_close('zero litter footprint',GasXferCoeffLitr2Atm_col(idg,1,1,1),0._r8)
  enddo
  SoilBefore=CondGasXSnowM_col(:,1,1,1)
  RawH2OLitr2CanAir_col=50._r8
  call UpdateSurfaceGasConductancesM(1,1,1,500._r8,0._r8)
  do idg=idg_beg,idg_NH3
    call check_close('bare independent of all litter resistance',CondGasXSnowM_col(idg,1,1,1),SoilBefore(idg))
  enddo

  ! With no litter pore resistance, isolate the explicit molecular film.
  call reset_surface()
  FracAsExposedSoil_col=0._r8; FracEffAsLitR_col=1._r8; FracSurfByLitR_col=1._r8
  call UpdateSurfaceGasConductancesM(1,1,1,0._r8,0._r8)
  do idg=idg_beg,idg_NH3
    call check_close('film-only species coefficient',GasXferCoeffLitr2Atm_col(idg,1,1,1), &
      0.25_r8/(0.03_r8+0.004_r8+0.0139_r8*FilmRatio(idg)))
    call check_close('same soil/litter interface without pores',CondGasXSnowM_col(idg,1,1,1), &
      GasXferCoeffLitr2Atm_col(idg,1,1,1))
  enddo
  call check_close('unchanged H2O total',1._r8/CondGasXSurf_col(1,1),0.03_r8+0.0179_r8)
  SoilBefore=CondGasXSnowM_col(:,1,1,1)
  RawTAtm2CanopySinkZ_col=0.13_r8
  call UpdateSurfaceGasConductancesM(1,1,1,0._r8,0._r8)
  do idg=idg_beg,idg_NH3
    call check_close('one species-independent atmospheric increment', &
      0.25_r8/CondGasXSnowM_col(idg,1,1,1)-0.25_r8/SoilBefore(idg),0.1_r8)
  enddo

  ! Litter pore resistance uses inverse diffusivity, separate from film scaling.
  RawTAtm2CanopySinkZ_col=0.03_r8
  call UpdateSurfaceGasConductancesM(1,1,1,0.5_r8,0._r8)
  do idg=idg_beg,idg_NH3
    call check_close('litter pores plus film',CondGasXSnowM_col(idg,1,1,1), &
      0.25_r8/(0.03_r8+0.004_r8+0.0139_r8*FilmRatio(idg)+0.5_r8*PoreRatio(idg)))
  enddo
  SoilBefore=CondGasXSnowM_col(:,1,1,1)
  call UpdateSurfaceGasConductancesM(1,1,1,1.5_r8,0._r8)
  do idg=idg_beg,idg_NH3
    call check_close('inverse-D litter pore increment', &
      0.25_r8/CondGasXSnowM_col(idg,1,1,1)-0.25_r8/SoilBefore(idg),PoreRatio(idg))
  enddo
  call UpdateSurfaceGasConductancesM(1,1,1,1.e12_r8,0._r8)
  call check_close('pore-dominated species limit', &
    CondGasXSnowM_col(idg_CO2,1,1,1)/CondGasXSnowM_col(idg_CH4,1,1,1),CGSG/CHSG)

  ! Four parallel patches must be mixed AFTER each gas's pore/film conversion.
  FracAsExposedSoil_col=0.5_r8; FracEffAsLitR_col=0.5_r8
  FracSurfByLitR_col=0.2_r8; FracSurfSnowCovered_col=0.5_r8
  call UpdateSurfaceGasConductancesM(1,1,1,0.5_r8,0.02_r8)
  do idg=idg_beg,idg_NH3
    GroundG=0.25_r8/0.004_r8 &
      +0.25_r8/(0.004_r8+0.0139_r8*FilmRatio(idg)+0.5_r8*PoreRatio(idg)) &
      +0.25_r8/(0.01_r8+0.04_r8*PoreRatio(idg)) &
      +0.25_r8/(0.01_r8+0.54_r8*PoreRatio(idg))
    LitterG=0.5_r8/(0.004_r8+0.0139_r8*FilmRatio(idg)) &
      +0.5_r8/(0.01_r8+0.04_r8*PoreRatio(idg))
    call check_close('species ground mosaic',CondGasXSnowM_col(idg,1,1,1),0.25_r8/(0.03_r8+1._r8/GroundG))
    call check_close('species litter mosaic and footprint',GasXferCoeffLitr2Atm_col(idg,1,1,1), &
      0.25_r8*0.2_r8/(0.03_r8+1._r8/LitterG))
  enddo
  FracAsExposedSoil_col=1._r8; FracEffAsLitR_col=0._r8
  call UpdateSurfaceGasConductancesM(1,1,1,0.5_r8,1.e12_r8)
  do idg=idg_beg,idg_NH3
    call check_close('impermeable snow leaves aerodynamic bypass',CondGasXSnowM_col(idg,1,1,1), &
      0.25_r8/(0.03_r8+0.008_r8))
  enddo

  ! Full snow removes snow-free films, but retains material diffusion below it.
  FracSurfSnowCovered_col=1._r8
  FracAsExposedSoil_col=0._r8; FracEffAsLitR_col=1._r8
  call UpdateSurfaceGasConductancesM(1,1,1,0.5_r8,0.04_r8)
  do idg=idg_beg,idg_NH3
    call check_close('full snow with subsnow litter pores',CondGasXSnowM_col(idg,1,1,1), &
      0.25_r8/(0.03_r8+0.01_r8+0.54_r8*PoreRatio(idg)))
    call check_close('snow path to litter water',GasXferCoeffLitr2Atm_col(idg,1,1,1), &
      0.25_r8*0.2_r8/(0.03_r8+0.01_r8+0.04_r8*PoreRatio(idg)))
  enddo
  SoilBefore=CondGasXSnowM_col(:,1,1,1); LitterBefore=GasXferCoeffLitr2Atm_col(:,1,1,1)
  RawH2OSoil2CanAir_col=10._r8; RawH2OLitr2CanAir_col=20._r8
  call UpdateSurfaceGasConductancesM(1,1,1,0.5_r8,0.04_r8)
  do idg=idg_beg,idg_NH3
    call check_close('full snow ignores bare film',CondGasXSnowM_col(idg,1,1,1),SoilBefore(idg))
    call check_close('full snow ignores litter film',GasXferCoeffLitr2Atm_col(idg,1,1,1),LitterBefore(idg))
  enddo
  RawH2OSnow2CanAir_col=0.03_r8
  call UpdateSurfaceGasConductancesM(1,1,1,0.5_r8,0.04_r8)
  do idg=idg_beg,idg_NH3
    call check_close('snow aerodynamic increment is common', &
      0.25_r8/CondGasXSnowM_col(idg,1,1,1)-0.25_r8/SoilBefore(idg),0.02_r8)
  enddo
  SoilBefore=CondGasXSnowM_col(:,1,1,1)
  call UpdateSurfaceGasConductancesM(1,1,1,0.5_r8,0.08_r8)
  do idg=idg_beg,idg_NH3
    call check_close('snow pore increment is species dependent', &
      0.25_r8/CondGasXSnowM_col(idg,1,1,1)-0.25_r8/SoilBefore(idg),0.04_r8*PoreRatio(idg))
  enddo

  ! Material temperature scaling is carried by the H2O pore resistance; its
  ! species ratio cancels the common TEFGASDIF factor at that same temperature.
  FracAsExposedSoil_col=1._r8; FracEffAsLitR_col=0._r8
  do it=1,2
    TFAC=TEFGASDIF(250._r8+60._r8*(it-1))
    call UpdateSurfaceGasConductancesM(1,1,1,0._r8,0.04_r8/TFAC)
    do idg=idg_beg,idg_NH3
      call check_close('temperature-scaled snow diffusion',CondGasXSnowM_col(idg,1,1,1), &
        0.25_r8/(0.03_r8+0.03_r8+0.04_r8*PoreRatio(idg)/TFAC))
    enddo
  enddo

  ! Both sides of litter dissolution must scale linearly with litter footprint.
  call reset_surface()
  FracSurfByLitR_col=1._r8
  call UpdateSurfaceGasConductancesM(1,1,1,0._r8,0._r8)
  FullLitterG=GasXferCoeffLitr2Atm_col(idg_CO2,1,1,1)
  FullExchange=FullLitterG*0.25_r8/(FullLitterG+0.25_r8)
  FracSurfByLitR_col=0.2_r8
  call UpdateSurfaceGasConductancesM(1,1,1,0._r8,0._r8)
  LitterG=GasXferCoeffLitr2Atm_col(idg_CO2,1,1,1)
  PartialExchange=LitterG*0.05_r8/(LitterG+0.05_r8)
  call check_close('air/liquid litter-area scaling',PartialExchange,0.2_r8*FullExchange)

  ! Coefficients include area/time and preserve every species at each M.
  SoilBefore=CondGasXSnowM_col(:,1,1,1); LitterBefore=GasXferCoeffLitr2Atm_col(:,1,1,1)
  AREA_3D=6._r8
  call UpdateSurfaceGasConductancesM(2,1,1,0._r8,0._r8)
  do idg=idg_beg,idg_NH3
    call check_close('sixfold soil area/time',CondGasXSnowM_col(idg,2,1,1),6._r8*SoilBefore(idg))
    call check_close('sixfold litter area/time',GasXferCoeffLitr2Atm_col(idg,2,1,1),6._r8*LitterBefore(idg))
  enddo
  FracSurfByLitR_col=0.4_r8; RawH2OLitr2CanAir_col=0.0358_r8
  call UpdateSurfaceGasConductancesM(2,1,1,0._r8,0._r8)
  do idg=idg_beg,idg_NH3
    call check_close('soil species M=1 retained',CondGasXSnowM_col(idg,1,1,1),SoilBefore(idg))
    call check_close('litter species M=1 retained',GasXferCoeffLitr2Atm_col(idg,1,1,1),LitterBefore(idg))
    call check_close('M=2 film with bounded-stability multiplier',GasXferCoeffLitr2Atm_col(idg,2,1,1), &
      1.5_r8*0.4_r8/(0.03_r8+2._r8*(0.004_r8+0.0139_r8*FilmRatio(idg))))
  enddo
  call check_close('M=1 litter footprint retained',FracSurfByLitRM_col(1,1,1),0.2_r8)
  call check_close('M=2 litter footprint retained',FracSurfByLitRM_col(2,1,1),0.4_r8)
  write(iulog,*) 'PASS: species molecular components, common aerodynamics, surface patches and substep storage'

contains
!------------------------------------------------------------------------------------------
  subroutine reset_surface()
    ! Restore an idealized cell with a known neutral aerodynamic/film partition.
    character(len=*), parameter :: subname='reset_surface'
    call PrintInfo('beg '//subname)
    NUM_col=1; AREA_3D=1._r8; dts_HeatWatTP=0.25_r8
    RawH2OSoil2CanAir_col=0.004_r8; RawIsoTSoil2CanAir_col=0.004_r8
    RawH2OLitr2CanAir_col=0.0179_r8; RawH2OSnow2CanAir_col=0.01_r8
    RawTAtm2CanopySinkZ_col=0.03_r8
    FracAsExposedSoil_col=1._r8; FracEffAsLitR_col=0._r8
    FracSurfSnowCovered_col=0._r8; FracSurfByLitR_col=0._r8
    call PrintInfo('end '//subname)
  end subroutine reset_surface
!------------------------------------------------------------------------------------------
  subroutine check_close(label,actual,expected)
    ! Fail when a transport invariant is violated or a coefficient is nonfinite.
    character(len=*), intent(in) :: label !description of the checked invariant, [-]
    real(r8), intent(in) :: actual !computed resistance [h/m] or conductance [m/h or m3]
    real(r8), intent(in) :: expected !reference value in the same units as actual
    character(len=*), parameter :: subname='check_close'
    call PrintInfo('beg '//subname)
    if(.not.ieee_is_finite(actual) .or. abs(actual-expected)>1.e-12_r8*max(1._r8,abs(expected)))then
      write(iulog,*) trim(label),actual,expected
      error stop 'ground resistance invariant failed'
    endif
    call PrintInfo('end '//subname)
  end subroutine check_close
end program GroundResistanceTest
