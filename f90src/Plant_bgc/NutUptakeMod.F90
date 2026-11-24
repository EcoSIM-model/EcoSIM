module NutUptakeMod
  ! Description:
  ! Compute plant nutrient/O2 uptake and root gas/solute exchange using accepted
  ! hydraulic fluxes, independently of the canopy water and energy solver.

  use data_kind_mod, only: r8 => DAT_KIND_R8,yearIJ_type
  use minimathmod,   only: safe_adb, vapsat, AZMAX1,dssign
  use TracerPropMod, only: gas_solubility
  use DebugToolMod, only : PrintInfo
  use PlantDebugMod, only : PrintRootTracer
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantPhotosynthesisAPIData, only : plt_photo
  use PlantRadiationAPIData, only : plt_rad
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantSoilChemistryAPIData, only : plt_soilchem
  use PlantAllometryAPIData, only : plt_allom
  use PlantBiomassAPIData, only : plt_biom
  use PlantEnergyWaterAPIData, only : plt_ew
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use PlantMathFuncMod
  use RootGasMod

  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public :: PlantNutientO2Uptake
  public :: ZeroNutrientUptake
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine PlantNutientO2Uptake(yearIJ,NZ,FDMP, PathLen_pvr,FineRootRadius,FracPRoot4Uptake,&
    FracMinRoot4Uptake_rpvr,FracSoiLayByPrimRoot,RootEffLen4Absorption_pvr)
  !
  !DESCRIPTION
  !doing plant population level nutrient, and O2 uptake
  implicit none
  type(yearIJ_type),intent(in) :: yearIJ
  integer, intent(in) :: NZ
  real(r8), intent(in):: FDMP
  real(r8), intent(in) :: PathLen_pvr(pltpar%jroots,JZ1),FineRootRadius(pltpar%jroots,JZ1),FracPRoot4Uptake(pltpar%jroots,JZ1,JP1)
  real(r8), intent(in) :: FracMinRoot4Uptake_rpvr(pltpar%jroots,JZ1,JP1)
  real(r8), intent(in) :: FracSoiLayByPrimRoot(JZ1,JP1)
  real(r8), intent(in) :: RootEffLen4Absorption_pvr(pltpar%jroots,JZ1)    ! [m]
  
  integer :: I,J
  real(r8)  :: PopPlantO2Uptake,PopPlantO2Demand

  associate(                                           &
    ZERO4Groth_pft    => plt_biom%ZERO4Groth_pft      ,& !input  :threshold zero for plang growth calculation, [-]
    PlantO2Stress_pft => plt_pheno%PlantO2Stress_pft   & !output :plant O2 stress indicator, [-]
  )
  I=yearIJ%I;J=yearIJ%J
  call FoliarNutrientInterception(NZ,FDMP)
!
!     ROOT(N=1) AD MYCORRHIZAL(N=2) O2 AND NUTRIENT UPTAKE
!
  call RootMycoO2NutrientUptake(I,J,NZ,PathLen_pvr,FineRootRadius,&
    FracPRoot4Uptake,FracMinRoot4Uptake_rpvr,FracSoiLayByPrimRoot,RootEffLen4Absorption_pvr &
    ,PopPlantO2Uptake,PopPlantO2Demand)

  call SumNutrientUptake(I,J,NZ)

  IF(PopPlantO2Demand.GT.ZERO4Groth_pft(NZ))THEN
    PlantO2Stress_pft(NZ)=AMAX1(AZMAX1(PopPlantO2Uptake)/PopPlantO2Demand,oscal_test)
  ELSE
    PlantO2Stress_pft(NZ)=1.0_r8
  ENDIF    
  end associate
  end subroutine PlantNutientO2Uptake

!----------------------------------------------------------------------------------------------------
  subroutine FoliarNutrientInterception(NZ,FDMP)

  implicit none
  integer, intent(in) :: NZ
  real(r8), intent(in):: FDMP

  character(len=*), parameter :: subname='FoliarNutrientInterception'
  real(r8) :: FNH3P
  real(r8) :: CNH3P
  real(r8) :: SNH3P
  real(r8) :: ZPOOLB
  integer :: NB

  associate(                                                          &
    TdegCCanopy_pft           => plt_ew%TdegCCanopy_pft              ,& !input  :canopy temperature, [oC]
    NU                        => plt_site%NU                         ,& !input  :current soil surface layer number, [-]
    AREA3                     => plt_site%AREA3                      ,& !input  :soil cross section area (vertical plane defined by its normal direction), [m2]
    ZERO4Groth_pft            => plt_biom%ZERO4Groth_pft             ,& !input  :threshold zero for plang growth calculation, [-]
    AtmGasc                   => plt_site%AtmGasc                    ,& !input  :atmospheric gas concentrations, [g m-3]
    CanopyLeafSheathC_brch    => plt_biom%CanopyLeafSheathC_brch     ,& !input  :plant branch leaf + sheath C, [g d-2]
    CanopyNonstElms_brch      => plt_biom%CanopyNonstElms_brch       ,& !input  :branch nonstructural element, [g d-2]
    LeafPetoNonstElmConc_brch => plt_biom%LeafPetoNonstElmConc_brch  ,& !input  :branch nonstructural C concentration, [g d-2]
    CanPStomaResistH2O_pft    => plt_photo%CanPStomaResistH2O_pft    ,& !input  :canopy stomatal resistance, [h m-1]
    RawH2OCanopy2Atm_pft         => plt_photo%RawH2OCanopy2Atm_pft         ,& !input  :canopy boundary layer resistance, [h m-1]
    LeafAreaLive_brch         => plt_morph%LeafAreaLive_brch         ,& !input  :branch leaf area, [m2 d-2]
    NumOfBranches_pft         => plt_morph%NumOfBranches_pft         ,& !input  :number of branches,[-]
    FracPARads2LiveCanopy_pft => plt_rad%FracPARads2LiveCanopy_pft   ,& !input  :fraction of incoming PAR absorbed by live canopy, [-]
    CanopyLeafArea_pft        => plt_morph%CanopyLeafArea_pft        ,& !input  :plant canopy leaf area, [m2 d-2]
    NH3Dep2Can_brch           => plt_rbgc%NH3Dep2Can_brch             & !output :gaseous NH3 flux fron root disturbance band, [g d-2 h-1]
  )

  call PrintInfo('beg '//subname)
  !
  !     NH3 EXCHANGE BETWEEN CANOPY AND ATMOSPHERE FROM NH3
  !     CONCENTRATION DIFFERENCES 'AtmGasc(idg_NH3)' (ATMOSPHERE FROM 'READS') AND
  !     'CNH3P' (CANOPY), AND FROM STOMATAL + BOUNDARY LAYER RESISTANCE
  !
  !     SNH3P =NH3 solubility at TdegCCanopy_pft
  !     TdegCCanopy_pft=canopy temperature (oC)
  !     FDMP,FNH3P=canopy dry matter content,NH3 concentration
  !     LeafAreaLive_brch,CanopyLeafArea_pft=branch,canopy leaf area
  !     CNH3P,AtmGasc(idg_NH3)=gaseous NH3 concentration in branch,atmosphere
  !     CZPOLB,ZPOOLB=nonstplt_rbgc%RUCtural N concentration,content in branch
  !     NH3Dep2Can_brch=NH3 flux between atmosphere and branch
  !     RA,CanPStomaResistH2O_pft=canopy boundary layer,stomatal resistance
  !     FracPARads2LiveCanopy_pft=fraction of radiation received by each live PFT canopy
  !
  SNH3P = gas_solubility(idg_NH3,TdegCCanopy_pft(NZ))
  FNH3P = 1.0E-04_r8*FDMP
  if(FracPARads2LiveCanopy_pft(NZ).GT.ZERO4Groth_pft(NZ))then
    D105: DO NB=1,NumOfBranches_pft(NZ)
      IF(CanopyLeafSheathC_brch(NB,NZ).GT.ZERO4Groth_pft(NZ) .AND. LeafAreaLive_brch(NB,NZ).GT.ZERO4Groth_pft(NZ) &
        .AND.CanopyLeafArea_pft(NZ).GT.ZERO4Groth_pft(NZ))THEN
        CNH3P                  = AZMAX1(FNH3P*LeafPetoNonstElmConc_brch(ielmn,NB,NZ)/SNH3P)
        ZPOOLB                 = AZMAX1(CanopyNonstElms_brch(ielmn,NB,NZ))
        NH3Dep2Can_brch(NB,NZ) = AMIN1(0.1_r8*ZPOOLB &
          ,AMAX1((AtmGasc(idg_NH3)-CNH3P)/(RawH2OCanopy2Atm_pft(NZ)+CanPStomaResistH2O_pft(NZ)) &
          *FracPARads2LiveCanopy_pft(NZ)*AREA3(NU)*LeafAreaLive_brch(NB,NZ)/CanopyLeafArea_pft(NZ),-0.1_r8*ZPOOLB))
        CanopyNonstElms_brch(ielmn,NB,NZ) = CanopyNonstElms_brch(ielmn,NB,NZ)+NH3Dep2Can_brch(NB,NZ)  
      ELSE
        NH3Dep2Can_brch(NB,NZ)=0.0_r8
      ENDIF
    ENDDO D105
  ELSE
    NH3Dep2Can_brch(1:NumOfBranches_pft(NZ),NZ)=0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine FoliarNutrientInterception

!----------------------------------------------------------------------------------------------------
  subroutine RootMycoO2NutrientUptake(I,J,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,&
    FracMinRoot4Uptake_rpvr,FracSoiLayByPrimRoot,RootEffLen4Absorption_pvr,PopPlantO2Uptake,PopPlantO2Demand)

  implicit none
  integer, intent(in) :: I,J,NZ
  real(r8), intent(in) :: PathLen_pvr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: FineRootRadius(pltpar%jroots,JZ1)
  real(r8), intent(in) :: FracPRoot4Uptake(pltpar%jroots,JZ1,JP1)
  real(r8), intent(in) :: FracMinRoot4Uptake_rpvr(pltpar%jroots,JZ1,JP1)   !minimum active fraction of roots for substance uptake from soil [none]
  real(r8), intent(in) :: FracSoiLayByPrimRoot(JZ1,JP1)
  real(r8), intent(in) :: RootEffLen4Absorption_pvr(pltpar%jroots,JZ1)
  real(r8), intent(out) :: PopPlantO2Uptake
  real(r8), intent(out) :: PopPlantO2Demand

  character(len=*), parameter :: subname='RootMycoO2NutrientUptake'
  real(r8) :: FCUP,FZUP,FPUP,FSatNutTransporter,PerPlantRootH2OUptake
  real(r8) :: dtPerPlantRootH2OUptake,FOXYX,PopPlantO2Uptake_vr
  real(r8) :: RootCO2Ar,RootCO2ArB
  real(r8) :: trc_solml_loc(idg_beg:idg_end)    !local copy of aqueous phase of the volatile tracers
  real(r8) :: trc_gasml_loc(idg_beg:idg_NH3)    !local copy of gaesous phase of the volatile tracers  
  real(r8) :: trc_solml_new(idg_beg:idg_end)    !local copy of aqueous phase of the volatile tracers
  real(r8) :: trc_gasml_new(idg_beg:idg_NH3)    !local copy of gaesous phase of the volatile tracers  
  real(r8) :: trc_solml_copy(idg_beg:idg_end)    !local copy of aqueous phase of the volatile tracers
  real(r8) :: trc_gasml_copy(idg_beg:idg_NH3)    !local copy of gaesous phase of the volatile tracers  
  real(r8) :: dmass0,RootMyMassC
  integer  :: N,L,idg

!     begin_execution
  associate(                                                                &
    THETW_vr                     => plt_soilchem%THETW_vr                  ,& !input  :volumetric water content, [m3 m-3]
    VLSoilPoreMicP_vr            => plt_soilchem%VLSoilPoreMicP_vr         ,& !input  :volume of soil layer, [m3 d-2]
    ZEROS2                       => plt_site%ZEROS2                        ,& !input  :threshold zero for numerical stability,[-]
    NU                           => plt_site%NU                            ,& !input  :current soil surface layer number, [-]
    NK                           => plt_site%NK                            ,& !input  :current hydrologically active layer, [-]
    ZERO                         => plt_site%ZERO                          ,& !input  :threshold zero for numerical stability, [-]
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                ,& !input  :threshold zero for plang growth calculation, [-]
    RAutoRootO2Limter_rpvr       => plt_rbgc%RAutoRootO2Limter_rpvr        ,& !input  :O2 constraint to root respiration (0-1), [-]
    RootRespPotent_pvr           => plt_rbgc%RootRespPotent_pvr            ,& !input  :root respiration unconstrained by O2, [g d-2 h-1]
    RootTotLenPerPlant_pvr       => plt_morph%RootTotLenPerPlant_pvr       ,& !input  :root layer length per plant, [m p-1]
    trcs_solml_vr                => plt_soilchem%trcs_solml_vr             ,& !input  :aqueous tracer, [g d-2]
    trcg_gasml_vr                => plt_soilchem%trcg_gasml_vr             ,& !input  :gas layer mass, [g d-2]
    Myco_pft                     => plt_morph%Myco_pft                     ,& !input  :mycorrhizal type (no or yes),[-]
    RootAbsorbLenPerPlant_pvr    => plt_morph%RootAbsorbLenPerPlant_pvr    ,& !input  :root layer length density, [m m-3]
    PopuRootMycoC_pvr            => plt_biom% PopuRootMycoC_pvr            ,& !input  :root layer C, [gC d-2]
    RootVH2O_pvr                 => plt_morph%RootVH2O_pvr                 ,& !input  :root layer volume water, [m2 d-2]
    Root1stActStructElms_rpvr    => plt_biom%Root1stActStructElms_rpvr     ,& !inoput :root layer active zone element in primary axes, [g d-2]
    RootMyco2ndStrutElms_rpvr    => plt_biom%RootMyco2ndStrutElms_rpvr     ,& !input  :root layer element secondary axes, [g d-2]
    RootMycoNonstElms_rpvr       => plt_biom%RootMycoNonstElms_rpvr        ,& !input  :root layer nonstructural element, [g d-2]
    NMaxRootBotLayer_pft         => plt_morph%NMaxRootBotLayer_pft         ,& !input  :maximum soil layer number for all root axes,[-]
    NumStructuralRootAxes_pft    => plt_morph%NumStructuralRootAxes_pft    ,& !input  :number of structural root axes,[-]
    trcs_deadroot2soil_pvr       => plt_rbgc%trcs_deadroot2soil_pvr        ,& !inoput :gases released to soil upong dying roots, [g d-2 h-1]
    trcg_rootml_pvr              => plt_rbgc%trcg_rootml_pvr               ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr              => plt_rbgc%trcs_rootml_pvr               ,& !inoput :root aqueous content, [g d-2]
    RootCO2Ar2Soil_pvr           => plt_rbgc%RootCO2Ar2Soil_pvr            ,& !inoput :root respiration released to soil, [gC d-2 h-1]
    RootO2Dmnd4Resp_pvr          => plt_rbgc%RootO2Dmnd4Resp_pvr            & !output :root O2 demand from respiration, [g d-2 h-1]
  )

  call PrintInfo('beg '//subname)

  PopPlantO2Uptake = 0._r8
  PopPlantO2Demand = 0._r8
  RootCO2Ar        = 0._r8
  RootCO2ArB       = 0._r8

!  dmass0=sum(trcs_solml_vr(idg_O2,NU:plt_site%NL))+sum(trcg_gasml_vr(idg_O2,NU:plt_site%NL))
!  write(904,*)I*1000+J/24.,dmass0,'beg'
  
  D950: DO L=NU,NK
    IF(VLSoilPoreMicP_vr(L).GT.ZEROS2 .AND. THETW_vr(L).GT.ZERO .AND. L<=NMaxRootBotLayer_pft(NZ)) then
      trc_solml_new  = 0._r8;trc_gasml_new = 0._r8
      
      D955: DO N  = 1, Myco_pft(NZ)      
        if(RootAbsorbLenPerPlant_pvr(N,L,NZ).GT.ZERO .AND. RootVH2O_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ) &
          .AND. PopuRootMycoC_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN

          call GetUptakeCapcity(I,J,N,L,NZ,FracPRoot4Uptake,FracMinRoot4Uptake_rpvr,FCUP,FZUP,FPUP,&
            FSatNutTransporter,PerPlantRootH2OUptake,dtPerPlantRootH2OUptake,FOXYX)

          !
          !     ROOT O2 DEMAND CALCULATED FROM O2 NON-LIMITED RESPIRATION RATE
          !
          !     RootO2Dmnd4Resp_pvr=O2 demand, g O2
          !     RootRespPotent_pvr=respiration unlimited by O2
          !     RootVH2O_pvr=root or myco aqueous volume
          !     FOXYX=fraction of total O2 demand from previous hour
          !
          RootO2Dmnd4Resp_pvr(N,L,NZ)=2.667_r8*RootRespPotent_pvr(N,L,NZ)

          !partition the gas for local uptake 
          DO idg=idg_beg,idg_end
            if(idg.NE.idg_O2)then
              trc_solml_loc(idg)=trcs_solml_vr(idg,L)*FracPRoot4Uptake(N,L,NZ)
            endif
          enddo
          trc_solml_loc(idg_CO2) = AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_CO2))
          trc_solml_loc(idg_O2)  = trcs_solml_vr(idg_O2,L)*FOXYX

          !  the two lines below may be redundant

          DO idg=idg_beg,idg_NH3
            if(idg/=idg_O2)then
              trc_gasml_loc(idg)=AZMAX1(trcg_gasml_vr(idg,L)*FracPRoot4Uptake(N,L,NZ))     
            endif          
          enddo
          trc_gasml_loc(idg_O2)  = AZMAX1(trcg_gasml_vr(idg_O2,L)*FOXYX)
          trc_gasml_loc(idg_CO2) = AZMAX1(trc_gasml_loc(idg_CO2))

          call RootSoilGasExchange(I,J,N,L,NZ,FineRootRadius,FracPRoot4Uptake,FracSoiLayByPrimRoot,&
            RootEffLen4Absorption_pvr,dtPerPlantRootH2OUptake,FOXYX,trc_gasml_loc,trc_solml_loc,PopPlantO2Uptake_vr)
        
          PopPlantO2Demand = PopPlantO2Demand+RootO2Dmnd4Resp_pvr(N,L,NZ)
          PopPlantO2Uptake = PopPlantO2Uptake+PopPlantO2Uptake_vr

          !update the soil gas
          DO idg=idg_beg,idg_NH3
            trc_gasml_new(idg)=trc_gasml_new(idg)+trc_gasml_loc(idg)
            trc_solml_new(idg)=trc_solml_new(idg)+trc_solml_loc(idg)
          ENDDO    

          idg=idg_NH3B
          trc_solml_new(idg)=trc_solml_new(idg)+trc_solml_loc(idg)

          call RootExudates(I,J,N,L,NZ)
          ! Carbon basis for NH4/NO3 capacity diagnostics: population secondary
          ! fine-root/mycorrhizal structural C plus nonstructural C in this layer.
          RootMyMassC=sum(RootMyco2ndStrutElms_rpvr(ielmc,N,L,1:NumStructuralRootAxes_pft(NZ),NZ)) &
            +RootMycoNonstElms_rpvr(ielmc,N,L,NZ)
          if(N==ipltroot)then
            ! Plant roots also include active primary-root structural C.
            RootMyMassC=RootMyMassC+sum(Root1stActStructElms_rpvr(ielmc,L,1:NumStructuralRootAxes_pft(NZ),NZ))
          endif
          !
          !     NUTRIENT UPTAKE
          !
          !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
          !     FCUP=limitation to active uptake respiration from CPOOLR
          !     FSatNutTransporter=protein concentration relative to 5%
          !     RootTotLenPerPlant_pvr=root,myco length per plant
          !
          IF(RAutoRootO2Limter_rpvr(N,L,NZ).GT.ZERO .AND. FCUP.GT.ZERO .AND. FSatNutTransporter.GT.ZERO &
            .AND. RootTotLenPerPlant_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
            !
            !     FZUP=limitn to active uptake respiration from CZPOLR
            !         
            call UptakeMineralNitrogen(I,J,N,L,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,FracMinRoot4Uptake_rpvr,&
              RootEffLen4Absorption_pvr,FCUP,FZUP,FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)
            !
            !     FPUP=limitn to active uptake respiration from CPPOLR
            !         
            call UptakeMineralPhosporhus(N,L,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,FracMinRoot4Uptake_rpvr,&
              RootEffLen4Absorption_pvr,FCUP,FPUP,FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)
          ENDIF
          RootCO2Ar=RootCO2Ar-plt_rbgc%RootCO2AutorX_pvr(N,L,NZ)      
        else
          ! Inactive absorptive roots must not compete for next-hour soil O2.
          RootO2Dmnd4Resp_pvr(N,L,NZ) = 0._r8
          RAutoRootO2Limter_rpvr(N,L,NZ) = 1._r8
          ! Absorptive geometry can disappear after respiration was calculated.
          ! Release previous-hour CO2 even though uptake is now inactive.
          RootCO2Ar2Soil_pvr(L,NZ) = RootCO2Ar2Soil_pvr(L,NZ)-plt_rbgc%RootCO2AutorX_pvr(N,L,NZ)
        ENDIF
      ENDDO D955      
    ELSE
      D956: DO N  = 1, Myco_pft(NZ)          
        ! Dry or out-of-root-zone layers contribute no current uptake demand.
        RootO2Dmnd4Resp_pvr(N,L,NZ) = 0._r8
        IF(L.EQ.NMaxRootBotLayer_pft(NZ)+1)THEN  
          RAutoRootO2Limter_rpvr(N,L,NZ) = RAutoRootO2Limter_rpvr(N,L-1,NZ)
        ELSE
          RAutoRootO2Limter_rpvr(N,L,NZ) = 0._r8
        ENDIF                      
        RootCO2ArB               = RootCO2ArB-plt_rbgc%RootCO2AutorX_pvr(N,L,NZ)
        RootCO2Ar2Soil_pvr(L,NZ) = RootCO2Ar2Soil_pvr(L,NZ)-plt_rbgc%RootCO2AutorX_pvr(N,L,NZ)
        DO idg=idg_beg,idg_NH3
          trcs_deadroot2soil_pvr(idg,L,NZ) = trcs_deadroot2soil_pvr(idg,L,NZ) + trcg_rootml_pvr(idg,N,L,NZ)
          trcs_deadroot2soil_pvr(idg,L,NZ) = trcs_deadroot2soil_pvr(idg,L,NZ) + trcs_rootml_pvr(idg,N,L,NZ)
          trcg_rootml_pvr(idg,N,L,NZ)      = 0._r8
          trcs_rootml_pvr(idg,N,L,NZ)      = 0._r8
        ENDDO
      ENDDO D956      
    ENDIF
  ENDDO D950
  
  call PrintInfo('end '//subname)
  end associate
  end subroutine RootMycoO2NutrientUptake

!----------------------------------------------------------------------------------------------------
  subroutine ZeroNutrientUptake()
    ! Description:
    ! Reset chemical uptake and gas-exchange outputs without changing accepted
    ! hydraulic or canopy energy fluxes. Retain previous O2 demand until each
    ! PFT reads its competition share; native handoff owns persistent CO2 debt.
    character(len=*), parameter :: subname='ZeroNutrientUptake'
    call PrintInfo('beg '//subname)
    associate( &
      RootSoilCO2Debt_vr => plt_rbgc%RootSoilCO2Debt_vr, & ! output: transient root/soil CO2 debt [gC/cell], layer
      RootCO2Ar2RootX_rpvr => plt_rbgc%RootCO2Ar2RootX_rpvr, & ! output: root/myco respiration released to root/myco, [gC /cell/hour], root type/layer/PFT
      RootCO2Ar2Soil_pvr => plt_rbgc%RootCO2Ar2Soil_pvr, & ! output: root respiration released to soil, [gC /cell/hour], layer/PFT
      RootO2_TotSink_pvr => plt_bgcr%RootO2_TotSink_pvr, & ! output: autotrophic respiration O2 sink [gO2/cell/hour], root type/layer/PFT
      RootCO2Ar2RootX_pvr => plt_rbgc%RootCO2Ar2RootX_pvr, & ! output: root respiration released to root, [gC /cell/hour], layer/PFT
      RCO2Emis2Root_rpvr => plt_rbgc%RCO2Emis2Root_rpvr, & ! output: aqueous CO2 flux from roots to root water, [g /cell/hour], root type/layer/PFT
      Soil2RootMycoExudE_pvr => plt_rbgc%Soil2RootMycoExudE_pvr, & ! output: root uptake (+ve) - exudation (-ve) of DOE, [g /cell/hour], element/root type/complex/layer/PFT
      RAutoRootO2Limter_rpvr => plt_rbgc%RAutoRootO2Limter_rpvr, & ! output: oxygen limitation of respiration [-], root type/layer/PFT
      RootNH4DmndSoil_pvr => plt_rbgc%RootNH4DmndSoil_pvr, & ! output: root uptake of NH4 non-band unconstrained by NH4, [g /cell/hour], root type/layer/PFT
      RootNutUptake_pvr => plt_rbgc%RootNutUptake_pvr, & ! output: root uptake of Nutrient band, [g /cell/hour], species/root type/layer/PFT
      VmaxNH4Root_pvr => plt_rbgc%VmaxNH4Root_pvr, & ! output: NH4 uptake capacity [gN/gC/hour], root type/layer/PFT
      VmaxNO3Root_pvr => plt_rbgc%VmaxNO3Root_pvr, & ! output: NO3 uptake capacity [gN/gC/hour], root type/layer/PFT
      Nutruptk_fClim_rpvr => plt_bgcr%Nutruptk_fClim_rpvr, & ! output: climate uptake modifier [-], root type/layer/PFT
      Nutruptk_fNlim_rpvr => plt_bgcr%Nutruptk_fNlim_rpvr, & ! output: nitrogen uptake modifier [-], root type/layer/PFT
      Nutruptk_fPlim_rpvr => plt_bgcr%Nutruptk_fPlim_rpvr, & ! output: phosphorus uptake modifier [-], root type/layer/PFT
      Nutruptk_fProtC_rpvr => plt_bgcr%Nutruptk_fProtC_rpvr, & ! output: protein uptake modifier [-], root type/layer/PFT
      RootOUlmNutUptake_pvr => plt_rbgc%RootOUlmNutUptake_pvr, & ! output: root uptake of NH4 band unconstrained by O2, [g /cell/hour], species/root type/layer/PFT
      RootNH4DmndBand_pvr => plt_rbgc%RootNH4DmndBand_pvr, & ! output: band NH4 demand [gN/cell/hour], root type/layer/PFT
      RootNO3DmndSoil_pvr => plt_rbgc%RootNO3DmndSoil_pvr, & ! output: non-band NO3 demand [gN/cell/hour], root type/layer/PFT
      RootNO3DmndBand_pvr => plt_rbgc%RootNO3DmndBand_pvr, & ! output: band NO3 demand [gN/cell/hour], root type/layer/PFT
      RootCUlmNutUptake_pvr => plt_rbgc%RootCUlmNutUptake_pvr, & ! output: root uptake of NH4 band unconstrained by root nonstructural C, [g /cell/hour], species/root type/layer/PFT
      RootH2PO4DmndSoil_pvr => plt_rbgc%RootH2PO4DmndSoil_pvr, & ! output: root uptake of H2PO4 non-band, [g /cell/hour], root type/layer/PFT
      RootH2PO4DmndBand_pvr => plt_rbgc%RootH2PO4DmndBand_pvr, & ! output: root uptake of H2PO4 band, [g /cell/hour], root type/layer/PFT
      RootH1PO4DmndSoil_pvr => plt_rbgc%RootH1PO4DmndSoil_pvr, & ! output: HPO4 demand in non-band by each root population, [g /cell/hour], root type/layer/PFT
      RootH1PO4DmndBand_pvr => plt_rbgc%RootH1PO4DmndBand_pvr, & ! output: HPO4 demand in band by each root population, [g /cell/hour], root type/layer/PFT
      RootN2Fix_pvr => plt_bgcr%RootN2Fix_pvr, & ! output: root N2 fixation, [gN /cell/hour], layer/PFT
      RootCO2Emis2Root_pvr => plt_bgcr%RootCO2Emis2Root_pvr, & ! output: root CO2 flux, [gC /cell/hour], layer/PFT
      RootCO2Emis2Root_vr => plt_bgcr%RootCO2Emis2Root_vr, & ! output: total root CO2 flux, [gC /cell/hour], layer
      trcs_Soil2plant_uptake_pvr => plt_rbgc%trcs_Soil2plant_uptake_pvr, & ! output: plant root-soil solute flux non-band, [g /cell/hour], species/layer/PFT
      Soil2RootMycoExudE_pft => plt_rbgc%Soil2RootMycoExudE_pft, & ! output: root uptake minus exudation [g element/cell/hour], element/PFT
      RootNH4Uptake_pft => plt_rbgc%RootNH4Uptake_pft, & ! output: total root uptake of NH4, [g /cell/hour], PFT
      RootNO3Uptake_pft => plt_rbgc%RootNO3Uptake_pft, & ! output: total root uptake of NO3, [g /cell/hour], PFT
      RootH2PO4Uptake_pft => plt_rbgc%RootH2PO4Uptake_pft, & ! output: total root uptake of PO4, [g /cell/hour], PFT
      RootHPO4Uptake_pft => plt_rbgc%RootHPO4Uptake_pft, & ! output: total root uptake of HPO4, [g /cell/hour], PFT
      RootN2Fix_pft => plt_rbgc%RootN2Fix_pft, & ! output: total root N2 fixation, [g /cell/hour], PFT
      RootO2Uptk_pvr => plt_rbgc%RootO2Uptk_pvr, & ! output: aqueous O2 flux from roots to root water, [g /cell/hour], root type/layer/PFT
      RootUptkSoiSol_pvr => plt_rbgc%RootUptkSoiSol_pvr, & ! output: aqueous CO2 flux from roots to soil water, [g /cell/hour], species/root type/layer/PFT
      trcg_air2root_flx_pvr => plt_rbgc%trcg_air2root_flx_pvr, & ! output: gaseous tracer flux through roots, [g /cell/hour], species/root type/layer/PFT
      trcg_Root_gas2aqu_flx_vr => plt_rbgc%trcg_Root_gas2aqu_flx_vr, & ! output: dissolution (+ve) - volatilization (-ve) gas flux in roots, [g /cell/hour], species/root type/layer/PFT
      trcs_Soil2plant_uptake_vr => plt_rbgc%trcs_Soil2plant_uptake_vr) ! output: total root-soil solute flux non-band, [g /cell/hour], species/layer
      RootSoilCO2Debt_vr=0._r8
      RootCO2Ar2RootX_rpvr=0._r8
      RootCO2Ar2Soil_pvr=0._r8
      RootO2_TotSink_pvr=0._r8
      RootCO2Ar2RootX_pvr=0._r8
      RCO2Emis2Root_rpvr=0._r8
      Soil2RootMycoExudE_pvr=0._r8
      RAutoRootO2Limter_rpvr=0._r8
      RootNH4DmndSoil_pvr=0._r8
      RootNutUptake_pvr=0._r8
      VmaxNH4Root_pvr=0._r8
      VmaxNO3Root_pvr=0._r8
      Nutruptk_fClim_rpvr=0._r8
      Nutruptk_fNlim_rpvr=0._r8
      Nutruptk_fPlim_rpvr=0._r8
      Nutruptk_fProtC_rpvr=0._r8
      RootOUlmNutUptake_pvr=0._r8
      RootNH4DmndBand_pvr=0._r8
      RootNO3DmndSoil_pvr=0._r8
      RootNO3DmndBand_pvr=0._r8
      RootCUlmNutUptake_pvr=0._r8
      RootH2PO4DmndSoil_pvr=0._r8
      RootH2PO4DmndBand_pvr=0._r8
      RootH1PO4DmndSoil_pvr=0._r8
      RootH1PO4DmndBand_pvr=0._r8
      RootN2Fix_pvr=0._r8
      RootCO2Emis2Root_pvr=0._r8
      RootCO2Emis2Root_vr=0._r8
      trcs_Soil2plant_uptake_pvr=0._r8
      Soil2RootMycoExudE_pft=0._r8
      RootNH4Uptake_pft=0._r8
      RootNO3Uptake_pft=0._r8
      RootH2PO4Uptake_pft=0._r8
      RootHPO4Uptake_pft=0._r8
      RootN2Fix_pft=0._r8
      RootO2Uptk_pvr=0._r8
      RootUptkSoiSol_pvr=0._r8
      trcg_air2root_flx_pvr=0._r8
      trcg_Root_gas2aqu_flx_vr=0._r8
      trcs_Soil2plant_uptake_vr=0._r8
    end associate
    call PrintInfo('end '//subname)
  end subroutine ZeroNutrientUptake


!----------------------------------------------------------------------------------------------------
  subroutine UptakeMineralPhosporhus(N,L,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,FracMinRoot4Uptake_rpvr,&
    RootEffLen4Absorption_pvr,FCUP,FPUP,FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)

  implicit none
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: PathLen_pvr(pltpar%jroots,JZ1),FineRootRadius(pltpar%jroots,JZ1),FracPRoot4Uptake(pltpar%jroots,JZ1,JP1)
  real(r8), intent(in) :: FracMinRoot4Uptake_rpvr(pltpar%jroots,JZ1,JP1)  !the minimum root fraction involved in substrate uptake
  real(r8), intent(in):: RootEffLen4Absorption_pvr(pltpar%jroots,JZ1)
  real(r8), intent(in):: FCUP
  real(r8), intent(in):: FPUP
  real(r8), intent(in):: FSatNutTransporter
  real(r8), intent(in):: PerPlantRootH2OUptake
  real(r8), intent(in):: RootMyMassC
  character(len=*), parameter :: subname='UptakeMineralPhosporhus'

  real(r8) :: TFPO4X,TFP14X,TFPOBX,TFP1BX
  real(r8) :: DIFFL
  real(r8) :: FP14X,FP1BX
  real(r8) :: FPO4X,FPOBX
  real(r8) :: POSGX
  real(r8) :: PATHL
  associate(                                                            &
    TortMicPM_vr                 => plt_site%TortMicPM_vr              ,& !input  :micropore soil tortuosity, [m3 m-3]
    ZEROS                        => plt_site%ZEROS                     ,& !input  :threshold zero for numerical stability,[-]
    ZERO2                        => plt_site%ZERO2                     ,& !input  :threshold zero for numerical stability,[-]
    RootH1PO4DmndBandPrev_pvr    => plt_rbgc%RootH1PO4DmndBandPrev_pvr ,& !input  :HPO4 demand in band by each root population, [g d-2 h-1]
    RootH1PO4DmndSoilPrev_pvr    => plt_rbgc%RootH1PO4DmndSoilPrev_pvr ,& !input  :HPO4 demand in non-band by each root population, [g d-2 h-1]
    RootH2PO4DmndSoilPrev_pvr    => plt_rbgc%RootH2PO4DmndSoilPrev_pvr ,& !input  :root uptake of H2PO4 non-band, [g d-2 h-1]
    RootH2PO4DmndBandPrev_pvr    => plt_rbgc%RootH2PO4DmndBandPrev_pvr ,& !input  :root uptake of H2PO4 band, [g d-2 h-1]
    SoluteDifusvtyT_vr            => plt_soilchem%SoluteDifusvtyT_vr     ,& !input  :aqueous diffusivity, [m2 h-1]
    RH1PO4EcoDmndBandPrev_vr     => plt_bgcr%RH1PO4EcoDmndBandPrev_vr  ,& !input  :HPO4 demand in band by all microbial, root, myco populations, [gP d-2 h-1]
    RH2PO4EcoDmndBandPrev_vr     => plt_bgcr%RH2PO4EcoDmndBandPrev_vr  ,& !input  :total root + microbial PO4 uptake band, [gP d-2 h-1]
    RH1PO4EcoDmndSoilPrev_vr     => plt_bgcr%RH1PO4EcoDmndSoilPrev_vr  ,& !input  :HPO4 demand in non-band by all microbial, root, myco populations, [gP d-2 h-1]
    RH2PO4EcoDmndSoilPrev_vr     => plt_bgcr%RH2PO4EcoDmndSoilPrev_vr   & !input  :total root + microbial PO4 uptake non-band, [gP d-2 h-1]
  )
  call PrintInfo('beg '//subname)
  TFPO4X=0.0_r8
  TFPOBX=0.0_r8
  TFP14X=0.0_r8
  TFP1BX=0.0_r8

!     begin_execution
  IF(RH2PO4EcoDmndSoilPrev_vr(L).GT.ZEROS)THEN
    FPO4X=AMAX1(FracMinRoot4Uptake_rpvr(N,L,NZ),RootH2PO4DmndSoilPrev_pvr(N,L,NZ)/RH2PO4EcoDmndSoilPrev_vr(L))
  ELSE
    FPO4X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RH2PO4EcoDmndBandPrev_vr(L).GT.ZEROS)THEN
    FPOBX=AMAX1(FracMinRoot4Uptake_rpvr(N,L,NZ),RootH2PO4DmndBandPrev_pvr(N,L,NZ)/RH2PO4EcoDmndBandPrev_vr(L))
  ELSE
    FPOBX=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RH1PO4EcoDmndSoilPrev_vr(L).GT.ZEROS)THEN
    FP14X=AMAX1(FracMinRoot4Uptake_rpvr(N,L,NZ),RootH1PO4DmndSoilPrev_pvr(N,L,NZ)/RH1PO4EcoDmndSoilPrev_vr(L))
  ELSE
    FP14X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RH1PO4EcoDmndBandPrev_vr(L).GT.ZEROS)THEN
    FP1BX=AMAX1(FracMinRoot4Uptake_rpvr(N,L,NZ),RootH1PO4DmndBandPrev_pvr(N,L,NZ)/RH1PO4EcoDmndBandPrev_vr(L))
  ELSE
    FP1BX=FracPRoot4Uptake(N,L,NZ)
  ENDIF

  TFPO4X=TFPO4X+FPO4X
  TFP14X=TFP14X+FP14X

  TFPOBX=TFPOBX+FPOBX
  TFP1BX=TFP1BX+FP1BX

  IF(FPUP.GT.ZERO2)THEN
    !
    !     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF H2PO4,HPO4
    !     FROM SOIL TO ROOT
    !
    !     POSGL=PO4 diffusivity
    !     TortMicPM_vr=soil tortuosity
    !     PathLen_pvr=path length of water and nutrient uptake
    !     FineRootRadius=root radius
    !     DIFFL=PO4 diffusion per plant
    !
    POSGX = SoluteDifusvtyT_vr(ids_H1PO4,L)*TortMicPM_vr(NPH,L)
    PATHL = AMIN1(PathLen_pvr(N,L),FineRootRadius(N,L)+SQRT(2.0_r8*POSGX))
    DIFFL = POSGX*safe_adb(RootEffLen4Absorption_pvr(N,L),LOG(PATHL/FineRootRadius(N,L)))

    call UptakeH2PO4(N,L,NZ,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)

    call UptakeHPO4(N,L,NZ,DIFFL,FP14X,FP1BX,FCUP,FPUP,FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)

  ENDIF
  call PrintInfo('end '//subname)
  end associate
  end subroutine UptakeMineralPhosporhus

!----------------------------------------------------------------------------------------------------
  subroutine UptakeNO3(I,J,N,L,NZ,FNO3X,FNOBX,PathLen_pvr,FineRootRadius,RootEffLen4Absorption_pvr,FCUP,FZUP,&
    FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  real(r8), intent(in):: FNO3X,FNOBX
  real(r8), intent(in) :: PathLen_pvr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: FineRootRadius(pltpar%jroots,JZ1)
  real(r8), intent(in) :: RootEffLen4Absorption_pvr(pltpar%jroots,JZ1)
  real(r8), intent(in):: FCUP
  real(r8), intent(in):: FZUP
  real(r8), intent(in):: FSatNutTransporter
  real(r8), intent(in):: PerPlantRootH2OUptake
  real(r8), intent(in):: RootMyMassC !population-total root/mycorrhizal carbon in this layer, [gC d-2]
  character(len=*), parameter :: subname='UptakeNO3'

  real(r8) :: DIFFL
  real(r8) :: DIFNO3,DIFNOB
  real(r8) :: PATHL
  real(r8) :: RMFNO3,RTKNO3,RTKNOP,RMFNOB,RTKNOB,RTKNPB
  real(r8) :: UPMX,UPMXP
  real(r8) :: VmaxNO3RootPerPlant !uptake capacity before O2 limitation, [gN plant-1 h-1]
  real(r8) :: ZOSGX,ZNO3M,ZNO3X,ZNOBM,ZNOBX
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig
  logical :: ldebug

! begin_execution
  associate(                                                    &
    PlantPopuLive_pft      => plt_site%PlantPopuLive_pft       ,& !input  :plant population, [d-2]
    TortMicPM_vr           => plt_site%TortMicPM_vr            ,& !input  :micropore soil tortuosity, [m3 m-3]
    ZERO                   => plt_site%ZERO                    ,& !input  :threshold zero for numerical stability, [-]
    RootSAreaPerPlant_pvr  => plt_morph%RootSAreaPerPlant_pvr  ,& !input  :root layer area per plant, [m2 p-1]
    fTgrowRootP_vr         => plt_pheno%fTgrowRootP_vr         ,& !input  :root layer temperature growth functiom, [-]
    RAutoRootO2Limter_rpvr => plt_rbgc%RAutoRootO2Limter_rpvr  ,& !input  :O2 constraint to root respiration (0-1), [-]
    VmaxNO3Root_pft        => plt_rbgc%VmaxNO3Root_pft         ,& !input  :maximum root NO3 uptake rate, [g m-2 h-1]
    CminNO3Root_pft        => plt_rbgc%CminNO3Root_pft         ,& !input  :minimum NO3 concentration for root NH4 uptake, [g m-3]
    KmNO3Root_pft          => plt_rbgc%KmNO3Root_pft           ,& !input  :Km for root NO3 uptake, [g m-3]
    RootCUlmNutUptake_pvr  => plt_rbgc%RootCUlmNutUptake_pvr   ,& !input  :root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
    RootNO3DmndSoil_pvr    => plt_rbgc%RootNO3DmndSoil_pvr     ,& !input  :root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr       ,& !input  :root uptake of Nutrient band, [g d-2 h-1]
    RootNO3DmndBand_pvr    => plt_rbgc%RootNO3DmndBand_pvr     ,& !input  :root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
    RootOUlmNutUptake_pvr  => plt_rbgc%RootOUlmNutUptake_pvr   ,& !input  :root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
    SoluteDifusvtyT_vr     => plt_soilchem%SoluteDifusvtyT_vr  ,& !input  :aqueous diffusivity, [m2 h-1]
    trcs_solml_vr          => plt_soilchem%trcs_solml_vr       ,& !input  :aqueous tracer, [g d-2]
    trcs_VLN_vr            => plt_soilchem%trcs_VLN_vr         ,& !input  :effective relative tracer volume, [-]
    trc_solcl_vr           => plt_soilchem%trc_solcl_vr        ,& !input  :aqueous tracer concentration, [g m-3]
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr        ,& !input  :soil micropore water content, [m3 d-2]
    VmaxNO3Root_pvr        => plt_rbgc%VmaxNO3Root_pvr          & !output :maximum NO3 uptake rate,  [gN h-1 (gC root)-1]
  )
  call PrintInfo('beg '//subname)
  ldebug=.false.
!
! PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NO3
! FROM SOIL TO ROOT
!
! ZOSGL=NO3 diffusivity
! TortMicPM_vr=soil tortuosity
! FineRootRadius=root radius
! PATH=path length of water and nutrient uptake
! DIFFL=NO3 diffusion per plant
!
  ZOSGX = SoluteDifusvtyT_vr(ids_NO3,L)*TortMicPM_vr(NPH,L)
  PATHL = AMIN1(PathLen_pvr(N,L),FineRootRadius(N,L)+SQRT(2.0_r8*ZOSGX))
  DIFFL = ZOSGX*safe_adb(RootEffLen4Absorption_pvr(N,L),LOG(PATHL/FineRootRadius(N,L)))
  !
  ! NO3 UPTAKE IN NON-BAND SOIL ZONE
  !
  !  VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3S=NO3 concentration in non-band
  !     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFNO3=soil-root convective NO3 flux per plant in non-band
  !     DIFNO3=soil-root NO3 diffusion per plant in non-band
  !
  VmaxNO3RootPerPlant=VmaxNO3Root_pft(N,NZ)*RootSAreaPerPlant_pvr(N,L,NZ) &
      *FSatNutTransporter*fTgrowRootP_vr(L,NZ)*AMIN1(FCUP,FZUP)

  IF(trcs_VLN_vr(ids_NO3,L).GT.ZERO.AND.trc_solcl_vr(ids_NO3,L).GT.CminNO3Root_pft(N,NZ))THEN
    RMFNO3 = PerPlantRootH2OUptake*trcs_VLN_vr(ids_NO3,L)
    DIFNO3 = DIFFL*trcs_VLN_vr(ids_NO3,L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max NO3 uptake in non-band unlimited,limited by O2
    !     RootSAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FSatNutTransporter=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=VmaxNO3RootPerPlant*trcs_VLN_vr(ids_NO3,L)

    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFNO3=soil-root convective N03 flux per plant in non-band
    !     DIFNO3=soil-root N03 diffusion per plant in non-band
    !     CNO3S=NO3 concentration in non-band
    !     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake,Km,min concn from PFT file
    !     RTKNO3,RTKNOP=NO3 uptake per plant in non-band lmtd,unlmtd by O2
    !     ZNO3M,ZNO3X=minimum,maximum NO3 available for uptake in non-band
    !     FNO3X=fraction of total NH4 uptake in non-band by root,myco populn
    !     RootNO3DmndSoil_pvr,RootNO3Uptake_pvr=NO3 uptake in non-band unlimited,limited by NO3
    !     RootOUlmNutUptake_pvr=NO3 uptake in non-band unlimited by O2
    !     RootCUlmNutUptake_pvr=NO3 uptake in non-band unlimited by nonstructural C
!
    ZNO3M = CminNO3Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_NO3,L)
    ZNO3X = AZMAX1(FNO3X*(trcs_solml_vr(ids_NO3,L)-ZNO3M))

    PlantSoluteUptakeConfig%SoluteConcMin   = CminNO3Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFNO3
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFNO3
    PlantSoluteUptakeConfig%UptakePerPlantRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_NO3,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = ZNO3X
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopuLive_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmNO3Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootNO3DmndSoil_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_NO3,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_NO3,N,L,NZ),RootNutUptake_pvr(ids_NO3,N,L,NZ),ldebug)

  ENDIF
  !
  !     NO3 UPTAKE IN BAND SOIL ZONE
  !
  !     VLNO3,VLNOB=fraction of soil volume in NO3 non-band,band
  !     CNO3B=NO3 concentration in band
  !     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFNOB=soil-root convective NO3 flux per plant in band
  !     DIFNOB=soil-root NO3 diffusion per plant in band
  !

  IF(trcs_VLN_vr(ids_NO3B,L).GT.ZERO.AND.trc_solcl_vr(ids_NO3B,L).GT.CminNO3Root_pft(N,NZ))THEN
    RMFNOB = PerPlantRootH2OUptake*trcs_VLN_vr(ids_NO3B,L)
    DIFNOB = DIFFL*trcs_VLN_vr(ids_NO3B,L)
    !
    !     NO3 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum NO3 uptake in band unlimited,limited by O2
    !     RootSAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FSatNutTransporter=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=VmaxNO3RootPerPlant*trcs_VLN_vr(ids_NO3B,L)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF NO3 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF NO3 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFNOB=soil-root convective NO3 flux per plant in band
    !     DIFNOB=soil-root NO3 diffusion per plant in band
    !     CNO3B=NH4 concentration in band
    !     VmaxNO3Root_pft,KmNO3Root_pft,CminNO3Root_pft=NO3 max uptake,Km,min concn from PFT file
    !     RTKNOB,RTKNPB=NO3 uptake per plant in band lmtd,unlmtd by O2
    !     ZNOBM,ZNOBX=minimum,maximum NO3 available for uptake in band
    !     FNOBX=fraction of total NO3 uptake in band by root,myco populn
    !     RootNO3DmndBand_pvr,RootNO3BUptake_pvr=NO3 uptake in band unlimited,limited by NH4
    !     RootOUlmNutUptake_pvr=NO3 uptake in band unlimited by O2
    !     RootCUlmNutUptake_pvr=NO3 uptake in band unlimited by nonstructural C
    !
    ZNOBM = CminNO3Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_NO3B,L)
    ZNOBX = AZMAX1(FNOBX*(trcs_solml_vr(ids_NO3B,L)-ZNOBM))

    PlantSoluteUptakeConfig%SoluteConcMin         = CminNO3Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx             = RMFNOB
    PlantSoluteUptakeConfig%SolDifusFlx           = DIFNOB
    PlantSoluteUptakeConfig%UptakePerPlantRateMax = UPMXP
    PlantSoluteUptakeConfig%O2Stress              = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc            = trc_solcl_vr(ids_NO3B,L)
    PlantSoluteUptakeConfig%SoluteMassMax         = ZNOBX
    PlantSoluteUptakeConfig%CAvailStress          = FCUP
    PlantSoluteUptakeConfig%PlantPopulation       = PlantPopuLive_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM              = KmNO3Root_pft(N,NZ)
  
    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootNO3DmndBand_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_NO3B,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_NO3B,N,L,NZ),RootNutUptake_pvr(ids_NO3B,N,L,NZ))

  ENDIF
  ! Normalize population capacity by population root/mycorrhizal carbon.
  if(RootMyMassC>1.e-4_r8)then
    VmaxNO3Root_pvr(N,L,NZ)=VmaxNO3RootPerPlant*PlantPopuLive_pft(NZ)/RootMyMassC
  else
    VmaxNO3Root_pvr(N,L,NZ)=0._r8
  endif
  call PrintInfo('end '//subname)
  end associate
  end subroutine UptakeNO3

!----------------------------------------------------------------------------------------------------
  subroutine UptakeNH4(I,J,N,L,NZ,FNH4X,FNHBX,PathLen_pvr,FineRootRadius,RootEffLen4Absorption_pvr,&
    FCUP,FZUP,FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)

  implicit none
  integer , intent(in) :: I,J
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: FNH4X,FNHBX,PathLen_pvr(pltpar%jroots,JZ1),FineRootRadius(pltpar%jroots,JZ1)
  real(r8), intent(in) :: RootEffLen4Absorption_pvr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: FCUP
  real(r8), intent(in) :: FZUP
  real(r8), intent(in) :: FSatNutTransporter
  real(r8), intent(in) :: PerPlantRootH2OUptake
  real(r8), intent(in) :: RootMyMassC !population-total root/mycorrhizal carbon in this layer, [gC d-2]
  real(r8) :: DIFFL
  real(r8) :: DIFNH4,DIFNHB
  real(r8) :: PATHL
  real(r8) :: RMFNH4,RTKNH4,RTKNHP,RMFNHB,RTKNHB,RTKNBP
  real(r8) :: UPMX,UPMXP
  real(r8) :: VmaxNH4RootPerPlant !uptake capacity before O2 limitation, [gN plant-1 h-1]
  real(r8) :: ZNHBX,ZNSGX,ZNH4M,ZNH4X,ZNHBM
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig
! begin_execution
  associate(                                                    &
    PlantPopuLive_pft      => plt_site%PlantPopuLive_pft       ,& !input  :plant population, [d-2]
    ZERO                   => plt_site%ZERO                    ,& !input  :threshold zero for numerical stability, [-]
    TortMicPM_vr           => plt_site%TortMicPM_vr            ,& !input  :micropore soil tortuosity, [m3 m-3]
    fTgrowRootP_vr         => plt_pheno%fTgrowRootP_vr         ,& !input  :root layer temperature growth functiom, [-]
    RAutoRootO2Limter_rpvr => plt_rbgc%RAutoRootO2Limter_rpvr  ,& !input  :O2 constraint to root respiration (0-1), [-]
    CMinNH4Root_pft        => plt_rbgc%CMinNH4Root_pft         ,& !input  :minimum NH4 concentration for root NH4 uptake, [g m-3]
    VmaxNH4Root_pft        => plt_rbgc%VmaxNH4Root_pft         ,& !input  :maximum root NH4 uptake rate, [g m-2 h-1]
    KmNH4Root_pft          => plt_rbgc%KmNH4Root_pft           ,& !input  :Km for root NH4 uptake, [g m-3]
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr       ,& !input  :root uptake of Nutrient band, [g d-2 h-1]
    RootCUlmNutUptake_pvr  => plt_rbgc%RootCUlmNutUptake_pvr   ,& !input  :root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
    RootOUlmNutUptake_pvr  => plt_rbgc%RootOUlmNutUptake_pvr   ,& !input  :root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
    RootNH4DmndSoil_pvr    => plt_rbgc%RootNH4DmndSoil_pvr     ,& !input  :root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
    RootNH4DmndBand_pvr    => plt_rbgc%RootNH4DmndBand_pvr     ,& !input  :root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
    RootSAreaPerPlant_pvr  => plt_morph%RootSAreaPerPlant_pvr  ,& !input  :root layer area per plant, [m p-1]
    SoluteDifusvtyT_vr     => plt_soilchem%SoluteDifusvtyT_vr  ,& !input  :aqueous diffusivity, [m2 h-1]
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr        ,& !input  :soil micropore water content, [m3 d-2]
    trcs_VLN_vr            => plt_soilchem%trcs_VLN_vr         ,& !input  :effective relative tracer volume, [-]
    trcs_solml_vr          => plt_soilchem%trcs_solml_vr       ,& !input  :aqueous tracer, [g d-2]
    trc_solcl_vr           => plt_soilchem%trc_solcl_vr        ,& !input  :aqueous tracer concentration, [g m-3]
    VmaxNH4Root_pvr        => plt_rbgc%VmaxNH4Root_pvr          & !output :maximum NH4 uptake rate, [gN h-1 (gC root)-1]
  )
! ZNSGL=NH4 diffusivity
! TortMicPM_vr=soil tortuosity
! PATH=path length of water and nutrient uptake
! FineRootRadius=root radius
! DIFFL=NH4 diffusion per plant
!
  ZNSGX = SoluteDifusvtyT_vr(idg_NH3,L)*TortMicPM_vr(NPH,L)
  PATHL = AMIN1(PathLen_pvr(N,L),FineRootRadius(N,L)+SQRT(2.0*ZNSGX))
  DIFFL = ZNSGX*safe_adb(RootEffLen4Absorption_pvr(N,L),LOG(PATHL/FineRootRadius(N,L)))
!
! NH4 UPTAKE IN NON-BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4S=NH4 concentration in non-band
! VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake,Km,min concn from PFT file
! PerPlantRootH2OUptake=root water uptake per plant
! RMFNH4=soil-root convective NH4 flux per plant in non-band
! DIFNH4=soil-root NH4 diffusion per plant in non-band
!
  VmaxNH4RootPerPlant=VmaxNH4Root_pft(N,NZ)*RootSAreaPerPlant_pvr(N,L,NZ) &
      *FSatNutTransporter*fTgrowRootP_vr(L,NZ)*AMIN1(FCUP,FZUP)
!  if(L<=3)write(1002,*)I*100+J,L,NZ,VmaxNH4RootPerPlant,VmaxNH4Root_pft(N,NZ),RootSAreaPerPlant_pvr(N,L,NZ), &
!      FSatNutTransporter,fTgrowRootP_vr(L,NZ),AMIN1(FCUP,FZUP)
  
  IF(trcs_VLN_vr(ids_NH4,L).GT.ZERO.AND.trc_solcl_vr(ids_NH4,L).GT.CMinNH4Root_pft(N,NZ))THEN
    RMFNH4=PerPlantRootH2OUptake*trcs_VLN_vr(ids_NH4,L)
    DIFNH4=DIFFL*trcs_VLN_vr(ids_NH4,L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=max NH4 uptake in non-band unlimited,limited by O2
!   RootSAreaPerPlant_pvr=root surface area per plant from grosub.f
!   FSatNutTransporter=protein concentration relative to 5%
!   fTgrowRootP_vr=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
!
    UPMXP=VmaxNH4RootPerPlant*trcs_VLN_vr(ids_NH4,L)

!
!   SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
!   SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
!   WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!   RMFNH4=soil-root convective NH4 flux per plant in non-band
!   DIFNH4=soil-root NH4 diffusion per plant in non-band
!   CNH4S=NH4 concentration in non-band
!   VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake,Km,min concn from PFT file
!   RTKNH4,RTKNHP=NH4 uptake per plant in non-band lmtd,unlmtd by O2
!   ZNH4M,ZNH4X=minimum,maximum NH4 available for uptake in non-band
!   FNH4X=fraction of total NH4 uptake in non-band by root,myco populn
!   RootNH4DmndSoil_pvr,RootNH4Uptake_pvr=NH4 uptake in non-band unlimited,limited by NH4
!   RootOUlmNutUptake_pvr=NH4 uptake in non-band unlimited by O2
!   RootCUlmNutUptake_pvr=NH4 uptake in non-band unlimited by nonstructural C
!
    ZNH4M = CMinNH4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_NH4,L)
    ZNH4X = AZMAX1(FNH4X*(trcs_solml_vr(ids_NH4,L)-ZNH4M))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinNH4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFNH4
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFNH4
    PlantSoluteUptakeConfig%UptakePerPlantRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_NH4,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = ZNH4X
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopuLive_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmNH4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootNH4DmndSoil_pvr(N,L,NZ),&
      RootOUlmNutUptake_pvr(ids_NH4,N,L,NZ),RootCUlmNutUptake_pvr(ids_NH4,N,L,NZ),&
      RootNutUptake_pvr(ids_NH4,N,L,NZ))

  ENDIF
!
! NH4 UPTAKE IN BAND SOIL ZONE
!
! VLNH4,VLNHB=fraction of soil volume in NH4 non-band,band
! CNH4B=NH4 concentration in band
! VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake,Km,min concn from PFT file
! PerPlantRootH2OUptake=root water uptake per plant
! RMFNHB=soil-root convective NH4 flux per plant in band
! DIFNHB=soil-root NH4 diffusion per plant in band
!

  IF(trcs_VLN_vr(ids_NH4B,L).GT.ZERO.AND.trc_solcl_vr(ids_NH4B,L).GT.CMinNH4Root_pft(N,NZ))THEN
    RMFNHB=PerPlantRootH2OUptake*trcs_VLN_vr(ids_NH4B,L)
    DIFNHB=DIFFL*trcs_VLN_vr(ids_NH4B,L)
!
!   NH4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
!   AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
!
!   UPMXP,UPMX=maximum NH4 uptake in band unlimited,limited by O2
!   RootSAreaPerPlant_pvr=root surface area per plant from grosub.f
!   FSatNutTransporter=protein concentration relative to 5%
!   fTgrowRootP_vr=temperature function for root growth
!   FCUP,FZUP=limitn to active uptake respiration from CCPOLR,CZPOLR
!   RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
!
    UPMXP=VmaxNH4RootPerPlant*trcs_VLN_vr(ids_NH4B,L)
!
!   SOLUTION FOR MASS FLOW + DIFFUSION OF NH4 IN AQUEOUS PHASE OF
!   SOIL = ACTIVE UPTAKE OF NH4 BY ROOT, CONSTRAINED BY COMPETITION
!   WITH OTHER ROOT AND MICROBIAL POPULATIONS
!
!   RMFNHB=soil-root convective NH4 flux per plant in band
!   DIFNHB=soil-root NH4 diffusion per plant in band
!   CNH4B=NH4 concentration in band
!   VmaxNH4Root_pft,KmNH4Root_pft,CMinNH4Root_pft=NH4 max uptake,Km,min concn from PFT file
!   RTKNHB,RTKNBP=NH4 uptake per plant in band lmtd,unlmtd by O2
!   ZNHBM,ZNHBX=minimum,maximum NH4 available for uptake in band
!   FNHBX=fraction of total NH4 uptake in band by root,myco populn
!   RootNH4DmndBand_pvr,RootNH4BUptake_pvr=NH4 uptake in band unlimited,limited by NH4
!   RootOUlmNutUptake_pvr=NH4 uptake in band unlimited by O2
!   RootCUlmNutUptake_pvr=NH4 uptake in band unlimited by nonstructural C
!

    ZNHBM=CMinNH4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_NH4B,L)
    ZNHBX=AZMAX1(FNHBX*(trcs_solml_vr(ids_NH4B,L)-ZNHBM))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinNH4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFNHB
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFNHB
    PlantSoluteUptakeConfig%UptakePerPlantRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_NH4B,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = ZNHBX
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopuLive_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmNH4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootNH4DmndBand_pvr(N,L,NZ),&
      RootOUlmNutUptake_pvr(ids_NH4B,N,L,NZ),RootCUlmNutUptake_pvr(ids_NH4B,N,L,NZ),&
      RootNutUptake_pvr(ids_NH4B,N,L,NZ))

  ENDIF
  ! Normalize population capacity by population root/mycorrhizal carbon.
  if(RootMyMassC>1.e-4_r8)then
    VmaxNH4Root_pvr(N,L,NZ)=VmaxNH4RootPerPlant*PlantPopuLive_pft(NZ)/RootMyMassC
  else
    VmaxNH4Root_pvr(N,L,NZ)=0._r8
  endif
  end associate
  end subroutine UptakeNH4

!----------------------------------------------------------------------------------------------------
  subroutine UptakeHPO4(N,L,NZ,DIFFL,FP14X,FP1BX,FCUP,FPUP,FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)

  implicit none
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: DIFFL
  real(r8), intent(in) :: FP14X,FP1BX
  real(r8), intent(in) :: FCUP,FPUP,FSatNutTransporter,PerPlantRootH2OUptake
  real(r8), intent(in) :: RootMyMassC
  real(r8) :: B,C
  real(r8) :: BP,CP
  real(r8) :: DIFH1P
  real(r8) :: DIFH1B
  real(r8) :: H1POM,H1POX,H1PXM,H1PXB
  real(r8) :: RMFH1P,RTKH1P,RTKHP1,RTKH1B,RTKHB1
  real(r8) :: RMFH2B
  real(r8) :: UPMX,UPMXP
  real(r8) :: X,Y
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig
  !     begin_execution
  associate(                                                    &
    PlantPopuLive_pft      => plt_site%PlantPopuLive_pft       ,& !input  :plant population, [d-2]
    ZERO                   => plt_site%ZERO                    ,& !input  :threshold zero for numerical stability, [-]
    fTgrowRootP_vr         => plt_pheno%fTgrowRootP_vr         ,& !input  :root layer temperature growth functiom, [-]
    RAutoRootO2Limter_rpvr => plt_rbgc%RAutoRootO2Limter_rpvr  ,& !input  :O2 constraint to root respiration (0-1), [-]
    CMinPO4Root_pft        => plt_rbgc%CMinPO4Root_pft         ,& !input  :minimum PO4 concentration for root NH4 uptake, [g m-3]
    KmPO4Root_pft          => plt_rbgc%KmPO4Root_pft           ,& !input  :Km for root PO4 uptake, [g m-3]
    VmaxPO4Root_pft        => plt_rbgc%VmaxPO4Root_pft         ,& !input  :maximum root PO4 uptake rate, [g m-2 h-1]
    RootCUlmNutUptake_pvr  => plt_rbgc%RootCUlmNutUptake_pvr   ,& !input  :root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
    RootH1PO4DmndSoil_pvr  => plt_rbgc%RootH1PO4DmndSoil_pvr   ,& !input  :HPO4 demand in non-band by each root population, [g d-2 h-1]
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr       ,& !input  :root uptake of Nutrient band, [g d-2 h-1]
    RootOUlmNutUptake_pvr  => plt_rbgc%RootOUlmNutUptake_pvr   ,& !input  :root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
    RootH1PO4DmndBand_pvr  => plt_rbgc%RootH1PO4DmndBand_pvr   ,& !input  :HPO4 demand in band by each root population, [g d-2 h-1]
    RootSAreaPerPlant_pvr  => plt_morph%RootSAreaPerPlant_pvr  ,& !input  :root layer area per plant, [m p-1]
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr        ,& !input  :soil micropore water content, [m3 d-2]
    trcs_VLN_vr            => plt_soilchem%trcs_VLN_vr         ,& !input  :effective relative tracer volume, [-]
    trc_solcl_vr           => plt_soilchem%trc_solcl_vr        ,& !input  :aqueous tracer concentration, [g m-3]
    trcs_solml_vr          => plt_soilchem%trcs_solml_vr        & !input  :aqueous tracer, [g d-2]
  )
  !
  !     HPO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4=HPO4 concentration in non-band
  !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant > 0._r8
  !     RMFH1P=soil-root convective HPO4 flux per plant in non-band
  !     DIFH1P=soil-root HPO4 diffusion per plant in non-band > 0._r8
!
    IF(trcs_VLN_vr(ids_H1PO4,L).GT.ZERO.AND.trc_solcl_vr(ids_H1PO4,L).GT.CMinPO4Root_pft(N,NZ))THEN
      RMFH1P = PerPlantRootH2OUptake*trcs_VLN_vr(ids_H1PO4,L)
      DIFH1P = DIFFL*trcs_VLN_vr(ids_H1PO4,L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=max HPO4 uptake in non-band unlimited,limited by O2
    !     RootSAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FSatNutTransporter=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1_r8*VmaxPO4Root_pft(N,NZ)*RootSAreaPerPlant_pvr(N,L,NZ) &
      *FSatNutTransporter*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_H1PO4,L)*AMIN1(FCUP,FPUP)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH1P=soil-root convective HPO4 flux per plant in non-band
    !     DIFH1P=soil-root HPO4 diffusion per plant in non-band
    !     CH1P4=HPO4 concentration in non-band
    !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH1P,RTKHP1=HPO4 uptake per plant in non-band lmtd,unlmtd by O2
    !     H1POM,H1POX=minimum,maximum HPO4 available for uptake in non-band
    !     FP14X=fraction of total HPO4 uptake in non-band by root,myco populn
    !     RootH1PO4DmndSoil_pvr,RootNutUptake_pvr=HPO4 uptake in non-band unlimited,limited by HPO4
    !     RootOUlmNutUptake_pvr(ids_H1PO4,=HPO4 uptake in non-band unlimited by O2
    !     RootCUlmNutUptake_pvr=HPO4 uptake in non-band unlimited by nonstructural C

    H1POM=CMinPO4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_H1PO4,L)
    H1POX=AZMAX1(FP14X*(trcs_solml_vr(ids_H1PO4,L)-H1POM))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinPO4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFH1P
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFH1P
    PlantSoluteUptakeConfig%UptakePerPlantRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_H1PO4,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = H1POX
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopuLive_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmPO4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootH1PO4DmndSoil_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_H1PO4,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_H1PO4,N,L,NZ),RootNutUptake_pvr(ids_H1PO4,N,L,NZ))

  ENDIF
  !
  !     HPO4 UPTAKE IN BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH1P4B=HPO4 concentration in band
  !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFH1B=soil-root convective HPO4 flux per plant in band
  !     DIFH1B=soil-root HPO4 diffusion per plant in band
  !
  IF(trcs_VLN_vr(ids_H1PO4B,L).GT.ZERO.AND.trc_solcl_vr(ids_H1PO4B,L).GT.CMinPO4Root_pft(N,NZ))THEN
    RMFH2B = PerPlantRootH2OUptake*trcs_VLN_vr(ids_H1PO4B,L)
    DIFH1B = DIFFL*trcs_VLN_vr(ids_H1PO4B,L)
    !
    !     HPO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum HPO4 uptake in band unlimited,limited by O2
    !     RootSAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FSatNutTransporter=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=0.1_r8*VmaxPO4Root_pft(N,NZ)*RootSAreaPerPlant_pvr(N,L,NZ) &
      *FSatNutTransporter*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_H1PO4B,L)*AMIN1(FCUP,FPUP)
    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF HPO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF HPO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH1B=soil-root convective HPO4 flux per plant in band
    !     DIFH1B=soil-root HPO4 diffusion per plant in band
    !     CH1P4B=HPO4 concentration in band
    !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH1B,RTKHB1=HPO4 uptake per plant in band lmtd,unlmtd by O2
    !     H1PXM,H1PXB=minimum,maximum HPO4 available for uptake in band
    !     FP1BX=fraction of total HPO4 uptake in band by root,myco populn
    !     RootH1PO4DmndBand_pvr,RootH1PO4BUptake_pvr=HPO4 uptake in band unlimited,limited by H2PO4
    !     RootOUlmNutUptake_pvr=HPO4 uptake in band unlimited by O2
    !     RootCUlmNutUptake_pvr=HPO4 uptake in band unlimited by nonstructural C

    H1PXM = CMinPO4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_H1PO4B,L)
    H1PXB = AZMAX1(FP1BX*(trcs_solml_vr(ids_H1PO4B,L)-H1PXM))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinPO4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFH2B
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFH1B
    PlantSoluteUptakeConfig%UptakePerPlantRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_H1PO4B,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = H1PXB
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopuLive_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmPO4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootH1PO4DmndBand_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_H1PO4B,N,L,NZ),RootNutUptake_pvr(ids_H1PO4B,N,L,NZ))

  ENDIF
  end associate
  end subroutine UptakeHPO4

!----------------------------------------------------------------------------------------------------
  subroutine UptakeH2PO4(N,L,NZ,DIFFL,FPO4X,FPOBX,FCUP,FPUP,FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)

  implicit none
  integer,  intent(in) :: N,L
  integer,  intent(in) :: NZ
  real(r8), intent(in) :: DIFFL
  real(r8), intent(in) :: FPO4X,FPOBX
  real(r8), intent(in) :: FCUP
  real(r8), intent(in) :: FPUP
  real(r8), intent(in) :: FSatNutTransporter
  real(r8), intent(in) :: PerPlantRootH2OUptake
  real(r8), intent(in) :: RootMyMassC

  real(r8) :: DIFH2P,DIFH2B
  real(r8) :: H2POM,H2POX,H2PXM,H2PXB
  real(r8) :: RTKHPB,RMFH2B,RMFH2P,RTKH2P,RTKHPP,RTKH2B
  real(r8) :: UPMX,UPMXP
  type(PlantSoluteUptakeConfig_type) :: PlantSoluteUptakeConfig
  !
  associate(                                                    &
    PlantPopuLive_pft      => plt_site%PlantPopuLive_pft       ,& !input  :plant population, [d-2]
    ZERO                   => plt_site%ZERO                    ,& !input  :threshold zero for numerical stability, [-]
    RAutoRootO2Limter_rpvr => plt_rbgc%RAutoRootO2Limter_rpvr  ,& !input  :O2 constraint to root respiration (0-1), [-]
    KmPO4Root_pft          => plt_rbgc%KmPO4Root_pft           ,& !input  :Km for root PO4 uptake, [g m-3]
    VmaxPO4Root_pft        => plt_rbgc%VmaxPO4Root_pft         ,& !input  :maximum root PO4 uptake rate, [g m-2 h-1]
    CMinPO4Root_pft        => plt_rbgc%CMinPO4Root_pft         ,& !input  :minimum PO4 concentration for root NH4 uptake, [g m-3]
    RootCUlmNutUptake_pvr  => plt_rbgc%RootCUlmNutUptake_pvr   ,& !input  :root uptake of NH4 band unconstrained by root nonstructural C, [g d-2 h-1]
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr       ,& !input  :root uptake of Nutrient band, [g d-2 h-1]
    RootOUlmNutUptake_pvr  => plt_rbgc%RootOUlmNutUptake_pvr   ,& !input  :root uptake of NH4 band unconstrained by O2, [g d-2 h-1]
    RootH2PO4DmndBand_pvr  => plt_rbgc%RootH2PO4DmndBand_pvr   ,& !input  :root uptake of H2PO4 band, [g d-2 h-1]
    RootH2PO4DmndSoil_pvr  => plt_rbgc%RootH2PO4DmndSoil_pvr   ,& !input  :root uptake of H2PO4 non-band, [g d-2 h-1]
    fTgrowRootP_vr         => plt_pheno%fTgrowRootP_vr         ,& !input  :root layer temperature growth functiom, [-]
    RootSAreaPerPlant_pvr  => plt_morph%RootSAreaPerPlant_pvr  ,& !input  :root layer area per plant, [m p-1]
    trcs_solml_vr          => plt_soilchem%trcs_solml_vr       ,& !input  :aqueous tracer, [g d-2]
    VLWatMicP_vr           => plt_soilchem%VLWatMicP_vr        ,& !input  :soil micropore water content, [m3 d-2]
    trc_solcl_vr           => plt_soilchem%trc_solcl_vr        ,& !input  :aqueous tracer concentration, [g m-3]
    trcs_VLN_vr            => plt_soilchem%trcs_VLN_vr          & !input  :effective relative tracer volume, [-]
  )
  !     H2PO4 UPTAKE IN NON-BAND SOIL ZONE
  !
  !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
  !     CH2P4=H2PO4 concentration in non-band
  !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
  !     PerPlantRootH2OUptake=root water uptake per plant
  !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
  !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
!
    IF(trcs_VLN_vr(ids_H1PO4,L).GT.ZERO.AND.trc_solcl_vr(ids_H2PO4,L).GT.CMinPO4Root_pft(N,NZ))THEN
      RMFH2P=PerPlantRootH2OUptake*trcs_VLN_vr(ids_H1PO4,L)
      DIFH2P=DIFFL*trcs_VLN_vr(ids_H1PO4,L)
      !
      !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
      !     AND FROM ROOT SURFACE AREA, C AND P CONSTRAINTS CALCULATED ABOVE
      !
      !     UPMXP,UPMX=max H2PO4 uptake in non-band unlimited,limited by O2
      !     RootSAreaPerPlant_pvr=root surface area per plant from grosub.f
      !     FSatNutTransporter=protein concentration relative to 5%
      !     fTgrowRootP_vr=temperature function for root growth
      !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
      !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
!
      UPMXP=VmaxPO4Root_pft(N,NZ)*RootSAreaPerPlant_pvr(N,L,NZ) &
        *FSatNutTransporter*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_H1PO4,L)*AMIN1(FCUP,FPUP)
      !
      !     SOLUTION FOR MASS FLOW + DIFFUSION OF H2PO4 IN AQUEOUS PHASE OF
      !     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY
      !     COMPETITION WITH OTHER ROOT AND MICROBIAL POPULATIONS
      !
      !     RMFH2P=soil-root convective H2PO4 flux per plant in non-band
      !     DIFH2P=soil-root H2PO4 diffusion per plant in non-band
      !     CH2P4=H2PO4 concentration in non-band
      !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
      !     RTKH2P,RTKHPP=H2PO4 uptake per plant in non-band lmtd,unlmtd by O2
      !     H2POM,H2POX=minimum,maximum H2PO4 available for uptake in non-band
      !     FPO4X=fraction of total H2PO4 uptake in non-band by root,myco populn
      !     RootH2PO4DmndSoil_pvr,RootH2PO4Uptake_pvr=H2PO4 uptake in non-band unlimited,limited by H2PO4
      !     RootOUlmNutUptake_pvr=H2PO4 uptake in non-band unlimited by O2
      !     RootCUlmNutUptake_pvr=H2PO4 uptake in non-band unlimited by nonstructural C

      H2POM = CMinPO4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_H1PO4,L)
      H2POX = AZMAX1(FPO4X*(trcs_solml_vr(ids_H2PO4,L)-H2POM))

      PlantSoluteUptakeConfig%SoluteConcMin   = CMinPO4Root_pft(N,NZ)
      PlantSoluteUptakeConfig%SolAdvFlx       = RMFH2P
      PlantSoluteUptakeConfig%SolDifusFlx     = DIFH2P
      PlantSoluteUptakeConfig%UptakePerPlantRateMax   = UPMXP
      PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
      PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_H2PO4,L)
      PlantSoluteUptakeConfig%SoluteMassMax   = H2POX
      PlantSoluteUptakeConfig%CAvailStress    = FCUP
      PlantSoluteUptakeConfig%PlantPopulation = PlantPopuLive_pft(NZ)
      PlantSoluteUptakeConfig%SoluteKM        = KmPO4Root_pft(N,NZ)

      call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootH2PO4DmndSoil_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_H2PO4,N,L,NZ),&
        RootCUlmNutUptake_pvr(ids_H2PO4,N,L,NZ),RootNutUptake_pvr(ids_H2PO4,N,L,NZ))

    ENDIF
    !
    !     H2PO4 UPTAKE IN BAND SOIL ZONE
    !
    !     VLPO4,VLPOB=fraction of soil volume in H2PO4 non-band,band
    !     CH2P4B=H2PO4 concentration in band
    !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
    !     PerPlantRootH2OUptake=root water uptake per plant
    !     RMFH2B=soil-root convective H2PO4 flux per plant in band
    !     DIFH2B=soil-root H2PO4 diffusion per plant in band
    !

  IF(trcs_VLN_vr(ids_H1PO4B,L).GT.ZERO.AND.trc_solcl_vr(ids_H2PO4B,L).GT.CMinPO4Root_pft(N,NZ))THEN
    RMFH2B=PerPlantRootH2OUptake*trcs_VLN_vr(ids_H1PO4B,L)
    DIFH2B=DIFFL*trcs_VLN_vr(ids_H1PO4B,L)
    !
    !     H2PO4 UPTAKE DEMAND FROM ROOT UPTAKE PARAMETERS ENTERED IN 'READQ'
    !     AND FROM ROOT SURFACE AREA, C AND N CONSTRAINTS CALCULATED ABOVE
    !
    !     UPMXP,UPMX=maximum H2PO4 uptake in band unlimited,limited by O2
    !     RootSAreaPerPlant_pvr=root surface area per plant from grosub.f
    !     FSatNutTransporter=protein concentration relative to 5%
    !     fTgrowRootP_vr=temperature function for root growth
    !     FCUP,FPUP=limitn to active uptake respiration from CCPOLR,CPPOLR
    !     RAutoRootO2Limter_rpvr=constraint by O2 consumption on all biological processes
    !
    UPMXP=VmaxPO4Root_pft(N,NZ)*RootSAreaPerPlant_pvr(N,L,NZ) &
      *FSatNutTransporter*fTgrowRootP_vr(L,NZ)*trcs_VLN_vr(ids_H1PO4B,L)*AMIN1(FCUP,FPUP)

    !
    !     SOLUTION FOR MASS FLOW + DIFFUSION OF PO4 IN AQUEOUS PHASE OF
    !     SOIL = ACTIVE UPTAKE OF H2PO4 BY ROOT, CONSTRAINED BY COMPETITION
    !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
    !
    !     RMFH2B=soil-root convective H2PO4 flux per plant in band
    !     DIFH2B=soil-root H2PO4 diffusion per plant in band
    !     CH2P4B=H2PO4 concentration in band
    !     VmaxPO4Root_pft,KmPO4Root_pft,CMinPO4Root_pft=H2PO4 max uptake,Km,min concn from PFT file
    !     RTKH2B,RTKHPB=H2PO4 uptake per plant in band lmtd,unlmtd by O2
    !     H2PXM,H2PXB=minimum,maximum H2PO4 available for uptake in band
    !     FPOBX=fraction of total H2PO4 uptake in band by root,myco populn
    !     RootH2PO4DmndBand_pvr,RootNutUptake_pvr=H2PO4 uptake in band unlimited,limited by H2PO4
    !     RootOUlmNutUptake_pvr=H2PO4 uptake in band unlimited by O2
    !     RootCUlmNutUptake_pvr=H2PO4 uptake in band unlimited by nonstructural C

    H2PXM=CMinPO4Root_pft(N,NZ)*VLWatMicP_vr(L)*trcs_VLN_vr(ids_H1PO4B,L)
    H2PXB=AZMAX1(FPOBX*(trcs_solml_vr(ids_H2PO4B,L)-H2PXM))

    PlantSoluteUptakeConfig%SoluteConcMin   = CMinPO4Root_pft(N,NZ)
    PlantSoluteUptakeConfig%SolAdvFlx       = RMFH2B
    PlantSoluteUptakeConfig%SolDifusFlx     = DIFH2B
    PlantSoluteUptakeConfig%UptakePerPlantRateMax   = UPMXP
    PlantSoluteUptakeConfig%O2Stress        = RAutoRootO2Limter_rpvr(N,L,NZ)
    PlantSoluteUptakeConfig%SoluteConc      = trc_solcl_vr(ids_H2PO4B,L)
    PlantSoluteUptakeConfig%SoluteMassMax   = H2PXB
    PlantSoluteUptakeConfig%CAvailStress    = FCUP
    PlantSoluteUptakeConfig%PlantPopulation = PlantPopuLive_pft(NZ)
    PlantSoluteUptakeConfig%SoluteKM        = KmPO4Root_pft(N,NZ)

    call SoluteUptakeByPlantRoots(PlantSoluteUptakeConfig,RootH2PO4DmndBand_pvr(N,L,NZ),RootOUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ),&
      RootCUlmNutUptake_pvr(ids_H2PO4B,N,L,NZ),RootNutUptake_pvr(ids_H2PO4B,N,L,NZ))

  ENDIF
  end associate
  end subroutine UptakeH2PO4

!----------------------------------------------------------------------------------------------------
  subroutine UptakeMineralNitrogen(I,J,N,L,NZ,PathLen_pvr,FineRootRadius,FracPRoot4Uptake,&
    FracMinRoot4Uptake_rpvr,RootEffLen4Absorption_pvr,FCUP,FZUP,FSatNutTransporter,PerPlantRootH2OUptake,RootMyMassC)

  implicit none
  integer , intent(in) :: I,J
  integer , intent(in) :: N,L
  integer , intent(in) :: NZ
  real(r8), intent(in) :: PathLen_pvr(pltpar%jroots,JZ1),FineRootRadius(pltpar%jroots,JZ1),FracPRoot4Uptake(pltpar%jroots,JZ1,JP1)
  real(r8), intent(in) :: FracMinRoot4Uptake_rpvr(pltpar%jroots,JZ1,JP1)
  real(r8), intent(in) :: RootEffLen4Absorption_pvr(pltpar%jroots,JZ1)
  real(r8), intent(in) :: FCUP
  real(r8), intent(in) :: FZUP
  real(r8), intent(in) :: FSatNutTransporter
  real(r8), intent(in) :: PerPlantRootH2OUptake
  real(r8), intent(in) :: RootMyMassC
  real(r8) :: FNO3X,FNOBX,FNH4X,FNHBX
  real(r8) :: TFNH4X,TFNO3X,TFNHBX,TFNOBX

!     begin_execution
  associate(                                                            &
    ZERO2                      => plt_site%ZERO2                       ,& !input  :threshold zero for numerical stability,[-]
    ZEROS                      => plt_site%ZEROS                       ,& !input  :threshold zero for numerical stability,[-]
    RootNO3DmndBandPrev_pvr    => plt_rbgc%RootNO3DmndBandPrev_pvr     ,& !input  :root uptake of NO3 non-band unconstrained by NO3, [g d-2 h-1]
    RootNH4DmndSoilPrev_pvr    => plt_rbgc%RootNH4DmndSoilPrev_pvr     ,& !input  :root uptake of NH4 non-band unconstrained by NH4, [g d-2 h-1]
    RootNH4DmndBandPrev_pvr    => plt_rbgc%RootNH4DmndBandPrev_pvr     ,& !input  :root uptake of NO3 band unconstrained by NO3, [g d-2 h-1]
    RootNO3DmndSoilPrev_pvr    => plt_rbgc%RootNO3DmndSoilPrev_pvr     ,& !input  :root uptake of NH4 band unconstrained by NH4, [g d-2 h-1]
    RNO3EcoDmndSoilPrev_vr     => plt_bgcr%RNO3EcoDmndSoilPrev_vr      ,& !input  :total root + microbial NO3 uptake non-band, [gN d-2 h-1]
    RNH4EcoDmndBandPrev_vr     => plt_bgcr%RNH4EcoDmndBandPrev_vr      ,& !input  :total root + microbial NH4 uptake band, [gN d-2 h-1]
    RNH4EcoDmndSoilPrev_vr     => plt_bgcr%RNH4EcoDmndSoilPrev_vr      ,& !input  :total root + microbial NH4 uptake non-band, [gN d-2 h-1]
    RNO3EcoDmndBandPrev_vr     => plt_bgcr%RNO3EcoDmndBandPrev_vr       & !input  :total root + microbial NO3 uptake band, [gN d-2 h-1]
  )
  TFNH4X=0.0_r8
  TFNHBX=0.0_r8
  TFNO3X=0.0_r8
  TFNOBX=0.0_r8

  IF(RNH4EcoDmndSoilPrev_vr(L).GT.ZEROS)THEN
    FNH4X=AMAX1(FracMinRoot4Uptake_rpvr(N,L,NZ),RootNH4DmndSoilPrev_pvr(N,L,NZ)/RNH4EcoDmndSoilPrev_vr(L))
  ELSE
    FNH4X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RNH4EcoDmndBandPrev_vr(L).GT.ZEROS)THEN
    FNHBX=AMAX1(FracMinRoot4Uptake_rpvr(N,L,NZ),RootNH4DmndBandPrev_pvr(N,L,NZ)/RNH4EcoDmndBandPrev_vr(L))
  ELSE
    FNHBX=FracPRoot4Uptake(N,L,NZ)
  ENDIF

  IF(RNO3EcoDmndSoilPrev_vr(L).GT.ZEROS)THEN
    FNO3X=AMAX1(FracMinRoot4Uptake_rpvr(N,L,NZ),RootNO3DmndSoilPrev_pvr(N,L,NZ)/RNO3EcoDmndSoilPrev_vr(L))
  ELSE
    FNO3X=FracPRoot4Uptake(N,L,NZ)
  ENDIF
  IF(RNO3EcoDmndBandPrev_vr(L).GT.ZEROS)THEN
    FNOBX=AMAX1(FracMinRoot4Uptake_rpvr(N,L,NZ),RootNO3DmndBandPrev_pvr(N,L,NZ)/RNO3EcoDmndBandPrev_vr(L))
  ELSE
    FNOBX=FracPRoot4Uptake(N,L,NZ)
  ENDIF

  TFNH4X=TFNH4X+FNH4X
  TFNO3X=TFNO3X+FNO3X
  TFNHBX=TFNHBX+FNHBX
  TFNOBX=TFNOBX+FNOBX
  
  IF(FZUP.GT.ZERO2)THEN
    !
    !     PARAMETERS FOR RADIAL MASS FLOW AND DIFFUSION OF NH4,NO3
    !     FROM SOIL TO ROOT
    !
    call UptakeNH4(I,J,N,L,NZ,FNH4X,FNHBX,PathLen_pvr,FineRootRadius,RootEffLen4Absorption_pvr,FCUP,FZUP,FSatNutTransporter,&
      PerPlantRootH2OUptake,RootMyMassC)

    call UptakeNO3(I,J,N,L,NZ,FNO3X,FNOBX,PathLen_pvr,FineRootRadius,RootEffLen4Absorption_pvr,FCUP,FZUP,FSatNutTransporter,&
      PerPlantRootH2OUptake,RootMyMassC)

  ENDIF
  end associate
  end subroutine UptakeMineralNitrogen

!----------------------------------------------------------------------------------------------------
  subroutine GetUptakeCapcity(I,J,N,L,NZ,FracPRoot4Uptake,FracMinRoot4Uptake_rpvr,FCUP,FZUP,FPUP,&
    FSatNutTransporter,PerPlantRootH2OUptake,dtPerPlantRootH2OUptake,FOXYX)
  !
  !Description:
  !calculate limitation factors for nutrient and gas uptake
  implicit none
  integer, intent(in)   :: I,J  
  integer, intent(in)   :: N,L
  integer, intent(in)   :: NZ
  REAL(R8), INTENT(IN)  :: FracPRoot4Uptake(pltpar%jroots,JZ1,JP1)
  real(r8), intent(in)  :: FracMinRoot4Uptake_rpvr(pltpar%jroots,JZ1,JP1)
  real(r8), intent(out) :: FCUP   !carbon-limitation factor for nutrient uptake, [0-1], greater value weaker limitation
  real(r8), intent(out) :: FZUP   !nitrogen-limitation factor for nutrient uptake, [0-1], greater value stronger limitation
  real(r8), intent(out) :: FPUP   !nitrogen-limitation factor for nutrient uptake, [0-1], greater value stronger limitation
  real(r8), intent(out) :: FSatNutTransporter  !normalized transporter density by reference to 5%
  real(r8), intent(out) :: PerPlantRootH2OUptake    !
  real(r8), intent(out) :: dtPerPlantRootH2OUptake  !per plant root water uptake within the time step of gas uptake, used for transpiration induced gas flow into roots.
  real(r8), intent(out) :: FOXYX  !effort fraction of plant NZ in total plant O2 demand.

  character(len=*), parameter :: subname='GetUptakeCapcity'
  real(r8) :: dss
  associate(                                                                   &
    RO2EcoDmndPrev_vr             => plt_bgcr%RO2EcoDmndPrev_vr               ,& !input  :total root + microbial O2 uptake, [g d-2 h-1]
    RootCO2EmisPot_pvr            => plt_rbgc%RootCO2EmisPot_pvr              ,& !input  :root CO2 efflux unconstrained by root nonstructural C, [g d-2 h-1]
    RootO2Dmnd4Resp_pvr           => plt_rbgc%RootO2Dmnd4Resp_pvr             ,& !input  :root O2 demand from respiration, [g d-2 h-1]
    PlantPopuLive_pft             => plt_site%PlantPopuLive_pft               ,& !input  :plant population, [d-2]
    ZEROS                         => plt_site%ZEROS                           ,& !input  :threshold zero for numerical stability,[-]
    ZERO                          => plt_site%ZERO                            ,& !input  :threshold zero for numerical stability, [-]
    RPlantRootH2OUptk_pvr         => plt_ew%RPlantRootH2OUptk_pvr             ,& !input  :root water uptake, [m2 d-2 h-1]
    RootProteinCMax_pft           => plt_allom%RootProteinCMax_pft            ,& !input  :reference root protein N, [gN g-1]
    ZERO4Groth_pft                => plt_biom%ZERO4Groth_pft                  ,& !input  :threshold zero for plang growth calculation, [-]
    RootProteinC_pvr              => plt_biom%RootProteinC_pvr                ,& !input  :root layer protein C, [gC d-2]
    RootMycoNonstElms_rpvr        => plt_biom%RootMycoNonstElms_rpvr          ,& !input  :root layer nonstructural element, [g d-2]
    RootNonstructElmConc_pvr     => plt_biom%RootNonstructElmConc_pvr       ,& !input  :root layer nonstructural C concentration, [g g-1]
    RootMycoActiveBiomC_pvr       => plt_biom%RootMycoActiveBiomC_pvr         ,& !input  :root layer structural C, [gC d-2]
    Nutruptk_fClim_rpvr           => plt_bgcr%Nutruptk_fClim_rpvr             ,& !output :carbon limitation for nutrient uptake,(0->1),stronger limitation, [-]
    Nutruptk_fNlim_rpvr           => plt_bgcr%Nutruptk_fNlim_rpvr             ,& !output :nitrogen limitation for nutrient uptake,(0->1),stronger limitation, [-]
    Nutruptk_fPlim_rpvr           => plt_bgcr%Nutruptk_fPlim_rpvr             ,& !output :phosphorus limitation for nutrient uptake,(0->1),stronger limitation, [-]
    Nutruptk_fProtC_rpvr          => plt_bgcr%Nutruptk_fProtC_rpvr            ,& !output :transporter scalar indicated by protein for nutrient uptake, greater value greater capacity, [-] 
    RootProteinConc_rpvr          => plt_biom%RootProteinConc_rpvr             & !output :root layer protein C concentration, [g g-1]
  )
  !
  !     UPTAKE CAPACITY 'FSatNutTransporter' DEPENDS ON ROOT,MYCORRHIZAL
  !     PROTEIN CONTENT RELATIVE TO 5% FOR WHICH ACTIVE UPTAKE
  !     PARAMETERS ARE DEFINED
  !
  !     RootProteinConc_rpvr,RootProteinCMax_pft=current,maximum protein concentration
  !     RootProteinC_pvr,WTRTL=protein content,mass
  !     FSatNutTransporter=protein concentration relative to 5%
  !
  call PrintInfo('beg '//subname)
  IF(RootMycoActiveBiomC_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
    RootProteinConc_rpvr(N,L,NZ) = AMIN1(RootProteinCMax_pft(NZ),RootProteinC_pvr(N,L,NZ)/RootMycoActiveBiomC_pvr(N,L,NZ))
    FSatNutTransporter           = RootProteinConc_rpvr(N,L,NZ)/0.05_r8
  ELSE
    RootProteinConc_rpvr(N,L,NZ)=RootProteinCMax_pft(NZ)
    FSatNutTransporter  =1.0_r8
  ENDIF
  Nutruptk_fProtC_rpvr(N,L,NZ)=FSatNutTransporter
  !
  !     RESPIRATION CONSTRAINT ON UPTAKE FROM NON-STRUCTURAL C
  !
  !     RootCO2EmisPot_pvr=total respiration from CPOOLR
  !     FCUP=measures C-limitation for nutrient uptake, using potential respiration as metric
  !     CPOOLR=nonstructural C content
  !C demand
  IF(RootCO2EmisPot_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
    FCUP=AZMAX1(AMIN1(1.0_r8,0.25_r8*RootMycoNonstElms_rpvr(ielmc,N,L,NZ)/RootCO2EmisPot_pvr(N,L,NZ)))
  ELSE
    FCUP=0.0_r8
  ENDIF
  Nutruptk_fClim_rpvr(N,L,NZ) = FCUP
  !
  !     FEEDBACK CONSTRAINT ON N UPTAKE FROM NON-STRUCTURAL N AND P
  !
  !     FZUP,FPUP=limitn to active uptake respiration from CZPOLR,CPPOLR
  !     CCPOLR,CZPOLR,CPPOLR=root non-structural C,N,P concentration
  !     ZCKI,PCKI,ZPKI,PZKI=N,P inhibition effect on N,P uptake
  !     dtPerPlantRootH2OUptake=water uptake at time step for gas flux calculations
  !
  IF(RootNonstructElmConc_pvr(ielmc,N,L,NZ).GT.ZERO)THEN
    FZUP=AMIN1(RootNonstructElmConc_pvr(ielmc,N,L,NZ)/(RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmn,N,L,NZ)/ZCKI) &
      ,safe_adb(RootNonstructElmConc_pvr(ielmp,N,L,NZ),(RootNonstructElmConc_pvr(ielmp,N,L,NZ)+RootNonstructElmConc_pvr(ielmn,N,L,NZ)/ZPKI)))
    FPUP=AMIN1(RootNonstructElmConc_pvr(ielmc,N,L,NZ)/(RootNonstructElmConc_pvr(ielmc,N,L,NZ)+RootNonstructElmConc_pvr(ielmp,N,L,NZ)/PCKI) &
      ,safe_adb(RootNonstructElmConc_pvr(ielmn,N,L,NZ),(RootNonstructElmConc_pvr(ielmn,N,L,NZ)+RootNonstructElmConc_pvr(ielmp,N,L,NZ)/PZKI)))
  ELSE
    FZUP=0.0_r8
    FPUP=0.0_r8
  ENDIF
  Nutruptk_fNlim_rpvr(N,L,NZ)=FZUP
  Nutruptk_fPlim_rpvr(N,L,NZ)=FPUP
  !NN=0
  if(PlantPopuLive_pft(NZ).GT.ZERO)then
    PerPlantRootH2OUptake  = AZMAX1(-RPlantRootH2OUptk_pvr(N,L,NZ))/PlantPopuLive_pft(NZ)
  else
    PerPlantRootH2OUptake  = 0._r8
  endif
  dtPerPlantRootH2OUptake = PerPlantRootH2OUptake*dts_gas
  !
  !     FACTORS CONSTRAINING O2 AND NUTRIENT UPTAKE AMONG
  !     COMPETING ROOT,MYCORRHIZAL AND MICROBIAL POPULATIONS
  !     IN BAND AND NON-BAND SOIL ZONES FROM DEMAND CALCULATED
  !     IN PREVIOUS HOUR
  !
  !
  IF(RO2EcoDmndPrev_vr(L).GT.ZEROS)THEN
    FOXYX=AMAX1(FracMinRoot4Uptake_rpvr(N,L,NZ),RootO2Dmnd4Resp_pvr(N,L,NZ)/RO2EcoDmndPrev_vr(L))
  ELSE
    FOXYX=FracPRoot4Uptake(N,L,NZ)
  ENDIF

  call PrintInfo('end '//subname)
  end associate
  end subroutine GetUptakeCapcity

!----------------------------------------------------------------------------------------------------
  subroutine RootExudates(I,J,N,L,NZ)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: N,L
  integer, intent(in) :: NZ
  character(len=*), parameter :: subname='RootExudates'

  real(r8) :: CPOOLX,CPOOLT
  real(r8) :: PPOOLX,ZPOOLX
  real(r8) :: VLWatMicPK,VLWatMicPT
  real(r8) :: XFRE(NumPlantChemElms)
  real(r8) :: RootExudE,scal
  integer :: K,NE
  !     begin_execution
  associate(                                                     &
    RootMycoNonstElms_rpvr  => plt_biom%RootMycoNonstElms_rpvr  ,& !input  :root layer nonstructural element, [g d-2]
    ZERO4Groth_pft          => plt_biom%ZERO4Groth_pft          ,& !input  :threshold zero for plang growth calculation, [-]
    ZEROS                   => plt_site%ZEROS                   ,& !input  :threshold zero for numerical stability,[-]
    ZEROS2                  => plt_site%ZEROS2                  ,& !input  :threshold zero for numerical stability,[-]
    VLWatMicPM_vr           => plt_site%VLWatMicPM_vr           ,& !input  :soil micropore water content, [m3 d-2]
    RootVH2O_pvr            => plt_morph%RootVH2O_pvr           ,& !input  :root layer volume water, [m2 d-2]
    FracBulkSOMC_vr         => plt_soilchem%FracBulkSOMC_vr     ,& !input  :fraction of total organic C in complex, [-]
    DOM_MicP_vr             => plt_soilchem%DOM_MicP_vr         ,& !input  :dissolved organic C micropore, [gC d-2]
    Soil2RootMycoExudE_pvr  => plt_rbgc%Soil2RootMycoExudE_pvr   & !inoput :soil to root/myco flux via exudation, (+) to plants, (-) to soil [g d-2 h-1]
  )
  call PrintInfo('beg '//subname)
  !
  !     ROOT EXUDATION OF C, N AND P DEPENDS ON CONCN DIFFERENCES
  !     BETWEEN ROOT NON-STRUCTURAL POOLS AND SOIL DISSOLVED POOLS
  !
  !     VLWatMicPMM=soil micropore water volume
  !     FracBulkSOMC_vr=fraction of total SOC in each substrate K from nitro.f
  !     RootVH2O_pvr=root aqueous volume
  !     CPOOLR,ZPOOLR,PPOOLR=non-structural C,N,P in root,myco
  !     XFRC,XFRN,XFRP=nonstructural C,N,P exchg at root-soil DOC equilibrium
  !     OQC=soil DOC
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     FEXUC,FEXUN,FEXUP=rate constant for root C,N,P exudation
  !     vapor pressure in watsub.f
  !      CanopyEvapTransLHeat_pft,SFLXC=canopylatent,sensible heat fluxes
  !      RA=canopy boundary layer resistance
  !     OSTR=O2 stress indicator
  !
  D195: DO K=1,jcplx
    VLWatMicPK=VLWatMicPM_vr(NPH,L)*FracBulkSOMC_vr(K,L)
    IF(VLWatMicPK.GT.ZEROS2.AND.RootVH2O_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
      VLWatMicPT=VLWatMicPK+RootVH2O_pvr(N,L,NZ)
      CPOOLX      = AMIN1(1.25E+03_r8*RootVH2O_pvr(N,L,NZ),RootMycoNonstElms_rpvr(ielmc,N,L,NZ))
      XFRE(ielmc) = (DOM_MicP_vr(idom_doc,K,L)*RootVH2O_pvr(N,L,NZ)-CPOOLX*VLWatMicPK)/VLWatMicPT

      !XFRC >0, plant absorbs dom 
      Soil2RootMycoExudE_pvr(ielmc,N,K,L,NZ)=FEXUDE(ielmc)*XFRE(ielmc)
      IF(DOM_MicP_vr(idom_doc,K,L).GT.ZEROS .AND. RootMycoNonstElms_rpvr(ielmc,N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
        CPOOLT                                = DOM_MicP_vr(idom_doc,K,L)+RootMycoNonstElms_rpvr(ielmc,N,L,NZ)
        ZPOOLX                                = 0.1_r8*RootMycoNonstElms_rpvr(ielmn,N,L,NZ)
        PPOOLX                                = 0.1_r8*RootMycoNonstElms_rpvr(ielmp,N,L,NZ)
        XFRE(ielmn)                           = (DOM_MicP_vr(idom_don,K,L)*RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-ZPOOLX*DOM_MicP_vr(idom_doc,K,L))/CPOOLT
        XFRE(ielmp)                           = (DOM_MicP_vr(idom_dop,K,L)*RootMycoNonstElms_rpvr(ielmc,N,L,NZ)-PPOOLX*DOM_MicP_vr(idom_doc,K,L))/CPOOLT
        Soil2RootMycoExudE_pvr(ielmn,N,K,L,NZ) = FEXUDE(ielmn)*XFRE(ielmn)    !>0 into roots
        Soil2RootMycoExudE_pvr(ielmp,N,K,L,NZ) = FEXUDE(ielmp)*XFRE(ielmp)
      ELSE
        Soil2RootMycoExudE_pvr(ielmn,N,K,L,NZ)=0.0_r8
        Soil2RootMycoExudE_pvr(ielmp,N,K,L,NZ)=0.0_r8
      ENDIF
    ELSE
      Soil2RootMycoExudE_pvr(1:NumPlantChemElms,N,K,L,NZ)=0.0_r8
    ENDIF
  ENDDO D195
  
  DO NE=1,NumPlantChemElms
    RootExudE=0._r8
    DO K=1,jcplx
      RootExudE=RootExudE+Soil2RootMycoExudE_pvr(NE,N,K,L,NZ)
    ENDDO
    !Avoid excessive exudation 
    if(RootExudE<0._r8)then
      if(RootMycoNonstElms_rpvr(NE,N,L,NZ)+RootExudE<0._r8)then
        scal=RootMycoNonstElms_rpvr(NE,N,L,NZ)/(-RootExudE)*0.999999_r8
        DO K=1,jcplx
          Soil2RootMycoExudE_pvr(NE,N,K,L,NZ)=Soil2RootMycoExudE_pvr(NE,N,K,L,NZ)*scal
        ENDDO
      endif
    endif    
  ENDDO      
  call PrintInfo('end '//subname)
  end associate
  end subroutine RootExudates

!----------------------------------------------------------------------------------------------------
  subroutine SumNutrientUptake(I,J,NZ)

  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), parameter :: subname='SumNutrientUptake'
  real(r8) :: DOM_uptk(1:NumPlantChemElms)
  real(r8) :: scal
  integer :: L  
  integer :: K,NE,N
  !     begin_execution
  associate(                                                    &
    NU                     => plt_site%NU                      ,& !input  :current soil surface layer number, [-]
    ZERO                   => plt_site%ZERO                    ,& !input  :threshold zero for numerical stability, [-]
    ZEROS2                 => plt_site%ZEROS2                  ,& !input  :threshold zero for numerical stability,[-]
    THETW_vr               => plt_soilchem%THETW_vr            ,& !input  :volumetric water content, [m3 m-3]
    VLSoilPoreMicP_vr      => plt_soilchem%VLSoilPoreMicP_vr   ,& !input  :volume of soil layer, [m3 d-2]
    Myco_pft               => plt_morph%Myco_pft               ,& !input  :mycorrhizal type (no or yes),[-]
    RootNutUptake_pvr      => plt_rbgc%RootNutUptake_pvr       ,& !input  :root uptake of Nutrient band, [g d-2 h-1]
    MaxSoilLays4Root_pft   => plt_morph%MaxSoilLays4Root_pft   ,& !input  :maximum soil layer number for all root axes,[-]
    REcoDOMProd_vr         => plt_bgcr%REcoDOMProd_vr          ,& !inoput :net microbial DOC flux, [gC d-2 h-1]
    RootMycoNonstElms_rpvr => plt_biom%RootMycoNonstElms_rpvr  ,& !inoput :root layer nonstructural element, [g d-2]
    Soil2RootMycoExudE_pvr => plt_rbgc%Soil2RootMycoExudE_pvr  ,& !inoput :root uptake (+ve) - exudation (-ve) of DOE, [g d-2 h-1]
    Soil2RootMycoExudE_pft => plt_rbgc%Soil2RootMycoExudE_pft  ,& !inoput :total root uptake (+ve) - exudation (-ve) of dissolved element, [g d-2 h-1]
    DOM_MicP_vr            => plt_soilchem%DOM_MicP_vr         ,& !inoput :dissolved organic matter in micropore, [g d-2]
    DOM_MicP_drib_vr       => plt_soilchem%DOM_MicP_drib_vr     & !inoput :dribbling flux for micropore dom,[g d-2]
  )
  call PrintInfo('beg '//subname)
  !
  !     TOTAL C,N,P EXCHANGE BETWEEN ROOTS AND SOIL
  !
  !     RDFOMC,RDFOMN,RDFOMP=nonstructl C,N,P exchange:-ve=exudn,+ve=uptake
  !     UPOMC,UPOMN,UPOMP=net PFT root-soil nonstructl C,N,P exchange
  !     RootNH4Uptake_pvr,RootNH4BUptake_pvr,RUPN03,RootNO3BUptake_pvr=uptake from non-band,band of NH4,NO3
  !     RootH2PO4Uptake_pvr,RootNutUptake_pvr,RootNutUptake_pvr,RootH1PO4BUptake_pvr=uptake from non-band,band of H2PO4,HPO4
  !     RootNH4Uptake_pft,RootNO3Uptake_pft,RootH2PO4Uptake_pft,RootHPO4Uptake_pft=PFT uptake of NH4,NO3,H2PO4,HPO4
  !

  D9501: DO L=MaxSoilLays4Root_pft(NZ),NU,-1
    IF(VLSoilPoreMicP_vr(L).GT.ZEROS2 .AND. THETW_vr(L).GT.ZERO)THEN  
      D295: DO K=1,jcplx
        DOM_uptk=0._r8
        DO N=1,Myco_pft(NZ)
          !positve means uptake from soil
          DO NE=1,NumPlantChemElms
            DOM_uptk(NE)=DOM_uptk(NE)+Soil2RootMycoExudE_pvr(NE,N,K,L,NZ)
          ENDDO
        ENDDO
        if(any(DOM_uptk/=0._r8))then
          DO NE=1,NumPlantChemElms
            call SubstrateDribbling(DOM_uptk(NE),DOM_MicP_drib_vr(NE,K,L),DOM_MicP_vr(NE,K,L))            
          ENDDO

          DO NE=1,NumPlantChemElms
            REcoDOMProd_vr(NE,K,L)=REcoDOMProd_vr(NE,K,L)-DOM_uptk(NE)
          ENDDO
        endif

        ! Opposing root/mycorrhizal exchanges can cancel in the soil budget.
        ! Apply each population's transfer even when the net soil flux is zero.
        DO N=1,Myco_pft(NZ)
          DO NE=1,NumPlantChemElms
            Soil2RootMycoExudE_pft(NE,NZ)     = Soil2RootMycoExudE_pft(NE,NZ)+Soil2RootMycoExudE_pvr(NE,N,K,L,NZ)
            RootMycoNonstElms_rpvr(NE,N,L,NZ) = RootMycoNonstElms_rpvr(NE,N,L,NZ)+Soil2RootMycoExudE_pvr(NE,N,K,L,NZ)
          ENDDO
        ENDDO
      ENDDO D295
      
    endif
  ENDDO  D9501
  !write(1114,*)I*1000+J/24.,1,2,0.0,plt_biom%RootMycoNonstElms_rpvr(ielmc,1,2,NZ),'af'//subname   
  call PrintInfo('end '//subname)
  end associate
  end subroutine SumNutrientUptake
  ![tail]
end module NutUptakeMod
