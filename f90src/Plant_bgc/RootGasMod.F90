module RootGasMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use minimathmod,   only: safe_adb, vapsat, AZMAX1, AZMIN1
  use minimathmod,   only: fixEXConsumpFlux,isclose, symmetric_flux_limiter, SubstrateDribbling
  use EcoSIMCtrlMod, only : lcoarseroot
  use DebugToolMod
  use EcosimConst
  use EcoSIMSolverPar
  use UptakePars
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantMorphologyAPIData, only : plt_morph
  use PlantPhenologyAPIData, only : plt_pheno
  use PlantSoilChemistryAPIData, only : plt_soilchem
  use PlantBiomassAPIData, only : plt_biom
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc
  use TracerIDMod
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public :: RootSoilGasExchange
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine RootSoilGasExchange(I,J,N,L,NZ,FineRootRadius,FracPRoot4Uptake,FracSoiLayByPrimRoot_pvr,&
    RootEffLen4Absorption_pvr,dtPerPlantRootH2OUptake,FOXYX,trc_gasml_loc,trc_solml_loc,PopPlantO2Uptake_vr)
  !
  !In the gas exchange between soil and roots, only O2 and CO2 involve active biochemical production/consumption
  !the other gases undergo physical exchange according to different water volumes inside/outside the roots.
  !O2 uptake includes that from inside the root through diffusion, and that directly from the soil through transpiration+diffusion.
  !CO2 is first released into roots, then diffuse into soil and atmosphere 
  !
  !The O2 transport for mycorrhizae needs some revision. Currently, they only get oxygen from soil. However, in saturated soil,
  !they get oxygen mostly from the roots that are linked parenchyma. Therefore, the current model tends to underestimate mycorrhizae
  !biomass for saturated soils.
  implicit none
  integer , intent(in) :: I,J,N,L,NZ
  real(r8), intent(in) :: FineRootRadius(jroots,JZ1),FracPRoot4Uptake(jroots,JZ1,JP1)
  real(r8), intent(in) :: FracSoiLayByPrimRoot_pvr(JZ1,JP1)
  real(r8), intent(in) :: RootEffLen4Absorption_pvr(jroots,JZ1)
  real(r8), intent(in) :: dtPerPlantRootH2OUptake   !root water uptake 
  real(r8), intent(in) :: FOXYX
  real(r8), intent(inout) :: trc_solml_loc(idg_beg:idg_end)    !local copy of aqueous phase of the volatile tracers
  real(r8), intent(inout) :: trc_gasml_loc(idg_beg:idg_NH3)    !local copy of gaesous phase of the volatile tracers  
  real(r8), intent(out):: PopPlantO2Uptake_vr

  character(len=*), parameter :: subname='RootSoilGasExchange'

  !local variables
  logical :: chkActiveUptk
  integer :: M,MX,NR
  real(r8) :: B,C
  real(r8) :: trcg_rootml_loc(idg_beg:idg_NH3)     !vol gaseous mass volatile in roots
  real(r8) :: trcs_rootml_loc(idg_beg:idg_NH3)     !vol aqueous mass volatile in roots
  real(r8) :: trcaqu_conc_soi_loc(idg_beg:idg_end)   !aqueous volatile tracer concentration in soil
  real(r8) :: trc_conc_root_loc(idg_beg:idg_NH3)   !aqueous tracer concentration in roots  [g m-3]
  real(r8) :: GasDifc_tscaled(idg_beg:idg_NH3)
  real(r8) :: SolDifc_tscaled(idg_beg:idg_NH3)
  real(r8) :: trcg_gcon_loc(idg_beg:idg_NH3)       !gaseous volatile concentration in roots 
  real(r8) :: COXYR                                !aqueous O2 concentration at root surface used for uptake [gO/m3]
  real(r8) :: trcs_maxRootml_loc(idg_beg:idg_NH3)
  real(r8) :: DIFOP                                !aqueous diffusivity of O2 within root
  real(r8) :: DifAqueVolatile(idg_beg:idg_end)     !aqueous diffusivity from soil to root
  real(r8) :: DIFOX,DiffusivitySolutEffP
  real(r8) :: DFAGas(idg_beg:idg_NH3)                  !maximum volatile flux from atmosphere to roots
  real(r8) :: DFGP,ROxySoil2Uptk

  real(r8) :: O2AquaDiffusvityP
  real(r8) :: RTVLWA  !root H2O vol for NH4
  real(r8) :: RTVLWB  !root H2O vol for NH4B
  real(r8) :: RGasTranspFlxPrev(idg_beg:idg_NH3)       !diagnosed increment gas flux
  real(r8) :: ROXYLX                         !prescribed dissolved O2 removal per substep (negative for supply)
  real(r8) :: O2AqueousForcing,O2GasForcing   !accepted forcing for the current substep
  real(r8) :: SoilCO2Debt                     !unfunded aqueous CO2 withdrawal in the local predictor, [gC d-2]
  real(r8) :: O2PhaseExchangeRate            !soil gas-water relaxation coefficient
  real(r8) :: O2SoilUptakeRequested          !soil uptake before the donor-budget limit
  real(r8) :: RGas_DisolvSoil_flx(idg_beg:idg_end)     !gas dissolution into aqueous concentration in soil
  real(r8) :: RTCR1,RTCR2,fMedium,SurvivingFineRootLength
  real(r8) :: RTCRA                                    !root conductance scalar for gas transport between atmosphere and root inner space [m]
  real(r8) :: RTARRX
  real(r8) :: RootCO2Prod_tscaled                      !root CO2 gas efflux due to root respiration at time step for gas flux calculations
  real(r8) :: RRADS
  real(r8) :: RSolUptkTransp(idg_beg:idg_end)          !volatile tracer uptake aided by transpiration
  real(r8) :: RootOxyUptakePerPlant
  real(r8) :: ROxySoil2UptkPerPlant                    !aqueous O2 uptake flux from soil due to diffusion and transpiration-aided-advection
  real(r8) :: ROxyRoot2UptkPerPlant                    !aqueous O2 diffusion uptake flux from inside roots
  real(r8) :: ROxyRoot2Uptk                            !oyxgen uptake by root, [gO/d2]
  real(r8) :: RSoilSolute2Roots(idg_beg:idg_end)       !root uptake of aqueous volatile from soil to inside roots
  real(r8) :: RDXAqueous(idg_beg:idg_end)              !maximum fluxes can be taken by roots
  real(r8) :: RDFAqueous(idg_beg:idg_end)
  real(r8) :: RUPNTX  !total uptake of NH4 from soil into roots
  real(r8) :: Root_gas2sol_flx(idg_beg:idg_NH3)    !gas dissolution into aqueous phase of the volatile tracers in roots
  real(r8) :: trcg_air2root_flx_loc(idg_beg:idg_NH3)   !diffusion flux of gas from atmosphere to inside roots
  real(r8) :: THETW1,THETM
  real(r8) :: RootOxyDemandPerPlant                    !O2 demand per plant at time step of flux calculation
  real(r8) :: DisolvedGasVolume(idg_beg:idg_NH3),VLWatMicPMO,VLWatMicPMM,VLsoiAirPMM
  real(r8) :: VOLWSP,VLWatMicPMA,VLWatMicPMB,VOLWSA,VOLWSB
  real(r8) :: VOLWAqueous(idg_beg:idg_end)  !solubility scaled aqueous volume 
  real(r8) :: VOLPNH3,VOLPNH3B
  real(r8) :: trcg_rootml_beg(idg_beg:idg_NH3)
  real(r8) :: trcs_rootml_beg(idg_beg:idg_NH3)
  real(r8) :: dtrc_err(idg_beg:idg_NH3)
  real(r8) :: X,trcs0,RTCRM,condC,condM
  real(r8) :: ZH3PA,ZH3PB,ZH3GA,ZH3GB
  real(r8) :: oscal=1._r8        !assuming O2 is close to half saturation constant
  integer  :: idg
  
!     begin_execution
  associate(                                                                 &
    RootSoilCO2Debt_vr           => plt_rbgc%RootSoilCO2Debt_vr            ,& !output :residual predictor debt for paired soil mass/debt handoff, [gC d-2]
    RootStrutElms_pft            => plt_biom%RootStrutElms_pft              ,& !input  :plant root structural element mass, [g d-2]
    ZERO4Groth_pft               => plt_biom%ZERO4Groth_pft                 ,& !input  :threshold zero for plang growth calculation, [-]
    PlantPopuLive_pft            => plt_site%PlantPopuLive_pft              ,& !input  :plant population, [d-2]
    CumSoilThickMidL_vr          => plt_site%CumSoilThickMidL_vr            ,& !input  :depth to middle of soil layer from surface of grid cell, [m]
    AtmGasc                      => plt_site%AtmGasc                        ,& !input  :atmospheric gas concentrations, [g m-3]
    ZEROS                        => plt_site%ZEROS                          ,& !input  :threshold zero for numerical stability,[-]
    ZERO                         => plt_site%ZERO                           ,& !input  :threshold zero for numerical stability, [-]
    VLWatMicPM_vr                => plt_site%VLWatMicPM_vr                  ,& !input  :soil micropore water content, [m3 d-2]
    VLsoiAirPM_vr                => plt_site%VLsoiAirPM_vr                  ,& !input  :soil air content, [m3 d-2]
    TortMicPM_vr                 => plt_site%TortMicPM_vr                   ,& !input  :micropore soil tortuosity, [m3 m-3]
    FILMM_vr                     => plt_site%FILMM_vr                       ,& !input  :soil water film thickness, [m]
    RGasTranspFlxPrev_vr         => plt_bgcr%RGasTranspFlxPrev_vr           ,& !input  :net gaseous flux, [g d-2 h-1]
    RO2AquaSourcePrev_vr         => plt_bgcr%RO2AquaSourcePrev_vr           ,& !input  :net aqueous O2 flux, [g d-2 h-1]
    ZERO4Uptk_pft                => plt_rbgc%ZERO4Uptk_pft                  ,& !input  :threshold zero for uptake calculation, [-]
    RootRespPotent_pvr           => plt_rbgc%RootRespPotent_pvr             ,& !input  :root respiration unconstrained by O2, [g d-2 h-1]
    RootO2Dmnd4Resp_pvr          => plt_rbgc%RootO2Dmnd4Resp_pvr            ,& !input  :root O2 demand from respiration, [g d-2 h-1]
    RootCO2AutorX_pvr            => plt_rbgc%RootCO2AutorX_pvr              ,& !input  :root respiration from previous time step, [g d-2 h-1]
    TScal4Difsvity_vr            => plt_soilchem%TScal4Difsvity_vr          ,& !input  :temperature effect on diffusivity,[-]
    trcs_VLN_vr                  => plt_soilchem%trcs_VLN_vr                ,& !input  :effective relative tracer volume, [-]
    GasDifcT_vr                   => plt_soilchem%GasDifcT_vr               ,& !input  :temperature-scaled gaseous diffusivity, [m2 h-1]
    SoluteDifusvtyT_vr            => plt_soilchem%SoluteDifusvtyT_vr        ,& !input  :temperature-scaled aqueous diffusivity, [m2 h-1]
    GasSolbility_vr              => plt_soilchem%GasSolbility_vr            ,& !input  :gas solubility, [m3 m-3]
    SoilWatAirDry_vr             => plt_soilchem%SoilWatAirDry_vr           ,& !input  :air-dry water content, [m3 m-3]
    VLSoilMicP_vr                => plt_soilchem%VLSoilMicP_vr              ,& !input  :total micropore volume in layer, [m3 d-2]
    FracAirFilledSoilPoreM_vr    => plt_soilchem%FracAirFilledSoilPoreM_vr  ,& !input  :soil air-filled porosity, [m3 m-3]
    DiffusivitySolutEffM_vr      => plt_soilchem%DiffusivitySolutEffM_vr    ,& !input  :coefficient for dissolution - volatilization, [-]
    iPlantCalendar_brch          => plt_pheno%iPlantCalendar_brch           ,& !input  :plant growth stage, [-]
    RootPoreTortu4Gas_pft        => plt_morph%RootPoreTortu4Gas_pft         ,& !input  :power function of root porosity used to calculate root gaseous diffusivity, [-]
    Root1stRadius_pvr            => plt_morph%Root1stRadius_pvr             ,& !input  :root layer diameter primary axes, [m]
    Root2ndEffLen4uptk_rpvr      => plt_morph%Root2ndEffLen4uptk_rpvr       ,& !input  :root layer average length, [m]
    RootPoreVol_pvr              => plt_morph%RootPoreVol_pvr               ,& !input  :root layer volume air, [m2 d-2]
    RootAbsorbLenPerPlant_pvr    => plt_morph%RootAbsorbLenPerPlant_pvr     ,& !input  :root layer length per plant, [m p-1]
    Root2ndXNumL_rpvr            => plt_morph%Root2ndXNumL_rpvr             ,& !input  :root layer number axes, [d-2]
    Root2ndLen_rpvr              => plt_morph%Root2ndLen_rpvr               ,& !input  :population fine-root length per structural axis, [m]
    RootMyco2ndStrutElms_rpvr    => plt_biom%RootMyco2ndStrutElms_rpvr      ,& !input  :fine-root structural C/N/P per axis, [g]
    NumStructuralRootAxes_pft    => plt_morph%NumStructuralRootAxes_pft     ,& !input  :number of active structural-root groups, [-]
    Root2ndRadius_rpvr           => plt_morph%Root2ndRadius_rpvr            ,& !input  :root layer diameter secondary axes, [m]
    RootRaidus_rpft              => plt_morph%RootRaidus_rpft               ,& !input  :root internal radius, [m]
    RootVH2O_pvr                 => plt_morph%RootVH2O_pvr                  ,& !input  :root space volume occupied by water in each layer, [m2 d-2]   
    RootPorosity_pft             => plt_morph%RootPorosity_pft              ,& !input  :root porosity, [m3 m-3]
    Root1stXNumL_pvr             => plt_morph%Root1stXNumL_pvr              ,& !input  :root layer number primary axes, [d-2]
    NGTopRootLayer_pft           => plt_morph%NGTopRootLayer_pft            ,& !input  :soil layer at planting depth, [-]
    RootMediumXNum_pvr           => plt_morph%RootMediumXNum_pvr            ,& !inoput :Number of medium size root axes in layer, [d-2]    
    Root1stTransptArea_pvr       => plt_morph%Root1stTransptArea_pvr        ,& !input  :transport area by 1st order root, [-] 
    RootFineFrac2Med_rpvr        => plt_morph%RootFineFrac2Med_rpvr         ,& !fine-axis-count-weighted medium-root routing by category, [-]
    RootMedTransptArea_pvr       => plt_morph%RootMedTransptArea_pvr        ,& !input  :transport area by medisum size roots, [-]
    MainBranchNum_pft            => plt_morph%MainBranchNum_pft             ,& !input  :number of main branch,[-]
    RootMediumLength_pvr         => plt_morph%RootMediumLength_pvr          ,& !input  :within layer mean length for medium size root axes, [m d-2]              
    RootO2Uptk_pvr               => plt_rbgc%RootO2Uptk_pvr                 ,& !inoput :aqueous O2 flux from roots to root water, [g d-2 h-1]
    RAutoRootO2Limter_rpvr       => plt_rbgc%RAutoRootO2Limter_rpvr         ,& !inoput :O2 constraint to root respiration (0-1), [-]
    REcoUptkSoilO2M_vr           => plt_rbgc%REcoUptkSoilO2M_vr             ,& !inoput :total O2 sink, [g d-2 t-1]
    RCO2Emis2Root_rpvr           => plt_rbgc%RCO2Emis2Root_rpvr             ,& !inoput :aqueous CO2 flux from roots to root water, [g d-2 h-1]
    trcg_air2root_flx_pvr        => plt_rbgc%trcg_air2root_flx_pvr          ,& !inoput :gaseous tracer flux through roots, [g d-2 h-1]
    trcg_Root_gas2aqu_flx_vr     => plt_rbgc%trcg_Root_gas2aqu_flx_vr       ,& !inoput :dissolution (+ve) - volatilization (-ve) gas flux in roots, [g d-2 h-1]
    RootUptkSoiSol_pvr           => plt_rbgc%RootUptkSoiSol_pvr             ,& !inoput :aqueous CO2 flux from roots to soil water, [g d-2 h-1]
    trcg_rootml_pvr              => plt_rbgc%trcg_rootml_pvr                ,& !inoput :root gas content, [g d-2]
    trcs_rootml_pvr              => plt_rbgc%trcs_rootml_pvr                ,& !inoput :root aqueous content, [g d-2]
    RootCO2Ar2Soil_pvr           => plt_rbgc%RootCO2Ar2Soil_pvr             ,& !inoput :root respiration released to soil, [gC d-2 h-1]
    RootCO2Ar2RootX_rpvr         => plt_rbgc%RootCO2Ar2RootX_rpvr           ,& !inoput :root respiration released to root, [gC d-2 h-1]
    RootO2_TotSink_pvr           => plt_bgcr%RootO2_TotSink_pvr             ,& !output :root O2 sink for autotrophic respiraiton, [gC d-2 h-1]
    RootAtmGasConductance_rpvr   => plt_rbgc%RootAtmGasConductance_rpvr      & !output :Conductance for gas diffusion between roots and atmosphere [m3 d-2 h-1]
  )
  
  call PrintInfo('beg '//subname)
  SoilCO2Debt=0._r8

  chkActiveUptk=RootRespPotent_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ) .AND. RootVH2O_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ) .AND. FOXYX.GT.ZERO4Uptk_pft(NZ)
!  if(N.eq.1)write(9340,*)I*1000+J/24.,L,chkActiveUptk,RootStrutElms_pft(ielmc,NZ),FracSoiLayByPrimRoot_pvr(L,NZ)
  IF(.not.chkActiveUptk)THEN    
    !root not active, add CO2 to soil
    RootCO2Ar2Soil_pvr(L,NZ)       = RootCO2Ar2Soil_pvr(L,NZ)-RootCO2AutorX_pvr(N,L,NZ)
    PopPlantO2Uptake_vr            = 0.0_r8
    RootO2_TotSink_pvr(N,L,NZ)     = 0._r8
    RAutoRootO2Limter_rpvr(N,L,NZ) = 1._r8
    RootCO2Ar2RootX_rpvr(N,L,NZ)   = 0._r8
  else
    !
    !     INITIALIZE VARIABLES USED IN ROOT GAS EXCHANGE
    !     (CO2, O2, CH4, N2, N2O, NH3, H2, Ar)
    !
    !     RTVLWA,RTVLWB=root aqueous volume in non-band,band
    !     dts_gas=time step of flux calculation (1/(NPH*NPT))
    !     RGasTranspFlxPrev(idg_O2)=diagnosed net O2 gas flux at time step of flux calculation
    !     RGasTranspFlxPrev(idg_CO2)=diagnosed net CO2 gas flux at time step of flux calculation
    !     ROXYLX=diagnosed net O2 aqueous flux at time step of flux calculation
    !

    DO idg=idg_beg,idg_NH3
      trcg_rootml_loc(idg) = AMAX1(ZERO4Groth_pft(NZ),trcg_rootml_pvr(idg,N,L,NZ))
      trcs_rootml_loc(idg) = AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_pvr(idg,N,L,NZ))
      trcg_rootml_beg(idg) = trcg_rootml_pvr(idg,N,L,NZ)
      trcs_rootml_beg(idg) = trcs_rootml_pvr(idg,N,L,NZ)
    ENDDO


    RTVLWA = RootVH2O_pvr(N,L,NZ)*trcs_VLN_vr(ids_NH4,L)
    RTVLWB = RootVH2O_pvr(N,L,NZ)*trcs_VLN_vr(ids_NH4B,L)

    DO idg=idg_beg,idg_NH3
      if(idg.eq.idg_CO2 .or. idg.eq.idg_O2)then
        RGasTranspFlxPrev(idg)  = RGasTranspFlxPrev_vr(idg,L)*dts_gas*FOXYX    
      else
        RGasTranspFlxPrev(idg)  = 0._r8
      endif
    ENDDO
    
    ROXYLX                = -RO2AquaSourcePrev_vr(L)*FOXYX*dts_gas   !positive removal from dissolved phase
    RootOxyDemandPerPlant = RootO2Dmnd4Resp_pvr(N,L,NZ)*dts_gas/PlantPopuLive_pft(NZ)*oscal
    !
    !     GASEOUS AND AQUEOUS DIFFUSIVITIES IN ROOT AND SOIL
    !
    IF(N.EQ.ipltroot .AND. iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).EQ.0)then
      do idg=idg_beg,idg_NH3    
        GasDifc_tscaled(idg) = GasDifcT_vr(idg,L)*dts_gas
      ENDDO
    else  
      do idg=idg_beg,idg_NH3
        GasDifc_tscaled(idg) = GasDifcT_vr(idg,L)*dts_gas*RootPoreTortu4Gas_pft(N,NZ)
      ENDDO
    endif

    do idg=idg_beg,idg_NH3
      SolDifc_tscaled(idg) = SoluteDifusvtyT_vr(idg,L)*dts_gas*FOXYX         !assuming oxygen stress affects uptake of other solutes
    enddo
    
    O2AquaDiffusvityP = SoluteDifusvtyT_vr(idg_O2,L)*dts_gas

    RGas_DisolvSoil_flx(idg_beg:idg_end) = 0.0_r8
    Root_gas2sol_flx(idg_beg:idg_NH3)    = 0.0_r8
    !
    !  ROOT CONDUCTANCE TO GAS TRANSFER
    !
    IF(RootStrutElms_pft(ielmc,NZ).GT.ZERO4Groth_pft(NZ) .AND. FracSoiLayByPrimRoot_pvr(L,NZ).GT.ZERO)THEN
      !primary roots conductance scalar[m]
      RTCR1 = AMAX1(PlantPopuLive_pft(NZ),Root1stXNumL_pvr(L,NZ))*Root1stTransptArea_pvr(N,L,NZ)/CumSoilThickMidL_vr(L)
      !
      RTCRM = 0._r8
      fMedium = 0._r8
      ! Absent or disabled medium roots leave the direct primary-fine pathway.
      if(lcoarseroot .and. RootMediumLength_pvr(L,NZ).GT.0._r8 .and. RootMediumXNum_pvr(L,NZ).GT.0._r8)then
        fMedium = RootFineFrac2Med_rpvr(N,L,NZ)
        if(RootMedTransptArea_pvr(N,L,NZ).GT.0._r8)then
          RTCRM=RootMediumXNum_pvr(L,NZ)*RootMedTransptArea_pvr(N,L,NZ)/RootMediumLength_pvr(L,NZ)
        endif
      endif

      SurvivingFineRootLength = 0._r8
      DO NR=1,NumStructuralRootAxes_pft(NZ)
        if(RootMyco2ndStrutElms_rpvr(ielmc,N,L,NR,NZ).GT.ZERO4Groth_pft(NZ))then
          SurvivingFineRootLength = SurvivingFineRootLength+AZMAX1(Root2ndLen_rpvr(N,L,NR,NZ))
        endif
      ENDDO

      RTCR2 = 0._r8
      if(Root2ndXNumL_rpvr(N,L,NZ).GT.0._r8 .and. Root2ndEffLen4uptk_rpvr(N,L,NZ).GT.0._r8)then
        RTCR2 = (Root2ndXNumL_rpvr(N,L,NZ)*PICON*Root2ndRadius_rpvr(N,L,NZ)**2) &
          /(FracSoiLayByPrimRoot_pvr(L,NZ)*Root2ndEffLen4uptk_rpvr(N,L,NZ))
      endif
      if(Root2ndXNumL_rpvr(N,L,NZ).GT.0._r8 .or. SurvivingFineRootLength.GT.0._r8)then
        ! A zero count alone does not prove that fine-root tissue is absent.
        ! Keep the fine-root pathway closed if surviving tissue has no conducting axes.
        condC = RTCR2*(1._r8-fMedium)
        condM = SeriesConductance(RTCR2*fMedium,RTCRM)
        RTCRA = SeriesConductance(RTCR1,condC+condM)
      elseif(RTCRM.GT.0._r8)then
        ! Neither fine-root axes nor surviving length: retain the primary-medium pathway.
        RTCRA = SeriesConductance(RTCR1,RTCRM)
      else
        ! Primary roots themselves retain a gas pathway when lateral roots are absent.
        RTCRA = RTCR1
      endif
    ELSE
      RTCRA=0.0_r8
    ENDIF

    !
    !     VARIABLES USED TO CALCULATE ROOT GAS TRANSFER
    !     BETWEEN AQUEOUS AND GASEOUS PHASES
    !
    IF(N.EQ.ipltroot .AND. iPlantCalendar_brch(ipltcal_Emerge,MainBranchNum_pft(NZ),NZ).GT.0 &
      .AND. RootAbsorbLenPerPlant_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
      RTARRX = RootEffLen4Absorption_pvr(N,L)/RootRaidus_rpft(N,NZ)       !length to radius ratio~tortuosity
      DIFOP  = O2AquaDiffusvityP*RTARRX
      DO idg  = idg_beg, idg_NH3
        DisolvedGasVolume(idg) = RootVH2O_pvr(N,L,NZ)*GasSolbility_vr(idg,L)
        DFAGas(idg)            = GasDifc_tscaled(idg)*RTCRA
      ENDDO
    ELSE
      !it is in the seed phase
      RTARRX                             = 0.0_r8
      DIFOP                              = 0.0_r8
      DisolvedGasVolume(idg_beg:idg_NH3) = 0.0_r8
      DFAGas(idg_beg:idg_NH3)            = 0.0_r8
      DFAGas(idg_O2)                     = 0._r8
    ENDIF

    DO idg=idg_beg,idg_NH3
      RootAtmGasConductance_rpvr(idg,N,L,NZ) = AMIN1(DFAGas(idg),RootPoreVol_pvr(N,L,NZ))
    ENDDO

    DFGP                = AMIN1(1.0,XNPD*SQRT(RootPorosity_pft(N,NZ))*TScal4Difsvity_vr(L))
    !The root respiration below is from the previous time step
    RootCO2Prod_tscaled = -RootCO2AutorX_pvr(N,L,NZ)*dts_gas

    !
    !     SOLVE FOR GAS EXCHANGE IN SOIL AND ROOTS DURING ROOT UPTAKE
    !     AT SMALLER TIME STEP NPH
    !
    D99: DO M=1,NPH
      !
      !     AQUEOUS GAS DIFFUSIVITY THROUGH SOIL WATER TO ROOT
      !
      if(isclose(VLWatMicPM_vr(M,L),0._r8))then
        ! Previous-hour respiration has already produced this CO2.
        ! Route this hydrological substep's share directly to soil.
        RootCO2Ar2Soil_pvr(L,NZ) = RootCO2Ar2Soil_pvr(L,NZ)+RootCO2Prod_tscaled*NPT
        cycle
      endif
      VLWatMicPMO = VLWatMicPM_vr(M,L)*FOXYX
      VLWatMicPMM = VLWatMicPM_vr(M,L)*FracPRoot4Uptake(N,L,NZ)
      VLsoiAirPMM = VLsoiAirPM_vr(M,L)*FracPRoot4Uptake(N,L,NZ)
      VOLWSP      = RootVH2O_pvr(N,L,NZ)+VLWatMicPMM
      VLWatMicPMA = VLWatMicPMM*trcs_VLN_vr(ids_NH4,L)
      VLWatMicPMB = VLWatMicPMM*trcs_VLN_vr(ids_NH4B,L)
      VOLWSA      = RTVLWA+VLWatMicPMA
      VOLWSB      = RTVLWB+VLWatMicPMB
      THETW1      = AZMAX1(VLWatMicPM_vr(M,L)/VLSoilMicP_vr(L),1.E-2_R8)

      IF(THETW1.GT.SoilWatAirDry_vr(L) .AND. FracPRoot4Uptake(N,L,NZ).GT.ZERO4Uptk_pft(NZ))THEN
        THETM  = TortMicPM_vr(M,L)*THETW1
        RRADS  = LOG((FILMM_vr(M,L)+FineRootRadius(N,L))/FineRootRadius(N,L))
        RTARRX = RootEffLen4Absorption_pvr(N,L)/RRADS

        DO idg=idg_beg,idg_end
          VOLWAqueous(idg)=VLWatMicPMM*GasSolbility_vr(idg,L)
        ENDDO
        VOLWAqueous(idg_NH3)  = VOLWAqueous(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
        VOLWAqueous(idg_NH3B) = VOLWAqueous(idg_NH3B)*trcs_VLN_vr(ids_NH4B,L)

        do idg = idg_beg, idg_NH3
          DifAqueVolatile(idg)=(THETM*SolDifc_tscaled(idg)+GasDifc_tscaled(idg) &
            *POROQ*FracAirFilledSoilPoreM_vr(M,L)**2/(VLSoilMicP_vr(L)*GasSolbility_vr(idg,L)))*RTARRX
        enddo
        DifAqueVolatile(idg_NH3B) = DifAqueVolatile(idg_NH3)
        DifAqueVolatile(idg_NH3)  = DifAqueVolatile(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
        DifAqueVolatile(idg_NH3B) = DifAqueVolatile(idg_NH3B)*trcs_VLN_vr(ids_NH4B,L)

        VOLPNH3  = VLsoiAirPMM*trcs_VLN_vr(ids_NH4,L)
        VOLPNH3B = VLsoiAirPMM*trcs_VLN_vr(ids_NH4B,L)
        !
        !     MASS FLOW OF GAS FROM SOIL TO ROOT AT SHORTER TIME STEP NPT
        !
        !     ROXYLX=prescribed soil aqueous O2 removal (negative for supply)
        !     VLWatMicPMM=micropore water volume
        !     RootVH2O_pvr,RootPoreVol_pvr=root aqueous,gaseous volume
        !     RMF*=soil convective solute flux:COS=CO2,OXS=O2,CHS=CH4,
        !     N2S=N2O,NHS=NH3 non-band,NHB=NH3 band,HGS=H2
        !
        D90: DO MX=1,NPT
          ! Apply external O2 forcing once, using fresh copies so a donor limit
          ! in one substep does not reduce the prescribed forcing in later ones.
          O2AqueousForcing=ROXYLX
          O2GasForcing=RGasTranspFlxPrev(idg_O2)
          call fixEXConsumpFlux(trc_solml_loc(idg_O2),O2AqueousForcing)
          call fixEXConsumpFlux(trc_gasml_loc(idg_O2),O2GasForcing,-1)
          do idg=idg_beg,idg_end
            if(idg.eq.idg_O2)then
              trcaqu_conc_soi_loc(idg_O2)=AMIN1(AtmGasc(idg_O2)*GasSolbility_vr(idg_O2,L),&
                AZMAX1(trc_solml_loc(idg_O2)/VLWatMicPMO))
            elseif(idg.eq.idg_NH3)then
              ! NH3 masses belong to separate nonband and fertilizer-band pools.
              ! Use their own water volumes; conductance and transpiration already
              ! carry the corresponding compartment fractions.
              trcaqu_conc_soi_loc(idg)=0._r8
              if(VLWatMicPMA.GT.0._r8) &
                trcaqu_conc_soi_loc(idg)=AZMAX1(trc_solml_loc(idg)/VLWatMicPMA)
            elseif(idg.eq.idg_NH3B)then
              trcaqu_conc_soi_loc(idg)=0._r8
              if(VLWatMicPMB.GT.0._r8) &
                trcaqu_conc_soi_loc(idg)=AZMAX1(trc_solml_loc(idg)/VLWatMicPMB)
            else
              trcaqu_conc_soi_loc(idg)=AZMAX1(trc_solml_loc(idg)/VLWatMicPMM)
            endif            
          enddo

          IF(RootPoreVol_pvr(N,L,NZ).GT.ZERO)THEN
            DO idg=idg_beg,idg_NH3
              trcg_gcon_loc(idg)=AZMAX1(trcg_rootml_loc(idg)/RootPoreVol_pvr(N,L,NZ))
            ENDDO
          ELSE
            trcg_gcon_loc(idg_beg:idg_NH3)=0.0_r8
          ENDIF

          do idg=idg_beg,idg_NH3
            trc_conc_root_loc(idg)=AZMAX1(trcs_rootml_loc(idg)/RootVH2O_pvr(N,L,NZ))
          enddo
          trc_conc_root_loc(idg_O2)=AMIN1(AtmGasc(idg_O2)*GasSolbility_vr(idg_O2,L),trc_conc_root_loc(idg_O2))
          DIFOX                    =DifAqueVolatile(idg_O2)+DIFOP
          do idg=idg_beg,idg_end
            RSolUptkTransp(idg)=dtPerPlantRootH2OUptake*trcaqu_conc_soi_loc(idg)
          enddo
          RSolUptkTransp(idg_NH3)  = RSolUptkTransp(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
          RSolUptkTransp(idg_NH3B) = RSolUptkTransp(idg_NH3B)*trcs_VLN_vr(ids_NH4B,L)
          !
          !     SOLUTION FOR MASS FLOW + DIFFUSION OF O2 IN AQUEOUS PHASES OF
          !     SOIL AND ROOT = ACTIVE UPTAKE OF O2 BY ROOT
          !     Oxygen is treated in its special way, assuming fluxes from soil and root are supporting uptake directory
          !     
          !     DifAqueVolatile(idg_O2)=O2 
          !     dtPerPlantRootH2OUptake=water uptake
          !     trcaqu_conc_soi_loc(idg_O2),trc_conc_root_loc(idg_O2)=soil,root aqueous O2 concentration
          !     RootOxyDemandPerPlant=O2 demand per plant
          !     RootOxyUptakePerPlant=root O2 uptake per plant
          !
          X=(DifAqueVolatile(idg_O2)+dtPerPlantRootH2OUptake)*trcaqu_conc_soi_loc(idg_O2)+DIFOP*trc_conc_root_loc(idg_O2)

          IF(X.GT.ZERO .AND. trc_solml_loc(idg_O2).GT.ZERO4Groth_pft(NZ))THEN
            !root take up O2
            B                     = -RootOxyDemandPerPlant-DIFOX*OXKM-X
            C                     = X*RootOxyDemandPerPlant
            RootOxyUptakePerPlant = (-B-SQRT(B*B-4.0_r8*C))/2.0_r8
            COXYR                 = (X-RootOxyUptakePerPlant)/DIFOX
            ROxySoil2UptkPerPlant = RSolUptkTransp(idg_O2)+DifAqueVolatile(idg_O2)*(trcaqu_conc_soi_loc(idg_O2)-COXYR)
            ROxyRoot2UptkPerPlant = DIFOP*(trc_conc_root_loc(idg_O2)-COXYR)
          ELSE
            !soil O2 concentration is too low, uptake is from inside root
            X=DIFOP*trc_conc_root_loc(idg_O2)
            IF(X.GT.ZERO .AND. trcs_rootml_loc(idg_O2).GT.ZERO4Groth_pft(NZ))THEN
              B                     = -RootOxyDemandPerPlant-DIFOP*OXKM-X
              C                     = X*RootOxyDemandPerPlant
              RootOxyUptakePerPlant = (-B-SQRT(B*B-4.0_r8*C))/2.0_r8
              COXYR                 = (X-RootOxyUptakePerPlant)/DIFOP
              ROxySoil2UptkPerPlant = 0.0_r8
              ROxyRoot2UptkPerPlant = DIFOP*(trc_conc_root_loc(idg_O2)-COXYR)
            ELSE
              RootOxyUptakePerPlant = 0.0_r8
              COXYR                 = 0.0_r8
              ROxySoil2UptkPerPlant = 0.0_r8
              ROxyRoot2UptkPerPlant = 0.0_r8
            ENDIF
          ENDIF
          !
          !     MASS FLOW + DIFFUSIVE EXCHANGE OF OTHER GASES
          !     BETWEEN ROOT AND SOIL, CONSTRAINED BY COMPETITION
          !     WITH OTHER ROOT AND MICROBIAL POPULATIONS
          !
          RSoilSolute2Roots(idg_beg:idg_end)=0.0_r8          
          ROxySoil2Uptk             = ROxySoil2UptkPerPlant*PlantPopuLive_pft(NZ)
          ROxyRoot2Uptk             = ROxyRoot2UptkPerPlant*PlantPopuLive_pft(NZ)
          RSoilSolute2Roots(idg_O2) = ROxySoil2Uptk
          RDFAqueous(idg_CO2)       = RSolUptkTransp(idg_CO2)+DifAqueVolatile(idg_CO2)*(trcaqu_conc_soi_loc(idg_CO2)-trc_conc_root_loc(idg_CO2))
          RDXAqueous(idg_CO2)       = (RootVH2O_pvr(N,L,NZ)*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_CO2)) &
            -VLWatMicPMM*AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_loc(idg_CO2)))/VOLWSP           !gradient flow fromm soil to roots

          IF(RDFAqueous(idg_CO2).GT.0.0_r8)THEN
            !CO2 from soil into roots          
            RSoilSolute2Roots(idg_CO2)=AMIN1(AZMAX1(RDXAqueous(idg_CO2)),RDFAqueous(idg_CO2)*PlantPopuLive_pft(NZ))
          ELSE
            !CO2 from roots into soil  
            RSoilSolute2Roots(idg_CO2)=AMAX1(AZMIN1(RDXAqueous(idg_CO2)),RDFAqueous(idg_CO2)*PlantPopuLive_pft(NZ))
          ENDIF

          IF(N.EQ.ipltroot)THEN
            !fluxes only involve roots
            DO idg=idg_beg,idg_NH3
              if(idg/=idg_CO2 .and. idg/=idg_NH3 .and. idg/=idg_O2)then
                !soil to root aqueous fluxes, transpiration + diffusion
                RDFAqueous(idg) = RSolUptkTransp(idg)+DifAqueVolatile(idg)*(trcaqu_conc_soi_loc(idg)-trc_conc_root_loc(idg))

                RDXAqueous(idg) = (RootVH2O_pvr(N,L,NZ)*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg)) &
                  -VLWatMicPMM*AMAX1(ZERO4Groth_pft(NZ),trcs_rootml_loc(idg)))/VOLWSP

                IF(RDFAqueous(idg).GT.0.0_r8)THEN
                  !flux to roots
                  RSoilSolute2Roots(idg)=AMIN1(AZMAX1(RDXAqueous(idg)),RDFAqueous(idg)*PlantPopuLive_pft(NZ))
                ELSE
                  !flux leaves roots
                  RSoilSolute2Roots(idg)=AMAX1(AZMIN1(RDXAqueous(idg)),RDFAqueous(idg)*PlantPopuLive_pft(NZ))
                ENDIF
              ENDIF
            ENDDO
            !NH3 flux from soil to roots
            RDFAqueous(idg_NH3)=RSolUptkTransp(idg_NH3)+DifAqueVolatile(idg_NH3)*(trcaqu_conc_soi_loc(idg_NH3)-trc_conc_root_loc(idg_NH3))

            IF(VOLWSA.GT.ZERO4Groth_pft(NZ))THEN
              ZH3PA              = trcs_rootml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
              RDXAqueous(idg_NH3) = (RTVLWA*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_NH3)) &
                -VLWatMicPMA*AMAX1(ZERO4Groth_pft(NZ),ZH3PA))/VOLWSA
            ELSE
              RDXAqueous(idg_NH3)=0.0_r8
            ENDIF


            IF(RDFAqueous(idg_NH3).GT.0.0_r8)THEN
              RSoilSolute2Roots(idg_NH3)=AMIN1(AZMAX1(RDXAqueous(idg_NH3)),RDFAqueous(idg_NH3)*PlantPopuLive_pft(NZ))
            ELSE
              RSoilSolute2Roots(idg_NH3)=AMAX1(AZMIN1(RDXAqueous(idg_NH3)),RDFAqueous(idg_NH3)*PlantPopuLive_pft(NZ))
            ENDIF

            RDFAqueous(idg_NH3B)=RSolUptkTransp(idg_NH3B)+DifAqueVolatile(idg_NH3B)*(trcaqu_conc_soi_loc(idg_NH3B)-trc_conc_root_loc(idg_NH3))
            IF(VOLWSB.GT.ZERO4Groth_pft(NZ))THEN
              ZH3PB               = trcs_rootml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)
              RDXAqueous(idg_NH3B) = (RTVLWB*AMAX1(ZERO4Groth_pft(NZ),trc_solml_loc(idg_NH3B)) &
                -VLWatMicPMB*AMAX1(ZERO4Groth_pft(NZ),ZH3PB))/VOLWSB
            ELSE
              RDXAqueous(idg_NH3B)=0.0_r8
            ENDIF

            IF(RDFAqueous(idg_NH3B).GT.0.0_r8)THEN
              RSoilSolute2Roots(idg_NH3B)=AMIN1(AZMAX1(RDXAqueous(idg_NH3B)),RDFAqueous(idg_NH3B)*PlantPopuLive_pft(NZ))
            ELSE
              RSoilSolute2Roots(idg_NH3B)=AMAX1(AZMIN1(RDXAqueous(idg_NH3B)),RDFAqueous(idg_NH3B)*PlantPopuLive_pft(NZ))
            ENDIF
          ENDIF

          ! Couple soil O2 uptake to donor-limited phase exchange. The pools
          ! already include external forcing, so it must not enter the budget again.
          O2PhaseExchangeRate=0._r8
          if(FracAirFilledSoilPoreM_vr(M,L).GT.AirFillPore_Min) &
            O2PhaseExchangeRate=FracPRoot4Uptake(N,L,NZ)*DiffusivitySolutEffM_vr(M,L)
          O2SoilUptakeRequested=ROxySoil2Uptk
          call LimitSoilOxygenUptake(trc_solml_loc(idg_O2),trc_gasml_loc(idg_O2), &
            VOLWAqueous(idg_O2),VLsoiAirPMM,O2PhaseExchangeRate,ROxySoil2Uptk,RGas_DisolvSoil_flx(idg_O2))

          if(O2SoilUptakeRequested.GT.0._r8 .and. ROxySoil2Uptk.LT.O2SoilUptakeRequested)then
            ! A negative inside-root uptake transfers soil O2 into root water.
            ! Reduce that transfer with its soil supply; retain independent root supply.
            if(ROxyRoot2Uptk.LT.0._r8)ROxyRoot2Uptk=ROxyRoot2Uptk*(ROxySoil2Uptk/O2SoilUptakeRequested)
          endif
          RSoilSolute2Roots(idg_O2)=ROxySoil2Uptk
          !
          !     GAS TRANSFER THROUGH ROOTS
          !
          IF(N.EQ.ipltroot .AND. RootPoreVol_pvr(N,L,NZ).GT.ZERO4Groth_pft(NZ))THEN
            RUPNTX=RSoilSolute2Roots(idg_NH3)+RSoilSolute2Roots(idg_NH3B)
            !
            !     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN ROOTS
            !     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
            !     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENT
            !
            trcs_maxRootml_loc(idg_CO2) = AZMAX1(trcs_rootml_loc(idg_CO2)+RootCO2Prod_tscaled)
            trcs_maxRootml_loc(idg_NH3) = AZMAX1(trcs_rootml_loc(idg_NH3)+RUPNTX)
            trcs_maxRootml_loc(idg_O2)  = AZMAX1(trcs_rootml_loc(idg_O2)-ROxyRoot2Uptk)
            DO idg=idg_beg,idg_NH3
              if(idg/=idg_CO2 .and. idg/=idg_NH3 .and. idg/=idg_O2)then
                trcs_maxRootml_loc(idg)=AZMAX1(trcs_rootml_loc(idg)+RSoilSolute2Roots(idg))
              endif

              Root_gas2sol_flx(idg)=AMAX1(-trcs_maxRootml_loc(idg),DFGP*(AMAX1(ZERO4Groth_pft(NZ),trcg_rootml_loc(idg))*DisolvedGasVolume(idg) &
                -trcs_maxRootml_loc(idg)*RootPoreVol_pvr(N,L,NZ))/(DisolvedGasVolume(idg)+RootPoreVol_pvr(N,L,NZ)))

              !>0._r8 atmosphere into roots, <0._r8 roots into atmosphere, assuming specific rate 1/hr
              trcg_air2root_flx_loc(idg)=RootAtmGasConductance_rpvr(idg,N,L,NZ)*(AtmGasc(idg)-trcg_gcon_loc(idg))
            enddo
          ELSE
            Root_gas2sol_flx(idg_beg:idg_NH3)      = 0.0_r8
            trcg_air2root_flx_loc(idg_beg:idg_NH3) = 0.0_r8
          ENDIF
          !
          !     UPDATE ROOT AQUEOUS, GASEOUS GAS CONTENTS AND CONCENTRATIONS'
          !     FOR ROOT AQUEOUS-GASEOUS, GASEOUS-ATMOSPHERE EXCHANGES
          !
          !Release autotrophic respiration CO2 into roots
          trcs_rootml_loc(idg_CO2) = trcs_rootml_loc(idg_CO2)+RootCO2Prod_tscaled

          DO idg=idg_beg,idg_NH3
            dtrc_err(idg)=trcg_rootml_loc(idg)+trcs_rootml_loc(idg)
            Root_gas2sol_flx(idg)=AMAX1(AMIN1(trcg_rootml_loc(idg),Root_gas2sol_flx(idg)),-trcs_rootml_loc(idg))

            trcg_rootml_loc(idg) = trcg_rootml_loc(idg)-Root_gas2sol_flx(idg)
            trcs_rootml_loc(idg) = trcs_rootml_loc(idg)+Root_gas2sol_flx(idg)

            call fixEXConsumpFlux(trcg_rootml_loc(idg),trcg_air2root_flx_loc(idg),-1)

            if(idg.eq.idg_O2)then
              ! Root withdrawal supports respiration and any reverse flow to soil.
              ! Limit both together before crediting the soil compartment.
              call LimitRootOxygenUptake(trcs_rootml_loc(idg),ROxyRoot2Uptk,ROxySoil2Uptk)
            elseif(idg.NE.idg_NH3)then
              !non-O2 gases are added to the root inside.
              call fixEXConsumpFlux(trcs_rootml_loc(idg),RSoilSolute2Roots(idg),-1)
            endif
          ENDDO

          ! Both soil NH3 compartments exchange with the same root aqueous pool.
          ! Include all inflows before sharing the donor budget among outflows.
          call LimitRootAmmoniaExchange(trcs_rootml_loc(idg_NH3), &
            RSoilSolute2Roots(idg_NH3),RSoilSolute2Roots(idg_NH3B))
          RSoilSolute2Roots(idg_O2)=ROxySoil2Uptk

          ! Soil phase exchange and pool updates use the accepted root transfers.
          !     GAS EXCHANGE BETWEEN GASEOUS AND AQUEOUS PHASES IN SOIL
          !     DURING ROOT UPTAKE DEPENDING ON CONCENTRATION DIFFERENCES
          !     CALCULATED FROM SOLUBILITIES, AND TRANSFER COEFFICIENTS
          !     FROM 'WATSUB'
          !     RSoilSolute2Roots=root aqueous gas uptake
          !     ROXYLX=soil net O2 aqueous flux

          !gas dissolution into soil
          IF(FracAirFilledSoilPoreM_vr(M,L).GT.AirFillPore_Min)THEN
            DiffusivitySolutEffP         = FracPRoot4Uptake(N,L,NZ)*DiffusivitySolutEffM_vr(M,L)
            RGas_DisolvSoil_flx(idg_CO2) = DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(idg_CO2))*VOLWAqueous(idg_CO2) &
              -(AMAX1(ZEROS,trc_solml_loc(idg_CO2))-RSoilSolute2Roots(idg_CO2))*VLsoiAirPMM)/(VOLWAqueous(idg_CO2)+VLsoiAirPMM)


            IF(N.EQ.ipltroot)THEN
              DO idg=idg_beg,idg_NH3-1
                if(idg/=idg_CO2 .and. idg/=idg_O2)then
                  RGas_DisolvSoil_flx(idg)=DiffusivitySolutEffP*(AMAX1(ZERO4Groth_pft(NZ),trc_gasml_loc(idg))*VOLWAqueous(idg) &
                    -(AMAX1(ZEROS,trc_solml_loc(idg))-RSoilSolute2Roots(idg))*VLsoiAirPMM)/(VOLWAqueous(idg)+VLsoiAirPMM)
                endif
              ENDDO

              IF(VOLWAqueous(idg_NH3)+VOLPNH3.GT.ZERO4Groth_pft(NZ))THEN
                ZH3GA               = trc_gasml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4,L)
                RGas_DisolvSoil_flx(idg_NH3) = DiffusivitySolutEffP &
                  *(AMAX1(ZERO4Groth_pft(NZ),ZH3GA)*VOLWAqueous(idg_NH3) &
                  -(AMAX1(ZEROS,trc_solml_loc(idg_NH3))-RSoilSolute2Roots(idg_NH3))*VOLPNH3) &
                  /(VOLWAqueous(idg_NH3)+VOLPNH3)
                ! Root efflux is negative; use its magnitude to bound phase exchange.
                RGas_DisolvSoil_flx(idg_NH3) = symmetric_flux_limiter( &
                  RGas_DisolvSoil_flx(idg_NH3),RSoilSolute2Roots(idg_NH3))
              ELSE
                RGas_DisolvSoil_flx(idg_NH3)=0.0_r8
              ENDIF

              IF(VOLWAqueous(idg_NH3B)+VOLPNH3B.GT.ZERO4Groth_pft(NZ))THEN
                ZH3GB                     = trc_gasml_loc(idg_NH3)*trcs_VLN_vr(ids_NH4B,L)
                RGas_DisolvSoil_flx(idg_NH3B) = DiffusivitySolutEffP &
                  *(AMAX1(ZERO4Groth_pft(NZ),ZH3GB)*VOLWAqueous(idg_NH3B) &
                  -(AMAX1(ZEROS,trc_solml_loc(idg_NH3B))-RSoilSolute2Roots(idg_NH3B))*VOLPNH3B) &
                  /(VOLWAqueous(idg_NH3B)+VOLPNH3B)
                ! Root efflux is negative; use its magnitude to bound phase exchange.
                RGas_DisolvSoil_flx(idg_NH3B) = symmetric_flux_limiter( &
                  RGas_DisolvSoil_flx(idg_NH3B),RSoilSolute2Roots(idg_NH3B))
              ELSE
                RGas_DisolvSoil_flx(idg_NH3B)=0.0_r8
              ENDIF
            ELSE
              DO idg=idg_beg,idg_NH3
                if(idg/=idg_CO2 .and. idg/=idg_O2)RGas_DisolvSoil_flx(idg)=0.0_r8
              ENDDO
            ENDIF
          ELSE
            RGas_DisolvSoil_flx(idg_beg:idg_end)=0.0_r8
          ENDIF
          ! Recalculate O2 phase exchange if the root donor reduced release to soil.
          ! Positive soil uptake was already limited before the root budget update.
          call LimitSoilOxygenUptake(trc_solml_loc(idg_O2),trc_gasml_loc(idg_O2), &
            VOLWAqueous(idg_O2),VLsoiAirPMM,O2PhaseExchangeRate,ROxySoil2Uptk,RGas_DisolvSoil_flx(idg_O2))
          RSoilSolute2Roots(idg_O2)=ROxySoil2Uptk
          !
          !     UPDATE GASEOUS, AQUEOUS GAS CONTENTS AND CONCENTRATIONS
          !     FROM GASEOUS-AQUEOUS EXCHANGE, SOIL GAS TRANSFERS
          idg=idg_CO2
          ! Keep transport forcing fixed; bound dissolution by the net gas supply.
          ! The accepted exchange also updates the aqueous pool below.
          RGas_DisolvSoil_flx(idg)=MIN(RGas_DisolvSoil_flx(idg), trc_gasml_loc(idg)+RGasTranspFlxPrev(idg))
          
          DO idg=idg_beg,idg_NH3
            if(idg.eq.idg_O2)then
              ! O2 transport forcing was applied before the coupled budget solve.
              trc_gasml_loc(idg)=trc_gasml_loc(idg)-RGas_DisolvSoil_flx(idg)
            else
              trc_gasml_loc(idg)=trc_gasml_loc(idg)-RGas_DisolvSoil_flx(idg)+RGasTranspFlxPrev(idg)
            endif
          ENDDO

          call fixEXConsumpFlux(trc_gasml_loc(idg_NH3),RGas_DisolvSoil_flx(idg_NH3B))

          !aqueous concentrations in soil
          DO idg=idg_beg,idg_end
            if(idg.eq.idg_CO2)then
              ! Fixed transport can force outgassing beyond dissolved supply.
              ! Carry that aqueous deficit through all M/MX substeps; subsequent
              ! dissolution or root release repays it before restoring free CO2.
              ! Aqueous mass minus debt retains the signed predictor budget.
              call SubstrateDribbling(RSoilSolute2Roots(idg)-RGas_DisolvSoil_flx(idg), &
                SoilCO2Debt,trc_solml_loc(idg))
            else
              trc_solml_loc(idg)=trc_solml_loc(idg)+RGas_DisolvSoil_flx(idg)-RSoilSolute2Roots(idg)
            endif
          enddo
          !
          !     ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE'
          !
          DO idg=idg_beg,idg_end
            RootUptkSoiSol_pvr(idg,N,L,NZ)=RootUptkSoiSol_pvr(idg,N,L,NZ)+RSoilSolute2Roots(idg)
          enddo
          !
          !     ACCUMULATE ROOT-ATMOSPHERE GAS EXCHANGE TO HOURLY TIME SCALE'
          !
          DO idg=idg_beg,idg_NH3
            trcg_Root_gas2aqu_flx_vr(idg,N,L,NZ) = trcg_Root_gas2aqu_flx_vr(idg,N,L,NZ)+Root_gas2sol_flx(idg)
            trcg_air2root_flx_pvr(idg,N,L,NZ)    = trcg_air2root_flx_pvr(idg,N,L,NZ)+trcg_air2root_flx_loc(idg)
          ENDDO
          !
          ! ACCUMULATE SOIL-ROOT GAS EXCHANGE TO HOURLY TIME SCALE'
          !
          ! RootO2Uptk_pvr=root O2 uptake from root internal
          ! REcoUptkSoilO2M_vr=total O2 uptake from soil by all microbial,root popns
          ! Root CO2 emission includes actual CO2 production from respiration and flux exchange with soil
          RCO2Emis2Root_rpvr(N,L,NZ)   = RCO2Emis2Root_rpvr(N,L,NZ)+RootCO2Prod_tscaled+RSoilSolute2Roots(idg_CO2)
          RootO2Uptk_pvr(N,L,NZ)       = RootO2Uptk_pvr(N,L,NZ)+ROxyRoot2Uptk    !uptake from O2 inside roots for root respiration
          REcoUptkSoilO2M_vr(M,L)      = REcoUptkSoilO2M_vr(M,L)+ROxySoil2Uptk      !uptake from soil O2 for root respiration
          RootCO2Ar2RootX_rpvr(N,L,NZ) = RootCO2Ar2RootX_rpvr(N,L,NZ)+RootCO2Prod_tscaled
        ENDDO D90       
      ELSE
        RootCO2Ar2Soil_pvr(L,NZ)       = RootCO2Ar2Soil_pvr(L,NZ)+RootCO2Prod_tscaled*NPT
      ENDIF
    ENDDO D99

    ! Export only outstanding debt, not every transient deficit. Each root/PFT
    ! owns a separate soil share; the column handoff sums those final debts.
    RootSoilCO2Debt_vr(L)=RootSoilCO2Debt_vr(L)+SoilCO2Debt

    DO idg=idg_beg,idg_NH3
      trcg_rootml_pvr(idg,N,L,NZ) = trcg_rootml_loc(idg)
      trcs_rootml_pvr(idg,N,L,NZ) = trcs_rootml_loc(idg)

      dtrc_err(idg)=trcg_rootml_beg(idg)+trcs_rootml_beg(idg)-trcg_rootml_loc(idg)-trcs_rootml_loc(idg) &
        +trcg_air2root_flx_pvr(idg,N,L,NZ)
      if(idg.eq.idg_O2)then
        dtrc_err(idg)=dtrc_err(idg)-RootO2Uptk_pvr(N,L,NZ)
      elseif(idg.eq.idg_CO2)then
        dtrc_err(idg)=dtrc_err(idg)+RCO2Emis2Root_rpvr(N,L,NZ)
      else
        dtrc_err(idg)=dtrc_err(idg)+RootUptkSoiSol_pvr(idg,N,L,NZ)
      endif 
    ENDDO    
    dtrc_err(idg_NH3)=dtrc_err(idg_NH3)+RootUptkSoiSol_pvr(idg_NH3B,N,L,NZ)

    !check mass conservation error
    !
    ! O2 CONSTRAINTS TO ROOT RESPIRATION DEPENDS UPON RATIO'
    ! O2 is taken from inside the root and directly from the aqueous soil O2.

    PopPlantO2Uptake_vr        = RootO2Uptk_pvr(N,L,NZ)+RootUptkSoiSol_pvr(idg_O2,N,L,NZ)
    RootO2_TotSink_pvr(N,L,NZ) = PopPlantO2Uptake_vr     !include O2 uptake from soil and from inside the roots
    !to be used in next iteration
    RAutoRootO2Limter_rpvr(N,L,NZ) = AMAX1(AMIN1(1.0_r8,AZMAX1(PopPlantO2Uptake_vr/RootO2Dmnd4Resp_pvr(N,L,NZ))),oscal_test)
  ENDIF

  call PrintInfo('end '//subname)
  end associate
  end subroutine RootSoilGasExchange
!----------------------------------------------------------------------------------------------------
  pure subroutine LimitRootOxygenUptake(aqueous_mass,root_uptake,soil_uptake)
  ! Positive uptake withdraws O2 from the named donor. Negative soil uptake
  ! is release from roots to soil; the sum of both uptakes is respiration.
  ! Phase exchange has already supplied/removed O2 from root water.
  implicit none
  real(r8), intent(inout) :: aqueous_mass,root_uptake,soil_uptake
  real(r8) :: requested

  requested=root_uptake
  if(requested.GT.0._r8)then
    root_uptake=MIN(requested,MAX(0._r8,aqueous_mass))
    if(soil_uptake.LT.0._r8)soil_uptake=soil_uptake*(root_uptake/requested)
  endif
  aqueous_mass=MAX(0._r8,aqueous_mass-root_uptake)
  end subroutine LimitRootOxygenUptake
!----------------------------------------------------------------------------------------------------
  pure subroutine LimitRootAmmoniaExchange(aqueous_mass,nonband_uptake,band_uptake)
  ! Positive flux enters roots. Both soil compartments share one aqueous donor
  ! after root phase exchange, including incoming NH3 from either compartment.
  ! Scale only outgoing transfers, preserving their relative shares and making
  ! the accepted transfers independent of the order of the soil compartments.
  implicit none
  real(r8), intent(inout) :: aqueous_mass,nonband_uptake,band_uptake
  real(r8) :: available,outgoing,scale

  available=MAX(0._r8,aqueous_mass)+MAX(0._r8,nonband_uptake)+MAX(0._r8,band_uptake)
  outgoing=-MIN(0._r8,nonband_uptake)-MIN(0._r8,band_uptake)
  if(outgoing.GT.available)then
    scale=available/outgoing
    if(nonband_uptake.LT.0._r8)nonband_uptake=nonband_uptake*scale
    if(band_uptake.LT.0._r8)band_uptake=band_uptake*scale
  endif
  aqueous_mass=MAX(0._r8,available+MIN(0._r8,nonband_uptake)+MIN(0._r8,band_uptake))
  end subroutine LimitRootAmmoniaExchange
!----------------------------------------------------------------------------------------------------
  pure subroutine LimitSoilOxygenUptake(aqueous_mass,gas_mass,soluble_volume,air_volume, &
    exchange_rate,uptake,exchange)
  ! Masses, uptake and exchange are population totals for one gas substep.
  ! Positive uptake consumes aqueous O2; positive exchange dissolves gaseous O2.
  ! External forcing has already been applied to both donor pools.
  implicit none
  real(r8), intent(in) :: aqueous_mass,gas_mass,soluble_volume,air_volume,exchange_rate
  real(r8), intent(inout) :: uptake
  real(r8), intent(out) :: exchange
  real(r8) :: aqueous,gas,a,b,total_volume,accessible_gas

  aqueous      = MAX(0._r8,aqueous_mass)
  gas          = MAX(0._r8,gas_mass)
  a            = 0._r8
  b            = 0._r8
  total_volume = soluble_volume+air_volume

  if(exchange_rate.GT.0._r8 .and. total_volume.GT.0._r8)then
    a = exchange_rate*soluble_volume/total_volume
    b = exchange_rate*air_volume/total_volume
  endif

  ! E(U)=a*gas-b*(aqueous-U). Requiring aqueous+E(U)-U >= 0
  ! gives U <= aqueous+a*gas/(1-b), also bounded by total donor mass.
  ! For b>=1 only the total-mass bound is needed for U>aqueous.
  accessible_gas=gas
  if(b.LT.1._r8)accessible_gas=gas*MIN(1._r8,a/(1._r8-b))
  if(uptake.GT.0._r8)uptake=MIN(uptake,aqueous+accessible_gas)
  exchange=a*gas-b*(aqueous-uptake)

  ! Limit both directions of phase transfer to the corresponding donor pool.
  exchange=MIN(gas,MAX(uptake-aqueous,exchange))

  ! Protect against roundoff at the exhausted-aqueous boundary, using the same
  ! addition order as the caller's aqueous+exchange-uptake update.
  if(uptake.GT.0._r8)uptake=MIN(uptake,MAX(0._r8,aqueous+exchange))
  end subroutine LimitSoilOxygenUptake
!----------------------------------------------------------------------------------------------------
  pure function SeriesConductance(a,b) result(g)
  implicit none
  real(r8), intent(in) :: a,b
  real(r8) :: g,small,large

  ! A closed segment closes its series pathway; avoid products and reciprocal overflow.
  g = 0._r8
  if(a.GT.0._r8 .and. b.GT.0._r8)then
    small = MIN(a,b)
    large = MAX(a,b)
    g = small/(1._r8+small/large)
  endif
  end function SeriesConductance
  ![tail]
end module RootGasMod
