module SoilPARAttenuationMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use EcoSimConst, only : OMCMassFrac
  implicit none

  private
  public :: CalcKSoilPAR

contains

  subroutine CalcKSoilPAR(sandMassLayer, clayMassLayer, drySoilMassLayer, orgCMassLayer, deltaZ, &
                          kSoil, tauSoil, parTransFrac, parAvgFrac, socMassFrac, &
                          betaSOCIn, KSOCIn, kSoilMineralIn, betaClayIn, betaSandIn, &
                          kSoilMinIn, kSoilMaxIn)
    implicit none

    ! Inputs
    real(r8), intent(in) :: sandMassLayer      ! sand mass in this layer, [Mg sand m-2]
    real(r8), intent(in) :: clayMassLayer      ! clay mass in this layer, [Mg clay m-2]
    real(r8), intent(in) :: drySoilMassLayer   ! dry soil mass, [Mg dry soil m-2]
    real(r8), intent(in) :: orgCMassLayer      ! organic C mass in this layer, [gC m-2]
    real(r8), intent(in) :: deltaZ             ! layer thickness, [m]
    
    ! Required output
    real(r8), intent(out) :: kSoil             ! PAR attenuation coefficient, [m-1]

    ! Optional diagnostic outputs
    real(r8), intent(out), optional :: tauSoil       ! layer optical depth, [-]
    real(r8), intent(out), optional :: parTransFrac  ! PAR_bottom / PAR_top, [-]
    real(r8), intent(out), optional :: parAvgFrac    ! PAR_layer_average / PAR_top, [-]
    real(r8), intent(out), optional :: socMassFrac   ! organic C / dry soil mass, [M C M-1 dry soil]

    ! Optional calibration parameters
    real(r8), intent(in), optional :: betaSOCIn        ! max SOC multiplier increment, default 2.0
    real(r8), intent(in), optional :: KSOCIn           ! half-saturation SOC fraction, default 0.01
    real(r8), intent(in), optional :: kSoilMineralIn   ! mineral-soil base k, [m-1], default 2000
    real(r8), intent(in), optional :: betaClayIn       ! clay texture multiplier slope, default 0.50
    real(r8), intent(in), optional :: betaSandIn       ! sand texture multiplier slope, default 0.0
    real(r8), intent(in), optional :: kSoilMinIn       ! lower bound for kSoil, [m-1], default 800
    real(r8), intent(in), optional :: kSoilMaxIn       ! upper bound for kSoil, [m-1], default 10000

    real(r8), parameter :: eps = 1.0e-12_r8
    real(r8), parameter :: tauMaxExp = 80.0_r8
    real(r8) :: sandMass
    real(r8) :: clayMass
    real(r8) :: orgCMassMg
    real(r8) :: orgMatterMassLayer
    real(r8) :: mineralMassLayer
    real(r8) :: sandFrac
    real(r8) :: clayFrac
    real(r8) :: sumSandClay
    real(r8) :: fSOC
    real(r8) :: fTexture
    real(r8) :: tauSoilCalc
    real(r8) :: parTransFracCalc
    real(r8) :: parAvgFracCalc
    real(r8) :: socMassFracCalc
    real(r8) :: betaSOC
    real(r8) :: KSOC
    real(r8) :: kSoilMineral
    real(r8) :: betaClay
    real(r8) :: betaSand
    real(r8) :: kSoilMin
    real(r8) :: kSoilMax
    real(r8) :: dz

    betaSOC      = 2.0_r8
    KSOC         = 0.01_r8
    kSoilMineral = 2000.0_r8
    betaClay     = 0.50_r8
    betaSand     = 0.00_r8
    kSoilMin     = 800.0_r8
    kSoilMax     = 10000.0_r8

    if (present(betaSOCIn))      betaSOC      = betaSOCIn
    if (present(KSOCIn))         KSOC         = KSOCIn
    if (present(kSoilMineralIn)) kSoilMineral = kSoilMineralIn
    if (present(betaClayIn))     betaClay     = betaClayIn
    if (present(betaSandIn))     betaSand     = betaSandIn
    if (present(kSoilMinIn))     kSoilMin     = kSoilMinIn
    if (present(kSoilMaxIn))     kSoilMax     = kSoilMaxIn

    sandMass = max(sandMassLayer, 0.0_r8)
    clayMass = max(clayMassLayer, 0.0_r8)
    orgCMassMg = max(orgCMassLayer, 0.0_r8)*1.e-6_r8
    orgMatterMassLayer = orgCMassMg/max(OMCMassFrac, eps)
    mineralMassLayer = max(drySoilMassLayer - orgMatterMassLayer, sandMass + clayMass, eps)

    sandFrac = min(max(sandMass/mineralMassLayer, 0.0_r8), 1.0_r8)
    clayFrac = min(max(clayMass/mineralMassLayer, 0.0_r8), 1.0_r8)
    sumSandClay = sandFrac + clayFrac
    if (sumSandClay > 1.0_r8) then
      sandFrac = sandFrac / sumSandClay
      clayFrac = clayFrac / sumSandClay
    endif

    dz = max(deltaZ, 0.0_r8)
    
    if (drySoilMassLayer > eps) then
      socMassFracCalc = orgCMassMg / drySoilMassLayer
    else
      socMassFracCalc = 0.0_r8
    endif

    ! Keep pathological inputs from forcing an optical coefficient far outside the prior range.
    socMassFracCalc = min(max(socMassFracCalc, 0.0_r8), 0.50_r8)

    ! SOC darkens soil. KSOC is the SOC mass fraction giving half the max increment.
    fSOC = 1.0_r8 + max(betaSOC, 0.0_r8) * socMassFracCalc / (socMassFracCalc + max(KSOC, eps))

    ! First-order texture prior: clay-rich soils tend to attenuate more.
    ! Sand has no default sign because measurements are site/mineralogy dependent.
    fTexture = 1.0_r8 + betaClay * clayFrac - betaSand * sandFrac
    fTexture = min(max(fTexture, 0.60_r8), 1.80_r8)

    if (kSoilMax < kSoilMin) then
      kSoilMax = kSoilMin
    endif

    kSoil = max(kSoilMineral, 0.0_r8) * fSOC * fTexture
    kSoil = min(max(kSoil, kSoilMin), kSoilMax)

    tauSoilCalc = kSoil * dz
    parTransFracCalc = exp(-min(tauSoilCalc, tauMaxExp))

    if (tauSoilCalc > eps) then
      parAvgFracCalc = (1.0_r8 - parTransFracCalc) / tauSoilCalc
    else
      parAvgFracCalc = 1.0_r8
    endif

    if (present(tauSoil))      tauSoil      = tauSoilCalc
    if (present(parTransFrac)) parTransFrac = parTransFracCalc
    if (present(parAvgFrac))   parAvgFrac   = parAvgFracCalc
    if (present(socMassFrac))  socMassFrac  = socMassFracCalc

  end subroutine CalcKSoilPAR

end module SoilPARAttenuationMod
