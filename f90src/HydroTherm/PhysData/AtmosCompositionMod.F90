module AtmosCompositionMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use ClimForcDataType
  use TracerIDMod
  use EcosimConst, only: Tref
  use DebugToolMod, only: PrintInfo
  implicit none
  private
  public :: UpdateAtmosGasMassConcentrations
contains
!------------------------------------------------------------------------------------------
  subroutine UpdateAtmosGasMassConcentrations(NHW,NHE,NVN,NVS,PressureRatio)
    ! Convert current mixing ratios [umol/mol] to elemental mass densities
    ! [g/m3] using mean air temperature [K] and atmospheric pressure [kPa].
    ! This conversion does not update dissolved boundary species or accumulators.
    integer, intent(in) :: NHW,NHE,NVN,NVS
    real(r8), optional, intent(in) :: PressureRatio !override P/101.325; ATS startup uses one
    integer :: NX,NY
    real(r8) :: tPBOT,tmp
    character(len=*), parameter :: subname='UpdateAtmosGasMassConcentrations'

    call PrintInfo('beg '//subname)
    DO NX=NHW,NHE
      DO NY=NVN,NVS
        if(present(PressureRatio))then
          tPBOT=PressureRatio
        else
          tPBOT=PBOT_col(NY,NX)/1.01325E+02_r8
        endif

        tmp                          = Tref/TairKClimMean_col(NY,NX)*tPBOT
        CCO2EI_gperm3_col(NY,NX)                = CO2EI_col(NY,NX)*5.36E-04_r8*tmp

        AtmGasCgperm3_col(idg_CO2,NY,NX) = CO2E_col(NY,NX)*5.36E-04_r8*tmp  !gC/m3
        AtmGasCgperm3_col(idg_CH4,NY,NX) = CH4E_col(NY,NX)*5.36E-04_r8*tmp  !gC/m3
        AtmGasCgperm3_col(idg_O2,NY,NX)  = OXYE_col(NY,NX)*1.43E-03_r8*tmp  !gO/m3
        AtmGasCgperm3_col(idg_N2,NY,NX)  = Z2GE_col(NY,NX)*1.25E-03_r8*tmp  !gN/m3
        AtmGasCgperm3_col(idg_N2O,NY,NX) = Z2OE_col(NY,NX)*1.25E-03_r8*tmp  !gN/m3
        AtmGasCgperm3_col(idg_AR,NY,NX)  = ARGE_col(NY,NX)*1.78E-03_r8*tmp  !gAr/m3
        AtmGasCgperm3_col(idg_NH3,NY,NX) = ZNH3E_col(NY,NX)*6.25E-04_r8*tmp !gN/m3
        AtmGasCgperm3_col(idg_H2,NY,NX)  = H2GE_col(NY,NX)*8.92E-05_r8*tmp  !gH/m3
        AtmGasCgperm3_col(idg_NH3B,NY,NX) = AtmGasCgperm3_col(idg_NH3,NY,NX)

      ENDDO
    ENDDO
    call PrintInfo('end '//subname)
  end subroutine UpdateAtmosGasMassConcentrations
!------------------------------------------------------------------------------------------
end module AtmosCompositionMod
