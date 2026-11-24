module SoilDiagsMod
  use data_kind_mod,      only: r8 => DAT_KIND_R8
  use CanopyDataType,     only: QVegET_col
  use GridDataType,       only: NU_col, NL_col
  use EcoSimConst,        only: RGasC
  use abortutils,         only: endrun
  use TracerPropMod,      only: MolecularWeight
  use ChemTranspDataType, only: GasSolbility_vr
  use EcoSIMCtrlDataType, only: ZEROS
  use DebugToolMod
  use SoilBGCDataType
  use GridDataType
  use SurfLitterDataType
  use CanopyDataType
  use SnowDataType
  use BalanceCheckDataType
  use SoilWaterDataType
  use ClimForcDataType
  use SurfSoilDataType
  use EcosimBGCFluxType
  use SoilHeatDataType
  use PlantDataRateType
  use EcoSimSumDataType
implicit none
  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: DiagSoilGasPressure
  contains

  subroutine DiagSoilGasPressure(I,J,NHW,NHE,NVN,NVS)  

  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS

  character(len=*), parameter :: subname='DiagSoilGasPressure'
  integer :: idg, NY,NX, L
  real(r8) :: GasMassSolubility(idg_beg:idg_end)    !g m3 /(mol d-2)
  real(r8) :: GasPres(idg_beg:idg_end)              !gas pressure for each tracer [Pa]

  call PrintInfo('beg '//subname)
  DO  NX=NHW,NHE
    DO  NY=NVN,NVS
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        if(VLWatMicP_vr(L,NY,NX).GT.ZEROS(NY,NX) .and. trcs_solml_vr(idg_CO2,L,NY,NX).GT.ZEROS(NY,NX))THEN
          do idg=idg_beg,idg_NH3
            GasMassSolubility(idg) =MolecularWeight(idg)*GasSolbility_vr(idg,L,NY,NX)*VLWatMicP_vr(L,NY,NX)  !conver into carbon g C/mol
          enddo  
          
          GasMassSolubility(idg_NH3B)=GasMassSolubility(idg_NH3)

          Soil_Gas_pressure_vr(L,NY,NX)=0._r8
          do idg=idg_beg,idg_end
            GasPres(idg)                  = trcs_solml_vr(idg,L,NY,NX)*RGasC*TKS_vr(L,NY,NX)/GasMassSolubility(idg)
            Soil_Gas_pressure_vr(L,NY,NX) = Soil_Gas_pressure_vr(L,NY,NX)+GasPres(idg)
          enddo
          do idg=idg_beg,idg_NH3
            Soil_Gas_Frac_vr(idg,L,NY,NX) = GasPres(idg)/Soil_Gas_pressure_vr(L,NY,NX)*1.E6_r8
          ENDDO
          Soil_Gas_Frac_vr(idg_NH3,L,NY,NX)=Soil_Gas_Frac_vr(idg_NH3,L,NY,NX)+GasPres(idg_NH3B)/Soil_Gas_pressure_vr(L,NY,NX)*1.E6_r8
        ELSE
          Soil_Gas_pressure_vr(L,NY,NX)             = 0._r8
          Soil_Gas_Frac_vr(idg_beg:idg_NH3,L,NY,NX) = 0._r8
        ENDIF
      ENDDO
    ENDDO
  ENDDO    
  call PrintInfo('end '//subname)
  end subroutine DiagSoilGasPressure
end module SoilDiagsMod
