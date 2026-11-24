module PlantAPI
  ! Public entry points. Keep calls inside their original PFT loops to preserve
  ! accumulation order; transfer implementations compile independently.
  use PlantColumnTransferMod, only : ReceivePlantColumns, &
    SendPlantColumnInputs, &
    SendPlantColumnState
  use PlantStateTransferMod, only : ReceivePlantState, &
    SendPlantState
  use PlantCanopyTransferMod, only : ReceivePlantCanopy, &
    SendPlantCanopy
  use PlantRootTransferMod, only : ReceivePlantRoots, &
    SendPlantRoots
  use PlantRootAxisTransferMod, only : ReceivePlantRootAxes, &
    SendPlantRootAxes
  use PlantTraitsTransferMod, only : SendPlantTraits
  use data_kind_mod, only : yearIJ_type
  use DebugToolMod, only : PrintInfo
  use GridDataType, only : NUI_col,NK_col
  use PlantMgmtDataType, only : NP0_col
  use EcoSIMCtrlDataType, only : DazCurrYear
  use PlantDataRateType, only : RootCO2Autor_col,RootCO2Autor_vr,RootCO2Ar2Root_col, &
    RootCO2Ar2Root_vr,RootCO2Ar2Soil_col,RootCO2Ar2Soil_vr,RootO2_TotSink_col,RootO2_TotSink_vr
  use EcoSIMSolverPar, only : NPH
  use SoilBGCDataType, only : REcoUptkSoilO2M_vr
  use SoilWaterDataType, only : FracAirFilledSoilPoreM_vr
  use PlantRootBGCAPIData, only : plt_rbgc
  use PlantSoilChemistryAPIData, only : plt_soilchem
implicit none

  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: PlantAPISend
  public :: PlantAPIRecv
  public :: RefreshPlantColumnForBGC

  contains

!------------------------------------------------------------------------------------------

  subroutine PlantAPIRecv(I,J,NY,NX)
  !
  !DESCRIPTION
  !
  implicit none
  integer, intent(in) :: I,J,NY,NX
  character(len=*), parameter :: subname='PlantAPIRecv'
  integer :: NZ,L,I1

  call PrintInfo('beg '//subname)
  I1=I+1;if(I1>DazCurrYear)I1=1
  call ReceivePlantColumns(I1,NY,NX)
  DO NZ=1,NP0_col(NY,NX)
    call ReceivePlantState(I,NY,NX,NZ)
    call ReceivePlantCanopy(NY,NX,NZ)
    call ReceivePlantRoots(NY,NX,NZ)
    call ReceivePlantRootAxes(NY,NX,NZ)
  ENDDO
  DO  L=NUI_col(NY,NX),NK_col(NY,NX)
    RootCO2Autor_col(NY,NX)   = RootCO2Autor_col(NY,NX)+RootCO2Autor_vr(L,NY,NX)
    RootCO2Ar2Root_col(NY,NX) = RootCO2Ar2Root_col(NY,NX)+ RootCO2Ar2Root_vr(L,NY,NX)
    RootCO2Ar2Soil_col(NY,NX) = RootCO2Ar2Soil_col(NY,NX)+RootCO2Ar2Soil_vr(L,NY,NX)
    RootO2_TotSink_col(NY,NX)    = RootO2_TotSink_col(NY,NX) + RootO2_TotSink_vr(L,NY,NX)
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine PlantAPIRecv


!------------------------------------------------------------------------------------------

  subroutine PlantAPISend(yearIJ,NY,NX)
  !
  !DESCRIPTION
  !Send data to plant model
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NY,NX
  integer :: L,M,NZ,I1,I
  character(len=*), parameter :: subname='PlantAPISend'

  call PrintInfo('beg '//subname)
  I=yearIJ%I
  call SendPlantColumnInputs(I,NY,NX,I1)
  DO NZ=1,NP0_col(NY,NX)
    call SendPlantTraits(NY,NX,NZ)
  ENDDO

  call SendPlantColumnState(I1,NY,NX)
  NZ100: DO NZ=1,NP0_col(NY,NX)
    call SendPlantState(I,NY,NX,NZ)
    call SendPlantCanopy(NY,NX,NZ)
    call SendPlantRoots(NY,NX,NZ)
    call SendPlantRootAxes(NY,NX,NZ)
  ENDDO NZ100

  DO L=1,NK_col(NY,NX)
    DO M=1,NPH
      plt_rbgc%REcoUptkSoilO2M_vr(M,L)           = REcoUptkSoilO2M_vr(M,L,NY,NX)
      plt_soilchem%FracAirFilledSoilPoreM_vr(M,L) = FracAirFilledSoilPoreM_vr(M,L,NY,NX)
    ENDDO
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine PlantAPISend

!------------------------------------------------------------------------------------------
  subroutine RefreshPlantColumnForBGC(yearIJ,NY,NX)
    ! Description:
    ! Refresh column soil, microbial and surface ledgers after transport and
    ! microbial BGC, retaining the parked plant phenology and hydraulic state.
    type(yearIJ_type), intent(in) :: yearIJ ! current model year/day/hour [-]
    integer, intent(in) :: NY,NX ! column row and column indices [-]
    character(len=*), parameter :: subname='RefreshPlantColumnForBGC'
    integer :: I1,L,M
    call PrintInfo('beg '//subname)
    associate( &
      RootO2Uptake => plt_rbgc%REcoUptkSoilO2M_vr, & ! output: microbial soil O2 uptake [gO2/cell/substep], substep/layer
      AirFraction => plt_soilchem%FracAirFilledSoilPoreM_vr) ! output: soil air-filled pore fraction [-], substep/layer
      call SendPlantColumnInputs(yearIJ%I,NY,NX,I1)
      call SendPlantColumnState(I1,NY,NX)
      do L=1,NK_col(NY,NX)
        do M=1,NPH
          RootO2Uptake(M,L)=REcoUptkSoilO2M_vr(M,L,NY,NX)
          AirFraction(M,L)=FracAirFilledSoilPoreM_vr(M,L,NY,NX)
        enddo
      enddo
    end associate
    call PrintInfo('end '//subname)
  end subroutine RefreshPlantColumnForBGC

end module PlantAPI
