module SurfaceAirExchangeMod
  ! Description:
  ! Select the ground-exchange air endpoint and assemble a current full-path
  ! resistance when no live or standing-dead canopy is exposed above snow.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use DebugToolMod, only: PrintInfo
  use MiniMathMod, only: RichardsonNumber
  use PhysPars, only: RAM
  use LandSurfDataType, only: canopy_exposed_col,RIB_col,RawIsoTAtm2CanopySinkZ_col
  use ClimForcDataType, only: TairK_col
  implicit none
  private
  public :: CanopyFreeSurface,SurfaceExchangeResistance
contains
!------------------------------------------------------------------------------------------
  function CanopyFreeSurface(NY,NX) result(canopy_free)
    ! Description:
    ! Use the same hourly exposure diagnostic as canopy aerodynamics, including
    ! fallow columns and canopies buried beneath snow. Hold it through the hour
    ! so resistance assembly and the selected air endpoint remain consistent.
    integer, intent(in) :: NY,NX       ! row and column indices [-]
    logical :: canopy_free             ! true when ground branches exchange with reference air [-]
    character(len=*), parameter :: subname='CanopyFreeSurface'

    call PrintInfo('beg '//subname)
    canopy_free=.not.canopy_exposed_col(NY,NX)
    call PrintInfo('end '//subname)
  end function CanopyFreeSurface
!------------------------------------------------------------------------------------------
  function SurfaceExchangeResistance(NY,NX,SurfaceTemperature,LocalResistance) result(resistance)
    ! Description:
    ! Return local resistance for a canopy-air endpoint, or add the atmospheric
    ! segment once for reference air. Each canopy-free patch uses its current
    ! temperature in the existing bounded atmospheric stability correction.
    integer, intent(in) :: NY,NX                ! row and column indices [-]
    real(r8), intent(in) :: SurfaceTemperature  ! current soil, litter or snow temperature [K]
    real(r8), intent(in) :: LocalResistance     ! local aerodynamic/film segment [h m-1]
    real(r8) :: resistance                      ! complete resistance to the selected air endpoint [h m-1]
    real(r8) :: neutral,ri,atmospheric
    character(len=*), parameter :: subname='SurfaceExchangeResistance'

    call PrintInfo('beg '//subname)
    resistance=LocalResistance
    if(CanopyFreeSurface(NY,NX))then
      neutral = RawIsoTAtm2CanopySinkZ_col(NY,NX)
      ri      = RichardsonNumber(RIB_col(NY,NX),TairK_col(NY,NX),SurfaceTemperature)
      ! The RawTAtm2CanopySinkZ segment, reevaluated for this patch's current
      ! temperature rather than the lagged canopy-air temperature.
      atmospheric = MAX(RAM,0.8_r8*neutral,MIN(1.2_r8*neutral,neutral/(1._r8-10._r8*ri)))
      resistance  = resistance+atmospheric
    endif
    call PrintInfo('end '//subname)
  end function SurfaceExchangeResistance
end module SurfaceAirExchangeMod
