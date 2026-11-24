module InitEcoSIM
!
!! DESCRIPTION
! initialize the data structure for EcoSIM
  use EcoSiMParDataMod, only: micpar, pltpar
  use abortutils,       only: endrun
  use EcoSIMCtrlMod
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: InitModules
  contains

  subroutine InitModules()

  use PlantDisturbsMod, only: InitPlantDisturbance
  use UptakesMod,       only: InitUptake
  use WatsubMod,        only: initWatsub
  use SoilBGCNLayMod,   only: InitNitro
  use RedistMod,        only: InitRedist
  use ErosionMod,       only: InitErosion
  use Hour1Mod,         only: InitHour1
  use HistDataType,     only: hist_ecosim
  use InitAllocMod,     only: InitAlloc
  use UnitMod,          only: units
  use TranspNoSaltMod,  only: InitTranspNoSalt
  use PlantInfoMod,     only: ReadPlantTraitTable
  use MicBGCAPI,        only: MicAPI_Init
  use GridConsts

  implicit  none

  call InitAlloc()

  !load plant trait table
  if(plant_model .and. .not. ats_cpl_mode) call ReadPlantTraitTable()

  call units%Initailize()

  call InitPlantDisturbance

  call InitUptake

  call initNitro

  call InitRedist

  call InitErosion

  call InitHour1(micpar%NumOfLitrCmplxs)

  call InitTranspNoSalt

  call hist_ecosim%Init(bounds)

  call MicAPI_Init

  end subroutine InitModules


end module InitEcoSIM
