module EcoSimSumDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use ElmIDMod
  use TracerIDMod
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8) :: TSoilH2G_lnd                             !total soil N2,[g d-2]
  real(r8), target, allocatable :: SurfGas_lnd(:)      !total soil gas, [g d-2]
  real(r8), target, allocatable :: PlantElemntStoreLandscape(:)    !total plant element (C,N,P, etc) balance,[g d-2]  
  real(r8) :: H2GOU                                    !cumulative H2 loss through lateral and lower boundaries,[g d-2]
  real(r8) :: TION                                     !total soil ion content	,[mol d-2]
  real(r8) :: TIONIN                                   !total surface ion flux	[mol d-2]
  real(r8) :: TIONOU                                   !total subsurface ion flux,	[mol d-2]
  real(r8) :: TSEDSO                                   !total soil sediment	,[Mg d-2]
  real(r8) :: TSedmErossLoss_lnds                      !total sediment subsurface flux,	[Mg d-2]
  real(r8) :: WatMassStore_lnd                         !total soil water content	,[m3 d-2]
  real(r8) :: HeatStore_lnd                            !total soil heat content	,[MJ d-2]
  real(r8) :: TSoilO2G_lnd                             !total soil O2 content	,[g d-2]
  real(r8) :: LitRMStoreLndscap(1:NumPlantChemElms)    !total soil litter OM content	,[g d-2]
  real(r8) :: POMHumStoreLndscap(NumPlantChemElms)     !total soil POM + humus C content	,[g d-2]
  real(r8) :: TDisolNH4_lnd                            !total soil NH4 content,	[g d-2]
  real(r8) :: tNO3_lnd                                 !total soil NO3 content,	[g d-2]
  real(r8) :: TDisolPi_lnd                             !total soil PO4 content,[g d-2]
  real(r8) :: CRAIN_lnd                                !total precipitation,[m3 d-2]
  real(r8) :: HEATIN_lnd                               !total surface heat flux,	[MJ d-2]
  real(r8) :: tAmendOrgC_lnd                           !total organic C amendment,	[g d-2]
  real(r8) :: TORGN                                    !total organic N amendment,	[g d-2]
  real(r8) :: TORGP                                    !total organic P amendment,	[g d-2]
  real(r8) :: QH2OLoss_lnds                            !total subsurface water flux,	[m3 d-2]
  real(r8) :: CEVAP                                    !total evaporation	,[m3 d-2]
  real(r8) :: CRUN                                     !total surface runoff	,[m3 d-2]
  real(r8) :: HeatOut_lnds                             !total subsurface heat flux,[MJ d-2]
  real(r8) :: OXYGOU                                   !total subsurface O2 flux	,[g d-2]
  real(r8) :: TOMOU_lnds(NumPlantChemElms)             !total subsurface C flux	,[g d-2]
  real(r8) :: TZIN                                     !total surface N flux	,[g d-2]
  real(r8) :: TPIN                                     !total surface P flux	,[g d-2]
  real(r8) :: Litrfall_lnds(NumPlantChemElms)          !total LitrFall C	,[g d-2]
  real(r8) :: TGasC_lnd                                !total soil CO2	,[g d-2]
  real(r8) :: TGasN_lnd                                !total soil N2	,[g d-2]

  private :: InitAllocate
  contains

!----------------------------------------------------------------------
  subroutine InitEcoSimSum
  implicit none

  call InitAllocate

  end subroutine InitEcoSimSum
!----------------------------------------------------------------------
  subroutine InitAllocate
  implicit none
  allocate(SurfGas_lnd(idg_beg:idg_NH3));SurfGas_lnd = 0._r8
  allocate(PlantElemntStoreLandscape(NumPlantChemElms)); PlantElemntStoreLandscape=0._r8
  end subroutine InitAllocate

!----------------------------------------------------------------------
  subroutine DestructEcoSimSum
  use abortutils, only : destroy
  implicit none
  call destroy(PlantElemntStoreLandscape)
  call destroy(SurfGas_lnd)
  end subroutine DestructEcoSimSum
end module EcoSimSumDataType
