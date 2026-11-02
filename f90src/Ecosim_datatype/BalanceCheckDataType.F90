module BalanceCheckDataType
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL     
  use GridConsts

  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8), target, allocatable :: WaterErr_col(:,:)     !water error
  real(r8), allocatable :: TillageLitr2SoilH2O_col(:,:) !internal water-equivalent transfer [m3/grid cell/step]
  real(r8), target, allocatable :: HeatErr_col(:,:)      !heat error
  real(r8), target, allocatable :: EcoElmErr_col(:,:)    !ecosystem element error
  real(r8), target, allocatable :: PlantElmErr_col(:,:)  !Plant chemical element error

  contains
!----------------------------------------------------------------------

  subroutine InitBalanceCheckData
  implicit none

  allocate(WaterErr_col(JY,JX)); WaterErr_col = 0._r8
  allocate(TillageLitr2SoilH2O_col(JY,JX)); TillageLitr2SoilH2O_col = 0._r8
  allocate(HeatErr_col(JY,JX)); HeatErr_col=0._r8
  allocate(EcoElmErr_col(JY,JX)); EcoElmErr_col=0._r8
  allocate(PlantElmErr_col(JY,JX)); PlantElmErr_col=0._r8
  end subroutine InitBalanceCheckData

!----------------------------------------------------------------------
  subroutine DestructBalanceCheckData
  use abortutils, only : destroy
  implicit none

  call destroy(WaterErr_col)
  call destroy(TillageLitr2SoilH2O_col)
  call destroy(HeatErr_col)
  call destroy(EcoElmErr_col)
  call destroy(PlantElmErr_col)
  end subroutine DestructBalanceCheckData
end module BalanceCheckDataType
