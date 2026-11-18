module SurfaceVaporExchangeMod
  ! Description:
  ! Evaluate a common equilibrium surface vapor density and partition signed
  ! atmospheric demand between liquid and ice, using physical volume fractions.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use EcoSimConst, only: TFice,LtHeatIceMelt,RGASC,DENSICE,EvapLHTC,SublmHTC
  use MiniMathMod, only: vapsat
  use DebugToolMod, only: PrintInfo
  implicit none
  private
  public :: SurfaceVaporDensity,PartitionSurfaceVapor
  character(len=*),private,parameter :: mod_filename= &
  __FILE__
contains
!------------------------------------------------------------------------------------------
  function SurfaceVaporDensity(Temperature,Potential,IceVolume) result(Density)
  !
  ! Description:
  ! Return vapor density in local phase equilibrium. Connected ice below
  ! freezing sets the chemical potential through constant-fusion-heat
  ! Clausius-Clapeyron, on the existing liquid vapsat reference. Do not apply
  ! the liquid matric/osmotic depression again on that ice-equilibrium branch.
  !
  real(r8),intent(in) :: Temperature ! surface temperature [K]
  real(r8),intent(in) :: Potential ! total liquid matric/osmotic potential [MPa]
  real(r8),intent(in) :: IceVolume ! vapor-connected physical ice volume [m3/cell]
  real(r8) :: Density ! equilibrium water vapor density [Mg/m3]
  real(r8) :: EquilibriumPotential
  character(len=*),parameter :: subname='SurfaceVaporDensity'

  call PrintInfo('beg '//subname)
  EquilibriumPotential=Potential
  if(IceVolume>0._r8.and.Temperature<TFice) &
    EquilibriumPotential=LtHeatIceMelt*(Temperature-TFice)/TFice
  Density=vapsat(Temperature)*EXP(18._r8*EquilibriumPotential/(RGASC*Temperature))
  call PrintInfo('end '//subname)
  end function SurfaceVaporDensity
!------------------------------------------------------------------------------------------
  subroutine PartitionSurfaceVapor(Demand,LiquidVolume,IceVolume,Temperature,DonorFraction, &
    LiquidFlux,IceFlux,LatentDerivative,IceDemand)
  !
  ! Description:
  ! Split one area/time-integrated demand by physical phase volume. Fluxes
  ! are positive into the surface and expressed as liquid-water equivalent.
  ! Bound each outward component by its own remaining donor; do not reassign
  ! an exhausted component to the other phase. Inward vapor forms frost only
  ! when supersaturated over ice; additional dry-matrix uptake stays liquid.
  ! Return d(latent heat)/d(demand) with stores and partition held fixed.
  !
  real(r8),intent(in) :: Demand ! potential signed vapor transfer [m3 water equivalent/cell/step]
  real(r8),intent(in) :: LiquidVolume ! remaining physical liquid [m3/cell]
  real(r8),intent(in) :: IceVolume ! remaining physical ice [m3/cell]
  real(r8),intent(in) :: Temperature ! surface temperature [K]
  real(r8),intent(in) :: DonorFraction ! maximum removable fraction this step [0..1]
  real(r8),intent(out) :: LiquidFlux ! accepted liquid transfer [m3 water equivalent/cell/step]
  real(r8),intent(out) :: IceFlux ! accepted ice transfer [m3 water equivalent/cell/step]
  real(r8),intent(out) :: LatentDerivative ! fixed-store latent derivative [MJ/Mg]
  real(r8),intent(in),optional :: IceDemand ! same conductance times air-minus-ice vapor density [m3 WE/cell/step]
  real(r8) :: Water,Ice,WaterFraction,IceFraction,LimitFraction
  character(len=*),parameter :: subname='PartitionSurfaceVapor'

  call PrintInfo('beg '//subname)
  LiquidFlux       = 0._r8
  IceFlux          = 0._r8
  LatentDerivative = 0._r8
  Water            = MAX(0._r8,LiquidVolume)
  Ice              = MAX(0._r8,IceVolume)
  LimitFraction    = MAX(0._r8,MIN(1._r8,DonorFraction))
  if(Demand>=0._r8)then
    LiquidFlux       = Demand
    LatentDerivative = EvapLHTC
    if(Temperature<TFice)then
      if(present(IceDemand))then
        IceFlux=MIN(Demand,MAX(0._r8,IceDemand))
        ! Both demands have the same derivative to air vapor density.
        if(IceDemand>0._r8)LatentDerivative=SublmHTC
      elseif(Ice>0._r8)then
        IceFlux          = Demand
        LatentDerivative = SublmHTC
      endif
      LiquidFlux=Demand-IceFlux
    endif
  elseif(Water+Ice>0._r8)then
    WaterFraction = Water/(Water+Ice)
    IceFraction   = Ice/(Water+Ice)
    LiquidFlux    = MAX(Demand*WaterFraction,-Water*LimitFraction)
    IceFlux       = MAX(Demand*IceFraction,-DENSICE*Ice*LimitFraction)
    if(Demand*WaterFraction>-Water*LimitFraction) &
      LatentDerivative=LatentDerivative+EvapLHTC*WaterFraction
    if(Demand*IceFraction>-DENSICE*Ice*LimitFraction) &
      LatentDerivative=LatentDerivative+SublmHTC*IceFraction
  endif
  call PrintInfo('end '//subname)
  end subroutine PartitionSurfaceVapor
end module SurfaceVaporExchangeMod
