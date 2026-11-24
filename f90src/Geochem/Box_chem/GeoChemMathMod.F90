module GeoChemMathMod
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SoluteParMod, only: TSLX
implicit none

  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME = &
  __FILE__

  contains


  function A2BC_ChemReact_Quad(A,B,C,DK)result(RAProd)
  !
  !using quadratic equation to compute reaction rates for chemical reaction type
  ! A <-> B+C 
  ! with DK=[B][C]/[A]
  ! where A, B, and C are chemical activity
  !
  ! (A-x) DK = (B+x)(C+x)
  ! x^2 + (B+C+DK)x +(B*C-A*DK)=0
  ! S0=(B+C+DK)
  ! Delta=S0*S0 - 4(B*C-A*DK)
  ! x=(-S0+sqrt(Delta))*tslx
  ! x> 0, consumption of A
  ! x'=(S0-sqrt(delta))*tslx, x'>0, production of A
  implicit none
  real(r8), intent(in) :: A   
  real(r8), intent(in) :: B
  real(r8), intent(in) :: C
  real(r8), intent(in) :: DK

  real(r8) :: RAProd  !>0, A increases
  real(r8) :: S0
  real(r8) :: S1

  S0     = B+C+DK
  ! Keep small positive discriminants: concentration cutoffs do not apply to squares.
  ! Use a sum of nonnegative terms and 4*(B*C-K*A)/(S0+sqrt(D)) to avoid cancellation.
  S1     = MAX(0._r8,(B-C)**2+2._r8*DK*(B+C)+DK**2+4._r8*DK*A)
  if(S0+SQRT(S1)>0._r8)then
    RAProd = AMAX1(AMIN1((4._r8*((B*C-A*DK)/(S0+SQRT(S1))))*tslx,B,C),-A)
  else
    RAProd = 0._r8
  endif

  end function A2BC_ChemReact_Quad

!------------------------------------------------------------------------------------------

  function A2BC_ChemReact_Grad(A,B,C,DK)result(RAProd)
  !using gradient equation to compute reaction rates for chemical reaction type
  ! A <-> B+C 
  ! with DK=[B][C]/[A]
  ! where A, B, and C are chemical activity

  implicit none
  real(r8), intent(in) :: A   
  real(r8), intent(in) :: B
  real(r8), intent(in) :: C
  real(r8), intent(in) :: DK

  real(r8) :: RAProd   !>0, A increases

  RAprod=tslx*(B*C-DK*A)/(B+DK)

  end function A2BC_ChemReact_Grad
end module GeoChemMathMod
