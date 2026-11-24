module NumericalAuxMod
!
! DESCRIPTION:
! Auxillary variable for implementing the numerical dribbling 
!
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use ElmIDMod
  use TracerIDMod
  use EcoSIMConfig, only : jcplx => jcplxc,jsken=>jskenc
  implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__

  real(r8),target,allocatable ::  trcs_solml_drib_vr(:,:,:,:)                    !dribbling flux for soil solute transport, [g d-2 h-1]
  real(r8),target,allocatable ::  DOM_MicP_drib_vr(:,:,:,:,:)                    !dribbling flux for dissolved organic matter in micropore,	[g d-2]
  contains


  subroutine InitNumericAux()

  implicit none

  allocate(trcs_solml_drib_vr(ids_beg:ids_end,0:JZ,JY,JX));    trcs_solml_drib_vr=0._r8
  allocate(DOM_MicP_drib_vr(idom_beg:idom_end,1:jcplx,0:JZ,JY,JX));DOM_MicP_drib_vr=0._r8
  
  end subroutine InitNumericAux



!----------------------------------------------------------------------
  subroutine DestructNumericAux
  use abortutils, only : destroy
  implicit none


  call destroy(trcs_solml_drib_vr)
  call destroy(DOM_MicP_drib_vr)

  end subroutine DestructNumericAux


end module NumericalAuxMod  