module MediumRootCohortMod
  ! Young tissue is a subset of the existing medium-root totals, not a new pool.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use DebugToolMod, only: PrintInfo
  implicit none
  private
  public :: MediumCohortGeometry, GrowMediumCohorts
contains
!------------------------------------------------------------------------------------------
  subroutine MediumCohortGeometry(number,length,elements,young_number,young_length,young_elements,volume,n,l,radius)
  !
  !! Description:
  ! All stocks describe one soil layer and parent structural-root group.
  ! Output cohort order: 1 = young, 2 = established (total minus young).
  real(r8),intent(in) :: number            !total population medium-axis count, including young axes, [axes per grid cell]
  real(r8),intent(in) :: length            !total population medium-root length, including young tissue, [m per grid cell]
  real(r8),intent(in) :: elements(:)       !total structural C/N/P, including young tissue; order C,N,P, [g element per grid cell]
  real(r8),intent(in) :: young_number      !young subset of the total population axis count, [axes per grid cell]
  real(r8),intent(in) :: young_length      !young subset of the total population root length, [m per grid cell]
  real(r8),intent(in) :: young_elements(:) !young subset of structural C/N/P; order C,N,P, [g element per grid cell]
  real(r8),intent(in) :: volume            !tissue volume per unit structural carbon used to infer radius, [m3 gC-1]
  real(r8),intent(out) :: n(2)             !population axis count in each output cohort, [axes per grid cell]
  real(r8),intent(out) :: l(2)             !population total length in each output cohort, [m per grid cell]
  real(r8),intent(out) :: radius(2)        !radius inferred from each cohort C stock and total length, [m]
  character(len=*),parameter :: subname='MediumCohortGeometry'
  real(r8) :: c(2)
  integer :: k

  call PrintInfo('beg '//subname)
  n(1)=MIN(MAX(young_number,0._r8),MAX(number,0._r8));n(2)=MAX(number-n(1),0._r8)
  l(1)=MIN(MAX(young_length,0._r8),MAX(length,0._r8));l(2)=MAX(length-l(1),0._r8)
  c(1)=MIN(MAX(young_elements(1),0._r8),MAX(elements(1),0._r8));c(2)=MAX(elements(1)-c(1),0._r8)
  radius=0._r8
  do k=1,2
    if(l(k)>0._r8.and.c(k)>0._r8)radius(k)=SQRT(c(k)*volume/(ACOS(-1._r8)*l(k)))
  enddo
  call PrintInfo('end '//subname)
  end subroutine MediumCohortGeometry
!------------------------------------------------------------------------------------------
  subroutine GrowMediumCohorts(growth,target,specific_length,thick_fraction,initial_length,mature_length, &
    extension_rate,number,length,elements,young_number,young_length,young_elements)
  !
  ! Description:
  ! The caller has already paid growth C/N/P and respiration. Partition that
  ! investment between established tips, young tips, and newly established axes.
  ! All stocks describe one soil layer and parent structural-root group.
  ! Young stocks are subsets of the totals; promotion changes only the subset.
  real(r8),intent(in) :: growth(:)            !paid net structural growth for this hourly call; order C,N,P, [g element per grid cell]
  real(r8),intent(in) :: target               !activated target for total population medium-axis count, [axes per grid cell]
  real(r8),intent(in) :: specific_length      !elongation length per unit structural C assigned to elongation, [m gC-1]
  real(r8),intent(in) :: thick_fraction       !fraction of growth C reserved for radial investment, [0-1]
  real(r8),intent(in) :: initial_length       !paid establishment length required per newly recruited axis, [m axis-1]
  real(r8),intent(in) :: mature_length        !length per axis at promotion to established; must exceed initial_length, [m axis-1]
  real(r8),intent(in) :: extension_rate       !potential extension rate per existing axis; assumes a one-hour call, [m axis-1 h-1]
  real(r8),intent(inout) :: number            !total population medium-axis count, including young axes, [axes per grid cell]
  real(r8),intent(inout) :: length            !total population medium-root length, including young tissue, [m per grid cell]
  real(r8),intent(inout) :: elements(:)       !total structural C/N/P, including young tissue; order C,N,P, [g element per grid cell]
  real(r8),intent(inout) :: young_number      !young subset of the total population axis count, [axes per grid cell]
  real(r8),intent(inout) :: young_length      !young subset of the total population root length, [m per grid cell]
  real(r8),intent(inout) :: young_elements(:) !young subset of structural C/N/P; order C,N,P, [g element per grid cell]

  character(len=*),parameter :: subname='GrowMediumCohorts'
  real(r8) :: deficit,demand_new,demand_old,budget,scale,new_axes,new_length,extension
  real(r8) :: young_extension,radial_c,young_c,young_fraction,old_young_fraction,promote,move_length,move_fraction

  call PrintInfo('beg '//subname)
  if(growth(1)>0._r8.and.specific_length>0._r8)then
    deficit    = MAX(target-number,0._r8)
    demand_new = deficit*initial_length

    ! A one-hour timestep; the configured extension rate is per existing tip.
    demand_old = MAX(number,0._r8)*extension_rate
    budget     = growth(1)*MAX(1._r8-thick_fraction,0._r8)*specific_length
    scale      = 0._r8

    if(demand_new+demand_old>0._r8)scale=MIN(1._r8,budget/(demand_new+demand_old))

    new_length     = scale*demand_new;new_axes = new_length/initial_length
    extension      = scale*demand_old
    young_fraction = 0._r8

    if(number>0._r8)young_fraction=MIN(1._r8,MAX(young_number/number,0._r8))

    young_extension=extension*young_fraction

    ! Unused developmental elongation capacity is explicitly radial investment.
    radial_c           = MAX(growth(1)-(new_length+extension)/specific_length,0._r8)
    old_young_fraction = 0._r8

    if(elements(1)>0._r8)old_young_fraction=MIN(1._r8,MAX(young_elements(1)/elements(1),0._r8))
    if(number<=0._r8)old_young_fraction=1._r8

    young_c        = (new_length+young_extension)/specific_length+radial_c*old_young_fraction
    young_fraction = MIN(1._r8,MAX(young_c/growth(1),0._r8))
    young_elements = young_elements+growth*young_fraction
    elements       = elements+growth
    young_length   = young_length+new_length+young_extension
    length         = length+new_length+extension
    young_number   = young_number+new_axes;number = number+new_axes

    ! Two-stage moment closure: extension beyond initial length supplies axes
    ! promoted at mature_length. Transfer length and its C/N/P together, without
    ! changing any aggregate stock. The residual young cohort retains its paid
    ! initial length; no old root is shortened by recruitment.
    if(young_number>0._r8.and.young_length>0._r8)then
      promote     = MIN(young_number,MAX(young_length-young_number*initial_length,0._r8)/(mature_length-initial_length))
      move_length = MIN(young_length,promote*mature_length)
      !
      if(promote>=young_number)move_length=young_length
      !
      move_fraction  = move_length/young_length
      young_elements = young_elements*(1._r8-move_fraction)
      young_length   = young_length-move_length
      young_number   = MAX(young_number-promote,0._r8)
    endif
  endif
  call PrintInfo('end '//subname)
  end subroutine GrowMediumCohorts
!------------------------------------------------------------------------------------------
end module MediumRootCohortMod
