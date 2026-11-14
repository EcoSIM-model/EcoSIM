program PlantPhaseWorkTest
  ! Description:
  ! Check ownership of intermediate plant state across two columns and phases,
  ! including pointer-backed arrays, scalar values and current shared soil inputs.
  use data_kind_mod, only: r8 => DAT_KIND_R8,yearIJ_type
  use EcoSiMParDataMod, only: micpar,pltpar
  use InitAllocMod, only: InitPlantMorphSize
  use PlantAPIData
  use PlantPhaseWorkMod
  use TracerIDMod, only: InitTracerIDs
  use DebugToolMod, only: PrintInfo
  use abortutils, only: iulog
  implicit none
  type(yearIJ_type) :: yearIJ
  integer :: NX
  call micpar%Init()
  call InitPlantMorphSize()
  pltpar%JP1=2;pltpar%JZ1=2;pltpar%jroots=2
  pltpar%NumCanopyLayers1=2;pltpar%NumOfLeafAzimuthSectors=1;pltpar%NumOfSkyAzimuthSects1=1
  pltpar%NumLeafInclinationClasses1=1;pltpar%MaxNodesPerBranch1=2
  call InitTracerIDs(.false.)
  call InitPlantAPIData()
  call InitPlantPhaseWork(1,2,1,1)
  yearIJ%year=2001;yearIJ%I=120;yearIJ%J=13
  plt_site%CO2E=399._r8;plt_ew%TKC_pft=270._r8
  do NX=1,2
    call SwapPlantPhaseColumn(1,NX)
    plt_site%CO2E=400._r8+NX
    plt_ew%TKC_pft=280._r8+NX
    plt_biom%TotBegVegE_pft=10._r8*NX
    plt_pheno%iPlantCalendar_brch=NX
    WaterDryMatter_pft(:,1,NX)=.1_r8*NX
    call SavePlantHydraulicColumn(yearIJ,1,NX)
    call check('shared scalar workspace restored',plt_site%CO2E==399._r8)
    call check('shared pointer targets restored',all(plt_ew%TKC_pft==270._r8))
  enddo
  ! Soil is supplied by the intervening physical/microbial updates, not parked.
  plt_soilchem%VLWatMicP_vr=.3_r8
  do NX=1,2
    call OpenPlantBGCColumn(yearIJ,1,NX)
    call check('column scalar retained',plt_site%CO2E==400._r8+NX)
    call check('column hydraulic array retained',all(plt_ew%TKC_pft==280._r8+NX))
    call check('plant balance baseline retained',all(plt_biom%TotBegVegE_pft==10._r8*NX))
    call check('column phenology retained',all(plt_pheno%iPlantCalendar_brch==NX))
    call check('column chemical inputs retained',all(WaterDryMatter_pft(:,1,NX)==.1_r8*NX))
    call check('current soil inputs retained',all(plt_soilchem%VLWatMicP_vr==.3_r8))
    plt_ew%TKC_pft=300._r8
    call ClosePlantBGCColumn(1,NX)
  enddo
  call DestroyPlantPhaseWork()
  call check('original API survives column cleanup',all(plt_ew%TKC_pft==270._r8))
  ! A new simulation can allocate/release the column work again without touching
  ! the original shared targets or reusing pointers to destroyed column arrays.
  call InitPlantPhaseWork(1,2,1,1)
  call DestroyPlantPhaseWork()
  call check('original API survives repeated lifecycle',all(plt_ew%TKC_pft==270._r8))
contains
!------------------------------------------------------------------------------------------
  subroutine check(Label,Condition)
    ! Description:
    ! Reject a violated column ownership or phase handoff contract.
    character(len=*), intent(in) :: Label ! diagnostic description [-]
    logical, intent(in) :: Condition ! true when the contract holds [-]
    character(len=*), parameter :: subname='check'
    call PrintInfo('beg '//subname)
    if(.not.Condition)then
      write(iulog,*)Label
      error stop 1
    endif
    call PrintInfo('end '//subname)
  end subroutine check
end program PlantPhaseWorkTest
