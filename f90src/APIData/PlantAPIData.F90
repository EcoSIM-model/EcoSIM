module PlantAPIData
  ! Compatibility facade; consumers should import their domain modules directly.
  use PlantAPICommonData
  use PlantSiteAPIData
  use PlantPhotosynthesisAPIData
  use PlantRadiationAPIData
  use PlantMorphologyAPIData
  use PlantPhenologyAPIData
  use PlantSoilChemistryAPIData
  use PlantAllometryAPIData
  use PlantBiomassAPIData
  use PlantEnergyWaterAPIData
  use PlantDisturbanceAPIData
  use PlantBGCRatesAPIData
  use PlantRootBGCAPIData
  implicit none
  public
contains

  subroutine InitPlantAPIData()

  implicit none

  JZ1                      => pltpar%JZ1
  NumCanopyLayers1         => pltpar%NumCanopyLayers1
  JP1                      => pltpar%JP1
  NumOfLeafAzimuthSectors1 => pltpar%NumOfLeafAzimuthSectors
  NumOfSkyAzimuthSects1    => pltpar%NumOfSkyAzimuthSects1
  NumLeafInclinationClasses1    => pltpar%NumLeafInclinationClasses1
  MaxNodesPerBranch1       => pltpar%MaxNodesPerBranch1

  !the following variable should be consistent with the soil bgc model
  jcplx                => pltpar%jcplx
  jsken                => pltpar%jsken
  NumLitterGroups      => pltpar%NumLitterGroups
  MaxNumBranches       => pltpar%MaxNumBranches
  MaxNumRootAxes       => pltpar%MaxNumRootAxes
  NumOfPlantMorphUnits => pltpar%NumOfPlantMorphUnits
  NumOfPlantLitrCmplxs => pltpar%NumOfPlantLitrCmplxs
  NumGrowthStages      => pltpar%NumGrowthStages
  jroots               => pltpar%jroots

  call plt_site%Init()

  call plt_rbgc%Init()

  call plt_bgcr%Init()

  call plt_pheno%Init()

  call plt_ew%Init()

  call plt_distb%Init()

  call plt_allom%Init()

  call plt_biom%Init()

  call plt_soilchem%Init()

  call plt_rad%Init()

  call plt_photo%Init()

  call plt_morph%Init()

  call InitAllocate()
  end subroutine InitPlantAPIData

  subroutine InitAllocate()
  implicit none


  end subroutine InitAllocate

  subroutine DestructPlantAPIData
  implicit none

  call plt_pheno%Destroy()

  call plt_bgcr%Destroy()

  call plt_ew%Destroy()

  call plt_distb%Destroy()

  call plt_allom%Destroy()

  call plt_biom%Destroy()

  call plt_soilchem%Destroy()

  call plt_rad%Destroy()

  call plt_photo%Destroy()

  call plt_morph%Destroy()

  call plt_site%Destroy()



  end subroutine DestructPlantAPIData
end module PlantAPIData
