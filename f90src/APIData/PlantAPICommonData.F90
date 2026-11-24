module PlantAPICommonData
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use data_const_mod, only : spval => DAT_CONST_SPVAL
  use ElmIDMod
  use abortutils, only : destroy
  use EcoSiMParDataMod, only : pltpar
  use TracerIDMod
implicit none
  save
  public

  integer, pointer :: NumGrowthStages              !number of growth stages
  integer, pointer :: MaxNumRootAxes               !maximum number of root layers
  integer, pointer :: MaxNumBranches               !maximum number of plant branches
  integer, pointer :: JP1                          !number of plants
  integer, pointer :: NumOfSkyAzimuthSects1        !number of sectors for the sky azimuth  [0,2*pi]
  integer, pointer :: jcplx                        !number of organo-microbial complexes
  integer, pointer :: NumOfLeafAzimuthSectors1     !number of sectors for the leaf azimuth, [0,pi]
  integer, pointer :: NumCanopyLayers1           !number of canopy layers
  integer, pointer :: JZ1                          !number of soil layers
  integer, pointer :: NumLeafInclinationClasses1      !number of sectors for the leaf zenith [0,pi/2]
  integer, pointer :: MaxNodesPerBranch1           !number of canopy nodes
  integer, pointer :: jsken                        !number of kinetic components in litter, PROTEIN(*,1),CH2O(*,2),CELLULOSE(*,3),LIGNIN(*,4) IN SOIL LITTER
  integer, pointer :: NumLitterGroups              !number of litter groups nonstructural(0,*),foliar(1,*),non-foliar(jroots,*),stalk(3,*),root(4,*), coarse woody (5,*)
  integer, pointer :: NumOfPlantMorphUnits         !number of organs involved in partition
  integer, pointer :: NumOfPlantLitrCmplxs         !number of plant litter complexes
  integer, pointer :: jroots                       !number of root types, root, myco

end module PlantAPICommonData
