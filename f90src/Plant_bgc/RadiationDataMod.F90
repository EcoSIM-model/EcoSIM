module RadiationDataMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use PlantAPICommonData,  only: JP1, NumCanopyLayers1
implicit none
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8),allocatable ::  LeafSWabsorptivityL_pft(:,:)           !Leaf shortwave radiation absorptivity profile
  real(r8),allocatable ::  LeafPARabsorptivityL_pft(:,:)          !Leaf PAR absorptivity profile
  real(r8),allocatable ::  RadSWLeafAlbedoL_pft(:,:)              !Leaf shortwave radiation scattering albedo profile
  real(r8),allocatable ::  RadPARLeafAlbedoL_pft(:,:)             !leaf PAR albedo profile
  real(r8),allocatable ::  RadSWLeafTransmitanceL_pft(:,:)        !leaf shortwave radiation transmittance profile
  real(r8),allocatable ::  RadPARLeafTransmitanceL_pft(:,:)       !leaf PAR transmittance profile
  real(r8),allocatable ::  StemSWabsorptivityL_pft(:,:)           !Stem shortwave radiation absorptivity profile
  real(r8),allocatable ::  StemPARabsorptivityL_pft(:,:)          !Stem PAR absorptivity profile
  real(r8),allocatable ::  RadSWStemAlbedoL_pft(:,:)              !Stem shortwave radiation scattering albedo profile
  real(r8),allocatable ::  RadPARStemAlbedoL_pft(:,:)             !Stem PAR albedo profile
  real(r8),allocatable ::  RadSWStemTransmitanceL_pft(:,:)        !Stem shortwave radiation transmittance profile
  real(r8),allocatable ::  RadPARStemTransmitanceL_pft(:,:)       !Stem PAR transmittance profile
  real(r8),allocatable ::  StemClumpFactor_pft(:)                 !stem clumping factor
  real(r8),allocatable ::  DStemSWabsorptivityL_pft(:,:)          !standing dead stem shortwave radiation absorptivity profile
  real(r8),allocatable ::  DStemPARabsorptivityL_pft(:,:)         !standing dead stem PAR absorptivity profile
  real(r8),allocatable ::  RadSWDStemAlbedoL_pft(:,:)             !standing dead stem shortwave radiation scattering albedo profile
  real(r8),allocatable ::  RadPARDstemAlbedoL_pft(:,:)            !standing dead stem PAR scattering albedo profile
  real(r8),allocatable ::  RadSWDstemTransmitanceL_pft(:,:)       !standing dead stem SW transmittance profile
  real(r8),allocatable ::  RadPARDstemTransmitanceL_pft(:,:)      !standing dead stem PAR transmittance profile

  public :: InitRadiationData
  public :: DestroyRadiationData
 contains
!------------------------------------------------------------------------------------------
  subroutine InitRadiationData
  implicit none
  allocate(LeafSWabsorptivityL_pft(NumCanopyLayers1,JP1));LeafSWabsorptivityL_pft=0._r8
  allocate(LeafPARabsorptivityL_pft(NumCanopyLayers1,JP1));LeafPARabsorptivityL_pft=0._r8
  allocate(RadSWLeafAlbedoL_pft(NumCanopyLayers1,JP1)); RadSWLeafAlbedoL_pft=0._r8
  allocate(RadPARLeafAlbedoL_pft(NumCanopyLayers1,JP1)); RadPARLeafAlbedoL_pft=0._r8
  allocate(RadSWLeafTransmitanceL_pft(NumCanopyLayers1,JP1)); RadSWLeafTransmitanceL_pft=0._r8
  allocate(RadPARLeafTransmitanceL_pft(NumCanopyLayers1,JP1)); RadPARLeafAlbedoL_pft=0._r8

  allocate(StemSWabsorptivityL_pft(NumCanopyLayers1,JP1)); StemSWabsorptivityL_pft=0._r8
  allocate(StemPARabsorptivityL_pft(NumCanopyLayers1,JP1)); StemPARAbsorptivityL_pft=0._r8
  allocate(RadSWStemAlbedoL_pft(NumCanopyLayers1,JP1)); RadSWStemAlbedoL_pft=0._r8
  allocate(RadPARStemAlbedoL_pft(NumCanopyLayers1,JP1)); RadPARStemAlbedoL_pft=0._r8
  allocate(RadSWStemTransmitanceL_pft(NumCanopyLayers1,JP1)); RadSWStemTransmitanceL_pft=0._r8
  allocate(RadPARStemTransmitanceL_pft(NumCanopyLayers1,JP1)); RadPARStemTransmitanceL_pft=0._r8
  allocate(StemClumpFactor_pft(JP1)); StemClumpFactor_pft=0._r8

  allocate(DStemSWabsorptivityL_pft(NumCanopyLayers1,JP1));      DStemSWabsorptivityL_pft=0._r8
  allocate(DStemPARabsorptivityL_pft(NumCanopyLayers1,JP1));    DStemPARabsorptivityL_pft=0._r8       
  allocate(RadSWDStemAlbedoL_pft(NumCanopyLayers1,JP1));      RadSWDStemAlbedoL_pft=0._r8         
  allocate(RadPARDstemAlbedoL_pft(NumCanopyLayers1,JP1));     RadPARDstemAlbedoL_pft=0._r8         
  allocate(RadSWDstemTransmitanceL_pft(NumCanopyLayers1,JP1));RadSWDstemTransmitanceL_pft=0._r8         
  allocate(RadPARDstemTransmitanceL_pft(NumCanopyLayers1,JP1)); RadPARDstemTransmitanceL_pft=0._r8       

  end subroutine InitRadiationData
!------------------------------------------------------------------------------------------
  subroutine DestroyRadiationData
  use abortutils,    only: destroy  
  implicit none
  
  call destroy(LeafSWabsorptivityL_pft)
  call destroy(LeafPARabsorptivityL_pft)
  call destroy(RadSWLeafAlbedoL_pft)
  call destroy(RadPARLeafAlbedoL_pft)
  call destroy(RadSWLeafTransmitanceL_pft)
  call destroy(RadPARLeafTransmitanceL_pft)

  call destroy(StemSWabsorptivityL_pft)
  call destroy(StemPARabsorptivityL_pft)
  call destroy(RadSWStemAlbedoL_pft)
  call destroy(RadPARStemAlbedoL_pft)
  call destroy(RadSWStemTransmitanceL_pft)
  call destroy(RadPARStemTransmitanceL_pft)
  call destroy(StemClumpFactor_pft)

  call destroy(DStemSWabsorptivityL_pft)      
  call destroy(DStemPARabsorptivityL_pft)     
  call destroy(RadSWDStemAlbedoL_pft)         
  call destroy(RadPARDstemAlbedoL_pft)        
  call destroy(RadSWDstemTransmitanceL_pft)   
  call destroy(RadPARDstemTransmitanceL_pft)  

  end subroutine DestroyRadiationData

end module RadiationDataMod