module LandSurfDataType
  ! Description:
  ! Store column surface geometry, aerodynamic exchange properties and air
  ! endpoints, including the hourly diagnosed exposure of the plant canopy.
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use GridConsts
  use TracerIDMod, only : idg_beg,idg_end
  implicit none
  public
  save
  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  real(r8) :: ALTIG                                                      !Altitude of landscape, [m]
  logical,allocatable :: canopy_exposed_col(:,:) ! canopy area protrudes above snow, [NY,NX], [-]
  real(r8),target,allocatable ::  SoilSurfRoughness_col(:,:)           !Initial soil surface roughness height, [m]
  real(r8),target,allocatable ::  ZeroPlaneDisplacem_col(:,:)            !Zero plane displacement height, [m]
  real(r8),target,allocatable ::  RoughnessLength_col(:,:)                   !Canopy surface roughness height, [m]
  real(r8),target,allocatable ::  SoiSurfRoughness(:,:)                  !Soil surface roughness height for calculating runoff velocity, [m]
  real(r8),target,allocatable ::  WindMesureHeight_col(:,:)              !Wind speed measurement height, [m]
  real(r8),target,allocatable ::  SurfElevation_col(:,:)                 !Reference surface elevation of grid cell, [m]
  real(r8),target,allocatable ::  RawIsoTSurf2CanopyHScal_col(:,:)       !scalar for isothermal aerodynamic resistance between zero-sink height and ground surface, [h m-1]
  real(r8),target,allocatable ::  RawIsoTAtm2CanopySinkZ_col(:,:)        !isothermal aerodynamic resistance between zero-sink height and wind ref height in atmosphere, [h m-1]
  real(r8),target,allocatable ::  RawTAtm2CanopySinkZ_col(:,:)           !aerodynamic resistance between zero-sink height and wind ref height in atmosphere, [h m-1]
  real(r8),target,allocatable ::  RawCanopyH2SinkZ_col(:,:)              !isothermal aerodynamic resistance bewtween canopy height and zero sink height, [h m-1]
  real(r8),target,allocatable ::  RIB_col(:,:)                           !Richardson number for calculating boundary layer resistance, [-]
  real(r8),target,allocatable ::  AltitudeLandscape_col(:,:)                          !Altitude of landscape, [m]
  real(r8),target,allocatable ::  SineGrndSlope_col(:,:)                 !Sine of slope, [-]
  real(r8),target,allocatable ::  CosineGrndSlope_col(:,:)               !Cosine of slope, [-]
  real(r8),target,allocatable ::  GroundSurfaceAzimuth_col(:,:)             !Azimuth of slope, [-]
  real(r8),target,allocatable ::  RefSurfElevation_col(:,:)                          !Altitude, [m]
  real(r8),target,allocatable ::  SL_col(:,:)                            !Slope, [degree]
  real(r8),target,allocatable ::  ASP_col(:,:)                           !Aspect , [degree]
  real(r8),target,allocatable ::  VPQ_col(:,:)                           ! air pressure at canopy sink height
  real(r8),target,allocatable ::  TKQ_col(:,:)                           ! air temperature at canopy sink height [K]
!----------------------------------------------------------------------

contains
  subroutine InitLandSurfData
  ! Description:
  ! Allocate surface fields; canopy exposure is rebuilt by hourly preparation.
  implicit none
  allocate(canopy_exposed_col(JY,JX)); canopy_exposed_col=.false.
  allocate(RawTAtm2CanopySinkZ_col(JY,JX)); RawTAtm2CanopySinkZ_col=0._r8
  allocate(RawCanopyH2SinkZ_col(JY,JX)); RawCanopyH2SinkZ_col=0._r8
  allocate(RawIsoTSurf2CanopyHScal_col(JY,JX)) ; RawIsoTSurf2CanopyHScal_col=0._r8
  allocate(SoilSurfRoughness_col(JY,JX));          SoilSurfRoughness_col=0._r8
  allocate(ZeroPlaneDisplacem_col(JY,JX));          ZeroPlaneDisplacem_col=0._r8
  allocate(RoughnessLength_col(JY,JX));          RoughnessLength_col=0._r8
  allocate(SoiSurfRoughness(JY,JX));          SoiSurfRoughness=0._r8
  allocate(WindMesureHeight_col(JY,JX));          WindMesureHeight_col=0._r8
  allocate(SurfElevation_col(JY,JX));         SurfElevation_col=0._r8
  allocate(RawIsoTAtm2CanopySinkZ_col(JY,JX));         RawIsoTAtm2CanopySinkZ_col=0._r8
  allocate(RIB_col(JY,JX));         RIB_col=0._r8
  allocate(AltitudeLandscape_col(JY,JX));        AltitudeLandscape_col=0._r8
  allocate(SineGrndSlope_col(JY,JX));        SineGrndSlope_col=0._r8
  allocate(CosineGrndSlope_col(JY,JX));        CosineGrndSlope_col=0._r8
  allocate(GroundSurfaceAzimuth_col(JY,JX));        GroundSurfaceAzimuth_col=0._r8
  allocate(RefSurfElevation_col(JY,JX));        RefSurfElevation_col=0._r8
  allocate(SL_col(JY,JX));          SL_col=0._r8
  allocate(ASP_col(JY,JX));         ASP_col=0._r8
  allocate(VPQ_col(JY,JX));         VPQ_col=0._r8  
  allocate(TKQ_col(JY,JX));         TKQ_col=0._r8  
  end subroutine InitLandSurfData

!----------------------------------------------------------------------
  subroutine DestructLandSurfData
  ! Description:
  ! Release the column surface geometry and exchange fields.
  use abortutils, only : destroy
  implicit none
  call destroy(canopy_exposed_col)
  call destroy(VPQ_col)    
  call destroy(TKQ_col)  
  call destroy(SoilSurfRoughness_col)
  call destroy(ZeroPlaneDisplacem_col)
  call destroy(RoughnessLength_col)
  call destroy(SoiSurfRoughness)
  call destroy(WindMesureHeight_col)
  call destroy(SurfElevation_col)
  call destroy(RawIsoTAtm2CanopySinkZ_col)
  call destroy(RIB_col)
  call destroy(AltitudeLandscape_col)
  call destroy(SineGrndSlope_col)
  call destroy(CosineGrndSlope_col)
  call destroy(GroundSurfaceAzimuth_col)
  call destroy(RefSurfElevation_col)
  call destroy(SL_col)
  call destroy(ASP_col)
  call destroy(RawCanopyH2SinkZ_col)
  call destroy(RawIsoTSurf2CanopyHScal_col)
  call destroy(RawTAtm2CanopySinkZ_col)
  end subroutine DestructLandSurfData

end module LandSurfDataType
