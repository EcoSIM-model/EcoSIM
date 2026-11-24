module InitVegBGC

  use EcosimConst
  use GridConsts
  use PlantBGCPars, only : ibackward,iforward
  implicit none
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
  public :: InitIrradianceGeometry
  contains
  ![header]
!------------------------------------------------------------------------------------------
  subroutine InitIrradianceGeometry(YSIN,YCOS,SkyAzimuthAngle)
  use CanopyRadDataType
  implicit none

  real(r8), intent(out) :: YSIN(NumOfSkyAzimuthSects)
  real(r8), intent(out) :: YCOS(NumOfSkyAzimuthSects)
  real(r8), intent(out) :: SkyAzimuthAngle(NumOfSkyAzimuthSects)
  real(r8) :: ZAZI(NumOfLeafAzimuthSectors)
  real(r8) :: YAGL    !sector angle of sky azimuth
  real(r8) :: OMEGZ,OMEGY
  real(r8) :: ZAGL,DAZI
  real(r8) :: DSkyAzimuthAngle
  integer :: L,M,N

  !     begin_execution
  !     SineLeafAngle,CosineLeafAngle=sine,cosine of leaf inclination class
  !     ZAZI=leaf azimuth class, it is pi because only on one side of the leaf
  !     SkyAzimuthAngle,YSIN,YCOS=sky azimuth,sine,cosine of sky azimuth
  !     OMEGA,OMEGX=incident aNGLe of diffuse radn at leaf,horizontal surface
  !     iScatteringDiffus:1=backscattering,2=forward scattering of sky radiation
  !

  D205: DO L=1,NumOfLeafAzimuthSectors
    ZAZI(L)=(L-0.5_r8)*PICON/real(NumOfLeafAzimuthSectors,r8)
  ENDDO D205

  !NumOfSkyAzimuthSects: number of sky azimuth sectors
  !NumOfLeafAzimuthSectors: number of leaf azimuth sectors
  !the following configuration is specific to NumOfSkyAzimuthSects equal to 4, so that the whole sky is separated into 4 zones.
  !in Grant 1989, six sky azimuth angles are used, so that the whole sky is separated into 10 zones.
  ! 
  DSkyAzimuthAngle=2._r8*PICON/real(NumOfSkyAzimuthSects,r8)
  YAGL           = PICON/real(NumOfSkyAzimuthSects,r8)
  D230: DO N=1,NumOfSkyAzimuthSects
    SkyAzimuthAngle(N)   = (N-0.5_r8)*DSkyAzimuthAngle    
    !The below is specifc for N=4, See Grant et al. (1989) for 10 sky zone.
    YSIN(N)              = SIN(YAGL)
    YCOS(N)              = COS(YAGL)
    TotSineSkyAngles_grd = TotSineSkyAngles_grd+YSIN(N)

    D225: DO L=1,NumOfLeafAzimuthSectors
      DAZI=COS(ZAZI(L)-SkyAzimuthAngle(N))
      DO  M=1,NumLeafInclinationClasses
        !cos(\sigma_{L,M,N}) in eq. (13) in Grant 1989
        OMEGY             = CosineLeafAngle(M)*YSIN(N)+SineLeafAngle(M)*YCOS(N)*DAZI
        OMEGA2Leaf(N,M,L) = ABS(OMEGY)
        OMEGX(N,M,L)      = OMEGA2Leaf(N,M,L)/YSIN(N)

        IF(CosineLeafAngle(M).GT.YSIN(N))THEN
          OMEGZ=ACOS(OMEGY)
        ELSE
          OMEGZ=-ACOS(OMEGY)
        ENDIF

        IF(OMEGZ.GT.-PICON2h)THEN
          ZAGL=YAGL+2.0_r8*OMEGZ
        ELSE
          ZAGL=YAGL-2.0_r8*(PICON+OMEGZ)
        ENDIF

        IF(ZAGL.GT.0.0_r8 .AND. ZAGL.LT.PICON)THEN
          iScatteringDiffus(N,M,L)=ibackward
        ELSE
          iScatteringDiffus(N,M,L)=iforward
        ENDIF
      ENDDO
    ENDDO D225
  ENDDO D230
  end subroutine InitIrradianceGeometry
  ![tail]
end module InitVegBGC
