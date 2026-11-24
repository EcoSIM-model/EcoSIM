module PlantDisturbMod
!
!! Description:
! code to apply distance to plants
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use SOMDataType
  use PlantBGCPars
  use PlantTraitDataType
  use GridConsts
  use FlagDataType
  use SoilWaterDataType
  use EcoSIMCtrlDataType
  use PlantDataRateType
  use FertilizerDataType
  use ClimForcDataType
  use EcosimConst
  use PlantMgmtDataType
  use RootDataType
  use CanopyDataType
  use EcoSimSumDataType
  use SoilBGCDataType
  use EcosimBGCFluxType
  use GridDataType
  implicit none
  private
  character(len=*),private, parameter :: mod_filename = &
  __FILE__
! end_include_section

  public :: PrepLandscapeGrazing
  contains

!------------------------------------------------------------------------------------------

  subroutine PrepLandscapeGrazing(I,J,NHW,NHE,NVN,NVS)
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NX,NY,NZ,nx1,NY1
  real(r8) :: TotalShootC4Grazing  !total canopy C that is subject to grazing/herbivory
  real(r8) :: GrazingArea !total area of active cells in the grazing section [m2]
!     begin_execution

  D2995: DO NX=NHW,NHE
    D2990: DO NY=NVN,NVS
      D2985: DO NZ=1,NP_col(NY,NX)
!
!     iHarvstType_pft=harvest type:0=none,1=grain,2=all above-ground
!                       ,3=pruning,4=grazing,5=fire,6=herbivory
!     LSG=landscape grazing section number
!     AvgCanopyBiomC2Graze_pft=section mean biomass density times local cell area
!       different pft can undergo different treatment/disturbance
        IF(iHarvstType_pft(NZ,I,NY,NX).EQ.iharvtyp_grazing .OR. iHarvstType_pft(NZ,I,NY,NX).EQ.iharvtyp_herbivo)THEN
          TotalShootC4Grazing = 0
          GrazingArea = 0._r8
          D1995: DO NX1=NHW,NHE
            D1990: DO NY1=NVN,NVS
              !when grazing occurs, it occurs to the same plant species over all grids of the landscape
              IF(LSG_pft(NZ,NY1,NX1).EQ.LSG_pft(NZ,NY,NX))THEN
                IF(IsPlantActive_pft(NZ,NY1,NX1).EQ.iTrue)THEN
                  TotalShootC4Grazing = TotalShootC4Grazing+ShootElms_pft(ielmc,NZ,NY1,NX1)
                  GrazingArea = GrazingArea+AREA_3D(3,NU_col(NY1,NX1),NY1,NX1)
                ENDIF
              ENDIF
            ENDDO D1990
          ENDDO D1995

          IF(GrazingArea.GT.0._r8)THEN
            ! Keep cell-total units: the downstream local/mean biomass ratio
            ! then compares densities, without counting local area twice.
            AvgCanopyBiomC2Graze_pft(NZ,NY,NX)=TotalShootC4Grazing/GrazingArea*AREA_3D(3,NU_col(NY,NX),NY,NX)
          ELSE
            AvgCanopyBiomC2Graze_pft(NZ,NY,NX)=ShootElms_pft(ielmc,NZ,NY,NX)
          ENDIF
        ENDIF
      ENDDO D2985
    ENDDO D2990
  ENDDO D2995
  end subroutine PrepLandscapeGrazing
end module PlantDisturbMod
