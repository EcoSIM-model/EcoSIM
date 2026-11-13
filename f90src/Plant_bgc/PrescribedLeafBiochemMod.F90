module PrescribedLeafBiochemMod
  ! Description:
  ! Supply LFPRTC and prescribed leaf geometry to StomatalDynamics.
  ! SP uses one diagnostic leaf node spanning its prescribed canopy layers
  ! for the potential fixation calculation that determines stomatal resistance.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use DebugToolMod, only: PrintInfo
  use ElmIDMod, only: iTrue, iFalse, ielmc, iphenotyp_evgreen
  use PlantAPICommonData
  use PlantSiteAPIData, only: plt_site
  use PlantMorphologyAPIData, only: plt_morph
  use PlantBiomassAPIData, only: plt_biom
  use PlantPhenologyAPIData, only: plt_pheno
  use PlantPhotosynthesisAPIData, only: plt_photo
  use PlantEnergyWaterAPIData, only: plt_ew
  use StomatesMod, only: PhotosynsDiag
  implicit none
  private
  public :: PreparePrescribedLeafBiochemistry
contains
!------------------------------------------------------------------------------------------
  subroutine PreparePrescribedLeafBiochemistry(I,J,NZ,ProteinDensity,LeafAreaByLayer)
    ! Description:
    ! Rebuild the SP diagnostic leaf node from LFPRTC and prescribed leaf area.
    ! Prescribed protein supplies nutrient status; no extra N/P or growth-phase
    ! inhibition is inferred from the separately prescribed osmotic concentration.
    integer, intent(in) :: I,J ! day of year and hour [-]
    integer, intent(in) :: NZ ! PFT slot in the current column [-]
    real(r8), intent(in) :: ProteinDensity ! LFPRTC [gC protein m-2 leaf]
    real(r8), intent(in) :: LeafAreaByLayer(NumCanopyLayers1) ! leaf area [m2/cell] by canopy layer
    character(len=*), parameter :: subname='PreparePrescribedLeafBiochemistry'
    real(r8) :: LeafArea
    integer :: L

    call PrintInfo('beg '//subname)
    associate( &
      NumBranches => plt_morph%NumOfBranches_pft, & !output: diagnostic branch count by PFT [-]
      NodeArea => plt_morph%LeafArea_node, & !output: leaf area by node/branch/PFT [m2/cell]
      LayerArea => plt_morph%CanopyLeafArea_lnode, & !output: leaf area by layer/node/branch/PFT [m2/cell]
      SectorArea => plt_morph%LeafAreaZsec_brch, & !output: leaf area by inclination/layer/node/branch/PFT [m2/cell]
      EffectiveArea => plt_photo%LeafEffArea_zsec, & !output: clumped sector leaf area [m2/cell]
      LeafAngles => plt_morph%LeafAngleClass_pft, & !input: inclination fractions by class/PFT [-]
      ClumpFactor => plt_morph%ClumpFactorNow_pft, & !input: current clumping fraction by PFT [-]
      Height => plt_morph%CanopyHeightLive_pft, & !input: live canopy height by PFT [m]
      SnowDepth => plt_ew%SnowDepth, & !input: current snow height [m]
      ZERO => plt_site%ZERO, & !input: numerical height tolerance [m]
      Population => plt_site%PlantPopuLive_pft, & !input: living plants by PFT [plants/cell]
      SLA => plt_morph%SLA1_pft, & !input: specific leaf area by PFT [m2 leaf/gC]
      LeafElements => plt_biom%LeafElmntNode_brch, & !output: diagnostic leaf elements by element/node/branch/PFT [g/cell]
      Protein => plt_biom%LeafProteinC_node, & !output: diagnostic protein C by node/branch/PFT [gC/cell]
      NonstructConc => plt_biom%LeafPetoNonstElmConc_brch, & !output: unused SP nutrient feedback by element/branch/PFT [g/gC]
      PhenologyType => plt_pheno%iPlantPhenolType_pft, & !output: SP activity follows prescribed LAI by PFT [-]
      BranchAlive => plt_pheno%isPlantBranchAlive_brch, & !output: diagnostic leaf branch eligibility by branch/PFT [-]
      GrainFailure => plt_pheno%HourFailGrainFill_brch, & !output: unused SP growth-phase hours by branch/PFT [h]
      NutrientControl => plt_photo%NutrientCtrlonC4Carboxy_node, & !output: prescribed protein nutrient sufficiency by node/branch/PFT [-]
      RubiscoCapacity => plt_photo%Vmax4RubiscoCarboxy_node, & !output: node/branch/PFT capacity [umol m-2 leaf s-1]
      PEPCapacity => plt_photo%Vmax4PEPCarboxy_node, & !output: node/branch/PFT capacity [umol m-2 leaf s-1]
      C4Pool => plt_photo%CPOOL4_node, & !output: unused SP mesophyll metabolite by node/branch/PFT [gC/cell]
      BundleCO2 => plt_photo%CMassCO2BundleSheath_node & !output: unused SP bundle-sheath CO2 by node/branch/PFT [gC/cell]
    )
      NumBranches(NZ)=0
      NodeArea(:,:,NZ)=0._r8
      LayerArea(:,:,:,NZ)=0._r8
      SectorArea(:,:,:,:,NZ)=0._r8
      EffectiveArea(:,:,:,:,NZ)=0._r8
      LeafElements(:,:,:,NZ)=0._r8
      Protein(:,:,NZ)=0._r8
      NonstructConc(:,:,NZ)=0._r8
      BranchAlive(:,NZ)=iFalse
      GrainFailure(:,NZ)=0._r8
      NutrientControl(:,:,NZ)=1._r8
      RubiscoCapacity(:,:,NZ)=0._r8
      PEPCapacity(:,:,NZ)=0._r8
      ! The stomatal potential-fixation calculation has no SP growth or C4
      ! metabolite history. Disable product inhibition from these unused pools
      ! rather than interpreting API fill values as accumulated carbon.
      C4Pool(:,:,NZ)=0._r8
      BundleCO2(:,:,NZ)=0._r8
      PhenologyType(NZ)=iphenotyp_evgreen
      LeafArea=sum(LeafAreaByLayer)
      if(Population(NZ)>0._r8 .and. LeafArea>0._r8)then
        NumBranches(NZ)=1
        BranchAlive(1,NZ)=iTrue
        NodeArea(1,1,NZ)=LeafArea
        LayerArea(:,1,1,NZ)=LeafAreaByLayer
        if(SLA(NZ)>0._r8)LeafElements(ielmc,1,1,NZ)=LeafArea/SLA(NZ)
        Protein(1,1,NZ)=ProteinDensity*LeafArea
        if(Height(NZ)>=SnowDepth-ZERO)then
          do L=1,NumCanopyLayers1
            SectorArea(:,L,1,1,NZ)=LeafAngles(:,NZ)*LeafAreaByLayer(L)/real(NumOfLeafAzimuthSectors1,r8)
          enddo
          EffectiveArea(:,:,1,1,NZ)=SectorArea(:,:,1,1,NZ)*ClumpFactor(NZ)
        endif
      endif
      call PhotosynsDiag(I,J,NZ)
    end associate
    call PrintInfo('end '//subname)
  end subroutine PreparePrescribedLeafBiochemistry
end module PrescribedLeafBiochemMod
