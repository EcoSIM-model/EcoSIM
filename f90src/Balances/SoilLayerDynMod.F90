module SoilLayerDynMod
! Description:
! subroutines to do soil relayering
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod,  only: micpar
  use PlantMgmtDataType, only: NP_col, NP0_col
  use minimathmod,       only: AZMAX1
  use UnitMod,           only: units
  use EcoSIMConfig,      only: ndbiomcp => NumDeadMicrbCompts,NumMicbAFunGrupsPerCmplx,NumMicbHFunGrupsPerCmplx
  use PlantTraitDataType, only : Myco_pft
  use SoilPhysParaMod   , only : IsPondLayer
  use NumericalAuxMod
  use DebugToolMod
  USE RedistDataMod      
  use RootDataType
  use PlantDataRateType, only: trcg_root_vr
  use GridDataType
  USE EcoSIMCtrlDataType
  USE SoilWaterDataType
  use SoilPropertyDataType
  use FlagDataType
  use SOMDataType
  use MicrobialDataType
  USE SoilBGCDataType
  use AqueChemDatatype
  use FertilizerDataType
  use SoilHeatDataType
  use EcoSIMCtrlMod
  use SoilPhysDataType
  use SurfSoilDataType
  use SurfLitterDataType
  use EcosimConst
implicit none


  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  real(r8) :: WDPOBDL,WDNOBD1,WDPOBD0,WDPOBD1
  real(r8) :: WDNHBDL,WDNHBD0,WDNHBD1,WDNOBDL,WDNOBD0
  real(r8), PARAMETER :: ZEROC=0.1E-03_r8

  integer, parameter :: isl_undef   = 0   !layer beneath is not defined
  integer, parameter :: isl_water   = 1   !layer beneath is water
  integer, parameter :: isl_soil    = 2   !layer beneath is soil
 
  integer, parameter :: ich_watlev    = 1
  integer, parameter :: iTopLayShrink = 2
  integer, parameter :: iTopLayGrow   = 3
  integer, parameter :: ich_frzthaw   = 4
  integer, parameter :: ich_erosion   = 5
  integer, parameter :: ich_sombgc    = 6

  integer, parameter :: iInnerLayer   = 1

  public :: UpdateSoilGrids
  contains


!------------------------------------------------------------------------------------------

  subroutine UpdateSoilGrids(I,J,NY,NX,DORGC_vr,DVLiceMicP_vr)
  ! Description:
  ! Relayer the soil profiles.
  !litter layer maybe mixed into the soil layers
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: DORGC_vr(JZ)          !change in organic matter, final-initial
  REAL(R8),INTENT(IN) :: DVLiceMicP_vr(JZ)  !change in ice volume, final-initial

  character(len=*), parameter :: subname = 'UpdateSoilGrids'
  real(r8) :: SoilLayBotEdgeLocNew_vr(0:JZ)   !edge location after mass update of different processes
  real(r8) :: SoilLayBotEdgeLocOld_vr(JZ)     !edge location before mass update
  integer  :: IFLGL(0:JZ,6)                !flag for soil thickness change by different processes
  real(r8) :: DDLYRX(3)
  integer :: NN,K,M,N,NR,NZ,L,idg
  integer :: L_rcept,L_donor,NUX,ICHKL,NGL,NU_old
  real(r8) :: FX,FY
  real(r8) :: FBO
  real(r8) :: FHO,FWO,FXVOLW
  real(r8) :: FO
  real(r8) :: DDLYRY(JZ)                   !change in layer thickness compared to the initial layer thickness
  real(r8) :: ZSOI(0:JZ)
  real(r8) :: ENGY0X,ENGY0,ENGY1X,ENGY1
  real(r8) :: XVOLWP   !volume of mobile water in litter layer
  real(r8) :: twat
  !     begin_execution
  !     SOIL SUBSIDENCE
  !
  NU_old=NU_col(NY,NX)
  
  call PrintInfo('beg '//subname)

  ZSOI=CumSoilThickness_vr(0:JZ,NY,NX)

  IF(.not. erosion_model)return
  !soil relayering can occur due to freeze-thaw, soc change, and erosion
  !set IFLGL(:,[1,4,5,6])
  call UpdateLayerEdges(I,J,NY,NX,DORGC_vr,DVLiceMicP_vr,SoilLayBotEdgeLocOld_vr,SoilLayBotEdgeLocNew_vr,IFLGL)


  !
  !     RECALCULATE SOIL LAYER THICKNESS
  !
  twat=0._r8
  DO L=NU_col(NY,NX),NL_col(NY,NX)
    twat=twat+VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)+(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))*DENSICE
  ENDDO


  ICHKL=0
  D245: DO L=NU_col(NY,NX),NL_col(NY,NX)-1
    !top -> down
    D230: DO NN=1,3

      If(NN.EQ.iInnerLayer)then
        call DiagInnerLayerThickChange(I,J,L,NY,NX,SoilLayBotEdgeLocOld_vr,IFLGL(:,ich_watlev),ICHKL,NUX,SoilLayBotEdgeLocNew_vr,DDLYRX(NN),DDLYRY,XVOLWP)        
      else !about pond surface change
        !IFLGL(:,2) set for drying pond, IFLGL(:,3) for increasing pond
        call DiagPondTopThickChange(I,J,L,NN,NY,NX,NUX,SoilLayBotEdgeLocNew_vr,IFLGL,DDLYRX,XVOLWP)
      endif
      !
      !     TRANSFER STATE VARIABLES BETWEEN LAYERS
      !      
      IF(ABS(DDLYRX(NN)).GT.ZERO)THEN   !layer L changed
        !L0,L1: source and target layers
        call ComputeLayerAllocFs(I,J,L,NN,NY,NX,NUX,DDLYRX,IFLGL,FX,FO,L_donor,L_rcept)
        
        IF(FX.GT.ZERO)THEN
          iResetSoilProf_col(NY,NX)=itrue
          FY=1.0_r8-FX
          IF(FY.LE.ZERO2)FY=0.0_r8

          ! Layer L0 is pond
          IF(IsPondLayer(SoilBulkDensity_vr(L_donor,NY,NX)))THEN !
            call UpdateLayerMaterials(L,L_donor,L_rcept,NY,NX,FX,FY,IFLGL)
  
            !layer L0 is soil  
          ELSE
            ! MOVE ALL MATTER WITH CHANGES IN LAYER DEPTHS
            !
            FHO=0._r8 ! Surface litter has no macropore transfer.
            IF(L_donor.NE.0)THEN
              call MoveFertMinerals(L,L_donor,L_rcept,NY,NX,FX,FO,FHO)
            ENDIF

            call MoveHeatWat(I,J,L,L_donor,L_rcept,NY,NX,FO,FX,XVOLWP,FHO,FWO)
            !
            !     SOIL FERTILIZER
            !
            call MoveFertSalt(L,L_donor,L_rcept,NY,NX,FX,FWO)
            !
            !     SOIL SOLUTES IN BAND
            !
            IF(L_donor.NE.0)THEN
              call MoveBandSolute(L,L_donor,L_rcept,NY,NX,FX,FWO,FO)

              call MoveDisolvGas(L_donor,L_rcept,NY,NX,FX,FWO)

              call MoveMacPoreSolute(L_donor,L_rcept,NY,NX,FHO)
            ENDIF
            ! SOIL ORGANIC MATTER
            call MoveSOM(L_donor,L_rcept,L,NY,NX,FO,IFLGL)

          ENDIF
          IF(NN.EQ.iInnerLayer)THEN
            IF(SoilBulkDensity_vr(L_donor,NY,NX).LE.ZERO .AND. SoilBulkDensity_vr(L_rcept,NY,NX).LE.ZERO &   !both donor and reception layers are both water
              .AND. VLWatMicP_vr(L_donor,NY,NX)+VLiceMicP_vr(L_donor,NY,NX).LE.ZEROS(NY,NX))THEN             !donor layer has not water, or no thickness
              CumDepz2LayBottom_vr(L_rcept,NY,NX) = CumDepz2LayBottom_vr(L_donor,NY,NX)
              SoilLayBotEdgeLocNew_vr(L_rcept)    = SoilLayBotEdgeLocNew_vr(L_donor)
            ENDIF
          ENDIF

          !
          !     RESET LOWER LAYER NUMBER WITH EROSION
          !
          IF(iErosionMode.EQ.ieros_frzthaweros .OR. iErosionMode.EQ.ieros_frzthawsomeros)THEN
            IF(L.EQ.NL_col(NY,NX) .AND. DLYR_3D(3,L,NY,NX).GT.DLYRI_3D(3,L,NY,NX))THEN
              NL_col(NY,NX)=MIN(NLI_col(NY,NX),NL_col(NY,NX)+1)
            ENDIF
            IF(L.EQ.NL_col(NY,NX)-1 .AND. CumDepz2LayBottom_vr(NL_col(NY,NX),NY,NX)-CumDepz2LayBottom_vr(L,NY,NX).LE.ZEROC)THEN
              CumDepz2LayBottom_vr(L,NY,NX)             = CumDepz2LayBottom_vr(L,NY,NX)+DLYR_3D(3,NL_col(NY,NX),NY,NX)
              SoilLayBotEdgeLocNew_vr(L)                = SoilLayBotEdgeLocNew_vr(L)+DLYR_3D(3,NL_col(NY,NX),NY,NX)
              CumDepz2LayBottom_vr(NL_col(NY,NX),NY,NX) = CumDepz2LayBottom_vr(L,NY,NX)
              SoilLayBotEdgeLocNew_vr(NL_col(NY,NX))    = SoilLayBotEdgeLocNew_vr(L)
              DLYR_3D(3,NL_col(NY,NX),NY,NX)            = 0.0_r8
              NL_col(NY,NX)                             = L
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D230
  ENDDO D245


  call RemoveNullPondLayer(I,J,NY,NX,IFLGL,NU_old)

  if(NU_col(NY,NX)/=NU_old)then
    ! Every active depth uses the current surface as its origin. Retired slots
    ! have no thickness; their coincident edges must also follow that surface.
    CumDepz2LayBottom_vr(1:NU_col(NY,NX)-1,NY,NX)=CumDepz2LayBottom_vr(0,NY,NX)
    CumSoilThickness_vr(0:NU_col(NY,NX)-1,NY,NX)=0._r8
    CumSoilThickMidL_vr(1:NU_col(NY,NX)-1,NY,NX)=0._r8
    DO L=NU_col(NY,NX),NL_col(NY,NX)
      CumSoilThickness_vr(L,NY,NX)=CumDepz2LayBottom_vr(L,NY,NX)-CumDepz2LayBottom_vr(0,NY,NX)
      CumSoilThickMidL_vr(L,NY,NX)=0.5_r8*(CumSoilThickness_vr(L-1,NY,NX)+CumSoilThickness_vr(L,NY,NX))
      SoilDepthMidLay_vr(L,NY,NX)=0.5_r8*(CumDepz2LayBottom_vr(L-1,NY,NX)+CumDepz2LayBottom_vr(L,NY,NX))
    ENDDO
  endif

  twat=0._r8
  DO L=NU_col(NY,NX),NL_col(NY,NX)
    twat=twat+VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)+(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))*DENSICE
    if(L+1<=NL_col(NY,NX))then
      if(CumSoilThickness_vr(L,NY,NX)>CumSoilThickness_vr(L+1,NY,NX))then
        print*,I*1000+J,'Ny,NX',NY,NX
        print*,'thick',CumSoilThickness_vr(1:NL_col(NY,NX),NY,NX),L,NU_col(NY,NX),NU_old
        print*,'cdl',CumDepz2LayBottom_vr(1:NL_col(NY,NX),NY,NX)
        print*,'DL',DLYR_3D(3,NU_col(NY,NX):NL_col(NY,NX),NY,NX)
!        call endrun('soil thickness update problematic in '//trim(mod_filename)//' at line',__LINE__)          
      endif
    endif
  ENDDO

!  if(I==348 .and. J==4)call endrun('soil thickness update problematic in '//trim(mod_filename)//' at line',__LINE__)          
  call PrintInfo('end '//subname)
  end subroutine UpdateSoilGrids
!------------------------------------------------------------------------------------------
  subroutine  RemoveNullPondLayer(I,J,NY,NX,IFLGL,NU_old)
  ! Description:
  ! Retire a thin pond layer only after conserving its remaining material pools;
  ! then shift overlying layers and advance the active surface index.
  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer, intent(in) :: IFLGL(0:JZ,6)                !flag for soil thickness change by different processes  
  integer, intent(in) :: NU_old
  character(len=*), parameter :: subname='RemoveNullPondLayer'
  integer :: L,L1,NZ,L2,ids,Ltarget
  integer :: MoveFlags(0:JZ,6) ! full-layer material transfer flags [-]
  real(r8) :: FX,FY
  
  call PrintInfo('beg '//subname)
  !identify missing layer 
  L1=0
  DO L=NL_col(NY,NX),NU_old,-1
    if(DLYR_3D(3,L,NY,NX)*cpw<VHCPNX_col(NY,NX) .and. IsPondLayer(SoilBulkDensity_vr(L,NY,NX)))then
      L1=L    
      EXIT
    endif
  ENDDO

  ! Layer retirement/index shifts move complete inventories, rather than the
  ! water-only transfer used when the surface pond grows over litter.
  IF(L1>=NU_old)THEN
    MoveFlags=IFLGL
    MoveFlags(:,iTopLayGrow)=isl_undef
    Ltarget=L1+1
    if(Ltarget>NL_col(NY,NX))Ltarget=L1-1
    if(Ltarget<NU_old)then
      write(iulog,*)subname,': no receiving layer for retired pond',NY,NX,L1
      call endrun('Cannot retire the only remaining pond layer',__LINE__)
    endif
    ! Hydrology already moved its residual water/heat. Gas, solute, organic
    ! and root pools must follow before this layer is excluded from budgets.
    call UpdateLayerMaterials(L1,L1,Ltarget,NY,NX,1._r8,0._r8,MoveFlags)
    L2=L1-1
    DO L=L1-1,NU_old,-1
       if(DLYR_3D(3,L,NY,NX).LE.ZERO2)CYCLE

       FX=1.0_r8;FY=0._r8
       !move one layer down
       call UpdateLayerMaterials(L2+1,L,L2+1,NY,NX,FX,FY,MoveFlags)
       
       DLYR_3D(3,L2+1,NY,NX)          = DLYR_3D(3,L2+1,NY,NX)+DLYR_3D(3,L,NY,NX)
       DLYR_3D(3,L,NY,NX)            = 0._r8
       CumDepz2LayBottom_vr(L,NY,NX) = CumDepz2LayBottom_vr(L-1,NY,NX)
       CumSoilThickness_vr(L,NY,NX)  = CumSoilThickness_vr(L-1,NY,NX)
       L2=L2-1
    ENDDO
    NU_col(NY,NX)=NU_old+1
    CumDepz2LayBottom_vr(0,NY,NX)=CumDepz2LayBottom_vr(NU_col(NY,NX),NY,NX)-DLYR_3D(3,NU_col(NY,NX),NY,NX)
  ENDIF

  if(NU_old<NU_col(NY,NX))then
    !find nonzero drip for zero mass layer
    DO L=NU_old,1,-1
      DO ids=ids_beg,ids_end
        if(trcs_solml_drib_vr(ids,L,NY,NX)>0._r8 .and. DLYR_3D(3,L,NY,NX)<ZERO2)then
          trcs_solml_drib_vr(ids,NU_col(NY,NX),NY,NX)=trcs_solml_drib_vr(ids,NU_col(NY,NX),NY,NX)+trcs_solml_drib_vr(ids,L,NY,NX)
          trcs_solml_drib_vr(ids,L,NY,NX)=0._r8        
        endif
      ENDDO
    ENDDO
  endif

  call PrintInfo('end '//subname)
  end subroutine RemoveNullPondLayer
!------------------------------------------------------------------------------------------
  subroutine DisplaceCurrLayerEdges(I,J,LX,NY,NX,DVLiceMicP_vr,DORGC_vr,DDLYTop_vr,DDLYBot_vr,IFLGL)

  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: LX,NY,NX               !current grid coordinates, Z, Y, X
  REAL(R8), INTENT(IN) :: DVLiceMicP_vr(JZ)      !change in ice volume, final - initial
  real(r8), intent(in) :: DORGC_vr(JZ)           !change in SOM, final - initial 
  real(r8), intent(inout) :: DDLYTop_vr(0:JZ,6)       !layer lower edge displacement
  real(r8), intent(inout) :: DDLYBot_vr(0:JZ,6)       !layer upper edge displacement
  integer , intent(inout) :: IFLGL(0:JZ,6)       !flag of change type
  character(len=*), parameter :: subname='DisplaceCurrLayerEdges'

    !     POND, from water to soil
    ! a layer is made up of soil micropores+macropores and solid materials (can be zero).
  call PrintInfo('beg '//subname)

  ! Current grid is pond water 
  IF(IsPondLayer(SoilBulkDensity_vr(LX,NY,NX)))THEN  
    !IFLGL can be set to isl_soil or isl_water
    call ComputePondEdgeChange(I,J,LX,NY,NX,DDLYTop_vr(:,ich_watlev),DDLYBot_vr(:,ich_watlev),IFLGL(:,ich_watlev))

    DDLYTop_vr(LX,ich_frzthaw) = DDLYTop_vr(LX+1,ich_frzthaw)
    DDLYBot_vr(LX,ich_frzthaw) = DDLYBot_vr(LX+1,ich_frzthaw)
    DDLYTop_vr(LX,ich_erosion) = DDLYTop_vr(LX+1,ich_erosion)
    DDLYBot_vr(LX,ich_erosion) = DDLYBot_vr(LX+1,ich_erosion)
    DDLYTop_vr(LX,ich_sombgc)  = DDLYTop_vr(LX+1,ich_sombgc)
    DDLYBot_vr(LX,ich_sombgc)  = DDLYBot_vr(LX+1,ich_sombgc)

    IF(LX.EQ.NU_col(NY,NX))then
      !set thickness change to litter layer
      DDLYTop_vr(LX-1,ich_frzthaw) = DDLYTop_vr(LX,ich_frzthaw)
      DDLYBot_vr(LX-1,ich_frzthaw) = DDLYBot_vr(LX,ich_frzthaw)
      DDLYTop_vr(LX-1,ich_erosion) = DDLYTop_vr(LX,ich_erosion)
      DDLYBot_vr(LX-1,ich_erosion) = DDLYBot_vr(LX,ich_erosion)
      DDLYTop_vr(LX-1,ich_sombgc)  = DDLYTop_vr(LX,ich_sombgc)
      DDLYBot_vr(LX-1,ich_sombgc)  = DDLYBot_vr(LX,ich_sombgc)
    endif
    !==================================================
    !Current layer is SOIL
    !==================================================
  ELSE
    !
    IFLGL(LX,ich_frzthaw) = isl_undef
    call FreezeThawSoilEdgeChange(I,J,LX,NY,NX,DVLiceMicP_vr,DDLYTop_vr(:,ich_frzthaw),DDLYBot_vr(:,ich_frzthaw),IFLGL(:,ich_frzthaw))

    IFLGL(LX,ich_erosion) = isl_undef !can be corrected into isl_water
    call ErosionSoilEdgeChange(I,J,LX,NY,NX,DDLYTop_vr(:,ich_erosion),DDLYBot_vr(:,ich_erosion),IFLGL(:,ich_erosion))

    IFLGL(LX,ich_sombgc)  = isl_undef !can be corrected into isl_soil or isl_water
    call SOMBGCSoilEdgeChange(I,J,LX,NY,NX,DORGC_vr,DDLYTop_vr(:,ich_sombgc),DDLYBot_vr(:,ich_sombgc),IFLGL(:,ich_sombgc))

  ENDIF
  call PrintInfo('end '//subname)
  end subroutine DisplaceCurrLayerEdges
!------------------------------------------------------------------------------------------

  subroutine SOMBGCSoilEdgeChange(I,J,LX,NY,NX,DORGC_vr,DDLYTop_vr,DDLYBot_vr,IFLGL)
  ! Soil layer thickness is changed due to addition/removal of organic matter
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: LX,NY,NX  
  real(r8), intent(in) :: DORGC_vr(JZ)              !change in SOM, final -initial 
  real(r8), intent(inout) :: DDLYTop_vr(0:JZ)       !displacement of layer upper edge
  real(r8), intent(inout) :: DDLYBot_vr(0:JZ)       !displacement of layer lower edge
  integer , intent(inout) :: IFLGL(0:JZ)       !flag of change type
  character(len=*), parameter :: subname='SOMBGCSoilEdgeChange'
  real(r8) :: DDLEqv_OrgC   !soil thickness added due to change in organic matter

  ! SOC GAIN OR LOSS, i.e. DORGC_vr(LX) significantly non-zero
  ! keeping macropore fraction
  ! SoiBulkDensityt0_vr: initial bulk density,

  IF((iErosionMode.EQ.ieros_frzthawsom .OR. iErosionMode.EQ.ieros_frzthawsomeros) &
    .AND. ABS(DORGC_vr(LX)).GT.ZEROS(NY,NX))THEN

    DDLEqv_OrgC = gC2MgOM*DORGC_vr(LX)/((1.0_r8-SoilFracAsMacP_vr(LX,NY,NX))*SoiBulkDensityt0_vr(LX,NY,NX)*AREA_3D(3,LX,NY,NX))
    !use the new definition, dVol/area
!    DDLEqv_OrgC = gC2MgOM*DORGC_vr(LX)/(SoiBulkDensityt0_vr(LX,NY,NX))/AREA_3D(3,LX,NY,NX)

    ! LX is bottom layer, or is litter layer
    IF(LX.EQ.NL_col(NY,NX) .OR. IsPondLayer(SoilBulkDensity_vr(LX+1,NY,NX)))THEN 
      DDLYTop_vr(LX) = DDLEqv_OrgC   !lower edge
      DDLYBot_vr(LX) = 0.0_r8        !upper edge
      IFLGL(LX)      = isl_water       !underneath layer is water
      !not bottom layer or litter layer
    ELSE          
      DDLYBot_vr(LX) = DDLYTop_vr(LX+1)
      DDLYTop_vr(LX) = DDLYBot_vr(LX)+DDLEqv_OrgC      
      IFLGL(LX)      = isl_soil

      ! top layer, the layer right above is water
      IF(LX.EQ.NU_col(NY,NX) .OR. IsPondLayer(SoilBulkDensity_vr(LX-1,NY,NX)))THEN
        DDLYTop_vr(LX-1) = DDLYTop_vr(LX)
        DDLYBot_vr(LX-1) = DDLYTop_vr(LX)
        if(IsPondLayer(SoilBulkDensity_vr(LX-1,NY,NX)))then
          IFLGL(LX-1)    = isl_water     !layer LX is water
        else
          IFLGL(LX-1) = isl_soil
        endif
      ENDIF
    ENDIF
    !erosion model if off or current layer som change cause no change  
  ELSE
    ! bottom layer
    IF(LX.EQ.NL_col(NY,NX))THEN
      DDLYTop_vr(LX) = 0.0_r8
      DDLYBot_vr(LX) = 0.0_r8
      !non-bottom layer  
    ELSE
      DDLYTop_vr(LX) = DDLYTop_vr(LX+1)
      DDLYBot_vr(LX) = DDLYTop_vr(LX+1)
    ENDIF
  ENDIF
  end subroutine SOMBGCSoilEdgeChange
!------------------------------------------------------------------------------------------  

  subroutine ErosionSoilEdgeChange(I,J,LX,NY,NX,DDLYTop_vr,DDLYBot_vr,IFLGL)
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: LX,NY,NX
  real(r8), intent(inout) :: DDLYTop_vr(0:JZ)       !displacement of layer upper edge
  real(r8), intent(inout) :: DDLYBot_vr(0:JZ)       !displacement of layer lower edge
  integer , intent(inout) :: IFLGL(0:JZ)       !flag of change type
  character(len=*), parameter :: subname='ErosionSoilEdgeChange'

  real(r8) :: DDLEqv_Erosion

  call PrintInfo('beg '//subname)
  !========================================================================================
  !     EROSION model is on
  !========================================================================================
  IF((iErosionMode.EQ.ieros_frzthaweros .OR. iErosionMode.EQ.ieros_frzthawsomeros) &
    .AND. ABS(tErosionSedmLoss_col(NY,NX)).GT.ZEROS(NY,NX))THEN
    !  5:  layer thickness change due to sediment erosion

    !Bottom layer
    IF(LX.EQ.NL_col(NY,NX))THEN
      ! total soil layer reduction due to erosion
      IF(SoilBulkDensity_vr(LX,NY,NX).GT.ZERO)THEN
        DDLEqv_Erosion = -tErosionSedmLoss_col(NY,NX)/(SoilMicPMassLayerMX(NY,NX)/VLSoilPoreMicP_vr(NU_col(NY,NX),NY,NX))
      ELSE
        DDLEqv_Erosion = 0._r8
      ENDIF
      DDLYTop_vr(LX)   = DDLEqv_Erosion     !move down
      DDLYBot_vr(LX)   = DDLEqv_Erosion     !move down, this makes sure that sediment thickness does not change
      IFLGL(LX)        = isl_water
      !do nothing to layers other than bottom  
    ELSE
      DDLYTop_vr(LX) = 0.0_r8
      DDLYBot_vr(LX) = 0.0_r8
    ENDIF
  ELSE
    DDLYTop_vr(LX) = 0.0_r8
    DDLYBot_vr(LX) = 0.0_r8
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine ErosionSoilEdgeChange

!------------------------------------------------------------------------------------------
  subroutine FreezeThawSoilEdgeChange(I,J,LX,NY,NX,DVLiceMicP_vr,DDLYTop_vr,DDLYBot_vr,IFLGL)
  !
  !Layer LX is soil
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: LX,NY,NX
  REAL(R8), INTENT(IN) :: DVLiceMicP_vr(JZ)       !change in ice volume, final - initial  
  real(r8), intent(inout) :: DDLYTop_vr(0:JZ)       !displacement of layer upper edge
  real(r8), intent(inout) :: DDLYBot_vr(0:JZ)       !displacement of layer lower edge
  integer , intent(inout) :: IFLGL(0:JZ)       !flag of change type

  character(len=*), parameter :: subname='FreezeThawSoilEdgeChange'
  real(r8) :: DENSJ
  real(r8) :: DDLWatEqv_IceMicP !thickness change due to freeze-thaw, >0 freeze, < 0 thaw

  call PrintInfo('beg '//subname)

  !     FREEZE-THAW
  !DVLiceMicP: change in ice volume, initial-final
  !DENSICE: mass density of ice, g/cm3
  !DENSJ:  void volume added per unit addition of ice (with respect to ice)
  !the change is put to layer thickness
  !DDLWatEqv_IceMicP: added thickness change

  ! 4: there is volume change due to freeze-thaw
  IF(ABS(DVLiceMicP_vr(LX)).GT.ZEROS(NY,NX))THEN
    DENSJ  = 1._r8-DENSICE

    !thickness change due to freeze-thaw
    DDLWatEqv_IceMicP = DVLiceMicP_vr(LX)*DENSJ/(AREA_3D(3,LX,NY,NX))

    !Layer LX is at the bottom
    IF(LX.EQ.NL_col(NY,NX))THEN
      !save lower edge displacement for layer LX
      DDLYTop_vr(LX) = DDLWatEqv_IceMicP       !expand/thrink at the layer bottom
      DDLYBot_vr(LX) = 0.0_r8                  !layer bottom is set to zero

      !inner layer or top layer
    ELSE
      !obtain the upper edge displacement for layer LX as the bottom edge displacement of layer LX+1
      DDLYBot_vr(LX) = DDLYTop_vr(LX+1)
      DDLYTop_vr(LX) = DDLYBot_vr(LX)+DDLWatEqv_IceMicP !lower edge displacement for layer LX

      !layer LX is top layer or it is the first soil layer beneath water
      IF(LX.EQ.NU_col(NY,NX) .OR. IsPondLayer(SoilBulkDensity_vr(LX-1,NY,NX)))THEN
        !initialize displacements of layer LX-1
        DDLYBot_vr(LX-1) = DDLYTop_vr(LX)
        DDLYTop_vr(LX-1) = DDLYTop_vr(LX)         
        IFLGL(LX-1)    = isl_undef
      ENDIF
    ENDIF
    ! no freezing-induced change in ice volume 
  ELSE        
    DDLWatEqv_IceMicP=0.0_r8
    IF(LX.EQ.NL_col(NY,NX))THEN
      DDLYTop_vr(LX) = 0.0_r8
      DDLYBot_vr(LX) = 0.0_r8
    ELSE
      DDLYTop_vr(LX) = DDLYTop_vr(LX+1)
      DDLYBot_vr(LX) = DDLYTop_vr(LX+1)
      
      IF(LX.EQ.NU_col(NY,NX))THEN
        DDLYTop_vr(LX-1) = DDLYTop_vr(LX)
        DDLYBot_vr(LX-1) = DDLYTop_vr(LX)
        IFLGL(LX-1)    = isl_undef
      ENDIF
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine FreezeThawSoilEdgeChange

!------------------------------------------------------------------------------------------
  subroutine ComputePondEdgeChange(I,J,LX,NY,NX,DDLYTop_vr,DDLYBot_vr,IFLGL)
  !
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: LX,NY,NX               !grid coordinates Z,Y,X
  real(r8), intent(inout) :: DDLYTop_vr(0:JZ)    !displacement of bottom edge of layer LX
  real(r8), intent(inout) :: DDLYBot_vr(0:JZ)    !displacement of top edge of layer LX
  integer , intent(inout) :: IFLGL(0:JZ)         !flag of next layer type, water or soil

  character(len=*), parameter :: subname='ComputePondEdgeChange'
  real(r8) :: DLYR_XMicP            !layer thickness change
  real(r8) :: DLYR_next

  call PrintInfo('beg '//subname)

  ! Current layer LX is supposed to be water
  ! Next layer is pond
  IF(IsPondLayer(SoilBulkDensity_vr(LX+1,NY,NX)))THEN   
    
    !DLYRI: initial water layer thickness, [m]
    !Current layer: Compute the excessive layer thickness minus the initial layer thickness  

    !DLYR_XMicP: layer thickness change
    DLYR_XMicP=DLYRI_3D(3,LX,NY,NX)-(VLWatMicP_vr(LX,NY,NX)+VLiceMicP_vr(LX,NY,NX))/AREA_3D(3,LX,NY,NX)

    !Next layer thickness DLYR_next: water+ice total thickness
    DLYR_next=(VLWatMicP_vr(LX+1,NY,NX)+VLiceMicP_vr(LX+1,NY,NX))/AREA_3D(3,LX,NY,NX)

    !Current layer is expanding (due to expansion or addition), 
    !Or current layer is non-shrinking, but the layer below has meaningful water thickness

    !DLYR_XMicP> 0, shrinking, DLYR_XMicP < 0, expanding, compared to initial layer thickness        
    !significant layer change
    IF(DLYR_XMicP.LT.-ZERO .OR. DLYR_next.GT.ZERO)THEN      !layer LX is expanding or layer LX+1 has significant water mass
      DDLYTop_vr(LX) = DDLYTop_vr(LX+1)+DLYR_XMicP          !Combine the excessive thickness to layer below
      DDLYBot_vr(LX) = AMIN1(DDLYTop_vr(LX+1),DLYR_next)    !lower displacement as the minimum, this limits the water disappearance rate

      !Layer below has significant thickness
      IF(DLYR_next.GT.ZERO)THEN
        IFLGL(LX)=isl_water
        !Layer below has no significant thickness
      ELSE
        IFLGL(LX)=isl_soil
      ENDIF

      !no significant change compared to initial thickness
    ELSE   
      !DLYR_XMicP > 0, shrinking, DLYR_XMicP < 0 expanding      
      DLYR_XMicP     = DLYR_3D(3,LX,NY,NX)-(VLWatMicP_vr(LX,NY,NX)+VLiceMicP_vr(LX,NY,NX))/AREA_3D(3,LX,NY,NX)      
      DDLYBot_vr(LX) = DDLYTop_vr(LX+1)               !copy the thickness of layer below        
      DDLYTop_vr(LX) = DDLYBot_vr(LX)+DLYR_XMicP      !Combine the excessive thickness to layer below
      IFLGL(LX)      = isl_soil
    ENDIF
    !next layer is soil
  ELSE
    !Compute layer thickness change, DLYR_3D thickness from previous time step
    DLYR_XMicP     = DLYR_3D(3,LX,NY,NX)-(VLWatMicP_vr(LX,NY,NX)+VLiceMicP_vr(LX,NY,NX))/AREA_3D(3,LX,NY,NX)
    IFLGL(LX)      = isl_soil                          !Mark the change type,  pond deepens
    DDLYBot_vr(LX) = DDLYTop_vr(LX+1)                  !update the top edge
    DDLYTop_vr(LX) = DDLYBot_vr(LX)+DLYR_XMicP         !accumulate the lower edge displacement, <0 (expansion)
  ENDIF

  !Reach the top layer or the layer above is soil or litter (it floats on water)
  IF(LX.EQ.NU_col(NY,NX))THEN              !current is top layer, layer above is litter
    DDLYTop_vr(LX-1) = DDLYTop_vr(LX)      !copy the current layer thickness to layer above
    DDLYBot_vr(LX-1) = DDLYTop_vr(LX)      !
    IFLGL(LX-1)      = isl_water           !layer LX-1 no change in thickness
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine ComputePondEdgeChange

!------------------------------------------------------------------------------------------
  subroutine UpdateLayerEdges(I,J,NY,NX,DORGC_vr,DVLiceMicP_vr,SoilLayBotEdgeLocOld_vr,SoilLayBotEdgeLocNew_vr,IFLGL)
! Description
! readjust layer thickness
! IFLGL: c1, ponding water, c2, pond disappear, c3, pond reappear, c4: freeze-thaw, c5: erosion, c6: som change
  implicit none
  integer,  intent(in) :: I,J
  integer,  intent(in) :: NY,NX                             !column location
  REAL(R8), INTENT(IN) :: DVLiceMicP_vr(JZ)                 !change in ice volume, final - initial
  real(r8), intent(in) :: DORGC_vr(JZ)                      !change in SOM, initial - final
  real(r8), intent(out) :: SoilLayBotEdgeLocOld_vr(JZ)      !location of layer lower edge before thickness update
  real(r8), intent(out) :: SoilLayBotEdgeLocNew_vr(0:JZ)    !cumulative depth after thickness update
  integer , intent(out) :: IFLGL(0:JZ,6)                    !flag of next layer type for different changes
  real(r8) :: DDLYTop_vr(0:JZ,6)                            !displacement of layer upper edge
  real(r8) :: DDLYBot_vr(0:JZ,6)                            !displacement of layer lower edge
  integer  :: LX,LY,LL,NN,L
  real(r8) :: DLYR
  character(len=*), parameter :: subname='UpdateLayerEdges'

! begin_execution
  call PrintInfo('beg '//subname)

  !Column surface is water layer, note SoilBulkDensity_vr<=0 is set at the model initialization
  !determine topy layer to decide if a column is pond or soil

  DDLYTop_vr      = 0._r8
  DDLYBot_vr      = 0._r8
  IFLGL(0:JZ,:) = isl_undef

  SoilLayBotEdgeLocNew_vr = 0._r8
  SoilLayBotEdgeLocOld_vr = 0._r8

  !loop grid by bottom-up
  D225: DO LX=NL_col(NY,NX),NU_col(NY,NX),-1
    !make a copy of the depth, bottom of the layer
    SoilLayBotEdgeLocOld_vr(LX)  = CumDepz2LayBottom_vr(LX,NY,NX)
    SoilLayBotEdgeLocNew_vr(LX)  = CumDepz2LayBottom_vr(LX,NY,NX)

    !==================================================
    !Compute edge displacement
    !==================================================
    call DisplaceCurrLayerEdges(I,J,LX,NY,NX,DVLiceMicP_vr,DORGC_vr,DDLYTop_vr,DDLYBot_vr,IFLGL)

    !========================================================================================
    !  Apply the change and RESET SOIL LAYER DEPTHS
    !========================================================================================
    D200: DO NN=1,6

      ! update edge location due to process NN
      IF(NN.NE.iTopLayShrink .AND. NN.NE.iTopLayGrow)THEN !not for LX=NU_col on pond grow/shrink

        IF(IsPondLayer(SoilBulkDensity_vr(LX,NY,NX)))THEN
          !========================================================================================
          !     POND
          !========================================================================================
          ! there are some change induced by process NN
          !next layer can be water for soil; when flag is undefined, the change is due to sediment or freeze-thaw
          !when flag is defined, depth is changed due to water (dis)appear; OM change also works here
          IF(IFLGL(LX,NN).NE.isl_undef)THEN 
            !move lower edge of layer LX
            CumDepz2LayBottom_vr(LX,NY,NX) = CumDepz2LayBottom_vr(LX,NY,NX)+DDLYBot_vr(LX,NN)
            SoilLayBotEdgeLocNew_vr(LX)    = SoilLayBotEdgeLocNew_vr(LX)+DDLYBot_vr(LX,NN)

            !layer LX is not top layer, and the next layer is soil; the column is pond
            IF(LX.NE.NU_col(NY,NX) .AND. IFLGL(LX,ich_watlev).EQ.isl_soil)THEN   !lowest water layer
              !move all layers above           
              DO LL=LX-1,0,-1
                CumDepz2LayBottom_vr(LL,NY,NX) = CumDepz2LayBottom_vr(LL,NY,NX)+DDLYBot_vr(LX,NN)
                SoilLayBotEdgeLocNew_vr(LL)    = SoilLayBotEdgeLocNew_vr(LL)+DDLYBot_vr(LX,NN)
              ENDDO
              !set top edge displacement to zero, this avoids double counting of displacement for shallower layers
              DDLYTop_vr(LX,NN)=0.0_r8
            ENDIF

            !LX is top layer
            IF(LX.EQ.NU_col(NY,NX))THEN
              !Thickness of layer LX
              DLYR                             = (VLWatMicP_vr(LX,NY,NX)+VLiceMicP_vr(LX,NY,NX))/AREA_3D(3,LX,NY,NX)
              CumDepz2LayBottom_vr(LX-1,NY,NX) = CumDepz2LayBottom_vr(LX,NY,NX)-DLYR
              SoilLayBotEdgeLocNew_vr(LX-1)    = SoilLayBotEdgeLocNew_vr(LX)-DLYR
            ENDIF
          ENDIF
        ELSE
          !========================================================================================
          !  SOIL
          !========================================================================================

          
          !================================================
          !     FREEZE-THAW process
          !================================================
          IF(NN.EQ.ich_frzthaw)THEN
            !update lower edge of layer LX
            CumDepz2LayBottom_vr(LX,NY,NX) = CumDepz2LayBottom_vr(LX,NY,NX)+DDLYBot_vr(LX,NN)
            SoilLayBotEdgeLocNew_vr(LX)    = SoilLayBotEdgeLocNew_vr(LX)+DDLYBot_vr(LX,NN)

            !LX is top layer
            IF(LX.EQ.NU_col(NY,NX))THEN
              !update lower edge of layer LX-1
              CumDepz2LayBottom_vr(LX-1,NY,NX) = CumDepz2LayBottom_vr(LX-1,NY,NX)+DDLYBot_vr(LX-1,NN)
              SoilLayBotEdgeLocNew_vr(LX-1)    = SoilLayBotEdgeLocNew_vr(LX-1)+DDLYBot_vr(LX-1,NN)
            ENDIF
          ENDIF

          !================================================
          ! SET SURFACE ELEVATION FOR SOIL EROSION
          !================================================
          IF(NN.EQ.ich_erosion .AND. IFLGL(LX,NN).EQ.isl_water)THEN
            CumDepz2LayBottom_vr(LX,NY,NX) = CumDepz2LayBottom_vr(LX,NY,NX)+DDLYBot_vr(LX,NN)
            SoilLayBotEdgeLocNew_vr(LX)    = SoilLayBotEdgeLocNew_vr(LX)+DDLYBot_vr(LX,NN)

            IF(LX.EQ.NU_col(NY,NX))THEN
              CumDepz2LayBottom_vr(LX-1,NY,NX) = CumDepz2LayBottom_vr(LX-1,NY,NX)+DDLYBot_vr(LX,NN)
              SoilLayBotEdgeLocNew_vr(LX-1)    = SoilLayBotEdgeLocNew_vr(LX-1)+DDLYBot_vr(LX,NN)
            ENDIF
          ENDIF

          !================================================
          ! SET SOIL LAYER DEPTHS FOR CHANGES IN SOC
          !================================================
          IF(NN.EQ.ich_sombgc .AND. IFLGL(LX,NN).NE.isl_undef)THEN
            CumDepz2LayBottom_vr(LX,NY,NX) = CumDepz2LayBottom_vr(LX,NY,NX)+DDLYBot_vr(LX,NN)
            SoilLayBotEdgeLocNew_vr(LX)    = SoilLayBotEdgeLocNew_vr(LX)+DDLYBot_vr(LX,NN)

            !Top layer or layer above is ponding water
            IF(LX.EQ.NU_col(NY,NX) .OR. IsPondLayer(SoilBulkDensity_vr(LX-1,NY,NX)))THEN
              CumDepz2LayBottom_vr(LX-1,NY,NX) = CumDepz2LayBottom_vr(LX-1,NY,NX)+DDLYBot_vr(LX-1,NN)
              SoilLayBotEdgeLocNew_vr(LX-1)    = SoilLayBotEdgeLocNew_vr(LX-1)+DDLYBot_vr(LX-1,NN)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO D200
    VLSoilMicP_vr(LX,NY,NX)=VLSoilPoreMicP_vr(LX,NY,NX)
  ENDDO D225
  VLSoilMicP_vr(0,NY,NX)=VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)

  call PrintInfo('end '//subname)
  end subroutine UpdateLayerEdges

!------------------------------------------------------------------------------------------

  subroutine ComputeLayerAllocFs(I,J,L,NN,NY,NX,NUX,DDLYRX,IFLGL,FX,FO,L_donor,L_rcept)

  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NN            !type of layer change
  integer, intent(in) :: NUX           !old top layer index
  integer, intent(in) :: L,NY,NX       !
  real(r8),intent(in) :: DDLYRX(3)     !layer thickness change
  integer, intent(in) :: IFLGL(0:JZ,6)
  real(r8), intent(out) :: FX        !fraction of soil or pond contents redistributed from L0 to L1
  real(r8), intent(out) :: FO        !allocation coefficient for organic material
  integer, intent(out) :: L_donor   !source
  integer, intent(out) :: L_rcept   !target
  character(len=*), parameter :: subname='ComputeLayerAllocFs'

  real(r8) :: DLEqv_MicP

! begin_execution
  call PrintInfo('beg '//subname)

  IF(DDLYRX(NN).GT.ZERO)THEN !moving down
    !
    !     LAYERS MOVE DOWN (FX>0.), pond drys up
    !
    IF(IFLGL(L,iTopLayShrink).EQ.isl_undef)THEN  !layer is sthrinking
      L_rcept = L
      L_donor = L+1
    ELSE
      L_rcept = NU_col(NY,NX)
      L_donor = NUX
    ENDIF

    IF((SoilBulkDensity_vr(L,NY,NX).LE.ZERO .AND. IFLGL(L,ich_watlev).EQ.isl_soil)      & !layer L is water and layer below is soil
        .OR. (DLYR_3D(3,L_donor,NY,NX).LE.ZEROC .AND. IFLGL(L,ich_sombgc).EQ.isl_water))THEN   !layer L0 disappear and layer L below is water
      FX = 1.0_r8
      FO = 1.0_r8
    ELSE
      IF(SoilBulkDensity_vr(L_donor,NY,NX).LE.ZERO)THEN !layer L0 is water
        !Thickness of layer L0
        DLEqv_MicP=(VLWatMicP_vr(L_donor,NY,NX)+VLiceMicP_vr(L_donor,NY,NX))/AREA_3D(3,L_donor,NY,NX)

        IF(DLEqv_MicP.GT.ZERO)THEN
          FX = AMIN1(1.0_r8,DDLYRX(NN)/DLEqv_MicP)
          FO = FX
        ELSE
          FX = 0.0_r8
          FO = 0.0_r8
        ENDIF
      ELSE
        IF(DLYR_3D(3,L,NY,NX).GT.ZERO .AND. DLYR_3D(3,L_donor,NY,NX).GT.ZERO)THEN
          FX = AMIN1(1.0,DDLYRX(NN)/DLYR_3D(3,L,NY,NX))
          FO = AMIN1(1.0,DDLYRX(NN)/DLYR_3D(3,L_donor,NY,NX))
        ELSE
          FX = 0.0_r8
          FO = 0.0_r8
        ENDIF
      ENDIF
    ENDIF
  ELSE
    !
    !     LAYERS MOVE UP (FX<=0.), pond builds up
    !
    IF(IFLGL(L,iTopLayGrow).EQ.isl_undef)THEN !
      L_rcept = L+1
      L_donor = L
    ELSE
      L_rcept =NU_col(NY,NX)
      L_donor = 0
    ENDIF

    !Layer L0 is pond water
    IF(SoilBulkDensity_vr(L_donor,NY,NX).LE.ZERO)THEN
      !Thickness of layer L0
      DLEqv_MicP=(VLWatMicP_vr(L_donor,NY,NX)+VLiceMicP_vr(L_donor,NY,NX))/AREA_3D(3,L_donor,NY,NX)
      IF(DLEqv_MicP.GT.ZERO)THEN
        FX = AMIN1(1.0,-DDLYRX(NN)/DLEqv_MicP)
        FO = FX
      ELSE
        FX = 0.0_r8
        FO = 0.0_r8
      ENDIF
    ELSE
      !
      IF(DLYR_3D(3,L,NY,NX).GT.ZERO .AND. DLYR_3D(3,L_donor,NY,NX).GT.ZERO)THEN
        FX = AMIN1(1.0,-DDLYRX(NN)/DLYR_3D(3,L,NY,NX))
        FO = AMIN1(1.0,-DDLYRX(NN)/DLYR_3D(3,L_donor,NY,NX))
      ELSE
        FX = 0.0_r8
        FO = 0.0_r8
      ENDIF
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine ComputeLayerAllocFs
!------------------------------------------------------------------------------------------

  subroutine DiagInnerLayerThickChange(I,J,L,NY,NX,SoilLayBotEdgeLocOld_vr,IFLGLWat,ICHKL,NUX,SoilLayBotEdgeLocNew_vr,DDLYRXNN,DDLYRY,XVOLWP)
  !
  !update layer L thickness, up-> down
  !
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: L                 !current layer index
  integer, intent(in) :: NY,NX             !current horizontal location
  real(r8), intent(in) :: SoilLayBotEdgeLocOld_vr(JZ)       !layer edge location before update 
  integer , intent(in) :: IFLGLWat(0:JZ)                !edge change flag due to different processes
  real(r8), intent(inout) :: SoilLayBotEdgeLocNew_vr(0:JZ)   !layer edge location after update
  integer,  intent(inout) :: ICHKL
  integer , intent(out) :: NUX                            !old top layer index
  real(r8), intent(out)   :: DDLYRXNN                    !thickness change of layer L due to all kinds of disturbances, >0 thrink/move down, <0 expand/move up
  real(r8), intent(inout) :: DDLYRY(JZ)                   !change in layer thickness required to restore the initial layer thickness
  real(r8), intent(out)   :: XVOLWP                       !volume of surface mobile water [m3/d2]
  character(len=*), parameter :: subname='DiagInnerLayerThickChange'

  real(r8) :: DLYR0
  real(r8) :: DLYRXX
  
  integer :: LL, loc

! begin_execution
  call PrintInfo('beg '//subname)

  NUX=0
  
  !layer L thickness
  DLYR_3D(3,L,NY,NX) = CumDepz2LayBottom_vr(L,NY,NX)-CumDepz2LayBottom_vr(L-1,NY,NX) 
  !make a copy of layer thickness
  DLYRXX             = DLYR_3D(3,L,NY,NX)                                           

  !layer L has no change, but L+1 has change
  IF(IFLGLWat(L).EQ.isl_undef .AND. IFLGLWat(L+1).NE.isl_undef)THEN
    DDLYRXNN = 0.0_r8
    loc      = 1
    IF(IsPondLayer(SoilBulkDensity_vr(L,NY,NX)))THEN    !layer L is ponding water
      DDLYRY(L)=DLYRI_3D(3,L,NY,NX)-DLYR_3D(3,L,NY,NX)  !required layer thickness for restoration
    ELSE
      DDLYRY(L)=0.0_r8                                  !no change for soil layer
    ENDIF
    !layer checked
    ICHKL=1

    !layer L has change, and (layer L+1 has no change, or layer L shrinked)  
  ELSEIF(IFLGLWat(L).EQ.isl_soil .AND. (IFLGLWat(L+1).EQ.isl_undef .OR. DLYR_3D(3,L,NY,NX).LE.DLYRI_3D(3,L,NY,NX)))THEN
    DDLYRXNN=0.0_r8
    loc=2
    !layer L is at top or layer L is not updated
    IF(L.EQ.NU_col(NY,NX) .OR. ICHKL.EQ.0)THEN
      DDLYRY(L)=0.0_r8
    ELSE
      DDLYRY(L)=DDLYRY(L-1)
    ENDIF

    IF(IFLGLWat(L).EQ.isl_soil .AND. IFLGLWat(L+1).EQ.isl_undef)ICHKL=0

    !layer L has change, or layer L+1 has no change
  ELSE
    loc=3
    IF(ICHKL.EQ.0)THEN
      DDLYRXNN = DLYRI_3D(3,L,NY,NX)-DLYR_3D(3,L,NY,NX)
      DDLYRY(L)  = DDLYRXNN
    ELSE
      DDLYRXNN = 0.0_r8
      DDLYRY(L)  = DDLYRY(L-1)
    ENDIF
  ENDIF

  !restore layer thickness
  CumDepz2LayBottom_vr(L,NY,NX) = CumDepz2LayBottom_vr(L,NY,NX)+DDLYRY(L)
  DLYR_3D(3,L,NY,NX)            = DLYR_3D(3,L,NY,NX)+DDLYRY(L)
  SoilDepthMidLay_vr(L,NY,NX)   = 0.5_r8*(CumDepz2LayBottom_vr(L,NY,NX)+CumDepz2LayBottom_vr(L-1,NY,NX))
  CumSoilThickness_vr(L,NY,NX)  = CumDepz2LayBottom_vr(L,NY,NX)-CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)

  !layer L is rigth above the bottom boundary
  IF(L.EQ.NL_col(NY,NX)-1)THEN
    DLYR_3D(3,L+1,NY,NX)           = CumDepz2LayBottom_vr(L+1,NY,NX)-CumDepz2LayBottom_vr(L,NY,NX)
    SoilDepthMidLay_vr(L+1,NY,NX)  = 0.5_r8*(CumDepz2LayBottom_vr(L+1,NY,NX)+CumDepz2LayBottom_vr(L,NY,NX))
    CumSoilThickness_vr(L+1,NY,NX) = CumDepz2LayBottom_vr(L+1,NY,NX)-CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)
  ENDIF

  !layer L is at the top boundary 
  IF(L.EQ.NU_col(NY,NX))THEN
    CumSoilThickMidL_vr(L,NY,NX)=0.5_r8*CumSoilThickness_vr(L,NY,NX)  !layer middle thickness
  ELSE
    CumSoilThickMidL_vr(L,NY,NX)=0.5_r8*(CumSoilThickness_vr(L,NY,NX)+CumSoilThickness_vr(L-1,NY,NX))
  ENDIF

  !current layer is soil 
  IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO)THEN
    DDLYRXNN=SoilLayBotEdgeLocNew_vr(L)-SoilLayBotEdgeLocOld_vr(L)
  ENDIF
      
  call PrintInfo('end '//subname)
  end subroutine DiagInnerLayerThickChange
!------------------------------------------------------------------------------------------
  subroutine DiagPondTopThickChange(I,J,L,NN,NY,NX,NUX,SoilLayBotEdgeLocNew_vr,IFLGL,DDLYRX,XVOLWP)
  ! Description:
  ! Diagnose pond surface recession and activate an adjacent layer on regrowth.
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NN
  integer, intent(in) :: L                 !current layer index
  integer, intent(in) :: NY,NX             !current horizontal location
  integer, intent(out):: NUX  
  real(r8), intent(inout) :: SoilLayBotEdgeLocNew_vr(0:JZ)   !layer edge location after update
  integer , intent(inout) :: IFLGL(0:JZ,6)                !edge change flag due to different processes
  real(r8), intent(out)   :: DDLYRX(3)                    !thickness change of layer L due to all kinds of disturbances, >0 thrink/move down, <0 expand/move up
  real(r8), intent(out)   :: XVOLWP                       !volume of surface mobile water [m3/d2]
  real(r8) :: DLYR0
  integer  :: LL
  character(len=*), parameter :: subname='DiagPondTopThickChange'

  call PrintInfo('beg '//subname)

  IF(L.NE.NU_col(NY,NX))then
    !return if it not top layer
    DDLYRX(NN)  = 0.0_r8
    IFLGL(L,NN) = isl_undef
    return
  endif    

  !below is for L==NU_col(NY,NX)
  !====================================================================================
  !     RESET POND SURFACE LAYER NUMBER IF LOST TO EVAPORATION
  !====================================================================================  
  IF(NN.EQ.iTopLayShrink)THEN !deal with top pond layer thrinking
    !layer L is lost
    IF(IsPondLayer(SoilBulkDensity_vr(NU_col(NY,NX),NY,NX))               & !current layer L=NU is pond water
      .AND. (VHeatCapacity_vr(NU_col(NY,NX),NY,NX).LE.VHCPNX_col(NY,NX)   & !current layer is dried out
      .OR. NUM_col(NY,NX).GT.NU_col(NY,NX)))THEN                            !current layer is moving down/dried out, note that 
                                                                            !NUM_col is copied from NU_col and updated in WatsubMod.F90,
                                                                            !and NU_col is to be updated
      NUX=NU_col(NY,NX)    
      DO LL=NUX+1,NL_col(NY,NX)
        IF(VLSoilPoreMicP_vr(LL,NY,NX).GT.ZEROS2(NY,NX))THEN !identify the first nonzero-thickness layer
          NU_col(NY,NX)         = LL
          DDLYRX(iTopLayShrink) = DLYR_3D(3,NUX,NY,NX)        !whole layer is gone > 0, moving down          
          IFLGL(L,iTopLayShrink) = isl_water                !next layer is water  
          
          DLYR_3D(3,NUX,NY,NX) = 0.0_r8      !it is dried out, so layer thickness is zero

          !layer NUX is in pond
          IF(IsPondLayer(SoilBulkDensity_vr(NUX,NY,NX)))THEN
            VGeomLayer_vr(NUX,NY,NX)     = AREA_3D(3,NUX,NY,NX)*DLYR_3D(3,NUX,NY,NX)
            VLSoilPoreMicP_vr(NUX,NY,NX) = VGeomLayer_vr(NUX,NY,NX)
          ENDIF
          exit
        ENDIF
      ENDDO
    ENDIF

    !====================================================================================
    !     RESET POND SURFACE LAYER NUMBER IF GAIN FROM PRECIPITATION
    !====================================================================================
  ELSEIF(NN.EQ.iTopLayGrow)THEN !work on top layer
    !CumLitRDepzInit_col is the initial litter layer bottom

    !obtain the free water volume defined as water beyond the water holding capacity of litter layer
    XVOLWP=AZMAX1(VLWatMicP_vr(0,NY,NX)-VLWatHeldCapSurf_col(NY,NX))
    IF(CumDepz2LayBottom_vr(0,NY,NX).GT.CumLitRDepzInit_col(NY,NX)              & !layer NU is lower than its intial depth, i.e. sinked
      .AND. XVOLWP.GT.(VLWatHeldCapSurf_col(NY,NX)+VHCPNX_col(NY,NX)/cpw))THEN    !sufficient free water for ponding

      !layer L is soil and current top layer is below the initial top layer
      IF(SoilBulkDensity_vr(L,NY,NX).GT.ZERO .AND. NU_col(NY,NX).GT.NUI_col(NY,NX))THEN
        
        ! Reactivate the adjacent retired pond slot. Jumping to NUI leaves
        ! zero-thickness layers between the new pond and the soil below.
        NU_col(NY,NX)  = L-1
        NUM_col(NY,NX) = NU_col(NY,NX)

        !assuming litter floats above water, compute the lower edge displacement
        DDLYRX(iTopLayGrow)  = -XVOLWP/AREA_3D(3,0,NY,NX)      !ponding water moving up <0
        IFLGL(L,iTopLayGrow) = isl_water

        !obtain litter layer + pond water thickness
        DLYR0 = (AZMAX1(VLWatMicP_vr(0,NY,NX)+VLiceMicP_vr(0,NY,NX)-VWatLitRHoldCapcity_col(NY,NX))+VLitR_col(NY,NX))/AREA_3D(3,0,NY,NX)

        !update litter layer thickness:=total thickness - free ponding water
        DLYR_3D(3,0,NY,NX) = DLYR0+DDLYRX(iTopLayGrow)

        !update the topsoil/water layer: add water
        DLYR_3D(3,NU_col(NY,NX),NY,NX) = DLYR_3D(3,NU_col(NY,NX),NY,NX)-DDLYRX(iTopLayGrow)

        !reset edges for inner layers 
        IF(L.GT.2)THEN    !L==NU_col      
          !set to upper edge of layer L. L varies from, for innter layers 
          DO LL=L-2,NU_col(NY,NX),-1
            CumDepz2LayBottom_vr(LL,NY,NX) = CumDepz2LayBottom_vr(L-1,NY,NX)
            SoilLayBotEdgeLocNew_vr(LL)    = SoilLayBotEdgeLocNew_vr(L-1)
          ENDDO
        ENDIF

        CumDepz2LayBottom_vr(0,NY,NX) = CumDepz2LayBottom_vr(NU_col(NY,NX),NY,NX)-DLYR_3D(3,NU_col(NY,NX),NY,NX)
        SoilLayBotEdgeLocNew_vr(0)    = SoilLayBotEdgeLocNew_vr(NU_col(NY,NX))-DLYR_3D(3,NU_col(NY,NX),NY,NX)

        SoilDepthMidLay_vr(NU_col(NY,NX),NY,NX)  = 0.5_r8*(CumDepz2LayBottom_vr(NU_col(NY,NX),NY,NX)+CumDepz2LayBottom_vr(0,NY,NX))
        CumSoilThickness_vr(NU_col(NY,NX),NY,NX) = DLYR_3D(3,NU_col(NY,NX),NY,NX)
        CumSoilThickMidL_vr(NU_col(NY,NX),NY,NX) = 0.5_r8*CumSoilThickness_vr(NU_col(NY,NX),NY,NX)
      ENDIF
    ENDIF
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine DiagPondTopThickChange
!------------------------------------------------------------------------------------------

  subroutine UpdateLayerMaterials(L,L0,L1,NY,NX,FX,FY,IFLGL)
  ! Description:
  ! Move material inventories from source layer L0 to receiving layer L1.
  implicit none
  integer, intent(in) :: L
  integer, intent(in) :: L0  !source layer
  integer, intent(in) :: L1  !target layer
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: FX,FY
  integer, intent(in) ::  IFLGL(0:JZ,6)
  integer :: N,M,NZ,K,NGL,NR,NE,NTU,NTSA,NTSAB,idg,NTP
  integer :: NTX,NTF,MID,idom,ids
  real(r8) :: ENGY0,ENGY1
  character(len=*), parameter :: subname='UpdateLayerMaterials'

! begin_execution
  call PrintInfo('beg '//subname)
  ! Receiving layer: pools present only in soil (no surface-litter pool).
  IF(L0.NE.0)THEN
    SAND_vr(L1,NY,NX)                  = SAND_vr(L1,NY,NX)+FX*SAND_vr(L0,NY,NX)
    SILT_vr(L1,NY,NX)                  = SILT_vr(L1,NY,NX)+FX*SILT_vr(L0,NY,NX)
    CLAY_vr(L1,NY,NX)                  = CLAY_vr(L1,NY,NX)+FX*CLAY_vr(L0,NY,NX)
    trcx_solml_vr(idx_CEC,L1,NY,NX) = trcx_solml_vr(idx_CEC,L1,NY,NX)+FX*trcx_solml_vr(idx_CEC,L0,NY,NX)
    trcx_solml_vr(idx_AEC,L1,NY,NX) = trcx_solml_vr(idx_AEC,L1,NY,NX)+FX*trcx_solml_vr(idx_AEC,L0,NY,NX)
    VLWatMacP_vr(L1,NY,NX) = VLWatMacP_vr(L1,NY,NX) +FX*VLWatMacP_vr(L0,NY,NX)
    VLiceMacP_vr(L1,NY,NX) = VLiceMacP_vr(L1,NY,NX) +FX*VLiceMacP_vr(L0,NY,NX)

    ! Banded fertilizer is allocated only for soil layers 1:JZ.
    DO NTF=ifertnb_beg,ifertnb_end
      FertN_mole_Band_vr(NTF,L1,NY,NX)=FertN_mole_Band_vr(NTF,L1,NY,NX)+FX*FertN_mole_Band_vr(NTF,L0,NY,NX)
    ENDDO
    FertP_mole_band_vr(L1,NY,NX)=FertP_mole_band_vr(L1,NY,NX)+FX*FertP_mole_band_vr(L0,NY,NX)
    trcs_solml_vr(ids_H1PO4B,L1,NY,NX)=trcs_solml_vr(ids_H1PO4B,L1,NY,NX)+FX*trcs_solml_vr(ids_H1PO4B,L0,NY,NX)
    trcs_solml_vr(ids_H2PO4B,L1,NY,NX)=trcs_solml_vr(ids_H2PO4B,L1,NY,NX)+FX*trcs_solml_vr(ids_H2PO4B,L0,NY,NX)

    IF(salt_model)THEN
      DO NTSAB=idsaltb_beg,idsaltb_end
        trcSalt_solml_vr(NTSAB,L1,NY,NX)=trcSalt_solml_vr(NTSAB,L1,NY,NX)+FX*trcSalt_solml_vr(NTSAB,L0,NY,NX)
      ENDDO
      DO NTSA=idsalt_beg,idsaltb_end
        trcSalt_soHml_vr(NTSA,L1,NY,NX)=trcSalt_soHml_vr(NTSA,L1,NY,NX)+FX*trcSalt_soHml_vr(NTSA,L0,NY,NX)
      ENDDO
    ENDIF

    DO NTX=idx_beg+1,idx_cation_end
      trcx_solml_vr(NTX,L1,NY,NX)=trcx_solml_vr(NTX,L1,NY,NX)+FX*trcx_solml_vr(NTX,L0,NY,NX)
    ENDDO

    DO NTX=idx_AEC+1,idx_end
      trcx_solml_vr(NTX,L1,NY,NX)=trcx_solml_vr(NTX,L1,NY,NX)+FX*trcx_solml_vr(NTX,L0,NY,NX)
    ENDDO

    DO NTP=idsp_beg,idsp_end
      trcp_saltpml_vr(NTP,L1,NY,NX)=trcp_saltpml_vr(NTP,L1,NY,NX)+FX*trcp_saltpml_vr(NTP,L0,NY,NX)
    ENDDO

    DO idg=idg_beg,idg_NH3
      trcg_gasml_vr(idg,L1,NY,NX)=trcg_gasml_vr(idg,L1,NY,NX)+FX*trcg_gasml_vr(idg,L0,NY,NX)
      trcg_root_vr(idg,L1,NY,NX)=trcg_root_vr(idg,L1,NY,NX)+FX*trcg_root_vr(idg,L0,NY,NX)
    ENDDO
    do ids=ids_beg,ids_end
      trcs_soHml_vr(ids,L1,NY,NX)=trcs_soHml_vr(ids,L1,NY,NX)+FX*trcs_soHml_vr(ids,L0,NY,NX)
    enddo

    ! Root and mycorrhizal pools.
    DO  NZ=1,NP0_col(NY,NX)
      ! Material inventories move independently of source/receiver activity.

      DO  NR=1,NumStructuralRootAxes_pft(NZ,NY,NX)
        DO NE=1,NumPlantChemElms
          RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ,NY,NX)=RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ,NY,NX) &
            +FX*RootMyco1stStrutElms_rpvr(NE,L0,NR,NZ,NY,NX)
        ENDDO
        Root1stLenPP_rpvr(L1,NR,NZ,NY,NX)  = Root1stLenPP_rpvr(L1,NR,NZ,NY,NX)+FX*Root1stLenPP_rpvr(L0,NR,NZ,NY,NX)
      ENDDO
      Root1stXNumL_pvr(L1,NZ,NY,NX)        = Root1stXNumL_pvr(L1,NZ,NY,NX)+FX*Root1stXNumL_pvr(L0,NZ,NY,NX)

      DO N=1,Myco_pft(NZ,NY,NX)
        DO idg=idg_beg,idg_NH3
          trcg_rootml_pvr(idg,N,L1,NZ,NY,NX) = trcg_rootml_pvr(idg,N,L1,NZ,NY,NX)+FX*trcg_rootml_pvr(idg,N,L0,NZ,NY,NX)
          trcs_rootml_pvr(idg,N,L1,NZ,NY,NX) = trcs_rootml_pvr(idg,N,L1,NZ,NY,NX)+FX*trcs_rootml_pvr(idg,N,L0,NZ,NY,NX)
        ENDDO

        DO  NR=1,NumStructuralRootAxes_pft(NZ,NY,NX)
          DO NE=1,NumPlantChemElms
            RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX)=RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX) &
              +FX*RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)
          ENDDO
          Root2ndLen_rpvr(N,L1,NR,NZ,NY,NX)   = Root2ndLen_rpvr(N,L1,NR,NZ,NY,NX)+FX*Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)
          Root2ndXNum_rpvr(N,L1,NR,NZ,NY,NX) = Root2ndXNum_rpvr(N,L1,NR,NZ,NY,NX)+FX*Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX)
        ENDDO

        DO NE=1,NumPlantChemElms
           RootMycoNonstElms_rpvr(NE,N,L1,NZ,NY,NX)=RootMycoNonstElms_rpvr(NE,N,L1,NZ,NY,NX)+FX* RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)
        ENDDO
        RootMycoActiveBiomC_pvr(N,L1,NZ,NY,NX) = RootMycoActiveBiomC_pvr(N,L1,NZ,NY,NX)+FX*RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX)
        PopuRootMycoC_pvr(N,L1,NZ,NY,NX)       = PopuRootMycoC_pvr(N,L1,NZ,NY,NX)+FX* PopuRootMycoC_pvr(N,L0,NZ,NY,NX)
        RootProteinC_pvr(N,L1,NZ,NY,NX)        = RootProteinC_pvr(N,L1,NZ,NY,NX)+FX*RootProteinC_pvr(N,L0,NZ,NY,NX)
        Root2ndXNumL_rpvr(N,L1,NZ,NY,NX)         = Root2ndXNumL_rpvr(N,L1,NZ,NY,NX)+FX*Root2ndXNumL_rpvr(N,L0,NZ,NY,NX)
        RootTotLenPerPlant_pvr(N,L1,NZ,NY,NX)     = RootTotLenPerPlant_pvr(N,L1,NZ,NY,NX)+FX*RootTotLenPerPlant_pvr(N,L0,NZ,NY,NX)
        RootLenDensPerPlant_pvr(N,L1,NZ,NY,NX) = RootLenDensPerPlant_pvr(N,L1,NZ,NY,NX)+FX*RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX)
        RootPoreVol_pvr(N,L1,NZ,NY,NX)         = RootPoreVol_pvr(N,L1,NZ,NY,NX)+FX*RootPoreVol_pvr(N,L0,NZ,NY,NX)
        RootVH2O_pvr(N,L1,NZ,NY,NX)            = RootVH2O_pvr(N,L1,NZ,NY,NX)+FX*RootVH2O_pvr(N,L0,NZ,NY,NX)
        Root1stRadius_pvr(N,L1,NZ,NY,NX)       = Root1stRadius_pvr(N,L1,NZ,NY,NX)+FX*Root1stRadius_pvr(N,L0,NZ,NY,NX)
        Root2ndRadius_rpvr(N,L1,NZ,NY,NX)       = Root2ndRadius_rpvr(N,L1,NZ,NY,NX)+FX*Root2ndRadius_rpvr(N,L0,NZ,NY,NX)
        RootSAreaPerPlant_pvr(N,L1,NZ,NY,NX)    = RootSAreaPerPlant_pvr(N,L1,NZ,NY,NX)+FX*RootSAreaPerPlant_pvr(N,L0,NZ,NY,NX)
        Root2ndEffLen4uptk_rpvr(N,L1,NZ,NY,NX)       = Root2ndEffLen4uptk_rpvr(N,L1,NZ,NY,NX)+FX*Root2ndEffLen4uptk_rpvr(N,L0,NZ,NY,NX)
      ENDDO
      DO NE=1,NumPlantChemElms
        RootNodulStrutElms_rpvr(NE,L1,NZ,NY,NX) = RootNodulStrutElms_rpvr(NE,L1,NZ,NY,NX)+FX*RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)
        RootNodulNonstElms_rpvr(NE,L1,NZ,NY,NX) = RootNodulNonstElms_rpvr(NE,L1,NZ,NY,NX)+FX*RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)
      ENDDO
    ENDDO
  ENDIF

  ! Receiving layer: pools also present in surface litter (L0 = 0).

  VLWatMicP_vr(L1,NY,NX)          = VLWatMicP_vr(L1,NY,NX) +FX*VLWatMicP_vr(L0,NY,NX)
  VLiceMicP_vr(L1,NY,NX)          = VLiceMicP_vr(L1,NY,NX) +FX*VLiceMicP_vr(L0,NY,NX)
  VLsoiAirP_vr(L1,NY,NX)          = VLsoiAirP_vr(L1,NY,NX) +FX*VLsoiAirP_vr(L0,NY,NX)
  VLMicP_vr(L1,NY,NX)             = VLMicP_vr(L1,NY,NX)    +FX*VLMicP_vr(L0,NY,NX)
  VLSoilMicP_vr(L1,NY,NX)         = VLSoilMicP_vr(L1,NY,NX)+FX*VLSoilMicP_vr(L0,NY,NX)
  VLWatMicPX_vr(L1,NY,NX)         = VLWatMicP_vr(L1,NY,NX)

  ! Reconstruct temperature after both macropore and micropore water transfers.
  ENGY1                           = VHeatCapacity_vr(L1,NY,NX)*TKS_vr(L1,NY,NX)
  ENGY0                           = VHeatCapacity_vr(L0,NY,NX)*TKS_vr(L0,NY,NX)
  ENGY1                           = ENGY1+FX*ENGY0
  VHeatCapSolidSoil_vr(L1,NY,NX) = VHeatCapSolidSoil_vr(L1,NY,NX)+FX*VHeatCapSolidSoil_vr(L0,NY,NX)
  VHeatCapacity_vr(L1,NY,NX)      = VHeatCapSolidSoil_vr(L1,NY,NX) &
    +cpw*(VLWatMicP_vr(L1,NY,NX)+VLWatMacP_vr(L1,NY,NX)) &
    +cpi*(VLiceMicP_vr(L1,NY,NX)+VLiceMacP_vr(L1,NY,NX))

  IF(VHeatCapacity_vr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS_vr(L1,NY,NX)=ENGY1/VHeatCapacity_vr(L1,NY,NX)
  ELSE
    TKS_vr(L1,NY,NX)=TKS_vr(L0,NY,NX)
  ENDIF
  TCS_vr(L1,NY,NX)=units%Kelvin2Celcius(TKS_vr(L1,NY,NX))

  DO NTF=ifertn_beg,ifertn_end
    FertN_mole_soil_vr(NTF,L1,NY,NX)=FertN_mole_soil_vr(NTF,L1,NY,NX)+FX*FertN_mole_soil_vr(NTF,L0,NY,NX)
  ENDDO
  FertP_mole_soil_vr(L1,NY,NX)=FertP_mole_soil_vr(L1,NY,NX) + FX*FertP_mole_soil_vr(L0,NY,NX)

  DO NTU=ids_nuts_beg,ids_nuts_end
    if(NTU/=ids_H2PO4B .and. NTU/=ids_H1PO4B)THEN
      trcs_solml_vr(NTU,L1,NY,NX)=trcs_solml_vr(NTU,L1,NY,NX)+FX*trcs_solml_vr(NTU,L0,NY,NX)
    ENDIF
  ENDDO

  IF(salt_model)THEN
    DO NTSA=idsalt_beg,idsalt_end
      trcSalt_solml_vr(NTSA,L1,NY,NX)=trcSalt_solml_vr(NTSA,L1,NY,NX)+FX*trcSalt_solml_vr(NTSA,L0,NY,NX)
    ENDDO
  ENDIF

  !Exclude NH3 and NH3B, which are accounted in nutrients
  DO idg=idg_beg,idg_end-2
    trcs_solml_vr(idg,L1,NY,NX)=trcs_solml_vr(idg,L1,NY,NX)+FX*trcs_solml_vr(idg,L0,NY,NX)
  ENDDO

  IF(IFLGL(L,iTopLayGrow).EQ.isl_undef)THEN
    DO  K=1,jcplx
      DO  N=1,NumMicbHFunGrupsPerCmplx
        DO  M=1,nlbiomcp
          DO NGL=JGniH(N),JGnfH(N)
            MID=micpar%get_micb_id(M,NGL)
            DO NE=1,NumPlantChemElms
              mBiomeHeter_vr(NE,MID,K,L1,NY,NX)=mBiomeHeter_vr(NE,MID,K,L1,NY,NX)+FX*mBiomeHeter_vr(NE,MID,K,L0,NY,NX)
            ENDDO
          enddo
        enddo
      enddo
    ENDDO
    DO  N=1,NumMicbAFunGrupsPerCmplx
      DO  M=1,nlbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeAutor_vr(NE,MID,L1,NY,NX)=mBiomeAutor_vr(NE,MID,L1,NY,NX)+FX*mBiomeAutor_vr(NE,MID,L0,NY,NX)
          ENDDO
        enddo
      enddo
    enddo

    DO K=1,jcplx
      DO  M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          OMBioResdu_vr(NE,M,K,L1,NY,NX)=OMBioResdu_vr(NE,M,K,L1,NY,NX)+FX*OMBioResdu_vr(NE,M,K,L0,NY,NX)
        ENDDO
      ENDDO
      DO idom=idom_beg,idom_end
        DOM_MicP_vr(idom,K,L1,NY,NX)      = DOM_MicP_vr(idom,K,L1,NY,NX)+FX*DOM_MicP_vr(idom,K,L0,NY,NX)
        DOM_MacP_vr(idom,K,L1,NY,NX) = DOM_MacP_vr(idom,K,L1,NY,NX)+FX*DOM_MacP_vr(idom,K,L0,NY,NX)
        SorbedOM_vr(idom,K,L1,NY,NX) = SorbedOM_vr(idom,K,L1,NY,NX)+FX*SorbedOM_vr(idom,K,L0,NY,NX)
      enddo

      DO M=1,jsken
        SolidOMAct_vr(M,K,L1,NY,NX)=SolidOMAct_vr(M,K,L1,NY,NX)+FX*SolidOMAct_vr(M,K,L0,NY,NX)
        DO NE=1,NumPlantChemElms
          SolidOM_vr(NE,M,K,L1,NY,NX)=SolidOM_vr(NE,M,K,L1,NY,NX)+FX*SolidOM_vr(NE,M,K,L0,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  DO ids=ids_beg,ids_end
    trcs_solml_drib_vr(ids,L1,NY,NX)=trcs_solml_drib_vr(ids,L1,NY,NX)+FX*trcs_solml_drib_vr(ids,L0,NY,NX)
  ENDDO
  !===================================================
  ! SOURCE LAYER
  !===================================================
  ! Source layer: pools present only in soil (no surface-litter pool).
  IF(L0.NE.0)THEN
    SAND_vr(L0,NY,NX)                  = FY*SAND_vr(L0,NY,NX)
    SILT_vr(L0,NY,NX)                  = FY*SILT_vr(L0,NY,NX)
    CLAY_vr(L0,NY,NX)                  = FY*CLAY_vr(L0,NY,NX)
    trcx_solml_vr(idx_CEC,L0,NY,NX) = FY*trcx_solml_vr(idx_CEC,L0,NY,NX)
    trcx_solml_vr(idx_AEC,L0,NY,NX) = FY*trcx_solml_vr(idx_AEC,L0,NY,NX)
    VLWatMacP_vr(L0,NY,NX) = FY*VLWatMacP_vr(L0,NY,NX)
    VLiceMacP_vr(L0,NY,NX) = FY*VLiceMacP_vr(L0,NY,NX)

    ! Banded fertilizer is allocated only for soil layers 1:JZ.
    DO NTF=ifertnb_beg,ifertnb_end
      FertN_mole_Band_vr(NTF,L0,NY,NX)=FY*FertN_mole_Band_vr(NTF,L0,NY,NX)
    ENDDO
    FertP_mole_band_vr(L0,NY,NX)=FY*FertP_mole_band_vr(L0,NY,NX)
    trcs_solml_vr(ids_H1PO4B,L0,NY,NX) = FY*trcs_solml_vr(ids_H1PO4B,L0,NY,NX)
    trcs_solml_vr(ids_H2PO4B,L0,NY,NX) = FY*trcs_solml_vr(ids_H2PO4B,L0,NY,NX)
    IF(salt_model)THEN
      DO NTSAB=idsaltb_beg,idsaltb_end
        trcSalt_solml_vr(NTSAB,L0,NY,NX)=FY*trcSalt_solml_vr(NTSAB,L0,NY,NX)
      ENDDO
      DO NTSA=idsalt_beg,idsaltb_end
        trcSalt_soHml_vr(NTSA,L0,NY,NX)=FY*trcSalt_soHml_vr(NTSA,L0,NY,NX)
      ENDDO
    ENDIF

    DO NTX=idx_beg+1,idx_cation_end
      trcx_solml_vr(NTX,L0,NY,NX)=FY*trcx_solml_vr(NTX,L0,NY,NX)
    ENDDO
    DO NTX=idx_AEC+1,idx_end
      trcx_solml_vr(NTX,L0,NY,NX)=FY*trcx_solml_vr(NTX,L0,NY,NX)
    ENDDO

    DO NTP=idsp_beg,idsp_end
      trcp_saltpml_vr(NTP,L0,NY,NX)=FY*trcp_saltpml_vr(NTP,L0,NY,NX)
    ENDDO

    DO idg=idg_beg,idg_NH3
      trcg_gasml_vr(idg,L0,NY,NX)=FY*trcg_gasml_vr(idg,L0,NY,NX)
      trcg_root_vr(idg,L0,NY,NX)=FY*trcg_root_vr(idg,L0,NY,NX)
    ENDDO
    do ids=ids_beg,ids_end
      trcs_soHml_vr(ids,L0,NY,NX)=FY*trcs_soHml_vr(ids,L0,NY,NX)
    enddo

    ! Root and mycorrhizal pools.
    DO  NZ=1,NP0_col(NY,NX)
      ! Material inventories move independently of source/receiver activity.
      DO NR=1,NumStructuralRootAxes_pft(NZ,NY,NX)
        DO NE=1,NumPlantChemElms
          RootMyco1stStrutElms_rpvr(NE,L0,NR,NZ,NY,NX) = FY*RootMyco1stStrutElms_rpvr(NE,L0,NR,NZ,NY,NX)
        ENDDO
        Root1stLenPP_rpvr(L0,NR,NZ,NY,NX)  = FY*Root1stLenPP_rpvr(L0,NR,NZ,NY,NX)
      ENDDO
      Root1stXNumL_pvr(L0,NZ,NY,NX)        = FY*Root1stXNumL_pvr(L0,NZ,NY,NX)
      DO  N=1,Myco_pft(NZ,NY,NX)
        DO idg=idg_beg,idg_NH3
          trcg_rootml_pvr(idg,N,L0,NZ,NY,NX) = FY*trcg_rootml_pvr(idg,N,L0,NZ,NY,NX)
          trcs_rootml_pvr(idg,N,L0,NZ,NY,NX) = FY*trcs_rootml_pvr(idg,N,L0,NZ,NY,NX)
        ENDDO
        DO NR=1,NumStructuralRootAxes_pft(NZ,NY,NX)
          DO NE=1,NumPlantChemElms
            RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX) = FY*RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)
          ENDDO
          Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)  = FY*Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)
          Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX) = FY*Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX)
        ENDDO
        DO NE=1,NumPlantChemElms
          RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)=FY* RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)
        ENDDO
        RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX) = FY*RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX)
        PopuRootMycoC_pvr(N,L0,NZ,NY,NX)       = FY* PopuRootMycoC_pvr(N,L0,NZ,NY,NX)
        RootProteinC_pvr(N,L0,NZ,NY,NX)        = FY*RootProteinC_pvr(N,L0,NZ,NY,NX)
        Root2ndXNumL_rpvr(N,L0,NZ,NY,NX)       = FY*Root2ndXNumL_rpvr(N,L0,NZ,NY,NX)
        RootTotLenPerPlant_pvr(N,L0,NZ,NY,NX)  = FY*RootTotLenPerPlant_pvr(N,L0,NZ,NY,NX)
        RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX) = FY*RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX)
        RootPoreVol_pvr(N,L0,NZ,NY,NX)         = FY*RootPoreVol_pvr(N,L0,NZ,NY,NX)
        RootVH2O_pvr(N,L0,NZ,NY,NX)            = FY*RootVH2O_pvr(N,L0,NZ,NY,NX)
        Root1stRadius_pvr(N,L0,NZ,NY,NX)       = FY*Root1stRadius_pvr(N,L0,NZ,NY,NX)
        Root2ndRadius_rpvr(N,L0,NZ,NY,NX)      = FY*Root2ndRadius_rpvr(N,L0,NZ,NY,NX)
        RootSAreaPerPlant_pvr(N,L0,NZ,NY,NX)   = FY*RootSAreaPerPlant_pvr(N,L0,NZ,NY,NX)
        Root2ndEffLen4uptk_rpvr(N,L0,NZ,NY,NX) = FY*Root2ndEffLen4uptk_rpvr(N,L0,NZ,NY,NX)
      ENDDO
      DO NE=1,NumPlantChemElms
        RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)=FY*RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)

        RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)=FY*RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)
      ENDDO
    ENDDO
  ENDIF

  ! Source layer: pools also present in surface litter (L0 = 0).
  VLWatMicP_vr(L0,NY,NX)          = FY*VLWatMicP_vr(L0,NY,NX)
  VLiceMicP_vr(L0,NY,NX)          = FY*VLiceMicP_vr(L0,NY,NX)
  VLsoiAirP_vr(L0,NY,NX)          = FY*VLsoiAirP_vr(L0,NY,NX)
  VLMicP_vr(L0,NY,NX)             = FY*VLMicP_vr(L0,NY,NX)
  VLSoilMicP_vr(L0,NY,NX)         = FY*VLSoilMicP_vr(L0,NY,NX)
  VLWatMicPX_vr(L0,NY,NX)         = VLWatMicP_vr(L0,NY,NX)
  ENGY0                           = FY*ENGY0
  VHeatCapSolidSoil_vr(L0,NY,NX) = FY*VHeatCapSolidSoil_vr(L0,NY,NX)

  ! Soil includes macropore heat capacity; surface litter does not.
  IF(L0.NE.0)THEN
    VHeatCapacity_vr(L0,NY,NX)=VHeatCapSolidSoil_vr(L0,NY,NX) &
      +cpw*(VLWatMicP_vr(L0,NY,NX)+VLWatMacP_vr(L0,NY,NX)) &
      +cpi*(VLiceMicP_vr(L0,NY,NX)+VLiceMacP_vr(L0,NY,NX))
  ELSE
    VHeatCapacity_vr(L0,NY,NX)=VHeatCapSolidSoil_vr(L0,NY,NX)+cpw*VLWatMicP_vr(L0,NY,NX)+cpi*VLiceMicP_vr(L0,NY,NX)
  ENDIF
  IF(VHeatCapacity_vr(L0,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS_vr(L0,NY,NX)=ENGY0/VHeatCapacity_vr(L0,NY,NX)
  ELSE
    TKS_vr(L0,NY,NX)=TKS_vr(L1,NY,NX)
  ENDIF
  TCS_vr(L0,NY,NX)=units%Kelvin2Celcius(TKS_vr(L0,NY,NX))

  DO NTF=ifertn_beg,ifertn_end
    FertN_mole_soil_vr(NTF,L0,NY,NX)=FY*FertN_mole_soil_vr(NTF,L0,NY,NX)
  ENDDO
  FertP_mole_soil_vr(L0,NY,NX)=FY*FertP_mole_soil_vr(L0,NY,NX)

  DO NTU=ids_nuts_beg,ids_nuts_end
    if(NTU/=ids_H1PO4B .and. NTU/=ids_H2PO4B)THEN
      trcs_solml_vr(NTU,L0,NY,NX)=FY*trcs_solml_vr(NTU,L0,NY,NX)
    ENDIF
  ENDDO
  IF(salt_model)THEN
    DO NTSA=idsalt_beg,idsalt_end
      trcSalt_solml_vr(NTSA,L0,NY,NX)=FY*trcSalt_solml_vr(NTSA,L0,NY,NX)
    ENDDO
  ENDIF

  !exclude NH3 and NH3B, which are accounted in nutrients
  DO idg=idg_beg,idg_end-2
    trcs_solml_vr(idg,L0,NY,NX)=FY*trcs_solml_vr(idg,L0,NY,NX)
  ENDDO

  DO ids=ids_beg,ids_end
    trcs_solml_drib_vr(ids,L0,NY,NX)=FY*trcs_solml_drib_vr(ids,L0,NY,NX)
  ENDDO

  IF(IFLGL(L,iTopLayGrow).EQ.isl_undef)THEN
    DO  K=1,jcplx
       DO N=1,NumMicbHFunGrupsPerCmplx
        DO M=1,nlbiomcp
          DO NGL=JGniH(N),JGnfH(N)
            MID=micpar%get_micb_id(M,NGL)
            DO NE=1,NumPlantChemElms
              mBiomeHeter_vr(NE,MID,K,L0,NY,NX)=FY*mBiomeHeter_vr(NE,MID,K,L0,NY,NX)
            ENDDO
          ENDDO
        enddo
      enddo
    ENDDO

    DO N=1,NumMicbAFunGrupsPerCmplx
      DO M=1,nlbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            mBiomeAutor_vr(NE,MID,L0,NY,NX)=FY*mBiomeAutor_vr(NE,MID,L0,NY,NX)
          ENDDO
        ENDDO
      enddo
    enddo

    DO K=1,jcplx
      DO  M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          OMBioResdu_vr(NE,M,K,L0,NY,NX)=FY*OMBioResdu_vr(NE,M,K,L0,NY,NX)
        ENDDO
      ENDDO
      do idom=idom_beg,idom_end
        DOM_MicP_vr(idom,K,L0,NY,NX)      = FY*DOM_MicP_vr(idom,K,L0,NY,NX)
        DOM_MacP_vr(idom,K,L0,NY,NX) = FY*DOM_MacP_vr(idom,K,L0,NY,NX)
        SorbedOM_vr(idom,K,L0,NY,NX) = FY*SorbedOM_vr(idom,K,L0,NY,NX)
      enddo
      DO  M=1,jsken
        SolidOMAct_vr(M,K,L0,NY,NX)=FY*SolidOMAct_vr(M,K,L0,NY,NX)
        do NE=1,NumPlantChemElms
          SolidOM_vr(NE,M,K,L0,NY,NX)=FY*SolidOM_vr(NE,M,K,L0,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  call PrintInfo('end '//subname)
  end subroutine UpdateLayerMaterials

!------------------------------------------------------------------------------------------

  subroutine MoveSOM(L0,L1,L,NY,NX,FO,IFLGL)
!
  implicit none
  integer, intent(in) :: L0,L1,L
  integer, intent(in) :: NY,NX
  real(r8),intent(in) :: FO
  integer, intent(in) :: IFLGL(0:JZ,6)

  integer :: K,N,M,NGL,NR,NZ,NE,idg,MID,idom
  real(r8) :: FXO,FRO
  real(r8) :: FXRTLG2,FXRTN2,FXEPOOLR,FXWTRTL
  real(r8) :: WTNDLE,FXEPOOLN
  real(r8) :: FXWTRT1E
  real(r8) :: FXWTRT2E,FXRTLG1
  real(r8) :: FXWTRTD,FXWSRTL,FNumAxesPerPrimaryRoot_pft,FXRTNL,FXRTLGP,FXRTDNP
  real(r8) :: FXRTVLP,FXRTVLW,FXRRAD1,FXRRAD2,FXRootSAreaPerPlant_pvr,FXRTLGA
  real(r8) :: FXOQN,FXOQP,FXOQA,FXOQCH,FXOQNH,FXOQPH,FXOQAH
  real(r8) :: FXOHC,FXOHN,FXOHP,FXOHA,FXOSC,FXOSA,FXOSN,FXOSP
  real(r8) :: FXOMC,FXOMN,FXOMP,FXORC,FXORN,FXORP,FXOQC
  real(r8) :: FXGA,FXGP
  
! begin_execution
  IF(IFLGL(L,3).EQ.0 .AND. L0.NE.0 &
    .AND. VLSoilPoreMicP_vr(L0,NY,NX).GT.ZEROS(NY,NX) &
    .AND. VLSoilPoreMicP_vr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    IF(L0.EQ.L.OR.CORGCI_vr(L0,NY,NX).LE.ZERO)THEN
      FXO=FO
    ELSE
      FXO=AMIN1(0.5,FO*AMIN1(10.0,CORGCI_vr(L1,NY,NX)/CORGCI_vr(L0,NY,NX)))
    ENDIF

    DO  K=1,jcplx
      DO  N=1,NumMicbHFunGrupsPerCmplx
        DO  M=1,nlbiomcp
          DO NGL=JGniH(N),JGnfH(N)
            MID=micpar%get_micb_id(M,NGL)          
            DO NE=1,NumPlantChemElms
              FXOMC                             = FXO*mBiomeHeter_vr(NE,MID,K,L0,NY,NX)
              mBiomeHeter_vr(NE,MID,K,L1,NY,NX) = mBiomeHeter_vr(NE,MID,K,L1,NY,NX)+FXOMC
              mBiomeHeter_vr(NE,MID,K,L0,NY,NX) = mBiomeHeter_vr(NE,MID,K,L0,NY,NX)-FXOMC
            ENDDO
          enddo
        enddo
      enddo
    ENDDO

    DO  N=1,NumMicbAFunGrupsPerCmplx
      DO  M=1,nlbiomcp
        DO NGL=JGniA(N),JGnfA(N)
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            FXOMC                           = FXO*mBiomeAutor_vr(NE,MID,L0,NY,NX)
            mBiomeAutor_vr(NE,MID,L1,NY,NX) = mBiomeAutor_vr(NE,MID,L1,NY,NX)+FXOMC
            mBiomeAutor_vr(NE,MID,L0,NY,NX) = mBiomeAutor_vr(NE,MID,L0,NY,NX)-FXOMC
          ENDDO
        enddo
      enddo
    enddo

    DO  K=1,jcplx
      DO  M=1,ndbiomcp
        DO NE=1,NumPlantChemElms
          FXORC                          = FXO*OMBioResdu_vr(NE,M,K,L0,NY,NX)
          OMBioResdu_vr(NE,M,K,L1,NY,NX) = OMBioResdu_vr(NE,M,K,L1,NY,NX)+FXORC
          OMBioResdu_vr(NE,M,K,L0,NY,NX) = OMBioResdu_vr(NE,M,K,L0,NY,NX)-FXORC
        ENDDO
      ENDDO

      DO idom=idom_beg,idom_end
        FXOQC                   = FXO*DOM_MicP_vr(idom,K,L0,NY,NX)
        DOM_MicP_vr(idom,K,L1,NY,NX) = DOM_MicP_vr(idom,K,L1,NY,NX)+FXOQC
        DOM_MicP_vr(idom,K,L0,NY,NX) = DOM_MicP_vr(idom,K,L0,NY,NX)-FXOQC
      enddo

      IF(SoilFracAsMacP_vr(L1,NY,NX).GT.ZERO.AND.SoilFracAsMacP_vr(L0,NY,NX).GT.ZERO)THEN
        DO idom=idom_beg,idom_end
          FXOQCH                       = FXO*DOM_MacP_vr(idom,K,L0,NY,NX)
          DOM_MacP_vr(idom,K,L1,NY,NX) = DOM_MacP_vr(idom,K,L1,NY,NX)+FXOQCH
          DOM_MacP_vr(idom,K,L0,NY,NX) = DOM_MacP_vr(idom,K,L0,NY,NX)-FXOQCH
        ENDDO
      ENDIF
       DO idom=idom_beg,idom_end
        FXOHC                        = FXO*SorbedOM_vr(idom,K,L0,NY,NX)
        SorbedOM_vr(idom,K,L1,NY,NX) = SorbedOM_vr(idom,K,L1,NY,NX)+FXOHC
        SorbedOM_vr(idom,K,L0,NY,NX) = SorbedOM_vr(idom,K,L0,NY,NX)-FXOHC
      ENDDO

      DO M=1,jsken
        DO NE=1,NumPlantChemElms
          FXOSC                       = FXO*SolidOM_vr(NE,M,K,L0,NY,NX)
          SolidOM_vr(NE,M,K,L1,NY,NX) = SolidOM_vr(NE,M,K,L1,NY,NX)+FXOSC
          SolidOM_vr(NE,M,K,L0,NY,NX) = SolidOM_vr(NE,M,K,L0,NY,NX)-FXOSC
        ENDDO
      ENDDO
    ENDDO
!
!     ROOTS (why return?) 05/21/2022, jyt
!
    return
    DO NZ=1,NP_col(NY,NX)
      IF(RootMycoActiveBiomC_pvr(ipltroot,L0,NZ,NY,NX).GT.ZERO4Groth_pft(NZ,NY,NX) &
        .AND.RootMycoActiveBiomC_pvr(ipltroot,L1,NZ,NY,NX).GT.ZERO4Groth_pft(NZ,NY,NX))THEN
        IF(L0.EQ.L.OR.CumSoilThickMidL_vr(L1,NY,NX).LE.ZERO)THEN
          FRO=FO
        ELSE
          FRO=AMIN1(0.5,FO*CumSoilThickMidL_vr(L0,NY,NX)/CumSoilThickMidL_vr(L1,NY,NX))
        ENDIF

        DO  NR=1,NumStructuralRootAxes_pft(NZ,NY,NX)
          DO NE=1,NumPlantChemElms
            FXWTRT1E                                     = FRO*RootMyco1stStrutElms_rpvr(NE,L0,NR,NZ,NY,NX)
            RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ,NY,NX) = RootMyco1stStrutElms_rpvr(NE,L1,NR,NZ,NY,NX)+FXWTRT1E
            RootMyco1stStrutElms_rpvr(NE,L0,NR,NZ,NY,NX) = RootMyco1stStrutElms_rpvr(NE,L0,NR,NZ,NY,NX)-FXWTRT1E
          ENDDO  
          FXRTLG1                         = FRO*Root1stLenPP_rpvr(L0,NR,NZ,NY,NX)
          Root1stLenPP_rpvr(L1,NR,NZ,NY,NX) = Root1stLenPP_rpvr(L1,NR,NZ,NY,NX)+FXRTLG1
          Root1stLenPP_rpvr(L0,NR,NZ,NY,NX) = Root1stLenPP_rpvr(L0,NR,NZ,NY,NX)-FXRTLG1
        ENDDO

        FNumAxesPerPrimaryRoot_pft           = FRO*Root1stXNumL_pvr(L0,NZ,NY,NX)
        Root1stXNumL_pvr(L1,NZ,NY,NX) = Root1stXNumL_pvr(L1,NZ,NY,NX)+FNumAxesPerPrimaryRoot_pft
        Root1stXNumL_pvr(L0,NZ,NY,NX) = Root1stXNumL_pvr(L0,NZ,NY,NX)-FNumAxesPerPrimaryRoot_pft

        DO  N=1,Myco_pft(NZ,NY,NX)
          DO idg=idg_beg,idg_NH3
            FXGA                               = FRO*trcg_rootml_pvr(idg,N,L0,NZ,NY,NX)
            trcg_rootml_pvr(idg,N,L1,NZ,NY,NX) = trcg_rootml_pvr(idg,N,L1,NZ,NY,NX)+FXGA
            trcg_rootml_pvr(idg,N,L0,NZ,NY,NX) = trcg_rootml_pvr(idg,N,L0,NZ,NY,NX)-FXGA

            FXGP                               = FRO*trcs_rootml_pvr(idg,N,L0,NZ,NY,NX)
            trcs_rootml_pvr(idg,N,L1,NZ,NY,NX) = trcs_rootml_pvr(idg,N,L1,NZ,NY,NX)+FXGP
            trcs_rootml_pvr(idg,N,L0,NZ,NY,NX) = trcs_rootml_pvr(idg,N,L0,NZ,NY,NX)-FXGP
          ENDDO

          DO  NR=1,NumStructuralRootAxes_pft(NZ,NY,NX)
            DO NE=1,NumPlantChemElms
              FXWTRT2E                                       = FRO*RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)
              RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX) = RootMyco2ndStrutElms_rpvr(NE,N,L1,NR,NZ,NY,NX)+FXWTRT2E
              RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX) = RootMyco2ndStrutElms_rpvr(NE,N,L0,NR,NZ,NY,NX)-FXWTRT2E
            ENDDO

            FXRTLG2                          = FRO*Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)
            Root2ndLen_rpvr(N,L1,NR,NZ,NY,NX) = Root2ndLen_rpvr(N,L1,NR,NZ,NY,NX)+FXRTLG2
            Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX) = Root2ndLen_rpvr(N,L0,NR,NZ,NY,NX)-FXRTLG2

            FXRTN2                             = FRO*Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX)
            Root2ndXNum_rpvr(N,L1,NR,NZ,NY,NX) = Root2ndXNum_rpvr(N,L1,NR,NZ,NY,NX)+FXRTN2
            Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX) = Root2ndXNum_rpvr(N,L0,NR,NZ,NY,NX)-FXRTN2
          ENDDO
          DO NE=1,NumPlantChemElms
             FXEPOOLR                                 = FRO*RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)
             RootMycoNonstElms_rpvr(NE,N,L1,NZ,NY,NX) = RootMycoNonstElms_rpvr(NE,N,L1,NZ,NY,NX)+FXEPOOLR
             RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX) = RootMycoNonstElms_rpvr(NE,N,L0,NZ,NY,NX)-FXEPOOLR
          ENDDO

          FXWTRTL                                = FRO*RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX)
          RootMycoActiveBiomC_pvr(N,L1,NZ,NY,NX) = RootMycoActiveBiomC_pvr(N,L1,NZ,NY,NX)+FXWTRTL
          RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX) = RootMycoActiveBiomC_pvr(N,L0,NZ,NY,NX)-FXWTRTL

          FXWTRTD                          = FRO*PopuRootMycoC_pvr(N,L0,NZ,NY,NX)
          PopuRootMycoC_pvr(N,L1,NZ,NY,NX) = PopuRootMycoC_pvr(N,L1,NZ,NY,NX)+FXWTRTD
          PopuRootMycoC_pvr(N,L0,NZ,NY,NX) = PopuRootMycoC_pvr(N,L0,NZ,NY,NX)-FXWTRTD

          FXWSRTL                         = FRO*RootProteinC_pvr(N,L0,NZ,NY,NX)
          RootProteinC_pvr(N,L1,NZ,NY,NX) = RootProteinC_pvr(N,L1,NZ,NY,NX)+FXWSRTL
          RootProteinC_pvr(N,L0,NZ,NY,NX) = RootProteinC_pvr(N,L0,NZ,NY,NX)-FXWSRTL

          FXRTNL                         = FRO*Root2ndXNumL_rpvr(N,L0,NZ,NY,NX)
          Root2ndXNumL_rpvr(N,L1,NZ,NY,NX) = Root2ndXNumL_rpvr(N,L1,NZ,NY,NX)+FXRTNL
          Root2ndXNumL_rpvr(N,L0,NZ,NY,NX) = Root2ndXNumL_rpvr(N,L0,NZ,NY,NX)-FXRTNL

          FXRTLGP                            = FRO*RootTotLenPerPlant_pvr(N,L0,NZ,NY,NX)
          RootTotLenPerPlant_pvr(N,L1,NZ,NY,NX) = RootTotLenPerPlant_pvr(N,L1,NZ,NY,NX)+FXRTLGP
          RootTotLenPerPlant_pvr(N,L0,NZ,NY,NX) = RootTotLenPerPlant_pvr(N,L0,NZ,NY,NX)-FXRTLGP

          FXRTDNP                                = FRO*RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX)
          RootLenDensPerPlant_pvr(N,L1,NZ,NY,NX) = RootLenDensPerPlant_pvr(N,L1,NZ,NY,NX)+FXRTDNP
          RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX) = RootLenDensPerPlant_pvr(N,L0,NZ,NY,NX)-FXRTDNP

          FXRTVLP                        = FRO*RootPoreVol_pvr(N,L0,NZ,NY,NX)
          RootPoreVol_pvr(N,L1,NZ,NY,NX) = RootPoreVol_pvr(N,L1,NZ,NY,NX)+FXRTVLP
          RootPoreVol_pvr(N,L0,NZ,NY,NX) = RootPoreVol_pvr(N,L0,NZ,NY,NX)-FXRTVLP

          FXRTVLW                     = FRO*RootVH2O_pvr(N,L0,NZ,NY,NX)
          RootVH2O_pvr(N,L1,NZ,NY,NX) = RootVH2O_pvr(N,L1,NZ,NY,NX)+FXRTVLW
          RootVH2O_pvr(N,L0,NZ,NY,NX) = RootVH2O_pvr(N,L0,NZ,NY,NX)-FXRTVLW

          FXRRAD1                          = FRO*Root1stRadius_pvr(N,L0,NZ,NY,NX)
          Root1stRadius_pvr(N,L1,NZ,NY,NX) = Root1stRadius_pvr(N,L1,NZ,NY,NX)+FXRRAD1
          Root1stRadius_pvr(N,L0,NZ,NY,NX) = Root1stRadius_pvr(N,L0,NZ,NY,NX)-FXRRAD1

          FXRRAD2                          = FRO*Root2ndRadius_rpvr(N,L0,NZ,NY,NX)
          Root2ndRadius_rpvr(N,L1,NZ,NY,NX) = Root2ndRadius_rpvr(N,L1,NZ,NY,NX)+FXRRAD2
          Root2ndRadius_rpvr(N,L0,NZ,NY,NX) = Root2ndRadius_rpvr(N,L0,NZ,NY,NX)-FXRRAD2

          FXRootSAreaPerPlant_pvr              = FRO*RootSAreaPerPlant_pvr(N,L0,NZ,NY,NX)
          RootSAreaPerPlant_pvr(N,L1,NZ,NY,NX) = RootSAreaPerPlant_pvr(N,L1,NZ,NY,NX)+FXRootSAreaPerPlant_pvr
          RootSAreaPerPlant_pvr(N,L0,NZ,NY,NX) = RootSAreaPerPlant_pvr(N,L0,NZ,NY,NX)-FXRootSAreaPerPlant_pvr

          FXRTLGA                          = FRO*Root2ndEffLen4uptk_rpvr(N,L0,NZ,NY,NX)
          Root2ndEffLen4uptk_rpvr(N,L1,NZ,NY,NX) = Root2ndEffLen4uptk_rpvr(N,L1,NZ,NY,NX)+FXRTLGA
          Root2ndEffLen4uptk_rpvr(N,L0,NZ,NY,NX) = Root2ndEffLen4uptk_rpvr(N,L0,NZ,NY,NX)-FXRTLGA
        ENDDO
!
!     ROOT NODULES
!
        DO NE=1,NumPlantChemElms
          WTNDLE                                 = FRO*RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)
          RootNodulStrutElms_rpvr(NE,L1,NZ,NY,NX) = RootNodulStrutElms_rpvr(NE,L1,NZ,NY,NX)+WTNDLE
          RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX) = RootNodulStrutElms_rpvr(NE,L0,NZ,NY,NX)-WTNDLE

          FXEPOOLN                               = FRO*RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)
          RootNodulNonstElms_rpvr(NE,L1,NZ,NY,NX) = RootNodulNonstElms_rpvr(NE,L1,NZ,NY,NX)+FXEPOOLN
          RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX) = RootNodulNonstElms_rpvr(NE,L0,NZ,NY,NX)-FXEPOOLN

        ENDDO
      ENDIF
    ENDDO
  ENDIF
  end subroutine MoveSOM
!------------------------------------------------------------------------------------------

  subroutine MoveMacPoreSolute(L0,L1,NY,NX,FHO)

  implicit none
  integer, intent(in) :: L0,L1,NY,NX
  real(r8), intent(in) :: FHO
  integer :: NTS,NTSA
  real(r8) :: FXSH,FXH

! begin_execution
!
!     SOIL MACROPORE N,P SOLUTES
!
  IF(SoilFracAsMacP_vr(L1,NY,NX).GT.ZERO.AND.SoilFracAsMacP_vr(L0,NY,NX).GT.ZERO)THEN

    DO NTS=ids_beg,ids_end
      FXSH                       = FHO*trcs_soHml_vr(NTS,L0,NY,NX)
      trcs_soHml_vr(NTS,L1,NY,NX) = trcs_soHml_vr(NTS,L1,NY,NX)+FXSH
      trcs_soHml_vr(NTS,L0,NY,NX) = trcs_soHml_vr(NTS,L0,NY,NX)-FXSH
    ENDDO

!
!     SOIL MACROPORE SOLUBLE SALTS
!
    IF(salt_model)THEN
      DO NTSA=idsalt_beg,idsaltb_end
        FXH                             = FHO*trcSalt_soHml_vr(NTSA,L0,NY,NX)
        trcSalt_soHml_vr(NTSA,L1,NY,NX) = trcSalt_soHml_vr(NTSA,L1,NY,NX)+FXH
        trcSalt_soHml_vr(NTSA,L0,NY,NX) = trcSalt_soHml_vr(NTSA,L0,NY,NX)-FXH
      ENDDO
    ENDIF
!
  ENDIF
  end subroutine MoveMacPoreSolute

!------------------------------------------------------------------------------------------

  subroutine MoveBandSolute(L,L0,L1,NY,NX,FX,FWO,FO)
  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8), intent(in) :: FX,FWO,FO
  real(r8) :: FAO,FCO
  integer  :: NTSAB,NTP,NTX
  real(r8) :: FXB,FXP,FXXC,FXXA
  real(r8) :: FXH1POB,FXH2POB

! only phosphrous is considered below because there is no gaseous phase.
  FXH1POB                           = FWO*trcs_solml_vr(ids_H1PO4B,L0,NY,NX)
  trcs_solml_vr(ids_H1PO4B,L1,NY,NX) = trcs_solml_vr(ids_H1PO4B,L1,NY,NX)+FXH1POB
  trcs_solml_vr(ids_H1PO4B,L0,NY,NX) = trcs_solml_vr(ids_H1PO4B,L0,NY,NX)-FXH1POB

  FXH2POB                           = FWO*trcs_solml_vr(ids_H2PO4B,L0,NY,NX)
  trcs_solml_vr(ids_H2PO4B,L1,NY,NX) = trcs_solml_vr(ids_H2PO4B,L1,NY,NX)+FXH2POB
  trcs_solml_vr(ids_H2PO4B,L0,NY,NX) = trcs_solml_vr(ids_H2PO4B,L0,NY,NX)-FXH2POB

  IF(salt_model)THEN
    DO NTSAB=idsaltb_beg,idsaltb_end
      FXB                              = FWO*trcSalt_solml_vr(NTSAB,L0,NY,NX)
      trcSalt_solml_vr(NTSAB,L1,NY,NX) = trcSalt_solml_vr(NTSAB,L1,NY,NX)+FXB
      trcSalt_solml_vr(NTSAB,L0,NY,NX) = trcSalt_solml_vr(NTSAB,L0,NY,NX)-FXB
    ENDDO
  ENDIF
!
!     SOIL ADSORBED CATIONS IN BAND, NON-BAND
!
  IF(L0.EQ.L .OR. CEC_vr(L0,NY,NX).LE.ZERO)THEN
    FCO=FO
  ELSE
    FCO=AMIN1(0.5_r8,FO*CEC_vr(L1,NY,NX)/CEC_vr(L0,NY,NX))
  ENDIF
  DO NTX=idx_beg,idx_cation_end
    FXXC                        = FCO*trcx_solml_vr(NTX,L0,NY,NX)
    trcx_solml_vr(NTX,L1,NY,NX) = trcx_solml_vr(NTX,L1,NY,NX)+FXXC
    trcx_solml_vr(NTX,L0,NY,NX) = trcx_solml_vr(NTX,L0,NY,NX)-FXXC
  ENDDO
!
!     SOIL ADSORBED ANIONS IN BAND, NON-BAND
!
  IF(L0.EQ.L.OR.AEC_vr(L0,NY,NX).LE.ZERO)THEN
    FAO=FO
  ELSE
    FAO=AMIN1(0.5_r8,FO*AEC_vr(L1,NY,NX)/AEC_vr(L0,NY,NX))
  ENDIF

  DO NTX=idx_cation_end+1,idx_end
    FXXA                        = FAO*trcx_solml_vr(NTX,L0,NY,NX)
    trcx_solml_vr(NTX,L1,NY,NX) = trcx_solml_vr(NTX,L1,NY,NX)+FXXA
    trcx_solml_vr(NTX,L0,NY,NX) = trcx_solml_vr(NTX,L0,NY,NX)-FXXA
  ENDDo
!
!     SOIL PRECIPITATES IN BAND, NON-BAND
!
  DO NTP=idsp_beg,idsp_end
    FXP                           = AMIN1(FX*trcp_saltpml_vr(NTP,L,NY,NX),trcp_saltpml_vr(NTP,L0,NY,NX))
    trcp_saltpml_vr(NTP,L1,NY,NX) = trcp_saltpml_vr(NTP,L1,NY,NX)+FXP
    trcp_saltpml_vr(NTP,L0,NY,NX) = trcp_saltpml_vr(NTP,L0,NY,NX)-FXP
  ENDDO

  end subroutine MoveBandSolute

!------------------------------------------------------------------------------------------


  Subroutine MoveDisolvGas(L0,L1,NY,NX,FX,FWO)
  implicit none
  integer, intent(in) :: L0,L1,NY,NX
  real(r8), intent(in) :: FX,FWO

  real(r8) :: FXG
  integer  :: idg
!
! Gas-phase NH3 moves here; dissolved NH3 and NH3B move in MoveFertSalt.
  DO idg=idg_beg,idg_NH3
    FXG                        = FWO*trcg_gasml_vr(idg,L0,NY,NX)
    trcg_gasml_vr(idg,L1,NY,NX) = trcg_gasml_vr(idg,L1,NY,NX)+FXG
    trcg_gasml_vr(idg,L0,NY,NX) = trcg_gasml_vr(idg,L0,NY,NX)-FXG
  ENDDO

  DO idg=idg_beg,idg_NH3-1
    FXG                        = FWO*trcs_solml_vr(idg,L0,NY,NX)
    trcs_solml_vr(idg,L1,NY,NX) = trcs_solml_vr(idg,L1,NY,NX)+FXG
    trcs_solml_vr(idg,L0,NY,NX) = trcs_solml_vr(idg,L0,NY,NX)-FXG
  ENDDO

  end Subroutine MoveDisolvGas

!------------------------------------------------------------------------------------------

  subroutine MoveFertSalt(L,L0,L1,NY,NX,FX,FWO)

  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8), intent(in):: FX,FWO
  real(r8) :: FXZ,FXNUT,FXZN
  integer  :: ids,idsalt,NTF

! begin_execution
  DO NTF=ifertn_beg,ifertn_end
    FXZN                             = AMIN1(FX*FertN_mole_soil_vr(NTF,L,NY,NX),FertN_mole_soil_vr(NTF,L0,NY,NX))
    FertN_mole_soil_vr(NTF,L1,NY,NX) = FertN_mole_soil_vr(NTF,L1,NY,NX)+FXZN
    FertN_mole_soil_vr(NTF,L0,NY,NX) = FertN_mole_soil_vr(NTF,L0,NY,NX)-FXZN
  ENDDO
  ! Broadcast P can be carried from the surface pool into soil.
  if(L0==0)then
    FXZN=FX*FertP_mole_soil_vr(0,NY,NX)
  else
    FXZN=AMIN1(FX*FertP_mole_soil_vr(L,NY,NX),FertP_mole_soil_vr(L0,NY,NX))
  endif
  FertP_mole_soil_vr(L1,NY,NX)=FertP_mole_soil_vr(L1,NY,NX)+FXZN
  FertP_mole_soil_vr(L0,NY,NX)=FertP_mole_soil_vr(L0,NY,NX)-FXZN
  IF (L0>0) then
    DO NTF=ifertnb_beg,ifertnb_end
      FXZN                        = AMIN1(FX*FertN_mole_Band_vr(NTF,L,NY,NX),FertN_mole_Band_vr(NTF,L0,NY,NX))
      FertN_mole_Band_vr(NTF,L1,NY,NX) = FertN_mole_Band_vr(NTF,L1,NY,NX)+FXZN
      FertN_mole_Band_vr(NTF,L0,NY,NX) = FertN_mole_Band_vr(NTF,L0,NY,NX)-FXZN
    ENDDO
    FXZN = AMIN1(FX*FertP_mole_band_vr(L,NY,NX),FertP_mole_band_vr(L0,NY,NX))    
    FertP_mole_band_vr(L1,NY,NX) = FertP_mole_band_vr(L1,NY,NX)+FXZN
    FertP_mole_band_vr(L0,NY,NX) = FertP_mole_band_vr(L0,NY,NX)-FXZN
  endif

!
!     SOIL N,P SOLUTES IN BAND, NON-BAND
!
  DO ids=ids_nuts_beg,ids_nuts_end
    if(ids/=ids_H2PO4B .and. ids/=ids_H1PO4B)THEN
      FXNUT                      = FWO*trcs_solml_vr(ids,L0,NY,NX)
      trcs_solml_vr(ids,L1,NY,NX) = trcs_solml_vr(ids,L1,NY,NX)+FXNUT
      trcs_solml_vr(ids,L0,NY,NX) = trcs_solml_vr(ids,L0,NY,NX)-FXNUT
    ENDIF
  ENDDO

!
!     SOIL SALT SOLUTES
!
  IF(salt_model)THEN
    DO idsalt=idsalt_beg,idsalt_end
      FXZ                             = FWO*trcSalt_solml_vr(idsalt,L0,NY,NX)
      trcSalt_solml_vr(idsalt,L1,NY,NX) = trcSalt_solml_vr(idsalt,L1,NY,NX)+FXZ
      trcSalt_solml_vr(idsalt,L0,NY,NX) = trcSalt_solml_vr(idsalt,L0,NY,NX)-FXZ
    ENDDO
  ENDIF
  end subroutine MoveFertSalt

!------------------------------------------------------------------------------------------

  subroutine MoveFertMinerals(L,L0,L1,NY,NX,FX,FO,FHO)
  implicit none
  integer, intent(in) :: L,L0,L1,NY,NX
  real(r8),intent(in) :: FX,FO
  real(r8),intent(out) :: FHO !donor macropore fraction; used for water/ice, heat, pore volume, and solutes
  real(r8) :: FBO

  real(r8) :: FXSAND,FXSILT,FXCLAY,FXROCK
  real(r8) :: FXWDNHB,FXDPNHB,FXWDNOB
  real(r8) :: FXDPNOB,FXWDPOB,FXDPPOB
  real(r8) :: FXVOLAH

  ! No macropore transfer unless the active branch below assigns a fraction.
  FHO=0._r8

  IF(DLYR_3D(3,L1,NY,NX).GT.ZERO.AND.DLYR_3D(3,L0,NY,NX).GT.ZERO)THEN
!
!     SOIL FERTILIZER BANDS
!
    IF(iFertNH4Band_col(NY,NX).EQ.ifert_on .AND. ROWSpaceNH4_col(NY,NX).GT.0.0_r8)THEN
      IF(L.EQ.NU_col(NY,NX) .OR. CumDepz2LayBottom_vr(L-1,NY,NX).LT.BandDepthNH4_col(NY,NX))THEN
        WDNHBDL                   = BandWidthNH4_vr(L,NY,NX)*DLYR_3D(3,L,NY,NX)
        WDNHBD0                   = BandWidthNH4_vr(L0,NY,NX)*DLYR_3D(3,L0,NY,NX)
        WDNHBD1                   = BandWidthNH4_vr(L1,NY,NX)*DLYR_3D(3,L1,NY,NX)
        FXWDNHB                   = AMIN1(FX*WDNHBDL,WDNHBD0)
        WDNHBD1                   = WDNHBD1+FXWDNHB
        WDNHBD0                   = WDNHBD0-FXWDNHB
        BandWidthNH4_vr(L1,NY,NX) = WDNHBD1/DLYR_3D(3,L1,NY,NX)
        BandWidthNH4_vr(L0,NY,NX) = WDNHBD0/DLYR_3D(3,L0,NY,NX)
        IF(CumDepz2LayBottom_vr(L,NY,NX).GE.BandDepthNH4_col(NY,NX))THEN
          FXDPNHB                       = AMIN1(FX*BandThicknessNH4_vr(L,NY,NX),BandThicknessNH4_vr(L0,NY,NX))
          BandThicknessNH4_vr(L1,NY,NX) = BandThicknessNH4_vr(L1,NY,NX)+FXDPNHB
          BandThicknessNH4_vr(L0,NY,NX) = BandThicknessNH4_vr(L0,NY,NX)-FXDPNHB
        ENDIF
        trcs_VLN_vr(ids_NH4B,L1,NY,NX)=AZMAX1(AMIN1(0.999_r8,BandWidthNH4_vr(L1,NY,NX) &
          /ROWSpaceNH4_col(NY,NX)*BandThicknessNH4_vr(L1,NY,NX)/DLYR_3D(3,L1,NY,NX)))
        trcs_VLN_vr(ids_NH4B,L0,NY,NX)=AZMAX1(AMIN1(0.999_r8,BandWidthNH4_vr(L0,NY,NX) &
          /ROWSpaceNH4_col(NY,NX)*BandThicknessNH4_vr(L0,NY,NX)/DLYR_3D(3,L0,NY,NX)))
        trcs_VLN_vr(ids_NH4,L1,NY,NX) = 1.0_r8-trcs_VLN_vr(ids_NH4B,L1,NY,NX)
        trcs_VLN_vr(ids_NH4,L0,NY,NX) = 1.0_r8-trcs_VLN_vr(ids_NH4B,L0,NY,NX)

        trcs_VLN_vr(idg_NH3B,L1,NY,NX) = trcs_VLN_vr(ids_NH4B,L1,NY,NX)
        trcs_VLN_vr(idg_NH3B,L0,NY,NX) = trcs_VLN_vr(ids_NH4B,L0,NY,NX)
        trcs_VLN_vr(idg_NH3,L1,NY,NX)  = trcs_VLN_vr(ids_NH4,L1,NY,NX)
        trcs_VLN_vr(idg_NH3,L0,NY,NX)  = trcs_VLN_vr(ids_NH4,L0,NY,NX)
      ENDIF
    ENDIF
    IF(iFertNO3Band_col(NY,NX).EQ.ifert_off .AND. ROWSpaceNO3Band_col(NY,NX).GT.0.0_r8)THEN
      IF(L.EQ.NU_col(NY,NX) .OR. CumDepz2LayBottom_vr(L-1,NY,NX).LT.BandDepthNO3_col(NY,NX))THEN
        WDNOBDL                   = BandWidthNO3_vr(L,NY,NX)*DLYR_3D(3,L,NY,NX)
        WDNOBD0                   = BandWidthNO3_vr(L0,NY,NX)*DLYR_3D(3,L0,NY,NX)
        WDNOBD1                   = BandWidthNO3_vr(L1,NY,NX)*DLYR_3D(3,L1,NY,NX)
        FXWDNOB                   = AMIN1(FX*WDNOBDL,WDNOBD0)
        WDNOBD1                   = WDNOBD1+FXWDNOB
        WDNOBD0                   = WDNOBD0-FXWDNOB
        BandWidthNO3_vr(L1,NY,NX) = WDNOBD1/DLYR_3D(3,L1,NY,NX)
        BandWidthNO3_vr(L0,NY,NX) = WDNOBD0/DLYR_3D(3,L0,NY,NX)
        IF(CumDepz2LayBottom_vr(L,NY,NX).GE.BandDepthNO3_col(NY,NX))THEN
          FXDPNOB                       = AMIN1(FX*BandThicknessNO3_vr(L,NY,NX),BandThicknessNO3_vr(L0,NY,NX))
          BandThicknessNO3_vr(L1,NY,NX) = BandThicknessNO3_vr(L1,NY,NX)+FXDPNOB
          BandThicknessNO3_vr(L0,NY,NX) = BandThicknessNO3_vr(L0,NY,NX)-FXDPNOB
        ENDIF
        trcs_VLN_vr(ids_NO3B,L1,NY,NX)=AZMAX1(AMIN1(0.999_r8,BandWidthNO3_vr(L1,NY,NX) &
          /ROWSpaceNO3Band_col(NY,NX)*BandThicknessNO3_vr(L1,NY,NX)/DLYR_3D(3,L1,NY,NX)))
        trcs_VLN_vr(ids_NO3B,L0,NY,NX)=AZMAX1(AMIN1(0.999_r8,BandWidthNO3_vr(L0,NY,NX) &
          /ROWSpaceNO3Band_col(NY,NX)*BandThicknessNO3_vr(L0,NY,NX)/DLYR_3D(3,L0,NY,NX)))
        trcs_VLN_vr(ids_NO3,L1,NY,NX)=1.0_r8-trcs_VLN_vr(ids_NO3B,L1,NY,NX)
        trcs_VLN_vr(ids_NO3,L0,NY,NX)=1.0_r8-trcs_VLN_vr(ids_NO3B,L0,NY,NX)

        trcs_VLN_vr(ids_NO2,L1,NY,NX)  = trcs_VLN_vr(ids_NO3,L1,NY,NX)
        trcs_VLN_vr(ids_NO2,L0,NY,NX)  = trcs_VLN_vr(ids_NO3,L0,NY,NX)
        trcs_VLN_vr(ids_NO2B,L1,NY,NX) = trcs_VLN_vr(ids_NO3B,L1,NY,NX)
        trcs_VLN_vr(ids_NO2B,L0,NY,NX) = trcs_VLN_vr(ids_NO3B,L0,NY,NX)
      ENDIF
    ENDIF
    IF(iFertPO4Band_col(NY,NX).EQ.ifert_on .AND. ROWSpacePO4_col(NY,NX).GT.0.0_r8)THEN
      IF(L.EQ.NU_col(NY,NX) .OR. CumDepz2LayBottom_vr(L-1,NY,NX).LT.BandDepthPO4_col(NY,NX))THEN
        WDPOBDL                   = BandWidthPO4_vr(L,NY,NX)*DLYR_3D(3,L,NY,NX)
        WDPOBD0                   = BandWidthPO4_vr(L0,NY,NX)*DLYR_3D(3,L0,NY,NX)
        WDPOBD1                   = BandWidthPO4_vr(L1,NY,NX)*DLYR_3D(3,L1,NY,NX)
        FXWDPOB                   = AMIN1(FX*WDPOBDL,WDPOBD0)
        WDPOBD1                   = WDPOBD1+FXWDPOB
        WDPOBD0                   = WDPOBD0-FXWDPOB
        BandWidthPO4_vr(L1,NY,NX) = WDPOBD1/DLYR_3D(3,L1,NY,NX)
        BandWidthPO4_vr(L0,NY,NX) = WDPOBD0/DLYR_3D(3,L0,NY,NX)
        IF(CumDepz2LayBottom_vr(L,NY,NX).GE.BandDepthPO4_col(NY,NX))THEN
          FXDPPOB                       = AMIN1(FX*BandThicknessPO4_vr(L,NY,NX),BandThicknessPO4_vr(L0,NY,NX))
          BandThicknessPO4_vr(L1,NY,NX) = BandThicknessPO4_vr(L1,NY,NX)+FXDPPOB
          BandThicknessPO4_vr(L0,NY,NX) = BandThicknessPO4_vr(L0,NY,NX)-FXDPPOB
        ENDIF
        trcs_VLN_vr(ids_H1PO4B,L1,NY,NX)=AZMAX1(AMIN1(0.999,BandWidthPO4_vr(L1,NY,NX) &
          /ROWSpacePO4_col(NY,NX)*BandThicknessPO4_vr(L1,NY,NX)/DLYR_3D(3,L1,NY,NX)))
        trcs_VLN_vr(ids_H1PO4B,L0,NY,NX)=AZMAX1(AMIN1(0.999,BandWidthPO4_vr(L0,NY,NX) &
          /ROWSpacePO4_col(NY,NX)*BandThicknessPO4_vr(L0,NY,NX)/DLYR_3D(3,L0,NY,NX)))
        trcs_VLN_vr(ids_H1PO4,L1,NY,NX)=1.0_r8-trcs_VLN_vr(ids_H1PO4B,L1,NY,NX)
        trcs_VLN_vr(ids_H1PO4,L0,NY,NX)=1.0_r8-trcs_VLN_vr(ids_H1PO4B,L0,NY,NX)

        trcs_VLN_vr(ids_H2PO4B,L1,NY,NX) = trcs_VLN_vr(ids_H1PO4B,L1,NY,NX)
        trcs_VLN_vr(ids_H2PO4B,L0,NY,NX) = trcs_VLN_vr(ids_H1PO4B,L0,NY,NX)
        trcs_VLN_vr(ids_H2PO4,L1,NY,NX)  = trcs_VLN_vr(ids_H1PO4,L1,NY,NX)
        trcs_VLN_vr(ids_H2PO4,L0,NY,NX)  = trcs_VLN_vr(ids_H1PO4,L0,NY,NX)
      ENDIF
    ENDIF
  ENDIF
!
!     SOIL MINERALS
!
  IF(L0.EQ.L.OR.SoiBulkDensityt0_vr(L0,NY,NX).LE.ZERO)THEN
    FBO=FX
  ELSE
    FBO=AMIN1(0.1,FX*SoiBulkDensityt0_vr(L1,NY,NX)/SoiBulkDensityt0_vr(L0,NY,NX))
  ENDIF
!     SoilBulkDensity_vr(L1,NY,NX)=(1.0-FO)*SoilBulkDensity_vr(L1,NY,NX)+FO*SoiBulkDensityt0_vr(L0,NY,NX)
  PH_vr(L1,NY,NX)   = (1.0_r8-FO)*PH_vr(L1,NY,NX)+FO*PH_vr(L0,NY,NX)
  FXSAND            = FBO*SAND_vr(L0,NY,NX)
  SAND_vr(L1,NY,NX) = SAND_vr(L1,NY,NX)+FXSAND
  SAND_vr(L0,NY,NX) = SAND_vr(L0,NY,NX)-FXSAND
  FXSILT            = FBO*SILT_vr(L0,NY,NX)
  SILT_vr(L1,NY,NX) = SILT_vr(L1,NY,NX)+FXSILT
  SILT_vr(L0,NY,NX) = SILT_vr(L0,NY,NX)-FXSILT
  FXCLAY            = FBO*CLAY_vr(L0,NY,NX)
  CLAY_vr(L1,NY,NX) = CLAY_vr(L1,NY,NX)+FXCLAY
  CLAY_vr(L0,NY,NX) = CLAY_vr(L0,NY,NX)-FXCLAY
  FXROCK            = FBO*ROCK_vr(L0,NY,NX)
  ROCK_vr(L1,NY,NX) = ROCK_vr(L1,NY,NX)+FXROCK
  ROCK_vr(L0,NY,NX) = ROCK_vr(L0,NY,NX)-FXROCK
!
!     MACROPORE FRACTION AND PORE VOLUME (WATER/ICE AND HEAT MOVE IN MoveHeatWat)
!
  IF(SoilFracAsMacP_vr(L1,NY,NX).GT.ZERO.AND.SoilFracAsMacP_vr(L0,NY,NX).GT.ZERO)THEN
    IF(L0.EQ.L.OR.SoilFracAsMacPt0_vr(L0,NY,NX).LE.ZERO)THEN
      FHO=FO
    ELSE
      FHO=AMIN1(0.5_r8,FO*SoilFracAsMacPt0_vr(L1,NY,NX)/SoilFracAsMacPt0_vr(L0,NY,NX))
    ENDIF
    SoilFracAsMacP_vr(L1,NY,NX) = (1.0_r8-FO)*SoilFracAsMacP_vr(L1,NY,NX)+FO*SoilFracAsMacP_vr(L0,NY,NX)
    FXVOLAH                     = FHO*VLMacroPores_vr(L0,NY,NX)
    VLMacroPores_vr(L1,NY,NX)         = VLMacroPores_vr(L1,NY,NX)+FXVOLAH
    VLMacroPores_vr(L0,NY,NX)         = VLMacroPores_vr(L0,NY,NX)-FXVOLAH
  ENDIF
  end subroutine MoveFertMinerals

!------------------------------------------------------------------------------------------

  subroutine MoveHeatWat(I,J,L,L0,L1,NY,NX,FO,FX,XVOLWP,FHO,FWO)
  implicit none
  integer, intent(in) :: I,J  
  integer, intent(in) :: L,NY,NX
  integer, intent(in) :: L0    !source layer
  integer, intent(in) :: L1    !target layer
  real(r8), intent(in) :: FO
  real(r8), intent(in) :: FX
  real(r8), intent(in) :: FHO      !donor macropore fraction from MoveFertMinerals; zero for litter
  real(r8), intent(in) :: XVOLWP
  real(r8), intent(out):: FWO      !mixing fraction of organic material
  real(r8) :: FXVOLW,FXVOLWH,FXVOLIH,FXVHC
  real(r8) :: FXENGY,FXVOLI,FXVLSoilMicP,FXVOLWX,FXVHCM
  real(r8) :: ENGY0,ENGY1

  IF(L0.EQ.L .OR. POROSI_vr(L0,NY,NX).LE.ZERO)THEN
    FWO=FO
  ELSE
    FWO=AMIN1(0.5,FO*POROSI_vr(L1,NY,NX)/POROSI_vr(L0,NY,NX))
  ENDIF
!             FXSCNV=FWO*SatHydroCondVert_vr(L0,NY,NX)
!             SatHydroCondVert_vr(L1,NY,NX)=SatHydroCondVert_vr(L1,NY,NX)+FXSCNV
!             SatHydroCondVert_vr(L0,NY,NX)=SatHydroCondVert_vr(L0,NY,NX)-FXSCNV
!             FXSCNH=FWO*SatHydroCondHrzn_vr(L0,NY,NX)
!             SatHydroCondHrzn_vr(L1,NY,NX)=SatHydroCondHrzn_vr(L1,NY,NX)+FXSCNH
!             SatHydroCondHrzn_vr(L0,NY,NX)=SatHydroCondHrzn_vr(L0,NY,NX)-FXSCNH

  IF(L0.EQ.0)THEN
    FXVOLW=FX*AZMAX1(XVOLWP+VLWatHeldCapSurf_col(NY,NX))
    Qinflx2Soil_col(NY,NX) = Qinflx2Soil_col(NY,NX)+FXVOLW
  ELSE
    FXVOLW=FWO*VLWatMicP_vr(L0,NY,NX)
  ENDIF
  VLWatMicP_vr(L1,NY,NX) = VLWatMicP_vr(L1,NY,NX)+FXVOLW
  VLWatMicP_vr(L0,NY,NX) = VLWatMicP_vr(L0,NY,NX)-FXVOLW
!     IF(VLiceMicP_vr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
  FXVOLI                 = FWO*VLiceMicP_vr(L0,NY,NX)
  VLiceMicP_vr(L1,NY,NX) = VLiceMicP_vr(L1,NY,NX)+FXVOLI
  VLiceMicP_vr(L0,NY,NX) = VLiceMicP_vr(L0,NY,NX)-FXVOLI
!     ENDIF
!     FXVOLA=FWO*VLMicP_vr(L0,NY,NX)
!     IF(L1.NE.NU_col(NY,NX))THEN
!     VLMicP_vr(L1,NY,NX)=VLMicP_vr(L1,NY,NX)+FXVOLA
!     ENDIF
!     IF(L0.NE.NU_col(NY,NX))THEN
!     VLMicP_vr(L0,NY,NX)=VLMicP_vr(L0,NY,NX)-FXVOLA
!     ENDIF
  FXVLSoilMicP                    = FWO*VLSoilMicP_vr(L0,NY,NX)
  VLSoilMicP_vr(L1,NY,NX)         = VLSoilMicP_vr(L1,NY,NX)+FXVLSoilMicP
  VLSoilMicP_vr(L0,NY,NX)         = VLSoilMicP_vr(L0,NY,NX)-FXVLSoilMicP
  FXVOLWX                         = FWO*VLWatMicPX_vr(L0,NY,NX)
  VLWatMicPX_vr(L1,NY,NX)         = VLWatMicPX_vr(L1,NY,NX)+FXVOLWX
  VLWatMicPX_vr(L0,NY,NX)         = VLWatMicPX_vr(L0,NY,NX)-FXVOLWX
  FXVHCM                          = FWO*VHeatCapSolidSoil_vr(L0,NY,NX)
  VHeatCapSolidSoil_vr(L1,NY,NX) = VHeatCapSolidSoil_vr(L1,NY,NX)+FXVHCM
  VHeatCapSolidSoil_vr(L0,NY,NX) = VHeatCapSolidSoil_vr(L0,NY,NX)-FXVHCM

  ! Transfer macropore water/ice and their sensible heat before updating temperatures.
  FXVOLWH=0._r8
  FXVOLIH=0._r8
  IF(L0.NE.0)THEN
    FXVOLWH                     = FHO*VLWatMacP_vr(L0,NY,NX)
    VLWatMacP_vr(L1,NY,NX)      = VLWatMacP_vr(L1,NY,NX)+FXVOLWH
    VLWatMacP_vr(L0,NY,NX)      = VLWatMacP_vr(L0,NY,NX)-FXVOLWH
    FXVOLIH                     = FHO*VLiceMacP_vr(L0,NY,NX)
    VLiceMacP_vr(L1,NY,NX)      = VLiceMacP_vr(L1,NY,NX)+FXVOLIH
    VLiceMacP_vr(L0,NY,NX)      = VLiceMacP_vr(L0,NY,NX)-FXVOLIH
  ENDIF
  FXVHC                           = FXVHCM+cpw*(FXVOLW+FXVOLWH)+cpi*(FXVOLI+FXVOLIH)
  FXENGY                          = TKS_vr(L0,NY,NX)*FXVHC
  ENGY1                           = VHeatCapacity_vr(L1,NY,NX)*TKS_vr(L1,NY,NX)+FXENGY
  ENGY0                           = VHeatCapacity_vr(L0,NY,NX)*TKS_vr(L0,NY,NX)-FXENGY
  VHeatCapacity_vr(L1,NY,NX)      = VHeatCapacity_vr(L1,NY,NX)+FXVHC
  VHeatCapacity_vr(L0,NY,NX)      = VHeatCapacity_vr(L0,NY,NX)-FXVHC
  IF(VHeatCapacity_vr(L1,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS_vr(L1,NY,NX)=ENGY1/VHeatCapacity_vr(L1,NY,NX)
  ELSE
    TKS_vr(L1,NY,NX)=TKS_vr(L,NY,NX)
  ENDIF
  TCS_vr(L1,NY,NX)=units%Kelvin2Celcius(TKS_vr(L1,NY,NX))
  IF(VHeatCapacity_vr(L0,NY,NX).GT.ZEROS(NY,NX))THEN
    TKS_vr(L0,NY,NX)=ENGY0/VHeatCapacity_vr(L0,NY,NX)
  ELSE
    TKS_vr(L0,NY,NX)=TKS_vr(L,NY,NX)
  ENDIF
  TCS_vr(L0,NY,NX)=units%Kelvin2Celcius(TKS_vr(L0,NY,NX))
  end subroutine MoveHeatWat

end module SoilLayerDynMod
