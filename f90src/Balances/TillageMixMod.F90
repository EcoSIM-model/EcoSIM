module TillageMixMod
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only: micpar
  use EcoSIMConfig,     only: ndbiomcp => NumDeadMicrbCompts
  use UnitMod,          only: units
  use minimathmod,      only: AZMAX1
  use MicrobialDiagMod, only: sumLitrOMLayL, sumORGMLayL
  use EcoSIMCtrlMod,    only: salt_model
  use MiniFuncMod,      only: gOC_to_m3_OM
  use DebugToolMod
  use EcosimConst
  use SOMDataType
  USE GridDataType
  USE GridConsts
  USE SoilBGCDataType
  use FertilizerDataType
  use AqueChemDatatype
  USE MicrobialDataType
  use SurfLitterDataType
  use SoilHeatDataType
  use SoilWaterDataType
  use BalanceCheckDataType, only: TillageLitr2SoilH2O_col
  use SoilPhysDataType
  USE SoilPropertyDataType
  USE ClimForcDataType
  USE FlagDataType
  use EcoSIMCtrlDataType
  USE EcoSimSumDataType
  implicit none
  character(len=*),private, parameter :: mod_filename =&
     __FILE__
  public :: ApplyTillageMixing
  contains

!------------------------------------------------------------------------------------------

  subroutine ApplyTillageMixing(I,J,NY,NX)
  implicit none
  integer, intent(in) :: I,J,NY,NX

  character(len=*), parameter :: subname='ApplyTillageMixing'
  integer :: K,N,M,L,LL,NGL,NTS
  real(r8) :: TDOMH(idom_beg:idom_end,0:micpar%NumOfLitrCmplxs)
  real(r8) :: FracTillCorp,CORP0,XCORP0,DCORPZ
  real(r8) :: FracMacPoreMixed,TransferredDOM
  REAL(R8) :: ZNHUX0,ZNHUXI,ZNFNX0
  real(r8) :: TZNFNG,TZNFN2
  real(r8) :: TNFNIH,TZNFNI
  real(r8) :: ENGYV(0:JZ),ENGYL,ENGYM(0:JZ)
  real(r8) :: ORGM(1:NumPlantChemElms),litrOM(1:NumPlantChemElms)
  real(r8) :: HFLXD
  real(r8) :: LitterWaterBeforeMix
  real(r8) :: TX,TL
  integer  :: NTX,NTP,NTG,NTSA,NTN,idom,MID,NE
  real(r8) :: FI(JZ),TI(JZ),TI1,FI1
  logical :: TillTest,tillTestL

!     begin_execution
!
  call PrintInfo('beg '//subname)
  !DCORP=soil mixing fraction with tillage
  !XCORP=factor for surface litter incorporation and soil mixing

  TillTest=J.EQ.INT(SolarNoonHour_col(NY,NX)) .AND. XTillCorp_col(NY,NX).LT.1.0_r8 .AND. DepzCorp_col(I,NY,NX).GT.0.0_r8
  if(.not. TillTest)return
!
  iResetSoilProf_col(NY,NX) = itrue
!   tillage depth cannot exceed whole soil column thickness
  DCORPZ = AMIN1(DepzCorp_col(I,NY,NX),CumSoilThickness_vr(NL_col(NY,NX),NY,NX))

!
!     ACCUMULATE STATE VARIABLES IN SURFACE RESIDUE FOR ADDITION
!     TO SOIL IN TILLAGE MIXING ZONE
!
!     IF(SoilOrgM_vr(ielmc,0,NY,NX).GT.ZEROS(NY,NX))THEN
!     XCORP0=AMAX1(XTillCorp_col(NY,NX),AMIN1(1.0,
!    2(VHeatCapLitRMin_col(NY,NX)/cpo)/SoilOrgM_vr(ielmc,0,NY,NX)))
  XCORP0=AMAX1(0.001_r8,XTillCorp_col(NY,NX))
!     ELSE
!     XCORP0=1.0
!     ENDIF
!   XCORP0: Fraction of material not being mixed

  call DeriveTillageProf(NY,NX,DCORPZ,LL,FI,TI)

  ! FI weights the tilled thickness in this column. Intensive properties
  ! are blended as properties; actual pools below retain mass redistribution.
  ! Conductivity/selectivity blending is an empirical arithmetic mixing rule,
  ! not an effective-conductivity or exchange-equilibrium calculation.
  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,SoiBulkDensityt0_vr(1:JZ,NY,NX))

  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,FieldCapacity_vr(1:JZ,NY,NX))

  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,WiltPoint_vr(1:JZ,NY,NX))

  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,SatHydroCondVert_vr(1:JZ,NY,NX))

  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,SatHydroCondHrzn_vr(1:JZ,NY,NX))

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,SAND_vr(1:JZ,NY,NX))

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,SILT_vr(1:JZ,NY,NX))

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,CLAY_vr(1:JZ,NY,NX))

  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,GKC4_vr(1:JZ,NY,NX))

  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,GKCA_vr(1:JZ,NY,NX))

  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,GKCM_vr(1:JZ,NY,NX))

  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,GKCN_vr(1:JZ,NY,NX))

  call MixIntensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,GKCK_vr(1:JZ,NY,NX))

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,FertN_mole_soil_vr(ifertn_beg:ifertn_end,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,FertN_mole_Band_vr(ifertnb_beg:ifertnb_end,1:JZ,NY,NX))

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,FertP_mole_soil_vr(0:JZ,NY,NX),XCORP0)

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,FertP_mole_band_vr(1:JZ,NY,NX))

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,trcs_solml_vr(ids_beg:idg_NH3,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,trcs_solml_vr(ids_nut_beg:ids_nuts_end,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,trcs_solml_vr(idg_NH3B:ids_nutb_end,1:JZ,NY,NX))

  if(salt_model)then
    CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcSalt_solml_vr(idsalt_beg:idsalt_end,0:JZ,NY,NX),XCORP0)

    CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcSalt_solml_vr(idsalt_pband_beg:idsalt_pband_end,1:JZ,NY,NX))
  endif

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcp_saltpml_vr(idsp_beg:idsp_CaSO4,1:JZ,NY,NX))

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcp_saltpml_vr(idsp_beg_band:idsp_end,1:JZ,NY,NX))

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcp_saltpml_vr(idsp_psoi_beg:idsp_psoi_end,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcx_solml_vr(idx_Hp:idx_cation_end,1:JZ,NY,NX))

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcx_solml_vr(idx_CEC:idx_CEC,1:JZ,NY,NX))

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcx_solml_vr(idx_AEC:idx_AEC,1:JZ,NY,NX))

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcx_solml_vr(idx_NH4:idx_NH4,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcx_solml_vr(idx_AEC+1:idx_anion_soil_end,0:JZ,NY,NX),XCORP0)

  ! Band-bound anions have no surface-litter pool.
  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcx_solml_vr(idx_anion_soil_end+1:idx_end,1:JZ,NY,NX))

  CALL MixExtensive2D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, trcg_gasml_vr(idg_beg:idg_NH3,1:JZ,NY,NX))

  CALL MixExtensive3D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, mBiomeAutor_vr(:,:,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive3D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, DOM_MicP_vr(:,1:micpar%NumOfLitrCmplxs,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive3D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, DOM_MicP_vr(:,micpar%NumOfLitrCmplxs+1:jcplx,1:JZ,NY,NX))

  CALL MixExtensive3D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, SorbedOM_vr(:,1:micpar%NumOfLitrCmplxs,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive3D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, SorbedOM_vr(:,micpar%NumOfLitrCmplxs+1:jcplx,1:JZ,NY,NX))

  CALL MixExtensive3D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, SolidOMAct_vr(:,1:micpar%NumOfLitrCmplxs,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive3D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, SolidOMAct_vr(:,micpar%NumOfLitrCmplxs+1:jcplx,1:JZ,NY,NX))

  CALL MixExtensive4D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, mBiomeHeter_vr(:,:,1:micpar%NumOfLitrCmplxs,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive4D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, mBiomeHeter_vr(:,:,micpar%NumOfLitrCmplxs+1:jcplx,1:JZ,NY,NX))

  CALL MixExtensive4D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, OMBioResdu_vr(:,:,1:micpar%NumOfLitrCmplxs,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive4D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, OMBioResdu_vr(:,:,micpar%NumOfLitrCmplxs+1:jcplx,1:JZ,NY,NX))

  CALL MixExtensive4D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, SolidOM_vr(:,:,1:micpar%NumOfLitrCmplxs,0:JZ,NY,NX),XCORP0)

  CALL MixExtensive4D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL, SolidOM_vr(:,:,micpar%NumOfLitrCmplxs+1:jcplx,1:JZ,NY,NX))

  !build the energy profile
  !ENGYV: water
  !ENGYM: material
  ENGYM(0) = cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,0,NY,NX))*TKS_vr(0,NY,NX)
  ENGYV(0) = (cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX))*TKS_vr(0,NY,NX)
  VHeatCapSolidSoil_vr(0,NY,NX)=cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,0,NY,NX))
  D3000: DO  L=NU_col(NY,NX),LL
    IF(DLYR_3D(3,L,NY,NX).GT.ZERO)THEN
      ENGYM(L) = VHeatCapSolidSoil_vr(L,NY,NX)*TKS_vr(L,NY,NX)
      ENGYV(L) = (cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)) &
          +cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX)))*TKS_vr(L,NY,NX)
    ENDIF    
  ENDDO D3000
  HFLXD=ENGYM(0)+ENGYV(0)

  !mix
  LitterWaterBeforeMix=VLWatMicP_vr(0,NY,NX)+DENSICE*VLiceMicP_vr(0,NY,NX)
  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,VLWatMicP_vr(0:JZ,NY,NX),XCORP0)

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,VLiceMicP_vr(0:JZ,NY,NX),XCORP0)
  ! Record incorporation for compartment checks; it is not infiltration or external water.
  TillageLitr2SoilH2O_col(NY,NX)=TillageLitr2SoilH2O_col(NY,NX)+LitterWaterBeforeMix &
    -VLWatMicP_vr(0,NY,NX)-DENSICE*VLiceMicP_vr(0,NY,NX)

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,ENGYV(0:JZ),XCORP0)

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,ENGYM(0:JZ),XCORP0)

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,VHeatCapSolidSoil_vr(0:JZ,NY,NX),XCORP0)

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,VLWatMacP_vr(1:JZ,NY,NX))

  call MixExtensive1D(XTillCorp_col(NY,NX),TI,FI,NU_col(NY,NX),LL,VLiceMacP_vr(1:JZ,NY,NX))

  !update 
  call sumLitrOMLayL(0,NY,NX,litrOM)

  SoilOrgM_vr(1:NumPlantChemElms,0,NY,NX)=litrOM(1:NumPlantChemElms)

  OMLitrC_vr(0,NY,NX)=litrOM(ielmc)

  VHeatCapacity_vr(0,NY,NX) = cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,0,NY,NX))+cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)

  VLitR_col(NY,NX)          = VLitR_col(NY,NX)*XCORP0
  VGeomLayer_vr(0,NY,NX)    = VGeomLayer_vr(0,NY,NX)*XCORP0

  ! Empty litter has no stored energy; retain its previous temperature.
  IF(VHeatCapacity_vr(0,NY,NX).GT.0._r8)THEN
    TKS_vr(0,NY,NX) = (ENGYM(0)+ENGYV(0))/VHeatCapacity_vr(0,NY,NX)
  ENDIF
  HFLXD                     = HFLXD-ENGYM(0)-ENGYV(0)
  DO  L=NU_col(NY,NX),LL
    IF(DLYR_3D(3,L,NY,NX).GT.ZERO)THEN
      VLWatMicPX_vr(L,NY,NX)    = VLWatMicP_vr(L,NY,NX)
      VHeatCapacity_vr(L,NY,NX) = VHeatCapSolidSoil_vr(L,NY,NX)+cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)) &
          +cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))          
      TKS_vr(L,NY,NX) = (ENGYM(L)+ENGYV(L))/VHeatCapacity_vr(L,NY,NX)
    ENDIF
  ENDDO

  HEATIN_lnd    = HEATIN_lnd-HFLXD
  HeatStore_lnd = HeatStore_lnd-HFLXD

  !reset
  EnergyImpact4Erosion_col(NY,NX) = 0.0_r8
!
!     EXTENT OF MIXING

    !
  FracTillCorp                = 1.0_r8-XTillCorp_col(NY,NX)

!
!     TEMPORARY ACCUMULATORS
!
  TNFNIH            = 0.0_r8
  TZNFN2       = 0.0_r8
  TZNFNI       = 0.0_r8
  ZNHUX0       = 0.0_r8
  ZNHUXI       = 0.0_r8
  ZNFNX0       = 0.0_r8

!   CORP0: Fraction of material being mixed
  CORP0=1.0_r8-XCORP0   

  DO K=1,micpar%NumOfLitrCmplxs
    DO idom=idom_beg,idom_end
      TDOMH(idom,K)  = DOM_MacP_vr(idom,K,0,NY,NX)*CORP0
!
!     REDUCE SURFACE RESIDUE STATE VARIABLES FOR INCORPORATION
!
      DOM_MacP_vr(idom,K,0,NY,NX) = DOM_MacP_vr(idom,K,0,NY,NX)*XCORP0
    enddo 
  ENDDO

  TZNFNG = ZNFNI_vr(0,NY,NX)*CORP0
  ZNHUX0 = AMAX1(ZNHUX0,ZNHU0_vr(0,NY,NX))
  ZNHUXI = AMAX1(ZNHUXI,ZNHUI_vr(0,NY,NX))
  ZNFNX0 = AMAX1(ZNFNX0,ZNFN0_vr(0,NY,NX))

!
!     REDISTRIBUTE SOIL STATE VARIABLES DURING TILLAGE
!
!     ACCUMULATE SOIL STATE VARIABLES WITHIN TILLAGE MIXING ZONE
!
  D100: DO L=NU_col(NY,NX),NL_col(NY,NX)
    !if current layer is within the tillage depth
    tillTestL=CumSoilThickness_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX).LT.DCORPZ .AND. DLYR_3D(3,L,NY,NX).GT.ZERO
    if(.not.tillTestL)exit
!     TL: fraction of layer L over the tillage zone
!     TI: fraction of layer L being tilled
    TL = AMIN1(DLYR_3D(3,L,NY,NX),DCORPZ-(CumSoilThickness_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)))
    FI1 = TL/DCORPZ
    TI1 = AMIN1(TL/DLYR_3D(3,L,NY,NX),1._r8)

    ZNHUX0 = AMAX1(ZNHUX0,ZNHU0_vr(L,NY,NX))
    ZNHUXI = AMAX1(ZNHUXI,ZNHUI_vr(L,NY,NX))
    ZNFNX0 = AMAX1(ZNFNX0,ZNFN0_vr(L,NY,NX))
    TZNFNI = TZNFNI+ZNFNI_vr(L,NY,NX)
    LL=L      
  ENDDO D100
!
!     CHANGE SOIL STATE VARIABLES IN TILLAGE MIXING ZONE
!     TO ACCOUNT FOR REDISTRIBUTION FROM MIXING
!   but non-P related precipitated species are not mixed, Jinyun Tang, Nov 30,2022
!   LL is the last layer tilled. 

  D2000: DO  L=NU_col(NY,NX),LL
    IF(DLYR_3D(3,L,NY,NX).GT.ZERO)THEN
      !TL: soil thickness in tillage zone
      !TI: fraction of current layer being tilled.
      !FI: fraction of current layer over the overall mixing zone
      !TX: fraction of current layer not tilled
      TL=AMIN1(DLYR_3D(3,L,NY,NX),DCORPZ-(CumSoilThickness_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)))
      FI1=TL/DCORPZ
      TI1=AMIN1(TL/DLYR_3D(3,L,NY,NX),1._r8)
      TX=AZMAX1(1.0_r8-TI1)
      
!     SoilFracAsMacP_vr(L,NY,NX)=XTillCorp_col(NY,NX)*SoilFracAsMacP_vr(L,NY,NX)

      ! Only the tilled, mixed portion of this layer changes pore domains.
      FracMacPoreMixed=TI1*FracTillCorp
      call MixSoluteMacpore(L,NY,NX,FracMacPoreMixed)

      DO  K=1,jcplx
        DO idom=idom_beg,idom_end
          TransferredDOM=FracMacPoreMixed*DOM_MacP_vr(idom,K,L,NY,NX)
          DOM_MicP_vr(idom,K,L,NY,NX)=DOM_MicP_vr(idom,K,L,NY,NX)+TransferredDOM
          DOM_MacP_vr(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)-TransferredDOM
        ENDDO  
      ENDDO
!
!     ADD STATE VARIABLES IN SURFACE RESIDUE INCORPORATED
!     WITHIN TILLAGE MIXING ZONE
!
      DO K=1,micpar%NumOfLitrCmplxs
        DO idom=idom_beg,idom_end
          DOM_MacP_vr(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)+FI1*TDOMH(idom,K)
        ENDDO
      ENDDO

      if(trcs_solml_vr(idg_O2,1,1,1)<0._r8)then
        print*,'tillage problem',trcs_solml_vr(idg_O2,1,1,1)
        stop
      endif

      ZNHU0_vr(L,NY,NX)=ZNHUX0
      ZNHUI_vr(L,NY,NX)=ZNHUXI
      ZNFN0_vr(L,NY,NX)=ZNFNX0
      ZNFNI_vr(L,NY,NX)=(TI1*ZNFNI_vr(L,NY,NX)+FracTillCorp*(FI1*TZNFNI-TI1*ZNFNI_vr(L,NY,NX))+TX*ZNFNI_vr(L,NY,NX)+FI1*TZNFNG)/FI1
      TZNFN2=TZNFN2+ZNFNI_vr(L,NY,NX)
    ENDIF
  ENDDO D2000

! nitrogen inhibitor
  ZNFN0_vr(0,NY,NX) = ZNFNX0
  ZNFNI_vr(0,NY,NX) = ZNFNI_vr(0,NY,NX)*XCORP0
  TZNFN2         = TZNFN2+TZNFNG
  TZNFNI         = TZNFNI+TZNFNG

  DO  L=NU_col(NY,NX),LL
    IF(TZNFN2.GT.ZERO)THEN
      ZNFNI_vr(L,NY,NX) = ZNFNI_vr(L,NY,NX)*TZNFNI/TZNFN2
      ZNFNI_vr(L,NY,NX) = ZNFNI_vr(L,NY,NX)+0.5_r8*(ZNFN0_vr(L,NY,NX)-ZNFNI_vr(L,NY,NX))
    ENDIF
  ENDDO

  call PrintInfo('end '//subname)
  end subroutine ApplyTillageMixing


!--------------------------------------------------------------------------------

  subroutine MixSoluteMacpore(L,NY,NX,FracMacPoreMixed)
  !
  ! Transfer macropore solutes to micropores in the tilled, mixed fraction.
  implicit none
  integer, intent(in) :: L,NY,NX
  real(r8), intent(in) :: FracMacPoreMixed !fraction of the whole layer affected [-]
  real(r8) :: TransferredSolute
  integer :: NTS,NTSA

  DO NTS=ids_beg,ids_end
    TransferredSolute=FracMacPoreMixed*trcs_soHml_vr(NTS,L,NY,NX)
    trcs_solml_vr(NTS,L,NY,NX)=trcs_solml_vr(NTS,L,NY,NX)+TransferredSolute
    trcs_soHml_vr(NTS,L,NY,NX)=trcs_soHml_vr(NTS,L,NY,NX)-TransferredSolute
  ENDDO

  if(salt_model)then
    DO NTSA=idsalt_beg,idsalt_end
      TransferredSolute=FracMacPoreMixed*trcSalt_soHml_vr(NTSA,L,NY,NX)
      trcSalt_solml_vr(NTSA,L,NY,NX)=trcSalt_solml_vr(NTSA,L,NY,NX)+TransferredSolute
      trcSalt_soHml_vr(NTSA,L,NY,NX)=trcSalt_soHml_vr(NTSA,L,NY,NX)-TransferredSolute
    ENDDO
  endif
  end subroutine MixSoluteMacpore
!--------------------------------------------------------------------------------

  subroutine MixIntensive1D(xdCorp,TI,FI,NU,LL,property1d)
  ! Mix material properties over the tilled volume, assuming the column's
  ! cross-sectional area is constant. A partially tilled layer retains its
  ! unmixed portion. No surface-litter property is incorporated here.
  implicit none
  real(r8), intent(in) :: xdCorp       !fraction of tilled material not mixed [-]
  real(r8), intent(in) :: TI(1:JZ)     !fraction of each layer inside the tillage zone [-]
  real(r8), intent(in) :: FI(1:JZ)     !each layer's share of total tilled thickness [-]
  integer, intent(in) :: NU,LL         !first and last layers in the tillage zone
  real(r8), intent(inout) :: property1d(:) !intensive property; element 1 is soil layer 1
  real(r8) :: MixedProperty,TotalWeight
  integer :: L

  MixedProperty=0._r8
  TotalWeight=0._r8
  do L=NU,LL
    MixedProperty=MixedProperty+FI(L)*property1d(L)
    TotalWeight=TotalWeight+FI(L)
  enddo
  if(TotalWeight<=0._r8)return
  MixedProperty=MixedProperty/TotalWeight

  do L=NU,LL
    property1d(L)=property1d(L)+(1._r8-xdCorp)*TI(L)*(MixedProperty-property1d(L))
  enddo
  end subroutine MixIntensive1D

!--------------------------------------------------------------------------------

  subroutine MixExtensive1D(xdCorp,TI,FI,NU,LL,mass1d,XCORP0)

  implicit none
  real(r8), intent(in) :: xdCorp    !fraction of layer not mixed  
  real(r8), intent(in) :: TI(1:JZ)  !fraction of layer L over the tillage zone
  real(r8), intent(in) :: FI(1:JZ)  !fraction of layer L being tilled
  integer , intent(in) :: NU,LL
  real(r8), dimension(:), intent(inout) :: mass1D
  real(r8), optional, intent(in) :: XCORP0
  real(r8)  :: tmass1D,tmass1DR,TX,CORP0,dCorp
  integer :: L,L1,L0

  L0=0  
  if(present(XCORP0))then
    L0         = 1
    CORP0      = 1.0_r8-XCORP0
    tmass1DR   = mass1D(L0)*CORP0
    mass1D(L0) = mass1D(L0)*XCORP0
  endif

  tmass1D = 0._r8
  dCorp   = 1.0_r8-xdCorp
  DO L=NU,LL    
    L1=L+L0
    tmass1D = tmass1D+TI(L)*mass1D(L1)*dCorp
  ENDDO

  D100: DO L=NU,LL
    L1=L+L0
    if(TI(L)<1._r8)then
      TX=1._r8-TI(L)
      mass1D(L1)=TI(L)*mass1D(L1)*xdCorp+FI(L)*tmass1D+TX*mass1D(L1)
    else
      mass1D(L1)=mass1D(L1)*xdCorp+FI(L)*tmass1D
    endif
    if(present(XCORP0))then
      mass1D(L1)=mass1D(L1)+FI(L)*tmass1DR
    endif
  ENDDO D100

  end subroutine MixExtensive1D

!--------------------------------------------------------------------------------

  subroutine MixExtensive2D(xdCorp,TI,FI,NU,LL,mass2d,XCORP0,name)

  implicit none
  real(r8), intent(in) :: xdCorp    !fraction of layer not mixed  
  real(r8), intent(in) :: TI(1:JZ)  !fraction of layer L over the tillage zone
  real(r8), intent(in) :: FI(1:JZ)  !fraction of layer L being tilled
  integer , intent(in) :: NU,LL

  real(r8), dimension(:,:), intent(inout) :: mass2D
  real(r8), optional, intent(in) :: XCORP0
  character(len=*), optional,intent(in) :: name
  real(r8)  :: tmass2DR(size(mass2D,1))
  real(r8)  :: tmass2D(size(mass2D,1))
  real(r8)  :: dCorp

  integer :: sz,ii,L,L0,L1
  real(r8) :: TX,CORP0   !fraction of layer L not tilled  

  if(present(name))print*,name
  sz=size(mass2D,1)
  L0=0
  if(present(XCORP0))then
    L0           = 1
    CORP0        = 1.0_r8-XCORP0
    tmass2DR     = mass2D(:,L0)*CORP0
    mass2D(:,L0) = mass2D(:,L0)*XCORP0
  endif

  tmass2D = 0._r8
  dCorp   = 1.0_r8-xdCorp
  DO L=NU,LL   
    L1=L+L0    
    DO ii=1,sz 
      tmass2D(ii) = tmass2D(ii)+TI(L)*mass2D(ii,L1)*dCorp
    enddo
  ENDDO

  D100: DO L=NU,LL
    L1=L+L0  
    do ii=1,sz
      if(TI(L)<1._r8)then
        TX=1._r8-TI(L)
        mass2D(ii,L1)=TI(L)*mass2D(ii,L1)*xdCorp+FI(L)*tmass2D(ii)+TX*mass2D(ii,L1)
      else
        mass2D(ii,L1)=mass2D(ii,L1)*xdCorp+FI(L)*tmass2D(ii)
      endif          
      if(present(XCORP0))then
        mass2D(ii,L1)=mass2D(ii,L1)+FI(L)*tmass2DR(ii)
      endif
    enddo
  ENDDO D100

  end subroutine MixExtensive2D

!--------------------------------------------------------------------------------

  subroutine MixExtensive3D(xdCorp,TI,FI,NU,LL,mass3d,XCORP0)

  implicit none
  real(r8), intent(in) :: xdCorp    !fraction of layer not mixed
  real(r8), intent(in) :: TI(1:JZ)  !fraction of layer L over the tillage zone
  real(r8), intent(in) :: FI(1:JZ)  !fraction of layer L being tilled
  integer , intent(in) :: NU,LL
  real(r8), dimension(:,:,:), intent(inout) :: mass3D
  real(r8), optional, intent(in) :: XCORP0
  real(r8)  :: tmass3D(size(mass3D,1),size(mass3D,2))
  real(r8)  :: tmass3DR(size(mass3D,1),size(mass3D,2))
  integer :: sz(3), ii, jj, L0,L1,L
  real(r8) :: TX,CORP0   !fraction of layer L not tilled  
  real(r8)  :: dCorp  

  sz=shape(mass3D)
  L0=0
  if(present(XCORP0))then
    L0             = 1
    CORP0          = 1.0_r8-XCORP0
    tmass3DR       = mass3D(:,:,L0)*CORP0
    mass3D(:,:,L0) = mass3D(:,:,L0)*XCORP0
  endif

  tmass3D = 0._r8
  dCorp   = 1.0_r8-xdCorp
  DO L=NU,LL  
    L1=L+L0  
    DO jj=1,sz(2) 
      DO ii=1,sz(1) 
        tmass3D(ii,jj) = tmass3D(ii,jj)+TI(L)*mass3D(ii,jj,L1)*dCorp
      enddo
    ENDDO
  ENDDO

  D100: DO L=NU,LL
    L1=L+L0  
    do jj=1,sz(2)
      do ii=1,sz(1)
          if(TI(L)<1._r8)then
            TX=1._r8-TI(L)
            mass3D(ii,jj,L1)=TI(L)*(mass3D(ii,jj,L1)*xdCorp+FI(L)*tmass3D(ii,jj))+TX*mass3D(ii,jj,L1)
          else
            mass3D(ii,jj,L1)=mass3D(ii,jj,L1)*xdCorp+FI(L)*tmass3D(ii,jj)
          endif          
          if(present(XCORP0))then
            mass3D(ii,jj,L1)=mass3D(ii,jj,L1)+FI(L1)*tmass3DR(ii,jj)
          endif
      enddo
    enddo
  ENDDO D100

  end subroutine MixExtensive3D

!--------------------------------------------------------------------------------

  subroutine MixExtensive4D(xdCorp,TI,FI,NU,LL,mass4d,XCORP0)

  implicit none
  real(r8), intent(in) :: xdCorp    !fraction of layer not mixed
  real(r8), intent(in) :: TI(1:JZ)  !fraction of layer L over the tillage zone
  real(r8), intent(in) :: FI(1:JZ)  !fraction of layer L being tilled
  integer,  intent(in) :: NU,LL
  real(r8), dimension(:,:,:,:), intent(inout) :: mass4D
  real(r8), optional, intent(in) :: XCORP0
  real(r8)  :: tmass4DR(size(mass4D,1),size(mass4D,2),size(mass4D,3))  
  real(r8)  :: tmass4D(size(mass4D,1),size(mass4D,2),size(mass4D,3))  
  real(r8) :: TX,CORP0   !fraction of layer L not tilled  
  integer :: sz(4), ii, jj,kk, L,L0,L1
  real(r8)  :: dCorp  
  
  sz=shape(mass4D)
  L0=0
  if(present(XCORP0))then
    L0               = 1
    CORP0            = 1.0_r8-XCORP0
    tmass4DR         = mass4D(:,:,:,L0)*CORP0
    mass4D(:,:,:,L0) = mass4D(:,:,:,L0)*XCORP0
  endif

  tmass4D = 0._r8
  dCorp   = 1.0_r8-xdCorp
  DO L=NU,LL  
    L1=L+L0  
    do kk=1,sz(3)
      DO jj=1,sz(2) 
        DO ii=1,sz(1) 
          tmass4D(ii,jj,kk) = tmass4D(ii,jj,kk)+TI(L)*mass4D(ii,jj,kk,L1)*dCorp
        enddo
      ENDDO
    enddo
  ENDDO

  D100: DO L=NU,LL
    L1=L+L0  
    do kk=1,sz(3)
      do jj=1,sz(2)
        do ii=1,sz(1)
          if(TI(L)<1._r8)then
            TX=1._r8-TI(L)
            mass4D(ii,jj,kk,L1)=TI(L)*mass4D(ii,jj,kk,L1)*xdCorp+FI(L)*tmass4D(ii,jj,kk)+TX*mass4D(ii,jj,kk,L1)
          else
            mass4D(ii,jj,kk,L1)=mass4D(ii,jj,kk,L1)*xdCorp+FI(L)*tmass4D(ii,jj,kk)
          endif     
          if(present(XCORP0))then
            mass4D(ii,jj,kk,L1)=mass4D(ii,jj,kk,L1)+FI(L)*tmass4DR(ii,jj,kk)
          endif
        enddo
      enddo
    enddo
  ENDDO D100

  end subroutine MixExtensive4D


 !--------------------------------------------------------------------------------
  subroutine DeriveTillageProf(NY,NX,DCORPZ,LL,FI,TI)
  !derive the mass allocation profile for tillage
  implicit none
  integer , intent(in) :: NY,NX
  real(r8), intent(in) :: DCORPZ
  integer , intent(out) :: LL
  real(r8), intent(out) :: FI(1:JZ)   !fraction of current layer over the overall mixing zone
  real(r8), intent(out) :: TI(1:JZ)   !fraction of layer L being tilled
  integer :: L
  logical :: tillTestL
  real(r8) :: TL  !fraction of layer L over the tillage zone

  LL=NU_col(NY,NX)
  D100: DO L=NU_col(NY,NX),NL_col(NY,NX)
    !if current layer is within the tillage depth
    tillTestL=CumSoilThickness_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX).LT.DCORPZ .AND. DLYR_3D(3,L,NY,NX).GT.ZERO
    if(.not.tillTestL)exit
    TL = AMIN1(DLYR_3D(3,L,NY,NX),DCORPZ-(CumSoilThickness_vr(L,NY,NX)-DLYR_3D(3,L,NY,NX)))
    FI(L) = TL/DCORPZ
    TI(L) = AMIN1(TL/DLYR_3D(3,L,NY,NX),1._r8)
    LL=L  
  ENDDO D100
  end subroutine DeriveTillageProf

end module TillageMixMod
