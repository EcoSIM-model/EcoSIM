module SoilDisturbMod
! Description:
! Remove soil and litter organic matter during fire or management, retaining
! accepted sensible energy and exporting the heat carried by withdrawn material.

! USES:
  use data_kind_mod, only : r8 => DAT_KIND_R8
  use abortutils  , only : endrun,iulog
  use minimathmod, only : safe_adb
  use EcoSiMParDataMod, only : micpar
  use MiniFuncMod,    only: gOC_to_m3_OM
  use EcoSIMCtrlMod,  only: plantOM4Heat
  use RootDataType,   only: RootMycoMassElm_vr
  use EcoSIMConfig,  only : NumMicbAFunGrupsPerCmplx,NumMicbHFunGrupsPerCmplx  
  use DebugToolMod
  use MicrobialDataType  
  use NitroPars
  use SOMDataType
  use ChemTranspDataType
  use FertilizerDataType
  use GridConsts
  use SoilPhysDataType
  use FlagDataType
  use SoilHeatDatatype
  use SoilWaterDataType
  use ClimForcDataType
  use EcosimConst
  use PlantMgmtDataType
  use EcoSimSumDataType
  use EcosimBGCFluxType
  use SoilBGCDataType
  use AqueChemDatatype
  use GridDataType
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: SOMRemovalByDisturbance
  contains
!------------------------------------------------------------------------------------------

  subroutine SOMRemovalByDisturbance(I,J,NY,NX)
  !
  ! Description:
  ! Remove organic pools for fire (22) or litter removal (21), transferring
  ! their elements and carried sensible heat to the corresponding ledgers.
  implicit none
  integer, intent(in) :: I,J,NY,NX
  character(len=*), parameter :: subname='SOMRemovalByDisturbance'

  integer :: L,K,M,N
  integer :: IFLGJ   !combustion flag, 0=yes, 1=no
  integer :: NLL,NGL,NTF,MID,ids,NE,idom
  real(r8) :: DC,DN,DP
  real(r8) :: DCORPC
  real(r8) :: FORGCX
  real(r8) :: HFLXD,rootC
  real(r8) :: EnergyBefore,EnergyRemaining ! accepted and retained sensible energy [MJ/cell]
  real(r8) :: ONX
  REAL(R8) :: OPX,OCA,OAH
  real(r8) :: ONL(4,1:jcplx),OPL(4,1:jcplx)
  real(r8) :: rmDOM(idom_beg:idom_end)
  real(r8) :: rmBiom(1:NumPlantChemElms )
  real(r8) :: OMelm(1:NumPlantChemElms)    !organic matter 
  real(r8) :: dOMelm(1:NumPlantChemElms)   !change of organic matter
  real(r8) :: DCORPC1

!     begin_execution
  call PrintInfo('beg '//subname)
  associate(                    &
   iprotein => micpar%iprotein, & ! input: protein kinetic-pool index [-]
   ilignin  => micpar%ilignin   & ! input: lignin kinetic-pool index [-]
  )
  IF(J.EQ.INT(SolarNoonHour_col(NY,NX)) .AND. &
    (iSoilDisturbType_col(I,NY,NX).EQ.itill_rmlitr .OR. iSoilDisturbType_col(I,NY,NX).EQ.itill_fire))THEN
    IF(iSoilDisturbType_col(I,NY,NX).EQ.itill_fire)THEN
      !fire
      iResetSoilProf_col(NY,NX) = itrue

      IFLGJ = 0
      NLL   = -1
      !identify burning depth
      D2945: DO L=0,NL_col(NY,NX)
        IF(L.EQ.0 .OR. L.GE.NUM_col(NY,NX))THEN
          IF(IFLGJ.EQ.1)THEN
            exit
          ELSEIF(THETW_vr(L,NY,NX).GT.VolMaxSoilMoist4Fire .OR. CSoilOrgM_vr(ielmc,L,NY,NX).LE.FORGC)THEN
            !combustion criteria not met: too much moisture, or too little organic matter
            IFLGJ=1
          ELSE
            NLL=L
          ENDIF
        ENDIF
      ENDDO D2945
    ELSE
      NLL=0
    ENDIF
    
    D2950: DO L=0,NLL
      IF(NLL.GE.0)THEN
        IF(iSoilDisturbType_col(I,NY,NX).EQ.itill_fire)THEN
          !fire
          IF(L.EQ.0)THEN
            FORGCX=0.0_r8
          ELSE
            FORGCX=FORGC
          ENDIF
          DCORPC=AMIN1(0.999_r8,DepzCorp_col(I,NY,NX))*(CSoilOrgM_vr(ielmc,L,NY,NX)-FORGCX) &
            /(AMAX1(CSoilOrgM_vr(ielmc,L,NY,NX),orgcden)-FORGCX)
        ELSE
          DCORPC=AMIN1(0.999_r8,DepzCorp_col(I,NY,NX))
        ENDIF
        !  QH2OLoss_lnds=QH2OLoss_lnds+DCORPC*VLWatMicP_vr(L,NY,NX)
        !  HeatOut_lnds=HeatOut_lnds+DCORPC*4.19*TKS_vr(L,NY,NX)*VLWatMicP_vr(L,NY,NX)
        !  VLWatMicP_vr(L,NY,NX)=VLWatMicP_vr(L,NY,NX)-DCORPC*VLWatMicP_vr(L,NY,NX)
        OMelm  = 0._r8
        dOMelm = 0._r8

        ONL(1:4,1:jcplx)=0._r8
        OPL(1:4,1:jcplx)=0._r8

        D2970: DO K=1,jcplx
          IF(L.NE.0.OR.(micpar%is_litter(K)))THEN
            !
            !     REMOVE heterotrophic MICROBIAL BIOMASS
            !
            D2960: DO N=1,NumMicbHFunGrupsPerCmplx
              DO NGL=JGniH(N),JGnfH(N)
                DO M=1,nlbiomcp
                  MID=micpar%get_micb_id(M,NGL)
                  DO NE=1,NumPlantChemElms 
                    rmBiom(NE)=DCORPC*mBiomeHeter_vr(NE,MID,K,L,NY,NX)
                  enddo

                  ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmn)
                  OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmp)

                  IF(micpar%is_litter(K))THEN
                    ONL(ilignin,K)=ONL(ilignin,K)+rmBiom(ielmn)-ONX
                    OPL(ilignin,K)=OPL(ilignin,K)+rmBiom(ielmp)-OPX
                  ELSE
                    ONL(iprotein,K)=ONL(iprotein,K)+rmBiom(ielmn)-ONX
                    OPL(iprotein,K)=OPL(iprotein,K)+rmBiom(ielmp)-OPX
                  ENDIF
                  DO  NE=1,NumPlantChemElms 
                    mBiomeHeter_vr(NE,MID,K,L,NY,NX) = mBiomeHeter_vr(NE,MID,K,L,NY,NX)-rmBiom(NE)
                    dOMelm(NE)                       = dOMelm(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
                  ENDDO
                  OMelm(ielmc)=OMelm(ielmc)+rmBiom(ielmc)
                  OMelm(ielmn)=OMelm(ielmn)+ONX
                  OMelm(ielmp)=OMelm(ielmp)+OPX
                enddo
              enddo
            ENDDO D2960
          ENDIF
        ENDDO D2970


        !          IF(L.NE.0.OR.(K.NE.3.AND.K.NE.4))THEN
        !
        !     REMOVE autotrophic MICROBIAL BIOMASS
        !
        DO  N=1,NumMicbAFunGrupsPerCmplx
          DO NGL=JGniA(N),JGnfA(N)
            DO M=1,nlbiomcp
              MID=micpar%get_micb_id(M,NGL)
              DO NE=1,NumPlantChemElms
                rmBiom(NE)=DCORPC*mBiomeAutor_vr(NE,MID,L,NY,NX)
              enddo

              ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmn)
              OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmp)

              ONL(ilignin,1)=ONL(ilignin,1)+rmBiom(ielmn)-ONX
              OPL(ilignin,1)=OPL(ilignin,1)+rmBiom(ielmp)-OPX

              DO NE=1,NumPlantChemElms
                mBiomeAutor_vr(NE,MID,L,NY,NX)=mBiomeAutor_vr(NE,MID,L,NY,NX)-rmBiom(NE)
                dOMelm(NE)=dOMelm(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
              enddo
              OMelm(ielmc)=OMelm(ielmc)+rmBiom(ielmc)
              OMelm(ielmn)=OMelm(ielmn)+ONX
              OMelm(ielmp)=OMelm(ielmp)+OPX              
            enddo
          enddo
        ENDDO
        !     ENDIF

        !
        !     REMOVE MICROBIAL RESIDUE
        !
        D2900: DO K=1,jcplx
          IF(L.NE.0 .OR. (micpar%is_litter(K)))THEN
            D2940: DO M=1,ndbiomcp
              DO NE=1,NumPlantChemElms
                rmBiom(NE)=DCORPC*OMBioResdu_vr(NE,M,K,L,NY,NX)
              enddo

              ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmn)
              OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmp)
              IF(micpar%is_litter(K))THEN
                ONL(ilignin,K)=ONL(ilignin,K)+rmBiom(ielmn)-ONX
                OPL(ilignin,K)=OPL(ilignin,K)+rmBiom(ielmp)-OPX
              ELSE
                ONL(iprotein,K)=ONL(iprotein,K)+rmBiom(ielmn)-ONX
                OPL(iprotein,K)=OPL(iprotein,K)+rmBiom(ielmp)-OPX
              ENDIF
              DO NE=1,NumPlantChemElms
                OMBioResdu_vr(NE,M,K,L,NY,NX) = OMBioResdu_vr(NE,M,K,L,NY,NX)-rmBiom(NE)
                dOMelm(NE)                    = dOMelm(NE)+OMBioResdu_vr(NE,M,K,L,NY,NX)
              enddo
              OMelm(ielmc)=OMelm(ielmc)+rmBiom(ielmc)
              OMelm(ielmn)=OMelm(ielmn)+ONX
              OMelm(ielmp)=OMelm(ielmp)+OPX                            
            ENDDO D2940
            !
            !     REMOVE DOC, DON, DOP
            !
            DO idom=idom_beg,idom_end
              rmDOM(idom)=DCORPC*DOM_MicP_vr(idom,K,L,NY,NX)
            enddo

            ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmn)
            OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmp)
            IF(micpar%is_litter(K))THEN
              ONL(ilignin,K)=ONL(ilignin,K)+rmDOM(ielmn)-ONX
              OPL(ilignin,K)=OPL(ilignin,K)+rmDOM(ielmp)-OPX
            ELSE
              ONL(iprotein,K)=ONL(iprotein,K)+rmDOM(ielmn)-ONX
              OPL(iprotein,K)=OPL(iprotein,K)+rmDOM(ielmp)-OPX
            ENDIF
            DO idom=idom_beg,idom_end
              DOM_MicP_vr(idom,K,L,NY,NX)=DOM_MicP_vr(idom,K,L,NY,NX)-rmDOM(idom)
            enddo

            OMelm(ielmc)=OMelm(ielmc)+rmDOM(idom_doc)+rmDOM(idom_acetate)            
            OMelm(ielmn)=OMelm(ielmn)+ONX
            OMelm(ielmp)=OMelm(ielmp)+OPX

            DO idom=idom_beg,idom_end
              rmDOM(idom)=DCORPC*DOM_MacP_vr(idom,K,L,NY,NX)
            enddo

            ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmn)
            OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmp)
            IF(micpar%is_litter(K))THEN
              ONL(ilignin,K)=ONL(ilignin,K)+rmDOM(ielmn)-ONX
              OPL(ilignin,K)=OPL(ilignin,K)+rmDOM(ielmp)-OPX
            ELSE
              ONL(iprotein,K)=ONL(iprotein,K)+rmDOM(ielmn)-ONX
              OPL(iprotein,K)=OPL(iprotein,K)+rmDOM(ielmp)-OPX
            ENDIF
            DO idom=idom_beg,idom_end
              DOM_MacP_vr(idom,K,L,NY,NX)=DOM_MacP_vr(idom,K,L,NY,NX)-rmDOM(idom)
            enddo

            OMelm(ielmc)=OMelm(ielmc)+rmDOM(idom_doc)+rmDOM(idom_acetate)            
            OMelm(ielmn)=OMelm(ielmn)+ONX
            OMelm(ielmp)=OMelm(ielmp)+OPX
            !
            !     REMOVE ADSORBED OM
            !
            DO idom=idom_beg,idom_end
              rmDOM(idom)=DCORPC*SorbedOM_vr(idom,K,L,NY,NX)
            enddo

            ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmn)
            OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmDOM(ielmp)
            IF(micpar%is_litter(K))THEN
              ONL(ilignin,K)=ONL(ilignin,K)+rmDOM(ielmn)-ONX
              OPL(ilignin,K)=OPL(ilignin,K)+rmDOM(ielmp)-OPX
            ELSE
              ONL(iprotein,K)=ONL(iprotein,K)+rmDOM(ielmn)-ONX
              OPL(iprotein,K)=OPL(iprotein,K)+rmDOM(ielmp)-OPX
            ENDIF
            DO idom=idom_beg,idom_end
              SorbedOM_vr(idom,K,L,NY,NX)=SorbedOM_vr(idom,K,L,NY,NX)-rmDOM(idom)
            enddo

            dOMelm(ielmc)=dOMelm(ielmc)+DOM_MicP_vr(idom_doc,K,L,NY,NX)+DOM_MacP_vr(idom_doc,K,L,NY,NX)+SorbedOM_vr(ielmc,K,L,NY,NX) &
              +DOM_MicP_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)
            dOMelm(ielmn)=dOMelm(ielmn)+DOM_MicP_vr(idom_don,K,L,NY,NX)+DOM_MacP_vr(idom_don,K,L,NY,NX)+SorbedOM_vr(ielmn,K,L,NY,NX)
            dOMelm(ielmp)=dOMelm(ielmp)+DOM_MicP_vr(idom_dop,K,L,NY,NX)+DOM_MacP_vr(idom_dop,K,L,NY,NX)+SorbedOM_vr(ielmp,K,L,NY,NX)

            OMelm(ielmc)=OMelm(ielmc)+rmDOM(idom_doc)+rmDOM(idom_acetate)
            OMelm(ielmn)=OMelm(ielmn)+ONX
            OMelm(ielmp)=OMelm(ielmp)+OPX!
            !     REMOVE RESIDUE
            !
            D2930: DO M=1,jsken
              do NE=1,NumPlantChemElms
                rmBiom(NE)=DCORPC*SolidOM_vr(NE,M,K,L,NY,NX)
              enddo

              OCA=DCORPC*SolidOMAct_vr(M,K,L,NY,NX)
              ONX=EFIRE(1,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmn)
              OPX=EFIRE(2,iSoilDisturbType_col(I,NY,NX))*rmBiom(ielmp)
              ONL(M,K)=ONL(M,K)+rmBiom(ielmn)-ONX
              OPL(M,K)=OPL(M,K)+rmBiom(ielmp)-OPX

              do NE=1,NumPlantChemElms
                SolidOM_vr(NE,M,K,L,NY,NX) = SolidOM_vr(NE,M,K,L,NY,NX)-rmBiom(NE)
                dOMelm(NE)                 = dOMelm(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
              enddo
                            
              SolidOMAct_vr(M,K,L,NY,NX)=SolidOMAct_vr(M,K,L,NY,NX)-OCA              
              OMelm(ielmc)=OMelm(ielmc)+rmBiom(ielmc)
              OMelm(ielmn)=OMelm(ielmn)+ONX
              OMelm(ielmp)=OMelm(ielmp)+OPX      
            ENDDO D2930
          ENDIF
        ENDDO D2900
        !
        !     ADD UNBURNED N,P TO ORG N, ORG P
        !
        D2905: DO K=1,jcplx
          DO  M=1,jsken
            SolidOM_vr(ielmn,M,K,L,NY,NX) = SolidOM_vr(ielmn,M,K,L,NY,NX)+ONL(M,K)
            SolidOM_vr(ielmp,M,K,L,NY,NX) = SolidOM_vr(ielmp,M,K,L,NY,NX)+OPL(M,K)
            DN                            = DN+ONL(M,K)
            DP                            = DP+OPL(M,K)
          enddo
        ENDDO D2905
        !
        !     REMOVE FERTILIZER IN RESIDUE
        !
        IF(iSoilDisturbType_col(I,NY,NX).EQ.itill_rmlitr)THEN
          OMelm(ielmn)=OMelm(ielmn)+DCORPC*(trcs_solml_vr(ids_NH4,L,NY,NX)+trcs_solml_vr(idg_NH3,L,NY,NX) &
            +trcs_solml_vr(ids_NO3,L,NY,NX)+trcs_solml_vr(ids_NO2,L,NY,NX))
          OMelm(ielmp)=OMelm(ielmp)+DCORPC*(trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H2PO4,L,NY,NX))
          ! Record the elemental mass exported from the molar pools below.
          OMelm(ielmn)=OMelm(ielmn)+DCORPC*natomw*(trcx_solml_vr(idx_NH4,L,NY,NX) &
            +SUM(FertN_mole_soil_vr(ifertn_beg:ifertn_end,L,NY,NX)))
          OMelm(ielmp)=OMelm(ielmp)+DCORPC*patomw*(trcp_saltpml_vr(idsp_AlPO4,L,NY,NX) &
            +trcp_saltpml_vr(idsp_FePO4,L,NY,NX)+trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX) &
            +2._r8*trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)+3._r8*trcp_saltpml_vr(idsp_Apatite,L,NY,NX))
          
          DCORPC1=1.0_r8-DCORPC

          trcs_solml_vr(idg_NH3,L,NY,NX)=DCORPC1*trcs_solml_vr(idg_NH3,L,NY,NX)
          do ids=ids_nut_beg,ids_nuts_end
            trcs_solml_vr(ids,L,NY,NX)=DCORPC1*trcs_solml_vr(ids,L,NY,NX)
          enddo
          trcx_solml_vr(idx_NH4,L,NY,NX)         = DCORPC1*trcx_solml_vr(idx_NH4,L,NY,NX)
          trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)    = DCORPC1*trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)
          trcp_saltpml_vr(idsp_FePO4,L,NY,NX)    = DCORPC1*trcp_saltpml_vr(idsp_FePO4,L,NY,NX)
          trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)   = DCORPC1*trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)
          trcp_saltpml_vr(idsp_Apatite,L,NY,NX)       = DCORPC1*trcp_saltpml_vr(idsp_Apatite,L,NY,NX)
          trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX) = DCORPC1*trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)

          DO NTF=ifertn_beg,ifertn_end
            FertN_mole_soil_vr(NTF,L,NY,NX)=DCORPC1*FertN_mole_soil_vr(NTF,L,NY,NX)
          ENDDO
        ENDIF
        SoilOrgM_vr(1:NumPlantChemElms,L,NY,NX)=dOMelm(1:NumPlantChemElms)
        
        ! Export sensible heat with the actual organic C withdrawn by this
        ! event, using the same capacity law as the remaining organic matter.
        ! ORGCX is an earlier hourly snapshot and includes other C changes.
        ! This accounts for carried sensible heat, not combustion heat release.
        EnergyBefore = VHeatCapacity_vr(L,NY,NX)*TKS_vr(L,NY,NX)
        HFLXD        = cpo*gOC_to_m3_OM(OMelm(ielmc))*TKS_vr(L,NY,NX)
        HeatOut_lnds = HeatOut_lnds+HFLXD
        ! This event follows hydrology. Book its export directly in the column
        ! ledger; a layer heat-source flux would be reapplied by thermal updates.
        HeatSource_col(NY,NX)=HeatSource_col(NY,NX)-HFLXD
        EnergyRemaining=EnergyBefore-HFLXD
        IF(L.EQ.0)then
          VHeatCapacity_vr(0,NY,NX)=cpw*VLWatMicP_vr(0,NY,NX)+cpi*VLiceMicP_vr(0,NY,NX)+cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,L,NY,NX))
        else
          VHeatCapSolidSoil_vr(L,NY,NX) = VHeatCapMineralSoil_vr(L,NY,NX)+cpo*gOC_to_m3_OM(SoilOrgM_vr(ielmc,L,NY,NX))
          VHeatCapacity_vr(L,NY,NX)     = VHeatCapSolidSoil_vr(L,NY,NX)+cpw*(VLWatMicP_vr(L,NY,NX)+VLWatMacP_vr(L,NY,NX)) &
            +cpi*(VLiceMicP_vr(L,NY,NX)+VLiceMacP_vr(L,NY,NX))
          ! Roots are separate from the soil organic pools removed above.
          if(plantOM4Heat .and. VHeatCapacity_vr(L,NY,NX)>0._r8)then
            rootC=sum(RootMycoMassElm_vr(ielmc,:,L,NY,NX))
            VHeatCapacity_vr(L,NY,NX)=VHeatCapacity_vr(L,NY,NX)+cpo*gOC_to_m3_OM(rootC)
          endif
        endif
        ! Updating the summarized organic capacity supplies no extra heat.
        ! Retain the accepted inventory minus the explicitly exported material
        ! enthalpy, including organic-pool changes since its earlier estimate.
        if(VHeatCapacity_vr(L,NY,NX)>0._r8)then
          TKS_vr(L,NY,NX)=EnergyRemaining/VHeatCapacity_vr(L,NY,NX)
          TCS_vr(L,NY,NX)=TKS_vr(L,NY,NX)-TC2K
        elseif(ABS(EnergyRemaining)>1.e-12_r8*MAX(1._r8,ABS(EnergyBefore)))then
          write(iulog,*)subname,': nonzero retained energy with zero capacity',I,J,L,NY,NX,EnergyRemaining
          call endrun('Disturbed soil/litter has no capacity for retained energy',__LINE__)
        endif

        IF(iSoilDisturbType_col(I,NY,NX).EQ.itill_rmlitr)THEN
          DO NE=1,NumPlantChemElms
            TOMOU_lnds(NE)=TOMOU_lnds(NE)+OMelm(NE)
          ENDDO

          HydroSufDOCFlx_col(NY,NX) = HydroSufDOCFlx_col(NY,NX)+OMelm(ielmc)
          HydroSufDONFlx_col(NY,NX) = HydroSufDONFlx_col(NY,NX)+OMelm(ielmn)
          HydroSufDOPFlx_col(NY,NX) = HydroSufDOPFlx_col(NY,NX)+OMelm(ielmp)
          Eco_NBP_CumYr_col(NY,NX)  = Eco_NBP_CumYr_col(NY,NX)-OMelm(ielmc)

        ELSEIF(iSoilDisturbType_col(I,NY,NX).EQ.itill_fire)THEN
          SurfGas_lnd(idg_CO2) = SurfGas_lnd(idg_CO2)-OMelm(ielmc)
          SurfGas_lnd(idg_O2)  = SurfGas_lnd(idg_O2)+2.667_r8*OMelm(ielmc)
          OXYGOU               = OXYGOU+2.667_r8*OMelm(ielmc)

          TOMOU_lnds(ielmn)          = TOMOU_lnds(ielmn)+OMelm(ielmn)
          TOMOU_lnds(ielmp)          = TOMOU_lnds(ielmp)+OMelm(ielmp)
          CO2byFire_CumYr_col(NY,NX) = CO2byFire_CumYr_col(NY,NX)-(1.0_r8-FrcAsCH4byFire)*OMelm(ielmc)
          CH4byFire_CumYr_col(NY,NX) = CH4byFire_CumYr_col(NY,NX)-FrcAsCH4byFire*OMelm(ielmc)
          O2byFire_CumYr_col(NY,NX)  = O2byFire_CumYr_col(NY,NX)+(1.0_r8-FrcAsCH4byFire)*2.667_r8*OMelm(ielmc)
          NH3byFire_CumYr_col(NY,NX) = NH3byFire_CumYr_col(NY,NX)-OMelm(ielmn)
          N2ObyFire_CumYr_col(NY,NX) = N2ObyFire_CumYr_col(NY,NX)-0.0_r8
          PO4byFire_CumYr_col(NY,NX) = PO4byFire_CumYr_col(NY,NX)-OMelm(ielmp)
          Eco_NBP_CumYr_col(NY,NX)   = Eco_NBP_CumYr_col(NY,NX)-OMelm(ielmc)
        ENDIF
      ENDIF
    ENDDO D2950
  ENDIF
  end associate
  call PrintInfo('end '//subname)
  end subroutine SOMRemovalByDisturbance
!------------------------------------------------------------------------------------------


end module SoilDisturbMod
