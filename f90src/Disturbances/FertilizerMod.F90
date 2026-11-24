module FertilizerMod
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use EcoSIMCtrlMod, only: lverb
  use DebugToolMod,  only: PrintInfo
  use abortutils,    only: iulog,endrun
  use EcoSIMConfig,  only : NumMicbAFunGrupsPerCmplx,NumMicbHFunGrupsPerCmplx
  use SurfLitterDataType
  use EcoSimSumDataType
  use FertilizerDataType
  use ClimForcDataType
  use SoilBGCDataType
  use FlagDataType
  use GridDataType
  use EcosimBGCFluxType
  use SOMDataType
  use EcoSiMParDataMod
  use EcoSIMCtrlDataType
  use SoilPropertyDataType
  use MicrobialDataType
  use AqueChemDatatype
  use MiniMathMod
  use EcosimConst
  use fileUtil
implicit none
  private

  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: ApplyFertilizerAtNoon

  contains
!------------------------------------------------------------------------------------------

  subroutine ApplyFertilizerAtNoon(I,J,NHW,NHE,NVN,NVS)
!
  implicit none
  integer, intent(in) :: I,J,NHW,NHE,NVN,NVS

  integer :: NX,NY
  real(r8) :: OFC(2),OFN(2),OFP(2)
  integer :: LFDPTH
  logical :: UreaLayers(0:JZ),NitrifierLayers(0:JZ)
  character(len=*), parameter :: subname='ApplyFertilizerAtNoon'
!     begin_execution

  call PrintInfo('beg '//subname)

  D8990: DO NX=NHW,NHE
    D8995: DO NY=NVN,NVS
      IF(J.EQ.INT(SolarNoonHour_col(NY,NX)))THEN

        call ApplyMineralFertilizer(I,J,NY,NX,LFDPTH,OFC,OFN,OFP,UreaLayers,NitrifierLayers)
!
!     SOIL LAYER NUMBER IN WHICH PLANT OR ANIMAL RESIDUES ARE APPLIED
!
       if(lverb)write(iulog,*) 'ApplyManure'
        call ApplyManure(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
!
!     FERTILIZER UREA, NITRIFICATION INHIBITORS
        if(lverb)write(iulog,*)'ApplyUreaNitrifierInhibitor'
        call ApplyUreaNitrifierInhibitor(I,J,NY,NX,UreaLayers,NitrifierLayers)

      ENDIF

    ENDDO D8995
  ENDDO D8990
  call PrintInfo('end '//subname)
  end subroutine ApplyFertilizerAtNoon
!------------------------------------------------------------------------------------------

  subroutine ApplyUreaNitrifierInhibitor(I,J,NY,NX,UreaLayers,NitrifierLayers)

  implicit none
  integer, intent(in) :: I,J,NY,NX
  logical, intent(in) :: UreaLayers(0:JZ),NitrifierLayers(0:JZ)
  character(len=*), parameter :: subname='ApplyUreaNitrifierInhibitor'
  integer :: L
!     begin_execution
  call PrintInfo('beg '//subname)
!
!     IYTYP=fertilizer release type from fertilizer input file
!     FERT=fertilizer type from fertilizer input file
!     iUreaHydInhibitorType_col=urea hydrolysis inhibitor type (1=no,2=yes)
!     ZNHU0,ZNHUI=initial,current urea hydrolysis inhibition activity
!     ZNFN0,ZNFNI=initial,current nitrification inhibition activity
!     Inhibitors are dimensionless activity factors with a maximum of 1,
!     not concentrations or conserved masses. Resetting receiving layers
!     to 1 and other layers to 0 for each treatment is an intentional
!     compromise; overlapping applications are not accumulated.
!     TODO: Revisit treatment overlap across dates and depths, including
!     the column-wide urea-inhibitor type, if a more detailed representation
!     is introduced. Preserving old activity would also require an explicit
!     assumption about how successive treatments interact.
! urea application
  IF(FERT(ifert_N_urea,I,NY,NX).GT.0._r8 .OR. FERT(ifert_N_urea_band,I,NY,NX).GT.0._r8)THEN  
    IF(IYTYP(iamendtyp_fert,I,NY,NX).EQ.0)THEN
      !default hydrolysis rate
      iUreaHydInhibitorType_col(NY,NX)=0
    ELSEIF(IYTYP(iamendtyp_fert,I,NY,NX).EQ.1 .OR. IYTYP(iamendtyp_fert,I,NY,NX).EQ.3)THEN
      !reduced hydrolysis rate
      iUreaHydInhibitorType_col(NY,NX)=1
    ELSE
      !urea hydrolysis is on    
      !more reduced hydrolysis rate
      iUreaHydInhibitorType_col(NY,NX)=2
    ENDIF

    ! Activity is intensive: initialize each receiving layer at full activity.
    ! The masks retain mineral placement even if manure changes LFDPTH.
    !urea inhibitor is on
    D9964: DO L=0,NL_col(NY,NX)
      IF(UreaLayers(L))THEN
        ZNHU0_vr(L,NY,NX)=1.0_r8  !inital inhibitor activity
        ZNHUI_vr(L,NY,NX)=1.0_r8  !current inhibitor activity
      ELSE
        ZNHU0_vr(L,NY,NX)=0._r8
        ZNHUI_vr(L,NY,NX)=0._r8
      ENDIF
    ENDDO D9964
  ENDIF

  IF(IYTYP(iamendtyp_fert,I,NY,NX).EQ.3 .OR. IYTYP(iamendtyp_fert,I,NY,NX).EQ.4)THEN
    !nitrification inhibitor is on
    D9965: DO L=0,NL_col(NY,NX)
      IF(NitrifierLayers(L))THEN
        ZNFN0_vr(L,NY,NX)=1.0_r8  !initial inhibitor activity
        ZNFNI_vr(L,NY,NX)=1.0_r8  !current inhibitor activity
      ELSE
        ZNFN0_vr(L,NY,NX)=0._r8
        ZNFNI_vr(L,NY,NX)=0._r8
      ENDIF
    ENDDO D9965
  ENDIF
  call PrintInfo('end '//subname)
  end subroutine ApplyUreaNitrifierInhibitor
!------------------------------------------------------------------------------------------

  subroutine ApplyManure(I,J,NY,NX,LFDPTH,OFC,OFN,OFP)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  integer, intent(inout) :: LFDPTH
  real(r8), intent(in) :: OFC(2),OFN(2),OFP(2)
  character(len=*), parameter :: subname='ApplyManure'
  real(r8) :: CORGCX
  real(r8) :: FDPTHM
  REAL(R8) :: OSCI,OSNI,OSPI
  REAL(R8) :: OSCX,OSNX,OSPX
  REAL(R8) :: OMC1,OMN1,OMP1
  real(r8) :: OQC1,OQN1,OQP1
  real(r8) :: OSC1,OSN1,OSP1
  REAL(R8) :: RNT,RPT
  integer  :: L,K,M,N,NN,NGL,amendment
  integer :: numGuildsH,numGuildsA
  real(r8) :: OMC1g,OMN1g,OMP1g,MicrobialShareFactor
!     begin_execution
  call PrintInfo('beg '//subname)

  associate(                           &
    k_fine_comp => micpar%k_fine_comp, &
    k_manure    => micpar%k_manure   , &
    ilignin     => micpar%ilignin    , &
    icellulos   => micpar%icellulos  , &
    icarbhyro   => micpar%icarbhyro  , &
    iprotein    => micpar%iprotein   , &
    OMCI        => micpar%OMCI      , &
    OMCF        => micpar%OMCF      , &
    OMCA        => micpar%OMCA      , &
    rNCOMC_ave  => micpar%rNCOMC_ave, &
    rPCOMC_ave  => micpar%rPCOMC_ave, &
    CNOFC       => micpar%CNOFC     , &
    CPOFC       => micpar%CPOFC       &
  )
!     LFDPTH=layer number
!
  IF(OFC(1)+OFC(2).GT.0._r8)THEN
    FDPTHM=FDPTH(I,NY,NX)+CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)
    IF(FDPTHM.LE.0._r8)THEN
      LFDPTH=0
    ELSE
      ! Organic-only records also require a destination in the active profile.
      LFDPTH=-1
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        IF(CumDepz2LayBottom_vr(L,NY,NX).GE.FDPTHM)THEN
          LFDPTH=L
          exit
        ENDIF
      ENDDO
      if(LFDPTH<0)call endrun('Organic amendment depth exceeds the soil profile in '//mod_filename,__LINE__)
    ENDIF
!
!     ALLOCATION OF PLANT RESIDUE APPLICATION TO
!     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     CFOSC=fraction of litter allocated to protein(1)
!     soluble CH2O(2), cellulose(3) and lignin(4)
!     ITYPE=litter type entered in fertilizer input file
!
!     MAIZE
!
    IF(IYTYP(iAmendtyp_plantRes,I,NY,NX).EQ.iPlantRes_maize)THEN
      CFOSC_vr(iprotein,k_fine_comp,LFDPTH,NY,NX)  = 0.080_r8
      CFOSC_vr(icarbhyro,k_fine_comp,LFDPTH,NY,NX) = 0.245_r8
      CFOSC_vr(icellulos,k_fine_comp,LFDPTH,NY,NX) = 0.613_r8
      CFOSC_vr(ilignin,k_fine_comp,LFDPTH,NY,NX)   = 0.062_r8
!
!     WHEAT
!
    ELSEIF(IYTYP(iAmendtyp_plantRes,I,NY,NX).EQ.iPlantRes_wheat)THEN
      CFOSC_vr(iprotein,k_fine_comp,LFDPTH,NY,NX)  = 0.125_r8
      CFOSC_vr(icarbhyro,k_fine_comp,LFDPTH,NY,NX) = 0.171_r8
      CFOSC_vr(icellulos,k_fine_comp,LFDPTH,NY,NX) = 0.560_r8
      CFOSC_vr(ilignin,k_fine_comp,LFDPTH,NY,NX)   = 0.144_r8
!
!     SOYBEAN
!
    ELSEIF(IYTYP(iAmendtyp_plantRes,I,NY,NX).EQ.iPlantRes_soybean)THEN
      CFOSC_vr(iprotein,k_fine_comp,LFDPTH,NY,NX)  = 0.138_r8
      CFOSC_vr(icarbhyro,k_fine_comp,LFDPTH,NY,NX) = 0.426_r8
      CFOSC_vr(icellulos,k_fine_comp,LFDPTH,NY,NX) = 0.316_r8
      CFOSC_vr(ilignin,k_fine_comp,LFDPTH,NY,NX)   = 0.120_r8
!
!     OLD STRAW
!
    ELSEIF(IYTYP(iAmendtyp_plantRes,I,NY,NX).EQ.iPlantRes_oldStraw)THEN
      CFOSC_vr(iprotein,k_fine_comp,LFDPTH,NY,NX)  = 0.075_r8
      CFOSC_vr(icarbhyro,k_fine_comp,LFDPTH,NY,NX) = 0.125_r8
      CFOSC_vr(icellulos,k_fine_comp,LFDPTH,NY,NX) = 0.550_r8
      CFOSC_vr(ilignin,k_fine_comp,LFDPTH,NY,NX)   = 0.250_r8
!
!     STRAW
!
    ELSEIF(IYTYP(iAmendtyp_plantRes,I,NY,NX).EQ.iPlantRes_Straw)THEN
      CFOSC_vr(iprotein,k_fine_comp,LFDPTH,NY,NX)  = 0.036_r8
      CFOSC_vr(icarbhyro,k_fine_comp,LFDPTH,NY,NX) = 0.044_r8
      CFOSC_vr(icellulos,k_fine_comp,LFDPTH,NY,NX) = 0.767_r8
      CFOSC_vr(ilignin,k_fine_comp,LFDPTH,NY,NX)   = 0.153_r8
!
!     COMPOST
!
    ELSEIF(IYTYP(iAmendtyp_plantRes,I,NY,NX).EQ.iPlantRes_compost)THEN
      CFOSC_vr(iprotein,k_fine_comp,LFDPTH,NY,NX)  = 0.143_r8
      CFOSC_vr(icarbhyro,k_fine_comp,LFDPTH,NY,NX) = 0.015_r8
      CFOSC_vr(icellulos,k_fine_comp,LFDPTH,NY,NX) = 0.640_r8
      CFOSC_vr(ilignin,k_fine_comp,LFDPTH,NY,NX)   = 0.202_r8
!
!     GREEN MANURE
!
    ELSEIF(IYTYP(iAmendtyp_plantRes,I,NY,NX).EQ.iPlantRes_GreeManure)THEN
      CFOSC_vr(iprotein,k_fine_comp,LFDPTH,NY,NX)  = 0.202_r8
      CFOSC_vr(icarbhyro,k_fine_comp,LFDPTH,NY,NX) = 0.013_r8
      CFOSC_vr(icellulos,k_fine_comp,LFDPTH,NY,NX) = 0.560_r8
      CFOSC_vr(ilignin,k_fine_comp,LFDPTH,NY,NX)   = 0.225_r8
!
!     SIMPLE SUBSTRATE
!
    ELSEIF(IYTYP(iAmendtyp_plantRes,I,NY,NX).EQ.iPlantRes_simple)THEN
      CFOSC_vr(iprotein,k_fine_comp,LFDPTH,NY,NX)  = 0.000_r8
      CFOSC_vr(icarbhyro,k_fine_comp,LFDPTH,NY,NX) = 1.000_r8
      CFOSC_vr(icellulos,k_fine_comp,LFDPTH,NY,NX) = 0.000_r8
      CFOSC_vr(ilignin,k_fine_comp,LFDPTH,NY,NX)   = 0.000_r8
    ELSE
      CFOSC_vr(iprotein,k_fine_comp,LFDPTH,NY,NX)  = 0.075_r8
      CFOSC_vr(icarbhyro,k_fine_comp,LFDPTH,NY,NX) = 0.125_r8
      CFOSC_vr(icellulos,k_fine_comp,LFDPTH,NY,NX) = 0.550_r8
      CFOSC_vr(ilignin,k_fine_comp,LFDPTH,NY,NX)   = 0.250_r8
    ENDIF
!
!     ALLOCATION OF ANIMAL MANURE APPLICATION TO
!     RESIDUE PROTEIN, CH2O, CELLULOSE, LIGNIN
!
!     RUMINANT
!
    IF(IYTYP(iAmendtyp_Manure,I,NY,NX).EQ.imanure_ruminant)THEN
      CFOSC_vr(iprotein,k_manure,LFDPTH,NY,NX)  = 0.036_r8
      CFOSC_vr(icarbhyro,k_manure,LFDPTH,NY,NX) = 0.044_r8
      CFOSC_vr(icellulos,k_manure,LFDPTH,NY,NX) = 0.630_r8
      CFOSC_vr(ilignin,k_manure,LFDPTH,NY,NX)   = 0.290_r8
!
!     NON-RUMINANT
!
    ELSEIF(IYTYP(iAmendtyp_Manure,I,NY,NX).EQ.imanure_nonruminant)THEN
      CFOSC_vr(iprotein,k_manure,LFDPTH,NY,NX)  = 0.138_r8
      CFOSC_vr(icarbhyro,k_manure,LFDPTH,NY,NX) = 0.401_r8
      CFOSC_vr(icellulos,k_manure,LFDPTH,NY,NX) = 0.316_r8
      CFOSC_vr(ilignin,k_manure,LFDPTH,NY,NX)   = 0.145_r8
!
!     GRAZING
!
    ELSEIF(IYTYP(iAmendtyp_Manure,I,NY,NX).EQ.imanure_grazing)THEN
      CFOSC_vr(iprotein,k_manure,LFDPTH,NY,NX)  = 0.036_r8
      CFOSC_vr(icarbhyro,k_manure,LFDPTH,NY,NX) = 0.044_r8
      CFOSC_vr(icellulos,k_manure,LFDPTH,NY,NX) = 0.630_r8
      CFOSC_vr(ilignin,k_manure,LFDPTH,NY,NX)   = 0.290_r8
!
!     OTHER
!
    ELSE
      CFOSC_vr(iprotein,k_manure,LFDPTH,NY,NX)  = 0.138_r8
      CFOSC_vr(icarbhyro,k_manure,LFDPTH,NY,NX) = 0.401_r8
      CFOSC_vr(icellulos,k_manure,LFDPTH,NY,NX) = 0.316_r8
      CFOSC_vr(ilignin,k_manure,LFDPTH,NY,NX)   = 0.145_r8
    ENDIF
!
!     DISTRIBUTE RESIDUE APPLICATION AMONG COMPONENTS OF RESIDUE complex
!
!     OFC,OFN,OFP=litter C,N,P application from fertilizer file
!

    ! Reserve the represented autotroph shares together with each heterotroph
    ! allocation, before debiting the common amendment C/N/P donor budgets.
    MicrobialShareFactor=1._r8
    DO NN=1,NumMicbAFunGrupsPerCmplx
      if(JGnfA(NN)>=JGniA(NN))MicrobialShareFactor=MicrobialShareFactor+OMCA(NN)
    ENDDO
    D2965: DO amendment=1,2
      ! Input slots are plant residue and animal manure; their destination
      ! complexes are fine litter and manure, respectively.
      IF(amendment==1)THEN
        K=k_fine_comp
      ELSE
        K=k_manure
      ENDIF
      OSCI=OFC(amendment)*AREA_3D(3,LFDPTH,NY,NX)
      OSNI=OFN(amendment)*AREA_3D(3,LFDPTH,NY,NX)
      OSPI=OFP(amendment)*AREA_3D(3,LFDPTH,NY,NX)
      IF(VLSoilMicPMass_vr(LFDPTH,NY,NX).GT.ZEROS(NY,NX))THEN
        CORGCX=OSCI/VLSoilMicPMass_vr(LFDPTH,NY,NX)
      ELSE
        CORGCX=orgcden
      ENDIF
      OSCX=0._r8
      OSNX=0._r8
      OSPX=0._r8
!
!     BIOMASSES OF MICROBIAL POPULATIONS IN RESIDUE
!
!     OMC,OMN,OMP=microbial biomass in litter application
!     OMCI=microbial biomass content in litter
!     OMCF,OMCA=hetero,autotrophic biomass composition in litter
!
      D2960: DO N=1,NumMicbHFunGrupsPerCmplx
        numGuildsH=JGnfH(N)-JGniH(N)+1
        ! Unrepresented groups leave their amendment share in the remainder.
        if(numGuildsH<=0)cycle D2960
        D2961: DO M=1,nlbiomcp
          OMC1=AZMAX1(AMIN1(OSCI*OMCI(M,K)*OMCF(N),(OSCI-OSCX)/MicrobialShareFactor))
          OMN1=AZMAX1(AMIN1(OMC1*rNCOMC_ave(M,N,K),(OSNI-OSNX)/MicrobialShareFactor))
          OMP1=AZMAX1(AMIN1(OMC1*rPCOMC_ave(M,N,K),(OSPI-OSPX)/MicrobialShareFactor))
          DO NGL=JGniH(N),JGnfH(N)
            associate(MID => micpar%get_micb_id(M,NGL))
              OMC1g=OMC1/numGuildsH
              OMN1g=OMN1/numGuildsH
              OMP1g=OMP1/numGuildsH
              mBiomeHeter_vr(ielmc,MID,K,LFDPTH,NY,NX)=mBiomeHeter_vr(ielmc,MID,K,LFDPTH,NY,NX)+OMC1g
              mBiomeHeter_vr(ielmn,MID,K,LFDPTH,NY,NX)=mBiomeHeter_vr(ielmn,MID,K,LFDPTH,NY,NX)+OMN1g
              mBiomeHeter_vr(ielmp,MID,K,LFDPTH,NY,NX)=mBiomeHeter_vr(ielmp,MID,K,LFDPTH,NY,NX)+OMP1g
            end associate
          ENDDO
          OSCX=OSCX+OMC1
          OSNX=OSNX+OMN1
          OSPX=OSPX+OMP1
          D2962: DO NN=1,NumMicbAFunGrupsPerCmplx
            numGuildsA=JGnfA(NN)-JGniA(NN)+1
            if(numGuildsA<=0)cycle D2962
            DO NGL=JGniA(NN),JGnfA(NN)
              associate(MID => micpar%get_micb_id(M,NGL))
                OMC1g=OMC1/numGuildsA
                OMN1g=OMN1/numGuildsA
                OMP1g=OMP1/numGuildsA
                mBiomeAutor_vr(ielmc,MID,LFDPTH,NY,NX)=mBiomeAutor_vr(ielmc,MID,LFDPTH,NY,NX)+OMC1g*OMCA(NN)
                mBiomeAutor_vr(ielmn,MID,LFDPTH,NY,NX)=mBiomeAutor_vr(ielmn,MID,LFDPTH,NY,NX)+OMN1g*OMCA(NN)
                mBiomeAutor_vr(ielmp,MID,LFDPTH,NY,NX)=mBiomeAutor_vr(ielmp,MID,LFDPTH,NY,NX)+OMP1g*OMCA(NN)
              end associate
            ENDDO
            OSCX=OSCX+OMC1*OMCA(NN)
            OSNX=OSNX+OMN1*OMCA(NN)
            OSPX=OSPX+OMP1*OMCA(NN)
          ENDDO D2962
        ENDDO D2961
      ENDDO D2960
!
!     DOC, DON AND DOP IN RESIDUE
!
!     OQC,OQN,OQP=DOC,DON,DOP in litter
!
      OQC1=AZMAX1(AMIN1(0.1_r8*OSCX,OSCI-OSCX))
      OQN1=AZMAX1(AMIN1(0.1_r8*OSNX,OSNI-OSNX))
      OQP1=AZMAX1(AMIN1(0.1_r8*OSPX,OSPI-OSPX))

      DOM_MicP_vr(idom_doc,K,LFDPTH,NY,NX)=DOM_MicP_vr(idom_doc,K,LFDPTH,NY,NX)+OQC1
      DOM_MicP_vr(idom_don,K,LFDPTH,NY,NX)=DOM_MicP_vr(idom_don,K,LFDPTH,NY,NX)+OQN1
      DOM_MicP_vr(idom_dop,K,LFDPTH,NY,NX)=DOM_MicP_vr(idom_dop,K,LFDPTH,NY,NX)+OQP1
!
!     REMAINDER DISTRIBUTED TO RESIDUE FRACTIONS
!
!     OSC,OSN,OSP,OSA=SOC,SON,SOP,colonized SOC in litter
!     VOLT=litter volume
!     AmendC_CumYr_flx_col,FertN_Flx_CumYr_col,FerP_Flx_CumYr_col=accumulated litter C,N,P application
!     Eco_NBP_CumYr_col=accumulated net biome productivity
!
      OSCX=OSCX+OQC1
      OSNX=OSNX+OQN1
      OSPX=OSPX+OQP1
      ! Normalize component weights independently of remaining carbon or
      ! nutrient mass. This is the same relative weighting without dividing
      ! by a remainder that can be zero after microbial/DOM allocation.
      RNT=0._r8
      RPT=0._r8
      DO M=1,jsken
        RNT=RNT+CFOSC_vr(M,K,LFDPTH,NY,NX)*CNOFC(M,K)
        RPT=RPT+CFOSC_vr(M,K,LFDPTH,NY,NX)*CPOFC(M,K)
      ENDDO
      D2970: DO M=1,jsken
        OSC1=CFOSC_vr(M,K,LFDPTH,NY,NX)*AZMAX1(OSCI-OSCX)
        IF(RNT.GT.ZERO)THEN
          OSN1=CFOSC_vr(M,K,LFDPTH,NY,NX)*CNOFC(M,K)/RNT*AZMAX1(OSNI-OSNX)
        ELSE
          OSN1=0._r8
        ENDIF
        IF(RPT.GT.ZERO)THEN
          OSP1=CFOSC_vr(M,K,LFDPTH,NY,NX)*CPOFC(M,K)/RPT*AZMAX1(OSPI-OSPX)
        ELSE
          OSP1=0._r8
        ENDIF
        SolidOM_vr(ielmc,M,K,LFDPTH,NY,NX)=SolidOM_vr(ielmc,M,K,LFDPTH,NY,NX)+OSC1
        SolidOM_vr(ielmn,M,K,LFDPTH,NY,NX)=SolidOM_vr(ielmn,M,K,LFDPTH,NY,NX)+OSN1
        SolidOM_vr(ielmp,M,K,LFDPTH,NY,NX)=SolidOM_vr(ielmp,M,K,LFDPTH,NY,NX)+OSP1
        SolidOMAct_vr(M,K,LFDPTH,NY,NX)=SolidOMAct_vr(M,K,LFDPTH,NY,NX)+OSC1*OMCI(1,K)
        IF(LFDPTH.EQ.0)THEN
          VGeomLayer_vr(LFDPTH,NY,NX)=VGeomLayer_vr(LFDPTH,NY,NX)+OSC1*ppmc/BulkDensLitR(k_fine_comp)
        ENDIF
      ENDDO D2970
      tAmendOrgC_lnd              = tAmendOrgC_lnd+OSCI
      TORGN                       = TORGN+OSNI
      TORGP                       = TORGP+OSPI
      AmendC_CumYr_flx_col(NY,NX) = AmendC_CumYr_flx_col(NY,NX)+OSCI
      FertN_Flx_CumYr_col(NY,NX)  = FertN_Flx_CumYr_col(NY,NX)+OSNI
      FerP_Flx_CumYr_col(NY,NX)    = FerP_Flx_CumYr_col(NY,NX)+OSPI
      IF(amendment==1)THEN
        ! Plant residue is external regardless of the animal-manure flag.
        Eco_NBP_CumYr_col(NY,NX)=Eco_NBP_CumYr_col(NY,NX)+OSCI
      ELSE
        !manure
        ! Grazing can merge internal returns into FERT and change its type.
        ! Preserve the independently read external share for application-day NBP.
        Eco_NBP_CumYr_col(NY,NX)=Eco_NBP_CumYr_col(NY,NX) &
          +ExternalManureC_col(I,NY,NX)*AREA_3D(3,LFDPTH,NY,NX)
      ENDIF
    ENDDO D2965
  ENDIF
  end associate
  call PrintInfo('end '//subname)
  end subroutine ApplyManure
!------------------------------------------------------------------------------------------

  subroutine ApplyMineralFertilizer(I,J,NY,NX,LFDPTH,OFC,OFN,OFP,UreaLayers,NitrifierLayers)
  implicit none
  integer, intent(in) :: I,J,NY,NX
  real(r8), intent(out) :: OFC(2),OFN(2),OFP(2)
  integer, intent(out) :: LFDPTH
  logical, intent(out) :: UreaLayers(0:JZ),NitrifierLayers(0:JZ)
  logical :: HasMineralFertilizer

  character(len=*), parameter :: subname='ApplyMineralFertilizer'
  real(r8) :: frcFert2Topsoil
  real(r8) :: frcFert2SurfLitr
  real(r8) :: CAC
  real(r8) :: CAS
  real(r8) :: CACX
  real(r8) :: CASX
  real(r8) :: FDPTHF
  real(r8) :: H0PO4T,H1PO4T,H2PO4T,H3PO4T
  real(r8) :: PMA,PMB,PHA
  real(r8) :: PALPOT,PFEPOT
  real(r8) :: PCAPDT,PCAPHT,PCAPMT
  real(r8) :: PMAX,PMBX,PHAX,POAX,POBX
  REAL(R8) :: XN4T
  real(r8) :: XOH0T,XOH1T,XOH2T,XH1PT,XH2PT
  real(r8) :: Z4A,Z3A,ZUA,ZOA,Z4B,Z3B
  REAL(R8) :: ZUB,ZOB
  real(r8) :: ZNH4T,ZNH3T,ZNO3T,ZNO2T
  real(r8) :: ZFE1PT,ZFE2PT
  real(r8) :: ZCA0PT,ZCA1PT,ZCA2PT
  real(r8) :: ZMG1PT,Z4AX,Z3AX,ZUAX,ZOAX
  real(r8) :: Z4BX,Z3BX,ZUBX,ZOBX,POA,POB
  integer :: L

!     begin_execution
!
!     NH4,NH3,UREA,NO3 FERTILIZER APPLICATION
!
!     *A,*B=broadcast,banded
!     Z4,Z3,ZU,ZO=NH4,NH3,urea,NO3
!
!   Z4A: NH4,Z3A:NH3,ZUA:urea,ZOA:NO3,Z4B:NH4Band,Z3B:NH3Band,ZUB:UreaBand,ZOB:NO3Band
!   PMA:CaH4P2O8,PMB:CaH4P2O8Band,PHA:apatite ((Ca5(PO4)3OH)),CAC:CaCO3,CAS:CaSO4
! 

  call PrintInfo('beg '//subname)
  Z4A=FERT(ifert_N_nh4,I,NY,NX)
  Z3A=FERT(ifert_N_nh3,I,NY,NX)
  ZUA=FERT(ifert_N_urea,I,NY,NX)
  ZOA=FERT(ifert_N_no3,I,NY,NX)
  Z4B=FERT(ifert_N_nh4_band,I,NY,NX)
  Z3B=FERT(ifert_N_nh3_band,I,NY,NX)
  ZUB=FERT(ifert_N_urea_band,I,NY,NX)
  ZOB=FERT(ifert_N_no3_band,I,NY,NX)
  POA=FERT(ifert_PO4_soil,I,NY,NX)
  POB=FERT(ifert_PO4_band,I,NY,NX)
!
!     MONOCALCIUM PHOSPHATE OR HYDROXYAPATITE
!
!     PM*,PH*=Ca(H2PO4)2,apatite
!
  PMA=FERT(ifert_P_Ca_H2PO4_2_soil,I,NY,NX)
  PMB=FERT(ifert_P_Ca_H2PO4_2_band,I,NY,NX)
  PHA=FERT(ifert_P_apatite,I,NY,NX)
!
!     LIME AND GYPSUM
!
!     CAC,CAS=CaCO3,CaSO4
!
  CAC=FERT(ifert_Ca_lime,I,NY,NX)
  CAS=FERT(ifert_Ca_gypsum,I,NY,NX)
!
!     PLANT(1) AND ANIMAL(2) RESIDUE C, N AND P
!
  OFC(1)=FERT(ifert_plant_resC,I,NY,NX)
  OFN(1)=FERT(ifert_plant_resN,I,NY,NX)
  OFP(1)=FERT(ifert_plant_resP,I,NY,NX)

  OFC(2)=FERT(ifert_plant_manuC,I,NY,NX)
  OFN(2)=FERT(ifert_plant_manuN,I,NY,NX)
  OFP(2)=FERT(ifert_plant_manuP,I,NY,NX)
!
!     SOIL LAYER NUMBER AT DEPTH OF FERTILIZER APPLICATION
!
!     LFDPTH=layer number
!     frcFert2SurfLitr=fraction of fertilizer applied to surface litter
!
  HasMineralFertilizer=Z4A+Z3A+ZUA+ZOA+Z4B+Z3B+ZUB+ZOB+PMA+PMB+PHA+CAC+CAS+POA+POB.GT.0._r8
  LFDPTH=0
  UreaLayers=.false.;NitrifierLayers=.false.
  ! Inhibitor-only records also need a placement determined from FDPTH.
  IF(HasMineralFertilizer .OR. IYTYP(iamendtyp_fert,I,NY,NX)==3 &
    .OR. IYTYP(iamendtyp_fert,I,NY,NX)==4)THEN
    FDPTHF=FDPTH(I,NY,NX)+CumDepz2LayBottom_vr(NU_col(NY,NX)-1,NY,NX)
    IF(FDPTHF.LE.0._r8.AND.isclose(Z4B+Z3B+ZUB+ZOB+PMB+POB,0._r8))THEN
      LFDPTH=0
      frcFert2SurfLitr=1.0_r8-EXP(-0.8E-02_r8*(SoilOrgM_vr(ielmc,0,NY,NX)/AREA_3D(3,0,NY,NX)))
    ELSE
      LFDPTH=-1
      D65: DO L=NUI_col(NY,NX),NL_col(NY,NX)
        IF(CumDepz2LayBottom_vr(L,NY,NX).GE.FDPTHF)THEN
          LFDPTH=L
          frcFert2SurfLitr=1.0_r8
          exit
        ENDIF
      ENDDO D65
      if(LFDPTH<0)call endrun('Fertilizer/inhibitor depth exceeds the soil profile in '//mod_filename,__LINE__)
    ENDIF

    frcFert2Topsoil=1.0_r8-frcFert2SurfLitr
    ! Record actual N destinations before any organic amendment changes LFDPTH.
    UreaLayers(LFDPTH)=(ZUA+ZUB)*frcFert2SurfLitr>0._r8
    NitrifierLayers(LFDPTH)=(Z4A+ZUA+ZOA+Z4B+ZUB+ZOB)*frcFert2SurfLitr>0._r8
    if(LFDPTH==0)then
      UreaLayers(NU_col(NY,NX))=(ZUA+ZUB)*frcFert2Topsoil>0._r8
      ! Ammonia fertilizer bypasses litter and enters topsoil directly.
      NitrifierLayers(NU_col(NY,NX))=(Z4A+ZUA+ZOA+Z4B+ZUB+ZOB)*frcFert2Topsoil+Z3A+Z3B>0._r8
    else
      NitrifierLayers(LFDPTH)=NitrifierLayers(LFDPTH) .or. (Z3A+Z3B>0._r8)
    endif
    if(.not.any(NitrifierLayers))then
      ! Without an N carrier, place the inhibitor using the application split.
      NitrifierLayers(LFDPTH)=frcFert2SurfLitr>0._r8
      if(LFDPTH==0)NitrifierLayers(NU_col(NY,NX))=frcFert2Topsoil>0._r8
    endif
  ENDIF

  IF(HasMineralFertilizer)THEN
!
!     RESET WIDTH AND DEPTH OF NH4 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWN=width of NH4 band row
!     DPNHB,WDNHB=depth,width of NH4 band
!     VLNHB,VLNH4=soil volume in NH4 band,non-band
!
    IF((Z4B+Z3B+ZUB.GT.0._r8) .OR. ((trcs_solml_vr(ids_NH4B,LFDPTH,NY,NX).GT.0._r8 &
      .OR. trcs_solml_vr(idg_NH3B,LFDPTH,NY,NX).GT.0._r8).AND.iFertNH4Band_col(NY,NX).EQ.ifert_off))THEN
      iFertNH4Band_col(NY,NX)=ifert_on
      ROWSpaceNH4_col(NY,NX)=ROWI(I,NY,NX)
      D50: DO L=NUI_col(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          BandThicknessNH4_vr(L,NY,NX)=DLYR_3D(3,L,NY,NX)
          BandWidthNH4_vr(L,NY,NX)=0._r8
        ELSEIF(L.EQ.LFDPTH)THEN
          BandThicknessNH4_vr(L,NY,NX)=AMAX1(0.025_r8,FDPTHF-CumDepz2LayBottom_vr(L-1,NY,NX))
          BandWidthNH4_vr(L,NY,NX)=AMIN1(0.025_r8,ROWSpaceNH4_col(NY,NX))
        ELSE
          BandThicknessNH4_vr(L,NY,NX)=0._r8
          BandWidthNH4_vr(L,NY,NX)=0._r8
        ENDIF
        IF(DLYR_3D(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN_vr(ids_NH4B,L,NY,NX)=AMIN1(0.999_r8,BandWidthNH4_vr(L,NY,NX)/ROWSpaceNH4_col(NY,NX) &
            *BandThicknessNH4_vr(L,NY,NX)/DLYR_3D(3,L,NY,NX))
        ELSE
          trcs_VLN_vr(ids_NH4B,L,NY,NX)=0._r8
        ENDIF
        trcs_VLN_vr(ids_NH4,L,NY,NX)    = 1.0_r8-trcs_VLN_vr(ids_NH4B,L,NY,NX)
        trcs_VLN_vr(idg_NH3B,L,NY,NX)   = trcs_VLN_vr(ids_NH4B,L,NY,NX)
        trcs_VLN_vr(idg_NH3,L,NY,NX)    = trcs_VLN_vr(ids_NH4,L,NY,NX)
        ZNH4T                           = trcs_solml_vr(ids_NH4,L,NY,NX)+trcs_solml_vr(ids_NH4B,L,NY,NX)
        ZNH3T                           = trcs_solml_vr(idg_NH3,L,NY,NX)+trcs_solml_vr(idg_NH3B,L,NY,NX)
        XN4T                            = trcx_solml_vr(idx_NH4,L,NY,NX)+trcx_solml_vr(idx_NH4B,L,NY,NX)
        trcs_solml_vr(ids_NH4,L,NY,NX)  = ZNH4T*trcs_VLN_vr(ids_NH4,L,NY,NX)
        trcs_solml_vr(idg_NH3,L,NY,NX)  = ZNH3T*trcs_VLN_vr(idg_NH3,L,NY,NX)
        trcs_solml_vr(ids_NH4B,L,NY,NX) = ZNH4T*trcs_VLN_vr(ids_NH4B,L,NY,NX)
        trcs_solml_vr(idg_NH3B,L,NY,NX) = ZNH3T*trcs_VLN_vr(idg_NH3B,L,NY,NX)
        trcx_solml_vr(idx_NH4,L,NY,NX)  = XN4T*trcs_VLN_vr(ids_NH4,L,NY,NX)
        trcx_solml_vr(idx_NH4B,L,NY,NX) = XN4T*trcs_VLN_vr(ids_NH4B,L,NY,NX)
      ENDDO D50
      BandDepthNH4_col(NY,NX)=BandThicknessNH4_vr(LFDPTH,NY,NX)+CumDepz2LayBottom_vr(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF NO3 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWO=width of NO3 band row
!     DPNOB,WDNOB=depth,width of NO3 band
!     VLNOB,VLNO3=soil volume in NO3 band,non-band
!
    IF((Z4B+Z3B+ZUB+ZOB.GT.0._r8) .OR. ((trcs_solml_vr(ids_NO3B,LFDPTH,NY,NX).GT.0._r8 &
      .OR. trcs_solml_vr(ids_NO2B,LFDPTH,NY,NX).GT.0._r8) .AND. iFertNO3Band_col(NY,NX).EQ.ifert_off))THEN
      iFertNO3Band_col(NY,NX)=ifert_on
      ROWSpaceNO3Band_col(NY,NX)=ROWI(I,NY,NX)
      D45: DO L=NUI_col(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          BandThicknessNO3_vr(L,NY,NX) = DLYR_3D(3,L,NY,NX)
          BandWidthNO3_vr(L,NY,NX)     = 0._r8
        ELSEIF(L.EQ.LFDPTH)THEN
          BandThicknessNO3_vr(L,NY,NX)=AMAX1(0.01_r8,FDPTHF-CumDepz2LayBottom_vr(L-1,NY,NX))
          BandWidthNO3_vr(L,NY,NX)=AMIN1(0.01_r8,ROWSpaceNO3Band_col(NY,NX))
        ELSE
          BandThicknessNO3_vr(L,NY,NX)=0._r8
          BandWidthNO3_vr(L,NY,NX)=0._r8
        ENDIF
        IF(DLYR_3D(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN_vr(ids_NO3B,L,NY,NX)=AMIN1(0.999_r8,BandWidthNO3_vr(L,NY,NX)/ROWSpaceNO3Band_col(NY,NX) &
            *BandThicknessNO3_vr(L,NY,NX)/DLYR_3D(3,L,NY,NX))
        ELSE
          trcs_VLN_vr(ids_NO3B,L,NY,NX)=0._r8
        ENDIF

        trcs_VLN_vr(ids_NO3,L,NY,NX)=1.0_r8-trcs_VLN_vr(ids_NO3B,L,NY,NX)
        trcs_VLN_vr(ids_NO2B,L,NY,NX)=trcs_VLN_vr(ids_NO3B,L,NY,NX)
        trcs_VLN_vr(ids_NO2,L,NY,NX)=trcs_VLN_vr(ids_NO3,L,NY,NX)
        ZNO3T=trcs_solml_vr(ids_NO3,L,NY,NX)+trcs_solml_vr(ids_NO3B,L,NY,NX)
        ZNO2T=trcs_solml_vr(ids_NO2,L,NY,NX)+trcs_solml_vr(ids_NO2B,L,NY,NX)

        trcs_solml_vr(ids_NO3,L,NY,NX)=ZNO3T*trcs_VLN_vr(ids_NO3,L,NY,NX)
        trcs_solml_vr(ids_NO2,L,NY,NX)=ZNO2T*trcs_VLN_vr(ids_NO2,L,NY,NX)
        trcs_solml_vr(ids_NO3B,L,NY,NX)=ZNO3T*trcs_VLN_vr(ids_NO3B,L,NY,NX)
        trcs_solml_vr(ids_NO2B,L,NY,NX)=ZNO2T*trcs_VLN_vr(ids_NO2B,L,NY,NX)
      ENDDO D45
      BandDepthNO3_col(NY,NX)=BandThicknessNO3_vr(LFDPTH,NY,NX)+CumDepz2LayBottom_vr(LFDPTH-1,NY,NX)
    ENDIF
!
!     RESET WIDTH AND DEPTH OF PO4 FERTILIZER BAND IF NEW BAND
!     AND ADD REMAINS OF ANY EXISTING FERTILIZER BAND TO NEW BAND
!
!     ROWP=width of H2PO4 band row
!     DPPOB,WDPOB=depth,width of H2PO4 band
!     VLPOB,VLPO4=soil volume in H2PO4 band,non-band
!
    IF((PMB+POB.GT.0.0_r8) .OR. (trcs_solml_vr(ids_H2PO4B,LFDPTH,NY,NX).GT.0._r8 .AND. iFertPO4Band_col(NY,NX).EQ.ifert_off))THEN
      iFertPO4Band_col(NY,NX)=ifert_on
      ROWSpacePO4_col(NY,NX)=ROWI(I,NY,NX)
      DO  L=NUI_col(NY,NX),JZ
        IF(L.LT.LFDPTH)THEN
          BandThicknessPO4_vr(L,NY,NX) = DLYR_3D(3,L,NY,NX)
          BandWidthPO4_vr(L,NY,NX)     = AMIN1(0.01,ROWSpacePO4_col(NY,NX))
        ELSEIF(L.EQ.LFDPTH)THEN
          BandThicknessPO4_vr(L,NY,NX) = AMAX1(0.01,FDPTHF-CumDepz2LayBottom_vr(L-1,NY,NX))
          BandWidthPO4_vr(L,NY,NX)     = AMIN1(0.01,ROWSpacePO4_col(NY,NX))
        ELSE
          BandThicknessPO4_vr(L,NY,NX) = 0._r8
          BandWidthPO4_vr(L,NY,NX)     = 0._r8
        ENDIF
        IF(DLYR_3D(3,L,NY,NX).GT.ZERO2)THEN
          trcs_VLN_vr(ids_H1PO4B,L,NY,NX)=AMIN1(0.999_r8,BandWidthPO4_vr(L,NY,NX)/ROWSpacePO4_col(NY,NX) &
          *BandThicknessPO4_vr(L,NY,NX)/DLYR_3D(3,L,NY,NX))
        ELSE
          trcs_VLN_vr(ids_H1PO4B,L,NY,NX)=0._r8
        ENDIF
        trcs_VLN_vr(ids_H1PO4,L,NY,NX)  = 1.0-trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcs_VLN_vr(ids_H2PO4B,L,NY,NX) = trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcs_VLN_vr(ids_H2PO4,L,NY,NX)  = trcs_VLN_vr(ids_H1PO4,L,NY,NX)

        H1PO4T = trcs_solml_vr(ids_H1PO4,L,NY,NX)+trcs_solml_vr(ids_H1PO4B,L,NY,NX)
        H2PO4T = trcs_solml_vr(ids_H2PO4,L,NY,NX)+trcs_solml_vr(ids_H2PO4B,L,NY,NX)
        XOH0T  = trcx_solml_vr(idx_OHe,L,NY,NX)+trcx_solml_vr(idx_OHeB,L,NY,NX)
        XOH1T  = trcx_solml_vr(idx_OH,L,NY,NX)+trcx_solml_vr(idx_OHB,L,NY,NX)
        XOH2T  = trcx_solml_vr(idx_OHp,L,NY,NX)+trcx_solml_vr(idx_OHpB,L,NY,NX)
        XH1PT  = trcx_solml_vr(idx_HPO4,L,NY,NX)+trcx_solml_vr(idx_HPO4B,L,NY,NX)
        XH2PT  = trcx_solml_vr(idx_H2PO4,L,NY,NX)+trcx_solml_vr(idx_H2PO4B,L,NY,NX)
        PALPOT = trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)+trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)
        PFEPOT = trcp_saltpml_vr(idsp_FePO4,L,NY,NX)+trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)
        PCAPDT = trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)+trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)
        PCAPHT = trcp_saltpml_vr(idsp_Apatite,L,NY,NX)+trcp_saltpml_vr(idsp_ApatiteBand,L,NY,NX)
        PCAPMT = trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)+trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX)

        trcs_solml_vr(ids_H1PO4,L,NY,NX)  = H1PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcs_solml_vr(ids_H2PO4,L,NY,NX)  = H2PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcs_solml_vr(ids_H1PO4B,L,NY,NX) = H1PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcs_solml_vr(ids_H2PO4B,L,NY,NX) = H2PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

        trcx_solml_vr(idx_OHe,L,NY,NX)          = XOH0T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_OH,L,NY,NX)           = XOH1T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_OHp,L,NY,NX)          = XOH2T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_HPO4,L,NY,NX)         = XH1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_H2PO4,L,NY,NX)        = XH2PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcx_solml_vr(idx_OHeB,L,NY,NX)         = XOH0T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_OHB,L,NY,NX)          = XOH1T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_OHpB,L,NY,NX)         = XOH2T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_HPO4B,L,NY,NX)        = XH1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcx_solml_vr(idx_H2PO4B,L,NY,NX)       = XH2PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_AlPO4,L,NY,NX)     = PALPOT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_FePO4,L,NY,NX)     = PFEPOT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_CaHPO4,L,NY,NX)    = PCAPDT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_Apatite,L,NY,NX)   = PCAPHT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_CaH4P2O8,L,NY,NX)  = PCAPMT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
        trcp_saltpml_vr(idsp_AlPO4B,L,NY,NX)    = PALPOT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_FePO4B,L,NY,NX)    = PFEPOT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_CaHPO4B,L,NY,NX)   = PCAPDT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_ApatiteBand,L,NY,NX)       = PCAPHT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        trcp_saltpml_vr(idsp_CaH4P2O8B,L,NY,NX) = PCAPMT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

        if(salt_model)then
          H0PO4T = trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)
          H3PO4T = trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)
          ZFE1PT = trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)
          ZFE2PT = trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)+trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)
          ZCA0PT = trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)
          ZCA1PT = trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)
          ZCA2PT = trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)+trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)
          ZMG1PT = trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)+trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)
          
          trcSalt_solml_vr(idsalt_H0PO4,L,NY,NX)=H0PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_H3PO4,L,NY,NX)=H3PO4T*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeHPO4,L,NY,NX)=ZFE1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeH2PO4,L,NY,NX)=ZFE2PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaPO4,L,NY,NX)=ZCA0PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaHPO4,L,NY,NX)=ZCA1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaH4P2O8,L,NY,NX)=ZCA2PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_MgHPO4,L,NY,NX)=ZMG1PT*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
          trcSalt_solml_vr(idsalt_H0PO4B,L,NY,NX)=H0PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)

          trcSalt_solml_vr(idsalt_H3PO4B,L,NY,NX)=H3PO4T*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeHPO4B,L,NY,NX)=ZFE1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_FeH2PO4B,L,NY,NX)=ZFE2PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaPO4B,L,NY,NX)=ZCA0PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaHPO4B,L,NY,NX)=ZCA1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_CaH4P2O8B,L,NY,NX)=ZCA2PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
          trcSalt_solml_vr(idsalt_MgHPO4B,L,NY,NX)=ZMG1PT*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
        endif

      ENDDO
      BandDepthPO4_col(NY,NX)=BandThicknessPO4_vr(LFDPTH,NY,NX)+CumDepz2LayBottom_vr(LFDPTH-1,NY,NX)
    ENDIF
!
!     UPDATE STATE VARIABLES FOR BROADCAST AND BANDED FERTILIZER
!     NH4, NH3, UREA, NO3, PO4, LIME AND GYPSUM IN SOIL
!     AND CONVERT FROM G TO MOLE
!
!     ZNH4FA,ZNH3FA,ZNHUFA,ZNO3FA=bdcast NH4,NH3,urea,NO3 fertilizer
!     ZNH4FB,ZNH3FB,ZNHUFB,ZNO3FB=banded NH4,NH3,urea,NO3 fertilizer
!     PCAPM1,PCAPD1,PCAPH1=concn of precip CaH4P2O8,CaHPO4,apatite in non-band
!     PCAPMB,PCAPDB,PCAPHB=concn of precip CaH4P2O8,CaHPO4,apatite in band
!     PCACO,PCASO=precipitated CaCO3,CaSO4

    Z4AX = Z4A*AREA_3D(3,LFDPTH,NY,NX)/natomw   !NH4 soil
    Z3AX = Z3A*AREA_3D(3,LFDPTH,NY,NX)/natomw   !NH3 soil
    ZUAX = ZUA*AREA_3D(3,LFDPTH,NY,NX)/natomw   !urea soil
    ZOAX = ZOA*AREA_3D(3,LFDPTH,NY,NX)/natomw   !NO3 soil
    Z4BX = Z4B*AREA_3D(3,LFDPTH,NY,NX)/natomw   !NH4 band
    Z3BX = Z3B*AREA_3D(3,LFDPTH,NY,NX)/natomw   !NH3 band
    ZUBX = ZUB*AREA_3D(3,LFDPTH,NY,NX)/natomw   !urea band
    ZOBX = ZOB*AREA_3D(3,LFDPTH,NY,NX)/natomw   !NO3 band
    POAX = POA*AREA_3D(3,LFDPTH,NY,NX)/patomw   !PO4 soil
    POBX = POB*AREA_3D(3,LFDPTH,NY,NX)/patomw   !PO4 band
    PMAX = PMA*AREA_3D(3,LFDPTH,NY,NX)/(2.0_r8*patomw)
    PMBX = PMB*AREA_3D(3,LFDPTH,NY,NX)/(2.0_r8*patomw)
    PHAX = PHA*AREA_3D(3,LFDPTH,NY,NX)/(3.0_r8*patomw)
    CACX = CAC*AREA_3D(3,LFDPTH,NY,NX)/40._r8   !lime CaCO3
    CASX = CAS*AREA_3D(3,LFDPTH,NY,NX)/40._r8   !gypsum CaSO4

    FertN_mole_soil_vr(ifert_N_nh4,LFDPTH,NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_N_nh4,LFDPTH,NY,NX)+Z4AX*frcFert2SurfLitr)
    FertN_mole_soil_vr(ifert_N_urea,LFDPTH,NY,NX) = AZMAX1(FertN_mole_soil_vr(ifert_N_urea,LFDPTH,NY,NX)+ZUAX*frcFert2SurfLitr)
    FertN_mole_soil_vr(ifert_N_no3,LFDPTH,NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_N_no3,LFDPTH,NY,NX)+ZOAX*frcFert2SurfLitr)

    ! Layer-zero applications are broadcast; banded pools exist only in soil.
    if(LFDPTH>0)then
      FertN_mole_Band_vr(ifert_N_nh4_band,LFDPTH,NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_N_nh4_band,LFDPTH,NY,NX)+Z4BX*frcFert2SurfLitr)
      FertN_mole_Band_vr(ifert_N_urea_band,LFDPTH,NY,NX) = AZMAX1(FertN_mole_Band_vr(ifert_N_urea_band,LFDPTH,NY,NX)+ZUBX*frcFert2SurfLitr)
      FertN_mole_Band_vr(ifert_N_no3_band,LFDPTH,NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_N_no3_band,LFDPTH,NY,NX)+ZOBX*frcFert2SurfLitr)
    endif

    ! Explicit PO4 enters dissolved pools as g P; mineral P remains in precipitates.
    trcs_solml_vr(ids_H2PO4,LFDPTH,NY,NX) = trcs_solml_vr(ids_H2PO4,LFDPTH,NY,NX) &
      +patomw*POAX*trcs_VLN_vr(ids_H1PO4,LFDPTH,NY,NX)*frcFert2SurfLitr
    trcs_solml_vr(ids_H2PO4B,LFDPTH,NY,NX) = trcs_solml_vr(ids_H2PO4B,LFDPTH,NY,NX) &
      +patomw*(POAX*trcs_VLN_vr(ids_H1PO4B,LFDPTH,NY,NX)+POBX)*frcFert2SurfLitr

    trcp_saltpml_vr(idsp_CaH4P2O8,LFDPTH,NY,NX)    = trcp_saltpml_vr(idsp_CaH4P2O8,LFDPTH,NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4,LFDPTH,NY,NX)*frcFert2SurfLitr
    trcp_saltpml_vr(idsp_CaH4P2O8B,LFDPTH,NY,NX)   = trcp_saltpml_vr(idsp_CaH4P2O8B,LFDPTH,NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4B,LFDPTH,NY,NX)*frcFert2SurfLitr+PMBX*frcFert2SurfLitr
    trcp_saltpml_vr(idsp_Apatite,LFDPTH,NY,NX)     = trcp_saltpml_vr(idsp_Apatite,LFDPTH,NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4,LFDPTH,NY,NX)*frcFert2SurfLitr
    trcp_saltpml_vr(idsp_ApatiteBand,LFDPTH,NY,NX) = trcp_saltpml_vr(idsp_ApatiteBand,LFDPTH,NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4B,LFDPTH,NY,NX)*frcFert2SurfLitr

    IF(LFDPTH.EQ.0)THEN
      FertN_mole_soil_vr(ifert_N_nh4,NU_col(NY,NX),NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_N_nh4,NU_col(NY,NX),NY,NX)+Z4AX*frcFert2Topsoil)
      FertN_mole_soil_vr(ifert_N_nh3,NU_col(NY,NX),NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_N_nh3,NU_col(NY,NX),NY,NX)+Z3AX)
      FertN_mole_soil_vr(ifert_N_urea,NU_col(NY,NX),NY,NX) = AZMAX1(FertN_mole_soil_vr(ifert_N_urea,NU_col(NY,NX),NY,NX)+ZUAX*frcFert2Topsoil)
      FertN_mole_soil_vr(ifert_N_no3,NU_col(NY,NX),NY,NX)  = AZMAX1(FertN_mole_soil_vr(ifert_N_no3,NU_col(NY,NX),NY,NX)+ZOAX*frcFert2Topsoil)

      FertN_mole_Band_vr(ifert_N_nh4_band,NU_col(NY,NX),NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_N_nh4_band,NU_col(NY,NX),NY,NX)+Z4BX*frcFert2Topsoil)
      FertN_mole_Band_vr(ifert_N_nh3_band,NU_col(NY,NX),NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_N_nh3_band,NU_col(NY,NX),NY,NX)+Z3BX)
      FertN_mole_Band_vr(ifert_N_urea_band,NU_col(NY,NX),NY,NX) = AZMAX1(FertN_mole_Band_vr(ifert_N_urea_band,NU_col(NY,NX),NY,NX)+ZUBX*frcFert2Topsoil)
      FertN_mole_Band_vr(ifert_N_no3_band,NU_col(NY,NX),NY,NX)  = AZMAX1(FertN_mole_Band_vr(ifert_N_no3_band,NU_col(NY,NX),NY,NX)+ZOBX*frcFert2Topsoil)

      ! Explicit PO4 enters dissolved pools as g P; mineral P remains in precipitates.
      trcs_solml_vr(ids_H2PO4,NU_col(NY,NX),NY,NX) = trcs_solml_vr(ids_H2PO4,NU_col(NY,NX),NY,NX) &
        +patomw*POAX*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)*frcFert2Topsoil
      trcs_solml_vr(ids_H2PO4B,NU_col(NY,NX),NY,NX) = trcs_solml_vr(ids_H2PO4B,NU_col(NY,NX),NY,NX) &
        +patomw*(POAX*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)+POBX)*frcFert2Topsoil

      trcp_saltpml_vr(idsp_CaH4P2O8,NU_col(NY,NX),NY,NX)    = trcp_saltpml_vr(idsp_CaH4P2O8,NU_col(NY,NX),NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)*frcFert2Topsoil
      trcp_saltpml_vr(idsp_CaH4P2O8B,NU_col(NY,NX),NY,NX)   = trcp_saltpml_vr(idsp_CaH4P2O8B,NU_col(NY,NX),NY,NX)+PMAX*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)*frcFert2Topsoil+PMBX*frcFert2Topsoil
      trcp_saltpml_vr(idsp_Apatite,NU_col(NY,NX),NY,NX)     = trcp_saltpml_vr(idsp_Apatite,NU_col(NY,NX),NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)*frcFert2Topsoil
      trcp_saltpml_vr(idsp_ApatiteBand,NU_col(NY,NX),NY,NX) = trcp_saltpml_vr(idsp_ApatiteBand,NU_col(NY,NX),NY,NX)+PHAX*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)*frcFert2Topsoil
    ELSE
      FertN_mole_soil_vr(ifert_N_nh3,LFDPTH,NY,NX)      = AZMAX1(FertN_mole_soil_vr(ifert_N_nh3,LFDPTH,NY,NX)+Z3AX*frcFert2SurfLitr)
      FertN_mole_Band_vr(ifert_N_nh3_band,LFDPTH,NY,NX) = AZMAX1(FertN_mole_Band_vr(ifert_N_nh3_band,LFDPTH,NY,NX)+Z3BX*frcFert2SurfLitr)
    ENDIF
    trcp_saltpml_vr(idsp_CaCO3,NU_col(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaCO3,NU_col(NY,NX),NY,NX)+CACX
    trcp_saltpml_vr(idsp_CaSO4,NU_col(NY,NX),NY,NX)=trcp_saltpml_vr(idsp_CaSO4,NU_col(NY,NX),NY,NX)+CASX
    TZIN                       = TZIN+natomw*(Z4AX+Z3AX+ZUAX+ZOAX+Z4BX+Z3BX+ZUBX+ZOBX)
    TPIN                       = TPIN+patomw*(POAX+POBX)+2._r8*patomw*(PMAX+PMBX)+3._r8*patomw*PHAX
    TIONIN                     = TIONIN+2.0_r8*(CACX+CASX)
    FertN_Flx_CumYr_col(NY,NX) = FertN_Flx_CumYr_col(NY,NX)+natomw*(Z4AX+Z4BX+Z3AX+Z3BX+ZUAX+ZUBX+ZOAX+ZOBX)
    FerP_Flx_CumYr_col(NY,NX)  = FerP_Flx_CumYr_col(NY,NX)+patomw*(POAX+POBX)+2._r8*patomw*(PMAX+PMBX)+3.0_r8*patomw*PHAX
  ENDIF

  call PrintInfo('end '//subname)  
  end subroutine ApplyMineralFertilizer

end module FertilizerMod
