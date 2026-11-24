module TranspSaltMod
!!
! Description:
!
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,      only: endrun,iulog
  use EcosimConst,      only: catomw, patomw
  use minimathmod  
  use DebugToolMod
  use SOMDataType
  use GridConsts
  use FlagDataType
  use SoilPhysDataType
  use EcoSIMSolverPar
  use EcoSIMCtrlDataType
  use ClimForcDataType
  use FertilizerDataType
  use SoilWaterDataType
  use SurfLitterDataType
  use SnowDataType
  use SurfSoilDataType
  use ChemTranspDataType
  use LandSurfDataType
  use SoilBGCDataType
  use RootDataType
  use AqueChemDatatype
  use SoilPropertyDataType
  use IrrigationDataType
  use PlantDataRateType
  use GridDataType
  use IngridTranspMod
  use TranspSaltDataMod
  implicit none

  private
  CHARACTER(LEN=*), PARAMETER :: MOD_FILENAME=&
  __FILE__

  public :: TranspSalt
  contains

!----------------------------------------------------------------------
  SUBROUTINE TranspSalt(I,J,NHW,NHE,NVN,NVS)
!
!     Description:
!
!     THIS SUBROUTINE CALCULATES 3-DIMENSIONAL FLUXES OF ALL SOIL
!     SALT SOLUTES
!
  implicit none
  integer, intent(in) :: I, J
  integer, intent(in) :: NHW,NHE,NVN,NVS

  integer :: NX,NY,L,M,ids
  integer, parameter :: max_salt_iterations=100
  integer :: iteration,stalled_iterations(idsalt_beg:idsaltb_end)
  character(len=*), parameter :: subname='TranspSalt'
  real(r8)  :: dpscal(idsalt_beg:idsaltb_end)
  real(r8)  :: pscal(idsalt_beg:idsaltb_end)
  real(r8) :: dpscal_max,remaining

!     execution begins here
!
!     TIME STEPS FOR SOLUTE FLUX CALCULATIONS
!
  call PrintInfo('beg '//subname)
  call InitSaltTransportModel(I,NHW,NHE,NVN,NVS)
!
!     TIME STEP USED IN GAS AND SOLUTE FLUX CALCULATIONS
!
  VFLWX=1._r8/NPH
  D30: DO M=1,NPH
   !
    pscal=1._r8;dpscal=1._r8
    dpscal_max=1._r8
    iteration=0
    stalled_iterations=0

    ! Apply ID-changing surface transfers once before independent soil retries.
    call ZeroFluxArraysM(NHW,NHE,NVN,NVS)
    call SurfaceSaltFluxM(M,NHW,NHE,NVN,NVS,band_transfer_only=.true.)

    do while(dpscal_max>1.e-2_r8)
      iteration=iteration+1
      if(iteration>max_salt_iterations)then
        write(iulog,*)subname,': retry limit at I,J,M,remaining=',I,J,M,dpscal_max
        call endrun(trim(mod_filename)//' at line',__LINE__)
      endif
      call ZeroFluxArraysM(NHW,NHE,NVN,NVS)

      call GetSaltTranspFlxM(I,M,NHW,NHE,NVN,NVS)
      !
      call UpdateSaltTranspM(M,NHW,NHE,NVN,NVS,dpscal,pscal)

      call AccumSaltFluxesM(M,NHW,NHE,NVN,NVS,dpscal,pscal)

      call copySaltStates(NHW,NHE,NVN,NVS)

      dpscal_max=0._r8
      DO ids=idsalt_beg,idsaltb_end
        remaining=dpscal(ids)*(1._r8-pscal(ids))
        if(dpscal(ids)>1.e-2_r8.and.remaining>=dpscal(ids))then
          stalled_iterations(ids)=stalled_iterations(ids)+1
        else
          stalled_iterations(ids)=0
        endif
        if(stalled_iterations(ids)>=3)then
          write(iulog,*)subname,': no accepted flux at I,J,M,tracer,remaining=',I,J,M,ids,remaining
          call endrun(trim(mod_filename)//' at line',__LINE__)
        endif
        dpscal(ids)=remaining
        dpscal_max=AMAX1(dpscal_max,dpscal(ids))
      enddo

    enddo
  ENDDO D30
  call PrintInfo('end '//subname)

  END subroutine TranspSalt
!------------------------------------------------------------------------------------------

  subroutine AccumSaltFluxesM(M,NHW,NHE,NVN,NVS,dpscal,pscal)
  implicit none
  integer, intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8), intent(in) :: dpscal(idsalt_beg:idsaltb_end)
  real(r8), intent(in):: pscal(idsalt_beg:idsaltb_end)

  real(r8) :: ppscal(idsalt_beg:idsaltb_end)
  integer  :: idsalt,NY,NX,L
  character(len=*), parameter :: subname='AccumSaltFluxesM'

  call PrintInfo('beg '//subname)
  ppscal=dpscal*pscal
  ! Replay every pool from the immutable staged baseline with the same
  ! accepted fraction used by the ledgers, even if only one tracer is limited.
  call UpdateSaltTranspM(M,NHW,NHE,NVN,NVS,ppscal)

  DO NX=NHW,NHE
    DO NY=NVN,NVS     
      DO idsalt=idsalt_beg,idsalt_end
        trcSalt_AquaADV_Snow2Litr_flx(idsalt,NY,NX)  = trcSalt_AquaADV_Snow2Litr_flx(idsalt,NY,NX)+ppscal(idsalt) &
          * trcSalt_AquaADV_Snow2Litr_flxM(idsalt,NY,NX)

        trcSalt_SnowDrift_flx_col(idsalt,NY,NX) = trcSalt_SnowDrift_flx_col(idsalt,NY,NX) + ppscal(idsalt) &
          * trcSalt_SnowDrift_flxM_col(idsalt,NY,NX)  
      ENDDO

      DO idsalt=idsalt_beg,idsaltb_end
        trcSalt_AquaADV_Snow2Soil_flx(idsalt,NY,NX)  = trcSalt_AquaADV_Snow2Soil_flx(idsalt,NY,NX)+ppscal(idsalt) &
          *trcSalt_AquaADV_Snow2Soil_flxM(idsalt,NY,NX)
      ENDDO

      call SummarizeSaltBoundaryFluxesM(NY,NX,ppscal)

      DO L=NU_col(NY,NX),NL_col(NY,NX)
        IF(VLSoilPoreMicP_vr(L,NY,NX)<=ZEROS2(NY,NX))cycle
        DO idsalt=idsalt_beg,idsaltb_end
          trcSalt_Mac2MicPore_flx_vr(idsalt,L,NY,NX)=trcSalt_Mac2MicPore_flx_vr(idsalt,L,NY,NX) &
            +ppscal(idsalt)*trcSalt_Mac2MicPore_flxM_vr(idsalt,L,NY,NX)
        ENDDO
      ENDDO

    ENDDO
  ENDDO  

  
  call PrintInfo('end '//subname)
  end subroutine AccumSaltFluxesM
!------------------------------------------------------------------------------------------

  subroutine SummarizeSaltBoundaryFluxesM(NY,NX,ppscal)
  ! Add accepted salt transport as elemental C/P exports, in grams.
  ! Internal snow/litter/soil and micro/macropore transfers are excluded.
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(in) :: ppscal(idsalt_beg:idsaltb_end)
  real(r8) :: carbon_weight(idsalt_beg:idsaltb_end), phosphorus_weight(idsalt_beg:idsaltb_end)
  real(r8) :: surface_export, subsurface_export
  integer :: idsalt,L
  character(len=*), parameter :: subname='SummarizeSaltBoundaryFluxesM'

  call PrintInfo('beg '//subname)
  carbon_weight=0._r8
  carbon_weight([idsalt_CO3,idsalt_HCO3,idsalt_CaCO3,idsalt_CaHCO3, &
    idsalt_MgCO3,idsalt_MgHCO3,idsalt_NaCO3])=catomw
  phosphorus_weight=0._r8
  ! Dissolved CaH4P2O8 IDs are legacy names for CaH2PO4(+), with one P atom.
  phosphorus_weight(idsalt_psoil_beg:idsalt_psoil_end)=patomw
  phosphorus_weight(idsalt_pband_beg:idsalt_pband_end)=patomw

  do idsalt=idsalt_beg,idsaltb_end
    if(carbon_weight(idsalt)==0._r8.and.phosphorus_weight(idsalt)==0._r8)cycle
    surface_export=0._r8
    if(idsalt<=idsalt_end)surface_export=-ppscal(idsalt)*(trcSalt_FloXSurRof_flxM_col(idsalt,NY,NX) &
      +trcSalt_SnowDrift_flxM_col(idsalt,NY,NX))
    ! Summed soil transport is top influx plus lateral/bottom net gain.
    ! Removing top influx excludes deposition and internal infiltration.
    L=NU_col(NY,NX)
    subsurface_export=trcSalt_MicpTranspFlxM_3D(idsalt,3,L,NY,NX) &
      +trcSalt_MacpTranspFlxM_3D(idsalt,3,L,NY,NX)
    do L=NU_col(NY,NX),NL_col(NY,NX)
      if(VLSoilPoreMicP_vr(L,NY,NX)<=ZEROS2(NY,NX))cycle
      subsurface_export=subsurface_export-trcSalt_Transp2Micp_flxM_vr(idsalt,L,NY,NX) &
        -trcSalt_Transp2Macp_flxM_vr(idsalt,L,NY,NX)
    enddo
    subsurface_export=subsurface_export*ppscal(idsalt)
    HydroSufDICFlx_col(NY,NX)=HydroSufDICFlx_col(NY,NX)+carbon_weight(idsalt)*surface_export
    HydroSubsDICFlx_col(NY,NX)=HydroSubsDICFlx_col(NY,NX)+carbon_weight(idsalt)*subsurface_export
    ! Each accepted fraction is accumulated once; the surface P total is annual.
    HydroSufDIPFlx_CumYr_col(NY,NX)=HydroSufDIPFlx_CumYr_col(NY,NX)+phosphorus_weight(idsalt)*surface_export
    HydroSubsDIPFlx_col(NY,NX)=HydroSubsDIPFlx_col(NY,NX)+phosphorus_weight(idsalt)*subsurface_export
  enddo
  call PrintInfo('end '//subname)
  end subroutine SummarizeSaltBoundaryFluxesM
!------------------------------------------------------------------------------------------

  subroutine ZeroAtmosSoluteFlux(NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: NY,NX
  integer :: idsalt

  ! Reset hourly atmospheric input before selecting current deposition.
  trcSalt_AquaAdv_flxM_snvr(:,1,NY,NX)=0._r8

  trcSalt_Precip2LitrM(:,NY,NX) = 0._r8

  trcSalt_Precip2MicpM(:,NY,NX) = 0._r8

  end subroutine ZeroAtmosSoluteFlux

!------------------------------------------------------------------------------------------

  subroutine AtmosSoluteFluxToTopsoil(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
  integer :: idsalt,ids
  character(len=*), parameter :: subname='AtmosSoluteFluxToTopsoil'


!     begin_execution
  call PrintInfo('beg '//subname)
  !
  DO idsalt=idsalt_beg,idsalt_end    
    trcSalt_Precip2LitrM(idsalt,NY,NX)   =(Rain2LitRSurf_col(NY,NX)*trcsalt_rain_mole_conc_col(idsalt,NY,NX) &
      +Irrig2LitRSurf_col(NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX))*dts_HeatWatTP  
  ENDDO


  DO idsalt=idsalt_beg,idsalt_KSO4
    trcSalt_Precip2MicpM(idsalt,NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcsalt_rain_mole_conc_col(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX))*dts_HeatWatTP
  ENDDO
  
  DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
    ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B  
    !soil 
    trcSalt_Precip2MicpM(idsalt,NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcsalt_rain_mole_conc_col(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX))*trcs_VLN_vr(ids_H1PO4,NU_col(NY,NX),NY,NX)*dts_HeatWatTP
    !band   
    trcSalt_Precip2MicpM(ids,NY,NX)=(Rain2SoilSurf_col(NY,NX)*trcsalt_rain_mole_conc_col(idsalt,NY,NX) &
      +Irrig2SoilSurf_col(NY,NX)*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX))*trcs_VLN_vr(ids_H1PO4B,NU_col(NY,NX),NY,NX)*dts_HeatWatTP    
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine AtmosSoluteFluxToTopsoil

!------------------------------------------------------------------------------------------

  subroutine GetSubHourFlux(I,NY,NX)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NY,NX
  integer :: L,idsalt,ids
  real(r8) :: IrrigWater


!     begin_execution
!
!     dts_HeatWatTP=1/no. of cycles h-1 for water, heat and solute flux calculations

  D10: DO L=NU_col(NY,NX),NL_col(NY,NX)
    DO idsalt=idsalt_beg,idsaltb_end
      trcSalt_RGeoChem_flxM_vr(idsalt,L,NY,NX)=-trcSalt_RGeoChem_flx_vr(idsalt,L,NY,NX)*dts_HeatWatTP
    ENDDO

    trcSalt_RGeoChem_flxM_vr(idsalt_Hp,L,NY,NX)=trcSalt_RGeoChem_flxM_vr(idsalt_Hp,L,NY,NX)-(RProd_Hp_vr(L,NY,NX))*dts_HeatWatTP
    !
    ! SOLUTE FLUXES FROM SUBSURFACE IRRIGATION
    ! placeholder: salt uptake by plant roots can be applied here. 
    !
    ! Keep trcSalt_Irrig_vr as the prescribed boundary concentration.
    ! Irrigation mass inputs belong only in the transport-substep flux array.
    IrrigWater=FWatIrrigate2MicP_vr(L,NY,NX)*dts_HeatWatTP
    trcSalt_Irrig_flxM_vr(:,L,NY,NX)=0._r8

    DO idsalt=idsalt_beg,idsalt_KSO4
      trcSalt_Irrig_flxM_vr(idsalt,L,NY,NX)=IrrigWater*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX)
    ENDDO

    DO idsalt=idsalt_H0PO4,idsalt_MgHPO4
      ids=idsalt-idsalt_H0PO4+idsalt_H0PO4B
      trcSalt_Irrig_flxM_vr(idsalt,L,NY,NX)=IrrigWater*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX)*trcs_VLN_vr(ids_H1PO4,L,NY,NX)
      trcSalt_Irrig_flxM_vr(ids,L,NY,NX)=IrrigWater*trcsalt_irrig_mole_conc_col(idsalt,I,NY,NX)*trcs_VLN_vr(ids_H1PO4B,L,NY,NX)
    ENDDO
    !
    !     SOLUTE DIFFUSIVITIES AT SUB-HOURLY TIME STEP
    !
    SoluteDifusivitytscaledM_vr(L,NY,NX)                 = SoluteDifusvtyT_vr(ids_H1PO4,L,NY,NX)*dts_HeatWatTP
    AquaIonDifusivty2_vr(idsalt_beg:idsalt_mend,L,NY,NX) = AquaIonDifusivty_vr(idsalt_beg:idsalt_mend,L,NY,NX)*dts_HeatWatTP

  ENDDO D10
  end subroutine GetSubHourFlux

!------------------------------------------------------------------------------------------
  subroutine copySaltStates(NHW,NHE,NVN,NVS)

  implicit none
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX
  integer :: L,idsalt

  DO NX=NHW,NHE
    DO NY=NVN,NVS
      D20: DO L=1,JS
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_ml2_snvr(idsalt,L,NY,NX)=trcSalt_ml_snvr(idsalt,L,NY,NX)
        ENDDO
      ENDDO D20

      L=0
      DO idsalt=idsalt_beg,idsaltb_end
        trcSalt_solml2_vr(idsalt,L,NY,NX)=trcSalt_solml_vr(idsalt,L,NY,NX)
      ENDDO
      DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO idsalt=idsalt_beg,idsaltb_end
          trcSalt_solml2_vr(idsalt,L,NY,NX)=trcSalt_solml_vr(idsalt,L,NY,NX)
          trcSalt_soHml2_vr(idsalt,L,NY,NX)=trcSalt_soHml_vr(idsalt,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  end subroutine copySaltStates
!------------------------------------------------------------------------------------------

  subroutine InitSaltTransportModel(I,NHW,NHE,NVN,NVS)
!
!     Description:
!
  implicit none
  integer, intent(in) :: I,NHW,NHE,NVN,NVS
  character(len=*), parameter :: subname='InitSaltTransportModel'
  integer :: NY,NX

  !     begin_execution
  call PrintInfo('beg '//subname)

  D9995: DO NX=NHW,NHE
    D9990: DO NY=NVN,NVS
      !
      call ZeroAtmosSoluteFlux(NY,NX)      

      ! Partial snow cover receives snow and exposed-ground deposition together.
      call AtmosSoluteFluxToTopsoil(I,NY,NX)
      call SaltFall2Snowpack(I,NY,NX)
      ! Stage soil sources independently of surface precipitation.
      call GetSubHourFlux(I,NY,NX)

    ENDDO D9990
  ENDDO D9995

  call copySaltStates(NHW,NHE,NVN,NVS)  
  call PrintInfo('end '//subname)
  end subroutine InitSaltTransportModel


!------------------------------------------------------------------------------------------

  subroutine UpdateSaltTranspM(M,NHW,NHE,NVN,NVS,dpscal,pscal)
!
!     Description:
!
  implicit none
  integer , intent(in) :: M,NHW,NHE,NVN,NVS
  real(r8), intent(in) :: dpscal(idsalt_beg:idsaltb_end)
  real(r8), optional, intent(out):: pscal(idsalt_beg:idsaltb_end)

  character(len=*), parameter :: subname='UpdateSaltTranspM'
  integer :: L,idsalt,NY,NX
  real(r8) :: pscal1(idsalt_beg:idsaltb_end)
  real(r8) :: flux
  !
  !     begin_execution
  !
  call PrintInfo('beg '//subname)

  ! *2 arrays are the immutable baseline; primary arrays hold trial results.
  ! Accumulate a common per-tracer donor limit across all pools and cells.
  ! The local array also makes replay without optional pscal safe.
  pscal1=1._r8
  D9695: DO NX=NHW,NHE
    D9690: DO NY=NVN,NVS
      !snow
      DO idsalt=idsalt_beg,idsalt_end
        flux=trcSalt_SnowDrift_flxM_col(idsalt,NY,NX)+trcSalt_Aqua_flxM_snvr(idsalt,1,NY,NX)
        flux=flux*dpscal(idsalt)
        call get_flux_scalar(trcSalt_ml2_snvr(idsalt,1,NY,NX),flux,trcSalt_ml_snvr(idsalt,1,NY,NX),pscal1(idsalt))
      ENDDO

      D9670: DO L=2,JS
        DO idsalt=idsalt_beg,idsalt_end
          flux=trcSalt_Aqua_flxM_snvr(idsalt,L,NY,NX)
          flux=flux*dpscal(idsalt)
          call get_flux_scalar(trcSalt_ml2_snvr(idsalt,L,NY,NX),flux,trcSalt_ml_snvr(idsalt,L,NY,NX),pscal1(idsalt))
        ENDDO
      ENDDO D9670
      !litter
      DO idsalt=idsalt_beg,idsalt_end        
        flux=trcSalt_FloXSurRof_flxM_col(idsalt,NY,NX)+trcSalt_MicpTranspFlxM_3D(idsalt,3,0,NY,NX)
        flux=flux*dpscal(idsalt)
        call get_flux_scalar(trcSalt_solml2_vr(idsalt,0,NY,NX),flux,trcSalt_solml_vr(idsalt,0,NY,NX),pscal1(idsalt))
      ENDDO

      D9685: DO L=NU_col(NY,NX),NL_col(NY,NX)
        IF(VLSoilPoreMicP_vr(L,NY,NX).GT.ZEROS2(NY,NX))THEN
          DO idsalt=idsalt_beg,idsaltb_end            
            flux=trcSalt_Transp2Micp_flxM_vr(idsalt,L,NY,NX)+trcSalt_Mac2MicPore_flxM_vr(idsalt,L,NY,NX) &
              +trcSalt_Irrig_flxM_vr(idsalt,L,NY,NX)-trcSalt_RGeoChem_flxM_vr(idsalt,L,NY,NX)
            flux=flux*dpscal(idsalt)
            call get_flux_scalar(trcSalt_solml2_vr(idsalt,L,NY,NX),flux,trcSalt_solml_vr(idsalt,L,NY,NX),pscal1(idsalt))

            flux=trcSalt_Transp2Macp_flxM_vr(idsalt,L,NY,NX)-trcSalt_Mac2MicPore_flxM_vr(idsalt,L,NY,NX)
            flux=flux*dpscal(idsalt)
            call get_flux_scalar(trcSalt_soHml2_vr(idsalt,L,NY,NX),flux,trcSalt_soHml_vr(idsalt,L,NY,NX),pscal1(idsalt))
            
          ENDDO
        ENDIF
      ENDDO D9685
    ENDDO D9690
  ENDDO D9695
  if(present(pscal))pscal=pscal1
  call PrintInfo('end '//subname)
  end subroutine UpdateSaltTranspM

!------------------------------------------------------------------------------------------

  subroutine ZeroFluxArraysM(NHW,NHE,NVN,NVS)

  implicit none
  integer , intent(in) :: NHW,NHE,NVN,NVS

  integer :: NY,NX
  integer :: L, idsalt

  ! A face may disappear after relayering. Clear all owner slots, including
  ! those above the destination's active surface or without a valid target.
  trcSalt_MicpTranspFlxM_3D=0._r8
  trcSalt_MacpTranspFlxM_3D=0._r8

  DO NX=NHW,NHE
    DO  NY=NVN,NVS
!
!     INITIALIZE SNOWPACK NET FLUX ACCUMULATORS
!
      D9855: DO L=1,JS
        DO idsalt=idsalt_beg,idsalt_end
          trcSalt_Aqua_flxM_snvr(idsalt,L,NY,NX)=0.0
        ENDDO
      ENDDO D9855
!
!     INITIALIZE SOIL SOLUTE NET FLUX ACCUMULATORS
!
      ! Include banded phosphate so retries start with no previous net flux.
      D9885: DO L=NU_col(NY,NX),NL_col(NY,NX)
        DO idsalt=idsalt_beg,idsaltb_end
          trcSalt_Transp2Micp_flxM_vr(idsalt,L,NY,NX)=0.0_r8
          trcSalt_Transp2Macp_flxM_vr(idsalt,L,NY,NX)=0.0_r8
        enddo
      ENDDO D9885      
    ENDDO
  ENDDO
  end subroutine ZeroFluxArraysM


end module TranspSaltMod
