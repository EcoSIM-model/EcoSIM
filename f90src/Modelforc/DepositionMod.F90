module DepositionMod
  use data_kind_mod,     only: r8 => DAT_KIND_R8
  use EcoSiMParDataMod, only: micpar
  use DebugToolMod,     only: PrintInfo
  use InitSOMBGCMOD,    only: MicrobeByLitterFall, InoculateCyanoBacter
  use ClimForcDataType
  use EcosimBGCFluxType
  use GridDataType
  use SOMDataType, only: DOM_MicP_vr, OMBioResdu_vr, SolidOM_vr, jsken
  use ElmIDMod
  use TracerIDMod, only: idom_doc, idom_don, idom_dop
  use EcoSIMHistMod, only: DATAP
  use PlantMgmtDataType, only: NP_col

implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: ApplyBioAerosol
contains


!------------------------------------------------------------------------------------------
  subroutine ApplyBioAerosol(I,J,NY,NX)
  !
  !apply bio-aerosol, made up by OM and
  !microbes
  !for simplicity, it is added to complex
  !fine litter.
  implicit none
  integer, intent(in) :: I,J
  integer, intent(in) :: NY,NX
  real(r8) :: OSCMK,OSCDep,ODMBC,ODNMB
  real(r8) :: OME_in(NumPlantChemElms)
  integer :: K,KL,NZ,M
  real(r8), parameter :: CFOSC(4)=(/0.075,0.125,0.550,0.250/)

  character(len=*), parameter :: subname='ApplyBioAerosol'

  call PrintInfo('beg '//subname)
  IF(DryDepoOMC(I)<=1.e-12_R8)return

  !add to fine litter group
  K=micpar%k_fine_comp

  !total organic C dry deposition
  OSCDep=DryDepoOMC(I)*AREA_3D(3,NU_col(NY,NX),NY,NX)*fAeroScalar_col(NY,NX)

  !total live microbial C deposition
  OSCMK =OSCDep*f_aerosol_LiveMB_col(NY,NX)

  !assuming 97.5% as non-cyanobacteria
  call MicrobeByLitterFall(I,J,K,NY,NX,OSCMK*0.975_r8,OME_in, mscal=1._r8)

  !assuming 2.5% as cyanobacteria
  KL=1
  call InoculateCyanoBacter(K,0,NY,NX,KL,CyanoInocC=OSCMK*0.025_r8,add_to_existing=.true.)

  !add DOM
  DOM_MicP_vr(idom_doc,K,0,NY,NX) = DOM_MicP_vr(idom_doc,K,0,NY,NX)+OSCDep*f_aerosol_DOM_col(NY,NX)
  DOM_MicP_vr(idom_don,K,0,NY,NX) = DOM_MicP_vr(idom_don,K,0,NY,NX)+OSCDep*f_aerosol_DOM_col(NY,NX)*0.05_r8
  DOM_MicP_vr(idom_dop,K,0,NY,NX) = DOM_MicP_vr(idom_dop,K,0,NY,NX)+OSCDep*f_aerosol_DOM_col(NY,NX)*0.005_r8

  !add dead microbial residual
  ODMBC=OSCDep*f_aerosol_DeadMB_col(NY,NX)
  OMBioResdu_vr(ielmc,iDbiom_labile,K,0,NY,NX) = OMBioResdu_vr(ielmc,iDbiom_labile,K,0,NY,NX)+ODMBC*micpar%FL(iLbiom_kinetic)
  OMBioResdu_vr(ielmc,iDbiom_recalc,K,0,NY,NX) = OMBioResdu_vr(ielmc,iDbiom_recalc,K,0,NY,NX)+ODMBC*micpar%FL(iLbiom_struct)
  OMBioResdu_vr(ielmn,iDbiom_labile,K,0,NY,NX) = OMBioResdu_vr(ielmn,iDbiom_labile,K,0,NY,NX)+ODMBC*micpar%FL(iLbiom_kinetic)*0.15_r8
  OMBioResdu_vr(ielmn,iDbiom_recalc,K,0,NY,NX) = OMBioResdu_vr(ielmn,iDbiom_recalc,K,0,NY,NX)+ODMBC*micpar%FL(iLbiom_struct) *0.09_r8
  OMBioResdu_vr(ielmp,iDbiom_labile,K,0,NY,NX) = OMBioResdu_vr(ielmp,iDbiom_labile,K,0,NY,NX)+ODMBC*micpar%FL(iLbiom_kinetic)*0.015_r8
  OMBioResdu_vr(ielmp,iDbiom_recalc,K,0,NY,NX) = OMBioResdu_vr(ielmp,iDbiom_recalc,K,0,NY,NX)+ODMBC*micpar%FL(iLbiom_struct) *0.009_r8

  !add dead nonmicroibal biomass, litter, dead moss/lichen
  !use the default partitioning into for kinetic components, CFOSC=(/0.075,0.125,0.550,0.250/)
  !the CNP ratio used below is very crude estimate
  ODNMB=OSCDep*f_aerosol_DeadNMB_col(NY,NX)
  DO M=1,jsken
    SolidOM_vr(ielmc,M,K,0,NY,NX)=  SolidOM_vr(ielmc,M,K,0,NY,NX)+ODNMB*CFOSC(M)
    SolidOM_vr(ielmn,M,K,0,NY,NX)=  SolidOM_vr(ielmn,M,K,0,NY,NX)+ODNMB*CFOSC(M)*0.02_r8
    SolidOM_vr(ielmp,M,K,0,NY,NX)=  SolidOM_vr(ielmp,M,K,0,NY,NX)+ODNMB*CFOSC(M)*0.002_r8
  ENDDO

  CumDryDepoOM_col(NY,NX)=CumDryDepoOM_col(NY,NX)+OSCDep

  !assign seed deposition for lichen and moss
  SeedCDeposition_pft(:,NY,NX)=0._r8;
  !each NZ corresponds to a different pft
  DO NZ=1,NP_col(NY,NX)
    if(DATAP(NZ,NY,NX)(1:4)=='lich')then
      SeedCDeposition_pft(NZ,NY,NX) = DryDepoOMC(I)*f_aerosol_LichB_col(NY,NX)
    elseif(DATAP(NZ,NY,NX)(1:4)=='moss')then
      SeedCDeposition_pft(NZ,NY,NX) = DryDepoOMC(I)*f_aerosol_MossB_col(NY,NX)
    endif
  ENDDO
  call PrintInfo('end '//subname)
  end subroutine ApplyBioAerosol

end module DepositionMod
