module PlantDebugMod
  use data_kind_mod, only : r8 => DAT_KIND_R8, yearIJ_type
  use abortutils, only : endrun
  use PlantBalMod
  use TracerIDMod
  use PlantAPICommonData
  use PlantSiteAPIData, only : plt_site
  use PlantMorphologyAPIData, only : plt_morph
  use PlantBiomassAPIData, only : plt_biom
  use PlantDisturbanceAPIData, only : plt_distb
  use PlantBGCRatesAPIData, only : plt_bgcr
  use PlantRootBGCAPIData, only : plt_rbgc

implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
  __FILE__
  public :: PrintRootTracer, PrintRootBiomass
  public :: RootMassBalCheck
  contains
  ![header]
!----------------------------------------------------------------------------------------------------
  subroutine PrintRootTracer(I,J,NZ,header)

  implicit none
  integer, intent(in) :: I,J,NZ
  character(len=*), intent(in) :: header
  associate(                                      &
    NU              => plt_site%NU               ,& !input  :current soil surface layer number, [-]
    NK              => plt_site%NK               ,& !input  :current hydrologically active layer, [-]
    Myco_pft        => plt_morph%Myco_pft        ,& !input  :mycorrhizal type (no or yes),[-]
    trcg_rootml_pvr => plt_rbgc%trcg_rootml_pvr   & !input  :root gas content, [g d-2]
  )
  return
  write(223,*)NZ,Myco_pft(NZ),NU,NK
  write(223,*)I*1000+J,trim(header),sum(trcg_rootml_pvr(idg_CO2,1:Myco_pft(NZ),NU:NK,NZ))+sum(trcg_rootml_pvr(idg_CO2,1:Myco_pft(NZ),NU:NK,NZ))
  end associate
  end subroutine PrintRootTracer

!----------------------------------------------------------------------------------------------------
  subroutine PrintRootBiomass(yearIJ,L,NR,NZ,header)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: L,NR,NZ
  character(len=*), intent(in) :: header
  real(r8) :: err

  associate(                                                          &
    RootMyco1stStrutElms_rpvr => plt_biom%RootMyco1stStrutElms_rpvr  ,& !inoput :root layer element primary axes, [g d-2]
    Root1stActStructElms_rpvr => plt_biom%Root1stActStructElms_rpvr  ,& !inoput :root layer active zone element in primary axes, [g d-2]
    Root1stLigStructElms_rpvr => plt_biom%Root1stLigStructElms_rpvr   & !inoput :root layer lignified zone element in primary axes, [g d-2]
  )

    err = RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ) &
      - Root1stActStructElms_rpvr(ielmc,L,NR,NZ)   &
      - Root1stLigStructElms_rpvr(ielmc,L,NR,NZ)
  
  write(999,*)yearIJ%I*1000+yearIJ%J/24._r8,err,L,NR,trim(header),RootMyco1stStrutElms_rpvr(ielmc,L,NR,NZ)
  if(abs(err)>1.e-8_r8)then    
    call endrun('C balance error test failure in '//trim(mod_filename)//' at line',__LINE__)                    
  endif
  end associate
  end subroutine PrintRootBiomass

!----------------------------------------------------------------------------------------------------
  subroutine RootMassBalCheck(yearIJ,NZ,opt,checktype,iut,info,id)
  implicit none
  type(yearIJ_type), intent(in) :: yearIJ
  integer, intent(in) :: NZ
  character(len=*), intent(in) :: opt         !'enter' or 'exit'
  character(len=1), intent(in) :: checktype   !'C','N' or 'P'
  integer, intent(in) :: iut                  !
  character(len=*), intent(in) :: info
  integer, optional, intent(in) :: id
  character(len=*), parameter :: subname='MassBalCheck'
  real(r8), save :: tmpval
  real(r8), save :: mass_inital(NumPlantChemElms)
  real(r8) :: mass_finale(NumPlantChemElms)
  real(r8) :: masserr
  integer :: id_loc

  if(opt=='enter')then
    call SumRootBiome(yearIJ,NZ,mass_inital)
    call SumRootAR(NZ);call SumLitfallBlg(NZ);
    tmpval=-plt_bgcr%RootAutoCO2_pft(NZ)+plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ)+ &
      plt_distb%RootLost2Fire_pft(ielmc,NZ)-plt_bgcr%Xfer2RootsC_pft(NZ)
  else
    id_loc=0
    if(present(id))id_loc=id

    call SumRootBiome(yearIJ,NZ,mass_finale)
    call SumRootAR(NZ);call SumLitfallBlg(NZ);
    masserr=mass_finale(ielmc)-mass_inital(ielmc)- &
        plt_bgcr%RootAutoCO2_pft(NZ)+plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ)+plt_distb%RootLost2Fire_pft(ielmc,NZ) &
        -tmpval-plt_bgcr%Xfer2RootsC_pft(NZ)

    write(iut,*)yearIJ%I*1000+yearIJ%J/24.,masserr,'mass',mass_finale(ielmc)-mass_inital(ielmc),mass_finale(ielmc),mass_inital(ielmc),'flux',&
        tmpval,-plt_bgcr%RootAutoCO2_pft(NZ)+plt_bgcr%LitrfallBlgrElms_pft(ielmc,NZ)+&
        plt_distb%RootLost2Fire_pft(ielmc,NZ)-plt_bgcr%Xfer2RootsC_pft(NZ),info,id_loc,&
        plt_morph%NumStructuralRootAxes_pft(NZ)    
!    write(iut,*)'plt_biom%SeasonalNonstElms_pft',plt_biom%SeasonalNonstElms_pft(:,NZ)
    write(iut,*)'plt_biom%RootElmsBeg_pft',plt_biom%RootElmsBeg_pft(:,NZ)
      
    if(abs(masserr)>1.e-8_r8)call endrun(trim(mod_filename)//' at line',__LINE__)               
  endif
  end subroutine RootMassBalCheck
  ![tail]
end module PlantDebugMod
