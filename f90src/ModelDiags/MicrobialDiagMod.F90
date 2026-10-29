module MicrobialDiagMod
!!
! DESCRIPTION:
! codes to do soil biological transformations
!
! USES:
  use data_kind_mod,    only: r8 => DAT_KIND_R8
  use abortutils,       only: endrun
  use minimathmod,      only: safe_adb, AZMAX1,AZERO
  use EcoSiMParDataMod, only: micpar
  use EcoSIMConfig , only : ndbiomcp => NumDeadMicrbCompts,NumMicbAFunGrupsPerCmplx,NumMicbHFunGrupsPerCmplx        
  use DebugToolMod
  use SoilWaterDataType
  use SurfLitterDataType
  use MicrobialDataType
  use NitroPars
  use SOMDataType
  use ChemTranspDataType
  use GridConsts
  use SoilBGCDataType
  use EcoSIMCtrlDataType
  use SoilPhysDataType
  use SoilPropertyDataType
  use GridDataType
  
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
  __FILE__

  public :: sumORGMLayL
  public :: sumLitrOMLayL
  public :: sumSurfOMCK
  public :: sumMicBiomLayL
  public :: sumHumOMLayL
  public :: SumMicbGroup
  public :: sumDOML
  public :: SumSolidOML

  contains
!------------------------------------------------------------------------------------------
  subroutine sumORGMLayL(L,NY,NX,ORGM,conly,info,clip_small_pools)
  !
  !sum up organic matter in layer L
  !including live microbial biomass, microbial residue, sorbed dom+acetate, and solid SOM
  
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  logical, optional, intent(in) :: conly
  character(len=*), optional, intent(in) :: info
  integer :: K,N,NGL,M,MID,NE,jcplx1,nelms,idom
  logical, optional, intent(in) :: clip_small_pools !default preserves the existing tiny-pool cleanup
  logical :: do_clip_small_pools
  logical :: conly_loc
  real(r8) :: orgm0(1:NumPlantChemElms)

  if(present(conly))then
    conly_loc=conly
  else
    conly_loc=.false.
  endif

  do_clip_small_pools=.true.
  if(present(clip_small_pools))do_clip_small_pools=clip_small_pools

  ORGM=0._r8
  !sumup heterotrophic microbes
  if(L==0)then
    jcplx1=micpar%NumOfLitrCmplxs
  else
    jcplx1=jcplx
  endif
  if(conly_loc)then
    nelms=1
  else
    nelms=NumPlantChemElms
  endif  
  
  !add autotrophic microbes
  DO  N=1,NumMicbAFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,nelms
          ORGM(NE)=ORGM(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      enddo
    enddo
  enddo

  DK100: DO K=1,jcplx1
    !add heterotrophic microbes
    DC100: DO  N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,nelms
            ORGM(NE)=ORGM(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo
    enddo DC100
  ENDDO DK100

  DK200: DO K=1,jcplx1
    !add microbial residual
    DO  M=1,ndbiomcp
      DO NE=1,nelms
        ORGM(NE)=AZERO(ORGM(NE)+OMBioResdu_vr(NE,M,K,L,NY,NX))
        if(ORGM(NE)<0._r8)then
        print*,'orgm3',ORGM(NE),OMBioResdu_vr(NE,M,K,L,NY,NX),L
        stop
        endif
      ENDDO    
    ENDDO

    !add dom
    if(do_clip_small_pools)then
      DO idom=idom_beg,idom_end
        if(abs(DOM_MicP_vr(idom,K,L,NY,NX))<1.e-12_r8)DOM_MicP_vr(idom,K,L,NY,NX)=0._r8
        if(abs(SorbedOM_vr(idom,K,L,NY,NX))<1.e-12_r8)SorbedOM_vr(idom,K,L,NY,NX)=0._r8
        if(abs(DOM_MacP_vr(idom,K,L,NY,NX))<1.e-12_r8)DOM_MacP_vr(idom,K,L,NY,NX)=0._r8
      ENDDO
    endif


    DO NE=1,nelms      
      ORGM(NE)=ORGM(NE)+DOM_MicP_vr(NE,K,L,NY,NX)+DOM_MacP_vr(NE,K,L,NY,NX)+SorbedOM_vr(NE,K,L,NY,NX)
        if(ORGM(NE)<0._r8)then
        print*,'orgm2',NE,ORGM,DOM_MicP_vr(NE,K,L,NY,NX),DOM_MacP_vr(NE,K,L,NY,NX),SorbedOM_vr(NE,K,L,NY,NX)
        stop
        endif
    ENDDO

    ORGM(ielmc)=ORGM(ielmc)+DOM_MicP_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    if(ORGM(ielmc)<0._r8)then
    print*,'orgmxxx',ORGM,DOM_MicP_vr(idom_acetate,K,L,NY,NX),DOM_MacP_vr(idom_acetate,K,L,NY,NX),SorbedOM_vr(idom_acetate,K,L,NY,NX)    
    print*,'L=',L
    if(present(info))print*,info
    stop
    endif

    !add solid organic matter, litter, manure etc
    DO  M=1,jsken
      DO NE=1,nelms
        orgm0(NE)=ORGM(NE)
        ORGM(NE)=ORGM(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
        if(ORGM(NE)<0._r8)then
        print*,K,NE,nelms,trim(info)//' orgm1 in '//trim(mod_filename)
        print*,'lay ',L,ORGM,SolidOM_vr(:,M,K,L,NY,NX),orgm0
        stop
        endif
      ENDDO  
    ENDDO  
  ENDDO DK200   

  end subroutine sumORGMLayL
!------------------------------------------------------------------------------------------
  subroutine SumSolidOML(ielm,L,NY,NX,SOMe)

  implicit none
  integer,intent(in) :: ielm
  integer,intent(in) :: L,NY,NX
  real(r8),intent(out) :: SOMe(1:jcplx)
  integer :: K,M,KL

  SOMe=0._r8
  if(L==0)then
    KL=micpar%NumOfLitrCmplxs
  else
    KL=jcplx
  endif  
  
  DO K=1,KL
    DO  M=1,jsken
      SOMe(K)=SOMe(K)+SolidOM_vr(ielm,M,K,L,NY,NX)
    ENDDO
  ENDDO
  
  end subroutine SumSolidOML
!------------------------------------------------------------------------------------------

  subroutine sumLitrOMLayL(L,NY,NX,ORGM,I,J)
  !
  !sum up litter OM (DOM + litter OM) in layer L
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  integer , optional, intent(in) :: I,J  
  integer :: K,N,NGL,M,MID,NE,idom
  real(r8) :: DOM_micp(idom_beg:idom_end)
  real(r8) :: DOM_macp(idom_beg:idom_end) 
  real(r8) :: BiomA(NumPlantChemElms),BiomH(NumPlantChemElms),BiomDead(NumPlantChemElms)
  real(r8) :: OMSorb(NumPlantChemElms), OMSolid(NumPlantChemElms)

  ORGM     = 0._r8
  DOM_micp = 0._r8
  DOM_macp = 0._r8
  BiomA    = 0._r8
  BiomH    = 0._r8
  BiomDead = 0._r8
  OMSolid  = 0._r8
  OMSorb   = 0._r8
  ! Attribute complex-independent autotrophic biomass to litter only.
  DO  N=1,NumMicbAFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          BiomA(NE)=BiomA(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      enddo
    enddo
  enddo

  DO K=1,micpar%NumOfLitrCmplxs
    !add live heterotrophic microbes
    DO  N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            BiomH(NE)=BiomH(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo
    enddo
    
   !add microbial residual
    DO  M=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        BiomDead(NE)=BiomDead(NE)+OMBioResdu_vr(NE,M,K,L,NY,NX)
      ENDDO    
    ENDDO

    !add dom
    DO idom=idom_beg,idom_end
      DOM_micp(idom) = DOM_micp(idom) + DOM_MicP_vr(idom,K,L,NY,NX)
      DOM_macp(idom) = DOM_macp(idom) + DOM_MacP_vr(idom,K,L,NY,NX)
    ENDDO

    !add dom
    DO NE=1,NumPlantChemElms
      OMSorb(NE)=OMSorb(NE)+SorbedOM_vr(NE,K,L,NY,NX)
    ENDDO

    OMSorb(ielmc)=OMSorb(ielmc)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    !add solid organic matter
    DO  M=1,jsken
      DO NE=1,NumPlantChemElms
        OMSolid(NE)=OMSolid(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
      ENDDO  
    ENDDO  
  ENDDO   
  
  DO NE=1,NumPlantChemElms
    ORGM(NE)=BiomA(NE)+BiomH(NE)+BiomDead(NE)+OMSorb(NE)+OMSolid(NE)+DOM_micp(NE)+DOM_macp(NE)
  ENDDO
  ORGM(ielmc)=ORGM(ielmc)+DOM_micp(idom_acetate) + DOM_macp(idom_acetate)
  
  end subroutine sumLitrOMLayL

!------------------------------------------------------------------------------------------

  subroutine sumHumOMLayL(L,NY,NX,ORGM)
  !
  !sum up humus OM in layer L; autotrophic biomass is included in litter
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)
  integer :: K,N,NGL,M,MID,NE,idom
  real(r8) :: DOM_micp(idom_beg:idom_end)
  real(r8) :: DOM_macp(idom_beg:idom_end) 

  ORGM=0._r8
  DOM_micp=0._r8
  DOM_macp=0._r8

  DO  K=micpar%NumOfLitrCmplxs+1,jcplx
   !sumup heterotrophic microbes
    DO  N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            ORGM(NE)=ORGM(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo
    enddo
  !add microbial residual
    DO  M=1,ndbiomcp
      DO NE=1,NumPlantChemElms
        ORGM(NE)=ORGM(NE)+OMBioResdu_vr(NE,M,K,L,NY,NX)
      ENDDO    
    ENDDO

    !add dom
    DO idom=idom_beg,idom_end
      DOM_micp(idom) = DOM_micp(idom) + DOM_MicP_vr(idom,K,L,NY,NX)
      DOM_macp(idom) = DOM_macp(idom) + DOM_MacP_vr(idom,K,L,NY,NX)
    ENDDO
    DO NE=1,NumPlantChemElms
      ORGM(NE)=ORGM(NE)+SorbedOM_vr(NE,K,L,NY,NX)
    ENDDO
    !add acetate
    ORGM(ielmc)=ORGM(ielmc)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    !add solid organic matter
    DO  M=1,jsken
      DO NE=1,NumPlantChemElms
        ORGM(NE)=ORGM(NE)+SolidOM_vr(NE,M,K,L,NY,NX)
      ENDDO  
    ENDDO  
  ENDDO    

  ORGM=ORGM+DOM_micp(1:NumPlantChemElms) + DOM_macp(1:NumPlantChemElms)

  ORGM(ielmc)=ORGM(ielmc) +DOM_micp(idom_acetate) + DOM_macp(idom_acetate)

  
  end subroutine sumHumOMLayL
  

!------------------------------------------------------------------------------------------

  subroutine sumMicBiomLayL(L,NY,NX,ORGM,I,J)
  !
  !sum up litter OM in layer L
  implicit none
  integer, intent(in) :: L, NY,NX
  real(r8), intent(out) :: ORGM(1:NumPlantChemElms)  !microbial biomass 
  integer,optional, intent(in) :: I,J
  character(len=*), parameter :: subname='sumMicBiomLayL'
  integer :: K,N,NGL,M,MID,NE,jcplx1
  real(r8) :: BiomHK(NumPlantChemElms,jcplx)

  call PrintInfo('beg '//subname)
  ORGM   = 0._r8
  BiomHK = 0._r8
  if(L==0)then
    jcplx1=micpar%NumOfLitrCmplxs
  else
    jcplx1=jcplx
  endif

  !add autotrophic microbes
  DO  N=1,NumMicbAFunGrupsPerCmplx
    DO NGL=JGniA(N),JGnfA(N)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          ORGM(NE)=ORGM(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      enddo
    enddo
  enddo

  !add heterotrophs
  DO K=1,jcplx1
    !add heterotrophic microbes
    DO  N=1,NumMicbHFunGrupsPerCmplx
      DO NGL=JGniH(N),JGnfH(N)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            BiomHK(NE,K)=BiomHK(NE,K)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo
    enddo
    DO NE=1,NumPlantChemElms
      ORGM(NE)=ORGM(NE)+BiomHK(NE,K)    
    enddo  
  enddo  

  call PrintInfo('end '//subname)
  end subroutine sumMicBiomLayL
!------------------------------------------------------------------------------------------

  subroutine sumSurfOMCK(NY,NX,SOMHeterKC,SOMAutorC)
  implicit none
  integer, intent(in) :: NY,NX
  real(r8), intent(out) :: SOMHeterKC(1:micpar%NumOfLitrCmplxs)  !total organic C in each litter complex
  real(r8), intent(out) :: SOMAutorC
  integer :: K,N,NGL,M,MID,NE,L

  SOMHeterKC = 0._r8
  SOMAutorC  = 0._r8

  L  = 0
  NE=ielmc
  !autotrophs
  DO  N=1,NumMicbAFunGrupsPerCmplx
    do NGL=JGniA(n),JGnfA(n)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        SOMAutorC=SOMAutorC+mBiomeAutor_vr(NE,MID,L,NY,NX)        
      ENDDO
    ENDDO
  enddo
  !live microbes
  DO K=1,micpar%NumOfLitrCmplxs
    DO  N=1,NumMicbHFunGrupsPerCmplx
      do NGL=JGniH(n),JGnfH(n)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          SOMHeterKC(K)=SOMHeterKC(K)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)          
        enddo
      enddo  
    enddo    

    !microbial residue
    DO  M=1,ndbiomcp
      SOMHeterKC(K)=SOMHeterKC(K)+OMBioResdu_vr(NE,M,K,L,NY,NX)              
    ENDDO  

    !add dom
    SOMHeterKC(K)=SOMHeterKC(K)+DOM_MicP_vr(NE,K,L,NY,NX)+DOM_MacP_vr(NE,K,L,NY,NX)+SorbedOM_vr(NE,K,L,NY,NX)
    
    !add acetate
    SOMHeterKC(K)=SOMHeterKC(K)+DOM_MicP_vr(idom_acetate,K,L,NY,NX)+DOM_MacP_vr(idom_acetate,K,L,NY,NX)+SorbedOM_vr(idom_acetate,K,L,NY,NX)    

    !add solid om
    DO M=1,jsken
      SOMHeterKC(K)=SOMHeterKC(K)+SolidOM_vr(NE,M,K,L,NY,NX)      
    ENDDO  
  ENDDO
  
  end subroutine sumSurfOMCK
!------------------------------------------------------------------------------------------

  subroutine SumMicbGroup(L,NY,NX,igroup,MicbE,isauto)
  !
  !Description
  !
  implicit none
  integer,  intent(in) :: L,NY,NX
  integer,  intent(in) :: igroup
  real(r8), intent(out):: MicbE(1:NumPlantChemElms)
  logical, optional, intent(in) :: isauto
  logical :: isauto_loc
  integer :: K,NE,M,MID,NGL
  character(len=*), parameter :: subname='SumMicbGroup'

  call PrintInfo('beg '//subname)
  if(present(isauto))then
    isauto_loc=isauto
  else
    isauto_loc=.false.
  endif  

  micBE=0._r8
  associate(mid_ANME2d => micpar%mid_AutoAMOANME2D, mid_NC10 => micpar%mid_AutoAMONC10)
  if(isauto_loc)then
    if(igroup /= micpar%mid_AutoAmmoniaOxidBacter        .and. &
       igroup /= micpar%mid_AutoNitriteOxidBacter        .and. & 
       igroup /= micpar%mid_AutoAeroCH4OxiBacter .and. &
       igroup /= micpar%mid_AutoH2GenoCH4GenArchea .and. &
       igroup /= mid_ANME2d .and. igroup /= mid_NC10) then
      call endrun('undefined autotroph group in '//trim(mod_filename),__LINE__)
    endif

    do NGL=JGniA(igroup),JGnfA(igroup)
      DO  M=1,nlbiomcp
        MID=micpar%get_micb_id(M,NGL)
        DO NE=1,NumPlantChemElms
          micBE(NE)=micBE(NE)+mBiomeAutor_vr(NE,MID,L,NY,NX)
        ENDDO
      ENDDO
    ENDDO

  else
    if(igroup /= micpar%mid_HeterAerobBacter        .and. &
       igroup /= micpar%mid_Facult_DenitBacter      .and. &
       igroup /= micpar%mid_Aerob_Fungi             .and. &
       igroup /= micpar%mid_fermentor               .and. &
       igroup /= micpar%mid_HeterAcetoCH4GenArchea  .and. &
       igroup /= micpar%mid_HeterAerobN2Fixer       .and. &
       igroup /= micpar%mid_HeterAnaerobN2Fixer     .and. &
       igroup /= micpar%mid_HeterMixtCynoBacter) then
      call endrun('undefined heterotroph group in '//trim(mod_filename),__LINE__)
    endif
    
    DO K=1,micpar%NumOfLitrCmplxs    
      do NGL=JGniH(igroup),JGnfH(igroup)
        DO  M=1,nlbiomcp
          MID=micpar%get_micb_id(M,NGL)
          DO NE=1,NumPlantChemElms
            micBE(NE)=micBE(NE)+mBiomeHeter_vr(NE,MID,K,L,NY,NX)
          ENDDO
        enddo
      enddo  
    enddo    
  endif

  end associate
  call PrintInfo('end '//subname)
  end subroutine SumMicbGroup

!------------------------------------------------------------------------------------------
  subroutine sumDOML(L,NY,NX,DOM)

  implicit none
  integer,  intent(in) :: L,NY,NX
  real(r8), intent(out):: DOM(idom_beg:idom_end)
  character(len=*), parameter :: subname='sumDOML'
  integer :: idom, K
  real(r8) :: DOM_micp(idom_beg:idom_end)
  real(r8) :: DOM_macp(idom_beg:idom_end) 
  call PrintInfo('beg '//subname)
  DOM      = 0._r8
  DOM_micp = 0._r8
  DOM_macp = 0._r8
  DO K = 1, micpar%NumOfLitrCmplxs
    DO idom=idom_beg,idom_end
      DOM_micp(idom) = DOM_micp(idom)+DOM_MicP_vr(idom,K,L,NY,NX)
      DOM_macp(idom) = DOM_macp(idom)+DOM_MacP_vr(idom,K,L,NY,NX)
    ENDDO
  ENDDO
  DOM=DOM_micp+DOM_macp
  call PrintInfo('end '//subname)
  end subroutine sumDOML

end module MicrobialDiagMod
