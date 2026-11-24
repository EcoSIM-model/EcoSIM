
module MLDataDiagMod

  use data_kind_mod,        only: r8 => DAT_KIND_R8,yearIJ_type

  use minimathmod,          only: AZMAX1,safe_adb,AZERO,AZERO1,real_truncate,isclose
  use EcoSiMParDataMod,     only: micpar
  use EcosimConst,          only: LtHeatIceMelt,Tref
  use DebugToolMod
  use abortutils,           only: endrun
  use EcoSIMSolverPar  
  use TracerIDMod
  use SoilBGCDataType
  USE PlantDataRateType
  USE SoilWaterDataType
  USE SurfLitterDataType
  USE EcosimBGCFluxType
  USE EcoSIMCtrlDataType
  USE ClimForcDataType
  USE GridDataType
  USE SoilPropertyDataType
  use SoilPhysDataType
  use SOMDataType
  use ChemTranspDataType
  use MicrobialDataType
  use IrrigationDataType
  use SoilHeatDataType
  use MLDataDiagType
  use MicrobialDiagMod
implicit none

  private
  character(len=*), parameter :: mod_filename = &
  __FILE__
  public :: MLDataDiag
  contains

!------------------------------------------------------------------------------------------
  subroutine MLDataDiag(yearIJ,NHW,NHE,NVN,NVS)
  !
  !write data for machine learning 
  !averaged from the top 60 cm
  use EcoSIMCtrlMod, only: etimer
  use abortutils,    only: iulog

  implicit none
  type(yearIJ_type), intent(in) :: yearIJ    
  integer, intent(in) :: NHW,NHE,NVN,NVS
  integer :: NY,NX
  INTEGER :: L,LF, stat
  REAL(r8) :: DPZ,DVOLL
  real(r8) :: DOC_vr(1:JZ)
  real(r8) :: Acet_vr(1:JZ)
  real(r8) :: MicbE(NumPlantChemElms)
  real(r8) :: Fermentor_vr(1:JZ)
  real(r8) :: AcetMGC_vr(1:JZ)
  real(r8) :: H2MGC_vr(1:JZ)
  real(r8) :: AeroMOC_vr(1:JZ)
  real(r8) :: AeroFungC_vr(1:JZ)
  real(r8) :: AeroBactC_vr(1:JZ)

  DOC_vr=0._R8
  Acet_vr=0._R8
  DO NX=NHW,NHE
    DO NY=NVN,NVS  
      !locate the lower boundary 60 cm
      DO L=1,NL_col(NY,NX)
        DVOLL=DLYR_3D(3,L,NY,NX)*AREA_3D(3,NU_col(NY,NX),NY,NX)
        DOC_vr(L)=sum(DOM_MicP_vr(idom_doc,:,L,NY,NX))
        Acet_vr(L)=sum(DOM_MicP_vr(idom_acetate,:,L,NY,NX))
        call SumMicbGroup(L,NY,NX,micpar%mid_fermentor,MicbE);Fermentor_vr(L)=MicbE(ielmc)
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAcetoCH4GenArchea,MicbE);AcetMGC_vr(L)=MicbE(ielmc)
        call SumMicbGroup(L,NY,NX,micpar%mid_AutoH2GenoCH4GenArchea,MicbE,isauto=.true.);H2MGC_vr(L)=MicbE(ielmc)
        call SumMicbGroup(L,NY,NX,micpar%mid_Aerob_Fungi,MicbE); AeroFungC_vr(L) = MicbE(ielmc)
        call SumMicbGroup(L,NY,NX,micpar%mid_HeterAerobBacter,MicbE); AeroBactC_vr(L)=MicbE(ielmc)
        call SumMicbGroup(L,NY,NX,micpar%mid_AutoAeroCH4OxiBacter,MicbE,isauto=.true.); AeroMOC_vr(L)=micBE(ielmc)

        IF(CumSoilThickness_vr(L,NY,NX)>=0.6_r8)THEN
          LF=L
          EXIT
        ENDIF
      ENDDO
          
      
      !obtain top 30 cm average 
      DPZ=0.3_r8
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), TKS_vr(1:LF,NY,NX), DPZ, TEMP30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), THETW_vr(1:LF,NY,NX),DPZ,THETW30cm_col(NY,NX)),__LINE__) !volume metric H2O 
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trc_solcl_vr(idg_O2,1:LF,NY,NX),DPZ,O2wConc30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trc_solcl_vr(idg_CO2,1:LF,NY,NX),DPZ,CO2wConc30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trc_solcl_vr(idg_H2,1:LF,NY,NX),DPZ,H2wConc30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trc_solcl_vr(idg_CH4,1:LF,NY,NX),DPZ,CH4wConc30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), DOC_vr(1:LF),DPZ,DOCConc30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), Acet_vr(1:LF),DPZ,AcetConc30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), RCH4ProdAcetcl_vr(1:LF,NY,NX),DPZ,RCH4ProdAcet30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), RCH4ProdHydrog_vr(1:LF,NY,NX),DPZ,RCH4ProdHG30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), RCH4Oxi_aero_vr(1:LF,NY,NX),DPZ,RAeroCH4Oxi30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trcs_RMicbUptake_vr(idg_CO2,1:LF,NY,NX),DPZ,RCO2Ht30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), RFerment_vr(1:LF,NY,NX),DPZ,RFerment30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), AcetMGC_vr(1:LF),DPZ,AcetMGC30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), H2MGC_vr(1:LF),DPZ,H2MGC30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), Fermentor_vr(1:LF),DPZ,FermC30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), AeroMOC_vr(1:LF),DPZ,AeroMOC30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), AeroFungC_vr(1:LF),DPZ,AeroHRFungC30cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), AeroBactC_vr(1:LF),DPZ,AeroHRBactC30cm_col(NY,NX)),__LINE__)  

      !obtain top 30 cm average 
      DPZ=0.6_r8
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), TKS_vr(1:LF,NY,NX), DPZ, TEMP60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), THETW_vr(1:LF,NY,NX),DPZ,THETW60cm_col(NY,NX)),__LINE__) !volume metric H2O 
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trc_solcl_vr(idg_O2,1:LF,NY,NX),DPZ,O2wConc60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trc_solcl_vr(idg_CO2,1:LF,NY,NX),DPZ,CO2wConc60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trc_solcl_vr(idg_H2,1:LF,NY,NX),DPZ,H2wConc60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fcon_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trc_solcl_vr(idg_CH4,1:LF,NY,NX),DPZ,CH4wConc60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), DOC_vr(1:LF),DPZ,DOCConc60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), Acet_vr(1:LF),DPZ,AcetConc60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), RCH4ProdAcetcl_vr(1:LF,NY,NX),DPZ,RCH4ProdAcet60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), RCH4ProdHydrog_vr(1:LF,NY,NX),DPZ,RCH4ProdHG60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), RCH4Oxi_aero_vr(1:LF,NY,NX),DPZ,RAeroCH4Oxi60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), trcs_RMicbUptake_vr(idg_CO2,1:LF,NY,NX),DPZ,RCO2Ht60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), RFerment_vr(1:LF,NY,NX),DPZ,RFerment60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), AcetMGC_vr(1:LF),DPZ,AcetMGC60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), H2MGC_vr(1:LF),DPZ,H2MGC60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), Fermentor_vr(1:LF),DPZ,FermC60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), AeroMOC_vr(1:LF),DPZ,AeroMOC60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), AeroFungC_vr(1:LF),DPZ,AeroHRFungC60cm_col(NY,NX)),__LINE__)
      call check_interp(compute_fmas_average(LF, CumSoilThickness_vr(1:LF,NY,NX), AeroBactC_vr(1:LF),DPZ,AeroHRBactC60cm_col(NY,NX)),__LINE__)  

    ENDDO
  ENDDO
  end subroutine MLDataDiag

!------------------------------------------------------------------------------------------

  function compute_fcon_average(n, x, f, x_eval, f_avg)result(stat)
  !! Computes the average of f(x) from 0 to x_eval using linear interpolation.
  !! Assumes the function starts at (0.0, 0.0) implicitly.
  !! assuming f is given as intensity
  integer, intent(in)  :: n          ! Number of strictly positive data points
  real(r8), intent(in) :: x(n)       ! x values: x(1) to x(n)
  real(r8), intent(in) :: f(n)       ! f(x) values: f(1) to f(n)
  real(r8), intent(in) :: x_eval     ! The upper limit of the average
  real(r8), intent(out):: f_avg      ! The computed average
  integer :: stat       ! Status flag: 0=Success, -1=Negative x, 1=Out of bounds

  integer :: i
  real(r8) :: integral
        
  ! Initialize defaults
  stat = 0
  f_avg = 0.0_r8
  integral = 0.0_r8

  ! Handle edge cases
  if (x_eval < 0.0_r8) then
    stat = -1  ! Error: Cannot evaluate average for negative x
    return
  end if
    
  if (x_eval == 0.0_r8) then
    f_avg = 0.0_r8  ! Limit of f(x) as x->0 is f(0) = 0
    return
  end if

  if (x_eval > x(n)) then
    stat = 1   ! Error: x_eval exceeds provided data bounds
    return
  end if

  integral=0._r8
  !within the first grid
  if(x_eval <= x(1))then
    f_avg=f(1)
    return
  else
    integral = f(1) * x(1)
  endif        

  ! Iterate through the intervals to compute the trapezoidal integral
  do i = 2, n    
    if (x_eval <= x(i)) then
      integral = integral + f(i) * (x_eval - x(i-1))        
      ! Exit loop as we have reached the evaluation point
      exit
    else
      ! Add the full trapezoidal area for this interval      
      integral = integral + f(i) * (x(i)- x(i-1))
    end if
  end do

  ! Compute the final average
  f_avg = integral / x_eval
    
  end function compute_fcon_average

!------------------------------------------------------------------------------------------

  function compute_fmas_average(n, x, f, x_eval, f_avg)result(stat)
  !! Computes the average of f(x) from 0 to x_eval using linear interpolation.
  !! Assumes the function starts at (0.0, 0.0) implicitly.
  !! assuming f is given as intensity
  integer, intent(in)  :: n          ! Number of strictly positive data points
  real(r8), intent(in) :: x(n)       ! x values: x(1) to x(n)
  real(r8), intent(in) :: f(n)       ! f(x) values: f(1) to f(n)
  real(r8), intent(in) :: x_eval     ! The upper limit of the average
  real(r8), intent(out):: f_avg      ! The computed average
  integer :: stat       ! Status flag: 0=Success, -1=Negative x, 1=Out of bounds

  integer :: i
  real(r8) :: integral
        
  ! Initialize defaults
  stat = 0
  f_avg = 0.0_r8
  integral = 0.0_r8

  ! Handle edge cases
  if (x_eval < 0.0_r8) then
    stat = -1  ! Error: Cannot evaluate average for negative x
    return
  end if
    
  if (x_eval == 0.0_r8) then
    f_avg = 0.0_r8  ! Limit of f(x) as x->0 is f(0) = 0
    return
  end if

  if (x_eval > x(n)) then
    stat = 1   ! Error: x_eval exceeds provided data bounds
    return
  end if

  !within the first grid
  if(x_eval <= x(1))then
    f_avg=f(1)
    return
  endif        

  ! Iterate through the intervals to compute the trapezoidal integral
  integral=0._r8
  do i = 1, n    
    if (x_eval <= x(i)) then
      integral = integral + f(i) * (x_eval - x(i-1))/(x(i)- x(i-1))        
      ! Exit loop as we have reached the evaluation point
      exit
    else
      ! Add the full trapezoidal area for this interval
      integral = integral + f(i)
    end if
  end do

  ! Compute the final average
  f_avg = integral / x_eval
    
  end function compute_fmas_average

!------------------------------------------------------------------------------------------
  subroutine check_interp(stat,line)
  implicit none
  integer, intent(in) :: stat
  integer, intent(in) :: line

  if (stat==-1)then
    call endrun('Error: Cannot evaluate average for negative x'//mod_filename,line)
  elseif(stat==1)then
    call endrun('Error: x_eval exceeds provided data bounds'//mod_filename,line)
  endif

  end subroutine check_interp
end module MLDataDiagMod