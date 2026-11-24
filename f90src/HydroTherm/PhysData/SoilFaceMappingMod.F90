module SoilFaceMappingMod
  ! Horizontal face storage uses the owner's layer and the destination column.
  ! The owner is the west/north cell; a negative flux reverses the actual donor.
  ! Vertical and external boundary faces retain their original layer indices.
  use data_kind_mod, only: r8 => DAT_KIND_R8
  use GridConsts, only: JZ,JY,JX
  use ElmIDMod, only: iWestEastDirection,iNorthSouthDirection,iVerticalDirection
  use GridDataType, only: NL_col,FlowDirIndicator_col,DLYR_3D,AREA_3D
  use SoilPropertyDataType, only: VLSoilPoreMicP_vr
  use EcoSIMCtrlDataType, only: ZEROS2
  use DebugToolMod, only: PrintInfo
  implicit none
  private
  integer, allocatable :: target_layer(:,:,:,:)
  integer :: west_edge,north_edge
  public :: GetSoilFaceGeometry
  public :: RebuildSoilFaceMap,DestroySoilFaceMap,SoilFaceSlot,SoilFaceIncomingRange
contains
!------------------------------------------------------------------------------------------
  subroutine RebuildSoilFaceMap(NHW,NHE,NVN,NVS)
    integer, intent(in) :: NHW,NHE,NVN,NVS
    integer :: NX,NY,L,N,TX,TY,LL
    character(len=*), parameter :: subname='RebuildSoilFaceMap'

    call PrintInfo('beg '//subname)
    if(.not.allocated(target_layer))allocate(target_layer(2,JZ,JY,JX))
    target_layer=0
    west_edge=NHW
    north_edge=NVN
    do NX=NHW,NHE
      do NY=NVN,NVS
        do N=FlowDirIndicator_col(NY,NX),2
          TX=NX;TY=NY
          if(N==iWestEastDirection)TX=NX+1
          if(N==iNorthSouthDirection)TY=NY+1
          if(TX>NHE.or.TY>NVS)cycle
          do L=1,NL_col(NY,NX)
            if(VLSoilPoreMicP_vr(L,NY,NX)<=ZEROS2(NY,NX))cycle
            do LL=L,NL_col(TY,TX)
              if(VLSoilPoreMicP_vr(LL,TY,TX)<=ZEROS2(TY,TX))cycle
              target_layer(N,L,NY,NX)=LL
              exit
            enddo
          enddo
        enddo
      enddo
    enddo
    call PrintInfo('end '//subname)
  end subroutine RebuildSoilFaceMap

!------------------------------------------------------------------------------------------
  subroutine DestroySoilFaceMap
    character(len=*), parameter :: subname='DestroySoilFaceMap'

    call PrintInfo('beg '//subname)
    if(allocated(target_layer))deallocate(target_layer)

    call PrintInfo('end '//subname)
  end subroutine DestroySoilFaceMap

!------------------------------------------------------------------------------------------
  pure integer function SoilFaceSlot(N,source_layer,destination_layer) result(slot)
    integer, intent(in) :: N,source_layer,destination_layer
    slot=destination_layer
    if(N/=3)slot=source_layer
  end function SoilFaceSlot

!------------------------------------------------------------------------------------------
  subroutine SoilFaceIncomingRange(N,L,NY,NX,first,last)
    ! All owners mapping to an active target form a contiguous layer range;
    ! inactive owner slots inside that range have zero flux. An empty range
    ! sums to zero. West/north external faces remain at the local layer L.
    integer, intent(in) :: N,L,NY,NX
    integer, intent(out) :: first,last
    integer :: SX,SY,K
    character(len=*), parameter :: subname='SoilFaceIncomingRange'

    call PrintInfo('beg '//subname)
    first=L;last=L
    if(N==iVerticalDirection .or. (N==iWestEastDirection.and.NX==west_edge) .or. &
      (N==iNorthSouthDirection .and. NY==north_edge))then
      call PrintInfo('end '//subname)
      return
    endif
    SX=NX;SY=NY
    if(N==iWestEastDirection)SX=NX-1
    if(N==iNorthSouthDirection)SY=NY-1
    first=1;last=0
    do K=1,NL_col(SY,SX)
      if(target_layer(N,K,SY,SX)/=L)cycle
      if(last==0)first=K
      last=K
    enddo
    call PrintInfo('end '//subname)
  end subroutine SoilFaceIncomingRange
!------------------------------------------------------------------------------------------
  subroutine GetSoilFaceGeometry(N,N1,N2,N3,N4,N5,N6,FaceArea,AreaOverDistance,Dispersivity)
    ! Use current dimensions of the actual layer pair. The face owner supplies
    ! the area, matching the convention used by water and heat transport.
    integer, intent(in) :: N,N1,N2,N3,N4,N5,N6
    real(r8), intent(out) :: FaceArea          ![m2]
    real(r8), intent(out) :: AreaOverDistance  ![m]
    real(r8), intent(out) :: Dispersivity      ![m]
    real(r8) :: Distance
    character(len=*), parameter :: subname='GetSoilFaceGeometry'

    call PrintInfo('beg '//subname)
    FaceArea=AREA_3D(N,N3,N2,N1)
    Distance=0.5_r8*(DLYR_3D(N,N3,N2,N1)+DLYR_3D(N,N6,N5,N4))
    AreaOverDistance=0._r8
    Dispersivity=0._r8
    if(Distance>0._r8)then
      AreaOverDistance=FaceArea/Distance
      Dispersivity=0.20_r8*Distance**1.07_r8
    endif
    call PrintInfo('end '//subname)
  end subroutine GetSoilFaceGeometry

!------------------------------------------------------------------------------------------
end module SoilFaceMappingMod
