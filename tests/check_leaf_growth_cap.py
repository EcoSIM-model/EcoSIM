"""Exercise production capped growth and canopy bookkeeping with occupied upper nodes.

Compile the actual Fortran procedures against small state fixtures. Check that a
retreating reproductive cap blocks growth without losing existing area or C/N/P,
that empty growth windows conserve allocation, and that reused slots are protected.
"""
import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def procedure(text, name):
    pattern = (r'^  (?:(?:pure|integer)\s+)?(?:subroutine|function) ' + name +
               r'\b.*?^  end (?:subroutine|function) ' + name + r'\b[^\n]*')
    return re.search(pattern, text, re.M | re.S | re.I)[0]


def fixture(branch, math, phenology):
    names = ['LimitNodeGrowthAllocation', 'GrowLeavesOnBranch',
             'GrowPetolShethOnBranch', 'GrowStalkOnBranch', 'AllocateLeaf2CanopyLayers']
    routines = '\n'.join(procedure(branch, name) for name in names)
    routines += '\n' + procedure(phenology, 'Emerged_plant_Phenology')
    files = dict(plt_pheno='PlantPhenologyAPIData', plt_morph='PlantMorphologyAPIData',
                 plt_biom='PlantBiomassAPIData', plt_allom='PlantAllometryAPIData',
                 plt_site='PlantSiteAPIData', plt_rad='PlantRadiationAPIData',
                 plt_ew='PlantEnergyWaterAPIData')
    types = []
    # Derive the fixture's field types and bounds from the production API allocation.
    for instance in sorted(set(re.findall(r'(plt_\w+)%', routines))):
        source = (ROOT / 'f90src/APIData' / (files[instance] + '.F90')).read_text()
        fields = []
        for field in sorted(set(re.findall(instance + r'%(\w+)', routines))):
            decl = re.search(r'^.*::\s*' + field + r'\b[^\n]*', source, re.M | re.I)[0]
            kind = 'integer' if re.search(r'^\s*integer', decl, re.I) else 'real(r8)'
            allocation = re.search(r'allocate\(this%' + field + r'\((.*?)\)\)', source, re.I)
            dims = '(' + allocation[1] + ')' if allocation else ''
            fields.append(f'  {kind} :: {field}{dims}=0')
        types += [f'type {instance}_type', *fields, 'end type',
                  f'type({instance}_type) :: {instance}']
    mini = (ROOT / 'f90src/Minimath/MiniMathMod.F90').read_text()
    helpers = '\n'.join(procedure(mini, name) for name in ['pMOD', 'AZMAX1_s', 'AZMAX1_d', 'isclose'])
    helpers += '\n' + '\n'.join(procedure(math, name) for name in
                               ['LeafDevelopmentNode', 'LeafNodeGrowthRanges', 'is_plant_woody_vascular'])
    ids = '\n'.join(re.findall(r'^  integer, parameter :: ibrch_.*$', branch, re.M))
    return '''module fixture
use iso_fortran_env, only: r8 => real64
use ElmIDMod
implicit none
integer, parameter :: JP1=1, MaxNumBranches=1, MaxNodesPerBranch1=25
integer, parameter :: NumCanopyLayers1=5, NumLeafInclinationClasses1=3, NumOfPlantMorphUnits=7
real(r8), parameter :: ppmc=1.e-6_r8, tiny_val=1.e-12_r8, PICON=3.141592653589793_r8
real(r8), parameter :: SLA2=-.33_r8, SSL2=-.5_r8, SNL2=-.67_r8
real(r8), parameter :: PSIMin4LeafOff(0:4)=[-200._r8,-2._r8,-2._r8,-2._r8,-2._r8]
real(r8), parameter :: FSTK=.05_r8, ZSTX=.001_r8
! The production procedure accepts this unused calendar argument.
type yearIJ_type
  integer :: I=1, J=1
end type
interface AZMAX1
  module procedure AZMAX1_s, AZMAX1_d
end interface
''' + ids + '\n' + '\n'.join(types) + '''
contains
subroutine PrintInfo(message)
  character(len=*), intent(in) :: message
end subroutine
subroutine DebugPrint(message,index)
  character(len=*), intent(in) :: message
  integer, intent(in) :: index
end subroutine
subroutine CalcPhenolEnvfactor(I,J,NZ,TFNP,WFNG,OFNG)
  integer, intent(in) :: I,J,NZ
  real(r8), intent(out) :: TFNP,WFNG,OFNG
  TFNP=1._r8; WFNG=1._r8; OFNG=1._r8
end subroutine
subroutine live_branch_phenology(I,J,NB,NZ,TFNP,WFNG,OFNG)
  integer, intent(in) :: I,J,NB,NZ
  real(r8), intent(in) :: TFNP,WFNG,OFNG
end subroutine
subroutine branch_specific_phenology(I,J,NB,NZ)
  integer, intent(in) :: I,J,NB,NZ
end subroutine
''' + helpers + '\n' + routines + '\nend module\n'


DRIVER = '''program check_leaf_growth_cap
use fixture
implicit none
integer :: lf,ll,sf,sl,k,hi,nco,lower,expected_leaf,expected_stalk
real(r8) :: part(7), original(7), before(3), after(3), upper(3,2)
real(r8) :: growth(3), alloc, stalk_tip
real(r8) :: layer_before(3,5,0:25,1,1), area_layers_before(5,0:25,1,1)
type(yearIJ_type) :: date
! Exhaustively compare ordinary, non-retreating windows with the previous code.
do hi=1,80
  do nco=1,25
    do lower=0,1
      call LeafNodeGrowthRanges(real(hi-1,r8)+.2_r8,0._r8,hi,nco,25,lower,lower,lf,ll,sf,sl)
      expected_leaf=max(lower,hi-nco+1)
      expected_stalk=max(min(lower,max(lower,hi-nco)),hi-23)
      if(lf/=expected_leaf.or.ll/=hi.or.sf/=expected_stalk.or.sl/=hi)error stop 'ordinary window changed'
    enddo
  enddo
enddo
! The final allowed index retains ECOSYS's INT(...)+1 convention.
if(LeafDevelopmentNode(12.4_r8,10.8_r8)/=11)error stop 'cap convention'
if(LeafDevelopmentNode(12.4_r8,0._r8)/=13)error stop 'unset cap'
if(LeafDevelopmentNode(12.4_r8,1.e-6_r8)/=13)error stop 'threshold'
if(LeafDevelopmentNode(12.4_r8,10._r8)/=11)error stop 'integer cap'
call LeafNodeGrowthRanges(12.4_r8,10.8_r8,13,3,25,1,1,lf,ll,sf,sl)
if(lf/=11.or.ll/=11.or.sl/=11)error stop 'cap intersection'
! A smaller cap cannot move the growth window back onto already completed leaves.
call LeafNodeGrowthRanges(12.4_r8,9.8_r8,13,3,25,1,1,lf,ll,sf,sl)
if(ll>=lf)error stop 'old leaves reopened'
original=[.3_r8,.2_r8,.2_r8,.1_r8,.05_r8,.05_r8,.1_r8]
part=original
call LimitNodeGrowthAllocation(part,lf,ll,sf,sl)
if(any(part(1:2)/=0._r8).or.abs(part(4)-.6_r8)>1.e-14_r8)error stop 'empty leaf sink'
if(abs(sum(part)-sum(original))>1.e-14_r8)error stop 'substrate allocation lost'
call LeafNodeGrowthRanges(35.4_r8,9.8_r8,36,3,25,1,1,lf,ll,sf,sl)
if(ll>=lf.or.sl>=sf)error stop 'circular slot owner violated'
part=original
call LimitNodeGrowthAllocation(part,lf,ll,sf,sl)
if(any(part(1:3)/=0._r8).or.abs(part(4)-.8_r8)>1.e-14_r8)error stop 'empty stalk sink'
! Set up actual tissue, including occupied nodes above the cap.
plt_site%PlantPopuLive_pft=10._r8
plt_site%ZERO=1.e-12_r8
plt_allom%FracGroth2Node_pft=1._r8
plt_allom%rProteinC2LeafN_pft=1._r8
plt_allom%rProteinC2LeafP_pft=1._r8
plt_morph%NumCogrowthNode_pft=3
plt_morph%SLA1_pft=.02_r8
plt_morph%PetolShethLen2Mass_pft=.125_r8
plt_morph%NodeLenPergC_pft=.1_r8
plt_morph%SinePetolShethAngle_pft=.7_r8
plt_morph%SineBranchAngle_pft=.7_r8
plt_morph%HypocotHeight_pft=.06_r8
plt_morph%SeedDepth_pft=.05_r8
plt_morph%MainBranchNum_pft=1
plt_morph%rLen2WidthLeaf_pft=2._r8
plt_morph%StemSpecVolume_pft=.001_r8
plt_morph%LeafAngleClass_pft=1._r8/3._r8
plt_rad%SineLeafAngle=[.2_r8,.6_r8,.9_r8]
plt_pheno%iPlantPhenolPattern_pft=iplt_annual
plt_pheno%KHiestGroLeafNode_brch=13
plt_biom%StalkStrutElms_brch(:,1,1)=[13._r8,1.3_r8,.13_r8]
do k=1,13
  plt_biom%LeafElmntNode_brch(:,k,1,1)=[1._r8,.1_r8,.01_r8]
  plt_biom%PetolShethElmntNode_brch(:,k,1,1)=[.5_r8,.05_r8,.005_r8]
  plt_biom%StructInternodeElms_brch(:,k,1,1)=[1._r8,.1_r8,.01_r8]
  plt_biom%LeafProteinC_node(k,1,1)=.1_r8
  plt_biom%PetoleProteinC_node(k,1,1)=.05_r8
  plt_morph%LeafArea_node(k,1,1)=.02_r8
  plt_morph%PetoleLength_node(k,1,1)=.01_r8
  plt_morph%StalkNodeVertLength_brch(k,1,1)=.01_r8
  plt_morph%StalkNodeHeight_brch(k,1,1)=real(k,r8)*.01_r8
enddo
do k=0,5
  plt_morph%CanopyHeightZ_col(k)=real(k,r8)*.1_r8
enddo
plt_biom%LeafStrutElms_brch(:,1,1)=sum(plt_biom%LeafElmntNode_brch(:,:,1,1),dim=2)
call AllocateLeaf2CanopyLayers(1,1,1,1,[.5_r8])
layer_before=plt_biom%LeafLayerElms_node
area_layers_before=plt_morph%CanopyLeafArea_lnode
! Exercise the actual phenology assignment, stubbing only unrelated stage rates.
plt_morph%NumOfBranches_pft=1
plt_morph%NumOfLeaves_brch=12.4_r8
plt_pheno%LeafNumberAtFloralInit_brch=10.8_r8
plt_pheno%isPlantBranchAlive_brch=itrue
plt_pheno%iEmbryophyteType_pft=4
plt_ew%PSICanopy_pft=-.5_r8
call Emerged_plant_Phenology(1,1,1)
if(plt_pheno%KHiestGroLeafNode_brch(1,1)/=13)error stop 'retained frontier retreated'
if(plt_morph%NumOfLeaves_brch(1,1)/=12.4_r8)error stop 'developmental clock changed'
if(plt_pheno%doRemobilization_brch(1,1)/=ifalse)error stop 'cap triggered node turnover'
! Applying the cap by itself cannot remove any tissue or its canopy projection.
call LeafNodeGrowthRanges(12.4_r8,10.8_r8,13,3,25,1,1,lf,ll,sf,sl)
plt_biom%LeafLayerElms_node=0._r8
plt_biom%CanopyLeafCLyr_pft=0._r8
plt_morph%CanopyLeafArea_lnode=0._r8
plt_morph%CanopyLeafAreaZ_pft=0._r8
call AllocateLeaf2CanopyLayers(1,1,1,1,[.5_r8])
if(any(layer_before/=plt_biom%LeafLayerElms_node))error stop 'cap changed existing layer CNP'
if(any(area_layers_before/=plt_morph%CanopyLeafArea_lnode))error stop 'cap changed existing area'
do k=1,13
  if(any(abs(sum(plt_biom%LeafLayerElms_node(:,:,k,1,1),dim=2)- &
    plt_biom%LeafElmntNode_brch(:,k,1,1))>1.e-14_r8))error stop 'node/layer CNP closure'
enddo
! Grow the one eligible leaf; upper leaf and petiole tissue must stay unchanged.
upper=plt_biom%LeafElmntNode_brch(:,12:13,1,1)
before=sum(plt_biom%LeafElmntNode_brch(:,:,1,1),dim=2)
growth=[.6_r8,.06_r8,.006_r8]
call GrowLeavesOnBranch(1,1,1,1,lf,ll,growth,1._r8,1._r8,alloc)
after=sum(plt_biom%LeafElmntNode_brch(:,:,1,1),dim=2)
if(any(abs(after-before-growth)>1.e-14_r8))error stop 'leaf CNP growth lost'
if(any(upper/=plt_biom%LeafElmntNode_brch(:,12:13,1,1)))error stop 'upper leaves grew'
upper=plt_biom%PetolShethElmntNode_brch(:,12:13,1,1)
before=sum(plt_biom%PetolShethElmntNode_brch(:,:,1,1),dim=2)
call GrowPetolShethOnBranch(date,1,1,lf,ll,growth,1._r8,1._r8,alloc)
after=sum(plt_biom%PetolShethElmntNode_brch(:,:,1,1),dim=2)
if(any(abs(after-before-growth)>1.e-14_r8))error stop 'petiole CNP growth lost'
if(any(upper/=plt_biom%PetolShethElmntNode_brch(:,12:13,1,1)))error stop 'upper petioles grew'
stalk_tip=plt_morph%StalkNodeHeight_brch(13,1,1)
upper=plt_biom%StructInternodeElms_brch(:,12:13,1,1)
before=sum(plt_biom%StructInternodeElms_brch(:,:,1,1),dim=2)
call GrowStalkOnBranch(1,1,1,1,sf,sl,growth,1._r8)
after=sum(plt_biom%StructInternodeElms_brch(:,:,1,1),dim=2)
if(any(abs(after-before-growth)>1.e-14_r8))error stop 'stalk CNP growth lost'
if(any(plt_biom%StructInternodeElms_brch(:,12:13,1,1)/=upper))error stop 'upper internodes grew'
if(plt_morph%StalkNodeHeight_brch(13,1,1)<=stalk_tip)error stop 'upper geometry not propagated'
if(abs(plt_morph%StalkNodeHeight_brch(13,1,1)-plt_morph%StalkNodeHeight_brch(11,1,1)-.02_r8)>1.e-14_r8) &
  error stop 'upper internode lengths changed'
! A new eligible empty node must still receive its first biomass.
call LeafNodeGrowthRanges(13.2_r8,0._r8,14,3,25,1,1,lf,ll,sf,sl)
call GrowLeavesOnBranch(1,1,1,1,lf,ll,growth,1._r8,1._r8,alloc)
if(plt_biom%LeafElmntNode_brch(1,14,1,1)<=0._r8)error stop 'empty eligible node trapped'
print *, 'PASS: 4000 ordinary windows; cap boundaries, no reopening, allocation, circular slots, canopy CNP/area and organ growth'
end program
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fc', default=os.environ.get('FC', 'gfortran'))
    parser.add_argument('--branch-source', type=Path, default=ROOT / 'f90src/Plant_bgc/PlantBranchMod.F90')
    parser.add_argument('--phenology-source', type=Path, default=ROOT / 'f90src/Plant_bgc/PlantPhenolMod.F90')
    parser.add_argument('--math-source', type=Path, default=ROOT / 'f90src/Plant_bgc/PlantMathFuncMod.F90')
    args = parser.parse_args()
    source = fixture(args.branch_source.read_text(), args.math_source.read_text(), args.phenology_source.read_text()) + DRIVER
    with tempfile.TemporaryDirectory(prefix='ecosim_leaf_cap_') as work:
        folder = Path(work)
        (folder / 'check.F90').write_text(source)
        subprocess.run(shlex.split(args.fc) + ['-ffree-line-length-none', '-fcheck=all',
                       '-ffpe-trap=invalid,zero,overflow',
                       str(ROOT / 'f90src/Modelconfig/ElmIDMod.F90'), str(folder / 'check.F90'),
                       '-o', str(folder / 'check')], cwd=folder, check=True)
        subprocess.run([str(folder / 'check')], cwd=folder, check=True)


if __name__ == '__main__':
    main()
