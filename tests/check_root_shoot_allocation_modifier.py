"""Check the ECOSYS PTRT sampling point using production allocation code.

Compile CalcPartitionCoeff and its scalar math helpers with isolated state.
Verify the exchange modifier separately from all seven normalized organ fractions.
Requires a GNU-compatible Fortran compiler (FC); no model inputs are needed.
"""

import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def procedure(text, kind, name):
    """Extract a complete production procedure without rewriting its equations."""
    match = re.search(
        rf"^\s*(?:pure\s+)?{kind} {name}\b.*?^\s*end {kind} {name}\b[^\n]*",
        text, re.IGNORECASE | re.MULTILINE | re.DOTALL,
    )
    if match is None:
        raise ValueError(f"Cannot find {name}")
    return match[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--source", type=Path,
                        default=ROOT / "f90src/Plant_bgc/PlantBranchMod.F90")
    args = parser.parse_args()
    text = args.source.read_text()
    allocation = procedure(text, "subroutine", "CalcPartitionCoeff")
    organ_ids = re.findall(r"^\s*integer, parameter ::\s*ibrch_.*$", text, re.MULTILINE)
    math = (ROOT / "f90src/Minimath/MiniMathMod.F90").read_text()
    helpers = procedure(math, "function", "AZMAX1_s") + "\n" + procedure(math, "function", "isclose")
    pars = (ROOT / "f90src/Modelpars/PlantBGCPars.F90").read_text()
    parameters = []
    for name in ("PART2LEAF_MIN", "PART2PETOL_MIN", "XFRX", "FracHour4LeafoffRemob"):
        parameters.append(re.search(rf"^\s*{name}\s*=.*$", pars, re.MULTILINE)[0])
    fixture = """
module allocation_fixture
  ! Description:
  ! Supply isolated state for production partitioning and scalar math procedures.
  use iso_fortran_env, only: r8 => real64
  use ElmIDMod
  implicit none
  integer, parameter :: NumOfPlantMorphUnits=7
  real(r8), parameter :: tiny_val=1.e-12_r8
  real(r8) :: PART2LEAF_MIN, PART2PETOL_MIN, XFRX, FracHour4LeafoffRemob(0:5)
  type configuration
    integer :: NumOfPlantMorphUnits=7
  end type
  type energy_state
    real(r8) :: TdegCCanopy_pft(1)=10._r8, PSICanopy_pft(1)=-0.5_r8
  end type
  type biomass_state
    real(r8) :: SapwoodBiomassC_brch(2,1)=1._r8
    real(r8) :: StalkRsrvElms_brch(3,2,1)=0.5_r8
  end type
  type site_state
    real(r8) :: ZERO=1.e-12_r8, AREA3(1)=1._r8
    integer :: NU=1
  end type
  type phenology_state
    real(r8) :: HourReq4LeafOff_brch(2,1)=1._r8, Hours4LeafOff_brch(2,1)=0._r8
    real(r8) :: TotReproNodeNumNormByMatrgrp_brch(2,1)=0._r8
    real(r8) :: TotalNodeNumNormByMatgrp_brch(2,1)=0._r8
    real(r8) :: TCChill4Seed_pft(1)=15._r8, HoursDoingRemob_brch(2,1)=0._r8
    integer :: iPlantCalendar_brch(10,2,1)=0, iPlantRootProfile_pft(1)=1
    integer :: iPlantDevelopPattern_pft(1)=inondeterminate, iPlantTurnoverPattern_pft(1)=0
    integer :: iPlantPhenolPattern_pft(1)=iplt_annual, iPlantPhenolType_pft(1)=iphenotyp_evgreen
  end type
  type morphology_state
    real(r8) :: NodeLenPergC_pft(1)=0.1_r8, CanopyLeafArea_pft(1)=1._r8
    real(r8) :: PARTS_brch(7,2,1)=0._r8
  end type
  type(configuration) :: pltpar
  type(energy_state) :: plt_ew
  type(biomass_state) :: plt_biom
  type(site_state) :: plt_site
  type(phenology_state) :: plt_pheno
  type(morphology_state) :: plt_morph
  interface AZMAX1
    module procedure AZMAX1_s
  end interface
""" + "\n".join(organ_ids) + """
contains
  subroutine PrintInfo(message)
    ! Description:
    ! Suppress production log messages during the isolated allocation test.
    character(len=*), intent(in) :: message !unused log message
  end subroutine PrintInfo
  subroutine DebugPrint(message, value)
    ! Description:
    ! Suppress production debug messages during the isolated allocation test.
    character(len=*), intent(in) :: message !unused log message
    integer, intent(in) :: value !unused PFT index
  end subroutine DebugPrint
""" + helpers + "\n" + allocation + """
end module allocation_fixture

program check_root_shoot_allocation_modifier
  ! Description:
  ! Exercise the ECOSYS sampling order independently of normalized organ growth.
  use allocation_fixture
  implicit none
  real(r8) :: expected(7), normal(7)
""" + "\n".join(parameters) + """

  expected=[0.725_r8,0.275_r8,0._r8,0._r8,0._r8,0._r8,0._r8]
  call check('pre-floral',1._r8,expected)

  plt_pheno%iPlantCalendar_brch(ipltcal_InitFloral,:,1)=150
  plt_pheno%TotalNodeNumNormByMatgrp_brch=0.25_r8
  normal=[0.475_r8,0.175_r8,0.21_r8,0.105_r8,0.0175_r8,0.0175_r8,0._r8]
  call check('vegetative',0.65_r8,normal)

  plt_biom%StalkRsrvElms_brch(ielmc,:,1)=0._r8
  expected=0.9_r8*normal
  expected(ibrch_resrv)=0.1945_r8
  call check('low reserves',0.585_r8,expected)

  plt_biom%StalkRsrvElms_brch(ielmc,:,1)=2._r8
  expected=normal
  expected(ibrch_stalk)=0.315_r8
  expected(ibrch_resrv)=0._r8
  call check('excess reserves',0.65_r8,expected)
  plt_biom%StalkRsrvElms_brch=0.5_r8

  plt_morph%CanopyLeafArea_pft=7.5_r8
  expected=[0.2375_r8,0.0875_r8,0.535_r8,0.105_r8,0.0175_r8,0.0175_r8,0._r8]
  call check('LAI redirection',0.325_r8,expected)
  plt_morph%CanopyLeafArea_pft=1._r8

  plt_morph%NodeLenPergC_pft=0._r8
  expected=normal
  expected(ibrch_stalk)=0._r8
  expected=expected/0.79_r8
  call check('SNL1 zero',0.65_r8,expected)
  plt_morph%NodeLenPergC_pft=0.1_r8

  plt_pheno%iPlantCalendar_brch(ipltcal_Anthesis,:,1)=170
  plt_pheno%iPlantCalendar_brch(ipltcal_BeginSeedFill,:,1)=180
  expected=[0.05_r8,0.02_r8,0.11625_r8,0._r8,0.11625_r8,0.11625_r8,0.58125_r8]
  call check('grain fill',0.07_r8,expected)
  plt_pheno%iPlantCalendar_brch(ipltcal_EndSeedFill,:,1)=250
  expected=[0.05_r8,0.02_r8,0._r8,0._r8,0.11625_r8,0.11625_r8,0._r8]/0.3025_r8
  call check('post-maturity',0.07_r8,expected)

  plt_pheno%iPlantDevelopPattern_pft=ideterminate
  expected=0._r8
  call check('zero total allocation',0._r8,expected)

  plt_pheno%iPlantCalendar_brch=0
  plt_pheno%iPlantDevelopPattern_pft=inondeterminate
  plt_pheno%iPlantPhenolPattern_pft=iplt_perennial
  plt_pheno%iPlantPhenolType_pft=iphenotyp_coldecid
  plt_pheno%iPlantTurnoverPattern_pft=1
  plt_pheno%Hours4LeafOff_brch=1._r8
  expected=[0._r8,0._r8,0.5_r8,0.5_r8,0._r8,0._r8,0._r8]
  call check('deciduous remobilization',1._r8,expected)
  print *, 'PASS: 10 allocation scenarios on two branches; exchange modifier and organ fractions'
contains
  subroutine check(label, modifier_expected, organ_expected)
    ! Description:
    ! Check the exchange modifier, normalized fractions and stored diagnostics.
    character(len=*), intent(in) :: label !scenario name
    real(r8), intent(in) :: modifier_expected !ECOSYS PTRT at the original sampling point, [-]
    real(r8), intent(in) :: organ_expected(7) !normalized organ allocation, [-]
    real(r8) :: part(7), modifier
    integer :: nb, eligible, active
    do nb=1,2
      call CalcPartitionCoeff(250,12,nb,1,part,eligible,active,modifier)
      if(abs(modifier-modifier_expected)>1.e-12_r8)then
        print *, label, ': wrong exchange modifier: ',modifier,modifier_expected
        error stop 1
      endif
      if(any(abs(part-organ_expected)>1.e-12_r8))then
        print *, label, ': wrong organ fractions: ',part
        error stop 2
      endif
      if(any(abs(plt_morph%PARTS_brch(:,nb,1)-part)>1.e-12_r8))error stop 'stored fractions differ'
    enddo
  end subroutine check
end program check_root_shoot_allocation_modifier
"""
    with tempfile.TemporaryDirectory(prefix="ecosim_root_shoot_modifier_") as work:
        folder = Path(work)
        source = folder / "check.F90"
        source.write_text(fixture)
        executable = folder / "check"
        subprocess.run(shlex.split(args.fc) + ["-ffree-line-length-none", "-fcheck=all",
                       "-fdefault-real-8", "-fdefault-double-8",
                       "-ffpe-trap=invalid,zero,overflow",
                       str(ROOT / "f90src/Modelconfig/ElmIDMod.F90"), str(source),
                       "-o", str(executable)], cwd=folder, check=True)
        subprocess.run([str(executable)], cwd=folder, check=True)


if __name__ == "__main__":
    main()
