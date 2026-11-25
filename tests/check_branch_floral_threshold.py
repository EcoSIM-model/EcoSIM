"""Compile the production branch-stage predicate and check ECOSYS floral thresholds.

Exercise main stems, annual side branches (including a zero maturity requirement),
perennials and nonzero initial node counts. No full-model inputs are required.
"""

import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--source", type=Path,
                        default=ROOT / "f90src/Plant_bgc/PlantPhenolMod.F90")
    args = parser.parse_args()
    text = args.source.read_text()
    match = re.search(
        r"^  function CheckBranchNodeState\b.*?^  end function CheckBranchNodeState\b",
        text, re.IGNORECASE | re.MULTILINE | re.DOTALL,
    )
    if match is None:
        raise ValueError(f"Cannot find CheckBranchNodeState in {args.source}")
    durations = re.findall(
        r"^.*parameter\s*::\s*GrothStageNorm4(?:Vegeta|Reprod)Pheno\s*=.*$",
        text, re.IGNORECASE | re.MULTILINE,
    )
    fixture = """
module branch_fixture
  ! Description:
  ! Supply minimal state for the unmodified production branch-stage predicate.
  use iso_fortran_env, only: r8 => real64
  use ElmIDMod
  implicit none
  type phenology_state
    real(r8) :: MatureGroup_pft(2)=12._r8
    real(r8) :: MatureGroup_brch(3,2)=12._r8
    real(r8) :: ReprodNodeNumNormByMatrgrp_brch(3,2)=0._r8
    real(r8) :: NodeNumNormByMatgrp_brch(3,2)=0._r8
  end type
  type morphology_state
    real(r8) :: ShootNodeNum_brch(3,2)=0._r8
    real(r8) :: NumOfLeaves_brch(3,2)=0._r8
    real(r8) :: ShootNodeNumAtInitFloral_brch(3,2)=0._r8
  end type
  type(phenology_state) :: plt_pheno
  type(morphology_state) :: plt_morph
""" + "\n".join(durations) + """
contains
  subroutine PrintInfo(message)
    ! Description:
    ! Suppress logging in this isolated predicate test.
    character(len=*), intent(in) :: message !unused production log message
  end subroutine PrintInfo
""" + match[0] + """
end module branch_fixture

program check_branch_floral_threshold
  ! Description:
  ! Check thresholds inherited from ECOSYS hfunc.f, including strict boundaries.
  use branch_fixture
  implicit none
  integer :: nb, nz, step
  ! Annual main stem requires 12 nodes; side branches require 9 and 0.
  ! Perennial branches retain 12, with a four-node initial offset in this fixture.
  integer, parameter :: threshold(3,2)=reshape([12,9,0,16,16,16],[3,2])
  logical :: actual, expected

  plt_pheno%MatureGroup_brch(:,1)=[12._r8,9._r8,0._r8]
  plt_morph%ShootNodeNumAtInitFloral_brch(:,2)=4._r8
  do nz=1,2
    do nb=1,3
      do step=0,36
        plt_morph%ShootNodeNum_brch(nb,nz)=real(step,r8)/2._r8
        expected=step>2*threshold(nb,nz)
        actual=CheckBranchNodeState(ipltcal_InitFloral,nb,nz,150,12)
        if(actual.neqv.expected)then
          print *, 'Wrong floral threshold: PFT, branch, nodes = ', nz,nb,real(step,r8)/2._r8
          error stop 1
        endif
      enddo
    enddo
  enddo
  print *, 'PASS: 222 main-stem/side-branch/perennial/offset floral-threshold checks'
end program check_branch_floral_threshold
"""
    with tempfile.TemporaryDirectory(prefix="ecosim_branch_floral_") as work:
        folder = Path(work)
        source = folder / "check.F90"
        source.write_text(fixture)
        executable = folder / "check"
        subprocess.run(shlex.split(args.fc) + ["-ffree-line-length-none", "-fcheck=all",
                       "-ffpe-trap=invalid,zero,overflow",
                       str(ROOT / "f90src/Modelconfig/ElmIDMod.F90"), str(source),
                       "-o", str(executable)], cwd=folder, check=True)
        subprocess.run([str(executable)], cwd=folder, check=True)


if __name__ == "__main__":
    main()
