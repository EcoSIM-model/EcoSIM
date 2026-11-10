"""Check NH3 gas delivery and budgets using the production source/sink procedure.

Compile the extracted procedure and its actual scalar sink limiter in isolation.
This avoids a test-only public model interface and isolates fertilizer NH3 delivery
from atmospheric exchange. Requires a GNU-compatible Fortran compiler (FC).
"""

import argparse
import os
from pathlib import Path
import re
import shlex
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def procedure(path, name):
    """Read one complete production procedure without duplicating its equations."""
    match = re.search(
        rf"^\s*subroutine {name}\b.*?^\s*end subroutine {name}\b[^\n]*",
        path.read_text(), re.IGNORECASE | re.MULTILINE | re.DOTALL,
    )
    if match is None:
        raise ValueError(f"Cannot find {name} in {path}")
    return match[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fc", default=os.environ.get("FC", "gfortran"))
    parser.add_argument("--source", type=Path,
                        default=ROOT / "f90src/Transport/Nonsalt/TranspNoSaltFastMod.F90")
    args = parser.parse_args()
    apply = procedure(args.source, "ApplyFastSourceSink")
    limiter = procedure(ROOT / "f90src/Minimath/MiniMathMod.F90", "SubstrateDribbling_scal")
    fixture = """
module source_fixture
  ! Description:
  ! Provide minimal state for the unmodified production source/sink procedures.
  use iso_fortran_env, only: r8 => real64
  implicit none
  integer, parameter :: idg_beg=1, idg_NH3=3
  integer :: NU_col(1,1)=1, NL_col(1,1)=2
  real(r8) :: trcs_solml2_vr(3,0:2,1,1)=0._r8
  real(r8) :: trcg_gasml2_vr(3,0:2,1,1)=0._r8
  real(r8) :: trcs_solml_drib_vr(3,0:2,1,1)=0._r8
  real(r8) :: RBGCSrceGasMM_vr(3,0:2,1,1)=0._r8
  real(r8) :: RBGCSinkGasMM_vr(3,0:2,1,1)=0._r8
  real(r8) :: trcg_NetPro_fast_col(3,1,1)=0._r8
  real(r8) :: trcs_netProd_lit_fast_col(3,1,1)=0._r8
  real(r8) :: RGasNetProdSoil_col(3,1,1)=0._r8
  interface SubstrateDribbling
    module procedure SubstrateDribbling_scal
  end interface
contains
""" + apply + "\n" + limiter + """
end module source_fixture

program check_nh3_source
  ! Description:
  ! Check delivered mass, gas versus water routing, budgets and sink carryover.
  use source_fixture
  implicit none
  integer :: step, scenario
  real(r8) :: broadcast, banded, expected
  real(r8), parameter :: tol=1.e-12_r8

  do scenario=1,3
    broadcast=0._r8
    banded=0._r8
    if(scenario/=2)broadcast=1._r8
    if(scenario/=1)banded=3._r8
    expected=broadcast+banded
    trcs_solml2_vr=7._r8
    trcg_gasml2_vr=0._r8
    trcs_solml_drib_vr=0._r8
    RBGCSrceGasMM_vr=0._r8
    RBGCSinkGasMM_vr=0._r8
    trcg_NetPro_fast_col=0._r8
    trcs_netProd_lit_fast_col=0._r8
    RGasNetProdSoil_col=0._r8
    ! Stage a one-hour source over four gas substeps and two soil layers.
    RBGCSrceGasMM_vr(idg_NH3,1,1,1)=broadcast/4._r8
    RBGCSrceGasMM_vr(idg_NH3,2,1,1)=banded/4._r8
    RBGCSrceGasMM_vr(1,1,1,1)=0.1_r8
    do step=1,4
      call ApplyFastSourceSink(1,12,1,1,1,1,1,step)
    enddo
    if(abs(sum(trcg_gasml2_vr(idg_NH3,:,1,1))-expected)>tol) &
      error stop 'NH3 source mass did not reach the gas pool'
    if(abs(trcg_gasml2_vr(idg_NH3,1,1,1)-broadcast)>tol) error stop 'wrong source layer'
    if(abs(trcg_gasml2_vr(idg_NH3,2,1,1)-banded)>tol) error stop 'wrong band source layer'
    if(any(abs(trcs_solml2_vr(idg_NH3,:,1,1)-7._r8)>tol)) &
      error stop 'NH3 gas source was also added to water'
    if(abs(trcg_NetPro_fast_col(idg_NH3,1,1)-expected)>tol) error stop 'NH3 column budget'
    if(abs(RGasNetProdSoil_col(idg_NH3,1,1)-expected)>tol) error stop 'NH3 soil budget'
    if(abs(trcs_solml2_vr(1,1,1,1)-7.4_r8)>tol) error stop 'other gas source changed'
    if(abs(trcg_NetPro_fast_col(1,1,1)-0.4_r8)>tol) error stop 'other gas budget changed'
  enddo

  ! A source must be available before a sink; unmet demand remains in carryover.
  trcg_gasml2_vr=0._r8
  trcg_gasml2_vr(idg_NH3,1,1,1)=0.25_r8
  RBGCSrceGasMM_vr=0._r8
  RBGCSinkGasMM_vr=0._r8
  trcs_solml_drib_vr=0._r8
  trcg_NetPro_fast_col=0._r8
  RGasNetProdSoil_col=0._r8
  RBGCSrceGasMM_vr(idg_NH3,1,1,1)=0.5_r8
  RBGCSinkGasMM_vr(idg_NH3,1,1,1)=1._r8
  call ApplyFastSourceSink(1,12,1,1,1,1,1,1)
  if(abs(trcg_gasml2_vr(idg_NH3,1,1,1))>tol) error stop 'negative/incorrect gas pool'
  if(abs(trcs_solml_drib_vr(idg_NH3,1,1,1)-0.25_r8)>tol) error stop 'wrong unmet sink'
  if(abs(trcg_NetPro_fast_col(idg_NH3,1,1)+0.5_r8)>tol) error stop 'net source/sink budget'
  if(abs(RGasNetProdSoil_col(idg_NH3,1,1)+0.5_r8)>tol) error stop 'net soil budget'
  RBGCSinkGasMM_vr=0._r8
  call ApplyFastSourceSink(1,12,1,1,1,1,1,2)
  if(abs(trcg_gasml2_vr(idg_NH3,1,1,1)-0.25_r8)>tol) error stop 'sink carryover not met'
  if(abs(trcs_solml_drib_vr(idg_NH3,1,1,1))>tol) error stop 'sink carryover remains'
  if(abs(trcg_NetPro_fast_col(idg_NH3,1,1))>tol) error stop 'cumulative column budget'
  if(abs(RGasNetProdSoil_col(idg_NH3,1,1))>tol) error stop 'cumulative soil budget'
  RBGCSrceGasMM_vr=0._r8
  call ApplyFastSourceSink(1,12,1,1,1,1,1,3)
  if(abs(trcg_gasml2_vr(idg_NH3,1,1,1)-0.25_r8)>tol) error stop 'zero-source state changed'
  print *, 'PASS: NH3 source delivery, phase/layer routing, budgets and sink carryover'
end program check_nh3_source
"""
    with tempfile.TemporaryDirectory(prefix="ecosim_nh3_source_") as work:
        folder = Path(work)
        source = folder / "check.F90"
        source.write_text(fixture)
        executable = folder / "check"
        subprocess.run(shlex.split(args.fc) + ["-ffree-line-length-none", "-fcheck=all",
                       "-ffpe-trap=invalid,zero,overflow", str(source), "-o", str(executable)],
                       cwd=folder, check=True)
        subprocess.run([str(executable)], cwd=folder, check=True)


if __name__ == "__main__":
    main()
