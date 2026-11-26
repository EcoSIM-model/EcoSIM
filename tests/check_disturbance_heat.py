"""Check fire/removal heat withdrawal against current and legacy capacities."""

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
    parser.add_argument("--fc", default=os.environ.get("FC", "agf13"))
    args = parser.parse_args()
    source = (ROOT / "f90src/Disturbances/SoilDisturbMod.F90").read_text()
    helper = re.search(
        r"^  pure function RemovedOrganicHeatCapacity\b.*?^  end function RemovedOrganicHeatCapacity\b",
        source, re.M | re.S | re.I,
    )
    assert helper
    mixing_source = (ROOT / "f90src/Balances/TillageMixMod.F90").read_text()
    mixing = re.search(r"^  subroutine MixExtensive1D\b.*?^  end subroutine MixExtensive1D\b",
                       mixing_source, re.M | re.S | re.I)
    assert mixing
    fixture = """module fixture
  use iso_fortran_env, only: r8=>real64
  implicit none
  integer, parameter :: JZ=2
contains
""" + helper[0] + "\n" + mixing[0] + """
end module
program check_heat
  use fixture
  implicit none
  real(r8) :: carbon,organic,other,removed,retained,fraction,capacity,temperature
  real(r8) :: loss,exported,energy,new_capacity,new_energy,new_temperature,nominal
  real(r8) :: solid(0:2),mineral(1:2),ti(2),fi(2)
  integer :: inventory,step
  temperature=290.73845698471968_r8
  do inventory=1,4
    carbon=5108576._r8
    organic=2.496e-6_r8*carbon/(.55_r8*1.3_r8)
    if(inventory==2)organic=1.5_r8 ! legacy inventory below current-C nominal capacity
    if(inventory==3)organic=0._r8 ! no represented organic sensible heat
    other=11.3_r8 ! retained mineral, water/ice and root capacities
    if(inventory==4)other=0._r8 ! organic-only dry litter
    do step=1,3
      fraction=.75_r8
      removed=carbon*fraction
      retained=carbon-removed
      capacity=other+organic
      energy=capacity*temperature
      loss=RemovedOrganicHeatCapacity(removed,retained,organic)
      if(abs(loss-fraction*organic)>1.e-12_r8)error stop 'wrong removal fraction'
      exported=loss*temperature
      new_capacity=capacity-loss
      new_energy=energy-exported
      if(new_capacity<=0._r8.or.new_energy<=0._r8)error stop 'negative retained inventory'
      new_temperature=new_energy/new_capacity
      if(abs(new_temperature-temperature)>1.e-10_r8)error stop 'donor temperature changed'
      if(abs(energy-new_capacity*new_temperature-exported)>1.e-9_r8)error stop 'heat lost or created'
      if(abs((new_capacity-other)-(organic-loss))>1.e-12_r8)error stop 'nonorganic capacity removed'
      organic=organic-loss
      carbon=retained
    enddo
  enddo
  if(RemovedOrganicHeatCapacity(0._r8,0._r8,0._r8)/=0._r8)error stop 'empty pool'
  if(RemovedOrganicHeatCapacity(0._r8,1._r8,2._r8)/=0._r8)error stop 'zero removal'
  if(RemovedOrganicHeatCapacity(1._r8,0._r8,2._r8)/=2._r8)error stop 'complete removal'
  ! Recorded SOAP failure: nominal removed-material heat exceeded the entire store.
  nominal=2.496e-6_r8*3831517.1367590865_r8/(.55_r8*1.3_r8)*temperature
  if(nominal<=3738.3793748893236_r8)error stop 'fixture no longer reproduces old overdraw'
  ! Tillage includes surface organic material but no surface mineral inventory.
  solid=[8._r8,10._r8,20._r8]
  mineral=[7._r8,15._r8]
  ti=1._r8
  fi=[.4_r8,.6_r8]
  call MixExtensive1D(.3_r8,ti,fi,1,2,solid,.5_r8)
  call MixExtensive1D(.3_r8,ti,fi,1,2,mineral)
  if(abs(sum(solid)-38._r8)>1.e-12_r8)error stop 'solid capacity lost in tillage'
  if(abs(sum(mineral)-22._r8)>1.e-12_r8)error stop 'mineral capacity lost in tillage'
  if(any(solid(1:2)<mineral))error stop 'negative organic capacity after tillage'
  print *, 'PASS: current/legacy/zero/organic-only capacities; repeated removals conserve heat and temperature'
end program
"""
    with tempfile.TemporaryDirectory(prefix="ecosim_disturbance_heat_") as work:
        folder = Path(work)
        (folder / "check.F90").write_text(fixture)
        executable = folder / "check"
        subprocess.run(shlex.split(args.fc) + [
            "-ffree-line-length-none", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
            str(folder / "check.F90"), "-o", str(executable),
        ], cwd=folder, check=True)
        subprocess.run([str(executable)], cwd=folder, check=True)


if __name__ == "__main__":
    main()
