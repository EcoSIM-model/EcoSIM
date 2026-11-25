"""Check production fire-schedule lookups across annual restart initialization."""

import argparse
import os
from pathlib import Path
import shlex
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fc", default=os.environ.get("FC", "agf13"))
    args = parser.parse_args()
    stubs = """
module data_kind_mod
  use iso_fortran_env, only: DAT_KIND_R8 => real64
end module
module EcoSIMCtrlMod
  character(len=80) :: soil_mgmt_in='fixture.nc'
end module
module abortutils
contains
  subroutine endrun(message,line)
    character(len=*), intent(in) :: message
    integer, optional, intent(in) :: line
    print *,message
    error stop 1
  end subroutine
end module
"""
    driver = """
program check_schedule
  use FireMod
  implicit none
  character(len=21) :: entry
  integer :: pass
  call config_fire('2021/a_pmf,b_smf;2023/c_pmf,d_smf')
  if(check_fire(2020,entry))error stop 'premature event'
  do pass=1,2
    if(.not.check_fire(2021,entry))error stop 'January reread lost fire'
    if(entry/='a_pmf,b_smf')error stop 'wrong 2021 event'
  enddo
  if(check_fire(2022,entry))error stop 'repeated event in a later year'
  do pass=1,2
    if(.not.check_fire(2023,entry))error stop 'lost second event'
    if(entry/='c_pmf,d_smf')error stop 'wrong 2023 event'
  enddo
  if(check_fire(2024,entry))error stop 'event after end of schedule'
  call config_fire('2021/a_pmf,b_smf;2023/c_pmf,d_smf')
  if(.not.check_fire(2023,entry))error stop 'failed to skip past event'
  if(entry/='c_pmf,d_smf')error stop 'wrong event after restart skip'
  call config_fire('')
  if(check_fire(2021,entry))error stop 'event in empty schedule'
  print *, 'PASS: fire schedule survives repeated annual reads and restart-year skips'
end program
"""
    with tempfile.TemporaryDirectory(prefix="ecosim_fire_schedule_") as work:
        folder = Path(work)
        (folder / "stubs.F90").write_text(stubs)
        (folder / "check.F90").write_text(driver)
        executable = folder / "check"
        subprocess.run(shlex.split(args.fc) + [
            "-ffree-line-length-none", "-fcheck=all", str(folder / "stubs.F90"),
            str(ROOT / "f90src/Disturbances/FireMod.F90"), str(folder / "check.F90"),
            "-o", str(executable),
        ], cwd=folder, check=True)
        subprocess.run([str(executable)], cwd=folder, check=True)


if __name__ == "__main__":
    main()
