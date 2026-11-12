"""Check SP root water/heat queues and exact continuation across the first hour.

Runs the SatePhenol example for two days, stopping/resuming after hour one.
Existing example outputs and inputs are read-only; --workdir must not exist.
"""

import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    work = args.workdir.resolve()
    subprocess.run([
        sys.executable, str(ROOT / "tests/check_restart_exactness.py"),
        "--executable", str(args.executable.resolve()),
        "--namelist", str(ROOT / "examples/run_dir/SatePhenol/SatePhenol_test.namelist"),
        "--workdir", str(work), "--days", "2", "--split-hours", "1",
    ], check=True)
    checkpoint, = (work / "first_segment").glob("*.ecosim.r.*.nc")
    with Dataset(checkpoint) as nc:
        nc.set_auto_maskandscale(False)
        water = nc["TWaterPlantRoot2Soil_vr"][:, 0]
        heat = nc["THeatLossRoot2Soil_vr"][1:len(water) + 1, 0]
        soil_temperature = nc["TKS"][1:len(water) + 1, 0]
        canopy_temperature = float(nc["TKCanopy_pft"][0])
        # This fixture has one living PFT. The queued flux is positive toward soil.
        assert np.any(water != 0.), "SP accepted root water never reached its queue"
        assert np.all(np.isfinite(water)) and np.max(np.abs(water)) < 1.
        expected_heat = water * 4.186 * np.where(water > 0., canopy_temperature, soil_temperature)
        np.testing.assert_allclose(heat, expected_heat, rtol=2.e-13, atol=1.e-14)
    result = {
        "first_hour_root_water_m3_per_cell": float(np.sum(water)),
        "first_hour_root_heat_MJ_per_cell": float(np.sum(heat)),
        "donor_enthalpy_matches": True,
        "exact_restart_hours": 48,
        "restart_split_hours": 1,
    }
    (work / "root_queue_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
