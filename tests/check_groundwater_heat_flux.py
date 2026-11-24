"""Check hourly recharge enthalpy and exact restart with the RiceUSTWT example.

Use a new --workdir. Source inputs and existing example outputs stay read-only.
The fixture reaches its first recharge at hour 12; the old heat accumulator
incorrectly repeated that input at hour 13. Hourly water/heat history must imply
a physically plausible donor temperature, independently checking heat units.
"""

import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

from check_restart_exactness import ROOT, setting


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    original = ROOT / "examples/run_dir/RiceUSTWT/riceUSTWT.namelist"
    text = original.read_text()

    def absolute_input(match):
        path = original.parent / match[1].rstrip()
        return "'" + str(path.resolve()) + "'" if path.is_file() else match[0]

    text = re.sub(r"'([^'\n]+)'", absolute_input, text)
    text = text.replace("hist_fincl1=", "hist_fincl1='HeatDischarge_FLX_col',", 1)
    text = setting(text, "hist_nhtfrq", "-1")
    fixture = work / "fixture.namelist"
    fixture.write_text(text)
    comparison = work / "comparison"
    subprocess.run([
        sys.executable, str(ROOT / "tests/check_restart_exactness.py"),
        "--executable", str(args.executable.resolve()), "--namelist", str(fixture),
        "--workdir", str(comparison), "--days", "2", "--split-hours", "12",
    ], check=True)

    history, = (comparison / "continuous").glob("*.ecosim.h0.*.nc")
    with Dataset(history) as nc:
        heat = np.ma.asarray(nc["HeatDischarge_FLX_col"][:, 0])
        water = np.ma.asarray(nc["QDischarge_FLX_col"][:, 0])
        assert nc["HeatDischarge_FLX_col"].units == "MJ/m2/hr"
        assert heat.shape == water.shape == (48,)
        assert not np.any(np.ma.getmaskarray(heat) | np.ma.getmaskarray(water))
        assert np.all(np.isfinite(heat)) and np.all(np.isfinite(water))
        wet = np.abs(water) > 1.e-12
        assert np.flatnonzero(wet)[0] == 11, "fixture no longer starts recharge at hour 12"
        assert np.all(heat[~wet] == 0.), "heat carried over into a zero-recharge hour"
        # Q is mm/hour and heat is MJ/m2/hour; cpw is 4.186 MJ/m3/K.
        donor_temperature = np.asarray(heat[wet] / water[wet]) * 1000. / 4.186
        assert np.all((donor_temperature > 240.) & (donor_temperature < 330.)), (
            "hourly heat includes an old flux or incorrect area/water conversion",
            donor_temperature,
        )
    result = {
        "hours": 48, "restart_split_hours": 12, "exact_restart": True,
        "first_recharge_hour": 12, "nonzero_recharge_hours": int(np.count_nonzero(wet)),
        "effective_donor_temperature_K": [float(donor_temperature.min()), float(donor_temperature.max())],
        "heat_units": "MJ/m2/hr",
    }
    (work / "groundwater_heat_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
