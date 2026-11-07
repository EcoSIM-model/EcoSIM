"""Check secondary-root C/N/P and mycorrhizal totals after checkpoint recovery.

Seed conflicting per-axis pools in a scratch radiation-only SP case. Verify
that prescribed live fine-root C remains authoritative after restart, regardless
of axis count, and does not imply root N/P or mycorrhizal C.
"""

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess

import numpy as np
from netCDF4 import Dataset

from check_annual_restart import ROOT, stage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = ROOT / "examples/run_dir/radiation_test/Radiation_test.namelist"
    text = re.sub(r"(?m)![^\n]*$", "", source.read_text())

    def absolute_input(match):
        candidate = source.parent / match[1].rstrip()
        return "'" + str(candidate.resolve()) + "'" if candidate.is_file() else match[0]

    text = re.sub(r"'([^'\n]+)'", absolute_input, text)
    for key, value, group in (
        ("finidat", "''", "ecosim"), ("case_name", "'root-aggregate'", "ecosim"),
        ("continue_run", ".false.", "ecosim"), ("hist_nhtfrq", "1", "ecosim"),
        ("hist_fincl1", "'RootC_2nd_pvr','RootN_2nd_pvr','RootP_2nd_pvr',"
         "'MycoBiomC_pvr','MycoArbuBiomC_pft'", "ecosim"),
        ("stop_option", "'nsteps'", "ecosim_time"), ("stop_n", "1", "ecosim_time"),
        ("rest_opt", "'never'", "ecosim_time"),
    ):
        text = stage.set_value(text, key, value, group)

    def run(folder, current, steps):
        (folder / "runfile.namelist").write_text(current)
        with (folder / "run.log").open("w") as log:
            result = subprocess.run([str(executable), "runfile.namelist"], cwd=folder,
                                    stdout=log, stderr=subprocess.STDOUT)
        body = (folder / "run.log").read_text(errors="replace")
        completed = re.findall(r"wrote out restart data at nstep\s*=\s*(\d+)", body)
        assert result.returncode == 0 and "ENDRUN" not in body, folder
        assert completed and int(completed[-1]) == steps, folder

    seed = work / "seed"
    seed.mkdir()
    run(seed, text, 1)
    checkpoint, = seed.glob("*.ecosim.r.*.nc")
    report = {}
    for count in (2, 0):
        folder = work / f"axes-{count}"
        folder.mkdir()
        for path in seed.iterdir():
            if path.suffix == ".nc" or path.name == "rpointer.esim":
                shutil.copy2(path, folder / path.name)
        with Dataset(folder / checkpoint.name, "r+") as nc:
            nc["NumStructuralRootAxes_pft"][0] = count
            pool = nc["RootMyco2ndStrutElms_rpvr"]

            def index(**indices):
                return tuple(indices.get(dim, slice(None)) for dim in pool.dimensions)

            pool[index(pft=0)] = 0
            for layer, masses in enumerate(((2., 3.), (4., 7.))):
                for axis, mass in enumerate(masses):
                    for element, ratio in enumerate((1., .1, .01)):
                        pool[index(pft=0, rootaxs=axis, levsoi=layer,
                                   rootyps=0, elmnts=element)] = mass * ratio
                    # ElmIDMod: imycorr_arbu=2 (zero-based NetCDF index 1).
                    pool[index(pft=0, rootaxs=axis, levsoi=layer,
                               rootyps=1, elmnts=0)] = .2 * mass
            # SP must ignore axis storage, including nonzero root N/P/myco C.
            pool[index(pft=0, rootaxs=2)] = 100.
        current = stage.set_value(text, "continue_run", ".true.", "ecosim")
        current = stage.set_value(current, "stop_n", "2", "ecosim_time")
        run(folder, current, 2)
        final = folder / "root-aggregate.ecosim.r.1980-01-01-007200.nc"
        with Dataset(final) as nc:
            volumes = np.array(nc["VGeomLayer_vr"][1:3, 0])
            area = volumes[0] / nc["DLYR_3D"][1, 2, 0]
            depth = np.asarray(nc["CumDepz2LayBottom_vr"][1:3, 0])
            # The radiation fixture uses IROOTJK=4: Jackson live biomass
            # 0.50 kg dry mass m-2, carbon fraction 0.488 and beta=0.980.
            fractions = -np.diff(np.r_[1., .980 ** (100. * depth)])
            totals = .50 * .488 * 1.e3 * fractions
        history, = folder.glob("*.ecosim.h0.*.nc")
        with Dataset(history) as nc:
            for name, expected in (("RootC_2nd_pvr", totals * area / volumes),
                                   ("RootN_2nd_pvr", np.zeros(2)),
                                   ("RootP_2nd_pvr", np.zeros(2)),
                                   ("MycoBiomC_pvr", np.zeros(2))):
                values = nc[name][-1, :2, 0]
                assert not np.any(np.ma.getmaskarray(values)), (count, name)
                np.testing.assert_allclose(values, expected, rtol=1.e-6, atol=1.e-6)
            myco_column = nc["MycoArbuBiomC_pft"][-1, 0]
            assert not np.any(np.ma.getmaskarray(myco_column)), count
            np.testing.assert_array_equal(myco_column, 0.)
        report[str(count)] = "PASS"
    (work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
