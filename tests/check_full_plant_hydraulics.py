"""Check same-hour soil root transfers from an active single-column checkpoint.

Run a two-day scratch branch with hourly history and checkpoints. Every hour's
soil-applied water must equal the current accepted total (roots plus mycorrhizae),
without reapplying the checkpoint's retained completed-hour transfer.
Check that accepted donor heat is retained for history/coupling, with the model's
water and energy checks enabled. --substeps exercises a different soil partition.
"""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

from check_annual_restart import stage
from check_restart_exactness import compare

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--seed", type=Path, required=True)
    parser.add_argument("--namelist", type=Path,
                        default=ROOT / "examples/run_dir/dryland/dryland.namelist")
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--substeps", type=int,
                        help="override each NPXS entry with this soil substep count")
    args = parser.parse_args()
    if args.substeps is not None and args.substeps < 1:
        parser.error("--substeps must be positive")
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = args.namelist.resolve()
    text = re.sub(r"(?m)![^\n]*$", "", source.read_text())

    def absolute_input(match):
        candidate = source.parent / match[1].rstrip()
        return "'" + str(candidate.resolve()) + "'" if candidate.is_file() else match[0]

    text = re.sub(r"'([^'\n]+)'", absolute_input, text)
    for key, value in (("hist_nhtfrq", "1"), ("hist_mfilt", "1000"),
                       ("hist_fincl1", "'RootH2OUptake_vr','RootH2OUptk_pvr','HeatUptk_vr'")):
        text = stage.set_value(text, key, value, "ecosim")
    if args.substeps is not None:
        text = stage.set_value(text, "NPXS", ",".join([str(args.substeps)] * 5), "ecosim")
    parent = work / "parent.namelist"
    parent.write_text(text)
    run = work / "branch"
    case_name = args.seed.name.split(".ecosim.r.", 1)[0] + ".hydrauliccheck"
    subprocess.run([sys.executable, str(ROOT / ".agents/skills/ecosim-restart-run/scripts/stage_restart_run.py"),
                    "--namelist", str(parent), "--restart", str(args.seed.resolve()), "--out", str(run),
                    "--case-name", case_name, "--stop-n", "2", "--stop-option", "ndays"],
                   stdout=(work / "stage.log").open("w"), check=True)
    namelist = run / "runfile.namelist"
    text = stage.set_value(namelist.read_text(), "rest_opt", "'nsteps'", "ecosim_time")
    text = stage.set_value(text, "rest_frq", "1", "ecosim_time")
    namelist.write_text(text)
    manifest = json.loads((run / "restart_run.json").read_text())
    with (run / "run.log").open("w") as log:
        subprocess.run([str(args.executable.resolve()), namelist.name], cwd=run,
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    body = (run / "run.log").read_text(errors="replace")
    steps = re.findall(r"wrote out restart data at nstep\s*=\s*(\d+)", body)
    assert "ENDRUN" not in body and steps and int(steps[-1]) == manifest["expected_nstep"], run
    with Dataset(manifest["inputs"]["grid_file_in"]) as grid:
        area = np.asarray(grid["DHI"][:]) * np.asarray(grid["DVI"][:])
        assert area.size == 1, "this focused test requires a single column"
        area = area.item()
    with Dataset(args.seed) as seed:
        assert seed["TWaterPlantRoot2Soil_vr"].shape[1] == 1
        retained = np.ma.filled(seed["TWaterPlantRoot2Soil_vr"][:, 0], 0.) * 1000. / area
    checkpoints = sorted(run.glob(case_name + ".ecosim.r.*.nc"))
    assert len(checkpoints) == 48, (len(checkpoints), run)
    expected = []
    expected_heat = []
    for path in checkpoints:
        assert not compare(path, path), path  # Valid-value finite checks.
        with Dataset(path) as checkpoint:
            current = np.ma.filled(checkpoint["TWaterPlantRoot2Soil_vr"][:, 0], 0.) * 1000. / area
            expected.append(current)
            # The heat profile includes layer zero (surface litter); root soil
            # history starts at layer one. Neither accepted array is subcycled
            # in place, because the external ATS soil solver also reads them.
            expected_heat.append(np.ma.filled(checkpoint["THeatLossRoot2Soil_vr"][1:, 0], 0.) / area)
    expected = np.array(expected)
    actual = []
    actual_heat = []
    for path in sorted(run.glob("*.ecosim.h0.*.nc")):
        with Dataset(path) as hist:
            actual.append(hist["RootH2OUptake_vr"][:, :, 0])
            actual_heat.append(hist["HeatUptk_vr"][:, :, 0])
    actual = np.ma.concatenate(actual)
    assert actual.shape == expected.shape
    valid = ~np.ma.getmaskarray(actual)
    assert np.isfinite(actual.data[valid]).all()
    assert np.max(np.abs(expected[1:])) > 1.e-6, "seed must exercise active root water exchange"
    np.testing.assert_allclose(actual.data[valid], expected[valid], rtol=2.e-6, atol=1.e-10)
    actual_heat = np.ma.concatenate(actual_heat)
    expected_heat = np.asarray(expected_heat)
    valid_heat = ~np.ma.getmaskarray(actual_heat)
    assert actual_heat.shape == expected_heat.shape
    assert np.isfinite(actual_heat.data[valid_heat]).all()
    np.testing.assert_allclose(actual_heat.data[valid_heat], expected_heat[valid_heat], rtol=2.e-6, atol=1.e-10)
    report = {"hours_checked": 48, "completed_hour_transfer_not_reapplied": True,
              "max_seed_retained_exchange_mm": float(np.max(np.abs(retained))),
              "max_exchange_mm": float(np.max(np.abs(expected[1:]))),
              "same_hour_total_root_and_myco_transfer": True,
              "accepted_hourly_donor_heat_preserved": True,
              "substeps_override": args.substeps}
    (work / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print("PASS: 48 hourly soil-applied transfers match current roots plus mycorrhizae;", work)


if __name__ == "__main__":
    main()
