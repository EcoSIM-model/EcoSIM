"""Check prescribed plant exchange under snow burial and exact restart recovery.

Runs two-day fixtures in a new scratch directory; source inputs remain read-only.
Requires NumPy and netCDF4. Includes a short buried PFT beneath an exposed one.
"""

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset

from check_annual_restart import ROOT, stage
from check_restart_exactness import compare


def check_case(folder, heights, areas):
    """Require buried PFT fluxes to vanish while exposed controls exchange."""
    continuous = folder / "continuous"
    checkpoints = list(continuous.glob("*.ecosim.r.*.nc"))
    assert len(checkpoints) == 1, checkpoints
    with Dataset(checkpoints[0]) as nc:
        snow_depth = float(nc["SnowDepth"][:].item())
        assert snow_depth > 0.01, "fixture lost its snow cover"
        exposed_heights = [h for h, a in zip(heights, areas) if h > 0.01 and a > 0.]
        if exposed_heights:
            assert snow_depth < min(exposed_heights), 'fixture buried the exposed control'
        total_area = nc['CanopyLeafArea_pft'][:len(heights)] + nc['CanopyStemSurfArea_pft'][:len(heights)]
        np.testing.assert_allclose(total_area, 1.5*np.asarray(areas))
        # This small one-m2 fixture must not accumulate fill-valued shoot heat.
        assert abs(float(nc['CanopyHeatStor_col'][:].item())) < 1.e3
        active_heights = nc['CanopyHeightLive_pft'][:len(heights)][np.asarray(areas) > 0.]
        np.testing.assert_allclose(nc['CanopyHeight_col'][:].item(), np.max(active_heights))
        for i, (height, area) in enumerate(zip(heights, areas)):
            if area == 0.:
                continue
            assert nc['PlantPopuLive_pft'][i] > 0.
            if height == 0.01:
                assert nc['LeafStalkAreaAct_pft'][i] == 0., 'buried PFT retained exposed area'
                for field in ("HeatXAir2PCan_pft", "Transpiration_pft"):
                    assert nc[field][i] == 0., (folder.name, field, i, nc[field][i])
                assert nc['WatHeldOnCanopy_pft'][i] == 0., 'buried PFT intercepted new rain'
                assert nc['SnowOnCanopy_pft'][i] == 0., 'buried PFT intercepted new snow'
    history, = continuous.glob("*.ecosim.h0.*.nc")
    assert not compare(history, history), history
    maxima = {}
    with Dataset(history) as nc:
        for field in ("CAN_H_pft", "CAN_LE_pft", "QVegTransp_pft"):
            values = np.ma.asarray(nc[field][:])
            assert values.shape[0] == 48, (field, values.shape)
            maxima[field] = []
            for i, (height, area) in enumerate(zip(heights, areas)):
                if area == 0.:
                    maxima[field].append(None)
                    continue
                assert not np.any(np.ma.getmaskarray(values[:, i])), field
                peak = float(np.max(np.abs(values[:, i])))
                maxima[field].append(peak)
                if height == 0.01:
                    assert peak == 0., (folder.name, field, i, peak)
                elif field == "CAN_H_pft":
                    assert peak > 1.e-6, "exposed control did not exchange sensible heat"
    if folder.name == 'shrinking':
        first, = (folder.parent/'shrinking'/'first_segment').glob('*.ecosim.r.*.nc')
        with Dataset(first) as nc:
            assert float(nc['CanopyHeight_col'][:].item()) > float(np.max(active_heights)), 'column height did not decrease'
    return {"heights_m": heights, "area_multipliers": areas,
            "final_snow_depth_m": snow_depth, "maximum_absolute_flux": maxima}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = ROOT / "examples/run_dir/SatePhenol/SatePhenol_test.namelist"
    template = re.sub(r"(?m)![^\n]*$", "", source.read_text())

    def absolute_input(match):
        candidate = source.parent / match[1].rstrip()
        return repr(str(candidate.resolve())) if candidate.is_file() else match[0]

    template = re.sub(r"'([^'\n]+)'", absolute_input, template)
    report = {}
    for name, heights, areas in (("buried", [0.01], [1.]), ("mixed", [17., 0.01], [1., 1.]),
                                ("exposed", [17.], [1.]), ("empty_tall", [17., 0.01], [0., 1.]),
                                ("shrinking", [17.], [1.])):
        inputs = work / (name + "-inputs")
        inputs.mkdir()
        phenology = inputs / "phenology.nc"
        shutil.copyfile(ROOT / "examples/inputs/blodgett/Blodget_sp_phenol_20261103.nc", phenology)
        with Dataset(phenology, "r+") as nc:
            nc["NZ"][:] = len(heights)
            first_type = nc["pft_type"][0, 0, 0, :].copy()
            for field in ("LAI", "SAI", "HTOP", "PlantPopuLive_pft"):
                nc[field][:] = 0.
            for i, (height, area) in enumerate(zip(heights, areas)):
                nc["pft_type"][0, 0, i, :] = first_type
                nc["LAI"][0, 0, i, :] = 1.2*area
                nc["SAI"][0, 0, i, :] = .3*area
                nc["HTOP"][0, 0, i, :] = height
                # The SP input schema requires positive population in named
                # slots even when their current prescribed canopy area is zero.
                nc["PlantPopuLive_pft"][0, 0, i] = .6
            if name == 'shrinking':
                nc['HTOP'][0, 0, 0, :] = 10.
                nc['HTOP'][0, 0, 0, 11] = 17.
        with Dataset(phenology) as nc:
            assert nc["NZ"][:].item() == len(heights)
            if name != 'shrinking':
                np.testing.assert_allclose(nc["HTOP"][0, 0, :len(heights), 0], heights)
        grid = inputs / "grid.nc"
        shutil.copyfile(ROOT / "examples/inputs/blodgett/Blodget_grid_20251115.nc", grid)
        with Dataset(grid, "r+") as nc:
            nc["DPTHSX"][:] = 3.
        with Dataset(grid) as nc:
            assert nc["DPTHSX"][:].item() == 3.
        text = template
        for key, value in (
            ("case_name", repr(name)), ("grid_file_in", repr(str(grid))),
            ("sp_phenol_in", repr(str(phenology))), ("pft_mgmt_in", "'NO'"),
            ("forc_periods", "2013,2014,1"), ("continue_run", ".false."),
            ("finidat", "''"), ("hist_nhtfrq", "1"), ("hist_mfilt", "8760"),
            ("hist_fincl1", "'CAN_H_pft','CAN_LE_pft','QVegTransp_pft','TKQ_col'"),
        ):
            text = stage.set_value(text, key, value, "ecosim")
        namelist = inputs / "case.namelist"
        namelist.write_text(text)
        with (inputs / "comparison.log").open("w") as log:
            subprocess.run([sys.executable, str(ROOT / "tests/check_restart_exactness.py"),
                            "--executable", str(executable), "--namelist", str(namelist),
                            "--workdir", str(work / name), "--days", "2", "--split-hours", "24"],
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        report[name] = check_case(work / name, heights, areas)
        print("PASS:", name, json.dumps(report[name]), flush=True)
    (work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
