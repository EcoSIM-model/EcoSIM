"""Check prescribed canopy totals and exchange densities in 1 and 6 m2 cells.

Runs two-day simulations in a new scratch directory. Requires numpy/netCDF4.
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
from check_restart_exactness import compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = ROOT / "examples/run_dir/SatePhenol/SatePhenol_test.namelist"
    template = re.sub(r"(?m)![^\n]*$", "", source.read_text())

    def absolute_input(match):
        candidate = source.parent / match[1].rstrip()
        return "'" + str(candidate.resolve()) + "'" if candidate.is_file() else match[0]

    template = re.sub(r"'([^'\n]+)'", absolute_input, template)
    phenology = work / "phenology.nc"
    shutil.copyfile(ROOT / "examples/inputs/blodgett/Blodget_sp_phenol_20261103.nc", phenology)
    with Dataset(phenology, "r+") as nc:
        for field, value in (("LAI", 1.2), ("SAI", .3), ("HTOP", 17.)):
            nc[field][:] = 0.
            nc[field][0, 0, 0, :] = value
    with Dataset(phenology) as nc:
        leaf_index = float(nc["LAI"][0, 0, 0, 0])
        stem_index = float(nc["SAI"][0, 0, 0, 0])
    results = {}
    for area, widths in ((1, (1., 1.)), (6, (2., 3.))):
        folder = work / str(area)
        folder.mkdir()
        grid = folder / "grid.nc"
        shutil.copyfile(ROOT / "examples/inputs/blodgett/Blodget_grid_20251115.nc", grid)
        with Dataset(grid, "r+") as nc:
            nc["DHI"][:] = widths[0]
            nc["DVI"][:] = widths[1]
        with Dataset(grid) as nc:
            assert nc["DHI"][:].item() * nc["DVI"][:].item() == area
        text = template
        for key, value, group in (
            ("case_name", "'canopy-area'", "ecosim"),
            ("grid_file_in", "'" + str(grid) + "'", "ecosim"),
            ("sp_phenol_in", "'" + str(phenology) + "'", "ecosim"),
            ("pft_mgmt_in", "'NO'", "ecosim"),
            ("forc_periods", "2013,2014,1", "ecosim"),
            ("continue_run", ".false.", "ecosim"), ("finidat", "''", "ecosim"),
            ("hist_nhtfrq", "1", "ecosim"), ("hist_mfilt", "8760", "ecosim"),
            ("hist_fincl1", "'BLYR_RS_H2O_pft','STOM_RSC_H2O_pft','CAN_RN_pft','CAN_H_pft','CAN_LE_pft'", "ecosim"),
            ("stop_option", "'nsteps'", "ecosim_time"), ("stop_n", "48", "ecosim_time"),
            ("rest_opt", "'never'", "ecosim_time"),
        ):
            text = stage.set_value(text, key, value, group)
        (folder / "case.namelist").write_text(text)
        with (folder / "run.log").open("w") as log:
            subprocess.run([str(executable), "case.namelist"], cwd=folder,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        body = (folder / "run.log").read_text(errors="replace")
        steps = re.findall(r"wrote out restart data at nstep\s*=\s*(\d+)", body)
        assert "ENDRUN" not in body and steps and int(steps[-1]) == 48, folder
        checkpoint, = folder.glob("*.ecosim.r.*.nc")
        history, = folder.glob("*.ecosim.h0.*.nc")
        assert not compare(checkpoint, checkpoint)
        assert not compare(history, history)
        results[area] = (checkpoint, history)
        with Dataset(checkpoint) as nc:
            np.testing.assert_allclose(nc["CanopyLeafArea_pft"][0], leaf_index * area, rtol=1.e-12)
            np.testing.assert_allclose(nc["CanopyStemSurfArea_pft"][0], stem_index * area, rtol=1.e-12)
            np.testing.assert_allclose(nc["CanopyLeafAreaZ_pft"][:, 0].sum(), leaf_index * area, rtol=1.e-12)
            np.testing.assert_allclose(nc["CanopyStemSurfAreaZ_pft"][:, 0].sum(), stem_index * area, rtol=1.e-12)
        with Dataset(history) as nc:
            np.testing.assert_allclose(nc["ECO_LAI_col"][:, 0], leaf_index, rtol=1.e-6)
            np.testing.assert_allclose(nc["ECO_SAI_col"][:, 0], stem_index, rtol=1.e-6)
    totals = ("CanopyLeafArea_pft", "CanopyStemSurfArea_pft", "LeafStalkAreaAct_pft",
              "CanopyLeafArea_col", "StemArea_col", "LeafStalkAreaAll_col",
              "CanopyLeafAreaZ_pft", "CanopyStemSurfAreaZ_pft",
              "CanopyLeafAareZ_col", "CanopyStemAareZ_col",
              "CanopyLeafSheathC_pft", "CanopySapwoodC_pft")
    geometry = ("RoughnessLength", "ZeroPlaneDisplacem_col", "RawIsoTAtm2CanopySinkZ_col")
    densities = ("ECO_LAI_col", "ECO_SAI_col", "BLYR_RS_H2O_pft", "STOM_RSC_H2O_pft")
    flux_fields = ("ECO_RADSW_col", "ECO_HeatS_col", "ECO_LE_col", "CAN_RN_pft", "CAN_H_pft", "CAN_LE_pft")
    flux_differences = {}
    with Dataset(results[1][0]) as a, Dataset(results[6][0]) as b:
        for name in totals:
            np.testing.assert_allclose(b[name][:], a[name][:] * 6., rtol=1.e-10, atol=1.e-12,
                                       err_msg=name)
        for name in geometry:
            np.testing.assert_allclose(b[name][:], a[name][:], rtol=1.e-10, atol=1.e-12,
                                       err_msg=name)
    with Dataset(results[1][1]) as a, Dataset(results[6][1]) as b:
        for name in densities:
            np.testing.assert_allclose(b[name][:], a[name][:], rtol=1.e-6, atol=1.e-6,
                                       err_msg=name)
        # Geometry/resistance invariance is the regression contract here.
        # Record coupled energy-flux differences separately: soil water/energy
        # and canopy-air feedback can expose other dependencies on cell size.
        # No tolerance is widened to declare those additional checks a pass.
        for name in flux_fields:
            x, y = a[name][:], b[name][:]
            peak = float(np.ma.abs(x).max())
            maximum = float(np.ma.abs(y - x).max())
            flux_differences[name] = {"max_abs_difference": maximum,
                                      "relative_to_reference_peak": maximum / peak if peak else 0.,
                                      "agrees_at_1e_6": bool(np.ma.allclose(x, y, rtol=1.e-6, atol=1.e-6))}
        assert np.any(np.abs(a["ECO_HeatS_col"][:]) > 0)
        assert np.any(np.abs(a["ECO_LE_col"][:]) > 0)
    report = {"result": "PASS", "areas_m2": [1, 6], "totals_scale_by_area": list(totals),
              "invariant_geometry": list(geometry), "invariant_history": list(densities),
              "coupled_flux_diagnostics": flux_differences}
    (work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
