"""Exercise standalone prescribed PFT inputs without a plant-management file.

Requires numpy and netCDF4. All generated inputs and short simulations stay in
a new --workdir. An annual branch checks changes to prescribed population density.
"""

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import numpy as np
from netCDF4 import Dataset, chartostring

from check_annual_restart import ROOT, stage
from check_restart_exactness import compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    source = ROOT / "examples/run_dir/SatePhenol/SatePhenol_test.namelist"
    text = re.sub(r"(?m)![^\n]*$", "", source.read_text())

    def absolute_input(match):
        candidate = source.parent / match[1].rstrip()
        return "'" + str(candidate.resolve()) + "'" if candidate.is_file() else match[0]

    text = re.sub(r"'([^'\n]+)'", absolute_input, text)
    # A non-unit cell area exposes a missing or duplicated area conversion
    # between the prescribed per-area stock and structural history output.
    grid_file = work / "grid-6m2.nc"
    shutil.copyfile(ROOT / "examples/inputs/blodgett/Blodget_grid_20251115.nc", grid_file)
    with Dataset(grid_file, "r+") as nc:
        nc["DHI"][:] = 2.
        nc["DVI"][:] = 3.
    with Dataset(grid_file) as nc:
        assert nc["DHI"][:].item() * nc["DVI"][:].item() == 6.
    for key, value, group in (
        ("grid_file_in", "'" + str(grid_file) + "'", "ecosim"),
        ("pft_mgmt_in", "'" + str(work / "does-not-exist.nc") + "'", "ecosim"),
        ("forc_periods", "2013,2014,1", "ecosim"),
        ("hist_nhtfrq", "1", "ecosim"),
        ("hist_fincl1", "'RootC_2nd_pvr','RootN_2nd_pvr','RootP_2nd_pvr',"
         "'Root2ndBiomC_pft','MycoBiomC_pvr','Root2nd_AxesNumL_pvr','Root2ndMeanLen_pvr'", "ecosim"),
        ("stop_option", "'ndays'", "ecosim_time"),
        ("stop_n", "2", "ecosim_time"),
    ):
        text = stage.set_value(text, key, value, group)
    with Dataset(ROOT / "input_data/ecosim_pftpar_20261021.nc") as nc:
        catalog = chartostring(nc["pfts"][:]).tolist()
        root_types = nc["IROOTJK"][:]
        mean_axis_lengths = np.array(nc["MFRUTL"][:], dtype=float)
        canopy_concentrations = np.array(nc["CCPOLT"][:], dtype=float)
    with Dataset(ROOT / "examples/inputs/blodgett/Blodget_sp_phenol_20261103.nc") as nc:
        monthly = {name: np.array(nc[name][0, 0, 0, :]) for name in ("LAI", "SAI", "HTOP")}

    def fixture(name, records, slots=2):
        path = work / (name + ".nc")
        with Dataset(path, "w", format="NETCDF3_CLASSIC") as nc:
            for dim, size in (("ntopou", 1), ("maxpfts", 2), ("nchar1", 10),
                              ("month", 12), ("year", None)):
                nc.createDimension(dim, size)
            for bound in ("NH1", "NV1", "NH2", "NV2"):
                nc.createVariable(bound, "i4", ("ntopou",))[:] = 1
            nc.createVariable("NZ", "i4", ("ntopou",))[:] = slots
            nc.createVariable("year", "i4", ("year",))[:] = [2013 + i for i in range(len(records))]
            names = nc.createVariable("pft_type", "S1", ("year", "ntopou", "maxpfts", "nchar1"))
            for i, record in enumerate(records):
                for j, pft in enumerate(record):
                    names[i, 0, j, :] = np.frombuffer(pft.encode().ljust(10, b"\0"), dtype="S1")
            population = nc.createVariable("PlantPopuLive_pft", "f8", ("year", "ntopou", "maxpfts"))
            population.units = "plants m-2"
            for i, record in enumerate(records):
                population[i, 0, :] = [((0.25 + .5 * i) if j == 0 else (3.5 + 2. * i)) if pft else 0.
                                       for j, pft in enumerate(record)]
            for field, values in monthly.items():
                var = nc.createVariable(field, "f4", ("year", "ntopou", "maxpfts", "month"))
                for i, record in enumerate(records):
                    for j, pft in enumerate(record):
                        var[i, 0, j, :] = values * (1 - 0.3 * j) if pft else 0
        # Validate dimensions, padded strings, slot counts and values after writing.
        with Dataset(path) as nc:
            assert nc["NZ"][:].tolist() == [slots]
            assert chartostring(nc["pft_type"][:])[:, 0, :].tolist() == records
            assert all(np.all(np.isfinite(nc[field][:])) for field in monthly)
        return path

    report = {}
    cases = (("two-pfts", [["ndlf00", "dfir00"]], 2, 2),
             ("mixed-root-types", [["ndlf00", "gr3a00"]], 2, 2),
             ("retained-slot", [["dfir00", ""]], 2, 1),
             ("bare", [["", ""]], 0, 0),
             ("too-many-slots", [["ndlf00", "dfir00"]], 3, None),
             ("name-after-blank", [["", "dfir00"]], 2, None),
             ("invalid-root-type-0", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-root-type-11", [["ndlf00", "dfir00"]], 2, None),
             ("minimum-axis-length", [["ndlf00", "dfir00"]], 2, 2),
             ("invalid-axis-length-zero", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-axis-length-negative", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-axis-length-short", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-axis-length-nan", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-axis-length-inf", [["ndlf00", "dfir00"]], 2, None),
             ("zero-canopy-concentration", [["ndlf00", "dfir00"]], 2, 2),
             ("invalid-canopy-concentration-negative", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-canopy-concentration-nan", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-canopy-concentration-inf", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-canopy-concentration-fill", [["ndlf00", "dfir00"]], 2, None),
             ("missing-population", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-population-zero", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-population-negative", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-population-nan", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-population-inf", [["ndlf00", "dfir00"]], 2, None),
             ("invalid-population-fill", [["ndlf00", "dfir00"]], 2, None))
    for name, records, slots, expected in cases:
        folder = work / name
        folder.mkdir()
        data = fixture(name, records, slots)
        if name == "missing-population" or name.startswith("invalid-population-"):
            with Dataset(data, "r+") as nc:
                if name == "missing-population":
                    nc.renameVariable("PlantPopuLive_pft", "omitted_population")
                else:
                    nc["PlantPopuLive_pft"][0, 0, 0] = {
                        "zero": 0., "negative": -1., "nan": np.nan,
                        "inf": np.inf, "fill": 9.969209968386869e36}[name.rsplit("-", 1)[1]]
        current = stage.set_value(text, "sp_phenol_in", "'" + str(data) + "'", "ecosim")
        axis_lengths = mean_axis_lengths.copy()
        concentrations = canopy_concentrations.copy()
        if name in ("mixed-root-types", "minimum-axis-length") or name.startswith("invalid-axis-length-"):
            traits = folder / "axis-length-traits.nc"
            shutil.copyfile(ROOT / "input_data/ecosim_pftpar_20261021.nc", traits)
            with Dataset(traits, "r+") as nc:
                if name == "mixed-root-types":
                    nc["MFRUTL"][catalog.index("ndlf34")] = .005
                    nc["MFRUTL"][catalog.index("gr3a34")] = .02
                    nc["CCPOLT"][catalog.index("ndlf34")] = .05
                    nc["CCPOLT"][catalog.index("gr3a34")] = .2
                else:
                    value = (.002 if name == "minimum-axis-length" else
                             {"zero": 0., "negative": -.01, "short": .001,
                              "nan": np.nan, "inf": np.inf}[name.rsplit("-", 1)[1]])
                    nc["MFRUTL"][catalog.index("ndlf34")] = value
            with Dataset(traits) as nc:
                nc.set_auto_mask(False)  # Invalid fixtures deliberately violate valid_min.
                axis_lengths = np.array(nc["MFRUTL"][:], dtype=float)
                concentrations = np.array(nc["CCPOLT"][:], dtype=float)
            current = stage.set_value(current, "pft_file_in", "'" + str(traits) + "'", "ecosim")
        if name == "zero-canopy-concentration" or name.startswith("invalid-canopy-concentration-"):
            traits = folder / "canopy-concentration-traits.nc"
            shutil.copyfile(ROOT / "input_data/ecosim_pftpar_20261021.nc", traits)
            value = (0. if name == "zero-canopy-concentration" else
                     {"negative": -.1, "nan": np.nan, "inf": np.inf,
                      "fill": -9999.}[name.rsplit("-", 1)[1]])
            with Dataset(traits, "r+") as nc:
                nc["CCPOLT"][catalog.index("ndlf34")] = value
            with Dataset(traits) as nc:
                nc.set_auto_mask(False)
                concentrations = np.array(nc["CCPOLT"][:], dtype=float)
            current = stage.set_value(current, "pft_file_in", "'" + str(traits) + "'", "ecosim")
        if name.startswith("invalid-root-type-"):
            traits = folder / "invalid-traits.nc"
            shutil.copyfile(ROOT / "input_data/ecosim_pftpar_20261021.nc", traits)
            invalid_type = int(name.rsplit("-", 1)[1])
            with Dataset(traits, "r+") as nc:
                nc["IROOTJK"][catalog.index("ndlf34")] = invalid_type
            with Dataset(traits) as nc:
                assert nc["IROOTJK"][catalog.index("ndlf34")] == invalid_type
            current = stage.set_value(current, "pft_file_in", "'" + str(traits) + "'", "ecosim")
        (folder / "case.namelist").write_text(current)
        with (folder / "run.log").open("w") as log:
            result = subprocess.run([str(executable), "case.namelist"], cwd=folder,
                                    stdout=log, stderr=subprocess.STDOUT)
        body = (folder / "run.log").read_text(errors="replace")
        if expected is None:
            message = ("missing PlantPopuLive_pft" if name == "missing-population"
                       else "zero-population SP slots require LAI=0" if name == "invalid-population-zero"
                       else "prescribed PlantPopuLive_pft must be finite and nonnegative" if name.startswith("invalid-population-")
                       else "prescribed PFT CCPOLT must be finite and nonnegative" if name.startswith("invalid-canopy-concentration-")
                       else "prescribed PFT MFRUTL must be finite and at least Root2ndTipLen4uptk" if name.startswith("invalid-axis-length-")
                       else "prescribed PFT IROOTJK must be between 1 and 10" if name.startswith("invalid-root-type-")
                       else "invalid PFT slot count" if name == "too-many-slots"
                       else "blank pft_type slots must follow")
            assert message in body, folder
            report[name] = "expected input rejection"
            continue
        assert result.returncode == 0 and "ENDRUN" not in body, folder
        steps = re.findall(r"wrote out restart data at nstep\s*=\s*(\d+)", body)
        assert steps and int(steps[-1]) == 48, folder
        checkpoint, = folder.glob("*.ecosim.r.1980-01-03-*.nc")
        assert not compare(checkpoint, checkpoint), folder  # Includes valid-value finite checks.
        with Dataset(checkpoint) as nc:
            assert nc["NP"][:].item() == expected
            assert nc["NP0"][:].item() == slots
            assert np.all(nc["PlantPopuLive_pft"][expected:slots] == 0)
            assert np.all(nc["IsPlantActive_pft"][expected:slots] == 0)
            assert np.all(nc["NumStructuralRootAxes_pft"][:expected] == 0)
            # Jan 2 uses the Dec/Jan interpolation; native area is the index
            # times cell area, while history returns the prescribed index.
            weight = 1.5 / 31.
            with Dataset(data) as prescribed:
                lai = np.array(prescribed["LAI"][0, 0, :expected, :])
                sai = np.array(prescribed["SAI"][0, 0, :expected, :])
            leaf_index = (.5 - weight) * lai[:, -1] + (.5 + weight) * lai[:, 0]
            stem_index = (.5 - weight) * sai[:, -1] + (.5 + weight) * sai[:, 0]
            np.testing.assert_allclose(nc["CanopyLeafArea_pft"][:expected], leaf_index * 6., rtol=1.e-7)
            np.testing.assert_allclose(nc["CanopyStemSurfArea_pft"][:expected], stem_index * 6., rtol=1.e-7)
            with Dataset(next(folder.glob("*.ecosim.h0.*.nc"))) as hist:
                np.testing.assert_allclose(hist["ECO_LAI_col"][-1, 0], leaf_index.sum(), rtol=1.e-6)
                np.testing.assert_allclose(hist["ECO_SAI_col"][-1, 0], stem_index.sum(), rtol=1.e-6)

            for j, name_in in enumerate(records[0][:expected]):
                # The grid's Koppen code replaces the deliberately different input suffix.
                loc = catalog.index(name_in[:4] + "34")
                assert nc["DATAPI"][j] == loc + 1
                kind = int(root_types[loc]) - 1
                # Jackson Table 1 live biomass/length and geometric profile, checked
                # against checkpoint values from the double-precision build.
                beta = [.943, .970, .950, .980, .967, .943, .982, .972, .972, .909][kind]
                population = .25 if j == 0 else 3.5
                length_per_area = [2.6, 4., 8.4, 6.1, 5.4, 112., 3.5, 4.1, 60.4, 7.4][kind] * 1.e3
                carbon_per_area = [.23, .13, .28, .50, .44, .95, .28, .33, .51, .34][kind] * .488 * 1.e3
                with Dataset(grid_file) as grid:
                    nl = int(grid["NJ"][:].item() + grid["NL1"][:].item())
                    nl += min(len(nc.dimensions["levsoi"]) - nl, int(grid["NL2"][:].item()))
                depth = nc["CumDepz2LayBottom_vr"][1:nl + 1, 0]
                fractions = -np.diff(np.r_[1., np.exp(np.log(beta) * depth * 100.)])
                np.testing.assert_allclose(nc["PlantPopuLive_pft"][j], population * 6., rtol=1.e-12)
                np.testing.assert_allclose(nc["RootTotLenPerPlant_pvr"][:nl, 0, j],
                                           length_per_area * fractions / population, rtol=1.e-10, atol=1.e-10)
                np.testing.assert_allclose(nc["PopuRootMycoC_pvr"][:nl, 0, j],
                                           carbon_per_area * fractions, rtol=1.e-10, atol=1.e-10)
                mean_length = axis_lengths[loc]
                axis_count = length_per_area * fractions * 6. / mean_length
                np.testing.assert_allclose(nc["Root2ndEffLen4uptk_rpvr"][:nl, 0, j], mean_length, rtol=1.e-12)
                np.testing.assert_allclose(nc["Root2ndXNumL_rpvr"][:nl, 0, j], axis_count, rtol=1.e-10)
                thickness = nc["DLYR_3D"][1:nl + 1, 2, 0]
                np.testing.assert_allclose(nc["DLYR_3D"][1:nl + 1, 0, 0] *
                                           nc["DLYR_3D"][1:nl + 1, 1, 0], 6.)
                history, = folder.glob("*.ecosim.h0.*.nc")
                with Dataset(history) as hist:
                    np.testing.assert_allclose(hist["POPN_pft"][-1, j], population, rtol=1.e-6)
                    assert hist["CanopyNonstOMConc_pft"].units == "g/gC"
                    canopy_concentration = hist["CanopyNonstOMConc_pft"][-1, j]
                    assert not np.any(np.ma.getmaskarray(canopy_concentration)), folder
                    np.testing.assert_allclose(canopy_concentration, concentrations[loc], rtol=1.e-6)
                    values = hist["RootC_2nd_pvr"][-1, :nl, j]
                    assert not np.any(np.ma.getmaskarray(values)), folder
                    np.testing.assert_allclose(values, carbon_per_area * fractions / thickness,
                                               rtol=1.e-6, atol=1.e-6)
                    np.testing.assert_allclose(hist["Root2ndBiomC_pft"][-1, j],
                                               carbon_per_area * fractions.sum(), rtol=1.e-6)
                    np.testing.assert_allclose(hist["Root2nd_AxesNumL_pvr"][-1, :nl, j], axis_count, rtol=1.e-6)
                    assert hist["Root2ndMeanLen_pvr"].units == "m"
                    values = hist["Root2ndMeanLen_pvr"][-1, :nl, j]
                    assert not np.any(np.ma.getmaskarray(values)), folder
                    np.testing.assert_allclose(values, mean_length, rtol=1.e-6)
                    for field in ("RootN_2nd_pvr", "RootP_2nd_pvr", "MycoBiomC_pvr"):
                        values = hist[field][-1, :nl, j]
                        assert not np.any(np.ma.getmaskarray(values)), (folder, field)
                        np.testing.assert_array_equal(values, 0.)
        report[name] = "PASS"

    # Change both PFT populations at January 1. A fresh-history
    # branch applies the new annual record after restoring the old year's state.
    annual_data = fixture("annual-population", [["ndlf00", "dfir00"], ["ndlf00", "dfir00"]])
    annual_text = stage.set_value(text, "sp_phenol_in", "'" + str(annual_data) + "'", "ecosim")
    for key, value, group in (("pft_mgmt_in", "'NO'", "ecosim"),
                              ("start_date", "'19801230000000'", "ecosim"),
                              ("hist_mfilt", "8760", "ecosim"),
                              ("stop_option", "'nsteps'", "ecosim_time"),
                              ("rest_opt", "'never'", "ecosim_time")):
        annual_text = stage.set_value(annual_text, key, value, group)

    def annual_segment(label, stop):
        folder = work / label
        folder.mkdir()
        (folder / "runfile.namelist").write_text(stage.set_value(annual_text, "stop_n", str(stop), "ecosim_time"))
        with (folder / "run.log").open("w") as log:
            subprocess.run([str(executable), "runfile.namelist"], cwd=folder,
                           stdout=log, stderr=subprocess.STDOUT, check=True)
        body = (folder / "run.log").read_text(errors="replace")
        steps = re.findall(r"wrote out restart data at nstep\s*=\s*(\d+)", body)
        assert "ENDRUN" not in body and steps and int(steps[-1]) == stop, folder
        return folder

    continuous = annual_segment("annual-continuous", 8976)  # Dec 30 + 10 days, leap year.
    first = annual_segment("annual-first", 8784)  # January 1, before reading the new year.
    branch = work / "annual-branch"
    with (work / "annual-branch.log").open("w") as log:
        subprocess.run([sys.executable, str(ROOT / ".agents/skills/ecosim-restart-run/scripts/stage_restart_run.py"),
                        "--namelist", str(first / "runfile.namelist"),
                        "--restart", str(next(first.glob("*.ecosim.r.*.nc"))),
                        "--out", str(branch), "--case-name", "population.branch", "--stop-n", "8",
                        "--stop-option", "ndays", "--exe", str(executable)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    final = sorted(continuous.glob("*.ecosim.r.*.nc"))[-1]
    differences = compare(final, sorted(branch.glob("population.branch.ecosim.r.*.nc"))[-1])
    differences = [d for d in differences if d.get("variable") not in
                   ("start_date", "timemgr_rst_start_ymd", "locfnh", "locfnhr")]
    # History time bounds use the model-year origin in days. Align the branch's
    # January origin to the reference's 1980 origin before exact comparison.
    reference_history = next(continuous.glob("*.ecosim.h0.1981-*.nc"))
    aligned_history = work / "annual-aligned-history.nc"
    shutil.copyfile(next(branch.glob("population.branch.ecosim.h0.*.nc")), aligned_history)
    with Dataset(aligned_history, "r+") as nc:
        nc["time_bounds"][:] = nc["time_bounds"][:] + 366.
    history_differences = compare(reference_history, aligned_history)
    (work / "annual-comparison.json").write_text(json.dumps(
        {"checkpoint": differences, "history": history_differences}, indent=2) + "\n")
    assert not differences and not history_differences
    with Dataset(final) as nc:
        np.testing.assert_allclose(nc["PlantPopuLive_pft"][:2], [.75 * 6., 5.5 * 6.], rtol=1.e-12)
        assert nc["NP"][:].item() == 2 and nc["NP0"][:].item() == 2
    with Dataset(next(first.glob("*.ecosim.h0.*.nc"))) as nc:
        np.testing.assert_allclose(nc["POPN_pft"][0, :2], [.25, 3.5], rtol=1.e-6)
    with Dataset(reference_history) as nc:
        np.testing.assert_allclose(nc["POPN_pft"][-1, :2], [.75, 5.5], rtol=1.e-6)
    report["annual-population-branch"] = "PASS"

    (work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
