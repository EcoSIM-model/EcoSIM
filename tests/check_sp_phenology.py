"""Exercise standalone prescribed PFT inputs without a plant-management file.

Requires numpy and netCDF4. All generated inputs and short simulations stay in
a new --workdir. Annual restart trajectories use check_annual_restart.py.
"""

import argparse
import json
from pathlib import Path
import re
import subprocess

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
    for key, value, group in (
        ("pft_mgmt_in", "'" + str(work / "does-not-exist.nc") + "'", "ecosim"),
        ("forc_periods", "2013,2014,1", "ecosim"),
        ("stop_option", "'ndays'", "ecosim_time"),
        ("stop_n", "2", "ecosim_time"),
    ):
        text = stage.set_value(text, key, value, group)
    with Dataset(ROOT / "input_data/ecosim_pftpar_20261021.nc") as nc:
        catalog = chartostring(nc["pfts"][:]).tolist()
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
             ("retained-slot", [["dfir00", ""]], 2, 1),
             ("bare", [["", ""]], 0, 0),
             ("too-many-slots", [["ndlf00", "dfir00"]], 3, None),
             ("name-after-blank", [["", "dfir00"]], 2, None))
    for name, records, slots, expected in cases:
        folder = work / name
        folder.mkdir()
        data = fixture(name, records, slots)
        current = stage.set_value(text, "sp_phenol_in", "'" + str(data) + "'", "ecosim")
        (folder / "case.namelist").write_text(current)
        with (folder / "run.log").open("w") as log:
            result = subprocess.run([str(executable), "case.namelist"], cwd=folder,
                                    stdout=log, stderr=subprocess.STDOUT)
        body = (folder / "run.log").read_text(errors="replace")
        if expected is None:
            message = "invalid PFT slot count" if name == "too-many-slots" else "blank pft_type slots must follow"
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
            for j, name_in in enumerate(records[0][:expected]):
                # The grid's Koppen code replaces the deliberately different input suffix.
                assert nc["DATAPI"][j] == catalog.index(name_in[:4] + "34") + 1
        report[name] = "PASS"

    (work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
