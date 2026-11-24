"""Check root C/N/P concentration aggregation and renamed checkpoint recovery.

Exercise full plant growth and a frozen radiation-only SP case in scratch.
The SP continuations carry distinct root and mycorrhizal concentrations through
both current and legacy checkpoint variable names, isolating cache reconstruction.
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
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    current_name = "RootNonstructElmConc_pvr"
    legacy_name = "RootNonstructElmConc_rpvr"

    def fixture(source, name, hours, periods):
        text = re.sub(r"(?m)![^\n]*$", "", source.read_text())

        def absolute_input(match):
            candidate = source.parent / match[1].rstrip()
            return "'" + str(candidate.resolve()) + "'" if candidate.is_file() else match[0]

        text = re.sub(r"'([^'\n]+)'", absolute_input, text)
        text = re.sub(r"(?mis)^\s*hist_fincl1\s*=.*?(?=^\s*[a-z]\w*\s*=|^\s*/)", "", text)
        for key, value, group in (
            ("case_name", "'" + name + "'", "ecosim"), ("finidat", "''", "ecosim"),
            ("continue_run", ".false.", "ecosim"), ("forc_periods", periods, "ecosim"),
            ("hist_nhtfrq", "1", "ecosim"), ("hist_mfilt", "8760", "ecosim"),
            ("hist_fincl1", "'RootNonstBConc_pvr'", "ecosim"),
            ("stop_option", "'nsteps'", "ecosim_time"), ("stop_n", str(hours), "ecosim_time"),
            ("rest_opt", "'never'", "ecosim_time"),
        ):
            text = stage.set_value(text, key, value, group)
        return text

    def run(folder, text, hours):
        folder.mkdir(exist_ok=True)
        (folder / "case.namelist").write_text(text)
        with (folder / "run.log").open("w") as log:
            result = subprocess.run([str(executable), "case.namelist"], cwd=folder,
                                    stdout=log, stderr=subprocess.STDOUT)
        body = (folder / "run.log").read_text(errors="replace")
        steps = re.findall(r"wrote out restart data at nstep\s*=\s*(\d+)", body)
        assert result.returncode == 0 and "ENDRUN" not in body, folder
        assert steps and int(steps[-1]) == hours, folder
        checkpoint = sorted(folder.glob("*.ecosim.r.*.nc"))[-1]
        history, = folder.glob("*.ecosim.h0.*.nc")
        assert not compare(checkpoint, checkpoint), checkpoint
        assert not compare(history, history), history
        return checkpoint, history

    full = work / "full"
    text = fixture(ROOT / "examples/run_dir/FireCA/FireCA.namelist", "root-nonst-full", 48, "1980,1989,1")
    checkpoint, history = run(full, text, 48)
    positive = 0
    with Dataset(checkpoint) as nc, Dataset(history) as hist:
        var = nc[current_name]
        concentrations = np.transpose(var[:], [var.dimensions.index(d) for d in
                                               ("pft", "rootyps", "levsoi", "elmnts")])
        assert hist["RootNonstBConc_pvr"].units == "g/gC"
        for pft, population in enumerate(nc["PlantPopuLive_pft"][:]):
            if np.ma.is_masked(population) or population <= 0:
                continue
            expected = concentrations[pft, 0, :, :].sum(axis=-1)
            values = hist["RootNonstBConc_pvr"][-1, :, pft]
            assert not np.any(np.ma.getmaskarray(values)), pft
            np.testing.assert_allclose(values, expected, rtol=1.e-6, atol=1.e-8)
            positive += int(np.count_nonzero(expected > 0))
    assert positive > 0, positive
    report = {"full": "PASS", "positive_full_layers": positive}

    seed = work / "seed"
    text = fixture(ROOT / "examples/run_dir/radiation_test/Radiation_test.namelist",
                   "root-nonst-rad", 1, "2013,2014,1")
    checkpoint, _ = run(seed, text, 1)
    finals = []
    for variant in ("current", "legacy"):
        folder = work / variant
        folder.mkdir()
        for path in seed.iterdir():
            if path.suffix == ".nc" or path.name == "rpointer.esim":
                shutil.copy2(path, folder / path.name)
        with Dataset(folder / checkpoint.name, "r+") as nc:
            var = nc[current_name]

            def index(**indices):
                return tuple(indices.get(d, slice(None)) for d in var.dimensions)

            var[index(pft=0)] = 0.
            for layer, concentrations in enumerate(((.02, .003, .001), (.05, .007, .002))):
                for element, value in enumerate(concentrations):
                    var[index(pft=0, rootyps=0, levsoi=layer, elmnts=element)] = value
                    var[index(pft=0, rootyps=1, levsoi=layer, elmnts=element)] = 10. * (element + 1)
            if variant == "legacy":
                nc.renameVariable(current_name, legacy_name)
        current = stage.set_value(text, "continue_run", ".true.", "ecosim")
        current = stage.set_value(current, "stop_n", "2", "ecosim_time")
        final, history = run(folder, current, 2)
        with Dataset(history) as hist, Dataset(final) as nc:
            assert current_name in nc.variables and legacy_name not in nc.variables
            values = hist["RootNonstBConc_pvr"][-1, :2, 0]
            assert not np.any(np.ma.getmaskarray(values)), variant
            np.testing.assert_allclose(values, [.024, .059], rtol=1.e-6, atol=1.e-8)
        finals.append((final, history))
        report[variant + "_recovery"] = "PASS"
    assert not compare(finals[0][0], finals[1][0]), "checkpoint name variants differ"
    assert not compare(finals[0][1], finals[1][1]), "history name variants differ"
    (work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
