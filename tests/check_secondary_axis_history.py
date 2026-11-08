"""Check full-model mean secondary-axis history against checkpoint geometry.

Run the FireCA example through its first growing season with hourly output.
Compare the final history sample to the final root length/count state, requiring
both positive means and layers without axes. SP values are checked separately
by check_sp_phenology.py.
"""

import argparse
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from netCDF4 import Dataset

from check_annual_restart import ROOT, stage
from check_restart_exactness import compare


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--workdir", type=Path, required=True)
    parser.add_argument("--days", type=int, default=180)
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    work = args.workdir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    steps = args.days * 24
    source = ROOT / "examples/run_dir/FireCA/FireCA.namelist"
    text = re.sub(r"(?m)![^\n]*$", "", source.read_text())

    def absolute_input(match):
        candidate = source.parent / match[1].rstrip()
        return "'" + str(candidate.resolve()) + "'" if candidate.is_file() else match[0]

    text = re.sub(r"'([^'\n]+)'", absolute_input, text)
    # This example has a multiline history-field list. Replace the whole list.
    text = re.sub(r"(?mis)^\s*hist_fincl1\s*=.*?(?=^\s*[a-z]\w*\s*=|^\s*/)", "", text)
    for key, value, group in (
        ("finidat", "''", "ecosim"), ("continue_run", ".false.", "ecosim"),
        ("forc_periods", "1980,1989,1", "ecosim"),
        ("hist_nhtfrq", "1", "ecosim"),
        ("hist_mfilt", str(steps + 1), "ecosim"),
        ("hist_fincl1", "'Root2ndMeanLen_pvr','Root2nd_AxesNumL_pvr'", "ecosim"),
        ("stop_option", "'nsteps'", "ecosim_time"),
        ("stop_n", str(steps), "ecosim_time"), ("rest_opt", "'never'", "ecosim_time"),
    ):
        text = stage.set_value(text, key, value, group)
    (work / "case.namelist").write_text(text)
    with (work / "run.log").open("w") as log:
        result = subprocess.run([str(executable), "case.namelist"], cwd=work,
                                stdout=log, stderr=subprocess.STDOUT)
    body = (work / "run.log").read_text(errors="replace")
    completed = re.findall(r"wrote out restart data at nstep\s*=\s*(\d+)", body)
    assert result.returncode == 0 and "ENDRUN" not in body, work
    assert completed and int(completed[-1]) == steps, work
    checkpoint, = work.glob("*.ecosim.r.*.nc")
    history, = work.glob("*.ecosim.h0.*.nc")
    assert not compare(checkpoint, checkpoint), checkpoint
    assert not compare(history, history), history
    positive = empty = 0
    with Dataset(checkpoint) as nc, Dataset(history) as hist:
        assert hist["Root2ndMeanLen_pvr"].units == "m"

        def ordered(name, dimensions):
            var = nc[name]
            return np.transpose(var[:], [var.dimensions.index(d) for d in dimensions])

        lengths = ordered("Root2ndLen_rpvr", ("pft", "rootyps", "levsoi", "rootaxs"))
        counts = ordered("Root2ndXNumL_rpvr", ("pft", "rootyps", "levsoi"))
        for pft, population in enumerate(nc["PlantPopuLive_pft"][:]):
            if np.ma.is_masked(population) or population <= 0:
                continue
            axes = int(nc["NumStructuralRootAxes_pft"][pft])
            total_length = np.asarray(lengths[pft, 0, :, :axes].sum(axis=-1))
            count = np.asarray(counts[pft, 0, :])
            expected = np.divide(total_length, count, out=np.zeros_like(total_length), where=count > 0)
            actual = hist["Root2ndMeanLen_pvr"][-1, :, pft]
            assert not np.any(np.ma.getmaskarray(actual)), pft
            np.testing.assert_allclose(actual, expected, rtol=1.e-6, atol=1.e-9)
            np.testing.assert_allclose(hist["Root2nd_AxesNumL_pvr"][-1, :, pft], count, rtol=1.e-6)
            positive += int(np.count_nonzero(expected > 0))
            empty += int(np.count_nonzero(count == 0))
    assert positive > 0 and empty > 0, (positive, empty)
    report = {"result": "PASS", "positive_layer_means": positive, "layers_without_axes": empty}
    (work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
