---
name: ecosim-restart-test
description: >-
  Test EcoSIM restart reproducibility using existing examples by comparing
  uninterrupted 10-day simulations with stopped-and-resumed runs. Use to
  validate restart changes or regression-screen restart capability.
---

# EcoSIM Restart Reproducibility Test

Use this for restart equivalence; use `ecosim-smoke-run` for a startup/finite-output screen.
The default comparison is 10 continuous days versus 5 days, stop, and 5 more days.

## Run

Build in the repository's configured terminal/compiler environment, using a matching
CMake cache. Run from the repository root:

```sh
cmake --build build/<configuration> --target ecosim.f90.x --parallel 6
.agents/skills/ecosim-restart-test/scripts/run_restart.sh -b build/<configuration>
```

Options follow the smoke-run interface:

- `-d DAYS`: total duration, default 10; the split defaults to halfway.
- `-b BUILD_DIR`: configured build; otherwise choose the newest available executable.
- `-o WORKDIR`: new or empty scratch directory; otherwise create a temporary directory.
- `-t CPU_SECS`: per-model-process CPU cap, default 1800.
- `-s HOURS`: override the split; `-s 125` tests history accumulation across an intraday restart.
- Positional arguments: example directories, such as `blodgett biocrust dryland`; default all.

Requires Python 3, NumPy, and netCDF4. The runner prefers the executable under
`drivers/ecosim/` and refuses a binary older than the Fortran sources.

## Comparison and outputs

The runner invokes [the shared comparator](../../../tests/check_restart_exactness.py)
for each eligible namelist. Inputs are resolved to absolute paths; simulations run
in separate scratch directories. Existing case outputs and restart pointers are
never overwritten. The continuation receives the checkpoint, matching restart
history (`.rh*`), and history files from the first segment.

Check all three runs reached their expected final timestep and checkpoint.
Compare every final checkpoint and history variable with exact numerical equality,
including dimensions/dtypes and valid-value NaN/Inf checks. Numeric fill values
are excluded from finite checks; do not filter NaN by magnitude. NetCDF byte layout
and changing metadata timestamps are not compared. No regression baseline is updated.

Read the `CASE / RC / DIFFS / RESULT` table, `summary.json`, `summary.txt`, and each
case's `comparison.json`. A missing checkpoint, incomplete run, process failure,
nonfinite valid value, or any unequal variable is a failure. Continue the other
cases and report failures with their logs; do not relax comparison tolerances.

Exclude `blodgett/Blodget.warming` (existing SIGILL failure) and
`dryland/dryland2` (requires a separately generated 2004 checkpoint), as in the
smoke workflow. These exclusions are recorded in the summary. Other unavailable
inputs remain reported failures rather than silent skips.

For a mature checkpoint, invoke the comparator directly with `--seed PATH` and
`--days 10`. Both branches recover the same seed before testing a newly written
checkpoint. Legacy checkpoints may omit state and cannot establish equivalence
to their original uninterrupted run.

Report the toolchain, duration/split, cases passed/failed/skipped, and report path.
A pass demonstrates equivalence only for the tested trajectories and build; it
is not evidence of scientific correctness or long-run stability. If inputs or
model failures block a case, retain its logs and report that limit.

## Recorded validation

See the [initial 10-day run record](references/validation-20261026.md) for
observed passes and unresolved failures. Do not treat these failures as
additional exclusions or assume a later build has the same result.
