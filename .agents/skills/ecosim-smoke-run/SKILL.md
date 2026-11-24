---
name: ecosim-smoke-run
description: "Run every EcoSIM example case for 20 simulated days by default in an isolated scratch mirror and report which ones fail. Use when asked to smoke-test, sanity-check, or regression-screen the example cases after a Fortran change, to verify a build still runs end to end, or to check example namelists for startup and early-timestep failures. Also use before opening a PR that touches f90src/."
---

# EcoSIM Example Smoke Run

## Use When

- A Fortran change under `f90src/` or `drivers/` needs a quick "does everything still run" check.
- You want to know which example cases break, before spending hours on a long simulation.
- You are screening example namelists after changing an input dataset or parameter file.

Do **not** use this to judge scientific correctness. See [Limits](#limits).

## What It Does

For each eligible namelist under `examples/run_dir/*/` (see exclusions below):

1. Stages it in a scratch mirror so committed outputs and restart pointers are never touched.
2. Forces a short run (`stop_option='ndays'`, `stop_n=DAYS`) from the case's own `start_date`.
3. Runs the model, capturing a per-case log.
4. Verifies exit code **and** that the run reached `nstep = DAYS x 24`.
5. Scans every history file for genuine NaN/Inf and implausible soil temperature.

## Workflow

Confirm the build is current first — a stale binary invalidates the whole exercise:

```bash
cmake --build build/<configuration> --parallel 8
```

Then run the smoke test from the repository root:

```bash
.agents/skills/ecosim-smoke-run/scripts/run_smoke.sh            # all cases, 20 days
.agents/skills/ecosim-smoke-run/scripts/run_smoke.sh -d 2       # faster screen
.agents/skills/ecosim-smoke-run/scripts/run_smoke.sh biocrust dryland   # subset
```

The runner prefers the executable under `drivers/ecosim/` over installed copies
and resolves the build directory before changing into scratch case directories.

Options: `-d DAYS` (default 20), `-b BUILD_DIR` (default: newest under `build/`),
`-o WORKDIR` (default: a `mktemp -d`), `-t CPU_SECS` per-case cap (default 1800).

The default run must complete 480 hourly steps (20 days x 24 hours). Use
`-d DAYS` when the user explicitly requests a different duration.

### Excluded namelists

**Exclude `blodgett/Blodget.warming` from routine smoke runs**, including runs
that explicitly select the `blodgett` directory. It has an existing SIGILL
failure and is reserved for separate investigation. The runner removes only
its scratch namelist; the original case and other Blodgett cases are preserved.

**Exclude `dryland/dryland2` from the smoke run.** It is not a standalone case:
it starts from `finidat='./dryland_maize.ecosim.r.2004-01-01-000000.nc'`, a
restart produced by running `dryland` itself out to 2004, and that file is not
tracked in the repository. Its result says nothing about the build either way.

The runner automatically excludes both namelists in its scratch mirror, including
when their directories are explicitly selected. It prints the exclusion reasons.
Original namelists and local restart files are preserved.

Exit status is 0 only if every case completed its days *and* the output scan is
clean. The script prints a `CASE / RC / NSTEP / RESULT` table, then the scan.

To re-scan existing output without re-running:

```bash
python3 .agents/skills/ecosim-smoke-run/scripts/check_history.py <dir-with-h0-files>
```

Historical 10-day reference run: all 17 namelists across 13 case directories pass at 10 days in
roughly 45 s total on an arm64 Mac. That count includes `dryland2`, which is
now excluded (see above), so expect 16 rows from a compliant run.

## Reading The Result

| Symptom | Meaning |
|---|---|
| `RC != 0` | crash or CPU cap hit — read the case log |
| `RC = 0` but `NSTEP` short | model stopped early without failing; treat as a failure |
| `FAIL` in the output scan | model produced NaN/Inf, or temperature left the plausible range |
| `(N all-fill vars)` | **not** a fault — plantless cases fill every PFT variable |
| `NOT CHECKED` lines | a variable's `_FillValue` is itself NaN, so it cannot be judged |

Benign log noise, not failures: `CHECK_VAR: variable ... is not on initial
dataset`, and `Do error checking on file` from reading a `finidat` restart.

## Limits

A short run exercises input reading, initialization, and the first `DAYS x 24`
timesteps. It says nothing about:

- **Scientific correctness.** Passing means "did not crash and produced finite
  numbers," nothing more. Report it that way.
- **Long-term behavior.** The run starts at each case's `start_date` with all
  spinup cycling skipped, and several cases start decades before their forcing
  period (biocrust and bare_soil at 1800, DaLake at 1930), so this is
  early-spinup behavior only.
- **Conservation.** C/N/P, water, and energy balance drift needs a multi-year
  run or `regression-tests/` (`mtest`, `rtest`).

Never update a regression baseline on the strength of a passing smoke run.

## Implementation Notes

Read `references/pitfalls.md` before modifying either script. It records the
non-obvious failures that this procedure is built around — in particular why the
scratch mirror must use symlinks at matching depth, and why a NaN scan that
filters on magnitude silently reports every file as clean.
