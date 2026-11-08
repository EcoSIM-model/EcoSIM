---
name: ecosim-restart-run
description: >-
  Start an EcoSIM simulation from an existing checkpoint (<case>.ecosim.r.*.nc)
  of a long run, either as a branch (new case name, fresh history files, extra
  history fields or changed inputs) or as a continuation of the parent case.
  Use for diagnostic or sensitivity runs over a few years of a mature
  simulation without repeating the spin-up, or to resume a stopped run.
---

# EcoSIM Restart Run

Use this to run from a mature checkpoint. Use `ecosim-restart-test` to test
restart reproducibility, and `ecosim-smoke-run` for short startup screens.

## Choose the mode

| | `branch` (default) | `continue` |
|---|---|---|
| Namelist | `continue_run=.false.`, `finidat='<checkpoint>'` | `continue_run=.true.`, `rpointer.esim` |
| Case name | must differ from the parent | parent's |
| History | new files from the checkpoint date; `hist_fincl*` may change | appended to the parent's history files; field list is read from the `.rh*` files and cannot change |
| Needs | checkpoint only | checkpoint, matching `.rh*` files, and the history files they name |
| Inputs | may differ (sensitivity runs) | should match the parent |

Branch is the right choice for diagnostics and experiments. In continue mode the
model reopens the history file named in the `.rh*` file and keeps writing to it,
so it must be a copy, and adding fields stops the run (`ntapes` and the field
list must match the checkpoint).

## Stage and run

From the repository root:

```sh
python3 .agents/skills/ecosim-restart-run/scripts/stage_restart_run.py \
  --namelist <case_dir>/<case>.namelist \
  --restart  <case_dir>/<case>.ecosim.r.2006-01-01-000000.nc \
  --out      <scratch>/branch2006 \
  --case-name <case>.diag --stop-n 2 \
  --add-fields CanopyNLim_pft,Root_C_pft,RootAR_pft \
  --copy-input pft_file_in \
  --exe build/<configuration>/bin/ecosim.f90.x
```

- `--stop-n` / `--stop-option` request additional time from the checkpoint.
  The stager converts that endpoint to an absolute `nsteps` limit: the driver
  replays its clock from January 1 before loading the checkpoint, and those
  skipped hours also increment its stop counter. Intraday checkpoint times are
  preserved; calendar month/year increments retain the day where possible and
  otherwise use the target month's last day.
- `--copy-input KEY` copies that input into the run directory so it can be
  edited there; every other input is referenced by absolute path and is never
  modified.
- `--set KEY=VALUE` writes extra `&ecosim` assignments verbatim.
- Without `--exe` the script only stages. With it, the script runs the model
  and fails unless the exit status is zero, the log contains
  `restFile_open Successfully`, and the log has no `ENDRUN` or Fortran runtime
  error. It also requires the expected final timestep and dated checkpoint;
  merely loading a checkpoint or stopping early is not a pass.
- The staged directory gets `runfile.namelist` and `restart_run.json`, which
  records the checkpoint, inputs and forcing-year mapping.

## What branch mode rewrites, and why

1. `start_date` is set to the checkpoint date. The driver otherwise starts from
   the parent `start_date`, because `is_restart()` is false for a branch.
2. `forc_periods`, `npxs` and `npys` are cut so that the first simulated year uses
   the forcing year the parent used for the checkpoint year. The script expands
   the parent periods year by year, as the driver loop does, and rebuilds them
   from the checkpoint year onward. Example: a parent starting in 1996 with
   `forc_periods = 2006,2015,1,2006,2024,1` and a 2006 checkpoint gives
   `forc_periods = 2006,2024,1`.
3. Transient PFT management (`pft_dflag = 1`) and soil management are selected
   by forcing year, so they follow the cut periods automatically. Planting
   records are skipped while the PFT is alive, so a checkpoint with live plants
   is not replanted.
4. The case name must differ; otherwise `restFile_getfile` stops with
   `Must change case name on branch run` unless `brnch_retain_casename` is set.

## Check the result

- In the log: `Reading restart dataset`, then `restFile_open Successfully`.
  For continue mode also `read restart/checkpoint info file: ecosim_rst <date>`.
- The new history file starts at the checkpoint date
  (`<case>.ecosim.h0.<YYYY>-01-01-00000.nc`).
- Compare state variables (for example LAI or soil C) on the first output days
  with the parent's history at the same dates. With unchanged inputs, they should
  match until the first input that differs takes effect.

## Pitfalls

- A continue-mode run without the `.rh*` file stops in `hist_restart_ncd` with
  `GETFIL: FAILED to get ./<case>.ecosim.rh0.<date>.nc`.
- An edited management NetCDF must keep blank-padded strings: the fertilizer
  reader stops on empty records filled with NUL bytes (`Bad real number in item 1
  of list input` at `ReadManagementMod.F90`).
- A continued or branched run from a legacy checkpoint is only as good as the
  state that checkpoint stores; use `ecosim-restart-test` to check restart
  equivalence when that matters.
