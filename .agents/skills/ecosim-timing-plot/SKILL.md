---
name: ecosim-timing-plot
description: "Visualize EcoSIM runtime-profiler output (timing/summary.csv and the per-step timing/time-<run-start>.txt) as a self-contained HTML report: total time per region, seconds per timestep, and per-region timestep traces. Use when asked to plot, chart, visualize, or summarize EcoSIM timing or profiling results, find which process dominates run time, or look for slow timesteps (history/restart write spikes)."
---

# EcoSIM Timing Plot

## Use When

- A run was made with `do_timing=.true.` (and optionally `do_timing_detail=.true.`) in `&ecosim`
  and you want to see where wall time goes.
- You need to compare process costs (`TranspNoSalt`, `WAT`, `NIT`, `PlantModel`, ...) or spot
  timesteps that are slower than the rest.

Not for scientific output; for history files use the ecosim-agent skills.

## Inputs

The standalone driver writes to `<cwd>/<case_name>_outputs/timing/`, where `<cwd>` is the
directory the executable was launched from:

| File | Written | Content |
|---|---|---|
| `summary.csv` | at normal shutdown only | per region: calls, total/mean/min/max seconds, % of whole run |
| `time-<YYYYMMDD-HHMMSS>.txt` | each step, if `do_timing_detail=.true.` | `step,timer,calls,seconds`; setup rows use step `-1` |

Older runs may have a single `time.txt`; it is read if no `time-*.txt` exists.

## Workflow

Run from the repository root. The script uses only the Python standard library.

```bash
# a <case>_outputs/ directory, its timing/ directory, or summary.csv
python3 .agents/skills/ecosim-timing-plot/scripts/plot_timing.py <case>_outputs

# choose a specific detail file, an output path, or a title
python3 .agents/skills/ecosim-timing-plot/scripts/plot_timing.py <case>_outputs/timing \
    --detail <case>_outputs/timing/time-20261008-110529.txt -o report.html --title "dryland 10 d"

# summary only
python3 .agents/skills/ecosim-timing-plot/scripts/plot_timing.py <case>_outputs --no-detail
```

By default the report is written to `timing/timing_report.html` next to `summary.csv`, and the
newest `time-*.txt` is used (the timestamp sorts chronologically). When several runs share a
`timing/` directory, pass `--detail` explicitly: `summary.csv` is replaced by every run, so it
belongs to the latest launch only.

To generate timing data in a scratch copy of a case, add to `&ecosim`:

```
do_timing=.true.
do_timing_detail=.true.
```

and follow the scratch-mirror approach of the `ecosim-smoke-run` skill so committed outputs are
not overwritten.

## Report contents

1. Headline tiles: whole-run seconds, timestep count, mean seconds per timestep, timestep share
   of the run, initialization time.
2. Horizontal bars of total seconds per region with % of run and a whole-run reference line,
   plus a table view.
3. Seconds per timestep (`Timestep` region) over the run.
4. Seconds per timestep for the 7 costliest regions other than `Timestep`, with legend,
   crosshair tooltip, and a table view. Remaining regions are named in the note.

Long runs are binned to at most about 1500 points per series. Bin widths are whole days (multiples
of 24 hourly steps) so daily events such as history writes do not alias. Each point is the mean
over the bin; the tooltip also shows the bin maximum, so isolated spikes stay visible.

## Interpreting

- Regions are inclusive and nested. Every call in `AdvanceModelOneYear` has its own region:
  - `AnnualSetup` contains `InitControlParms`, `STARTS`, `ReadPlantInfo`, `ReadClimSoilForcing`,
    `STARTQ`, `STARTE` and `ReadSoilWarmingTref` (the conditional ones appear only when they run).
  - `Timestep` contains `WTHR`, `EcoSIMOneStep`, `HistUpdate`, `HistUpdateHbuf`, `ClockUpdate`,
    `HistoryIO` and `RestartWrite`.
  - `EcoSIMOneStep` contains the process regions (`HOUR1`, `WAT`, `NIT`, `PlantModel`,
    `soluteModel`, `TranspNoSalt`, `TranspSalt`, `EROSION`, `REDIST`, `Diagnostics`).
  - Daily calls (`DAY`, `SetAnnualAccumlators`, `WriteBBGCForc`) sit outside `Timestep`.
  Percentages therefore overlap and do not sum to 100.
- Whole-run time minus `Initialization`, `AnnualSetup`, `Timestep`, the daily regions,
  `RestartRead` and `SummarizeTracerMass` is uninstrumented driver work (year loop, cleanup).
- Daily regions run outside any step; in the detail file they are recorded with the next step
  (`DAY` with the day's first hour, `WriteBBGCForc` with the following day's first hour).
- A region that runs only on some steps (e.g. `RestartWrite`) counts as zero on other steps
  in the per-step plot.
- `RestartRead` is recorded in step 0's detail rows, not under setup.
- Wall-clock times depend on machine load; compare runs made on the same machine and build.

## Limits

- An aborted run leaves `summary.csv` empty; the script stops with an error. The detail file
  still holds every completed step and can be inspected directly.
- Only the standalone `ecosim.f90.x` driver writes timing files.
