# utils

|File             | Description                                  |
|-----------------|----------------------------------------------|
|abortutils.F90   | tools for safe abortion of the model         |
|data_kind_mod.F90| setting byte size for different data types   |
|TestMod.F90      | code to support regression tests             |
|EcoSimConst.F90  | module defines constant ecosim parameters    |
|fileUtil.F90     | module to handle files                       |
|getfilename.c    | C code to handle filenames                   |
|MiniMathMod.F90  | module with functions to handle some simple maths safely|
|timings.F90      | Code to do timing                            |

### Runtime profiling (standalone EcoSIM driver)

Set `do_timing=.true.` in `&ecosim` to write
`<case>_outputs/timing/summary.csv` at normal shutdown. It reports call counts,
wall seconds (total, mean, minimum, maximum), and percentages of elapsed run time.
The run clock starts after namelist reading and output-directory setup, and ends
following model cleanup. Named regions include initialization, annual setup,
every subroutine call in `AdvanceModelOneYear` (named after the routine), the
processes in `Run_EcoSIM_one_step`, diagnostics, history I/O, restart
reads/writes, and timesteps.
Uninstrumented driver work still contributes to whole-run elapsed time.

Set `do_timing_detail=.true.` as well for optional `timing/time-<run-start>.txt`
CSV records (`step,timer,calls,seconds`), where `<run-start>` is the local launch
time as `YYYYMMDD-HHMMSS`. Steps are zero based and continue across years;
setup records use step `-1`. The detail format replaces the previous wide table.
The summary is replaced for each new run. Each detail-enabled run writes its own
file, so earlier detail files are kept; two launches within the same second
share a name and the later one replaces the earlier. Copy the summary before
repeating a case.

For new regions, declare `type(timer_stamp_type) :: stamp`, then call
`start_timer(stamp)` and `end_timer('RegionName',stamp)`. Each simultaneously
active region needs its own stamp. `end_timer_loop()` writes and resets timestep
records without resetting run totals. Timer names retain their full length.
Nested timings are inclusive, so totals and percentages overlap. Accumulation is
serial and is not thread safe; instrument outside parallel regions. Runs that
abort before normal finalization may not produce a complete summary.
