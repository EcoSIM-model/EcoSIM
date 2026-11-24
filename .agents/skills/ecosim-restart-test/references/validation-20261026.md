# Initial 10-day restart validation: 2026-10-26

Built with VSCode's bash environment and agf13 (GCC 13), Darwin arm64,
double precision, Release. Ran the default suite: 240 continuous steps versus
120 steps followed by a continuation to step 240. All 18 eligible cases
completed all three model executions; 11 comparisons passed, 7 failed.

| Case | Result | Unequal variables/findings |
|---|---|---:|
| DaLake/DaLake | PASS | 0 |
| Fen/FenStordIsland | PASS | 0 |
| FireCA/FireCA | FAIL | 2 |
| Pond/Pond | FAIL | 175 |
| RiceUSTWT/riceUSTWT | PASS | 0 |
| SatePhenol/SatePhenol_test | FAIL | 127 |
| bare_soil/BareSoil | FAIL | 1 |
| biocrust/biocrust | FAIL | 1 |
| biocrust_hillslope/biocrust_hillslope | FAIL | 1 |
| blodgett/Blodget.dfir | PASS | 0 |
| blodgett/Blodget.ndlf.lig | PASS | 0 |
| blodgett/Blodget.ndlf.nolig | PASS | 0 |
| climeConst/climeConst | PASS | 0 |
| dryland/dryland | PASS | 0 |
| lake/lake | FAIL | 1 |
| macpore_test/Blodget.macpore | PASS | 0 |
| macpore_test/Blodget.nomacpore | PASS | 0 |
| radiation_test/Radiation_test | PASS | 0 |

Skipped: Blodget.warming (existing SIGILL failure) and dryland2 (external
2004 checkpoint dependency). The additional dryland comparison with an
intraday split at hour 125 passed with zero differences.

FireCA differs in cumulative gas-balance diagnostics. BareSoil, biocrust and
lake differ in the cumulative gas-balance checkpoint field. Biocrust hillslope
differs in history cumulative organic-matter dry deposition. Pond and
SatePhenol differ in multiple physical/chemical states and history outputs.
These are observed differences, not established root causes. No tolerances
or baselines were changed; these failures were unresolved at this initial run.

The skill validator passed. Comparator checks confirmed equal/fill data pass,
a small numeric perturbation fails, and a NaN valid value fails.

Initial generated summary/logs (local, not versioned):
`ecosim-agent/result/Blodget_dfir_20261026/restart_validation_10day/`.
Raw NetCDF outputs remain at `/private/tmp/ecosim-restart-10day-20261026`.
The validation record is historical; rerun before drawing conclusions about
another build. Scientific correctness and long-term behavior are untested.

## Follow-up fixes and final validation

After retaining geometry and deposition accumulators, correcting the prescribed
phenology active-PFT count, and restoring the active snow-layer count, the same
10-day Release suite passes **18/18 eligible cases with zero checkpoint/history
differences**. The two documented exclusions are unchanged.

Pond, Biocrust hillslope, and SatePhenol passed targeted Debug/Release checks;
FireCA, BareSoil, biocrust, and lake passed targeted Debug checks. Intraday
restarts at hour 125 passed for Pond, hillslope, SatePhenol and FireCA.
A FireCA checkpoint with the new snow-layer-count field removed also recovered
exactly through the compatibility reconstruction. Validation used VSCode's bash
and agf13 (GCC 13), Darwin arm64, double precision. No tolerances or baselines
were relaxed. This establishes equivalence for these short trajectories only.

The shared cause of the last four diagnostic failures was omitted snow-layer
count: restored snow gases were left out of the initial gas-mass sum.
Final generated reports (local, not versioned) are under
`ecosim-agent/result/Blodget_dfir_20261026/restart_fix_fireca_20261026/`.
