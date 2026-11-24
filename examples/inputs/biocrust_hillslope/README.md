# biocrust_hillslope input data

Extends the single-cell `examples/inputs/biocrust` bare-soil/biocrust case into a
5-position hillslope transect, using the multi-topo-unit grid/pft_mgmt layout
demonstrated in `examples/inputs/DaLake`. Generated 2026-10-01.

## Files

- `biocrust_hillslope_grid_20261001.nc` (+ `.cdl`): extends
  `examples/inputs/biocrust/BareSoil_grid_20260831.nc` from `ntopou=nrow=1` to
  `ntopou=nrow=5`.
- `biocrust_hillslope_pft_mgmt_20261001.nc` (+ `.cdl`): extends
  `examples/inputs/biocrust/biocrust_pft_mgmt.nc` from `ntopou=1` to `ntopou=5`.
- Climate forcing and GHG forcing are **not duplicated** here; the run
  namelist (`examples/run_dir/biocrust_hillslope/biocrust_hillslope.namelist`)
  points back at `examples/inputs/biocrust/BareSoil_clim_20250203.nc` and the
  shared `input_data/fatm_hist_GHGs_1750-2025.nc`, since atmospheric forcing
  is `ngrid`-dimensioned (one weather station for the whole grid) and is
  unaffected by the number of topo units.

## Hillslope layout

Single column (`NH1=NH2=1`), five rows, one grid-file topo-unit entry per
row (`NV1=NV2=1..5`) because `readimod.F90` writes each entry's soil data to
the single point `(NV1,NH1)` rather than broadcasting across a row range —
confirmed by tracing `f90src/IOutils/readimod.F90` lines ~393-460. Positions,
north (row 1) to south (row 5):

| Row | Position   | SL0 (deg) |
|-----|------------|-----------|
| 1   | summit     | 2         |
| 2   | shoulder   | 5         |
| 3   | backslope  | 10        |
| 4   | footslope  | 5         |
| 5   | toeslope   | 2         |

`DVI` = 1 m per row (5 m transect length), `DHI` = 1 m (unchanged column
width), `ASPX` = 90 deg at every position (one planar slope line, not a 2D
catchment).

## What was copied unchanged vs. what was invented

- **Copied unchanged** across all 5 positions: every `ntopou x nlevs` soil
  profile (`CDPTH` ... `RSPmL`), the `ntopou`-scalar soil/litter chemistry
  (`PSIFC`, `PSIWP`, `ALBS`, `PH0`, `RSCf`/`RSNf`/`RSPf`/..., `IXTYP1`,
  `IXTYP2`, `NUI`, `NJ`, `NL1`, `NL2`, `ISOILR`), and every `ngrid`-level
  scalar (`ALATG`, `ALTIG`, `ATCAG`, ..., `fAeroMicb`, `fAeroMoss`,
  `fAeroLich`). **This is a simplifying assumption, not a measured hillslope
  soil survey**: the same parent material, litter chemistry, and atmosphere
  are assumed at every position; only slope differs.
- **Varied intentionally**: `SL0` only, to give the lateral water/heat
  redistribution routines (which key off slope) something to act on.
- **pft_mgmt**: one management entry (`NV1=1,NV2=5`) applies the same
  `lich22` + `moss22` community to the whole hillslope, reusing the
  grouping convention `DaLake_pft_20250902.nc` uses for its entries that
  span multiple rows (confirmed by tracing the `NH1,NH2,NV1,NV2` read loops
  in `f90src/IOutils/PlantInfoMod.F90` and `ReadManagementMod.F90`, which
  loop `DO NX=NH1,NH2; DO NY=NV1,NV2` when assigning PFT/management data,
  unlike the grid file's single-point assignment).

If real topographic/soil survey data for a specific hillslope site becomes
available (differing soil depth, texture, or litter inputs by position), it
should replace the uniform profile replication above — this file only
establishes the structural (ntopou=5) extension the user asked for.

## Validation performed

- Generated via a one-off script
  (not checked into the repo) using `netCDF4`, reading the real `.nc`
  binaries as ground truth (the committed `BareSoil_grid_20260831.nc.cdl`
  text dump was found to be stale relative to its `.nc`: it lists
  `fAeroLiveMB`/`fAeroDeadMB`/`fAeroDeadNMB`/`fAeroDOM`, which are not
  actually present; the real file has a single `fAeroMicb` instead).
- Read back with `netCDF4` to confirm dimensions, the `NV1/NV2/SL0` ladder,
  and that unpopulated `pft_mgmt` entries (rows 2-5) resolve to the masked
  fill value, which the read loops in `PlantInfoMod.F90` use as the
  stop condition.
- Ran `examples/run_dir/biocrust_hillslope/` through the `ecosim-smoke-run`
  skill against the current build; see that run's log for the outcome.
