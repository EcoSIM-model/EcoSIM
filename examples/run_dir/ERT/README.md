# ERT full-model example

Run `ERT.namelist` from this directory (prefer a scratch mirror for testing).
It is an eight-column full plant/microbial/geochemical configuration using
`../../inputs/ERT/grid_SOAP_8grid.nc` and
`../../inputs/ERT/pft_SOAP_v1.nc`. The management file prescribes constant
vegetation (`pft_dflag=0`); its first record is reused by the native reader.
The SOAP grid/management files retain their original values from the local
`benchmark/smallset/ERT/` workspace.

Hourly forcing is `../../inputs/ERT/nldas2.ert.1979-2024.hour.nc`. The current
shared PFT parameter table and atmospheric concentrations are selected from
`../../../input_data/`. The supplied namelist cold-starts case `ERT` in 1930
for 40 years with the original forcing schedule.

The six-column SP alternative is
[../ERT_SatePhenol/ERT6.namelist](../ERT_SatePhenol/ERT6.namelist). Shared input
provenance is in [../../inputs/ERT/README.md](../../inputs/ERT/README.md).

Validation: paths, grid/management bounds, selected PFT identifiers and input
schemas were checked. The shared preflight currently supports annual plant
management only and rejects this valid constant-management format; a runtime
smoke test of this full-model example has not been performed for this addition.
