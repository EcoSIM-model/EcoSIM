# Standing-dead-only SatePhenol example

This synthetic thermal/water-exchange example reuses SatePhenol's Blodgett grid,
meteorological forcing, atmospheric forcing and PFT parameters. It cold-starts
with one `ndlf34` slot, zero living population and LAI, constant dead stem area
index 0.3 m²/m², and dead canopy height 17 m. The first simulation year (1980)
uses the 2013 meteorological forcing, following the parent's forcing convention.
The namelist defaults to two calendar years (17,544 hours from 1980-01-01
to 1982-01-01), with hourly output. The smoke command below uses its own
20-day default.

For a named prescribed slot with `PlantPopuLive_pft=0`, `LAI` must be zero,
`SAI` specifies standing-dead stem area, and `HTOP` specifies dead height.
For positive living population, these fields retain their existing live-canopy
meaning. Initial dead C comes from `WTSTDI=5000 gC/m²` in the selected shared PFT
file; initial N/P follow its woody-litter composition. The example does not add
dead biomass each day or simulate a mortality event. Prescribed live root C and
root length are zero. No dead-root stock is required for passive canopy exchange.

Exposed dead stems absorb radiation, exchange sensible heat, intercept rain/snow,
and evaporate/sublimate intercepted stores. They do not photosynthesize or
transpire. Heat capacity uses the outer 2.5 mm dead-stem shell (limited by existing
dead carbon and total stem area), plus water/snow stores. Interior dead carbon
remains in the biomass pools. Fully snow-buried material retains its temperature and capacity with zero
direct atmospheric exchange; tissue-to-snow conduction is not modeled.

The model reads `examples/inputs/SatePhenolStandingDead/standing_dead.nc`.
Its editable source is `standing_dead.nc.cdl` in the same directory. From the
repository root, regenerate it with:

```sh
ncgen -o examples/inputs/SatePhenolStandingDead/standing_dead.nc \
  examples/inputs/SatePhenolStandingDead/standing_dead.nc.cdl
```

To run this example without writing into its source directory, use the repository
smoke runner (replace the build directory for your platform):

```sh
.agents/skills/ecosim-smoke-run/scripts/run_smoke.sh \
  -b build/Darwin-arm64-double-Debug SatePhenolStandingDead
```

`CAN_H_pft`, `CAN_RN_pft`, and `CAN_LE_pft` show passive exchange;
`QVegTransp_pft` and `RootC_2nd_pvr` must remain zero. Checkpoints retain the
existing canopy temperature and intercepted-water/snow stores. This is a model
capability example, not a calibrated standing-dead forest simulation.

Check a completed run with:

```sh
python3 tests/check_sp_standing_dead.py --run-dir /path/to/scratch/case --hours 480
```

For a completed run using the namelist’s two-year duration, use `--hours 17544`.
