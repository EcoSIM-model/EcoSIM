# Dryland NH3Band + NO3Band example

This case derives from `../dryland/dryland.namelist`. It retains the maize–soybean
rotation, weather, soil/grid, plant traits, planting and harvest schedules. It
shares those inputs with dryland and selects its own soil-management file:

`../../inputs/dryland_band_only/dryland_band_only_soilmgmt_20261110.nc`

Only `NH3Band` and `NO3Band` supply fertilizer N. `NO3Band` retains its original
dryland dose; `NH3Band` supplies the remaining N so annual N totals match the
updated dryland example, including the additional 2 g N/m² requested there.

| Date | NH3Band (g N/m²) | NO3Band (g N/m²) | Total N (kg N/ha) |
| --- | ---: | ---: | ---: |
| 16 May 2001 | 10.54 | 4.27 | 148.1 |
| 16 May 2003 | 7.98 | 2.99 | 109.7 |
| 7 April 2005 | 9.84 | 3.92 | 137.6 |

There is no fertilizer application in 2002, 2004 or 2006. Application depth
remains 0.02 m and row spacing remains 0.76 m. Banded phosphate, tillage,
residue/manure settings and irrigation selectors retain the dryland values.
This is an experimental fertilizer-form comparison, not an observed management
record. Compared with dryland, equal total N does not imply equal crop response
or equal gaseous losses.

Editable JSON and CDL representations accompany the NetCDF. After editing them,
regenerate and validate the NetCDF: that is the file the model reads.

The case adds NH3 gas, dissolved NH3, NH3 surface-emission and cumulative NH3
budget-error fields to the original dryland history selections. The configured
run is six years; a short January smoke run does not reach the first fertilizer
application on day 136. Use a scratch copy for validation and run at least
140 days to cross that event. Use a binary containing the NH3 gas-source fix
in `ApplyFastSourceSink`.

From this directory, launch the configured full case with the executable located
in your matching build directory, for example:

```bash
../../../build/Darwin-arm64-double-Release/drivers/ecosim/ecosim.f90.x dryland_band_only.namelist
```

This writes simulation output and restart files in this directory; use a copy
when preserving previous results.

Validation on 10 November 2026: GNU Fortran 13.4.0 Release build; a scratch
140-day run completed 3360 hourly steps through the first fertilizer event.
History contained finite valid values, and cumulative amended N was 14.81 g N/m².
This is a runtime check, not validation of crop-yield accuracy.
