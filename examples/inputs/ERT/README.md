# ERT prescribed phenology inputs

`../../run_dir/ERT_SatePhenol/ERT6.namelist` selects these files for the
six-column SP example. The model reads the `.nc` files; paired `.nc.cdl`
files are editable text representations (`ncgen -o file.nc file.nc.cdl`).

| Columns | PFT code | Living population (plants/m²) | Canopy height (m) |
| --- | --- | --- | --- |
| 1–4 | `ndlf43` | 0.1 | 17 |
| 5–6 | `gr3s43` (C3 grass) | 300 | **0.5, provisional** |

The population densities and PFT codes retain the corresponding entries from
`ERT_pft_20250819.nc`. Each column has one active
PFT. Grid `IETYPG=0` preserves the supplied PFT suffix. Traits are read from
`../../../input_data/ecosim_pftpar_20261021.nc`; shared traits are unchanged.

`ERT6_sp_phenol_20261114.nc` supplies monthly LAI, SAI, height and population.
Its single year record repeats as a climatology for every forcing year. Tree
LAI, SAI and height are copied from the input selected by SatePhenol,
`../blodgett/Blodget_sp_phenol_20261103.nc`. Grass LAI and SAI were supplied by
the user on 2026-11-14; grass height remains a provisional assumption.

| Month | Tree LAI | Tree SAI | Grass LAI | Grass SAI |
| --- | ---: | ---: | ---: | ---: |
| Jan | 1.1852 | 0.3190 | 0.05 | 0.25 |
| Feb | 1.1821 | 0.3058 | 0.05 | 0.20 |
| Mar | 1.1554 | 0.3058 | 0.05 | 0.15 |
| Apr | 1.2433 | 0.3032 | 0.10 | 0.15 |
| May | 1.2922 | 0.3058 | 0.65 | 0.10 |
| Jun | 1.3341 | 0.3117 | 1.80 | 0.15 |
| Jul | 1.2296 | 0.3433 | 2.50 | 0.25 |
| Aug | 1.4118 | 0.3032 | 2.20 | 0.35 |
| Sep | 1.4343 | 0.3084 | 1.10 | 0.45 |
| Oct | 1.3941 | 0.3292 | 0.25 | 0.40 |
| Nov | 1.2721 | 0.3656 | 0.05 | 0.35 |
| Dec | 1.2218 | 0.3249 | 0.05 | 0.30 |

`ERT_grid_20250826.nc` retains the existing six-column ERT grid. Climate
forcing has been moved from the run folder into this directory without changing
its contents. `ERT6.namelist` selects `nldas2.ert.1979-2024.hour.nc`; the
alternatives `ERT_ERA5clim_1979-2024.hour.nc` and
`ERT.1979-2021.hour.clim.nc` are also stored here. The six topographic units share one meteorological grid group. The
forcing filename is inherited; its provenance and timestamp convention have
not been independently established. Existing soil and boundary settings have
been preserved, rather than recalibrated.

The namelist enables `ldo_sp_mode`, disables microbial/geochemical cycling,
and obtains PFT identity and population from `sp_phenol_in` instead of plant
management. It cold-starts case `ERT6_SP` and retains the existing 1930 start,
4-year duration and original forcing schedule. These model years use the
1980–1989 forcing cycle. Atmospheric concentrations use shared
`fatm_hist_GHGs_1750-2023.nc` on the model-year axis.

Run from a scratch copy with the same relative path depth (or absolute input
paths) to protect existing outputs and restart pointers. The separate
`../../run_dir/ERT/ERT.namelist` selects the eight-column full-model example.
It uses `grid_SOAP_8grid.nc` and `pft_SOAP_v1.nc`, copied unchanged from the
local `benchmark/smallset/ERT/` workspace, together with the shared current PFT
table. Both examples use paths relative to their own run directories.
