# ERT_SatePhenol

`ERT6.namelist` is the six-column prescribed-phenology case: tree canopy in
columns 1–4, C3 grass in columns 5–6. Inputs and their provenance are documented
in [../../inputs/ERT/README.md](../../inputs/ERT/README.md).

The default is a 4-year cold-start run named `ERT6_SP`. Run a scratch copy
when testing; existing ERT outputs, checkpoints and restart pointers have been
retained. `../ERT/ERT.namelist` is a separate eight-column SOAP setup and has not been
converted to SP; its climate path also points to `../../inputs/ERT/`.
