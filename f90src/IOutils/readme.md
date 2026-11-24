# ioutils

|File                 | Description                                                                         |
|---------------------|-------------------------------------------------------------------------------------|
|BLOCKDATA001.f       |Set up data for doy calculation                                                      |
|foutp.f              |Set up label for plant variable output                                               |
|fouts.f              |Set up label for soil variable output                                                |
|outpd.f              |Daily output for plant C, N, P, water and heat                                       |
|outph.f              |Hourly output for plant C, N, P, water and heat                                      |
|outsd.f              |Daily output for soil C, N, P, water and heat                                        |
|outsh.f              |Hourly output for soil C, N, P, water and heat                                       |
|readimod.f              |Read soil and topographic inputs                                                     |
|readqmod.f              |Read pft and plant management inputs                                                 |
|readsmod.f              |Read soil and soil management inputs                                                 |
|routp.f              |Read pft checkpoint files                                                            |
|routq.f              |Read plant and management checkpoint files                                           |
|routs.f              |Read soil checkpoint files                                                           |
|split.f              |Splits output for C, N, P, water and heat to grids                                   |
|splitc.f             |Caller of splits                                                                     |
|woutp.f              |Write checkpoint file for plant state variables                                      |
|woutq.f              |Write pft species names for checkpoint                                               |
|wouts.f              |Write soil variables for checkpoint                                                  |
|splitp.c             |Split plant data file into separate files determined by grid cell indices and species|
|splits.c             |Split soil data file into separate files determined by grid cell indices             |
|blk10.h              |Soil water related variables                                                         |
|blk11a.h             |Soil heat and water related variables                                                |
|blk11b.h             |Soil gas related variables                                                           |
|blk12a.h             |Soil gas related variables                                                           |
|blk12b.h             |Nutrient uptake related variables                                                    |
|blk13a.h             |soil C, N, and P related variables                                                   |
|blk13b.h             |Soil gas and C, N, and P related variables                                           |
|blk13c.h             |Nutrient fertilizer related soil variables                                           |
|blk13d.h             |DOM variables                                                                        |
|blk14.h              |C, N and P balance related variables                                                 |
|blk15a.h             |Water and heat flux related variables                                                |
|blk15b.h             |Runoff related lateral gas exchanges                                                 |
|blk16.h              |Uptake related variables                                                             |
|blk17.h              |Doy related variables                                                                |
|blk18a.h             |Plant-soil exchange related variables                                                |
|blk18b.h             |Gas and nutrients related rate variables                                             |
|blk19a.h             |Salt model related variables                                                         |
|blk19b.h             |Salt model related variables                                                         |
|blk19c.h             |Salt model related variables                                                         |
|blk19d.h             |Salt model related variables                                                         |
|blk1cp.h             |Plant C, N, and P variables                                                          |
|blk1cr.h             |Plant root related variables                                                         |
|blk1g.h              |Photosynthesis related variables                                                     |
|blk1n.h              |Plant N variables                                                                    |
|blk1p.h              |Plant P varaibles                                                                    |
|blk1s.h              |Plant nutrient recycling related variables                                           |
|blk1u.h              |Canopy water and heat related variables                                              |
|blk20a.h             |Runoff chemistry related variables                                                   |
|blk20b.h             |Runoff chemistry related variables                                                   |
|blk20c.h             |Runoff chemistry related variables                                                   |
|blk20d.h             |Runoff C, N and P fluxes related variables                                           |
|blk20e.h             |Runoff C, N and P fluxes related variables                                           |
|blk20f.h             |Erosion related variables                                                            |
|blk21a.h             |Precipitation model related variables                                                |
|blk21b.h             |Fertilzer related variables                                                          |
|blk22a.h             |Irrigation related variables                                                         |
|blk22b.h             |Irrigation related chemical variables                                                |
|blk22c.h             |Deposition related variables                                                         |
|blk2a.h              |Climate related variables                                                            |
|blk2b.h              |Wet deposition related variables                                                     |
|blk2c.h              |Irrigation related tracer variables                                                  |
|blk3.h               |Plant phenology related variables                                                    |
|blk5.h               |Radiation related plant variables                                                    |
|blk6.h               |Radiation related variables                                                          |
|blk8a.h              |Soil water related variables                                                         |
|blk8b.h              |Soil phyiscal variables (include seeding depth)                                      |
|blk9a.h              |Canopy related variables                                                             |
|blk9b.h              |Plant biochemistry related variables                                                 |
|blk9c.h              |Harvest related variables                                                            |
|blkc.h               |Plant biochemistry related variables                                                 |
|blktest.h            |Soil warming related variables                                                       |
|filec.h              |Char variables for file management                                                   |
|files.h              |File management files                                                                |
|parameters.h         |Basic size parameters                                                                |
