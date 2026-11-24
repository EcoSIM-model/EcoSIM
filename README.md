# EcoSIM

A biogeochemical modeling library spins off the ecosys model.

# Pull requests (for contributors & collaborators)

Before submitting your PR:
- make sure smoke-run-test and restart-run-test are passed (via coding agent).
- Search for existing PRs to prevent duplicating efforts
- Rebase your change with main branch & Test your changes:
-Create separate PRs for each feature or fix:
 - Avoid combining unrelated changes in a single PR
 - For intricate features, consider opening a feature request first to discuss and align expectations

 - Consider allowing write access to your branch for faster reviews, as reviewers can push commits directly
- If you are a new contributor
  - Limit your open PRs to 1
  - Do not submit trivial fixes (e.g. typos, formatting changes)
- submit the PR to next, avoid main

## Download

git clone --recursive git@github.com:jinyun1tang/EcoSIM.git

## Building

Use the script in the root directory to build the code:

sh build_EcoSIM.sh

or run it as a bash script:

./build_EcoSIM.sh

This should do everything including building and linking all the required tpls, and configuring, building and linking EcoSIM itself.

The executable will be found in ./build/bin/ecosim.f90.x

This has been tested with gcc on multiple systems, but not intel compilers. See the .github/workflow/ecosim-ci.yml file for examples. You will need cmake to use the build_EcoSIM.sh script.  

## Additional options

The build_EcoSIM.sh script lets you set some opttional parameters:

debug=0
mpi=0
shared=0
verbose=0
sanitize=0
regression_test=0

precision="double"
prefix=""
systype=""

#Leave empty to just use the environment variable compiler
CC=""
CXX=""
FC=""

You can also set options via the command line, for example:

./build_EcoSIM.sh CC=/path/to/cc CXX=/path/to/cxx FC=/path/to/fortran --debug --regression_test

Run

./build_EcoSIM.sh --help 

for a full list of optional arguments

## using docker on mac 
1. make sure docker is installed on your machine, please see https://www.docker.com, and ask for help from whatever chatbot you have access to.

2. to turn on rocky-linux docker container, run

./run_docker.sh --fresh --rocky

where --fresh means the first time run docker, otherwise, use 

./run_docker.sh --rocky

One can also try ubuntu container. However, it seems buggy on my apple machine.

## Repository layout

- `f90src/` — Fortran source of the EcoSIM model (process modules, drivers, I/O, data types)
- `calibration_tool/` — agent-based calibration framework for EcoSIM; wraps parameter perturbation, simulation, and target extraction into a reset/step environment. See `calibration_tool/usage.html`
- `examples/` — example model inputs and ready-to-run case directories (e.g. `examples/run_dir/biocrust`), each with namelist and case-specific forcing
- `input_data/` — shared model input datasets, including PFT parameter NetCDF files, atmospheric GHG forcing, and climate inputs referenced by multiple cases
- `ecosim-agent/` — Python preprocessing and analysis bridge (git submodule): climate/soil/management data preparation, PFT parameter editing, output target extraction; canonical skills live under `ecosim-agent/.agents/skills/`
- `ecosim_knowledge/` — HTML notes documenting model science: process formulations such as canopy water potential, transpiration feedbacks, leaf respiration, and microbial complex kinetics
- `tests/` and `regression-tests/` — unit and regression test resources




