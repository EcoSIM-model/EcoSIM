# EcoSIM repository guidance

## Scope and purpose

This repository contains the EcoSIM biogeochemical model and a biocrust simulation
workspace. Use these instructions throughout the repository; follow more specific
`AGENTS.md` files in subdirectories where applicable. Explicit user instructions
take precedence over repository guidance.

## Repository map

- `f90src/`: Fortran model source. `Main/` and `Ecosim_mods/` coordinate model
  execution; `Modelforc/` handles forcing; `IOutils/` handles inputs, history,
  and restarts.
- `f90src/Plant_bgc/`, `Microbial_bgc/`, `Geochem/`, `HydroTherm/`, and
  `Transport/`: biological, chemical, hydrological, thermal, and transport processes.
- `f90src/Ecosim_datatype/`, `APIData/`, `Modelconfig/`, and `Modelpars/`:
  state/flux definitions, API data, identifiers, and parameters.
- `f90src/Balances/` and `ModelDiags/`: accounting and diagnostics.
- `drivers/`, `f90src/APIs/`, and `f90src/ATSUtils/`: drivers and coupling interfaces.
- `input_data/`: shared model inputs and PFT parameter NetCDF/CDL files.
- `examples/inputs/`: case-specific inputs. `examples/run_dir/`: case namelists
  and simulation working directories, including `biocrust/`.
- `ecosim-agent/`: a separate Git submodule for preprocessing and analysis.
  Read its `AGENTS.md` before working there. Its canonical skills are under
  `ecosim-agent/.agents/skills/`; generated preprocessing outputs belong under
  `ecosim-agent/result/` unless the user specifies another destination.
- `tests/` and `regression-tests/`: validation resources and regression tooling.
- `build/` and `local/`: generated build/install products.
- `3rd-partylibs/`: third-party dependency submodules.

## Knowledge documents and user preferences

- When asked to commit a fix or change, commit locally only. Do not push or
  otherwise update a remote unless the user explicitly says "commit to remote".
- Save generated EcoSIM process explanations, code walkthroughs, and similar
  knowledge documents under `ecosim_knowledge/` by default, unless the user
  specifies another destination. This is a persistent user preference for
  future requests.
- Write HTML files to `ecosim_knowledge/` only when the user explicitly asks
  for an HTML format.
- Inside `ecosim_knowledge/`, stage and commit only HTML files and the content
  they reference (images, CSS, data files they load). Leave everything else
  there untracked, including `.DS_Store` and intermediate or scratch files.
- Keep relative source links valid from that directory when creating or moving
  a document.

## Working with model code

- When a code review identifies an issue, trace the same underlying problem
  through related equations, variables, callers, and downstream consumers
  before reporting it. Group confirmed instances sharing a root cause or
  requiring a coordinated fix into one finding, listing the affected sites
  and the full scope of the recommended correction. Keep independent issues
  separate and distinguish confirmed effects from unverified concerns.
- Inspect `git status` before editing. Preserve existing user changes, untracked
  inputs, simulation outputs, and restart files; do not clean or reset them.
- Trace callers and CMake source lists to establish the active implementation.
  Similar routines occur in initialization, process, and coupling modules; do not
  assume one is obsolete or change every copy without checking its role.
- When multiple correct implementation options are available, choose the simpler
  one. Prioritize clarity in code, equations, naming, and comments.
- Follow nearby Fortran conventions, precision kinds, array ordering, and module
  interfaces. Update allocation, initialization, API mapping, history, and restart
  handling as needed when adding persistent state or flux variables.
- Every Fortran module must include a header comment describing its purpose and
  responsibilities.
- Every subroutine and function must include a description comment explaining
  its purpose.
- Each module or procedure description must start with a standalone comment line
  `! Description:`, followed by the actual description on subsequent comment
  lines. For example:
  ```fortran
  ! Description:
  ! Compute soil matric potential from liquid water content.
  ```
- For new subroutines, declare `intent(in)`, `intent(out)`, or `intent(inout)`
  for each data argument and annotate its declaration with a comment explaining
  its meaning, units, and array dimensions where applicable.
- New Fortran procedures should declare
  `character(len=*), parameter :: subname='procedure_name'`, using the actual
  procedure name as the string value. Include `call PrintInfo('beg '//subname)`
  at entry and `call PrintInfo('end '//subname)` at normal exit. Separate
  procedures with the following comment line:
  `!------------------------------------------------------------------------------------------`
- When a diagnostic `write` statement is needed before calling `endrun`, write
  to `iulog` rather than the default output unit `*` (use `write(iulog,*)`).
- Use `associate` blocks only when creating a new subroutine or modifying a
  subroutine that already contains an `associate` block. In those subroutines,
  use aliases for derived-type components in executable statements, adding
  aliases to the existing block or introducing a block in the new subroutine
  as needed. When modifying an existing subroutine without an `associate`
  block, keep direct derived-type component references such as
  `plt_morph%Root1stDepz_raxes(NR,NZ)`; do not introduce an `associate` block.
  Exception: when adding history-output variables,
  avoid `associate` aliases and use direct history components such as
  `this%h1D_variable_col` in allocation, registration, and update code, following
  the surrounding history-output conventions.
- If a procedure uses `associate`, place `associate (...)` immediately after
  `call PrintInfo('beg '//subname)`, except that safety checks for an early stop
  may appear between this call and `associate`. Place `end associate` immediately
  before `call PrintInfo('end '//subname)`.
- Annotate every alias in an `associate` block with a comment describing the
  variable's meaning, units, array dimensions where applicable, and whether it
  is an input, output, or input/output within the procedure. Apply this rule to
  new blocks and to aliases added or modified in existing blocks.
- Keep edits focused. Do not hand-edit generated build files or modify dependency
  submodules unless that is part of the requested work.
- Treat `ecosim-agent` changes as submodule changes, separate from root-repository
  changes. Do not automatically commit, update submodule pointers, or publish work.
- Preserve staged changes inside `ecosim-agent`. Do not edit files with staged
  changes or alter their index entries (including staging, unstaging, or committing)
  unless the user explicitly requests updates to those staged contents.

## Scientific and input conventions

- When asked to evaluate a plant-trait `.desc` file (including
  `*.plant_trait.desc.*`), read `ecosim-agent/AGENTS.md` and use the relevant
  skills under `ecosim-agent/.agents/skills/`, starting with
  `ecosim-plant-trait-sanity-check/SKILL.md`. Treat the `.desc` file as read-only.
- Track units explicitly: distinguish carbon mass from dry mass, per-individual
  quantities from per-area quantities, and per-area values from grid-cell totals.
  Check the actual equations when comments and metadata disagree.
- Preserve C/N/P, water, and energy accounting for affected processes. Explain
  newly introduced sources, sinks, transfers, or numerical clipping.
- Read enum values from `f90src/Modelconfig/ElmIDMod.F90`; do not rely solely on
  NetCDF flag descriptions. Validate PFT codes against the parameter file selected
  by the case, and check whether grid Koppen settings replace the PFT suffix.
- EcoSIM treats lichen and moss through nonvascular plant pathways. Interpret
  seed/planting terminology according to the implemented pools and equations;
  do not equate establishment reserves with an initialized mature community.
- Planting strings use `DDMMYYYY`, population density in units per square meter,
  and depth in meters. Honor user-provided values without silently substituting
  biological defaults. For a single grid cell, use
  `NH1 = NV1 = NH2 = NV2 = 1`; `NZ` is the active PFT count.
- Use the relevant Python-tools skill for management or other preprocessing.
  Validate NetCDF dimensions, strings, units, and active entries by reading the
  generated file back. Keep editable source inputs alongside generated artifacts
  when practical. Keep paired CDL/NetCDF representations consistent when changing
  an input dataset, and identify which file the model actually reads.

## Site run requests

- For requests to create, set up, or run an EcoSIM case for a site (for
  example "create an ecosim run for site US-UMB"), follow the
  `ecosim-agent` run-folder workflow
  (`ecosim-agent/.agents/skills/ameriflux-namelist-generator/references/runfolder_workflow.md`)
  and track progress with `ecosim-agent/.agents/skills/ecosim-site-run-tracker/`.
  Initialize the tracker with the request text first, record the executable,
  test run, full run, and evaluation results with its `set` command, and
  report stage status from `ecosim-agent/result/<SITE_ID>/PROGRESS.md`.

## Build and validation

- When the agent runs inside or is linked to VSCode, build and test in the
  environment of VSCode's integrated terminal rather than the agent's default
  shell, which may resolve a different compiler toolchain. Before building,
  determine that environment from the VSCode settings: workspace
  `.vscode/settings.json` overrides the user `settings.json` (location depends
  on the OS). Read `terminal.integrated.defaultProfile.<os>` and the matching
  `terminal.integrated.profiles.<os>` entry for the shell and its arguments, and
  `terminal.integrated.env.<os>` for environment variables (expand references
  such as `${env:PATH}`). Run build and test commands through that shell with
  those variables exported. Do not record the resolved values in the repository;
  re-read the settings for each session.
- Reuse a build directory whose `CMakeCache.txt` compilers match that
  environment; do not mix toolchains in one build directory. If none matches,
  configure a new one. Report the toolchain used along with validation results.
  If VSCode is not involved or no terminal settings exist, use the default shell
  and say so.
- Run build commands from the repository root. Inspect supported options with
  `bash build_EcoSIM.sh --help`; the standard build is
  `bash build_EcoSIM.sh`. Debug/regression configuration is available through
  `bash build_EcoSIM.sh --debug --regression_test`.
- For an existing configured build, use
  `cmake --build build/<configuration> --parallel <jobs>`, selecting an actual
  configuration directory and checking its compiler/cache first. Locate the
  resulting executable rather than assuming a fixed path from older documentation.
- `build_EcoSIM.sh --clean` deletes both `build/` and `local/`; do not use it as
  routine validation.
- Select validation appropriate to the change. Fortran changes generally need a
  build and a focused test or short simulation exercising the modified process.
  Input-only changes need schema/value checks; documentation-only changes do not
  require compiling the model.
- To screen the example cases after a Fortran change, use the
  `ecosim-smoke-run` skill (`.agents/skills/ecosim-smoke-run/`, also exposed at
  `.claude/skills/` and `.codex/skills/`). It runs each example for a few
  simulated days in a scratch mirror and checks exit status, steps completed,
  and output for NaN/Inf. A pass means the model ran and produced finite
  numbers; it is not evidence of scientific correctness.
- To test restart reproducibility, use the `ecosim-restart-test` skill
  (`.agents/skills/ecosim-restart-test/`, also exposed at `.claude/skills/`
  and `.codex/skills/`). It compares 10-day uninterrupted example runs with
  five-day runs followed by five-day continuations in scratch directories,
  requiring exact checkpoint/history equality and finite valid values.
  Report failed cases; do not relax tolerances or update baselines to pass.
- To run from an existing checkpoint of a long simulation (branch with new
  history fields or changed inputs, or continue the parent case), use the
  `ecosim-restart-run` skill (`.agents/skills/ecosim-restart-run/`, also
  exposed at `.claude/skills/` and `.codex/skills/`).
- Inspect `regression-tests/Makefile` before running its targets. It provides
  `mtest` and `rtest`; verify the executable path and required data first. Do not
  update regression baselines merely to make a changed result pass.
- Run simulations in a scratch/copy of the case when existing outputs or restart
  pointers could be overwritten. Inspect the namelist, input paths, run duration,
  and restart mode before launching.
- Report what changed, validation actually performed, and any unresolved limits.
  Distinguish a successful build from demonstrated scientific correctness.
