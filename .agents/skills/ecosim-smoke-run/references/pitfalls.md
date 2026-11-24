# Pitfalls behind the smoke-run procedure

Each item below caused a wrong result or a wasted run during development. They
are the reason the scripts look the way they do.

## 1. A magnitude filter silently hides every NaN

The output scan must exclude EcoSIM's numeric `_FillValue` (~1e36) so padding in
inactive PFT slots and soil layers is not mistaken for data. The obvious way to
write that is fatal:

```python
keep = np.abs(a) < 1e30        # WRONG: abs(nan) < 1e30 is False
```

Any comparison against NaN is `False`, so this drops NaN along with the
sentinels and the scan reports every file clean. Correct form:

```python
keep = ~np.isfinite(a) | (np.abs(a) < 1e30)
```

Non-finite values are kept precisely because they are what the scan is hunting.
**Verify any change to this logic by injecting a NaN into a copy of a real
output file and confirming the checker fails.** The clean-file case passing is
not evidence the checker works.

## 2. netCDF4 auto-masking manufactures NaN

`Dataset(...)` masks `_FillValue` by default. Filling that masked array with
NaN (or letting `np.ma` propagate) makes nearly every file look broken — the
first scan during development reported 19 of 21 files bad, all spurious. Always
`nc.set_auto_mask(False)` and handle fill explicitly.

Check the fill value itself, too: if `_FillValue` is NaN, model-produced NaN is
indistinguishable from padding and the variable cannot be judged. In current
output all 221 float variables carry a numeric fill, so nothing is skipped — but
that is a property of today's history writer, not a guarantee.

## 3. Relative input paths pin the scratch directory depth

Namelists reference inputs as `../../inputs/<case>/` and `../../../input_data/`.
A case therefore only works at depth `<root>/examples/run_dir/<case>/`. Copying
cases to an arbitrary temp directory breaks every path.

The mirror reproduces that depth and symlinks the two input trees:

```
$WORK/eco/input_data      -> $ROOT/input_data
$WORK/eco/examples/inputs -> $ROOT/examples/inputs
$WORK/eco/examples/run_dir/<case>/   (namelists only)
```

This is also what keeps committed outputs and `rpointer.esim` untouched — worth
preserving, since `examples/run_dir/` holds a large volume of untracked results.

## 4. macOS has no `timeout`

`timeout 1800 ...` returns 127 on macOS, which looks exactly like a model
failure — the first full sweep reported all 17 cases failing for this reason.
Use `( ulimit -t SECS; ... )` instead; it is portable and caps CPU time.

## 5. Exit code 0 does not mean the run finished

The model can exit cleanly having done less than asked. Check that the log's
final `wrote out restart data at nstep = N` matches `DAYS * 24`. That identity
holds because every example namelist sets `delta_time=3600.`; if a case ever
uses a different timestep, this check needs to read `delta_time` per case.

## 6. `sed -i` is not portable

BSD `sed` requires `-i ''`, GNU `sed` rejects it. Write to a temp file and
`mv` instead.

## 7. Cases with local restart dependencies

`dryland2` sets `finidat='./dryland_maize.ecosim.r.2004-01-01-000000.nc'`, a
file living in the run directory rather than under `examples/inputs/`. The
staging step copies anything referenced by a `'./...'` path; a case added later
with the same pattern is handled automatically.

That restart is *not* tracked in the repository (`git ls-files
examples/run_dir/dryland/` lists only the two namelists) -- it only exists if
someone has run `dryland` out to 2004 locally. So `dryland2` is excluded from
the smoke run by policy: see "Excluded namelists" in `SKILL.md`. Whether it
passes or fails depends on the operator's untracked working tree, not on the
build, so it is not evidence either way.

## 8. Expected benign output

- `CHECK_VAR: variable <name> is not on initial dataset` — normal when a
  `finidat` restart predates a newly added state variable.
- `Do error checking on file` — a status line, not an error.
- High all-fill variable counts (79) in `bare_soil`, `climeConst`, and `lake`:
  these cases have no plants, so every PFT variable is fill.
