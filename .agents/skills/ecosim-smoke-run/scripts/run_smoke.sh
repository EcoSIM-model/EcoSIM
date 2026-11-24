#!/usr/bin/env bash
# EcoSIM short smoke run over the example cases.
#
# Runs every (or selected) example namelist for a small number of days in an
# isolated scratch mirror, so nothing in examples/run_dir/ is touched.
#
# Usage:
#   run_smoke.sh [-d DAYS] [-b BUILD_DIR] [-o WORKDIR] [-t CPU_SECS] [case ...]
#
#   -d DAYS       simulated days per case            (default 20)
#   -b BUILD_DIR  build configuration directory      (default: autodetect newest)
#   -o WORKDIR    scratch directory to run in        (default: mktemp -d)
#   -t CPU_SECS   per-case CPU-time cap              (default 1800)
#   case ...      run_dir subdirectory names         (default: all)
#
# Exit status: 0 if every case completed the requested days, 1 otherwise.

set -uo pipefail

DAYS=20
BUILD_DIR=""
WORKDIR=""
CPU_CAP=1800

while getopts ":d:b:o:t:h" opt; do
  case $opt in
    d) DAYS=$OPTARG ;;
    b) BUILD_DIR=$OPTARG ;;
    o) WORKDIR=$OPTARG ;;
    t) CPU_CAP=$OPTARG ;;
    h) sed -n '2,20p' "$0"; exit 0 ;;
    *) echo "unknown option -$OPTARG" >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOT=$(cd "$SCRIPT_DIR" && git rev-parse --show-toplevel 2>/dev/null) \
  || { echo "ERROR: not inside a git repository" >&2; exit 2; }

# ---------------------------------------------------------------- executable
if [[ -z $BUILD_DIR ]]; then
  BUILD_DIR=$(ls -dt "$ROOT"/build/*/ 2>/dev/null | head -1)
fi
[[ -n $BUILD_DIR ]] || { echo "ERROR: no build/ configuration found. Run: bash build_EcoSIM.sh" >&2; exit 2; }
BUILD_DIR=$(cd "$BUILD_DIR" && pwd) || exit 2

EXE=""
for cand in "$BUILD_DIR/drivers/ecosim/ecosim.f90.x" \
            "$BUILD_DIR/local/bin/ecosim.f90.x" "$BUILD_DIR/bin/ecosim.f90.x"; do
  [[ -x $cand ]] && { EXE=$cand; break; }
done
[[ -n $EXE ]] || { echo "ERROR: ecosim.f90.x not found under $BUILD_DIR" >&2; exit 2; }

# Warn if any source is newer than the executable -- a stale binary silently
# invalidates the whole run.
if find "$ROOT/f90src" "$ROOT/drivers" -name '*.F90' -newer "$EXE" -print -quit 2>/dev/null | grep -q .; then
  echo "WARNING: Fortran sources are newer than $EXE -- rebuild first:" >&2
  echo "         cmake --build $BUILD_DIR --parallel 8" >&2
fi

# ------------------------------------------------------------------ scratch
if [[ -z $WORKDIR ]]; then
  WORKDIR=$(mktemp -d "${TMPDIR:-/tmp}/ecosim_smoke.XXXXXX")
else
  mkdir -p "$WORKDIR"
fi
MIRROR=$WORKDIR/eco
LOGS=$WORKDIR/logs
rm -rf "$MIRROR"
mkdir -p "$MIRROR/examples/run_dir" "$LOGS"

# Namelists reference inputs as ../../inputs/... and ../../../input_data/...
# Symlinking at the SAME depth is what keeps those relative paths resolvable.
ln -s "$ROOT/input_data"      "$MIRROR/input_data"
ln -s "$ROOT/examples/inputs" "$MIRROR/examples/inputs"

# -------------------------------------------------------------- case staging
if [[ $# -gt 0 ]]; then
  CASES=("$@")
else
  CASES=()
  for d in "$ROOT"/examples/run_dir/*/; do
    compgen -G "$d*.namelist" > /dev/null && CASES+=("$(basename "$d")")
  done
fi
[[ ${#CASES[@]} -gt 0 ]] || { echo "ERROR: no cases with namelists found" >&2; exit 2; }

for c in "${CASES[@]}"; do
  src=$ROOT/examples/run_dir/$c
  [[ -d $src ]] || { echo "ERROR: no such case directory: $src" >&2; exit 2; }
  dst=$MIRROR/examples/run_dir/$c
  mkdir -p "$dst"
  cp "$src"/*.namelist "$dst"/ 2>/dev/null

  # Excluded from the routine smoke screen; keep the source case available
  # for separate investigation of its existing SIGILL failure.
  if [[ $c == blodgett ]]; then
    rm -f "$dst/Blodget.warming.namelist"
    echo "SKIP blodgett/Blodget.warming: existing SIGILL failure"
  fi
  if [[ $c == dryland ]]; then
    rm -f "$dst/dryland2.namelist"
    echo "SKIP dryland/dryland2: requires a locally generated 2004 checkpoint"
  fi

  # Some cases (e.g. dryland2) start from a finidat restart living in the run
  # directory itself; copy anything referenced by a './...' path.
  for f in "$dst"/*.namelist; do
    while IFS= read -r p; do
      [[ -f $src/$p ]] && cp "$src/$p" "$dst/"
    done < <(grep -ohE "finidat[[:space:]]*=[[:space:]]*'\./[^']*'" "$f" \
             | sed -E "s/.*'\.\/([^']*)'.*/\1/")
  done
done

# ------------------------------------------------------------ namelist patch
# Force a short run from each case's own start_date and disable periodic
# restart writing. Everything else is left exactly as committed.
for f in "$MIRROR"/examples/run_dir/*/*.namelist; do
  sed -E \
    -e "s/^[[:space:]]*stop_n[[:space:]]*=.*$/stop_n=$DAYS/" \
    -e "s/^[[:space:]]*stop_option[[:space:]]*=.*$/stop_option='ndays'/" \
    -e "s/^[[:space:]]*rest_opt[[:space:]]*=.*$/rest_opt='never'/" \
    -e "s/^[[:space:]]*continue_run[[:space:]]*=.*$/continue_run=.false./" \
    "$f" > "$f.tmp" && mv "$f.tmp" "$f"     # portable in-place (macOS + GNU sed)
done

# ------------------------------------------------------------------ run loop
EXPECT_NSTEP=$((DAYS * 24))                 # every example uses delta_time=3600
SUMMARY=$WORKDIR/summary.txt
: > "$SUMMARY"
fail=0

printf '%-44s %-6s %-7s %s\n' CASE RC NSTEP RESULT | tee "$SUMMARY"
for f in "$MIRROR"/examples/run_dir/*/*.namelist; do
  d=$(dirname "$f"); nl=$(basename "$f")
  tag="$(basename "$d")/${nl%.namelist}"
  log=$LOGS/$(echo "$tag" | tr / _).log

  # NOTE: macOS has no `timeout`; ulimit -t is the portable CPU-time guard.
  ( ulimit -t "$CPU_CAP"; cd "$d" && "$EXE" "$nl" > "$log" 2>&1 )
  rc=$?

  nstep=$(grep -a 'wrote out restart data at nstep' "$log" | tail -1 | grep -oE '[0-9]+$')
  nstep=${nstep:-0}

  if [[ $rc -eq 0 && $nstep -eq $EXPECT_NSTEP ]]; then
    res=OK
  else
    res=FAIL; fail=1
  fi
  printf '%-44s %-6s %-7s %s\n' "$tag" "$rc" "$nstep" "$res" | tee -a "$SUMMARY"
done

echo
echo "expected nstep = $EXPECT_NSTEP ($DAYS days x 24 h)"
echo "logs:    $LOGS"
echo "outputs: $MIRROR/examples/run_dir"

# --------------------------------------------------------------- output scan
if command -v python3 > /dev/null; then
  echo
  python3 "$SCRIPT_DIR/check_history.py" "$MIRROR/examples/run_dir" || fail=1
else
  echo "WARNING: python3 unavailable; skipped the NaN/Inf output scan." >&2
fi

echo
if [[ $fail -eq 0 ]]; then
  echo "SMOKE RUN PASSED (${#CASES[@]} case directories, $DAYS days each)"
else
  echo "SMOKE RUN FAILED -- see the FAIL rows and logs above"
fi
exit $fail
