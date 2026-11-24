#!/usr/bin/env bash
# Same case-directory interface as ecosim-smoke-run; defaults to 10 days.
set -eu
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$SCRIPT_DIR/run_restart.py" "$@"
