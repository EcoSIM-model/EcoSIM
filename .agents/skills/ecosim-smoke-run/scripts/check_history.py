#!/usr/bin/env python3
"""Scan EcoSIM history output for genuine NaN/Inf and implausible values.

Usage:  check_history.py RUN_DIR [--tmin -90] [--tmax 80]

RUN_DIR is a directory tree containing *.ecosim.h*.nc files.
Exit status is 1 if anything suspicious is found, else 0.

The important subtlety: EcoSIM writes a numeric _FillValue (~1e36) into
inactive PFT slots and inactive soil layers.  netCDF4 auto-masking turns those
into NaN, so a naive `isnan` scan reports almost every file as broken.  This
script disables auto-masking and excludes fill explicitly, so only NaN/Inf that
the model actually produced is reported.
"""

import argparse
import glob
import os
import sys

import numpy as np
from netCDF4 import Dataset

FILL_CUTOFF = 1e30          # values at/above this are sentinels, not data


def real_values(var):
    """Return the finite-candidate data of `var` with fill sentinels removed.

    Returns (values, skipped_reason).  `skipped_reason` is non-None when the
    variable cannot be judged -- notably when _FillValue is itself NaN, which
    makes model-produced NaN indistinguishable from padding.
    """
    a = np.asarray(var[:], dtype=float)
    fv = getattr(var, "_FillValue", None)
    if fv is None:
        fv = getattr(var, "missing_value", None)

    if fv is not None and np.isnan(float(fv)):
        return None, "_FillValue is NaN"

    # Keep every non-finite value -- those are exactly what we are hunting.
    # Only FINITE values are subject to the sentinel cutoff; writing this as a
    # bare `abs(a) < CUTOFF` silently discards NaN, because any comparison
    # against NaN is False, and the scan then reports everything as clean.
    keep = ~np.isfinite(a) | (np.abs(a) < FILL_CUTOFF)
    if fv is not None:
        keep &= a != float(fv)
    return a[keep], None


def scan_file(path, tmin, tmax):
    issues, skipped, allfill = [], [], []

    with Dataset(path) as nc:
        nc.set_auto_mask(False)                 # critical -- see module docstring
        for name, var in nc.variables.items():
            if not np.issubdtype(var.dtype, np.floating):
                continue

            vals, reason = real_values(var)
            if reason:
                skipped.append(f"{name} ({reason})")
                continue

            if vals.size == 0:
                allfill.append(name)            # normal for PFT vars in bare cases
                continue

            n_nan = int(np.isnan(vals).sum())
            n_inf = int(np.isinf(vals).sum())
            if n_nan or n_inf:
                issues.append(f"{name}: NaN={n_nan} Inf={n_inf}")
                continue

            # Soil/air temperature is the cheapest physical tripwire we have.
            if name in ("TEMP_vr", "TKS_vr") and vals.size:
                lo, hi = float(vals.min()), float(vals.max())
                if lo < tmin or hi > tmax:
                    issues.append(f"{name}: range [{lo:.4g}, {hi:.4g}] degC "
                                  f"outside [{tmin}, {tmax}]")

    return issues, skipped, allfill


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("--tmin", type=float, default=-90.0,
                    help="min plausible soil temperature, degC (default -90)")
    ap.add_argument("--tmax", type=float, default=80.0,
                    help="max plausible soil temperature, degC (default 80)")
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.run_dir, "**", "*.ecosim.h*.nc"),
                             recursive=True))
    if not files:
        print(f"OUTPUT SCAN: no history files under {args.run_dir}", file=sys.stderr)
        return 1

    bad = 0
    for path in files:
        rel = os.path.relpath(path, args.run_dir)
        try:
            issues, skipped, allfill = scan_file(path, args.tmin, args.tmax)
        except Exception as exc:                # unreadable file is itself a failure
            print(f"  FAIL {rel}: cannot read ({exc})")
            bad += 1
            continue

        if issues:
            bad += 1
            print(f"  FAIL {rel}")
            for line in issues[:12]:
                print(f"         {line}")
            if len(issues) > 12:
                print(f"         ... +{len(issues) - 12} more")
        else:
            note = f"  ({len(allfill)} all-fill vars)" if allfill else ""
            print(f"  ok   {rel}{note}")

        if skipped:
            print(f"         NOT CHECKED: {len(skipped)} var(s) with NaN _FillValue")
            for line in skipped[:5]:
                print(f"           {line}")

    print(f"\nOUTPUT SCAN: {len(files) - bad}/{len(files)} history files clean")
    if bad:
        print("Note: 'all-fill vars' alone is not a fault -- plantless cases "
              "(bare_soil, lake, climeConst) legitimately fill every PFT variable.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
