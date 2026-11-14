"""Compare the compiled VG/Mualem functions with 200-digit Decimal equations.

Uses only Python's standard library. The Fortran driver emits varied shapes,
pressure scales, and wet/dry contents; no model baseline is rewritten.
"""

import argparse
from decimal import Decimal, localcontext
from pathlib import Path
import subprocess
import sys


def exact(value):
    """Interpret each printed value as its actual binary double input/output."""
    return Decimal.from_float(float(value))


def check(executable):
    result = subprocess.run([str(executable), "dump"], check=True, capture_output=True, text=True)
    samples = [line.split()[1:] for line in result.stdout.splitlines() if line.startswith("SAMPLE ")]
    if len(samples) != 360:
        raise AssertionError(f"Expected 360 reference samples, received {len(samples)}")
    largest = Decimal(0)
    with localcontext() as context:
        context.prec = 200
        residual, porosity = exact(.05), exact(.45)
        dry = Decimal(-2500)
        tiny = exact(sys.float_info.min)
        tolerance = Decimal("3e-12")
        for row in samples:
            n, alpha, ell, theta, pressure, recovered, conductivity = map(exact, row)
            m = (n - 1) / n
            se = (theta - residual) / (porosity - residual)
            expected_pressure = max(dry, -((se ** (-1 / m) - 1) ** (1 / n)) / alpha)
            expected_water = residual + (porosity - residual) / (
                1 + (alpha * abs(max(pressure, dry))) ** n
            ) ** m
            expected_k = se ** ell * (1 - (1 - se ** (1 / m)) ** m) ** 2
            # Production intentionally rounds conductivities below normal double
            # range to zero; do not mask any larger discrepancies with an atol.
            if expected_k < tiny:
                expected_k = Decimal(0)
            for label, actual, expected in (
                ("pressure", pressure, expected_pressure),
                ("inverse", recovered, expected_water),
                ("conductivity", conductivity, expected_k),
            ):
                if expected == 0:
                    if actual != 0:
                        raise AssertionError(f"{label}: expected zero, received {actual}; {row}")
                    continue
                relative = abs((actual - expected) / expected)
                largest = max(largest, relative)
                if relative > tolerance:
                    raise AssertionError(f"{label}: relative difference {relative}; {row}")
    print(f"PASS: {len(samples)} VG samples, maximum relative difference {largest:.3E}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    args = parser.parse_args()
    check(args.executable.resolve(strict=True))


if __name__ == "__main__":
    main()
