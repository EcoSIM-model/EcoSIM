"""Compare 10 continuous EcoSIM days with a 5-day run plus 5-day continuation.

Select example directories as with ecosim-smoke-run. Requires numpy/netCDF4.
"""
import argparse
import json
from pathlib import Path
import resource
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[4]
EXCLUSIONS = {
    "blodgett/Blodget.warming": "existing SIGILL failure; investigate separately",
    "dryland/dryland2": "requires a locally generated 2004 checkpoint",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-d", "--days", type=int, default=10)
    parser.add_argument("-b", "--build-dir", type=Path, help="configured build (default: newest executable)")
    parser.add_argument("-o", "--workdir", type=Path, help="new/empty directory (default: temporary directory)")
    parser.add_argument("-t", "--cpu-secs", type=int, default=1800, help="CPU limit per model process")
    parser.add_argument("-s", "--split-hours", type=int, help="restart hour (default: half the run)")
    parser.add_argument("cases", nargs="*", help="example case directories (default: all)")
    args = parser.parse_args()
    split = args.split_hours if args.split_hours is not None else args.days * 12
    if not 0 < split < args.days * 24 or args.cpu_secs <= 0:
        parser.error("require positive CPU limit and restart split inside the run")
    try:
        import numpy
        import netCDF4
    except ImportError as error:
        parser.error(str(error))
    builds = [args.build_dir.resolve()] if args.build_dir else list((ROOT / "build").glob("*/"))
    executables = []
    for build in builds:
        # Prefer the build product because installed copies can be stale.
        found = next((build / suffix for suffix in (
            "drivers/ecosim/ecosim.f90.x", "local/bin/ecosim.f90.x", "bin/ecosim.f90.x")
            if (build / suffix).is_file()), None)
        if found:
            executables.append(found)
    if not executables:
        parser.error("no executable found; build EcoSIM first or supply -b BUILD_DIR")
    executable = max(executables, key=lambda path: path.stat().st_mtime)
    newer = next((p for tree in (ROOT / "f90src", ROOT / "drivers")
                  for p in tree.rglob("*.F90") if p.stat().st_mtime > executable.stat().st_mtime), None)
    if newer:
        parser.error(f"executable is older than {newer}; rebuild first")
    examples = ROOT / "examples/run_dir"
    cases = args.cases or sorted(p.name for p in examples.iterdir()
                                 if p.is_dir() and any(p.glob("*.namelist")))
    selected, skipped = [], []
    for case in dict.fromkeys(cases):
        if Path(case).name != case or not (examples / case).is_dir():
            parser.error(f"not an example case directory: {case}")
        namelists = sorted((examples / case).glob("*.namelist"))
        if not namelists:
            parser.error(f"no namelists in {case}")
        for namelist in namelists:
            tag = f"{case}/{namelist.stem}"
            if tag in EXCLUSIONS:
                skipped.append({"case": tag, "reason": EXCLUSIONS[tag]})
            else:
                selected.append((tag, namelist))
    if not selected:
        parser.error("no eligible namelists selected")
    work = args.workdir.resolve() if args.workdir else Path(tempfile.mkdtemp(prefix="ecosim_restart_"))
    if work.exists() and (not work.is_dir() or any(work.iterdir())):
        parser.error(f"output directory must be new or empty: {work}")
    work.mkdir(parents=True, exist_ok=True)
    (work / "logs").mkdir()
    report = {"executable": str(executable), "days": args.days, "split_hours": split,
              "cpu_secs": args.cpu_secs, "skipped": skipped, "results": []}
    lines = [f"{'CASE':44} {'RC':>4} {'DIFFS':>7} RESULT"]
    print(f"Executable: {executable}\nOutputs: {work}", flush=True)
    for entry in skipped:
        print(f"SKIP {entry['case']}: {entry['reason']}", flush=True)
    print(lines[0], flush=True)

    def cpu_limit():
        _, hard = resource.getrlimit(resource.RLIMIT_CPU)
        limit = args.cpu_secs if hard == resource.RLIM_INFINITY else min(args.cpu_secs, hard)
        resource.setrlimit(resource.RLIMIT_CPU, (limit, hard))

    for tag, namelist in selected:
        folder = work / "cases" / tag
        log = work / "logs" / (tag.replace("/", "_") + ".log")
        command = [sys.executable, str(ROOT / "tests/check_restart_exactness.py"),
                   "--executable", str(executable), "--namelist", str(namelist),
                   "--workdir", str(folder), "--days", str(args.days), "--split-hours", str(split)]
        with log.open("w") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, preexec_fn=cpu_limit)
        comparison = folder / "comparison.json"
        differences = None
        if comparison.exists():
            details = json.loads(comparison.read_text())
            differences = sum(len(items) for items in details["comparisons"].values())
        passed = result.returncode == 0 and differences == 0
        entry = {"case": tag, "status": "PASS" if passed else "FAIL", "returncode": result.returncode,
                 "differences": differences, "comparison": str(comparison), "log": str(log)}
        report["results"].append(entry)
        count = str(differences) if differences is not None else "-"
        line = f"{tag:44} {result.returncode:4} {count:>7} {entry['status']}"
        lines.append(line)
        print(line, flush=True)
        (work / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
        (work / "summary.txt").write_text("\n".join(lines) + "\n")
    failed = sum(entry["status"] == "FAIL" for entry in report["results"])
    print(f"\n{len(selected) - failed}/{len(selected)} cases passed; {work / 'summary.json'}", flush=True)
    return int(failed != 0)


if __name__ == "__main__":
    raise SystemExit(main())
