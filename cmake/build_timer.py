#!/usr/bin/env python3
"""Time build commands without changing their arguments, output, or exit status.

The append-only TSV log is shared safely by parallel Make/Ninja jobs on Unix.
ECOSIM_BUILD_TIMING_LOG and ECOSIM_BUILD_RUN_ID group launcher records with the
configure/build/install phases recorded by build_EcoSIM.sh.
"""

import argparse
import csv
from datetime import datetime, timezone
import fcntl
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

FIELDS = ("started_utc", "run_id", "phase", "elapsed_seconds", "exit_code", "step")


def step_name(command, phase):
    if phase == "compile":
        output = command[command.index("-o") + 1] if "-o" in command and command.index("-o") + 1 < len(command) else None
        sources = [arg for arg in command if arg != output and arg.lower().endswith(
            (".f", ".f90", ".f95", ".f03", ".f08", ".for", ".c", ".cc", ".cpp", ".cxx"))]
        if sources:
            return sources[-1]
    if phase == "link":
        if "-o" in command and command.index("-o") + 1 < len(command):
            return command[command.index("-o") + 1]
        # Static library rules run the archiver and (on some platforms) ranlib.
        archives = [arg for arg in command if arg.endswith((".a", ".lib"))]
        if archives:
            return Path(command[0]).name + " " + archives[0]
    return " ".join(shlex.quote(arg) for arg in command)


def append_record(path, record):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a+", newline="", encoding="utf-8") as stream:
        # Lock only the append, never the timed command. Compilation stays parallel.
        fcntl.flock(stream, fcntl.LOCK_EX)
        stream.seek(0, os.SEEK_END)
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t")
        if stream.tell() == 0:
            writer.writeheader()
        writer.writerow(record)
        stream.flush()
        fcntl.flock(stream, fcntl.LOCK_UN)


def run_command(args):
    command = args.command
    if command[:1] == ["--"]:
        command = command[1:]
    if not command:
        raise ValueError("a command is required after --")
    log = os.environ.get("ECOSIM_BUILD_TIMING_LOG") or args.log
    if not log:
        raise ValueError("--log or ECOSIM_BUILD_TIMING_LOG is required")
    started = datetime.now(timezone.utc).isoformat(timespec="milliseconds")
    start = time.monotonic()
    try:
        status = subprocess.call(command)
        if status < 0:
            status = 128 - status
    except OSError as error:
        print("build timer: " + str(error), file=sys.stderr)
        status = 127
    except KeyboardInterrupt:
        status = 130
    elapsed = time.monotonic() - start
    append_record(log, dict(zip(FIELDS, (
        started, os.environ.get("ECOSIM_BUILD_RUN_ID", "direct"), args.phase,
        f"{elapsed:.6f}", status, step_name(command, args.phase)))))
    return status


def report(args):
    with open(args.log, newline="", encoding="utf-8") as stream:
        fcntl.flock(stream, fcntl.LOCK_SH)
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if args.run_id:
        rows = [row for row in rows if row["run_id"] == args.run_id]
    print("EcoSIM build timing summary")
    print("Log: " + str(Path(args.log).resolve()))
    if args.run_id:
        print("Run: " + args.run_id)
    print("Durations are wall seconds. Parallel steps overlap; their sum is not build wall time.")
    print(f"Recorded commands: {len(rows)}; failed: {sum(int(r['exit_code']) != 0 for r in rows)}")
    print("\nBuild phases:")
    print(f"{'SECONDS':>10}  {'EXIT':>4}  PHASE")
    for row in sorted(rows, key=lambda r: r["started_utc"]):
        if row["phase"] not in ("compile", "link"):
            print(f"{float(row['elapsed_seconds']):10.3f}  {row['exit_code']:>4}  {row['phase']}")
    print("\nCompile/link steps, slowest first:")
    print(f"{'SECONDS':>10}  {'EXIT':>4}  {'PHASE':<8}  STEP")
    for row in sorted(rows, key=lambda r: float(r["elapsed_seconds"]), reverse=True):
        if row["phase"] in ("compile", "link"):
            print(f"{float(row['elapsed_seconds']):10.3f}  {row['exit_code']:>4}  {row['phase']:<8}  {row['step']}")
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="action")
    run = sub.add_parser("run", help="record one command's elapsed time and exit code")
    run.add_argument("--log")
    run.add_argument("--phase", required=True)
    run.add_argument("command", nargs=argparse.REMAINDER)
    summary = sub.add_parser("report", help="summarize an existing timing log")
    summary.add_argument("--log", required=True)
    summary.add_argument("--run-id")
    args = parser.parse_args()
    if not args.action:
        parser.error("an action is required")
    try:
        return run_command(args) if args.action == "run" else report(args)
    except (OSError, ValueError) as error:
        print("build timer: " + str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
