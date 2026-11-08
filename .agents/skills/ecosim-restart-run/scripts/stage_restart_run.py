#!/usr/bin/env python3
"""Stage (and optionally launch) an EcoSIM run that starts from a checkpoint.

Modes
  branch    New case from a checkpoint: fresh history files, history fields may
            change, inputs may differ.  Sets continue_run=.false. and finidat.
  continue  Resume the parent case: appends to the parent's history files;
            history fields cannot change.  Sets continue_run=.true. and writes
            rpointer.esim.

The parent case directory is never written to.  Input paths are rewritten to
absolute paths; files named with --copy-input are copied so they can be edited.
"""
import argparse
import calendar
from datetime import datetime, timedelta
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

INPUT_KEYS = ("pft_file_in", "grid_file_in", "pft_mgmt_in", "clm_hour_file_in",
              "clm_day_file_in", "soil_mgmt_in", "atm_ghg_in", "clm_factor_in",
              "sp_phenol_in")
TIME_KEYS = ("stop_n", "stop_option", "rest_opt", "rest_frq", "delta_time")
RESTART_NAME = re.compile(r"^(?P<case>.+)\.ecosim\.r\.(?P<y>\d{4})-(?P<m>\d{2})-(?P<d>\d{2})-(?P<s>\d+)\.nc$")


def active(line):
    return not line.lstrip().startswith("!")


def strip_inline_comment(line):
    """Remove namelist comments while retaining ! and escaped quotes in strings."""
    quote = None
    i = 0
    while i < len(line):
        char = line[i]
        if quote:
            if char == quote:
                if i + 1 < len(line) and line[i + 1] == quote:
                    i += 1
                else:
                    quote = None
        elif char in "'\"":
            quote = char
        elif char == "!":
            return line[:i]
        i += 1
    return line


def get_value(text, key):
    """Last active assignment wins, as in Fortran namelist reads."""
    val = None
    for line in text.splitlines():
        line = strip_inline_comment(line)
        m = re.match(rf"^\s*{key}\s*=\s*(.*?)\s*$", line, re.IGNORECASE)
        if m and active(line):
            val = m.group(1)
    return val


def set_value(text, key, value, group):
    """Replace every active assignment of key, or insert it at the top of group."""
    lines, hit = text.splitlines(), False
    for i, line in enumerate(lines):
        m = re.match(rf"^(\s*){key}\s*=", line, re.IGNORECASE)
        if m and active(line):
            lines[i] = f"{m.group(1)}{key} = {value}"
            hit = True
    if not hit:
        for i, line in enumerate(lines):
            if line.strip().lower() == f"&{group}":
                lines.insert(i + 1, f"    {key} = {value}")
                break
        else:
            sys.exit(f"namelist group &{group} not found")
    return "\n".join(lines) + "\n"


def ints(value):
    return [int(v) for v in re.split(r"[,\s]+", value.strip()) if v]


def expand_periods(fp, npxs, npys):
    """Per simulated year: (forcing year, npx, npy), following the driver loop."""
    out = []
    for p in range(len(fp) // 3):
        y1, y2, ncyc = fp[3 * p:3 * p + 3]
        step = -1 if y1 > y2 else 1
        for _ in range(ncyc):
            out += [(y, npxs[p], npys[p]) for y in range(y1, y2 + step, step)]
    return out


def compress_periods(seq):
    """Rebuild forc_periods/npxs/npys from per-year entries (one cycle each)."""
    fp, px, py = [], [], []
    start = 0
    for i in range(1, len(seq) + 1):
        if i < len(seq):
            step = seq[i][0] - seq[i - 1][0]
            first_step = seq[start + 1][0] - seq[start][0] if i - start >= 2 else step
            same = abs(step) == 1 and step == first_step and seq[i][1:] == seq[start][1:]
            if same:
                continue
        fp += [seq[start][0], seq[i - 1][0], 1]
        px.append(seq[start][1])
        py.append(seq[start][2])
        start = i
    return fp, px, py


def history_files_for_continue(restart):
    """Read history filenames from the master checkpoint, not its .rh* tapes."""
    try:
        import netCDF4
    except ImportError:
        sys.exit("continue mode needs netCDF4 to read history file names from .rh* files")
    names = []
    with netCDF4.Dataset(restart) as nc:
        if "locfnh" not in nc.variables:
            sys.exit(f"master checkpoint lacks history filenames (locfnh): {restart}")
        for row in nc.variables["locfnh"][:]:
            name = b"".join(c for c in row.filled(b"") if c).decode().strip(" \x00")
            if name:
                names.append(Path(name).name)
    return names


def run_end_time(checkpoint, count, unit, step_seconds):
    """Advance from the checkpoint, retaining time of day for calendar units."""
    if unit == "nsteps":
        return checkpoint + timedelta(seconds=count * step_seconds)
    if unit == "ndays":
        return checkpoint + timedelta(days=count)
    months = count * (12 if unit == "nyears" else 1)
    year, month0 = divmod(checkpoint.year * 12 + checkpoint.month - 1 + months, 12)
    month = month0 + 1
    return checkpoint.replace(year=year, month=month,
                              day=min(checkpoint.day, calendar.monthrange(year, month)[1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--namelist", required=True, type=Path, help="parent case namelist")
    ap.add_argument("--restart", required=True, type=Path, help="parent checkpoint <case>.ecosim.r.YYYY-MM-DD-SSSSS.nc")
    ap.add_argument("--out", required=True, type=Path, help="new or empty run directory")
    ap.add_argument("--mode", choices=("branch", "continue"), default="branch")
    ap.add_argument("--case-name", help="branch case name (default <parent>.branch)")
    ap.add_argument("--stop-n", type=int, required=True, help="length counted from the checkpoint")
    ap.add_argument("--stop-option", default="nyears", choices=("nsteps", "ndays", "nmonths", "nyears"))
    ap.add_argument("--add-fields", default="", help="comma list appended to hist_fincl1 (branch only)")
    ap.add_argument("--copy-input", action="append", default=[], metavar="KEY",
                    help="copy this input (e.g. pft_file_in) instead of referencing it")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=VALUE",
                    help="extra &ecosim assignment, value written verbatim (e.g. disp_planttrait=.true.)")
    ap.add_argument("--exe", type=Path, help="run this executable after staging")
    args = ap.parse_args()
    if args.stop_n <= 0:
        ap.error("--stop-n must be positive")

    nml_path, rst = args.namelist.resolve(), args.restart.resolve()
    case_dir, out = nml_path.parent, args.out.resolve()
    m = RESTART_NAME.match(rst.name)
    if not rst.is_file():
        sys.exit(f"checkpoint file not found: {rst}")
    if not m:
        sys.exit(f"checkpoint name not <case>.ecosim.r.YYYY-MM-DD-SSSSS.nc: {rst.name}")
    if out.exists() and any(out.iterdir()):
        sys.exit(f"output directory is not empty: {out}")
    if out == case_dir:
        sys.exit("refusing to stage into the parent case directory")
    out.mkdir(parents=True, exist_ok=True)

    text = nml_path.read_text()
    parent_case = get_value(text, "case_name").strip("'\"")
    rst_year = int(m["y"])
    if not 0 <= int(m["s"]) < 86400:
        ap.error("checkpoint seconds must lie within its day")
    checkpoint_time = datetime(int(m["y"]), int(m["m"]), int(m["d"])) + timedelta(seconds=int(m["s"]))
    rst_date = checkpoint_time.strftime("%Y%m%d%H%M%S")

    # inputs: absolute paths, or copies the user can edit
    record = {"mode": args.mode, "parent_namelist": str(nml_path), "checkpoint": str(rst), "inputs": {}}
    for key in INPUT_KEYS:
        raw = get_value(text, key)
        if raw is None:
            continue
        val = raw.strip("'\"")
        if val.upper() in ("NO", ""):
            continue
        src = (case_dir / val).resolve()
        if not src.exists():
            sys.exit(f"{key} not found: {src}")
        if key in args.copy_input:
            shutil.copy2(src, out / src.name)
            dst = f"./{src.name}"
        else:
            dst = str(src)
        text = set_value(text, key, f"'{dst}'", "ecosim")
        record["inputs"][key] = dst

    if args.mode == "branch":
        case = args.case_name or f"{parent_case}.branch"
        if case == parent_case:
            sys.exit("branch case name must differ from the parent case name")
        shutil.copy2(rst, out / rst.name)
        start_year = int(get_value(text, "start_date").strip("'\"")[:4])
        seq = expand_periods(ints(get_value(text, "forc_periods")),
                             ints(get_value(text, "npxs")), ints(get_value(text, "npys")))
        idx = rst_year - start_year
        if not 0 <= idx < len(seq):
            sys.exit(f"checkpoint year {rst_year} outside the parent forcing sequence")
        fp, px, py = compress_periods(seq[idx:])
        for key, value in (("case_name", f"'{case}'"), ("continue_run", ".false."),
                           ("finidat", f"'./{rst.name}'"), ("start_date", f"'{rst_date}'"),
                           ("forc_periods", ",".join(map(str, fp))),
                           ("npxs", ", ".join(map(str, px))), ("npys", ", ".join(map(str, py)))):
            text = set_value(text, key, value, "ecosim")
        fields = [f for f in (get_value(text, "hist_fincl1") or "").replace("'", "").replace('"', "").split(",") if f.strip()]
        fields += [f.strip() for f in args.add_fields.split(",") if f.strip() and f.strip() not in fields]
        if fields:
            text = set_value(text, "hist_fincl1", ",".join(f"'{f}'" for f in fields), "ecosim")
        record.update(case_name=case, start_date=rst_date, forcing_years=[s[0] for s in seq[idx:]][:12])
    else:
        if args.add_fields:
            sys.exit("continue mode cannot change history fields; use --mode branch")
        stem = rst.name.replace(".ecosim.r.", ".ecosim.rh")
        tag = rst.name.split(".ecosim.r.")[1]
        rh_files = sorted(rst.parent.glob(f"{parent_case}.ecosim.rh*.{tag}"))
        if not rh_files:
            sys.exit(f"no history restart files ({stem.split('.rh')[0]}.rh*.{tag}) next to the checkpoint")
        for f in [rst, *rh_files]:
            shutil.copy2(f, out / f.name)
        for name in history_files_for_continue(rst):
            shutil.copy2(rst.parent / name, out / name)
        (out / "rpointer.esim").write_text(f"./{rst.name}\n")
        text = set_value(text, "continue_run", ".true.", "ecosim")
        record["copied"] = [f.name for f in rh_files]

    # The standalone driver replays its clock from January 1 before loading
    # the checkpoint. Its stop counter includes those skipped hours. Express
    # the requested additional duration as an absolute model-step endpoint.
    raw_step = get_value(text, "delta_time") or "1800."
    step_seconds = float(raw_step.split("!", 1)[0].strip().replace("D", "E").replace("d", "e"))
    if step_seconds <= 0:
        ap.error("delta_time must be positive")
    target_time = run_end_time(checkpoint_time, args.stop_n, args.stop_option, step_seconds)
    origin = datetime(int(get_value(text, "start_date").strip("'\"")[:4]), 1, 1)
    target_steps = (target_time - origin).total_seconds() / step_seconds
    if not target_steps.is_integer() or target_steps <= 0:
        ap.error("requested endpoint must align with the model timestep after its clock origin")
    text = set_value(text, "stop_option", "'nsteps'", "ecosim_time")
    text = set_value(text, "stop_n", str(int(target_steps)), "ecosim_time")
    record.update(requested_stop_n=args.stop_n, requested_stop_option=args.stop_option,
                  expected_end=target_time.isoformat(), expected_nstep=int(target_steps))
    for item in args.set:
        key, _, value = item.partition("=")
        text = set_value(text, key.strip(), value.strip(), "ecosim")

    (out / "runfile.namelist").write_text(text)
    (out / "restart_run.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))

    if args.exe:
        log = out / "log"
        with open(log, "w") as fh:
            rc = subprocess.call([str(args.exe.resolve()), "runfile.namelist"], cwd=out, stdout=fh, stderr=subprocess.STDOUT)
        body = log.read_text(errors="replace")
        loaded = "restFile_open Successfully" in body
        steps = re.findall(r"wrote out restart data at nstep\s*=\s*(\d+)", body)
        completed = bool(steps) and int(steps[-1]) == int(target_steps)
        stamp = target_time.strftime("%Y-%m-%d-") + f"{target_time.hour * 3600 + target_time.minute * 60 + target_time.second:06d}"
        completed = completed and (out / f"{parent_case if args.mode == 'continue' else case}.ecosim.r.{stamp}.nc").is_file()
        print(f"exit={rc} checkpoint_loaded={loaded} completed={completed} expected_nstep={int(target_steps)} log={log}")
        if rc != 0 or not loaded or not completed or "ENDRUN" in body or "Fortran runtime error" in body:
            sys.exit(1)


if __name__ == "__main__":
    main()
