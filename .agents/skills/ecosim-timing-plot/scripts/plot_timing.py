#!/usr/bin/env python3
"""Render EcoSIM runtime-profiler output as a self-contained HTML report.

Reads <case>_outputs/timing/summary.csv and, when present, a per-step detail
file (time-<run-start>.txt, or the legacy time.txt). Standard library only.
"""
import argparse
import csv
import html
import json
import math
import sys
from pathlib import Path

MAX_POINTS = 1500   # per-series points after binning
MAX_SERIES = 7      # categorical slots; extra regions are listed in the table only
STEPS_PER_DAY = 24  # the standalone driver advances one hour per step


def read_summary(path):
    """Return (whole_run_seconds, rows) from summary.csv."""
    run_seconds = None
    lines = []
    with path.open(newline='') as stream:
        for line in stream:
            if line.startswith('#'):
                if 'Whole run seconds:' in line:
                    run_seconds = float(line.split(':', 1)[1])
                continue
            lines.append(line)
    if not lines:
        sys.exit(f'{path}: no data rows (did the run reach finalize_timer?)')
    rows = []
    for row in csv.DictReader(lines):
        rows.append({
            'timer': row['timer'],
            'calls': int(row['calls']),
            'total': float(row['total_s']),
            'mean': float(row['mean_s']),
            'min': float(row['min_s']),
            'max': float(row['max_s']),
            'percent': float(row['run_percent']),
        })
    return run_seconds, rows


def read_detail(path):
    """Return ({timer: {step: seconds}}, setup totals, sorted steps)."""
    per_step, setup, steps = {}, {}, set()
    with path.open(newline='') as stream:
        for row in csv.DictReader(stream):
            step, name, seconds = int(row['step']), row['timer'], float(row['seconds'])
            if step < 0:
                setup[name] = setup.get(name, 0.0) + seconds
                continue
            steps.add(step)
            per_step.setdefault(name, {})[step] = per_step.get(name, {}).get(step, 0.0) + seconds
    return per_step, setup, sorted(steps)


def bin_series(per_step, steps):
    """Bin every timer onto a common step axis; absent steps count as zero."""
    if not steps:
        return 1, [], {}
    first, last = steps[0], steps[-1]
    width = max(1, math.ceil((last - first + 1) / MAX_POINTS))
    if width > 1:
        # Whole days of hourly steps, so daily events (history writes) do not alias.
        width = STEPS_PER_DAY * math.ceil(width / STEPS_PER_DAY)
    nbins = (last - first) // width + 1
    starts = [first + i * width for i in range(nbins)]
    counts = [0] * nbins
    for step in steps:
        counts[(step - first) // width] += 1
    series = {}
    for name, values in per_step.items():
        sums, peaks = [0.0] * nbins, [0.0] * nbins
        for step, seconds in values.items():
            idx = (step - first) // width
            sums[idx] += seconds
            peaks[idx] = max(peaks[idx], seconds)
        series[name] = {
            'mean': [s / c if c else None for s, c in zip(sums, counts)],
            'max': [p if c else None for p, c in zip(peaks, counts)],
        }
    return width, starts, series


def pick_detail_file(timing_dir):
    candidates = sorted(timing_dir.glob('time-*.txt'))
    if candidates:
        return candidates[-1]   # YYYYMMDD-HHMMSS sorts chronologically
    legacy = timing_dir / 'time.txt'
    return legacy if legacy.exists() else None


def build_payload(summary_path, detail_path, title):
    run_seconds, rows = read_summary(summary_path)
    payload = {
        'title': title,
        'summaryFile': str(summary_path),
        'detailFile': str(detail_path) if detail_path else None,
        'runSeconds': run_seconds,
        'summary': rows,
        'detail': None,
    }
    if detail_path:
        per_step, setup, steps = read_detail(detail_path)
        width, starts, series = bin_series(per_step, steps)
        totals = {name: sum(v.values()) for name, v in per_step.items()}
        ranked = sorted((n for n in series if n != 'Timestep'), key=lambda n: -totals[n])
        payload['detail'] = {
            'steps': len(steps),
            'firstStep': steps[0] if steps else None,
            'lastStep': steps[-1] if steps else None,
            'binWidth': width,
            'binStarts': starts,
            'timestep': series.get('Timestep'),
            'regions': [{'name': n, **series[n]} for n in ranked[:MAX_SERIES]],
            'omitted': ranked[MAX_SERIES:],
            'setup': setup,
        }
    return payload


def render(payload):
    template = (Path(__file__).with_name('report_template.html')).read_text()
    data = json.dumps(payload, separators=(',', ':')).replace('</', '<\\/')
    return (template.replace('__TITLE__', html.escape(payload['title']))
                    .replace('__DATA__', data))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('timing', type=Path,
                        help='timing/ directory, a <case>_outputs/ directory, or summary.csv')
    parser.add_argument('--detail', type=Path,
                        help='detail file (default: newest time-*.txt next to summary.csv)')
    parser.add_argument('--no-detail', action='store_true', help='plot the summary only')
    parser.add_argument('-o', '--output', type=Path,
                        help='output HTML (default: <timing dir>/timing_report.html)')
    parser.add_argument('--title', help='report title (default: case directory name)')
    args = parser.parse_args()

    target = args.timing
    if target.is_dir() and (target / 'timing').is_dir():
        target = target / 'timing'
    summary = target if target.is_file() else target / 'summary.csv'
    if not summary.is_file():
        sys.exit(f'{summary}: not found')
    timing_dir = summary.parent
    detail = None if args.no_detail else (args.detail or pick_detail_file(timing_dir))
    if detail and not detail.is_file():
        sys.exit(f'{detail}: not found')
    title = args.title or timing_dir.resolve().parent.name.removesuffix('_outputs') or 'EcoSIM'
    output = args.output or timing_dir / 'timing_report.html'
    output.write_text(render(build_payload(summary, detail, title)))
    print(f'summary: {summary}')
    print(f'detail:  {detail or "none"}')
    print(f'report:  {output}')


if __name__ == '__main__':
    main()
