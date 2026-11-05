"""Explicit M19 selection. No recursive milestone runners or physical devices.

Use --groups for development only. Qualification requires every REQUIRED group
plus the separately recorded native silent UI check and cross-platform review.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import subprocess
import time

REQUIRED = ['core', 'stereo', 'surround', 'rates', 'failures', 'inputs', 'gating', 'real',
            'measurement', 'section', 'requirements']
p = argparse.ArgumentParser()
p.add_argument('binaries', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--cinema', type=Path, required=True)
p.add_argument('--groups', nargs='+', choices=REQUIRED, default=REQUIRED)
a = p.parse_args()
assert len(a.groups) == len(set(a.groups)) and not a.output.exists()
a.output.mkdir(parents=True)
bins = a.binaries.resolve()

def find(name):
    paths = list(bins.glob('**/' + name))
    paths = [v for v in paths if v.is_file()]
    assert len(paths) == 1, (name, paths)
    return paths[0]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

record = dict(required=REQUIRED, selected=a.groups, machine=platform.machine(),
              os=subprocess.check_output(['sw_vers'], text=True),
              cpu=subprocess.check_output(['sysctl', '-n', 'machdep.cpu.brand_string'], text=True).strip(),
              cinema_sha256=digest(a.cinema), cases=[], physical_devices_opened=False, video_loaded=False)
start = time.monotonic()
for group in a.groups:
    report = a.output.resolve() / (group + '.json')
    root = a.output.resolve() / group
    if group == 'measurement':
        binary = find('LevelScopeMeasurementTests')
        args = ['--report', str(report)]
    else:
        binary = find({'section': 'LevelScopeSectionProposalTests',
                       'requirements': 'LevelScopeJobRequirementsTests'}.get(group, 'LevelScopeLraProposalTests'))
        args = ['--root', str(root), '--output', str(report)]
        if group == 'section': args += ['--smoke']
        elif group != 'requirements': args += ['--group', group]
        if group == 'real': args += ['--input', str(a.cinema.resolve())]
    command = ['/usr/bin/time', '-l', str(binary), *args]
    began = time.monotonic()
    print('START', group, flush=True)
    with (a.output / (group + '.log')).open('w') as log:
        done = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    duration = time.monotonic() - began
    evidence = json.loads(report.read_text()) if report.exists() else {}
    raw = (a.output / (group + '.log')).read_text()
    rss = re.search(r'(\d+)\s+maximum resident set size', raw)
    row = dict(group=group, command=command, executable_sha256=digest(binary), exit_code=done.returncode,
               seconds=duration, peak_rss_bytes=int(rss.group(1)) if rss else None,
               passed=done.returncode == 0 and evidence.get('passed') is True,
               build=evidence.get('build'))
    record['cases'].append(row)
    record['seconds'] = time.monotonic() - start
    record['selected_passed'] = all(v['passed'] for v in record['cases'])
    record['required_automated_complete'] = (set(a.groups) == set(REQUIRED) and
        len(record['cases']) == len(REQUIRED) and record['selected_passed'])
    record['qualified'] = False  # Native UI, frozen identity and cross-platform review remain separate.
    (a.output / 'run.json').write_text(json.dumps(record, indent=2) + '\n')
    print('PASS' if row['passed'] else 'FAIL', group, round(duration, 3), flush=True)
raise SystemExit(0 if record['selected_passed'] else 1)
