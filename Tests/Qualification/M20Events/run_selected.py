"""Explicit silent Logic-stage manifest. No historical runner or physical device."""
import argparse
import json
from pathlib import Path
import subprocess
import time


def main():
    p = argparse.ArgumentParser()
    p.add_argument('package', type=Path)
    p.add_argument('root', type=Path)
    p.add_argument('--music', type=Path, required=True)
    p.add_argument('--cinema', type=Path, required=True)
    p.add_argument('--gating', type=Path, required=True)
    args = p.parse_args()
    package, root = args.package.resolve(), args.root.resolve()
    root.mkdir(parents=True, exist_ok=False)
    receipts, fixtures = root / 'receipts', root / 'fixtures'
    receipts.mkdir()
    fixtures.mkdir()
    selected = []

    def grouped(name, target, group, extra=()):
        selected.append((name, [str(package / 'bin' / target), '--group', group,
                               '--root', str(fixtures / name), '--output', str(receipts / (name + '.json')), *map(str, extra)]))

    for group in ('core', 'features', 'stereo', 'surround', 'heterogeneous', 'fallback', 'population'):
        grouped('event-' + group, 'LevelScopeEventPlanningTests', group)
    for mode, label in enumerate(('legacy20', 'event20', 'margin20', 'margin40', 'margin60')):
        grouped('music-' + label, 'LevelScopeEventPlanningTests', 'external',
                ('--input', args.music.resolve(), '--experiment', mode))
    for label, file in (('gating', args.gating), ('cinema', args.cinema)):
        grouped('event-' + label, 'LevelScopeEventPlanningTests', 'external',
                ('--input', file.resolve(), '--experiment', 4))
    for group in ('core', 'stereo', 'surround', 'heterogeneous', 'failures', 'audio', 'population'):
        grouped('region-' + group, 'LevelScopeRegionPlanningTests', group)
    for group in ('core', 'stereo', 'inputs'):
        grouped('lra-' + group, 'LevelScopeLraProposalTests', group)
    selected += [('region-ui', [str(package / 'bin/LevelScopeRegionUITests'), str(fixtures / 'region-ui')]),
                 ('timeline', [str(package / 'bin/LevelScopeTimelineProgressTests'), '--root', str(fixtures / 'timeline'),
                               '--output', str(receipts / 'timeline.json')])]
    (root / 'manifest.json').write_text(json.dumps(selected, indent=2) + '\n')
    runs = []
    for name, command in selected:
        print('START ' + name, flush=True)
        start = time.monotonic()
        with (root / (name + '.log')).open('w') as log:
            try:
                result = subprocess.run(['/usr/bin/time', '-l', *command], stdout=log, stderr=subprocess.STDOUT, timeout=600)
                code = result.returncode
            except subprocess.TimeoutExpired:
                code = 124
        if name == 'region-ui' and (fixtures / 'region-ui/result.json').is_file():
            (receipts / 'region-ui.json').write_bytes((fixtures / 'region-ui/result.json').read_bytes())
        runs.append(dict(name=name, exit_code=code, wall_seconds=time.monotonic()-start))
        (root / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        print(('PASS ' if code == 0 else 'FAIL ') + name, flush=True)
    return int(any(r['exit_code'] != 0 for r in runs))


if __name__ == '__main__':
    raise SystemExit(main())
