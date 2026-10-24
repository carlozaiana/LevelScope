#!/usr/bin/env python3
"""Read M6 physical reports; never start a device, modify evidence or relax gates."""
import argparse
import json
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('evidence', type=Path)
p.add_argument('--output', type=Path)
a = p.parse_args()
root = a.evidence
identity = json.loads((root / 'candidate-01-identity.json').read_text())
expected_binary = identity['artifacts_sha256']['LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification']
states = ('source', 'current')
rates = ('24000-1001', '24-1', '25-1', '30000-1001', '30-1')
matrix = [profile + '-' + rate for profile in ('proxy', 'lt', '422') for rate in rates] + ['51_SyncTest']
summary = dict(source_digest=identity['source_digest'], binary_sha256=expected_binary,
               scope='Physical audio-clock/UI-paint reports only; no optical/acoustic synchronization claim',
               complete=False, passed=False, missing=[], failures=[], platforms={})

def phase_identity(path, source_digest=False):
    if not path.exists():
        summary['missing'].append(str(path.relative_to(root)))
        return
    data = json.loads(path.read_text())
    failed = []
    if data.get('binary_sha256') != expected_binary:
        failed.append('binary_hash')
    if source_digest and data.get('source_digest') != identity['source_digest']:
        failed.append('source_digest')
    if failed:
        summary['failures'].append(dict(report=str(path.relative_to(root)), checks=failed))

def report(path, domain, platform, seconds):
    if not path.exists():
        summary['missing'].append(str(path.relative_to(root)))
        return None
    data = json.loads(path.read_text())
    expected_arch, expected_rate, expected_device = (
        ('arm64', 48000, 'Mac Studio Speakers') if platform == 'sequoia'
        else ('x86_64', 44100, 'Built-in Output'))
    checks = dict(passed=data.get('passed') is True,
                  native_arch=data.get('architecture') == expected_arch,
                  rate=data.get('sample_rate') == expected_rate,
                  device=data.get('device') == expected_device,
                  muted=data.get('all_routes_muted') is True,
                  domain=data.get('audition_state') == ('Source' if domain == 'source' else 'Current State'),
                  gain=data.get('processing_gain_db') == (0 if domain == 'source' else -3.305074118944),
                  measurement_build=data.get('measurement_build') == 'afe72a024daea1a94ff298567bb185d2c054dcdaea61a38bcdfc36b9a1df5e06',
                  duration=data.get('requested_seconds') == seconds and data.get('source_seconds', 0) >= seconds,
                  eligible=data.get('qualification_eligible') is True,
                  memory_method=data.get('memory_method') == 'prores-memory-v2-user-approved-2026-10-20',
                  memory=data.get('memory_passed') is True,
                  zero_skips=data.get('actual_paint_skipped_frames') == 0,
                  zero_blanks=data.get('blank_paints') == 0,
                  zero_bounds=data.get('bound_failures') == 0,
                  zero_deadlines=data.get('paint_deadline_failures') == 0,
                  zero_underruns=data.get('underrun_blocks') == 0)
    failed = [k for k, v in checks.items() if not v]
    if failed:
        summary['failures'].append(dict(report=str(path.relative_to(root)), checks=failed))
    return data

for platform in ('sequoia', 'high-sierra'):
    resumed = root / (platform + '-video-resume-2026-10-24-01')
    entry = dict(matrix={}, h264={}, endurance={})
    summary['platforms'][platform] = entry
    phase_identity(resumed / 'identity.json', source_digest=True)
    phase_identity(resumed / 'prores-endurance-current' / 'identity.json')
    for state in states:
        folder = (root / ('sequoia-prores-matrix-' + state + '-01') if platform == 'sequoia'
                  else resumed / ('prores-matrix-' + state))
        count = 0
        phase_identity(folder / 'identity.json')
        for case in matrix:
            if report(folder / (case + '.json'), state, platform, 100 if case == '51_SyncTest' else 18):
                count += 1
        entry['matrix'][state] = count
        for kind in ('cfr', 'vfr'):
            case = kind + '-' + state
            d = report(resumed / ('h264-' + case + '.json'), state, platform, 18)
            entry['h264'][case] = d is not None and d.get('passed') is True
    for profile in ('proxy', 'lt', '422'):
        d = report(resumed / 'prores-endurance-current' / (profile + '-1800.json'), 'current', platform, 1800)
        if d is not None:
            entry['endurance'][profile] = {k: d.get(k) for k in (
                'passed', 'source_seconds', 'actual_paint_skipped_frames', 'blank_paints',
                'bound_failures', 'paint_deadline_failures', 'underrun_blocks',
                'peak_sampled_playback_rss_bytes', 'first_minute_mean_non_reusable_resident_bytes',
                'last_minute_mean_non_reusable_resident_bytes', 'first_minute_mean_physical_footprint_bytes',
                'last_minute_mean_physical_footprint_bytes', 'memory_passed')}
summary['complete'] = not summary['missing']
summary['passed'] = summary['complete'] and not summary['failures']
if a.output:
    a.output.write_text(json.dumps(summary, indent=2) + '\n')
print(json.dumps(summary, indent=2))
raise SystemExit(0 if summary['passed'] else 1)
