#!/usr/bin/env python3
"""Actual-platform M8 qualification, using injected devices only. UI is separate."""
import argparse
import json
import pathlib
import shutil
import subprocess
import sys
import time

from collect_m1_environment import sha256


def validate_results(core, playback, identity, environment):
    arch = environment['execution_arch']
    prefix = '15.' if arch == 'arm64' else '10.13.'
    if (arch not in ('arm64', 'x86_64') or not environment['platform_guard_passed']
            or not environment['product_version']['output'].startswith(prefix)
            or environment['artifact_integrity_errors']
            or environment['source_digest'] != identity['source_digest']):
        raise ValueError('Actual platform / source / artifact guard failed')
    for report in (core, playback):
        if (not report.get('passed') or report['build'] != identity['measurement_build']
                or report['architecture'] != arch):
            raise ValueError('Failed or mismatched M8 suite')
    p = core['performance']
    for key, cap in [('capture_ms', 100), ('ui_adoption_ms', 100), ('save_ms', 5000),
                     ('parse_ms', 5000), ('active_cancel_ms', 1000),
                     ('rss_growth_bytes', 256 * 1024 * 1024),
                     ('session_bytes', 128 * 1024 * 1024), ('ui_callback_p95_ms', 100)]:
        if not 0 <= p[key] <= cap:
            raise ValueError('Persistence budget failed: ' + key)
    if p['ui_callback_count'] < 5:
        raise ValueError('Missing responsiveness evidence')
    media = core['long_media']
    if (media['frames'] != 518400000 or media['rate'] != 48000 or media['channels'] != 8
            or not 0 <= media['active_cancel_ms'] <= 1000 or media['cancel_progress'] <= 0):
        raise ValueError('Incomplete real three-hour verification/cancellation')
    expected = {(rate, channels, gain) for rate in (44100, 48000, 96000, 192000)
        for channels in (2, 6, 8) for gain in (0., -6.020599913279624, 6.020599913279624, 24.)}
    cases = playback['sample_cases']
    if (len(cases) != 48 or {(r['rate'], r['channels'], r['gain_db']) for r in cases} != expected
            or playback['reference_cases'] != 5 or playback['physical_devices_opened']):
        raise ValueError('Incomplete restored audition/reference matrix')
    if not playback['checks'] or not all(r['passed'] for r in playback['checks']):
        raise ValueError('Failed playback assertions')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--platform', required=True, choices=('sequoia', 'high-sierra'))
    p.add_argument('--output', required=True, type=pathlib.Path)
    args = p.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    identity = json.loads(pathlib.Path('identity.json').read_text())
    runs = []

    def run(name, argv, timeout=7200):
        start = time.monotonic()
        with (output / (name + '.log')).open('w') as log:
            code = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, timeout=timeout).returncode
        runs.append(dict(name=name, argv=argv, exit_code=code, seconds=time.monotonic()-start,
                         log_sha256=sha256(output / (name + '.log'))))
        (output / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        print(name, 'PASS' if code == 0 else 'FAIL', flush=True)
        return code == 0

    passed = run('regressions', [sys.executable, '-B', 'Tests/run_m7_qualification.py',
        '--platform', args.platform, '--output', str(output / 'regressions'), '--completion-regressions'])
    environment = json.loads((output / 'regressions/environment.json').read_text())
    # Never continue a suite after a failed actual-platform/artifact guard.
    if not environment['platform_guard_passed'] or environment['artifact_integrity_errors']:
        raise SystemExit('Actual-platform/artifact guard failed')
    for name in ('environment.json', 'identity.json'):
        shutil.copy2(output / 'regressions' / name, output / name)
    arch = 'arm64' if args.platform == 'sequoia' else 'x86_64'
    def executable(name):
        return ['/usr/bin/arch', '-' + arch, str(pathlib.Path('bin', name).resolve())]
    passed &= run('sessions', executable('LevelScopeSessionTests') + ['--output', str(output / 'sessions.json'),
        '--exchange-out', str(output / 'exchange-out'), '--long-media', str(pathlib.Path('m7-fixtures/three-hour-71-48000.flac').resolve())])
    passed &= run('restored-playback', executable('LevelScopeSessionPlaybackTests') + ['--report', str(output / 'restored-playback.json'),
        '--samples', str(output / 'restored-samples'), '--video', str(pathlib.Path('video-fixtures').resolve()),
        '--prores', str(pathlib.Path('m7-video-fixtures').resolve())])
    passed &= run('m8-runner', [sys.executable, '-B', 'Tests/M8QualificationRunnerTests.py'])
    error = ''
    try:
        validate_results(json.loads((output/'sessions.json').read_text()),
            json.loads((output/'restored-playback.json').read_text()), identity, environment)
    except (ValueError, KeyError, OSError, TypeError) as e:
        error = str(e); passed = False
    (output/'summary.json').write_text(json.dumps(dict(automated_passed=bool(passed), error=error,
        qualification_complete=False, interactive_ui='Separate evidence required',
        cross_mac_reopen='Bidirectional exchange is separate required evidence',
        physical_devices_opened=False, source_digest=identity['source_digest'],
        identity_sha256=sha256(pathlib.Path('identity.json'))), indent=2)+'\n')
    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())
