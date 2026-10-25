#!/usr/bin/env python3
"""Frozen M9 actual-platform regression runner. Physical/UI/exchange gates are separate."""
import argparse
import json
import pathlib
import shutil
import subprocess
import sys
import time
from collect_m1_environment import sha256
from compare_m9_extended_oracles import validate as validate_extended_oracles


def validate_results(report, long_report, identity, environment):
    arch = environment['execution_arch']
    prefix = '15.' if arch == 'arm64' else '10.13.'
    if (arch not in ('arm64', 'x86_64') or not environment['platform_guard_passed']
            or not environment['product_version']['output'].startswith(prefix)
            or environment['artifact_integrity_errors'] or environment['source_digest'] != identity['source_digest']):
        raise ValueError('M9 actual-platform / artifact guard failed')
    for item in (report, long_report):
        if (not item.get('passed') or item['build'] != identity['measurement_build']
                or item['architecture'] != arch or item['physical_devices_opened']
                or not item['checks'] or not all(row['passed'] for row in item['checks'])):
            raise ValueError('M9 suite failed, unguarded or mismatched')
    expected = {(r, c, g, 'M9 mono group Current') for r in (44100, 48000, 96000, 192000)
                for c in (2, 6) for g in (0., 6.020599913279624, -6.020599913279624, -3.305074118944, -60., 24.)}
    expected |= {(r, c, 0., 'M9 mono group Source') for r in (44100, 48000, 96000, 192000) for c in (2, 6)}
    rows = report['cases']
    if len(rows) != 56 or {(r['rate'], r['channels'], r['gain_db'], r['kind']) for r in rows} != expected or report['reference_cases'] != 10:
        raise ValueError('M9 incomplete capture/reference matrix')
    migrations = report['m8_migrations']
    if len(migrations) != 3 or {r['channels'] for r in migrations} != {2, 6, 8} or not all(r['original_unchanged'] for r in migrations):
        raise ValueError('M9 missing preserved M8 migration matrix')
    media = long_report['long_programme']
    if media['frames'] != 518400000 or media['rate'] != 48000 or media['channels'] != 6:
        raise ValueError('M9 missing actual three-hour six-file evidence')
    for key, cap in [('save_ms', 5000), ('parse_ms', 5000), ('cancel_ms', 1000),
                     ('preparation_cancel_ms', 1000), ('max_seek_prefill_ms', 2000),
                     ('ui_callback_p95_ms', 100), ('rss_growth_bytes', 256 * 1024 * 1024),
                     ('session_bytes', 128 * 1024 * 1024)]:
        if not 0 <= media[key] <= cap:
            raise ValueError('M9 resource budget exceeded: ' + key)
    if media['ui_callback_count'] < 5 or media['owned_payload_bytes'] != media['original_payload_bytes']:
        raise ValueError('M9 missing responsiveness or bounded copy evidence')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--platform', required=True, choices=('sequoia', 'high-sierra'))
    parser.add_argument('--output', required=True, type=pathlib.Path)
    args = parser.parse_args()
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=False)
    identity = json.loads(pathlib.Path('identity.json').read_text()); runs = []
    def run(name, argv, timeout=14400):
        start = time.monotonic()
        with (output / (name + '.log')).open('w') as log:
            code = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, timeout=timeout).returncode
        runs.append(dict(name=name, argv=argv, exit_code=code, seconds=time.monotonic()-start,
                         log_sha256=sha256(output / (name + '.log'))))
        (output / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        print(name, 'PASS' if code == 0 else 'FAIL', flush=True)
        return code == 0
    passed = run('regressions', [sys.executable, '-B', 'Tests/run_m8_qualification.py', '--platform', args.platform, '--output', str(output / 'regressions')])
    environment = json.loads((output / 'regressions/environment.json').read_text())
    if not environment['platform_guard_passed'] or environment['artifact_integrity_errors']:
        raise SystemExit('M9 actual-platform / artifact guard failed')
    for name in ('environment.json', 'identity.json'):
        shutil.copy2(output / 'regressions' / name, output / name)
    arch = 'arm64' if args.platform == 'sequoia' else 'x86_64'
    executable = ['/usr/bin/arch', '-' + arch, str(pathlib.Path('bin/LevelScopeProgrammeTests').resolve())]
    passed &= run('programmes', executable + ['--output', str(output / 'programmes.json'), '--samples', str(output / 'samples'),
        '--video', str(pathlib.Path('video-fixtures').resolve()), '--prores', str(pathlib.Path('m7-video-fixtures').resolve()),
        '--exchange-out', str(output / 'exchange-out'), '--m8-session', str(pathlib.Path('m9-fixtures/m8-sessions').resolve())])
    passed &= run('extended-oracles', executable + ['--only', '--extended-oracle', '--output', str(output / 'extended-oracles.json')])
    passed &= run('long-programme', executable + ['--only', '--long-programme', str(pathlib.Path('m9-fixtures/long-48000/programme.levelscope').resolve()), '--output', str(output / 'long-programme.json')])
    passed &= run('m9-runner', [sys.executable, '-B', 'Tests/M9QualificationRunnerTests.py'])
    error = ''
    try:
        validate_results(json.loads((output/'programmes.json').read_text()), json.loads((output/'long-programme.json').read_text()), identity, environment)
        extended = json.loads((output/'extended-oracles.json').read_text())
        validate_extended_oracles(extended)
        if extended['build'] != identity['measurement_build'] or extended['architecture'] != arch:
            raise ValueError('Extended oracle identity mismatch')
    except (ValueError, KeyError, OSError, TypeError) as e:
        error = str(e); passed = False
    (output/'summary.json').write_text(json.dumps(dict(automated_passed=bool(passed), error=error, qualification_complete=False,
        interactive_ui='Separate required evidence', cross_mac_exchange='Separate bidirectional required evidence',
        physical_video='Separate coordinated required evidence', physical_multichannel_listening='Separate, not claimed',
        source_digest=identity['source_digest'], identity_sha256=sha256(pathlib.Path('identity.json'))), indent=2)+'\n')
    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())
