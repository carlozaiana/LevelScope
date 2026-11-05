"""Run only the affected hardware-independent monitoring/device/audition/preparation selection."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('candidate', type=Path)
    parser.add_argument('results', type=Path)
    args = parser.parse_args()
    candidate, results = args.candidate.resolve(), args.results.resolve()
    assert not results.exists(), 'Use fresh result paths'
    results.mkdir(parents=True)
    identity = json.loads((candidate / 'identity.json').read_text())
    for rel, digest in identity['artifacts_sha256'].items():
        assert sha(candidate / rel) == digest, rel
    commands = [
        ('devices', 'LevelScopeDeviceSettingsTests', [str(results / 'device-files')], results / 'device-files/result.json'),
        ('audition', 'LevelScopeAuditionTests', ['--report', str(results / 'audition.json')], results / 'audition.json'),
        ('preparation', 'LevelScopePreparationFeedbackTests', [str(results / 'preparation-files')], results / 'preparation-files/result.json'),
        ('session-playback', 'LevelScopeSessionPlaybackTests', ['--report', str(results / 'session-playback.json')], results / 'session-playback.json'),
    ]
    report = {'passed': False, 'os': platform.mac_ver()[0], 'architecture': platform.machine(),
              'candidate_identity_sha256': sha(candidate / 'identity.json'), 'build': identity['measurement_build'],
              'source_commit': identity['source_commit'], 'physical_devices_opened': False, 'groups': []}
    try:
        for name, executable, arguments, path in commands:
            command = [str(candidate / 'bin' / executable), *arguments]
            started = time.monotonic()
            with (results / (name + '.log')).open('w') as log:
                outcome = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=300)
            data = json.loads(path.read_text()) if path.exists() else {}
            passed = outcome.returncode == 0 and data.get('passed') and data.get('build') == identity['measurement_build'] and all(c['passed'] for c in data.get('checks', []))
            report['groups'].append({'name': name, 'command': command, 'returncode': outcome.returncode,
                                     'seconds': time.monotonic() - started, 'passed': bool(passed), 'checks': len(data.get('checks', []))})
            print(name, 'PASS' if passed else 'FAIL', flush=True)
            assert passed, name
        for rel, digest in identity['artifacts_sha256'].items():
            assert sha(candidate / rel) == digest, rel
        report['passed'] = True
    finally:
        (results / 'summary.json').write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
