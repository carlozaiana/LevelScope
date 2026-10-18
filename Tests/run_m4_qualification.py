#!/usr/bin/env python3
"""Run an isolated M4 bundle on actual High Sierra or arm64 Sequoia; no installs.

Run from the bundle root with a Python 3.8+ runtime. Silent interactive UI checks
are separate. Every run requires a fresh output directory and intact artifacts.
"""
import argparse
import datetime
import json
import pathlib
import platform
import subprocess
import sys
import time

from collect_m1_environment import command, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--platform', required=True, choices=['high-sierra', 'sequoia'])
    parser.add_argument('--output', required=True, type=pathlib.Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    identity_path = pathlib.Path('identity.json').resolve()
    identity = json.loads(identity_path.read_text())
    os_version = command('/usr/bin/sw_vers', '-productVersion')
    machine = platform.machine()
    arch = 'x86_64' if args.platform == 'high-sierra' else 'arm64'
    valid = machine == arch and os_version['exit_code'] == 0 and os_version['output'].startswith(
        '10.13.' if args.platform == 'high-sierra' else '15.')
    errors = [p for p, digest in identity['artifacts_sha256'].items()
              if not pathlib.Path(p).is_file() or sha256(pathlib.Path(p)) != digest]
    binaries = list(pathlib.Path('bin').glob('LevelScope*'))
    binaries.append(pathlib.Path('LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone'))
    environment = dict(run_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        requested_platform=args.platform, product_version=os_version,
        os_build=command('/usr/bin/sw_vers', '-buildVersion'), machine=machine, execution_arch=arch,
        cpu=command('/usr/sbin/sysctl', '-n', 'machdep.cpu.brand_string'),
        hardware_model=command('/usr/sbin/sysctl', '-n', 'hw.model'),
        python_version=sys.version, python_executable=sys.executable,
        platform_guard_passed=valid, artifact_integrity_errors=errors,
        identity_sha256=sha256(identity_path), source_digest=identity['source_digest'],
        artifacts={str(p): dict(sha256=sha256(p), bytes=p.stat().st_size,
                                file=command('/usr/bin/file', str(p)),
                                dependencies=command('/usr/bin/otool', '-L', str(p))) for p in binaries})
    (output / 'environment.json').write_text(json.dumps(environment, indent=2) + '\n')
    if not valid or errors:
        raise SystemExit('Platform or artifact integrity guard failed; see environment.json')
    (output / 'identity.json').write_text(identity_path.read_text())
    runs = []

    def run(name, argv):
        start = time.monotonic()
        with (output / (name + '.log')).open('w') as log:
            completed = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, timeout=1800)
        runs.append(dict(name=name, argv=argv, exit_code=completed.returncode,
                         seconds=time.monotonic() - start, log_sha256=sha256(output / (name + '.log'))))
        (output / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
        print(name, 'PASS' if completed.returncode == 0 else 'FAIL', flush=True)
        return completed.returncode == 0

    def executable(target):
        return ['/usr/bin/arch', '-' + arch, str(pathlib.Path('bin', target).resolve())]

    passed = run('diagnostics', executable('LevelScopeDiagnosticsTests') + ['--report', str(output / 'diagnostics.json')])
    passed &= run('proposal', executable('LevelScopeProposalTests') + ['--report', str(output / 'proposal.json')])
    passed &= run('compliance', executable('LevelScopeComplianceTests') + ['--report', str(output / 'compliance.json')])
    passed &= run('measurement', executable('LevelScopeMeasurementTests'))
    passed &= run('scaffold', executable('LevelScopeScaffoldTests'))
    fixtures = ('cfr-exact.mov', 'vfr-exact.mov', 'two-streams-audio.mov', 'timecode.mov',
                'cfr-exact.mp4', 'cfr-exact.m4v', 'half-hour.mov')
    passed &= run('video', executable('LevelScopeVideoTests') + [str(pathlib.Path('video-fixtures', p).resolve()) for p in fixtures])
    passed &= run('runner', [sys.executable, '-B', 'Tests/EbuReferenceRunnerTests.py'])
    passed &= run('m3-runner', [sys.executable, '-B', 'Tests/M3QualificationRunnerTests.py'])
    passed &= run('m4-runner', [sys.executable, '-B', 'Tests/M4QualificationRunnerTests.py'])
    for suite, script in [('itu', 'run_reference_measurements.py'), ('ebu', 'run_ebu_measurements.py')]:
        passed &= run(suite, [sys.executable, '-B', 'Tests/' + script, '--binary',
                              str(pathlib.Path('bin/LevelScopeMeasurementTests').resolve()),
                              '--vectors', str(pathlib.Path('vectors', suite).resolve()),
                              '--report', str(output / (suite + '.json')), '--compare',
                              str(pathlib.Path('baseline', suite + '.json').resolve()), '--arch', arch])
    (output / 'summary.json').write_text(json.dumps(dict(automated_passed=bool(passed),
        interactive_ui='Separate evidence required', platform=args.platform,
        source_digest=identity['source_digest'], identity_sha256=sha256(identity_path)), indent=2) + '\n')
    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())
