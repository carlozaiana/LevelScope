#!/usr/bin/env python3
"""Record OS/CPU/runtime and executable hashes without opening audio devices."""
import argparse
import datetime
import hashlib
import json
import pathlib
import platform
import subprocess
import sys


def command(*args):
    completed = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True, timeout=30)
    return dict(command=list(args), exit_code=completed.returncode, output=completed.stdout.strip())


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=pathlib.Path)
    parser.add_argument('--identity', required=True, type=pathlib.Path)
    parser.add_argument('--measurement', required=True, type=pathlib.Path)
    parser.add_argument('--scaffold', required=True, type=pathlib.Path)
    parser.add_argument('--app', required=True, type=pathlib.Path, help='App executable, not bundle directory')
    parser.add_argument('--require-os', required=True, choices=['high-sierra', 'sequoia'])
    args = parser.parse_args()
    product = command('/usr/bin/sw_vers', '-productVersion')
    machine = platform.machine()
    valid = product['exit_code'] == 0 and ((args.require_os == 'high-sierra' and product['output'].startswith('10.13.') and machine == 'x86_64')
              or (args.require_os == 'sequoia' and product['output'].startswith('15.') and machine == 'arm64'))
    identity = json.loads(args.identity.read_text())
    integrity_errors = []
    for relative, expected in identity['artifacts_sha256'].items():
        path = args.identity.parent / relative
        if not path.is_file() or sha256(path) != expected:
            integrity_errors.append(relative)
    report = dict(run_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  expected_platform=args.require_os, platform_guard_passed=valid,
                  product_version=product, os_build=command('/usr/bin/sw_vers', '-buildVersion'),
                  machine=machine, cpu=command('/usr/sbin/sysctl', '-n', 'machdep.cpu.brand_string'),
                  hardware_model=command('/usr/sbin/sysctl', '-n', 'hw.model'),
                  python_version=sys.version, python_executable=sys.executable,
                  source_build_identity=identity, identity_sha256=sha256(args.identity),
                  artifact_integrity_passed=not integrity_errors,
                  artifact_integrity_errors=integrity_errors, artifacts={})
    for name in ('measurement', 'scaffold', 'app'):
        path = getattr(args, name)
        report['artifacts'][name] = dict(sha256=sha256(path), bytes=path.stat().st_size,
                                         file=command('/usr/bin/file', str(path)),
                                         dependencies=command('/usr/bin/otool', '-L', str(path)))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print('Platform guard:', 'PASS' if valid else 'FAIL', product['output'], machine)
    return 0 if valid and not integrity_errors else 1


if __name__ == '__main__':
    sys.exit(main())
