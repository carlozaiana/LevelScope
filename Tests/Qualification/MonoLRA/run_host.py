"""Run the finite, device-free mono qualification on either actual macOS host."""
import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as file:
        for chunk in iter(lambda: file.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def verify(package, expected):
    manifest = package / 'identity.json'
    assert sha(manifest) == expected, 'Unexpected frozen package'
    identity = json.loads(manifest.read_text())
    for name, digest in identity['artifacts_sha256'].items():
        assert sha(package / name) == digest, name
    return identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('package', type=Path)
    parser.add_argument('identity_sha256')
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    package = args.package.resolve()
    output = args.output.resolve()
    identity = verify(package, args.identity_sha256)
    output.mkdir(parents=True, exist_ok=False)
    groups = []
    for name, target in [('mono', 'LraProposal'), ('mono-readiness', 'LraProposal'), ('mono-ui', 'RegionUI'),
                         ('inner-stereo', 'EventPlanning'), ('inner-surround', 'EventPlanning')]:
        group_root = output / name
        report = output / (name + '.json')
        if name == 'mono-ui':
            command = [group_root, '--mono-lra']
            report = group_root / 'result.json'
        else:
            command = ['--group', name, '--root', group_root, '--output', report]
        binary = package / 'bin' / ('LevelScope' + target + 'Tests')
        started = time.monotonic()
        with (output / (name + '.log')).open('w') as log:
            proc = subprocess.run([str(value) for value in [binary] + command], stdout=log, stderr=subprocess.STDOUT, timeout=360)
        data = json.loads(report.read_text()) if report.exists() else {}
        row = dict(group=name, exit=proc.returncode, assertions=len(data.get('checks', [])),
                   passed=proc.returncode == 0 and data.get('passed') is True and data.get('build') == identity['measurement_build'],
                   report_sha256=sha(report) if report.exists() else None, seconds=time.monotonic() - started)
        groups.append(row)
        print(name, row['passed'], row['assertions'], flush=True)
    verify(package, args.identity_sha256)
    summary = dict(passed=all(group['passed'] for group in groups), groups=groups, os=platform.mac_ver()[0],
                   os_build=subprocess.check_output(['sw_vers', '-buildVersion'], universal_newlines=True).strip(),
                   architecture=platform.machine(), implementation=identity['source_commit'], package=args.identity_sha256,
                   verified_entries_before_and_after=len(identity['artifacts_sha256']), physical_outputs_opened=False)
    (output/'receipt.json').write_text(json.dumps(summary, indent=2)+'\n')
    return 0 if summary['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
