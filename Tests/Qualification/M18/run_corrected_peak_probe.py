"""Repeat only the Intel service probes affected by the diagnostic policy fix."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def signature(data):
    decisions = {
        name: {k: v for k, v in row.items() if k not in ('main_preflight_ms', 'main_accept_ms')}
        for name, row in data['outcomes'].items()
    }
    return (data['identities'], decisions, data['file_sha256'],
            data['source_lufs'], data['current_lufs'], data['base_build'])


def main():
    parser = argparse.ArgumentParser()
    for name in ('candidate', 'original', 'probe', 'output'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--probe-sha256', required=True)
    args = parser.parse_args()
    root, original, probe, out = (getattr(args, n).resolve()
                                  for n in ('candidate', 'original', 'probe', 'output'))
    out.mkdir(parents=True, exist_ok=False)
    identity = json.loads((root / 'identity.json').read_text())
    baseline = {r['name']: r for r in json.loads((original / 'commands.json').read_text())}
    rows, failure = [], ''
    try:
        assert platform.system() == 'Darwin' and platform.machine() == 'x86_64'
        assert sha(probe) == args.probe_sha256
        for channels in (2, 6):
            expected = None
            origins = []
            for mode in ('serial', 'app-mode'):
                name = 'peaks-' + str(channels) + '-' + mode
                report = original / (name + '.json')
                assert baseline[name]['exit_code'] == 0
                assert sha(report) == baseline[name]['report_sha256']
                data = json.loads(report.read_text())
                assert data['base_build'] == identity['measurement_build'] and len(data['outcomes']) == 7
                current = signature(data)
                if expected is None:
                    expected = current
                assert current == expected
                origins.append(dict(mode=mode, report_sha256=sha(report),
                                    seconds=baseline[name]['seconds']))
            name = 'peaks-' + str(channels) + '-bounded-corrected'
            report = out / (name + '.json')
            fixture = root / 'peak-fixtures' / ('programme-' + str(channels) + '.wav')
            assert sha(fixture) == identity['artifacts_sha256'][str(fixture.relative_to(root))]
            command = ['/usr/bin/arch', '-x86_64', str(probe), str(fixture), str(report), '--decisions']
            print('START', name, flush=True)
            started = time.monotonic()
            with (out / (name + '.log')).open('w') as log:
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            data = json.loads(report.read_text())
            exact = result.returncode == 0 and signature(data) == expected
            rows.append(dict(channels=channels, command=command, exit_code=result.returncode,
                             seconds=time.monotonic() - started, exact_match=exact,
                             report_sha256=sha(report), baseline_receipts=origins))
            assert exact, name + ': exact serial/app comparison failed'
            print('PASS', name, round(rows[-1]['seconds'], 2), flush=True)
        assert len(rows) == 2
    except Exception as error:
        failure = repr(error)
    summary = dict(passed=not failure, failure=failure, architecture=platform.machine(),
                   build=identity['measurement_build'], identity_sha256=sha(root / 'identity.json'),
                   corrected_probe_sha256=sha(probe), checker_sha256=sha(__file__), cases=rows,
                   execution='Private diagnostic forces bounded service on Intel; production app remains serial',
                   supersedes='Only the originally mislabeled Intel bounded-service receipts',
                   physical_devices_opened=False)
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))
    return 1 if failure else 0


if __name__ == '__main__':
    raise SystemExit(main())
