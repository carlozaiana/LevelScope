#!/usr/bin/env python3
"""Run pinned official ITU vectors; no downloads, audio devices or installations.

Extract the archive identified by Reference/itu-bs2217-manifest.json locally.
Audio is not part of the repository or shipped application. A missing/mismatched
fixture fails the run. Use --compare with another platform's saved report.
"""
import argparse
import concurrent.futures
import hashlib
import json
import pathlib
import platform
import subprocess
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--binary', required=True, type=pathlib.Path)
parser.add_argument('--vectors', required=True, type=pathlib.Path)
parser.add_argument('--report', required=True, type=pathlib.Path)
parser.add_argument('--compare', type=pathlib.Path)
parser.add_argument('--arch', choices=['arm64', 'x86_64'], help='macOS: execute universal test binary in this architecture')
args = parser.parse_args()
manifest_path = pathlib.Path(__file__).parent / 'Reference/itu-bs2217-manifest.json'
manifest = json.loads(manifest_path.read_text())


def check(case):
    path = args.vectors / case['file']
    if hashlib.sha256(path.read_bytes()).hexdigest() != case['sha256']:
        raise ValueError('Reference hash mismatch: ' + str(path))
    command = [str(args.binary.resolve()), '--file', str(path.resolve())]
    if args.arch:
        command = ['/usr/bin/arch', '-' + args.arch] + command
    output = subprocess.check_output(command, text=True)
    metrics = json.loads(output)
    expected, measured = case['integrated_lufs'], metrics['I']
    passed = measured is None if expected is None else measured is not None and abs(measured - expected) <= case['tolerance_lu']
    return dict(file=case['file'], sha256=case['sha256'], expected=expected, measured=metrics, passed=passed)


with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    results = list(pool.map(check, manifest['files']))
for result in results:
    print(('PASS' if result['passed'] else 'FAIL'), result['file'], result['measured']['I'], flush=True)
comparison_errors = []
if args.compare:
    previous = json.loads(args.compare.read_text())
    previous_results = {r['file']: r for r in previous['results']}
    for current in results:
        old = previous_results.get(current['file'])
        if old is None or old['sha256'] != current['sha256']:
            comparison_errors.append(current['file'] + ': missing/different reference')
            continue
        for metric, tolerance in [('I', .01), ('M', .01), ('S', .01), ('TP', .01), ('LRA', .1)]:
            a, b = current['measured'][metric], old['measured'][metric]
            if (a is None) != (b is None) or (a is not None and abs(a-b) > tolerance):
                comparison_errors.append(current['file'] + ': ' + metric)
report = dict(platform=platform.platform(), machine=platform.machine(), execution_arch=args.arch or platform.machine(),
              binary_sha256=hashlib.sha256(args.binary.read_bytes()).hexdigest(),
              manifest_sha256=hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
              results=results, comparison_errors=comparison_errors)
args.report.write_text(json.dumps(report, indent=2) + '\n')
print(f"{sum(r['passed'] for r in results)}/{len(results)} official vectors passed; {len(comparison_errors)} comparison errors")
sys.exit(0 if all(r['passed'] for r in results) and not comparison_errors else 1)
