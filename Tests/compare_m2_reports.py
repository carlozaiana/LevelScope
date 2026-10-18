#!/usr/bin/env python3
"""Compare two completed M2 bundles; profile allowances are never widened."""
import argparse
import json
import math

from pathlib import Path
from run_ebu_measurements import compare_results


def compare(left, right):
    report = dict(left=str(left), right=str(right), errors=[], maximum_differences={},
                  tolerances=dict(I=.01, M=.01, S=.01, TP=.01, LRA=.1))

    def read(root, name):
        return json.loads((root / (name + '.json')).read_text())

    for name, keys in [('environment', ('identity_sha256', 'source_digest')),
                       ('compliance', ('profile_digest', 'checks', 'engine', 'evaluator', 'build'))]:
        a, b = read(left, name), read(right, name)
        for key in keys:
            if key not in a or key not in b or a[key] != b[key]:
                report['errors'].append(name + ': ' + key)
    for root in (left, right):
        env, m2, summary = read(root, 'environment'), read(root, 'compliance'), read(root, 'summary')
        if (not env['platform_guard_passed'] or env['artifact_integrity_errors']
                or not m2['passed'] or not summary['automated_passed']
                or m2['architecture'] != env['execution_arch']):
            report['errors'].append(str(root) + ': qualification or architecture failed')
    platforms = {read(root, 'environment')['requested_platform'] for root in (left, right)}
    if platforms != {'high-sierra', 'sequoia'}:
        report['errors'].append('Requires actual High Sierra and arm64 Sequoia reports')
    for suite in ('itu', 'ebu'):
        a, b = read(left, suite), read(right, suite)
        if suite == 'ebu':
            report['errors'].extend(compare_results(a['results'], b))
        previous = {r['file']: r for r in b['results']}
        if set(previous) != {r['file'] for r in a['results']}:
            report['errors'].append(suite + ': file set mismatch')
        maximum = {m: 0 for m in report['tolerances']}
        for row in a['results']:
            old = previous.get(row['file'])
            if old is None or old['sha256'] != row['sha256']:
                report['errors'].append(suite + ': fixture mismatch ' + row['file'])
                continue
            for metric, tolerance in report['tolerances'].items():
                x, y = row['measured'][metric], old['measured'][metric]
                if (x is None) != (y is None):
                    report['errors'].append(suite + ': availability ' + row['file'] + ' ' + metric)
                elif x is not None:
                    if not math.isfinite(x) or not math.isfinite(y) or abs(x-y) > tolerance:
                        report['errors'].append(suite + ': tolerance ' + row['file'] + ' ' + metric)
                    maximum[metric] = max(maximum[metric], abs(x-y))
        report['maximum_differences'][suite] = maximum
    report['passed'] = not report['errors']
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('left', type=Path)
    parser.add_argument('right', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    result = compare(args.left, args.right)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('Cross-platform comparison:', 'PASS' if result['passed'] else 'FAIL', result['errors'])
    raise SystemExit(0 if result['passed'] else 1)
