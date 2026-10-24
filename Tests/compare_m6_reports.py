#!/usr/bin/env python3
"""Compare exact captured M6 streams; metric budgets remain in the M1–M5 comparators."""
import argparse
import array
import hashlib
import json
import math
import pathlib


def compare_capture(left, right, zero_gain):
    if len(left) != len(right) or len(left) % 4:
        raise ValueError('Missing or different sample counts')
    a, b = array.array('I'), array.array('I')
    a.frombytes(left)
    b.frombytes(right)
    maximum = 0
    for x, y in zip(a, b):
        if x & 0x7f800000 == 0x7f800000 or y & 0x7f800000 == 0x7f800000:
            raise ValueError('Non-finite capture')
        if zero_gain and x != y:
            raise ValueError('Zero gain must be bit exact across platforms')
        # Monotonic IEEE754 key, including sign and signed zero.
        key = lambda n: 0xffffffff - n if n & 0x80000000 else n + 0x80000000
        maximum = max(maximum, abs(key(x) - key(y)))
        if maximum > (0 if zero_gain else 1):
            raise ValueError('Nonzero gain exceeds one float32 ULP')
    return maximum


def compare(left, right):
    a = json.loads((left / 'audition.json').read_text())
    b = json.loads((right / 'audition.json').read_text())
    if not a.get('passed') or not b.get('passed') or a['build'] != b['build']:
        raise ValueError('Failed suite or different source/build identity')
    if {a['architecture'], b['architecture']} != {'arm64', 'x86_64'}:
        raise ValueError('Two actual executable architectures required; also inspect OS guard')
    environments = [json.loads((root / 'environment.json').read_text()) for root in (left, right)]
    for report, environment in zip((a, b), environments):
        version = environment['product_version']['output']
        prefix = '15.' if report['architecture'] == 'arm64' else '10.13.'
        if not environment['platform_guard_passed'] or not version.startswith(prefix) or environment['execution_arch'] != report['architecture']:
            raise ValueError('Actual platform guard failed')
    if environments[0]['source_digest'] != environments[1]['source_digest']:
        raise ValueError('Different source manifests')
    rows = a['sample_cases']
    if len(rows) != 124 or len(b['sample_cases']) != len(rows):
        raise ValueError('Incomplete capture matrix')
    results = []
    for first, second in zip(rows, b['sample_cases']):
        for key in ('kind', 'rate', 'channels', 'gain_db', 'frames', 'file'):
            if first[key] != second[key]:
                raise ValueError('Different case identity: ' + key)
        if first['lra_stable'] != second['lra_stable']:
            raise ValueError('Different LRA stability')
        metric_errors = {}
        for name, metric in first['metrics'].items():
            other = second['metrics'][name]
            if metric['status'] != other['status'] or (metric['value'] is None) != (other['value'] is None):
                raise ValueError('Different metric availability/status')
            if metric['value'] is not None:
                error = abs(metric['value'] - other['value'])
                if not math.isfinite(error) or error > (0.1 if name == 'lra' else 0.01):
                    raise ValueError('Captured metric exceeds M1 cross-platform budget: ' + name)
                metric_errors[name] = error
        data = []
        for root, row in ((left, first), (right, second)):
            if pathlib.Path(row['file']).name != row['file']:
                raise ValueError('Invalid capture filename')
            raw = (root / 'samples' / row['file']).read_bytes()
            if len(raw) != row['frames'] * row['channels'] * 4 or hashlib.sha256(raw).hexdigest() != row['sha256']:
                raise ValueError('Capture integrity mismatch')
            data.append(raw)
        maximum = compare_capture(*data, first['gain_db'] == 0)
        results.append(dict(case=first['file'], kind=first['kind'], maximum_float32_ulp=maximum, metric_errors=metric_errors))
    return dict(passed=True, source_build=a['build'], sample_cases=results,
                maximum_float32_ulp=max(r['maximum_float32_ulp'] for r in results),
                scope='Captured sample parity only; OS, metrics, UI and physical video require separate evidence')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('left', type=pathlib.Path)
    p.add_argument('right', type=pathlib.Path)
    p.add_argument('--output', type=pathlib.Path, required=True)
    args = p.parse_args()
    try:
        result = compare(args.left, args.right)
    except (ValueError, KeyError, OSError) as error:
        result = dict(passed=False, error=str(error))
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
