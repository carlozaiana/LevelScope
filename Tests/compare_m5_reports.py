#!/usr/bin/env python3
"""Compare M5 actual-platform qualification and enforce three-hour budgets."""
import argparse
import json
import math
from pathlib import Path
from compare_m4_reports import compare_m4

LIMITS = dict(aggregation_ms=2000, additional_section_storage_bytes=16*1024*1024,
              interaction_p95_ms=100, active_cancel_ms=1000)
VALUE_TOLERANCES = dict(minimum=1e-9, median=1e-9, maximum=1e-9, share=1e-12)


def compare_values(left, right):
    if not isinstance(left, list) or not isinstance(right, list) or len(left) != 24 or len(right) != 24:
        return False
    try:
        for a, b in zip(left, right):
            for key in ('sample_rate', 'signal_kind', 'begin', 'end'):
                if not isinstance(a[key], int) or a[key] != b[key]: return False
            for family in ('m', 's', 'lra'):
                x, y = a[family], b[family]
                for key in ('contained', 'crossing', 'padding', 'included', 'missing'):
                    if not isinstance(x[key], int) or x[key] < 0 or x[key] != y[key]: return False
                for key, tolerance in VALUE_TOLERANCES.items():
                    v, w = x[key], y[key]
                    if v is None or w is None:
                        if v is not None or w is not None: return False
                    elif (not isinstance(v, (int, float)) or not isinstance(w, (int, float))
                          or not math.isfinite(v) or not math.isfinite(w) or abs(v-w) > tolerance):
                        return False
    except (KeyError, TypeError):
        return False
    return True


def compare_m5(left, right):
    report = compare_m4(left, right)
    records = [json.loads((root / 'sections.json').read_text()) for root in (left, right)]
    if not compare_values(*(row.get('measured_section_values') for row in records)):
        report['errors'].append('sections: measured values, availability or membership mismatch')
    for key in ('method', 'checks', 'build', 'engine'):
        if key not in records[0] or key not in records[1] or records[0][key] != records[1][key]:
            report['errors'].append('sections: ' + key)
    for root, row in zip((left, right), records):
        env = json.loads((root / 'environment.json').read_text())
        if (not row.get('passed') or row.get('architecture') != env['execution_arch']
                or not row.get('checks') or not all(c.get('passed') for c in row['checks'])):
            report['errors'].append(str(root) + ': sections failed or wrong architecture')
        performance = row.get('performance', {})
        if performance.get('generated_fixture_seconds') != 10800 or performance.get('section_count') != 1000:
            report['errors'].append(str(root) + ': three-hour fixture missing')
        for key, limit in LIMITS.items():
            value = performance.get(key)
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0 or value > limit:
                report['errors'].append(str(root) + ': performance ' + key)
    report['section_performance_limits'] = LIMITS
    report['section_value_tolerances'] = VALUE_TOLERANCES
    report['passed'] = not report['errors']
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('left', type=Path); p.add_argument('right', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); result = compare_m5(args.left, args.right)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('M5 cross-platform comparison:', 'PASS' if result['passed'] else 'FAIL', result['errors'])
    raise SystemExit(0 if result['passed'] else 1)
