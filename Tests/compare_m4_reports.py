#!/usr/bin/env python3
"""Compare M4 actual-platform qualification and enforce three-hour budgets."""
import argparse
import json
import math
from pathlib import Path
from compare_m3_reports import compare_m3

LIMITS = dict(analysis_ms=2000, worker_including_identity_ms=2000,
              analysis_storage_bytes=64*1024*1024, viewport_prepare_p95_ms=50,
              offscreen_timeline_paint_p95_ms=100, diagnostics_peak_rss_growth_bytes=128*1024*1024,
              active_cancel_ms=1000)


def compare_m4(left, right):
    report = compare_m3(left, right)
    records = [json.loads((root / 'diagnostics.json').read_text()) for root in (left, right)]
    for key in ('method', 'checks', 'build', 'engine'):
        if key not in records[0] or key not in records[1] or records[0][key] != records[1][key]:
            report['errors'].append('diagnostics: ' + key)
    for root, row in zip((left, right), records):
        env = json.loads((root / 'environment.json').read_text())
        if (not row.get('passed') or row.get('architecture') != env['execution_arch']
                or not row.get('checks') or not all(c.get('passed') for c in row['checks'])):
            report['errors'].append(str(root) + ': diagnostics failed or wrong architecture')
        performance = row.get('performance', {})
        if performance.get('generated_fixture_seconds') != 10800:
            report['errors'].append(str(root) + ': three-hour fixture missing')
        for key, limit in LIMITS.items():
            value = performance.get(key)
            if not isinstance(value, (int, float)) or not math.isfinite(value) or value > limit:
                report['errors'].append(str(root) + ': performance ' + key)
    report['performance_limits'] = LIMITS
    report['passed'] = not report['errors']
    return report


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('left', type=Path); p.add_argument('right', type=Path)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args(); result = compare_m4(args.left, args.right)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('M4 cross-platform comparison:', 'PASS' if result['passed'] else 'FAIL', result['errors'])
    raise SystemExit(0 if result['passed'] else 1)
