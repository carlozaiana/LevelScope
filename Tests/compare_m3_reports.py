#!/usr/bin/env python3
"""Compare completed actual-platform M3 results without widening profile rules."""
import argparse
import json
from pathlib import Path
from compare_m2_reports import compare


def compare_m3(left, right):
    report = compare(left, right)
    a = json.loads((left / 'proposal.json').read_text())
    b = json.loads((right / 'proposal.json').read_text())
    for key in ('method', 'checks', 'profile_digest', 'build', 'engine', 'scalar_outcomes'):
        if key not in a or key not in b or a[key] != b[key]:
            report['errors'].append('proposal: ' + key)
    for root, proposal in ((left, a), (right, b)):
        environment = json.loads((root / 'environment.json').read_text())
        if not proposal['passed'] or proposal['architecture'] != environment['execution_arch']:
            report['errors'].append('proposal: execution failed or wrong architecture')
    report['passed'] = not report['errors']
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('left', type=Path)
    parser.add_argument('right', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = compare_m3(args.left, args.right)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('M3 cross-platform comparison:', 'PASS' if result['passed'] else 'FAIL', result['errors'])
    raise SystemExit(0 if result['passed'] else 1)
