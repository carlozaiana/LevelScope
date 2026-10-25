#!/usr/bin/env python3
"""Compare restored samples, exact availability, and actual-platform M9 guards."""
import argparse
import hashlib
import json
import math
import pathlib

from compare_m7_reports import compare_capture
from run_m9_qualification import validate_results


def compare(left, right):
    roots = (left, right)
    reports = [json.loads((p/'programmes.json').read_text()) for p in roots]
    cores = [json.loads((p/'long-programme.json').read_text()) for p in roots]
    identities = [json.loads((p/'identity.json').read_text()) for p in roots]
    environments = [json.loads((p/'environment.json').read_text()) for p in roots]
    for c, r, i, e in zip(cores, reports, identities, environments):
        validate_results(r, c, i, e)
    if ({r['architecture'] for r in reports} != {'arm64', 'x86_64'}
            or identities[0] != identities[1]):
        raise ValueError('Different candidate or missing actual architecture')
    if cores[0]['long_programme']['programme_identity'] != cores[1]['long_programme']['programme_identity']:
        raise ValueError('Different real long-media fixture')
    rows = []
    for a, b in zip(reports[0]['cases'], reports[1]['cases']):
        for key in ('kind', 'rate', 'channels', 'gain_db', 'frames', 'file', 'lra_stable', 'sample_peak_frames', 'true_peak_eighth_frames', 'ordered_roles', 'integrated_included', 'lra_included'):
            if a[key] != b[key]:
                raise ValueError('Different capture case identity: ' + key)
        values = []
        for root, row in zip(roots, (a, b)):
            if pathlib.Path(row['file']).name != row['file']:
                raise ValueError('Invalid sample path')
            data = (root/'samples'/row['file']).read_bytes()
            if (len(data) != row['frames'] * row['channels'] * 4
                    or hashlib.sha256(data).hexdigest() != row['sha256']):
                raise ValueError('Capture integrity mismatch')
            values.append(data)
        ulp = compare_capture(*values, a['gain_db'] == 0)
        errors = {}
        for name, metric in a['metrics'].items():
            other = b['metrics'][name]
            if metric['status'] != other['status'] or (metric['value'] is None) != (other['value'] is None):
                raise ValueError('Metric availability mismatch')
            if metric['value'] is not None:
                difference = abs(metric['value']-other['value'])
                if not math.isfinite(difference) or difference > (.1 if name == 'lra' else .01):
                    raise ValueError('M1 cross-platform metric budget exceeded')
                errors[name] = difference
        rows.append(dict(case=a['file'], maximum_float32_ulp=ulp, metric_errors=errors))
    return dict(passed=True, source_digest=identities[0]['source_digest'], source_build=reports[0]['build'],
                maximum_float32_ulp=max(r['maximum_float32_ulp'] for r in rows), sample_cases=rows,
                scope='M9 native/long-media/resource/sample comparison; UI and bidirectional session exchange are separate')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('left', type=pathlib.Path)
    p.add_argument('right', type=pathlib.Path)
    p.add_argument('--output', required=True, type=pathlib.Path)
    args = p.parse_args()
    try:
        result = compare(args.left, args.right)
    except (ValueError, KeyError, TypeError, OSError) as error:
        result = dict(passed=False, error=str(error))
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
