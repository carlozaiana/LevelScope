#!/usr/bin/env python3
"""Add exact long-file channel locations to the frozen aggregate metric comparison."""
import hashlib
import json
import math
import pathlib
import sys

r = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(r / 'candidate-03/Tests'))
from compare_m10_reports import compare_peaks

paths = [r / (p + '-candidate-03-long/long-export/report.json') for p in ('sequoia', 'high-sierra')]
a, b = [json.loads(p.read_text()) for p in paths]
errors, rows = [], []
for domain in ('source', 'current', 'artifact'):
    x, y = a[domain], b[domain]
    try:
        assert x['frames'] == y['frames'] == '518400000'
        assert x['sample_rate'] == y['sample_rate'] == 48000
        assert len(x['per_channel_peaks']) == len(y['per_channel_peaks']) == 6
        assert x['measured_content_sha256'] == y['measured_content_sha256'] or domain == 'artifact'
        # Exact output locations are required for identical artifact bytes.
        # The approved policy permits quantized output differences; never turn
        # this additional report into a new identical-hash acceptance rule.
        identical_audio = domain != 'artifact' or x['sha256'] == y['sha256']
        maximum = compare_peaks(x['per_channel_peaks'], y['per_channel_peaks'], identical_audio)
        assert x['evaluation']['verdict'] == y['evaluation']['verdict']
        for name, metric in x['metrics'].items():
            other = y['metrics'][name]
            assert metric['status'] == other['status'] and metric['unit'] == other['unit']
            if metric['value'] is None:
                assert other['value'] is None
            else:
                assert math.isfinite(metric['value']) and math.isfinite(other['value'])
                assert abs(metric['value'] - other['value']) <= (.1 if name == 'lra' else .01)
        rows.append(dict(domain=domain, passed=True, maximum_per_channel_difference_db=maximum,
            roles_statuses_units_exact=True, exact_peak_locations_required=identical_audio,
            peak_locations_equal=all(p[k]==q[k] for p,q in zip(x['per_channel_peaks'],y['per_channel_peaks']) for k in ('sample_peak_frame','true_peak_eighth_frame'))))
    except (AssertionError, KeyError, TypeError) as exc:
        errors.append(domain + ': ' + repr(exc))
result = dict(passed=not errors, errors=errors, domains=rows,
    evidence=[dict(path=str(p.relative_to(r)), sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in paths],
    scope='Additional long-file metric and exact per-channel location comparison; comparison tolerances do not alter profile allowances.')
(r / 'm10-long-detail-comparison-candidate-03.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
raise SystemExit(0 if result['passed'] else 1)
