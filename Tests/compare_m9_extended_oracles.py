#!/usr/bin/env python3
"""Full-resolution stable-LRA/gating and silence oracle comparison; no devices."""
import argparse
import json
import math
from pathlib import Path


def validate(report):
    arch = report['architecture']
    prefix = '15.' if arch == 'arm64' else '10.13.'
    if (arch not in ('arm64', 'x86_64') or not report['os'].split()[-1].startswith(prefix)
            or not report['passed'] or report['physical_devices_opened']
            or not report['checks'] or not all(c['passed'] for c in report['checks'])):
        raise ValueError('Extended oracle failed or not executed on an actual qualified platform')
    rows = report['extended_oracles']
    expected = {(c, r, s, d) for c in (2, 6) for r in (44100, 48000, 96000, 192000)
                for s in (False, True) for d in ('Source', 'Current State')}
    if len(rows) != 32 or {(r['channels'], r['rate'], r['silence'], r['domain']) for r in rows} != expected:
        raise ValueError('Incomplete extended oracle matrix')
    for row in rows:
        if (not row['same_platform_full_result_exact']
                or row['frames'] != row['rate'] * (4 if row['silence'] else 64) + 17
                or row['gain_db'] != (0 if row['domain'] == 'Source' else -3.305074118944)):
            raise ValueError('Extended oracle identity/precision differs')
        for name in ('sample_peak_frames', 'true_peak_eighth_frames', 'sample_peaks', 'true_peaks'):
            if len(row[name]) != row['channels']:
                raise ValueError('Missing per-channel extended evidence')
        windows = row['windows']
        if {w[0] for w in windows} != {0, 1, 2} or any(not math.isfinite(w[2]) or w[2] < 0 for w in windows):
            raise ValueError('Missing/non-finite full-resolution windows')
        step = row['rate'] // 10
        for family, start, end in ((0, 4 * step, row['frames']), (1, 30 * step, row['frames']),
                                   (2, 30 * step, row['frames'] + 15 * step)):
            if [w[1] for w in windows if w[0] == family] != list(range(start, end + 1, step)):
                raise ValueError('Incomplete or displaced full-resolution window grid')
        metrics = row['metrics']
        if set(metrics) != {'integrated', 'lra', 'momentary', 'short_term', 'sample_peak', 'true_peak'}:
            raise ValueError('Missing metric availability')
        if any(m['value'] is not None and not math.isfinite(m['value']) for m in metrics.values()):
            raise ValueError('Non-finite metric')
        if row['silence']:
            if (metrics['integrated']['value'] is not None or metrics['lra']['value'] is not None
                    or any(m['status'] != 1 or m['value'] is not None for m in metrics.values())
                    or row['integrated_included'] or row['lra_included'] or any(w[2] for w in windows)):
                raise ValueError('Silence became a measured loudness or included window')
        elif (any(m['status'] != 0 or m['value'] is None for m in metrics.values())
                or not row['lra_stable'] or metrics['lra']['value'] is None or metrics['lra']['value'] <= 1
                or not 0 < row['integrated_included'] < sum(w[0] == 0 for w in windows)
                or not 0 < row['lra_included'] < sum(w[0] == 2 for w in windows)):
            raise ValueError('Missing stable nonzero LRA or excluded-window evidence')


def compare(left, right):
    for report in (left, right):
        validate(report)
    if left['build'] != right['build'] or {left['architecture'], right['architecture']} != {'arm64', 'x86_64'}:
        raise ValueError('Different implementation or missing actual platform')
    results = []
    for a, b in zip(left['extended_oracles'], right['extended_oracles']):
        for key in ('channels', 'rate', 'frames', 'silence', 'domain', 'gain_db', 'lra_stable',
                    'integrated_included', 'lra_included', 'sample_peak_frames', 'true_peak_eighth_frames'):
            if a[key] != b[key]:
                raise ValueError('Extended identity/availability/location mismatch: ' + key)
        errors = {}
        for name, metric in a['metrics'].items():
            other = b['metrics'][name]
            if metric['status'] != other['status'] or (metric['value'] is None) != (other['value'] is None):
                raise ValueError('Different metric status: ' + name)
            if metric['value'] is not None:
                error = abs(metric['value'] - other['value'])
                if error > (.1 if name == 'lra' else .01):
                    raise ValueError('M1 metric comparison budget exceeded: ' + name)
                errors[name] = error
        for key in ('integrated_gate', 'lra_gate', 'p10', 'p95'):
            x, y = a[key], b[key]
            if (x is None) != (y is None) or x is not None and (not math.isfinite(x+y) or abs(x-y) > .01):
                raise ValueError('Gate/percentile mismatch: ' + key)
        maximum = 0.
        if len(a['windows']) != len(b['windows']):
            raise ValueError('Window coverage mismatch')
        for x, y in zip(a['windows'], b['windows']):
            if x[:2] != y[:2] or (x[2] == 0) != (y[2] == 0):
                raise ValueError('Window position/zero-energy mismatch')
            if x[2] > 0:
                maximum = max(maximum, abs(10 * math.log10(x[2] / y[2])))
        if maximum > .01:
            raise ValueError('Full-resolution window comparison budget exceeded')
        for key in ('sample_peaks', 'true_peaks'):
            for x, y in zip(a[key], b[key]):
                if (not math.isfinite(x+y) or x < 0 or y < 0 or (x == 0) != (y == 0)
                        or x > 0 and abs(20 * math.log10(x/y)) > .01):
                    raise ValueError('Per-channel peak comparison budget exceeded')
        results.append(dict(channels=a['channels'],rate=a['rate'],silence=a['silence'],domain=a['domain'],
                            metric_errors=errors,maximum_window_lu_difference=maximum))
    return dict(passed=True,build=left['build'],cases=results,
                scope='32 full-resolution Source/Current dynamic/silence interleaved oracles; no playback/device claim')


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('left',type=Path);p.add_argument('right',type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    try:
        result=compare(json.loads(a.left.read_text()),json.loads(a.right.read_text()))
    except (ValueError,KeyError,TypeError,OSError,ZeroDivisionError) as error:
        result=dict(passed=False,error=str(error))
    a.output.write_text(json.dumps(result,indent=2)+'\n');return 0 if result['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
