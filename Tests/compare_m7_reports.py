#!/usr/bin/env python3
"""Compare exact captured M7 streams; metric budgets remain in the M1–M5 comparators."""
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
    a = json.loads((left / 'multichannel.json').read_text())
    b = json.loads((right / 'multichannel.json').read_text())
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
    if a['profile_digest'] != 'ef6f5cc071a8c9d5a347da0dd61fcc46e9f24da0aed787bc9f3ea48cb3db7e66' or a['profile_digest'] != b['profile_digest'] or a['layout_capability'] != 'pcm-7.1-bs2051-system-i-v1' or a['layout_capability'] != b['layout_capability'] or a['mapping_version'] != b['mapping_version']:
        raise ValueError('Different profile or layout capability')
    rows = a['sample_cases']
    if len(rows) != 218 or len(b['sample_cases']) != len(rows):
        raise ValueError('Incomplete capture matrix')
    required = {(f'M7 format {fmt}', rate, gain) for fmt in range(9)
                for rate in (44100, 48000, 96000, 192000)
                for gain in (0., -6.020599913279624, 6.020599913279624, -3.305074118944, -60., 24.)}
    actual = {(r['kind'], r['rate'], r['gain_db']) for r in rows if r['kind'].startswith('M7 format ')}
    muted = {r['kind'] for r in rows if r['kind'].startswith('muted highest peak ')}
    if actual != required or muted != {'muted highest peak 3', 'muted highest peak 4'}:
        raise ValueError('Incomplete M7 format/rate/gain or muted-peak coverage')
    if any(r['channels'] != 8 or r['frames'] <= 0 for r in rows):
        raise ValueError('M7 captures must contain all eight channels')
    if len({r['file'] for r in rows}) != len(rows):
        raise ValueError('Duplicate capture identity')
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
            raw = (root / 'm7-samples' / row['file']).read_bytes()
            if len(raw) != row['frames'] * row['channels'] * 4 or hashlib.sha256(raw).hexdigest() != row['sha256']:
                raise ValueError('Capture integrity mismatch')
            data.append(raw)
        maximum = compare_capture(*data, first['gain_db'] == 0)
        results.append(dict(case=first['file'], kind=first['kind'], maximum_float32_ulp=maximum, metric_errors=metric_errors))
    return dict(passed=True, source_build=a['build'], sample_cases=results,
                maximum_float32_ulp=max(r['maximum_float32_ulp'] for r in results),
                scope='Captured sample parity only; OS, metrics, UI and physical video require separate evidence')


def validate_long_seeks(rows):
    required = {(frame, gain) for frame in (0, 123456789, 518352000, 518399999)
                for gain in (0., -3.305074118944)}
    if len(rows) != 8 or {(r['frame'], r['gain_db']) for r in rows} != required:
        raise ValueError('Incomplete long Source/Current seek matrix')
    if any(not r['all_eight_exact'] or not math.isfinite(r['buffer_ms'])
           or not 0 <= r['buffer_ms'] <= 2000 for r in rows):
        raise ValueError('Long seek parity or responsiveness failed')


def compare_long(left, right):
    reports = [json.loads((root / 'multichannel-long.json').read_text()) for root in (left, right)]
    if not all(r['passed'] for r in reports) or reports[0]['build'] != reports[1]['build']:
        raise ValueError('Long-programme suite failed or differs in build')
    a, b = [r['long_programme'] for r in reports]
    for data in (a, b): validate_long_seeks(data['seek_cases'])
    if a['source_sha256'] != b['source_sha256'] or a['decoded_frames'] != 518400000 or b['decoded_frames'] != 518400000:
        raise ValueError('Different or incomplete real three-hour PCM fixture')
    errors = {}
    for domain in ('source', 'current'):
        first, second = a[domain], b[domain]
        if first['gain_db'] != second['gain_db'] or len(first['peaks']) != 8 or len(second['peaks']) != 8:
            raise ValueError('Different long processing plan or missing channel')
        values = {}
        for key in ('integrated_lufs', 'true_peak_dbtp', 'lra_lu'):
            x, y = first[key], second[key]
            if (x is None) != (y is None): raise ValueError('Different long metric availability')
            if x is not None:
                error = abs(x-y)
                if not math.isfinite(error) or error > (0.1 if key == 'lra_lu' else 0.01):
                    raise ValueError('Long metric outside cross-platform budget: '+key)
                values[key] = error
        for x, y in zip(first['peaks'], second['peaks']):
            for key in ('channel', 'sample_frame', 'true_peak_eighth_frame'):
                if x[key] != y[key]: raise ValueError('Different long channel/peak position')
            for key in ('sample_peak', 'true_peak'):
                if x[key] <= 0 or y[key] <= 0 or not math.isfinite(x[key]+y[key]) or abs(20*math.log10(x[key]/y[key])) > 0.01:
                    raise ValueError('Long per-channel peak outside budget')
        errors[domain] = values
    return dict(passed=True, source_sha256=a['source_sha256'], frames=a['decoded_frames'], metric_errors=errors)


def compare_video(left, right):
    reports = [json.loads((root / 'multichannel-video.json').read_text()) for root in (left, right)]
    if not all(r['passed'] for r in reports) or reports[0]['build'] != reports[1]['build']:
        raise ValueError('Video isolation suite failed or differs in build')
    expected = {'two-streams-audio.mov', 'vfr-exact.mov', 'proxy-with-tone.mov', 'lt-with-tone.mov', '422-with-tone.mov'}
    for report in reports:
        rows = report['video_cases']
        if len(rows) != 5 or {r['file'] for r in rows} != expected:
            raise ValueError('Incomplete M7 native video matrix')
        if any(r['selected_channels'] != 8 or r['observed_audio_outputs'] != 0 or r['observed_video_outputs'] <= 0
               or not r['source_and_current_exact_capture'] for r in rows):
            raise ValueError('Video source isolation failed')
    if {r['file']:r['sha256'] for r in reports[0]['video_cases']} != {r['file']:r['sha256'] for r in reports[1]['video_cases']}:
        raise ValueError('Different native video fixtures')
    return dict(passed=True, native_cases=5, selected_channels=8, observed_embedded_audio_outputs=0)

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('left', type=pathlib.Path)
    p.add_argument('right', type=pathlib.Path)
    p.add_argument('--output', type=pathlib.Path, required=True)
    args = p.parse_args()
    try:
        result = compare(args.left, args.right)
        result["long_programme"] = compare_long(args.left, args.right)
        result["native_video"] = compare_video(args.left, args.right)
        result["scope"] = "M7 automated cross-platform comparison; UI and physical timing remain separate"
    except (ValueError, KeyError, OSError) as error:
        result = dict(passed=False, error=str(error))
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    return 0 if result['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
