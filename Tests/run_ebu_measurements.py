#!/usr/bin/env python3
"""Qualify the offline kernel against the locally supplied, pinned EBU v05 WAVs.

No downloads, playback, device access or installation. Official minimum-requirements
checks and project sanity checks on supplementary alignment files stay separate.
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import json
import math
import pathlib
import platform
import subprocess
import sys


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def evaluate(rule, metrics):
    metric = rule['metric']
    kind = rule['kind']
    if kind == 'scalar':
        actual, expected = [metrics[metric]], [rule['expected']]
    elif kind == 'curve':
        curve = [p for p in metrics[metric + '_curve']
                 if rule['start_frame'] <= p['end_frame'] <= rule['end_frame']]
        frames = list(range(rule['start_frame'], rule['end_frame'] + 1, rule['step_frames']))
        if not frames or [p['end_frame'] for p in curve] != frames:
            raise ValueError('Missing, repeated or unordered ' + metric + ' curve points')
        actual = [p['value'] for p in curve]
        expected = [rule['expected']] * len(actual)
    elif kind == 'checkpoints':
        points = metrics['checkpoints']
        expected = rule['expected']
        frames = [rule['period_frames'] * (i + 1) for i in range(len(expected))]
        if not frames or [p['end_frame'] for p in points] != frames:
            raise ValueError('Missing, repeated or unordered checkpoints')
        actual = [p[metric] for p in points]
    else:
        raise ValueError('Unknown check kind: ' + kind)
    valid = bool(actual) and len(actual) == len(expected) and all(finite(v) for v in actual + expected)
    errors = [a - e for a, e in zip(actual, expected)] if valid else []
    passed = valid and all(rule['lower_error'] <= error <= rule['upper_error'] for error in errors)
    return dict(kind=kind, metric=metric, expected=rule['expected'],
                lower_error=rule['lower_error'], upper_error=rule['upper_error'],
                readings=len(actual), measured_min=min(actual) if valid else None,
                measured_max=max(actual) if valid else None,
                error_min=min(errors) if errors else None, error_max=max(errors) if errors else None,
                passed=passed)


def run_case(case, binary, vectors, arch):
    result = dict(file=case['file'], sha256=case['sha256'], group=case['group'], passed=False)
    try:
        path = vectors / case['file']
        if sha256(path) != case['sha256']:
            raise ValueError('Reference SHA-256 mismatch')
        command = [str(binary.resolve()), '--file', str(path.resolve()), '--roles', ','.join(case['roles'])]
        if any(r['kind'] == 'curve' for r in case['checks']):
            command.append('--curves')
        periods = {r['period_frames'] for r in case['checks'] if r['kind'] == 'checkpoints'}
        if len(periods) > 1:
            raise ValueError('Conflicting checkpoint intervals')
        if periods:
            command += ['--checkpoint-frames', str(next(iter(periods)))]
        if arch:
            command = ['/usr/bin/arch', '-' + arch] + command
        completed = subprocess.run(command, capture_output=True, text=True, timeout=300, check=True)
        metrics = json.loads(completed.stdout)
        if sha256(path) != case['sha256']:
            raise ValueError('Reference changed during measurement')
        for field in ('sample_rate', 'channels', 'frames', 'bits_per_sample'):
            if metrics[field] != case[field]:
                raise ValueError('Decoded metadata mismatch: ' + field)
        if metrics['roles'] != ','.join(case['roles']):
            raise ValueError('Channel-role mapping mismatch')
        if arch and metrics['execution_arch'] != arch:
            raise ValueError('Executable architecture mismatch')
        if not case['checks']:
            raise ValueError('No acceptance rules')
        checks = [evaluate(rule, metrics) for rule in case['checks']]
        result.update(measured=metrics,
                      checks=checks, passed=all(c['passed'] for c in checks),
                      measurement_json_sha256=hashlib.sha256(completed.stdout.encode()).hexdigest())
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as error:
        result['error'] = str(error)
        if isinstance(error, subprocess.CalledProcessError):
            result['stderr'] = error.stderr
    return result


def compare_results(results, previous):
    errors = []
    old_by_file = {r['file']: r for r in previous['results']}
    if set(old_by_file) != {r['file'] for r in results}:
        errors.append('Comparison file inventory differs')
    for current in results:
        old = old_by_file.get(current['file'])
        if not old or old['sha256'] != current['sha256'] or not old['passed'] or not current['passed']:
            errors.append(current['file'] + ': missing, changed or failed reference result')
            continue
        for metric, tolerance in [('I', .01), ('M', .01), ('S', .01), ('TP', .01), ('SP', .01), ('LRA', .1)]:
            a, b = current['measured'][metric], old['measured'][metric]
            if (a is None) != (b is None) or (a is not None and
                    (not finite(a) or not finite(b) or abs(a - b) > tolerance)):
                errors.append(current['file'] + ': ' + metric)
        # Include time-dependent criteria, not just whole-file summary metrics.
        if len(current['checks']) != len(old['checks']):
            errors.append(current['file'] + ': acceptance rules differ')
            continue
        for a, b in zip(current['checks'], old['checks']):
            if (a['kind'], a['metric'], a['readings']) != (b['kind'], b['metric'], b['readings']):
                errors.append(current['file'] + ': criterion coverage differs')
                continue
            tolerance = .1 if a['metric'] == 'LRA' else .01
            for key in ('measured_min', 'measured_max'):
                if not finite(a[key]) or not finite(b[key]) or abs(a[key] - b[key]) > tolerance:
                    errors.append(current['file'] + ': ' + a['metric'] + ' ' + key)
        for key in ('M_curve', 'S_curve', 'checkpoints'):
            a, b = current['measured'].get(key, []), old['measured'].get(key, [])
            if [p['end_frame'] for p in a] != [p['end_frame'] for p in b]:
                errors.append(current['file'] + ': ' + key + ' time positions differ')
                continue
            fields = [('value', .01)] if key.endswith('_curve') else [('M', .01), ('S', .01)]
            for x, y in zip(a, b):
                for metric, tolerance in fields:
                    v, w = x[metric], y[metric]
                    if (v is None) != (w is None) or (v is not None and
                            (not finite(v) or not finite(w) or abs(v - w) > tolerance)):
                        errors.append(current['file'] + ': ' + key + ' at ' + str(x['end_frame']))
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', required=True, type=pathlib.Path)
    parser.add_argument('--vectors', required=True, type=pathlib.Path, help='Directory containing the 70 WAVs and readme.txt')
    parser.add_argument('--report', required=True, type=pathlib.Path)
    parser.add_argument('--compare', type=pathlib.Path)
    parser.add_argument('--arch', choices=['arm64', 'x86_64'])
    args = parser.parse_args()
    manifest_path = pathlib.Path(__file__).parent / 'Reference/ebu-v05-manifest.json'
    manifest = json.loads(manifest_path.read_text())
    failures = []
    expected_names = {case['file'] for case in manifest['files']}
    actual_names = {p.name for p in args.vectors.glob('*.wav')}
    if len(expected_names) != 70 or len(manifest['files']) != 70 or actual_names != expected_names:
        failures.append('Fixture inventory mismatch: missing=' + str(sorted(expected_names - actual_names))
                        + ', extra=' + str(sorted(actual_names - expected_names)))
    try:
        if sha256(args.vectors / 'readme.txt') != manifest['readme_sha256']:
            failures.append('README hash mismatch')
    except OSError as error:
        failures.append(str(error))
    archive = args.vectors.parent / manifest['archive']
    archive_verified = archive.is_file() and sha256(archive) == manifest['archive_sha256']
    if archive.is_file() and not archive_verified:
        failures.append('Archive hash mismatch')
    binary_hash = sha256(args.binary)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        results = list(pool.map(lambda case: run_case(case, args.binary, args.vectors, args.arch), manifest['files']))
    for result in results:
        print(('PASS' if result['passed'] else 'FAIL'), result['file'], result.get('error', ''), flush=True)
        for check in result.get('checks', []):
            if not check['passed']:
                print(' ', json.dumps(check), flush=True)
    identities = {(r['measured']['engine_version'], r['measured']['build_identity'], r['measured']['execution_arch'])
                  for r in results if 'measured' in r}
    if len(identities) != 1 or sha256(args.binary) != binary_hash:
        failures.append('Measurement build changed or identities disagree')
    if args.compare:
        previous = json.loads(args.compare.read_text())
        if previous['manifest_sha256'] != sha256(manifest_path):
            failures.append('Comparison manifest differs')
        failures += compare_results(results, previous)
    groups = {}
    for group in ('minimum_requirements', 'supplementary_alignment'):
        subset = [r for r in results if r['group'] == group]
        groups[group] = dict(files=len(subset), passed=sum(r['passed'] for r in subset),
                             rules=sum(len(r.get('checks', [])) for r in subset),
                             readings=sum(c['readings'] for r in subset for c in r.get('checks', [])))
    report = dict(schema_version=1, suite=manifest['suite'],
                  run_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  platform=platform.platform(), machine=platform.machine(),
                  execution_arch=next(iter(identities))[2] if len(identities) == 1 else None,
                  binary_sha256=binary_hash, manifest_sha256=sha256(manifest_path),
                  runner_sha256=sha256(pathlib.Path(__file__)), archive_verified=archive_verified,
                  groups=groups, results=results, errors=failures,
                  passed=all(r['passed'] for r in results) and not failures)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps(groups), flush=True)
    for failure in failures:
        print('FAIL', failure, flush=True)
    print('EBU qualification:', 'PASS' if report['passed'] else 'FAIL', flush=True)
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
