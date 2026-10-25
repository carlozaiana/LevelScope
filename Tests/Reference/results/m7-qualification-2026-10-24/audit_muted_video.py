#!/usr/bin/env python3
"""Post-run M7 evidence audit; reads reports/media identities, never opens devices.

This audit is qualification documentation, not a frozen candidate build input.
It retains the existing ProRes memory/paint gates and adds M7 channel identity.
"""
import argparse
import ast
import hashlib
import json
import math
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--root', type=Path, required=True)
p.add_argument('--candidate', type=Path, required=True)
p.add_argument('--fixtures', type=Path, required=True)
p.add_argument('--actual', type=Path, required=True)
p.add_argument('--vfr', type=Path, required=True)
p.add_argument('--memory-verifier', type=Path, required=True)
p.add_argument('--sequoia-physical', type=Path,
               help='Explicit retained attempt to audit; original reports are never replaced')
p.add_argument('--high-sierra-physical', type=Path)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
if a.output.exists():
    raise SystemExit('Use a fresh audit output; previous evidence is preserved.')

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def read(path):
    return json.loads(path.read_text())

# Reuse the independently implemented, existing primary memory verifier without
# invoking its old candidate-04 command-line driver or altering its thresholds.
tree = ast.parse(a.memory_verifier.read_text())
function = next(n for n in tree.body if isinstance(n, ast.FunctionDef)
                and n.name == 'verify_memory')
namespace = {'math': math}
exec(compile(ast.Module(body=[function], type_ignores=[]), str(a.memory_verifier), 'exec'), namespace)
verify_memory = namespace['verify_memory']

identity = read(a.candidate / 'identity.json')
manifest = read(a.fixtures / 'manifest.json')
binary_key = 'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
binary_hash = identity['artifacts_sha256'][binary_key]
build = read(a.root / 'sequoia-02/multichannel.json')['build']
errors, rows = [], []

def require(condition, label):
    if not condition:
        errors.append(label)

def inspect(path, arch, device, rate, domain, movie_hash, duration, fps=None):
    if not path.exists():
        errors.append(str(path) + ': missing required report')
        return
    data = read(path)
    gain = 0.0 if domain == 'source' else -3.305074118944
    expected = dict(passed=True, architecture=arch, device=device, sample_rate=rate,
        audition_state='Source' if domain == 'source' else 'Current State',
        processing_gain_db=gain, measurement_build=build, video_sha256=movie_hash,
        audio_sha256=identity['artifacts_sha256'][f'm7-fixtures/clock-71-{rate}.flac'],
        requested_seconds=duration, all_routes_muted=True, selected_audio_channels=8,
        layout_capability='pcm-7.1-bs2051-system-i-v1',
        pre_monitor_nonzero_channel_mask=255, pre_monitor_channel_count_valid=True,
        actual_paint_skipped_frames=0, skipped_frame_limit=0, blank_paints=0,
        bound_failures=0, paint_deadline_failures=0, underrun_blocks=0,
        attachment_within_budget=True, memory_passed=True, hal_output_timestamp=True)
    failures = [key for key, value in expected.items() if data.get(key) != value]
    failures += verify_memory(data)
    if data.get('source_seconds', -1) < duration:
        failures.append('incomplete duration')
    if not data.get('audition_identity') or data.get('pre_monitor_frames_observed', 0) <= 0:
        failures.append('missing selected-audio identity/capture')
    if fps and data.get('unique_painted_frames', 0) < math.floor((duration - 2) * fps) - 1:
        failures.append('painted coverage')
    if duration >= 120 and data.get('memory_growth_assessed') is not True:
        failures.append('missing memory growth assessment')
    errors.extend(str(path) + ': ' + item for item in failures)
    rows.append(dict(report=str(path.relative_to(a.root)), sha256=sha(path),
                     passed=not failures, duration=duration,
                     maximum_frame_pts_skew_seconds=data.get('maximum_frame_pts_skew_seconds'),
                     actual_paint_skipped_frames=data.get('actual_paint_skipped_frames'),
                     peak_rss_bytes=data.get('peak_sampled_playback_rss_bytes')))

for platform, arch, device, rate in (
    ('sequoia', 'arm64', 'Mac Studio Speakers', 48000),
    ('high-sierra', 'x86_64', 'Built-in Output', 44100),
):
    root = (a.sequoia_physical if platform == 'sequoia' else a.high_sierra_physical) \
        or a.root / (platform + '-muted-video-02')
    for phase, domains in (('matrix', ('source', 'current')), ('endurance', ('current',))):
        for domain in domains:
            directory = root / ('prores-' + phase + '-' + domain)
            try:
                run_identity, summary = read(directory / 'identity.json'), read(directory / 'summary.json')
                expected_identity = dict(binary_sha256=binary_hash, architecture=arch,
                    fixture_manifest_sha256=sha(a.fixtures / 'manifest.json'),
                    actual_sha256=sha(a.actual), device=device, sample_rate=rate,
                    all_routes_muted=True, audition_state=domain,
                    selected_audio_sha256=identity['artifacts_sha256'][f'm7-fixtures/clock-71-{rate}.flac'],
                    processing_gain_db=None if domain == 'source' else -3.305074118944)
                for key, value in expected_identity.items():
                    require(run_identity.get(key) == value, str(directory) + ': identity ' + key)
                require(summary.get('passed') is True, str(directory) + ': summary failed')
            except (OSError, ValueError) as error:
                errors.append(str(directory) + ': ' + str(error))
                continue
            cases = [(name, 18) for name, row in manifest.items() if row.get('required')]
            if phase == 'endurance':
                cases = [(profile + '-1800.mov', 1800) for profile in ('proxy', 'lt', '422')]
            expected_names = {Path(name).stem for name, _ in cases}
            if phase == 'matrix':
                expected_names.add(a.actual.stem)
            actual_names = {path.stem for path in directory.glob('*.json')
                            if path.stem not in ('identity', 'runs', 'summary')}
            require(actual_names == expected_names, str(directory) + ': required case set')
            require(summary.get('test_count') == len(expected_names), str(directory) + ': case count')
            for name, duration in cases:
                fixture = manifest[name]
                inspect(directory / (Path(name).stem + '.json'), arch, device, rate,
                        domain, fixture['sha256'], duration, fixture['scale'] / fixture['step'])
            if phase == 'matrix':
                inspect(directory / (a.actual.stem + '.json'), arch, device, rate,
                        domain, sha(a.actual), 100, 25)
    for domain in ('source', 'current'):
        inspect(root / ('h264-cfr-' + domain + '.json'), arch, device, rate, domain,
                identity['artifacts_sha256']['video-fixtures/half-hour.mov'], 18)
        inspect(root / ('h264-vfr-' + domain + '.json'), arch, device, rate, domain,
                sha(a.vfr), 18)
    try:
        require(read(root / 'summary.json').get('passed') is True, str(root) + ': overall run failed')
    except (OSError, ValueError) as error:
        errors.append(str(root) + ': ' + str(error))

result = dict(passed=not errors, source_digest=identity['source_digest'],
    binary_sha256=binary_hash, memory_verifier_sha256=sha(a.memory_verifier),
    sequoia_physical=str(a.sequoia_physical or a.root / 'sequoia-muted-video-02'),
    high_sierra_physical=str(a.high_sierra_physical or a.root / 'high-sierra-muted-video-02'),
    scope='Muted physical audio-clock/UI-paint evidence only; not optical/acoustic sync or physical 7.1 listening',
    errors=errors, reports=rows)
a.output.write_text(json.dumps(result, indent=2) + '\n')
print('M7 muted video audit:', 'PASS' if result['passed'] else 'FAIL', len(rows), 'reports', len(errors), 'errors')
raise SystemExit(0 if result['passed'] else 1)
