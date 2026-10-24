#!/usr/bin/env python3
"""Resume only outstanding M6 candidate-01 physical checks; all routes muted."""
import argparse, hashlib, json, pathlib, platform, subprocess, sys, time

p = argparse.ArgumentParser()
for name in ('candidate', 'fixtures', 'actual', 'vfr', 'device', 'sample-rate', 'output'):
    p.add_argument('--' + name, required=True)
p.add_argument('--include-prores-matrices', action='store_true')
a = p.parse_args()
candidate = pathlib.Path(a.candidate).resolve()
fixtures = pathlib.Path(a.fixtures).resolve()
out = pathlib.Path(a.output).resolve()
out.mkdir(parents=True, exist_ok=False)
binary = candidate / 'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
def sha(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()
expected = json.loads((candidate / 'identity.json').read_text())
assert sha(binary) == expected['artifacts_sha256'][str(binary.relative_to(candidate))]
identity = dict(source_digest=expected['source_digest'], binary_sha256=sha(binary),
                architecture=platform.machine(), os=subprocess.check_output(['sw_vers'], text=True),
                device=a.device, sample_rate=int(a.sample_rate), all_routes_muted=True,
                monitor_gain_db=-30, helper_sha256=sha(__file__),
                vfr_sha256=sha(a.vfr), fixture_manifest_sha256=sha(fixtures / 'manifest.json'))
(out / 'identity.json').write_text(json.dumps(identity, indent=2) + '\n')
runs = []
def run(name, cmd, summary):
    print('START', name, flush=True)
    start = time.monotonic()
    with (out / (name + '.log')).open('x') as log:
        rc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT).returncode
    data = json.loads(summary.read_text()) if summary.exists() else {}
    passed = rc == 0 and data.get('passed') is True
    runs.append(dict(name=name, exit_code=rc, passed=passed,
                     elapsed_seconds=time.monotonic()-start, argv=cmd))
    (out / 'runs.json').write_text(json.dumps(runs, indent=2) + '\n')
    print('FINISH', name, 'PASS' if passed else 'FAIL', flush=True)
    if not passed:
        (out / 'summary.json').write_text(json.dumps(dict(passed=False, identity=identity, runs=runs), indent=2) + '\n')
        raise SystemExit(1)
def prores(domain, phase):
    name = 'prores-' + phase + '-' + domain
    result = out / name
    cmd = [sys.executable, '-B', str(candidate / 'Tests/run_prores_playback_qualification.py'),
           '--binary', str(binary), '--fixtures', str(fixtures), '--actual', a.actual,
           '--output', str(result), '--device', a.device, '--sample-rate', a.sample_rate,
           '--phase', phase, '--audition-state', domain]
    if domain == 'current':
        cmd += ['--processing-gain', '-3.305074118944']
    run(name, cmd, result / 'summary.json')
if a.include_prores_matrices:
    for domain in ('source', 'current'):
        prores(domain, 'matrix')
for domain in ('source', 'current'):
    for kind, movie in [('cfr', candidate / 'video-fixtures/half-hour.mov'), ('vfr', pathlib.Path(a.vfr))]:
        name = 'h264-' + kind + '-' + domain
        report = out / (name + '.json')
        cmd = [str(binary), '--allow-audio', '--muted', '--device', a.device,
               '--audio', str(fixtures / ('selected-clock-' + a.sample_rate + '.wav')),
               '--video', str(movie), '--seconds', '18', '--report', str(report), '--audition-state', domain]
        if domain == 'current':
            cmd += ['--processing-gain', '-3.305074118944']
        run(name, cmd, report)
prores('current', 'endurance')
(out / 'summary.json').write_text(json.dumps(dict(passed=True, identity=identity, runs=runs), indent=2) + '\n')
