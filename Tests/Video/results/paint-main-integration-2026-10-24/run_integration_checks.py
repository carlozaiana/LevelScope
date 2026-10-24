#!/usr/bin/env python3
"""Device-closed checks for the exact diagnostic integration; no playback."""
import datetime
import hashlib
import json
import os
import pathlib
import platform
import subprocess
import sys

ROOT = pathlib.Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation')
BUILD = ROOT / 'build/paint-main-integration-20261024-01'
CMAKE = BUILD / 'cmake'
OUT = BUILD / 'checks'
OUT.mkdir()
ENV = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
summary = {'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
           'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
           'native_machine': platform.machine(), 'physical_playback': False,
           'audio_devices_opened': False, 'checks': [], 'all_passed': False}

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save():
    (OUT / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')

def run(name, argv, expected=0, timeout=600):
    print('START', name, flush=True)
    with (OUT / (name + '.log')).open('x') as log:
        result = subprocess.run([str(x) for x in argv], cwd=ROOT, env=ENV,
                                stdout=log, stderr=subprocess.STDOUT, timeout=timeout)
    summary['checks'].append({'name': name, 'argv': [str(x) for x in argv],
                              'expected_exit_code': expected, 'exit_code': result.returncode})
    save()
    if result.returncode != expected:
        print((OUT / (name + '.log')).read_text()[-4000:], flush=True)
        raise RuntimeError(name + ' returned an unexpected exit code')
    print('PASS', name, flush=True)

try:
    candidate = json.loads((ROOT / 'build/paint-diagnostic-20261024-01/candidate-03/identity.json').read_text())
    for path, digest in candidate['source_files_sha256'].items():
        assert sha(ROOT / path) == digest, path
    summary['tested_source_files_identical'] = len(candidate['source_files_sha256'])
    summary['reused_candidate_source_digest'] = candidate['source_digest']
    listed = json.loads(subprocess.check_output(['ctest', '--test-dir', str(CMAKE), '--show-only=json-v1'], text=True))
    expected_tests = {'video-timing', 'offline-measurement', 'measured-compliance',
                      'global-trim-proposals', 'measured-diagnostics', 'programme-sections',
                      'verified-audition', 'standalone-scaffold'}
    summary['ctest_names'] = [t['name'] for t in listed['tests']]
    assert set(summary['ctest_names']) == expected_tests
    run('ctest', ['ctest', '--test-dir', CMAKE, '--output-on-failure', '--parallel', '1'])
    run('trace-buffer', [CMAKE / 'LevelScopeVideoPaintTraceTests'])
    run('runner-guards', [sys.executable, '-B', ROOT / 'Tests/VideoPaint/RunnerTests.py'])
    apps = {'LevelScopeStandalone': 'LevelScope Standalone',
            'LevelScopeVideoQualification': 'LevelScope Video Qualification',
            'LevelScopeVideoPaintDiagnostic': 'LevelScope Video Paint Diagnostic'}
    executables = {}
    for target, name in apps.items():
        executables[target] = CMAKE / (target + '_artefacts') / 'Release' / (name + '.app') / 'Contents/MacOS' / name
    diag = executables['LevelScopeVideoPaintDiagnostic']
    trace = OUT / 'offline-trace.json'
    run('diagnostic-offline-startup', [diag, '--diagnostic-offline-check', trace])
    before = sha(trace)
    run('exclusive-output-guard', [diag, '--diagnostic-offline-check', trace], expected=1)
    assert sha(trace) == before
    data = json.loads(trace.read_text())
    assert data['schema'] == 'levelscope-paint-trace-v1'
    assert data['valid'] is False  # Deliberately no display observer in the offline export check.
    assert data['ui_dropped'] == 0 and data['display_dropped'] == 0
    summary['offline_trace_intentionally_not_qualification_evidence'] = True
    summary['instrumentation_target_isolation'] = {}
    for target, exe in executables.items():
        flags = (CMAKE / 'CMakeFiles' / (target + '.dir') / 'flags.make').read_text()
        expected = target == 'LevelScopeVideoPaintDiagnostic'
        assert ('-DLEVELSCOPE_VIDEO_PAINT_DIAGNOSTICS=1' in flags) == expected
        for arch in ('arm64', 'x86_64'):
            symbols = subprocess.check_output(['nm', '-arch', arch, str(exe)], text=True)
            assert ('paintdiagnostic' in symbols) == expected, (target, arch)
        summary['instrumentation_target_isolation'][target] = expected
    invocation = json.loads((BUILD / 'build-invocation.json').read_text())
    for target in invocation['named_targets']:
        if target in executables:
            continue
        candidates = [p for p in [CMAKE / target, CMAKE / (target + '_artefacts') / 'Release' / target] if p.is_file()]
        assert len(candidates) == 1, (target, candidates)
        executables[target] = candidates[0]
    identities = {}
    for target, exe in executables.items():
        archs = subprocess.check_output(['lipo', '-archs', str(exe)], text=True).strip().split()
        assert set(archs) == {'arm64', 'x86_64'}, target
        load = subprocess.check_output(['otool', '-arch', 'x86_64', '-l', str(exe)], text=True)
        lines = load.splitlines()
        versions = [lines[i + 2].strip() for i, line in enumerate(lines) if 'cmd LC_VERSION_MIN_MACOSX' in line]
        assert versions == ['version 10.13'], (target, versions)
        identities[target] = {'path': str(exe), 'sha256': sha(exe), 'architectures': archs,
                              'x86_64_deployment_target': '10.13'}
    (OUT / 'binary-identities.json').write_text(json.dumps(identities, indent=2) + '\n')
    summary['fresh_binary_count'] = len(identities)
    plugin_executables = [str(p) for bundle in CMAKE.rglob('*.vst3')
                          for p in bundle.glob('Contents/MacOS/*') if p.is_file()]
    assert not plugin_executables
    summary['vst3_executables_built'] = plugin_executables
    summary['vst3_configure_metadata_only'] = True
    summary['all_passed'] = True
except Exception as error:
    summary['error'] = str(error)
    raise
finally:
    summary['completed_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    save()
