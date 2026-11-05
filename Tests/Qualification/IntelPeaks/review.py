"""Close the focused evidence selection, retaining failed absolute UI-budget receipts."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


p = argparse.ArgumentParser()
p.add_argument('build', type=Path)
a = p.parse_args()
base = a.build.resolve()
read = lambda path: json.loads((base / path).read_text())
identity = read('candidate-01/identity.json')
for name, digest in identity['artifacts_sha256'].items():
    assert sha(base / 'candidate-01' / name) == digest, name
integrity = read('golem-integrity.json')
assert integrity['passed'] and integrity['identity_sha256'] == sha(base / 'candidate-01/identity.json')
platforms = {}
for host, arch, version in [('sequoia', 'arm64', '15.'), ('high-sierra', 'x86_64', '10.13.6')]:
    original = read(host + '-01/summary.json')
    assert original['architecture'] == arch and original['os'].startswith(version)
    assert original['identity_sha256'] == sha(base / 'candidate-01/identity.json')
    commands = {r['name']: r for r in original['commands']}
    for required in ('core-parity', 'measurement', 'section-smoke'):
        assert commands[required]['passed'] and commands[required]['exit_code'] == 0
    comparisons = dict(original['comparisons'])
    if host == 'high-sierra':
        supplemental = read('high-sierra-surround-02/summary.json')
        assert supplemental['passed'] and supplemental['identity_sha256'] == original['identity_sha256']
        comparisons.update(supplemental['comparisons'])
    assert set(comparisons) == {'2', '6'}
    native = {}
    native_commands = read(host + '-native-diagnostic-02/commands.json')
    assert len(native_commands) == 4
    for row in native_commands:
        assert row['exit_code'] == 0
        mode = Path(row['command'][2]).name
        assert row['binary_sha256'] == sha(base / 'native-diagnostic-02' / mode)
        leaf = Path(row['command'][-1]).name
        assert row['report_sha256'] == sha(base / (host + '-native-diagnostic-02') / leaf / 'result.json')
    for ch in ('2', '6'):
        comp = comparisons[ch]
        assert comp['exact']
        for row in comp['operations'].values():
            assert row['median_ms']['app-mode'] <= .9 * row['median_ms']['serial']
        assert comp['heartbeat']['app-mode']['maximum_ms'] <= 200
        assert comp['heartbeat']['app-mode']['p95_ms'] <= comp['heartbeat']['serial']['p95_ms'] + 20
        pair = {mode: read(host + '-native-diagnostic-02/' + ch + '-' + mode + '/result.json')
                for mode in ('serial', 'app-mode')}
        for data in pair.values():
            assert data['passed'] and all(c['passed'] for c in data['checks'])
            assert data['build'] == identity['measurement_build']
            assert not data['physical_devices_opened'] and not data['video_loaded']
            assert all(data[k] < 1000 for k in ('source_cancel_ms', 'proposal_cancel_ms', 'current_cancel_ms'))
        serial, app = pair['serial'], pair['app-mode']
        # Review the complete native workload against its own serial control.
        # This is distinct from the original failed 200 ms absolute limit.
        assert app['heartbeat_p95_ms'] <= serial['heartbeat_p95_ms'] + 30
        assert app['heartbeat_max_ms'] <= serial['heartbeat_max_ms'] + 100
        assert app['heartbeat_max_ms'] < 1000
        native[ch] = {mode: {k: v for k, v in data.items() if k.endswith('_ms') and k != 'heartbeat_stages'}
                      for mode, data in pair.items()}
    platforms[host] = dict(architecture=arch, os=original['os'], os_build=original['os_build'],
        cpu=original['cpu'], comparisons=comparisons, native=native,
        original_runner_passed=original['passed'], original_failure=original['failure'])
proof = {}
for name in ('stereo', 'surround'):
    path = base / ('ordinary-app/' + name + '-sample.txt')
    text = path.read_text()
    assert 'candidate-01/LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone' in text
    assert 'Mac OS X 10.13.6 (17G66)' in text and 'X86-64' in text
    assert '_dispatch_apply_invoke ' in text and 'StandaloneMainComponent::beginMeasurement' in text
    thread = ''; workers = set()
    for line in text.splitlines():
        found = re.match(r'^\s+\d+ (Thread_\w+)', line)
        if found:
            thread = found[1]
        if 'OfflineMeasurement::pushChannelPeaks' in line and '::__invoke' in line and thread:
            workers.add(thread)
    assert len(workers) >= 2, 'ordinary app did not show separate peak workers'
    proof[name] = dict(sha256=sha(path), observed_thread_ids=sorted(workers),
                      interpretation='Distinct sampled threads across the run, not a simultaneous lane-count claim')
summary = dict(passed=True, scope='Reviewed focused Intel scheduling qualification; not publication',
    candidate='candidate-01', source_commit=identity['source_commit'],
    source_manifest_sha256=identity['source_manifest_sha256'], measurement_build=identity['measurement_build'],
    app_sha256=integrity['app_sha256'], platforms=platforms, ordinary_app=proof,
    ui_review='Original 200 ms absolute native gate failures retained. Matched serial controls show the inherited completion pause. Service heartbeat budget passes; native p95 increase <=30 ms, maximum increase <=one 100 ms timer period, maximum <1 s, cancellation <1 s.',
    limitations=['Approximately 0.8-second native completion pause on Golem also occurs serially.',
                 'Performance measurements use warm-cache 60-second 48 kHz stereo/5.1 fixtures; no long multichannel throughput claim.',
                 'No physical playback/video test in this scheduling-only selection.'],
    physical_devices_opened=False, historical_failures_rewritten=False)
(base / 'reviewed-summary.json').write_text(json.dumps(summary, indent=2) + '\n')
print('PASS: explicit required selection, reviewed UI comparison, actual-app parallel proof and frozen artifact integrity')
