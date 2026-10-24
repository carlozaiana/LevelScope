#!/usr/bin/env python3
"""Verify/preview by default. Physical playback requires fresh coordination and --run.
No settings writes, retries, default-device opens, or changes to frozen packages.
"""
import argparse
import datetime
import hashlib
import json
import pathlib
import plistlib
import platform
import subprocess
import sys
import time

PROFILES = {
    'sequoia': {'device': 'Mac Studio Speakers', 'rate': 48000, 'arch': 'arm64', 'os_prefix': '15.'},
    'high-sierra': {'device': 'Built-in Output', 'rate': 44100, 'arch': 'x86_64', 'os_prefix': '10.13.'},
}
ORDER = ['control', 'diagnostic', 'control', 'diagnostic', 'control', 'diagnostic']
def sha256(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def write_new(path, value):
    with pathlib.Path(path).open('x') as f: json.dump(value, f, indent=2, default=str); f.write('\n')

def verify_package(package):
    identity = json.loads((package / 'identity.json').read_text())
    for rel, expected in identity['artifacts_sha256'].items():
        path = package / rel
        if not path.is_file() or sha256(path) != expected:
            raise RuntimeError('Package integrity mismatch: ' + rel)
    return identity

def audio_snapshot():
    raw = subprocess.check_output(['/usr/sbin/system_profiler', '-xml', 'SPAudioDataType'], timeout=45)
    # Exclude profiler execution timestamps; retain the complete device configuration.
    return [entry["_items"] for entry in plistlib.loads(raw) if "_items" in entry]

def find_devices(value, name):
    found = []
    if isinstance(value, dict):
        if value.get('_name') == name and 'coreaudio_device_srate' in value: found.append(value)
        for item in value.values(): found += find_devices(item, name)
    elif isinstance(value, list):
        for item in value: found += find_devices(item, name)
    return found

def check_device(snapshot, profile):
    matches = find_devices(snapshot, profile['device'])
    if len(matches) != 1: raise RuntimeError('Requested device missing or ambiguous')
    if float(matches[0]['coreaudio_device_srate']) != profile['rate']:
        raise RuntimeError('Device rate changed; will not change it or launch playback')
    if int(matches[0].get('coreaudio_device_output', 0)) != 2:
        raise RuntimeError('Expected exactly two output channels')

def command(package, identity, profile, role, output):
    exe = package / identity['executables'][role]
    args = [str(exe), '--allow-audio', '--muted', '--device', profile['device'], '--seconds', '100',
        '--audio', str(package / 'fixtures' / ('selected-clock-%d.wav' % profile['rate'])),
        '--video', str(package / 'fixtures/51_SyncTest.mov'), '--report', str(output / 'report.json'),
        '--audition-state', 'source']
    if role == 'diagnostic': args += ['--trace', str(output / 'trace.json')]
    return args

def trace_summary(path):
    data = json.loads(path.read_text())
    if data.get('schema') != 'levelscope-paint-trace-v1' or not data.get('valid'):
        raise RuntimeError('Invalid or incomplete diagnostic trace')
    if data['ui_dropped'] or data['display_dropped']:
        raise RuntimeError('Trace overflow; diagnostic is invalid')
    ui, native = data['ui_events'], data['display_events']
    if not ui or not native: raise RuntimeError('Missing trace stream')
    for stream in [ui, native]:
        if any(len(e) != 14 or e[4] <= 0 for e in stream): raise RuntimeError('Invalid event schema')
        if any(a[0] > b[0] for a, b in zip(stream, stream[1:])): raise RuntimeError('Nonmonotonic trace')
    if ui[-1][1] != 10: raise RuntimeError('Trace has no closing event')
    summary = {'trace_sha256': sha256(path), 'ui_events': len(ui), 'display_events': len(native),
        'allocated_buffer_bytes': data['allocated_buffer_bytes'], 'nominal_hz': data['nominal_hz'],
        'display_id': data['display_id'], 'scope': 'Callback timing, not scanout; independent display-link observer is not JUCE internal callback'}
    for label, events in [('timer', [e for e in ui if e[1] == 1]), ('paint', [e for e in ui if e[1] == 5]),
                          ('main_vblank', [e for e in ui if e[1] == 7]), ('independent_display', native)]:
        gaps = [(b[0] - a[0]) / 1e6 for a, b in zip(events, events[1:])]
        summary[label + '_max_interval_ms'] = max(gaps, default=0)
    for label, begin, end in [('full_paint', 5, 6), ('full_timer', 1, 2)]:
        start = None; durations = []
        for event in ui:
            if event[1] == begin: start = event[0]
            if event[1] == end and start is not None: durations.append((event[0] - start) / 1e6); start = None
        summary[label + '_max_ms'] = max(durations, default=0)
    summary['window_transitions'] = [e for e in ui if e[1] == 9]
    return summary

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package', type=pathlib.Path, required=True)
    parser.add_argument('--platform', choices=PROFILES, required=True)
    parser.add_argument('--output', type=pathlib.Path)
    parser.add_argument('--run', action='store_true')
    parser.add_argument('--coordination-note', default='')
    args = parser.parse_args(argv)
    package = args.package.resolve(); profile = PROFILES[args.platform]
    identity = verify_package(package)
    version = subprocess.check_output(['/usr/bin/sw_vers', '-productVersion'], universal_newlines=True).strip()
    if platform.machine() != profile['arch'] or not version.startswith(profile['os_prefix']):
        raise RuntimeError('Host does not match the selected comparison platform')
    snapshot = audio_snapshot(); check_device(snapshot, profile)
    preview = {'mode': 'preview; no test launched', 'profile': profile, 'order': ORDER,
               'seconds_per_run': 100, 'identity_sha256': sha256(package / 'identity.json'),
               'source_digest': identity['source_digest'], 'audio_configuration': snapshot}
    if not args.run:
        print(json.dumps(preview, indent=2, default=str)); return 0
    if not args.coordination_note.strip() or args.output is None:
        raise RuntimeError('--run requires the fresh user coordination note and a new output directory')
    output = args.output.resolve()
    if package == output or package in output.parents:
        raise RuntimeError('Results must be outside the immutable package')
    output.mkdir(parents=True, exist_ok=False)
    preview['coordination_note'] = args.coordination_note
    preview['recorded_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    preview['display_snapshot'] = plistlib.loads(subprocess.check_output(['/usr/sbin/system_profiler', '-xml', 'SPDisplaysDataType'], timeout=45))
    preview['power_settings'] = subprocess.check_output(['/usr/bin/pmset', '-g', 'custom'], universal_newlines=True)
    preview['process_inventory'] = subprocess.check_output(['/bin/ps', '-axo', 'pid=,comm='], universal_newlines=True)
    write_new(output / 'preflight.json', preview)
    results = []
    for i, role in enumerate(ORDER, 1):
        current = audio_snapshot(); check_device(current, profile)
        folder = output / ('%02d-%s' % (i, role)); folder.mkdir()
        argv = command(package, identity, profile, role, folder)
        write_new(folder / 'invocation.json', {'argv': argv, 'audio_before': current, 'binary_sha256': sha256(argv[0]),
                 'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat()})
        print('Starting fixed run %d/6: %s, muted, 100 seconds' % (i, role), flush=True)
        started = time.monotonic()
        with (folder / 'process.log').open('x') as log:
            # subprocess.run timeout terminates only the child it started.
            completed = subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, timeout=240)
        report = json.loads((folder / 'report.json').read_text())
        if not report.get('all_routes_muted') or report['device'] != profile['device'] or report['sample_rate'] != profile['rate']:
            raise RuntimeError('Run disagrees with coordinated routing/device/rate')
        if report.get('audition_state') != 'Source' or report.get('processing_gain_db') != 0:
            raise RuntimeError('Unexpected selected-audio processing identity')
        for key, rel in [('audio_sha256', 'fixtures/selected-clock-%d.wav' % profile['rate']), ('video_sha256', 'fixtures/51_SyncTest.mov')]:
            if report[key] != identity['artifacts_sha256'][rel]: raise RuntimeError('Playback fixture identity mismatch')
        result = {'number': i, 'role': role, 'exit_code': completed.returncode, 'wall_seconds': time.monotonic() - started,
            'gate_passed': report['passed'], 'paint_skips': report['actual_paint_skipped_frames'],
            'report_sha256': sha256(folder / 'report.json')}
        if role == 'diagnostic':
            if not report.get('diagnostic_trace_valid'): raise RuntimeError('Diagnostic reported invalid tracing')
            result['trace'] = trace_summary(folder / 'trace.json')
        write_new(folder / 'outcome.json', result); results.append(result)
        # Continue the fixed comparison after a recorded timing failure, never retry it.
    after = audio_snapshot(); check_device(after, profile); verify_package(package)
    write_new(output / 'summary.json', {'scope': 'Fixed diagnostic comparison; no milestone acceptance or fix claim',
        'runs': results, 'all_existing_gates_passed': all(x['gate_passed'] for x in results),
        'audio_after': after, 'audio_configuration_unchanged': after == snapshot,
        'identity_sha256': sha256(package / 'identity.json')})
    return 0 if all(x['gate_passed'] for x in results) else 1

if __name__ == '__main__':
    try: sys.exit(main())
    except Exception as error: print('STOP: ' + str(error), file=sys.stderr); sys.exit(2)
