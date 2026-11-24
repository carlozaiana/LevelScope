"""Run only explicitly selected M20/UI groups from an immutable candidate."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time


def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def main():
    region_groups = ('core', 'stereo', 'surround', 'heterogeneous', 'failures', 'population', 'audio', 'gating', 'cinema')
    other_groups = ('monitor', 'control-variants', 'ui', 'event-ui', 'inner-ui', 'applied-ui', 'operation-ui', 'playback-ui', 'graph-ui', 'gain-input-ui', 'current-sections-ui', 'applied-sections', 'processing-smoke', 'graph-boundaries', 'gain-entry', 'devices', 'preparation', 'timeline', 'requirements', 'section-smoke', 'legacy-core', 'legacy-stereo', 'legacy-failures', 'legacy-gating')
    event_groups = ('event-stereo', 'event-surround', 'event-fallback', 'inner-stereo', 'inner-surround')
    parser = argparse.ArgumentParser()
    parser.add_argument('candidate', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--case', action='append', required=True, choices=region_groups + other_groups + event_groups)
    parser.add_argument('--gating-input', type=Path)
    parser.add_argument('--cinema-input', type=Path)
    parser.add_argument('--inner-session', type=Path)
    args = parser.parse_args()
    candidate, destination = args.candidate.resolve(), args.destination.resolve()
    assert not destination.exists(), 'Use fresh receipt paths'
    assert len(set(args.case)) == len(args.case), 'Do not silently duplicate requested cases'
    identity = json.loads((candidate / 'identity.json').read_text())
    assert sha(candidate / 'identity.json') == (candidate / 'identity.sha256').read_text().split()[0], 'Frozen identity checksum'
    inputs = {}
    for name, path in (('gating', args.gating_input), ('cinema', args.cinema_input), ('inner-ui', args.inner_session), ('applied-ui', args.inner_session), ('operation-ui', args.inner_session), ('graph-ui', args.inner_session)):
        if name in args.case:
            assert path and path.is_file(), 'Explicit bounded input required for ' + name
            inputs[name] = {'path': str(path.resolve()), 'sha256': sha(path)}

    def verify():
        for name, digest in identity['artifacts_sha256'].items():
            assert sha(candidate / name) == digest, name

    verify()
    destination.mkdir(parents=True)
    receipt = {'passed': False, 'os': platform.mac_ver()[0], 'architecture': platform.machine(),
               'source_commit': identity['source_commit'], 'build': identity['measurement_build'],
               'candidate_identity_sha256': sha(candidate / 'identity.json'), 'requested_cases': args.case,
               'physical_devices_opened': False, 'combined_qualification_complete': False, 'inputs': inputs, 'groups': []}
    try:
        for case in args.case:
            root, output = destination / (case + '-files'), destination / (case + '.json')
            if case in event_groups:
                target = 'LevelScopeEventPlanningTests'
                arguments = ['--group', case if case.startswith('inner-') else case[6:], '--root', str(root), '--output', str(output)]
            elif case in region_groups:
                target = 'LevelScopeRegionPlanningTests'
                arguments = ['--group', case, '--root', str(root), '--output', str(output)]
                if case in ('gating', 'cinema'):
                    supplied = args.gating_input if case == 'gating' else args.cinema_input
                    assert supplied and supplied.is_file(), 'Explicit bounded input required for ' + case
                    arguments += ['--input', str(supplied.resolve())]
            elif case == 'monitor':
                target = 'LevelScopeAuditionTests'
                arguments = ['--monitor-only', '--report', str(output)]
            elif case.startswith('legacy-'):
                target = 'LevelScopeLraProposalTests'
                arguments = ['--group', case[7:], '--root', str(root), '--output', str(output)]
            elif case in ('gain-input-ui', 'current-sections-ui', 'ui', 'event-ui', 'inner-ui', 'applied-ui', 'operation-ui', 'playback-ui', 'graph-ui', 'control-variants', 'devices', 'preparation'):
                target = {'gain-input-ui': 'LevelScopeRegionUITests', 'current-sections-ui': 'LevelScopeRegionUITests', 'graph-ui': 'LevelScopeRegionUITests', 'playback-ui': 'LevelScopeRegionUITests', 'operation-ui': 'LevelScopeRegionUITests', 'applied-ui': 'LevelScopeRegionUITests', 'inner-ui': 'LevelScopeRegionUITests', 'event-ui': 'LevelScopeRegionUITests', 'control-variants': 'LevelScopeRegionUITests', 'ui': 'LevelScopeRegionUITests', 'devices': 'LevelScopeDeviceSettingsTests',
                          'preparation': 'LevelScopePreparationFeedbackTests'}[case]
                arguments, output = [str(root)], root / 'result.json'
                if case == 'control-variants':
                    arguments += ['--variants']
                elif case in ('event-ui', 'playback-ui', 'gain-input-ui', 'current-sections-ui'):
                    arguments += ['--'+case]
                elif case in ('inner-ui', 'applied-ui', 'operation-ui', 'graph-ui'):
                    arguments += ['--' + case, str(args.inner_session.resolve())]
            elif case in ('applied-sections', 'processing-smoke'):
                target = 'LevelScopeSectionProcessingTests'
                arguments = ['--applied-sections' if case == 'applied-sections' else '--smoke', '--root', str(root), '--output', str(output)]
            else:
                target = {'gain-entry': 'LevelScopeTimelineProgressTests', 'graph-boundaries': 'LevelScopeTimelineProgressTests', 'timeline': 'LevelScopeTimelineProgressTests', 'requirements': 'LevelScopeJobRequirementsTests',
                          'section-smoke': 'LevelScopeSectionProposalTests'}[case]
                arguments = ['--root', str(root), '--output', str(output)]
                if case in ('graph-boundaries', 'gain-entry'):
                    arguments += ['--'+case]
                if case == 'section-smoke':
                    arguments += ['--smoke']
            command = [str(candidate / 'bin' / target), *arguments]
            started = time.monotonic()
            with (destination / (case + '.log')).open('w') as log:
                outcome = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=900 if case in ('inner-ui', 'applied-ui', 'operation-ui', 'graph-ui') else 300)
            data = json.loads(output.read_text()) if output.exists() else {}
            checks = data.get('checks', [])
            passed = outcome.returncode == 0 and data.get('passed') and data.get('build') == identity['measurement_build'] and bool(checks) and all(c['passed'] for c in checks)
            receipt['groups'].append({'name': case, 'command': command, 'returncode': outcome.returncode,
                                      'seconds': time.monotonic() - started, 'passed': bool(passed),
                                      'checks': len(checks), 'receipt': str(output.relative_to(destination)),
                                      'receipt_sha256': sha(output) if output.exists() else None})
            print(case, 'PASS' if passed else 'FAIL', flush=True)
            assert passed, case
        verify()
        for info in inputs.values():
            assert sha(Path(info['path'])) == info['sha256'], 'Original bounded fixture changed'
        receipt['passed'] = True
    finally:
        (destination / 'summary.json').write_text(json.dumps(receipt, indent=2) + '\n')


if __name__ == '__main__':
    main()
