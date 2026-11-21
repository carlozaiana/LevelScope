"""Finite competing-boundary native Sequoia evidence; not combined qualification."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('candidate', type=Path)
    parser.add_argument('run', type=Path)
    parser.add_argument('logic', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists(), 'Use a fresh evidence record'
    identity = json.loads((args.candidate / 'identity.json').read_text())
    assert sha(args.candidate / 'identity.json') == (args.candidate / 'identity.sha256').read_text().split()[0]
    for name, digest in identity['artifacts_sha256'].items():
        assert sha(args.candidate / name) == digest, name
    logic = json.loads((args.logic / 'identity.json').read_text())
    assert sha(args.logic / 'identity.json') == 'da714861bab2c40fd05ab350feead3a06b84b29884993298519b1ada2925c5ed'
    inherited = {name: digest for name, digest in logic['source_files_sha256'].items()
                 if name.startswith('Source/Core/') or name.startswith('Source/Standalone/SessionDocument')
                 or name.startswith('Source/Standalone/StandaloneSectionProposal')
                 or name.startswith('Source/Standalone/StandaloneLraProposal')
                 or name.startswith('Source/Standalone/StandaloneEventRegionProposal')
                 or name.startswith('Source/Standalone/RegionAssessmentReceipt')
                 or name.startswith('Source/Standalone/StandaloneRegionProposal')}
    assert inherited and any('InnerEventTracker' in name for name in inherited)
    for name, digest in inherited.items():
        assert identity['source_sha256'].get(name) == digest, 'Shared Logic changed: ' + name
    run = json.loads((args.run / 'summary.json').read_text())
    expected = {'inner-ui', 'event-ui', 'ui', 'timeline', 'inner-stereo', 'inner-surround'}
    assert run['passed'] and set(run['requested_cases']) == expected
    assert run['source_commit'] == identity['source_commit'] and run['build'] == identity['measurement_build']
    assert run['candidate_identity_sha256'] == sha(args.candidate / 'identity.json')
    assert run['os'].startswith('15.') and run['architecture'] == 'arm64'
    assert not run['physical_devices_opened'] and not run['combined_qualification_complete']
    assert len(run['groups']) == len(expected)
    for group in run['groups']:
        path = args.run / group['receipt']
        assert group['passed'] and group['returncode'] == 0 and sha(path) == group['receipt_sha256']
        data = json.loads(path.read_text())
        assert data['passed'] and data['build'] == identity['measurement_build']
        assert len(data['checks']) == group['checks'] and all(c['passed'] for c in data['checks'])
    result = {
        'scope': 'Competing-boundary native UI integration on Sequoia only',
        'source_commit': identity['source_commit'], 'measurement_build': identity['measurement_build'],
        'candidate_identity_sha256': sha(args.candidate / 'identity.json'),
        'logic_checkpoint': identity['logic_handover'], 'unchanged_shared_logic_files': inherited,
        'groups': run['groups'], 'focused_checks': sum(g['checks'] for g in run['groups']),
        'passed': True, 'combined_qualification_complete': False, 'physical_devices_opened': False,
        'pending': ['actual High Sierra / new-method cross-platform execution',
                    'independent programme validation', 'coordinated listening',
                    'separately coordinated broad video/endurance stage', 'authorized main integration'],
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('PASS:', result['focused_checks'], 'focused Sequoia assertions; combined qualification remains pending')


if __name__ == '__main__':
    main()
