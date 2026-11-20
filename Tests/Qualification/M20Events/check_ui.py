"""Check the finite Sequoia integration evidence, not general M20 readiness."""
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
    parser.add_argument('gate', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    assert not args.output.exists(), 'Use a new evidence record'
    identity = json.loads((args.candidate / 'identity.json').read_text())
    assert sha(args.candidate / 'identity.json') == (args.candidate / 'identity.sha256').read_text().split()[0]
    for name, digest in identity['artifacts_sha256'].items():
        assert sha(args.candidate / name) == digest, name
    gate = json.loads((args.gate / 'identity.json').read_text())
    assert not gate['product_alternative_selection_enabled']
    for name, digest in gate['artifacts_sha256'].items():
        assert sha(args.gate / name) == digest, 'Authority gate ' + name
    authority = json.loads((args.gate / 'result.json').read_text())
    assert authority['passed'] and authority['build'] == gate['measurement_build']
    assert authority['checks'] and all(c['passed'] for c in authority['checks'])
    run = json.loads((args.run / 'summary.json').read_text())
    expected = {'event-ui', 'ui', 'timeline', 'event-stereo', 'event-surround', 'event-fallback'}
    assert run['passed'] and set(run['requested_cases']) == expected
    assert run['source_commit'] == identity['source_commit'] and run['build'] == identity['measurement_build']
    assert run['candidate_identity_sha256'] == sha(args.candidate / 'identity.json')
    assert run['os'].startswith('15.') and run['architecture'] == 'arm64'
    assert not run['physical_devices_opened'] and not run['combined_qualification_complete']
    assert len(run['groups']) == len(expected) and {g['name'] for g in run['groups']} == expected
    for group in run['groups']:
        path = args.run / group['receipt']
        assert group['passed'] and group['returncode'] == 0 and sha(path) == group['receipt_sha256']
        data = json.loads(path.read_text())
        assert data['passed'] and data['build'] == identity['measurement_build']
        assert len(data['checks']) == group['checks'] and all(c['passed'] for c in data['checks'])
    result = {
        'scope': 'Complete-event native UI integration on Sequoia only',
        'source_commit': identity['source_commit'],
        'measurement_build': identity['measurement_build'],
        'candidate_identity_sha256': sha(args.candidate / 'identity.json'),
        'gate_identity_sha256': sha(args.gate / 'identity.json'),
        'gate_checks': len(authority['checks']),
        'groups': run['groups'],
        'focused_checks': sum(g['checks'] for g in run['groups']),
        'passed': True, 'combined_m20_qualification_complete': False,
        'physical_devices_opened': False,
        'pending': ['actual High Sierra and new-method cross-platform execution',
                    'broader detector controls and independent holdouts',
                    'coordinated listening and musical quality review'],
    }
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('PASS:', result['focused_checks'], 'focused Sequoia assertions; combined qualification remains pending')


if __name__ == '__main__':
    main()
