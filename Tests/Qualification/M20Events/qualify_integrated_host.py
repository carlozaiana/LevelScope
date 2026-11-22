"""Finite silent Golem run from the immutable integrated candidate (Python 3.6+)."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

EXPECTED_IDENTITY = '8f1497d12d4d5a1c3a0e8cd4aa3e864cc70a540023a399041d5dbd808ac93e79'

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('candidate', type=Path)
    parser.add_argument('session', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--group', action='append', choices=['inner-ui', 'event-ui', 'ui', 'inner-stereo', 'inner-surround', 'inner-population'], required=True)
    parser.add_argument('--expected-identity', default=EXPECTED_IDENTITY, help='Explicit preserved replacement package identity when qualifying a documented fix')
    args = parser.parse_args()
    candidate, root = args.candidate.resolve(), args.output.resolve()
    assert platform.mac_ver()[0] == '10.13.6' and platform.machine() == 'x86_64', 'Actual Intel High Sierra required'
    assert not root.exists() and len(args.group) == len(set(args.group)), 'Fresh destinations and explicit finite groups'
    identity = json.loads((candidate/'identity.json').read_text())
    assert sha(candidate/'identity.json') == args.expected_identity
    def verify():
        for name, digest in identity['artifacts_sha256'].items():
            assert sha(candidate/name) == digest, name
    verify()
    root.mkdir(parents=True)
    result = dict(passed=False, os=platform.mac_ver()[0], architecture=platform.machine(), host=platform.node(),
                  candidate_identity_sha256=args.expected_identity, build=identity['measurement_build'],
                  source_commit=identity['source_commit'], planned=args.group, groups=[],
                  physical_devices_opened=False, video_or_endurance_run=False)
    try:
        for group in args.group:
            files = root/(group+'-files')
            receipt = root/(group+'.json')
            if group in ('inner-ui', 'event-ui', 'ui'):
                command = [str(candidate/'bin/LevelScopeRegionUITests'), str(files)]
                if group != 'ui': command.append('--'+group)
                if group == 'inner-ui': command.append(str(args.session.resolve()))
                receipt = files/'result.json'
            else:
                command = [str(candidate/'bin/LevelScopeEventPlanningTests'), '--group', group, '--root', str(files), '--output', str(receipt)]
            print('START '+group, flush=True)
            start = time.monotonic()
            with (root/(group+'.log')).open('w') as log:
                process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=900)
            data = json.loads(receipt.read_text()) if receipt.exists() else {}
            passed = process.returncode == 0 and data.get('passed') and data.get('checks') and all(c['passed'] for c in data['checks']) and data.get('build') == identity['measurement_build']
            if passed and group == 'inner-population':
                row = next(c for c in data['cases'] if '64_proxy_seconds' in c)
                passed = row['64_proxy_seconds'] <= 20 and row['actual_peak_heap_bytes'] <= 64*1024*1024 and row['maximum_feature_reservation'] <= 10*1024*1024
            result['groups'].append(dict(name=group, command=command, seconds=time.monotonic()-start,
                returncode=process.returncode, passed=bool(passed), checks=len(data.get('checks', [])),
                receipt=str(receipt.relative_to(root)), receipt_sha256=sha(receipt) if receipt.exists() else None))
            (root/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
            print(('PASS ' if passed else 'FAIL ')+group, flush=True)
            if not passed: raise RuntimeError('Failed focused group: '+group)
        verify()
        result['passed'] = True
    finally:
        (root/'summary.json').write_text(json.dumps(result, indent=2)+'\n')

if __name__ == '__main__':
    main()
