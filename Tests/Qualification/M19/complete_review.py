"""Read-only M19 completion audit; does not run tests, rebuild, or open devices.

Original numerical-only and frozen-package receipts retain their original status.
This separate decision requires explicit listening review and every selected group.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('root', type=Path)
parser.add_argument('--output', type=Path, required=True)
args = parser.parse_args()
root = args.root.resolve()
repo = Path(__file__).resolve().parents[3]
result = {'date': '2026-11-05', 'qualified': False, 'ready_for_checkpoint_review': False,
          'publication_authorized': False, 'checks': []}

def read(path):
    return json.loads(path.read_text())

def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def git(*arguments):
    return subprocess.check_output(['git', *arguments], cwd=repo, text=True).strip()

def check(ok, label):
    result['checks'].append({'name': label, 'passed': bool(ok)})
    if not ok:
        raise AssertionError(label)

def assertions_passed(report):
    checks = report.get('checks')
    # The inherited measurement runner records an assertion count, while the
    # newer runners retain named assertion records. Both also require passed.
    return (checks > 0 if isinstance(checks, int) else
            bool(checks) and all(v['passed'] for v in checks))

try:
    identities = {}
    for candidate in ['candidate-01', 'candidate-02']:
        folder = root / candidate
        identity = read(folder / 'identity.json')
        identities[candidate] = identity
        check(all((folder / name).is_file() and sha(folder / name) == digest
                  for name, digest in identity['artifacts_sha256'].items()),
              candidate + ': every frozen artifact unchanged')
        digest = hashlib.sha256(json.dumps(identity['source_files_sha256'],
                    sort_keys=True, separators=(',', ':')).encode()).hexdigest()
        check(digest == identity['source_manifest_sha256'], candidate + ': manifest digest')
    accepted = identities['candidate-02']
    check(all((repo / name).is_file() and sha(repo / name) == digest
              for name, digest in accepted['source_files_sha256'].items()),
          'Current checkout still matches every candidate-02 manifest file')
    changes = git('diff', '--name-only', identities['candidate-01']['source_commit'],
                  accepted['source_commit'], '--', 'Source', 'CMakeLists.txt').splitlines()
    check(changes == ['Source/Standalone/LraPassagesPanel.cpp',
                      'Source/Standalone/LraPassagesPanel.h'], 'Reuse restricted to UI legend correction')
    remote_integrity = read(root / 'high-sierra-final-integrity.json')
    check(remote_integrity['passed'] and remote_integrity['machine'] == 'x86_64'
          and '10.13.6' in remote_integrity['os'], 'Actual High Sierra package integrity')
    for item in remote_integrity['packages']:
        check(item['passed'] and not item['mismatches']
              and item['identity_sha256'] == sha(root / item['candidate'] / 'identity.json'),
              'Identical local/Golem frozen identity ' + item['candidate'])

    required = {'core', 'stereo', 'surround', 'rates', 'failures', 'inputs', 'gating',
                'real', 'measurement', 'section', 'requirements'}
    for host, architecture, version in [('sequoia', 'arm64', '15.7.9'),
                                         ('high-sierra', 'x86_64', '10.13.6')]:
        run = read(root / (host + '-01/run.json'))
        check(run['machine'] == architecture and version in run['os'], host + ': actual platform')
        check(run['required_automated_complete'] and len(run['cases']) == len(required)
              and {v['group'] for v in run['cases']} == required, host + ': all eleven groups')
        for row in run['cases']:
            report = read(root / (host + '-01') / (row['group'] + '.json'))
            check(row['exit_code'] == 0 and report['passed']
                  and report['build'] == identities['candidate-01']['measurement_build']
                  and assertions_passed(report),
                  host + ': passing ' + row['group'] + ' original receipt')
        native = read(root / (host + '-native-02/result.json'))
        check(native['passed'] and native['build'] == accepted['measurement_build']
              and all(v['passed'] for v in native['checks'])
              and native['passages_p95_ms'] <= 100 and native['cancel_ms'] < 1000
              and not native['physical_devices_opened'] and not native['video_loaded'],
              host + ': fresh candidate-02 native UI, cancellation and measured workflow')
        rates = read(root / (host + '-surround-rates.json'))
        check(rates['passed'] and rates['build'] == accepted['measurement_build']
              and all(v['passed'] for v in rates['checks'])
              and {v['rate'] for v in rates['cases']} == {44100, 96000, 192000},
              host + ': supplemental 5.1 rate cells')
        resources = read(root / ('resources-sequoia-02.json' if host == 'sequoia'
                                  else 'resources-high-sierra.json'))
        check(resources['passed'] and resources['cancelled_cleanly']
              and resources['grouping_seconds'] + resources['planning_seconds'] <= resources['budget_seconds']
              and resources['resident_bytes'] <= 64 * 1024 * 1024
              and not resources['full_duration_audio_decoded'], host + ': bounded population resources')
    check(read(root / 'numerical-ui-review.json')['numerical_and_native_complete'],
          'Original numerical/native cross-platform review complete')
    comparison = read(root / 'surround-rates-comparison.json')
    check(comparison['passed'] and len(comparison['comparisons']) == 3
          and all(all(v['checks'].values()) and all(v['numeric_fields_exactly_equal'].values())
                  for v in comparison['comparisons']), 'Supplemental cross-platform exact agreement')
    check(remote_integrity['supplemental_surround_rates_binary_sha256']
          == sha(root / 'supplemental/surround-rates'), 'Same supplemental executable on both Macs')
    for probe in read(root / 'supplemental-provenance.json')['probes']:
        check(sha(repo / probe['source']) == probe['source_sha256']
              and sha(repo / probe['executable']) == probe['executable_sha256'],
              'Supplemental source/binary identity ' + probe['source'])
    for name, count, monitor in [('listening-stereo', 4, -5), ('listening-stereo-repeat', 4, -5),
                                 ('listening-cinema', 4, -5), ('listening-cinema-reference-zero', 5, 0)]:
        receipt = read(root / (name + '.json'))
        check(receipt['passed'] and receipt['device_closed'] and receipt['monitor_db'] == monitor
              and len(receipt['clips']) == count and all(v['underruns'] == 0 for v in receipt['clips']),
              name + ': approved sequence complete and device closed')
    observation = read(root / 'listening-observations.json')['final_user_inspection']
    check(observation['accepted_for_bounded_m19_boundary_review']
          and observation['verbatim'] == 'I checked all 4 stereo and 5.1 files and they are flawless.',
          'Explicit user acceptance of bounded supplied-file inspection')
    before = [line.split() for line in (root / 'checkpoint-tags-before.txt').read_text().splitlines()]
    refs = dict(line.split()[::-1] for line in (root / 'remote-final-refs.txt').read_text().splitlines())
    for name, tag, target in before:
        check(git('rev-parse', name) == tag and git('rev-parse', name + '^{}') == target
              and refs['refs/tags/' + name] == tag and refs['refs/tags/' + name + '^{}'] == target,
              'Preserved local/remote tag object and target ' + name)
    main = '51e47f300fcc38a386a5a92a7277029175cfab0c'
    check(git('rev-parse', 'main') == main and refs['refs/heads/main'] == main, 'Main unchanged at M18')
    result.update(qualified=True, ready_for_checkpoint_review=True,
                  reviewer_sha256=sha(Path(__file__)),
                  branch='codex/m19-lra-passages', tested_source_commit=accepted['source_commit'],
                  source_manifest_sha256=accepted['source_manifest_sha256'],
                  measurement_build=accepted['measurement_build'],
                  main=main, existing_tags_preserved=len(before),
                  scope='Approved impact-selected stereo/5.1 M19; no new full-duration/video or physical multichannel routing qualification',
                  listening='Bounded coordinated playback plus user-reported supplied stereo/5.1 file inspection; earlier observations retained',
                  limitations=['Cinema 15 LU attempt fails honestly and cannot Apply',
                               'At most three private full-programme passes; no arbitrary LRA attainment guarantee',
                               'Mono/7.1 new automatic LRA method outside scope',
                               'Full physical multichannel routing/listening remains separate'])
except Exception as error:
    result['failure'] = str(error)
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'checks'}, indent=2))
raise SystemExit(0 if result['qualified'] else 1)
