"""Run only the named offline boundary checks. Each run owns a fresh output directory."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
MUSIC = Path('/Users/fuzz/N15/Mixdown/MusicInBackAndForegroundSomeDialog.wav')
OLD = json.loads((HERE / 'results/boundary-investigation-01/provenance.json').read_text())['saved_sessions']
EVENT = 'LevelScopeEventPlanningTests'


def sha(path):
    with path.open('rb') as stream:
        result = hashlib.sha256()
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('binary_root', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--only', nargs='+', default=[])
    parser.add_argument('--compile', action='store_true')
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    bins = args.binary_root.resolve()
    identity = {} if args.compile else json.loads((bins/'identity.json').read_text())
    selected = {}

    def executable(target):
        return bins / (f'{target}_artefacts/Release' if args.compile else 'bin') / target

    def group(label, name, *extra, target=EVENT):
        selected[label] = [executable(target), '--group', name, '--root', root/label,
                           '--output', root/(label+'.json'), *extra]

    for mode, label in [(2, 'prior-margin20'), (6, 'inner-margin20'), (1, 'prior-strict20'),
                        (5, 'inner-strict20'), (4, 'prior-margin60'), (7, 'inner-margin60')]:
        group(label, 'external', '--input', MUSIC, '--experiment', str(mode), *(['--listening'] if mode in (2,6) else []))
    group('inner-evidence', 'inner-evidence', '--input', MUSIC)
    for name in ['inner-controls','inner-features','inner-population','inner-stereo','inner-surround',
                 'inner-fallback','core','features','stereo','fallback']:
        group(name, name)
    for name in ['m19','m20']:
        group('saved-'+name, 'saved-listening', '--input', Path(OLD[name]['path']))
    for label, file in [('gate', 'build/m20-logic-2026-11-06/fixtures-01/legacy-gating/gating.wav'),
                        ('cinema', 'build/m19-qualification-2026-11-05/media/cinema-27m-35m.wav')]:
        group(label, 'external', '--input', REPO/file, '--experiment', '7')
    for name in ['core','failures','stereo']:
        group('region-'+name,name,target='LevelScopeRegionPlanningTests')
    group('lra-core','core',target='LevelScopeLraProposalTests')
    for label, flags in [('ui',[]), ('event-ui',['--event-ui'])]:
        selected[label] = [executable('LevelScopeRegionUITests'),root/label,*flags]
    selected['timeline'] = [executable('LevelScopeTimelineProgressTests'),'--root',root/'timeline',
                            '--output',root/'timeline.json']
    wanted = args.only or list(selected)
    for label in wanted:
        if label not in selected:
            raise ValueError(label)
    if not args.compile:
        for target in {selected[label][0] for label in wanted}:
            assert sha(target)==identity['artifacts_sha256'][str(target.relative_to(bins))]
    output = dict(source=identity.get('source_base_commit'),measurement_build=identity.get('measurement_build'),
                  package_identity_sha256=(sha(bins/'identity.json') if identity else None),
                  system=platform.mac_ver()[0],architecture=platform.machine(),devices_opened=False,
                  planned=wanted,runs={})
    (root/'commands.json').write_text(json.dumps({k:list(map(str,selected[k])) for k in wanted},indent=2)+'\n')
    for label in wanted:
        print('START',label,flush=True)
        if label.startswith('saved-'):
            original=OLD[label.removeprefix('saved-')]
            assert sha(Path(original['path']))==original['sha256']
        start=time.monotonic()
        with (root/(label+'.log')).open('w') as log:
            proc=subprocess.run(list(map(str,selected[label])),stdout=log,stderr=subprocess.STDOUT,timeout=600)
        if label in ('ui','event-ui') and (root/label/'result.json').exists():
            shutil.copyfile(root/label/'result.json',root/(label+'.json'))
        receipt=root/(label+'.json')
        data=json.loads(receipt.read_text()) if receipt.exists() else {}
        passed=(proc.returncode==0 and data.get('passed') and data.get('checks')
                and all(check.get('passed') for check in data['checks'])
                and (args.compile or data.get('build')==identity['measurement_build']))
        output['runs'][label]=dict(passed=bool(passed),code=proc.returncode,wall_s=round(time.monotonic()-start,3),
                                   checks=len(data.get('checks',[])),failure=data.get('failure',''),
                                   receipt_sha256=sha(receipt) if receipt.exists() else None)
        (root/'run.json').write_text(json.dumps(output,indent=2)+'\n')
        print('PASS' if passed else 'FAIL',label,data.get('failure',''),flush=True)
    return int(not all(r['passed'] for r in output['runs'].values()))


if __name__=='__main__':
    raise SystemExit(main())
