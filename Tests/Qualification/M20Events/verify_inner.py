"""Verify the frozen run and assemble local, source-aligned listening copies."""
import argparse
import json
import random
import shutil
from pathlib import Path

from run_inner import MUSIC, OLD, sha


def main():
    parser=argparse.ArgumentParser()
    for name in ('package','run','destination'):
        parser.add_argument(name,type=Path)
    args=parser.parse_args()
    package,run,dest=(p.resolve() for p in (args.package,args.run,args.destination))
    identity=json.loads((package/'identity.json').read_text())
    execution=json.loads((run/'run.json').read_text())
    assert execution['source']==identity['source_base_commit']
    assert execution['package_identity_sha256']==sha(package/'identity.json')
    for section,folder in [('artifacts_sha256',package),('source_files_sha256',package/'source')]:
        for path,digest in identity[section].items():
            assert sha(folder/path)==digest,path
    expected=set('prior-margin20 inner-margin20 prior-strict20 inner-strict20 prior-margin60 inner-margin60 inner-evidence inner-controls inner-features inner-population inner-stereo inner-surround inner-fallback core features stereo fallback saved-m19 saved-m20 gate cinema region-core region-failures region-stereo lra-core ui event-ui timeline'.split())
    assert set(execution['runs'])==expected and len(execution['planned'])==len(expected)
    cases={}
    for label,receipt in execution['runs'].items():
        path=run/(label+'.json')
        data=json.loads(path.read_text())
        assert data['passed'] and receipt['passed'] and receipt['code']==0 and receipt['receipt_sha256']==sha(path),label
        assert data['build']==identity['measurement_build'] and all(x['passed'] for x in data['checks']),label
        cases[label]=data.get('cases',[])
    assert sha(MUSIC)=='7ba6404041c94c909df290a6829a0fe92ce594329a1d575490c2bb4a5b6cc8f0'
    for label,rows in cases.items():
        if label.startswith(('prior-','inner-margin','inner-strict')) or label in ('gate','cinema'):
            row=rows[0]
            assert row['passes']==len(row['attempts'])<=3 and row['proxies']<=64 and row['planning_bytes']<=64*1024**2
            assert row['actual_assessment_peak_heap_bytes']<=64*1024**2
            assert row['feature_high_water']<=10*1024**2
            if label not in ('gate','cinema'):
                assert row['content_sha256']==sha(MUSIC) and row['frames']==13228800
    for control in ['margin20','strict20','margin60']:
        prior=cases['prior-'+control][0]
        new=cases['inner-'+control][0]
        assert prior['settings']['schema']==2 and new['settings']['schema']==3
        assert prior['current_snapshot']==new['current_snapshot'] and prior['source_snapshot']==new['source_snapshot']
        for key in ('mode','leftMs','rightMs','global','automaticMargins','input','splits','protected','cuts','minimumMs','separationMs'):
            assert prior['settings'][key]==new['settings'][key],key
    primary=cases['inner-margin20'][0]
    assert any(a['family'].startswith('joint') and a['pass'] and a.get('listening_session') for a in primary['attempts'])
    assert all(a['coverage']!='Supported complete' for a in primary['attempts'])
    assert 'Review required' in primary['outcome']
    for original in OLD.values():
        assert sha(Path(original['path']))==original['sha256']
    dest.mkdir(exist_ok=False)
    sessions=dest/'sessions';sessions.mkdir()
    blind=dest/'blind';blind.mkdir()
    shutil.copyfile(MUSIC,dest/'Source.wav')
    copied={}
    for name,original in OLD.items():
        file=sessions/(name+'.levelscope');shutil.copyfile(original['path'],file)
        assert sha(file)==original['sha256']
        copied[str(file.relative_to(dest))]=sha(file)
    neutral=sessions/'Source-reference.levelscope'
    shutil.copyfile(run/'inner-margin20/source-reference.levelscope',neutral)
    copied[str(neutral.relative_to(dest))]=sha(neutral)
    candidates={}
    def add(path,entry):
        digest=sha(path);candidates.setdefault(digest,dict(path=path,origins=[]))['origins'].append(entry)
    for name in OLD:
        add(run/('saved-'+name)/'programme.wav',dict(origin='saved-'+name,classification='preserved listening session',session=f'sessions/{name}.levelscope'))
    for label in ('prior-margin20','inner-margin20'):
        for i,a in enumerate(cases[label][0]['attempts'],1):
            entry=dict(origin=label,attempt=i,receipt=a['receipt'],plan=a['plan'],coverage=a['coverage'],
                       selected_rules_pass=a['pass'],il=a['il'],lra=a['lra'],tp=a['tp'],global_db=a['global_db'])
            if 'listening_session' in a:
                file=sessions/(label+'-'+a['listening_session'])
                shutil.copyfile(run/label/a['listening_session'],file)
                copied[str(file.relative_to(dest))]=sha(file)
                entry['session']=str(file.relative_to(dest))
            add(run/label/a['listening_pcm'],entry)
    order=list(candidates);random.Random(20261120).shuffle(order)
    key={}
    for n,digest in enumerate(order,1):
        label=f'Comparison-{n:02}.wav'
        shutil.copyfile(candidates[digest]['path'],blind/label)
        assert sha(blind/label)==digest
        key[label]=dict(sha256=digest,origins=candidates[digest]['origins'])
    manifest=dict(implementation_commit=identity['source_base_commit'],measurement_build=identity['measurement_build'],
                  source_sha256=sha(dest/'Source.wav'),frames=13228800,sample_rate=48000,
                  intended_gains=True,listening_performed=False,comparisons=key,sessions=copied,
                  controls={name:cases[name][0] for name in ('prior-strict20','inner-strict20','prior-margin60','inner-margin60','gate','cinema')})
    (dest/'comparison-key.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (dest/'Start here.md').write_text('''# Listening comparison

Use Source.wav as the common reference. Every WAV in the blind folder starts at the same source frame and runs the full 4 min 35.6 s. The files use the intended delivery gains. Both channels belong to the programme. Disable player normalization, crossfades and other gain changes.

Start with a continuous pass. Then compare 9–33 s and 48–64 s, keeping enough lead-in and tail to hear the bed. Mark the first relevant sound, the last comparable hit and the point where you want gain to return. Check the trough between the two louder parts of the later phrase. Record any lifted tail, late entrance, pumping or unnecessary attenuation of speech or music elsewhere. Keep a per-passage preference and a whole-programme preference. A separately logged monitor-level-matched pass can help identify a global-gain preference.

The WAVs include measured failures as listening controls. Identical renders appear once. Open comparison-key.json after recording independent impressions to see the method, pass/fail, coverage and exact plan. The sessions folder contains preserved originals as copies and independently remeasured passing alternatives. The revised sessions require the supplied local reader. Use the WAVs for a self-contained comparison.
''')
    print(f"Verified {len(expected)} groups, {sum(r['checks'] for r in execution['runs'].values())} assertions, {len(key)} distinct blinded renders.")


if __name__=='__main__':
    main()
