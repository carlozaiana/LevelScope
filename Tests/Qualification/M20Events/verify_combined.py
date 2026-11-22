"""Finite acceptance and preservation receipt for the 22-November UI qualification."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

repo = Path(__file__).resolve().parents[3]
evidence = repo/'Tests/Qualification/M20Events/results/combined-2026-11-22'
work = repo/'build/m20-combined-qualification-2026-11-22'
candidate = repo/'build/m20-competing-boundary-ui-2026-11-21/candidate-03'
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as stream:
        for b in iter(lambda: stream.read(1024*1024), b''): h.update(b)
    return h.hexdigest()
def read(p): return json.loads(p.read_text())
def git(*args): return subprocess.check_output(['git',*args],cwd=repo,text=True).strip()

before=read(evidence/'preservation-before.json')
result={'passed':False,'checks':{},'qualification_scope':'Focused actual High Sierra, native Sequoia, cross-platform and independent programmes; video/endurance deferred'}
checks=result['checks']
def require(name, condition):
    checks[name]=bool(condition)
    assert condition,name
try:
    identity=read(candidate/'identity.json')
    require('final_identity',sha(candidate/'identity.json')==(candidate/'identity.sha256').read_text().split()[0])
    result['candidate']={'path':str(candidate),'identity_sha256':sha(candidate/'identity.json'),'source_commit':identity['source_commit'],'measurement_build':identity['measurement_build']}
    original=read(Path(before['candidate_path'])/'identity.json')
    original_sources={k:v for k,v in original['source_sha256'].items() if k.startswith('Source/')}
    final_sources={k:v for k,v in identity['source_sha256'].items() if k.startswith('Source/')}
    changed=[k for k in original_sources if original_sources[k]!=final_sources.get(k)]
    expected=['Source/Standalone/RegionProposalReviewPanel.cpp','Source/Standalone/RegionProposalReviewPanel.h','Source/Standalone/StandaloneProposalPreferences.cpp','Source/Standalone/StandaloneProposalPreferences.h']
    require('only_four_ui_source_files_changed',sorted(changed)==expected and original_sources.keys()==final_sources.keys())
    result['source_changes']=changed
    result['unchanged_production_sources']={k:v for k,v in final_sources.items() if k not in changed}
    packages=list(before['packages'])+[
        {'path':str(candidate.parent/'candidate-02'),'identity_sha256':'486adef440a6237a466cf0384b87dc7b0ccc02b0ef9261101640a6514c6c6b1c'},
        {'path':str(candidate),'identity_sha256':sha(candidate/'identity.json')}]
    for package in packages:
        p=Path(package['path']);require('package_identity:'+str(p),sha(p/'identity.json')==package['identity_sha256'])
        manifest=read(p/'identity.json')
        for k,v in manifest['artifacts_sha256'].items(): require('artifact:'+str(p/k),sha(p/k)==v)
    result['preserved_package_count']=len(packages)
    for group in ['manuals','listening_files']:
        for p,v in before[group].items():require(group+':'+p,sha(Path(p))==v)
    for p,v in before['untracked'].items():require('unrelated:'+p,sha(repo/p)==v)
    require('main_unchanged',git('rev-parse','main')==before['main'])
    require('tracking_main_unchanged',git('rev-parse','origin/main')==before['origin_main'])
    tags=git('for-each-ref','--format=%(refname) %(objectname) %(*objectname)','refs/tags').splitlines()
    require('all_tag_objects_targets_unchanged',tags==before['tags'])
    def worktree_records(text):
        return {block.splitlines()[0]:block for block in text.strip().split('\n\n')}
    old_trees=worktree_records(before['worktrees']);new_trees=worktree_records(git('worktree','list','--porcelain'))
    require('other_worktrees_unchanged',all(new_trees.get(k)==v for k,v in old_trees.items() if k!='worktree '+str(repo)) and old_trees.keys()==new_trees.keys())
    result['platform_groups']={}
    for platform,folder in [('sequoia','final-sequoia-01'),('high-sierra','golem-return/final-high-sierra-01')]:
        root=work/folder;summary=read(root/'summary.json')
        require(platform+':summary',summary['passed'] and summary['candidate_identity_sha256']==result['candidate']['identity_sha256'])
        require(platform+':actual_os_architecture',(summary['os'],summary['architecture'])==(('10.13.6','x86_64') if platform=='high-sierra' else ('15.7.9','arm64')))
        require(platform+':exact_groups',{g['name'] for g in summary['groups']}=={'inner-ui','event-ui','ui'})
        for g in summary['groups']:
            receipt=root/g['receipt'];d=read(receipt)
            require(platform+':'+g['name'],sha(receipt)==g['receipt_sha256'] and d['passed'] and all(c['passed'] for c in d['checks']) and d['build']==identity['measurement_build'])
        result['platform_groups'][platform]=summary['groups']
    numeric=read(work/'final-golem-numerics.json')
    require('replacement_cross_platform_exact',numeric['comparison']['passed'])
    require('replacement_numeric_engine_unchanged',read(work/'replacement-v-original-numerics.json')['passed'])
    result['cross_platform_numeric_receipt_sha256']=sha(work/'final-golem-numerics.json')
    for name in ['inner-stereo','inner-surround','inner-population']:
        d=read(evidence/('candidate-01/high-sierra-03/'+name+'.json'))
        require('inherited_golem:'+name,d['passed'] and all(c['passed'] for c in d['checks']))
    for name in ['DialogAmbiance-02','5.1_Busy-01']:
        d=read(evidence/('independent/'+name+'/result.json'))
        require('independent_rejection:'+name,d['passed'] and d['measured_passing_attempts']==0 and len(d['attempts'])==3 and not any(a['pass'] for a in d['attempts']))
    for name in ['independent/exchange-local-01','independent/exchange-from-golem-01','candidate-01/exchange-01']:
        d=read(evidence/(name+'/result.json'))
        require('exchange:'+name,d['passed'])
    exchange=read(evidence/'independent/exchange-comparison.json')
    require('exchange_numerics_and_cache_authority',exchange['opposite_platforms'] and not exchange['foreign_caches_used'] and exchange['all_numerical_fields_and_peak_positions_identical_in_json'])
    result['inherited_evidence']={
        'original_integrated_sequoia':'6 groups / 3004 assertions, inner-ui-01; affected UI replaced by new runs',
        'unchanged_logic':'28 groups / 3556 assertions, inner-01, by existing source reuse map',
        'actual_high_sierra_service_resource':'candidate-01 high-sierra-03: inner-stereo60, inner-surround57, inner-population17',
        'independent_holdouts':'fixed inputs/settings; 27 stereo + 24 surround checks, no passing correction',
        'exchange':'34 same-platform control + 35 checks each direction; schema11/programme9/history4; no metric deltas',
        'note':'Four UI presentation files differ. Underlying authority, engine, codec, export and test/harness sources match; Source-wide measurement build identity deliberately differs.'}
    require('user_development_listening_acceptance',read(evidence/'user-observations.json')['development_listening']['accepted_for_current_development_stage'])
    result['passed']=True
finally:
    (evidence/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
print('Combined focused qualification PASS; all preservation checks passed.')
