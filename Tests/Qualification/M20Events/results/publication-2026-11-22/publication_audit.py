import hashlib
import json
from pathlib import Path
import subprocess
import sys

work=Path(__file__).resolve().parent
root=work.parent.parent
qualification=root/'Tests/Qualification/M20Events/results/combined-2026-11-22'
branch='codex/m20-competing-boundary-ui-2026-11-21'
tag='m20-integrated-competing-boundaries-2026-11-22'
checkpoint='7dafaa8b0f7bde7dd42a6d95f6cec6dbbd3914da'
tested='2f90dab6eaa25f9d2ce9f033341c948b2994d9da'
integration=root/'build/m18-publication-2026-11-03/integration'
def read(p): return json.loads(p.read_text())
def git(*args,cwd=root): return subprocess.check_output(['git',*args],cwd=cwd,text=True).strip()
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as stream:
        for data in iter(lambda:stream.read(1024*1024),b''):h.update(data)
    return h.hexdigest()
def refs(path):
    return {line.split('\t')[1]:line.split('\t')[0] for line in path.read_text().splitlines() if '\t' in line and not line.startswith('ref:')}
def save(path,data):
    assert not path.exists(),path
    path.write_text(json.dumps(data,indent=2)+'\n')
def preflight():
    previous=read(qualification/'preservation-before.json')
    assert read(qualification/'acceptance.json')['passed']
    protected={};packages=[]
    definitions=previous['packages']+[
        {'path':str(root/'build/m20-competing-boundary-ui-2026-11-21/candidate-02'),'identity_sha256':'486adef440a6237a466cf0384b87dc7b0ccc02b0ef9261101640a6514c6c6b1c'},
        {'path':str(root/'build/m20-competing-boundary-ui-2026-11-21/candidate-03'),'identity_sha256':'f44c1fea69a41416c3c1dd0284a3a5b9cd52757bc484f9623cbb7fb13aa48dd3'}]
    for package in definitions:
        p=Path(package['path']);manifest=read(p/'identity.json')
        assert sha(p/'identity.json')==package['identity_sha256']
        for name,digest in manifest['artifacts_sha256'].items():
            assert sha(p/name)==digest,str(p/name);protected[str(p/name)]=digest
        for name in ['identity.json','identity.sha256']:
            if (p/name).exists():protected[str(p/name)]=sha(p/name)
        packages.append(dict(package,indexed_artifacts=len(manifest['artifacts_sha256'])))
    for group in ['manuals','listening_files']:
        for name,digest in previous[group].items():
            assert sha(Path(name))==digest,name;protected[name]=digest
    for folder in [qualification,root/'Tests/Qualification/M20Events/results/inner-ui-01',root/'Tests/Qualification/M20Events/results/inner-01']:
        for p in folder.rglob('*'):
            if p.is_file():protected[str(p)]=sha(p)
    media=[('DialogAmbiance.wav','602485eef95774826a575b0069d358e6723a31f6928dbd6e04dd0df44c98eb90'),('5.1_Busy.wav','4fa30bfd9a8fbb64e68fc00eaf11f8adf81ead8502d05ce639b22cc07645d7b6'),('MusicInBackAndForegroundSomeDialog.wav','7ba6404041c94c909df290a6829a0fe92ce594329a1d575490c2bb4a5b6cc8f0')]
    for name,digest in media:
        p=Path('/Users/fuzz/N15/Mixdown')/name;assert sha(p)==digest;protected[str(p)]=digest
    remote=refs(work/'remote-before.txt')
    tags=[line.split(' ') for line in (work/'local-tags-before.txt').read_text().splitlines()];assert len(tags)==26
    for ref,obj,peeled in tags:assert remote[ref]==obj and remote[ref+'^{}']==peeled
    assert remote['refs/heads/main']==git('rev-parse','main')==git('rev-parse','origin/main')=='22961f03de7919e1de631f8b13eb4b852108ac0b'
    assert git('branch','--show-current')==branch and git('rev-parse','HEAD')==checkpoint
    assert 'refs/heads/'+branch not in remote and 'refs/tags/'+tag not in remote
    assert not git('diff','HEAD','--','Source','Tests','CMakeLists.txt')
    subprocess.run(['git','merge-base','--is-ancestor',remote['refs/heads/main'],checkpoint],cwd=root,check=True)
    assert git('rev-parse',tested+':Source')==git('rev-parse',checkpoint+':Source')
    assert not git('diff',tested,checkpoint,'--','Source','CMakeLists.txt')
    manifest=read(root/'build/m20-competing-boundary-ui-2026-11-21/candidate-03/identity.json')
    for name,digest in manifest['source_sha256'].items():
        assert hashlib.sha256(subprocess.check_output(['git','show',checkpoint+':'+name],cwd=root)).hexdigest()==digest,name
    worktrees=[]
    for block in git('worktree','list','--porcelain').split('\n\n'):
        path=Path(block.splitlines()[0][9:]);tracked=git('diff','HEAD','--',cwd=path);untracked={}
        for name in subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=path).decode().split('\0'):
            if name and (path/name).is_file():untracked[name]=sha(path/name)
        worktrees.append({'path':str(path),'head':git('rev-parse','HEAD',cwd=path),'record':block,'tracked_diff_sha256':hashlib.sha256(tracked.encode()).hexdigest(),'status':git('status','--porcelain=v1',cwd=path),'untracked':untracked})
    assert not git('diff','HEAD','--',cwd=integration)
    result={'authorization':{'source_chat':'01a13c8f-90a9-7150-813f-cb057626ef7e','user_message':'01a203b7-041f-7793-8ca1-938e1fd5877e','text':'OK. Please tell the UI chat directly to publish the integrated M20 checkpoint to GitHub and merge it into main, preserving all existing tags, the tested version and the agreed qualification deferrals.'},'repository':git('remote','get-url','origin'),'branch':branch,'checkpoint':checkpoint,'tested':tested,'source_tree':git('rev-parse',checkpoint+':Source'),'tag':tag,'main_before':remote['refs/heads/main'],'remote_before':remote,'local_tags_before':tags,'packages':packages,'protected':protected,'worktrees':worktrees,'candidate_source_entries_verified':len(manifest['source_sha256']),'qualification_receipt_sha256':sha(qualification/'acceptance.json'),'passed':True}
    save(work/'preflight.json',result)
    print('Preflight PASS:',len(packages),'packages;',len(protected),'protected files;',len(tags),'existing annotated tags;',len(manifest['source_sha256']),'candidate source entries. Fast-forward available.')
def verify(label,expected):
    before=read(work/'preflight.json');remote=refs(work/('remote-'+label+'.txt'));checks={}
    def check(name,ok):checks[name]=bool(ok);assert ok,name
    for name,digest in before['protected'].items():check('protected:'+name,sha(Path(name))==digest)
    local_tags={row.split(' ')[0]:row.split(' ')[1:] for row in git('for-each-ref','--format=%(refname) %(objectname) %(*objectname)','refs/tags').splitlines()}
    for ref,obj,peeled in before['local_tags_before']:
        check('old_local_tag:'+ref,local_tags.get(ref)==[obj,peeled]);check('old_remote_tag:'+ref,remote.get(ref)==obj and remote.get(ref+'^{}')==peeled)
    new_ref='refs/tags/'+tag;tag_object=git('rev-parse',new_ref)
    check('new_tag_annotated',git('cat-file','-t',tag_object)=='tag')
    check('new_tag_target',git('rev-parse',new_ref+'^{}')==checkpoint and remote.get(new_ref)==tag_object and remote.get(new_ref+'^{}')==checkpoint)
    check('exact_tag_inventory',set(local_tags)=={x[0] for x in before['local_tags_before']}|{new_ref})
    for ref,value in before['remote_before'].items():
        if ref not in ['HEAD','refs/heads/main']:check('unrelated_remote_ref:'+ref,remote.get(ref)==value)
    check('exact_remote_inventory',set(remote)==set(before['remote_before'])|{'refs/heads/'+branch,new_ref,new_ref+'^{}'})
    for ref in ['refs/heads/main','refs/heads/'+branch]:check('published:'+ref,remote.get(ref)==expected)
    check('remote_default_head',remote.get('HEAD')==expected)
    check('local_heads',git('rev-parse','main')==git('rev-parse',branch)==git('rev-parse','origin/main')==expected)
    check('tested_source',git('rev-parse',expected+':Source')==before['source_tree'] and not git('diff',tested,expected,'--','Source','CMakeLists.txt'))
    subprocess.run(['git','merge-base','--is-ancestor',checkpoint,expected],cwd=root,check=True)
    for name in ['Docs/ProjectIndex/89_M20CombinedFocusedQualification.md','Tests/Qualification/M20Events/results/combined-2026-11-22/acceptance.json']:
        check('qualified_evidence:'+name,subprocess.check_output(['git','show',expected+':'+name],cwd=root)==subprocess.check_output(['git','show',checkpoint+':'+name],cwd=root))
    finder_changes=[]
    for state in before['worktrees']:
        path=Path(state['path']);mutable=path in [root,integration]
        check('worktree_head:'+str(path),git('rev-parse','HEAD',cwd=path)==(expected if mutable else state['head']))
        check('worktree_tracked_changes:'+str(path),hashlib.sha256(git('diff','HEAD','--',cwd=path).encode()).hexdigest()==state['tracked_diff_sha256'])
        for name,digest in state['untracked'].items():
            actual=sha(path/name) if (path/name).is_file() else None
            if actual!=digest and Path(name).name=='.DS_Store':finder_changes.append(str(path/name))
            else:check('unrelated_work:'+str(path/name),actual==digest)
    result={'passed':all(checks.values()),'stage':label,'main':expected,'checkpoint':checkpoint,'tag':tag,'tag_object':tag_object,'branch':branch,'source_tree':before['source_tree'],'protected_files':len(before['protected']),'old_tags_preserved':len(before['local_tags_before']),'candidate_source_entries_verified':before['candidate_source_entries_verified'],'finder_metadata_changes':finder_changes,'checks':checks,'runtime_tests_repeated':False,'rebuild_performed':False}
    save(work/(label+'-verification.json'),result);print(label,'PASS:',len(checks),'checks; main',expected,'tag object',tag_object)
if sys.argv[1]=='preflight':preflight()
else:verify(sys.argv[1],sys.argv[2])
