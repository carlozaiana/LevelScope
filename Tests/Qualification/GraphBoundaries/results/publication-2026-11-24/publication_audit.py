import hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path('/Users/fuzz/.codex/worktrees/a411/LevelScope')
WORK=ROOT/'build/shared-source-graph-publication-2026-11-24'
INTEGRATION=ROOT/'build/m18-publication-2026-11-03/integration'
BRANCH='codex/graph-region-boundaries-2026-11-24'
TAG='m20-shared-source-graph-boundaries-2026-11-24'
CHECKPOINT='2590d3fbdaaf437e313f8646021039d3467c4ef2'
TESTED='f801ac06307f25f423a953a7e35b59c28013ccc7'
SOURCE='31be38bd890fee3667a6a976fbc3524a2b26e609'
MAIN='18adeb520006b7fe7e0a22300cc26abd28eec51c'
QUAL=ROOT/'Tests/Qualification/GraphBoundaries/results/combined-2026-11-24'
def git(*args,cwd=ROOT):return subprocess.check_output(['git',*args],cwd=cwd,text=True).strip()
def read(p):return json.loads(p.read_text())
def save(p,v):
 with p.open('x') as f:json.dump(v,f,indent=2);f.write('\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def refs(p):return {line.split('\t')[1]:line.split('\t')[0] for line in p.read_text().splitlines() if '\t' in line and not line.startswith('ref:')}
def tags():return {line.split(' ')[0]:line.split(' ')[1:] for line in git('for-each-ref','--format=%(refname) %(objectname) %(*objectname)','refs/tags').splitlines()}
def ancestor(a,b):subprocess.run(['git','merge-base','--is-ancestor',a,b],cwd=ROOT,check=True)
def preflight():
 assert git('remote','get-url','origin')=='https://github.com/carlozaiana/LevelScope.git'
 assert git('branch','--show-current')==BRANCH and git('rev-parse','HEAD')==CHECKPOINT
 assert not git('diff','HEAD','--') and not git('diff','HEAD','--',cwd=INTEGRATION)
 assert git('rev-parse','HEAD',cwd=INTEGRATION)==MAIN
 remote=refs(WORK/'remote-before.txt');local=tags()
 assert remote['refs/heads/main']==git('rev-parse','main')==git('rev-parse','origin/main')==MAIN
 assert {k for k in remote if k.startswith('refs/tags/') and not k.endswith('^{}')}==set(local)
 for ref,(obj,peeled) in local.items():assert remote[ref]==obj and remote.get(ref+'^{}')==peeled
 assert 'refs/tags/'+TAG not in remote and 'refs/heads/'+BRANCH not in remote
 for a in [MAIN,'f9eeae554421b4ef4f8b48d53a69785aa3f6af46','614377f4e573bdd0fc34726584418b9333d092c2',TESTED]:ancestor(a,CHECKPOINT)
 assert git('rev-parse','HEAD:Source')==git('rev-parse',TESTED+':Source')==SOURCE
 assert not git('diff',TESTED,CHECKPOINT,'--','Source','CMakeLists.txt')
 qual=read(QUAL/'acceptance.json');assert qual['focused_both_platform_qualification_passed'] and len(qual['actual_platforms'])==2
 assert all(p['passed'] and p['assertions']==1506 for p in qual['actual_platforms'])
 prior=read(ROOT/'build/boundaries-workflow-publication-2026-11-23/preflight.json');protected=dict(prior['protected'])
 for name,digest in protected.items():assert sha(Path(name))==digest,name
 definitions=prior['packages']+[{'path':str(Path('/Users/fuzz/Developer/LevelScope/build/session-source-2026-11-23')/name)} for name in ['candidate-01','candidate-02','candidate-03']]+[{'path':str(ROOT/'build/graph-region-boundaries-2026-11-24/candidate-01'),'identity_sha256':'cecac62f3f98a5c92e1b87d938ab94f813b32f7cd2984cf8e7156a7bf7af43cf'}]
 packages=[]
 for definition in definitions:
  p=Path(definition['path']);digest=sha(p/'identity.json')
  if 'identity_sha256' in definition:assert digest==definition['identity_sha256']
  if p==Path('/Users/fuzz/Developer/LevelScope/build/session-source-2026-11-23/candidate-03'):assert digest=='94c42f1eca24de3adc671723b4b01cd4bcc986fe602b4c6c6d25484379079843'
  manifest=read(p/'identity.json')
  for name,h in manifest['artifacts_sha256'].items():assert sha(p/name)==h,str(p/name);protected[str(p/name)]=h
  for name in ['identity.json','identity.sha256']:
   if (p/name).exists():protected[str(p/name)]=sha(p/name)
  packages.append({'path':str(p),'identity_sha256':digest,'artifacts':len(manifest['artifacts_sha256'])})
 candidate=Path(definitions[-1]['path']);manifest=read(candidate/'identity.json');entries=0
 for name,digest in manifest['artifacts_sha256'].items():
  if name.startswith('source/'):
   assert hashlib.sha256(subprocess.check_output(['git','show',CHECKPOINT+':'+name[7:]],cwd=ROOT)).hexdigest()==digest,name
   entries+=1
 for folder in ['Tests/Qualification/SharedSource/results','Tests/Qualification/GraphBoundaries/results','Tests/Qualification/PlaybackPreparation/results']:
  for p in (ROOT/folder).rglob('*'):
   if p.is_file():protected[str(p)]=sha(p)
 worktrees=[]
 for block in git('worktree','list','--porcelain').split('\n\n'):
  p=Path(block.splitlines()[0][9:]);untracked={}
  for name in subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=p).decode().split('\0'):
   if name and (p/name).is_file():untracked[name]=sha(p/name)
  worktrees.append({'path':str(p),'head':git('rev-parse','HEAD',cwd=p),'branch':git('branch','--show-current',cwd=p),'tracked_diff_sha256':hashlib.sha256(subprocess.check_output(['git','diff','HEAD','--binary'],cwd=p)).hexdigest(),'untracked':untracked})
 save(WORK/'preflight.json',{'passed':True,'repository':git('remote','get-url','origin'),'branch':BRANCH,'checkpoint':CHECKPOINT,'implementation':TESTED,'source_tree':SOURCE,'main_before':MAIN,'tag':TAG,'tags':local,'remote':remote,'packages':packages,'protected':protected,'worktrees':worktrees,'manifest_entries_verified':entries,'qualification_receipt_sha256':sha(QUAL/'acceptance.json'),'user_acceptance':read(WORK/'user-acceptance.json')})
 print('Preflight PASS:',len(packages),'packages,',len(protected),'protected files,',len(local),'tags,',entries,'source entries; fast-forward available.')
def verify(label,expected):
 before=read(WORK/'preflight.json');remote=refs(WORK/('remote-'+label+'.txt'));local=tags();checks={}
 def check(name,ok):checks[name]=bool(ok);assert ok,name
 for name,digest in before['protected'].items():check('protected:'+name,sha(Path(name))==digest)
 for ref,pair in before['tags'].items():check('tag:'+ref,local.get(ref)==pair and remote.get(ref)==pair[0] and remote.get(ref+'^{}')==pair[1])
 new='refs/tags/'+TAG;obj=git('rev-parse',new)
 check('annotated',git('cat-file','-t',obj)=='tag')
 check('checkpoint_tag',local.get(new)==[obj,CHECKPOINT] and remote.get(new)==obj and remote.get(new+'^{}')==CHECKPOINT)
 check('tag_inventory',set(local)==set(before['tags'])|{new})
 for ref,value in before['remote'].items():
  if ref not in ['HEAD','refs/heads/main']:check('unrelated_remote:'+ref,remote.get(ref)==value)
 check('remote_inventory',set(remote)==set(before['remote'])|{new,new+'^{}','refs/heads/'+BRANCH})
 for ref in ['HEAD','refs/heads/main','refs/heads/'+BRANCH]:check('remote:'+ref,remote.get(ref)==expected)
 check('local_heads',git('rev-parse','main')==git('rev-parse',BRANCH)==git('rev-parse','origin/main')==expected)
 check('exact_source',git('rev-parse',expected+':Source')==SOURCE and not git('diff',TESTED,expected,'--','Source','CMakeLists.txt'))
 ancestor(CHECKPOINT,expected)
 for name in ['Docs/ProjectIndex/94_SharedVerifiedSourceAndPlaybackPreparation.md','Docs/ProjectIndex/95_GraphAppliedRegionBoundaries.md','Tests/Qualification/GraphBoundaries/results/combined-2026-11-24/acceptance.json']:
  check('handover:'+name,subprocess.check_output(['git','show',CHECKPOINT+':'+name],cwd=ROOT)==subprocess.check_output(['git','show',expected+':'+name],cwd=ROOT))
 finder=[]
 for state in before['worktrees']:
  p=Path(state['path']);check('head:'+str(p),git('rev-parse','HEAD',cwd=p)==(expected if p in [ROOT,INTEGRATION] else state['head']))
  check('branch:'+str(p),git('branch','--show-current',cwd=p)==state['branch'])
  check('edits:'+str(p),hashlib.sha256(subprocess.check_output(['git','diff','HEAD','--binary'],cwd=p)).hexdigest()==state['tracked_diff_sha256'])
  for name,digest in state['untracked'].items():
   actual=sha(p/name) if (p/name).is_file() else None
   if actual!=digest and Path(name).name=='.DS_Store':finder.append(str(p/name))
   else:check('untracked:'+str(p/name),actual==digest)
 save(WORK/(label+'-verification.json'),{'passed':all(checks.values()),'stage':label,'main':expected,'checkpoint':CHECKPOINT,'tag':TAG,'tag_object':obj,'branch':BRANCH,'source_tree':SOURCE,'packages_preserved':len(before['packages']),'protected_files':len(before['protected']),'old_tags_preserved':len(before['tags']),'manifest_entries_verified':before['manifest_entries_verified'],'finder_metadata_changes':finder,'checks':checks,'runtime_tests_repeated':False,'rebuild_performed':False,'physical_devices_opened':False,'audible_tests_run':False})
 print(label,'PASS:',len(checks),'checks; main',expected,'; annotated tag',obj)
if __name__=='__main__':
 if sys.argv[1]=='preflight':preflight()
 else:verify(sys.argv[1],sys.argv[2])
