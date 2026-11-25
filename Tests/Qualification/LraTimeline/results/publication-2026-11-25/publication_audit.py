import hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path('/Users/fuzz/.codex/worktrees/a411/LevelScope')
WORK=ROOT/'build/lra-navigation-publication-2026-11-25'
INTEGRATION=ROOT/'build/m18-publication-2026-11-03/integration'
BRANCH='codex/lra-lanes-resizing-2026-11-24'
TAG='m20-navigation-lra-lanes-resizing-2026-11-25'
HANDOVER='670e987fc4bc84eaa4d5f222ea8b150d5cb9ed7b'
TESTED='61c93ae86d49e4d35bfec798f6568344f5cb8fc5'
SOURCE='97492cc70c3ef845dbfcbe65e946413380e7f4f6'
MAIN='b5ff52b9c87422b28f27a089c66b746fd2447a05'
CANDIDATE=ROOT/'build/lra-lanes-resizing-2026-11-24/candidate-02'
QUAL=ROOT/'Tests/Qualification/LraTimeline/results/2026-11-24'
def git(*args,cwd=ROOT):return subprocess.check_output(['git',*args],cwd=cwd,text=True).strip()
def read(p):return json.loads(p.read_text())
def save(p,v):
 with p.open('x') as f:json.dump(v,f,indent=2);f.write('\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def refs(p):return {line.split('\t')[1]:line.split('\t')[0] for line in p.read_text().splitlines() if '\t' in line}
def tags():return {line.split(' ')[0]:line.split(' ')[1:] for line in git('for-each-ref','--format=%(refname) %(objectname) %(*objectname)','refs/tags').splitlines()}
def ancestor(a,b):subprocess.run(['git','merge-base','--is-ancestor',a,b],cwd=ROOT,check=True)
def preflight():
 assert git('remote','get-url','origin')=='https://github.com/carlozaiana/LevelScope.git'
 assert git('branch','--show-current')==BRANCH and git('rev-parse','HEAD')==HANDOVER
 assert not git('diff','HEAD','--') and not git('diff','HEAD','--',cwd=INTEGRATION)
 assert git('branch','--show-current',cwd=INTEGRATION)=='main' and git('rev-parse','HEAD',cwd=INTEGRATION)==MAIN
 remote=refs(WORK/'remote-before.txt');local=tags()
 assert remote['refs/heads/main']==git('rev-parse','main')==git('rev-parse','origin/main')==MAIN
 assert {k for k in remote if k.startswith('refs/tags/') and not k.endswith('^{}')}==set(local)
 for ref,(obj,peeled) in local.items():assert remote[ref]==obj and remote.get(ref+'^{}')==peeled
 assert 'refs/tags/'+TAG not in remote and 'refs/heads/'+BRANCH not in remote
 for a in [MAIN,'583bb6114728483f28a11a62b5c77a669d0433c5','68a2cd1f50ed8d651957e484636421dcea301f4e','fbd650cf78e8dc8b36bfaee6d8626a9aec816e1f',TESTED]:ancestor(a,HANDOVER)
 assert git('rev-parse','HEAD:Source')==git('rev-parse',TESTED+':Source')==SOURCE
 assert not git('diff',TESTED,HANDOVER,'--','Source','CMakeLists.txt')
 q=read(QUAL/'combined.json');assert q['passed'] and q['per_host_assertions']==2762
 for platform in ['sequoia-02','high-sierra-02']:
  s=read(QUAL/platform/'summary.json');assert s['passed'] and sum(g['checks'] for g in s['groups'])==2762
  assert s['candidate_identity_sha256']=='843f513b6f0afed08b13eda0ca13692721b2506aa65df5dfaddd500a54b5cab6'
 prior=read(ROOT/'build/lra-lanes-resizing-2026-11-24/preservation-before.json');protected=dict(prior['protected'])
 for name,digest in protected.items():assert sha(Path(name))==digest,name
 # Extend the inherited inventory to every local preserved candidate directory,
 # including failed candidates; never rewrite their manifests or frozen files.
 packages=[]
 for identity in sorted((ROOT/'build').glob('*/candidate-*/identity.json')):
  p=identity.parent;d=read(identity)
  if 'artifacts_sha256' in d:
   for name,h in d['artifacts_sha256'].items():assert sha(p/name)==h,str(p/name)
  for f in p.rglob('*'):
   if f.is_file():protected[str(f)]=sha(f)
  packages.append({'path':str(p),'identity_sha256':sha(identity),'inventoried_files':sum(f.is_file() for f in p.rglob('*'))})
 for folder in ['Tests/Qualification/LraTimeline','Tests/Qualification/GainNavigation','Tests/Qualification/CurrentSections']:
  for p in (ROOT/folder).rglob('*'):
   if p.is_file():protected[str(p)]=sha(p)
 for name in ['Docs/ProjectIndex/100_GainEntryAndBoundaryNavigation.md','Docs/ProjectIndex/101_MeasuredLRALanesAndIndependentHeights.md']:
  p=ROOT/name;protected[str(p)]=sha(p)
 d=read(CANDIDATE/'identity.json');assert sha(CANDIDATE/'identity.json')=='843f513b6f0afed08b13eda0ca13692721b2506aa65df5dfaddd500a54b5cab6'
 entries=0
 for name,h in d['artifacts_sha256'].items():
  if name.startswith('source/'):
   assert hashlib.sha256(subprocess.check_output(['git','show',TESTED+':'+name[7:]],cwd=ROOT)).hexdigest()==h,name
   entries+=1
 worktrees=[]
 for block in git('worktree','list','--porcelain').split('\n\n'):
  p=Path(block.splitlines()[0][9:]);untracked={}
  for name in subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=p).decode().split('\0'):
   if name and (p/name).is_file():untracked[name]=sha(p/name)
  worktrees.append({'path':str(p),'head':git('rev-parse','HEAD',cwd=p),'branch':git('branch','--show-current',cwd=p),'tracked_diff_sha256':hashlib.sha256(subprocess.check_output(['git','diff','HEAD','--binary'],cwd=p)).hexdigest(),'untracked':untracked})
 save(WORK/'preflight.json',{'passed':True,'repository':git('remote','get-url','origin'),'branch':BRANCH,'qualification_handover':HANDOVER,'tested_checkpoint':TESTED,'source_tree':SOURCE,'main_before':MAIN,'tag':TAG,'tags':local,'remote':remote,'local_candidates':packages,'inherited_protected_files':len(prior['protected']),'protected':protected,'worktrees':worktrees,'manifest_entries_verified':entries,'qualification_receipt_sha256':sha(QUAL/'combined.json'),'user_acceptance':read(WORK/'user-acceptance.json')})
 print('Preflight PASS:',len(packages),'local candidates,',len(protected),'protected files,',len(local),'tags,',entries,'manifest source entries; fast-forward available.')
def verify(label,expected):
 before=read(WORK/'preflight.json');remote=refs(WORK/('remote-'+label+'.txt'));local=tags();checks={}
 checkpoint=(WORK/'checkpoint.txt').read_text().strip()
 def check(name,ok):checks[name]=bool(ok);assert ok,name
 for name,digest in before['protected'].items():check('protected:'+name,sha(Path(name))==digest)
 for ref,pair in before['tags'].items():check('tag:'+ref,local.get(ref)==pair and remote.get(ref)==pair[0] and remote.get(ref+'^{}')==pair[1])
 new='refs/tags/'+TAG;obj=git('rev-parse',new)
 check('annotated',git('cat-file','-t',obj)=='tag')
 check('checkpoint_tag',local.get(new)==[obj,checkpoint] and remote.get(new)==obj and remote.get(new+'^{}')==checkpoint)
 check('tag_inventory',set(local)==set(before['tags'])|{new})
 for ref,value in before['remote'].items():
  if ref not in ['HEAD','refs/heads/main']:check('unrelated_remote:'+ref,remote.get(ref)==value)
 check('remote_inventory',set(remote)==set(before['remote'])|{new,new+'^{}','refs/heads/'+BRANCH})
 for ref in ['HEAD','refs/heads/main','refs/heads/'+BRANCH]:check('remote:'+ref,remote.get(ref)==expected)
 check('local_heads',git('rev-parse','main')==git('rev-parse',BRANCH)==git('rev-parse','origin/main')==expected)
 check('exact_source',git('rev-parse',expected+':Source')==SOURCE and not git('diff',TESTED,expected,'--','Source','CMakeLists.txt'))
 ancestor(HANDOVER,checkpoint);ancestor(checkpoint,expected)
 for name in ['Docs/ProjectIndex/100_GainEntryAndBoundaryNavigation.md','Docs/ProjectIndex/101_MeasuredLRALanesAndIndependentHeights.md','Tests/Qualification/LraTimeline/results/2026-11-24/combined.json']:
  check('handover:'+name,subprocess.check_output(['git','show',HANDOVER+':'+name],cwd=ROOT)==subprocess.check_output(['git','show',expected+':'+name],cwd=ROOT))
 finder=[]
 for state in before['worktrees']:
  p=Path(state['path']);check('head:'+str(p),git('rev-parse','HEAD',cwd=p)==(expected if p in [ROOT,INTEGRATION] else state['head']))
  check('branch:'+str(p),git('branch','--show-current',cwd=p)==state['branch'])
  check('edits:'+str(p),hashlib.sha256(subprocess.check_output(['git','diff','HEAD','--binary'],cwd=p)).hexdigest()==state['tracked_diff_sha256'])
  for name,digest in state['untracked'].items():
   actual=sha(p/name) if (p/name).is_file() else None
   if actual!=digest and Path(name).name=='.DS_Store':finder.append(str(p/name))
   else:check('untracked:'+str(p/name),actual==digest)
 save(WORK/(label+'-verification.json'),{'passed':all(checks.values()),'stage':label,'main':expected,'checkpoint':checkpoint,'tag':TAG,'tag_object':obj,'branch':BRANCH,'source_tree':SOURCE,'local_candidates_preserved':len(before['local_candidates']),'protected_files':len(before['protected']),'old_tags_preserved':len(before['tags']),'manifest_entries_verified':before['manifest_entries_verified'],'finder_metadata_changes':finder,'checks':checks,'runtime_tests_repeated':False,'rebuild_performed':False,'physical_devices_opened':False,'audible_tests_run':False,'golem_accessed':False})
 print(label,'PASS:',len(checks),'checks; main',expected,'; annotated tag',obj)
if __name__=='__main__':
 if sys.argv[1]=='preflight':preflight()
 else:verify(sys.argv[1],sys.argv[2])
