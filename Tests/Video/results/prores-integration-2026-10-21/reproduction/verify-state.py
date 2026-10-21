from pathlib import Path
import json,subprocess,hashlib,sys
R=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation');O=Path('/Users/fuzz/Developer/LevelScope');B=R/'build/prores-integration-20261021';P=R/'build/prores-20261020/candidate-04'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4194304),b''):h.update(b)
 return h.hexdigest()
def git(*a,cwd=R):return subprocess.check_output(['git',*a],cwd=cwd,text=True).strip()
before=json.loads((B/'before.json').read_text());i=json.loads((P/'identity.json').read_text())
remote={v.split()[1]:v.split()[0] for v in git('ls-remote','origin','refs/heads/main','refs/heads/codex/video-implementation','refs/tags/*').splitlines()}
tags=[]
for row in before['tags']:
 n,o,c=row['ref'],row['object'],row['commit']
 tags.append(dict(row,unchanged=(git('rev-parse',n)==o==remote.get(n) and git('rev-parse',n+'^{commit}')==c==remote.get(n+'^{}'))))
original_bad=[n for n,h in before['original_inventory'].items() if not (O/n).is_file() or sha(O/n)!=h]
untracked_bad=[n for n,h in before['video_untracked'].items() if not (R/n).is_file() or sha(R/n)!=h]
source_bad=[n for n,h in i['source_files_sha256'].items() if sha(R/n)!=h]
art_bad=[n for n,h in i['artifacts_sha256'].items() if sha(P/n)!=h]
source_diff=git('diff','--name-only',i['base_commit'],'--','CMakeLists.txt','Source','Tests',':!Tests/Video/results',':!Tests/Video/README.md')
report=dict(main_before=before['main_before'],local_main=git('rev-parse','main'),remote_main=remote['refs/heads/main'],video_head=git('rev-parse','HEAD'),remote_video=remote['refs/heads/codex/video-implementation'],checkpoint_tags=tags,original_head=git('rev-parse','HEAD',cwd=O),original_branch=git('branch','--show-current',cwd=O),original_status=git('status','--porcelain',cwd=O),original_inventory_changes=original_bad,video_untracked_changes=untracked_bad,tested_source=i['base_commit'],source_count=len(i['source_files_sha256']),source_mismatches=source_bad,source_diff_since_candidate=source_diff,artifact_count=len(i['artifacts_sha256']),artifact_mismatches=art_bad,candidate_identity_sha256=sha(P/'identity.json'),worktrees=git('worktree','list','--porcelain'))
report['passed']=all(x['unchanged'] for x in tags) and not any([original_bad,untracked_bad,source_bad,source_diff,art_bad]) and report['original_head']==before['original_head'] and report['original_branch']==before['original_branch'] and report['original_status']==before['original_status']
out=B/sys.argv[1];assert not out.exists();out.write_text(json.dumps(report,indent=2)+'\n')
print('Preservation/source verification:',report['passed']);assert report['passed']
