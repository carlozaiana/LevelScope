import hashlib,json,subprocess
from pathlib import Path
root=Path(__file__).resolve().parent
repo=root.parent.parent
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=repo,text=True).strip()
before=json.loads((root/'preservation-before.json').read_text())
prior=json.loads((repo/'build/boundaries-workflow-publication-2026-11-23/preflight.json').read_text())
checks={}
for name,digest in prior['protected'].items():checks[name]=sha(Path(name))==digest
candidate=Path(before['candidate03'])
checks[str(candidate/'identity.json')]=sha(candidate/'identity.json')==before['identity_sha256']
for name,digest in before['artifacts'].items():checks[str(candidate/name)]=sha(candidate/name)==digest
checks['main']=git('rev-parse','main')==before['main']
checks['tags']=git('for-each-ref','--format=%(refname) %(objectname) %(*objectname)','refs/tags')==before['tags'].strip()
failed=[k for k,v in checks.items() if not v]
result={'passed':not failed,'previous_packages':len(prior['packages']),'previous_protected_files':len(prior['protected']),'shared_source_candidate03_artifacts':len(before['artifacts']),'checks':len(checks),'failed':failed,'main':git('rev-parse','main'),'tags':len(before['tags'].splitlines()),'scope':'Read-only; no remote actions. Current graph work only; original packages, evidence, media, sessions and manuals preserved.'}
(root/'preservation-after.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2));assert not failed
