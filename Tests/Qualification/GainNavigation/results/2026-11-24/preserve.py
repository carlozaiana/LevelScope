import hashlib,json,subprocess,sys
from pathlib import Path
root=Path(__file__).resolve().parent
prior=Path('/Users/fuzz/Developer/LevelScope/build/mono-lra-2026-11-24/publication-2026-11-24/preflight.json')
def run(*args):return subprocess.check_output(args,text=True).strip()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def refs():return {r:[run('git','rev-parse',r),run('git','rev-parse',r+'^{}')] for r in run('git','tag').splitlines()}
before=root/'preservation-before.json'
if not before.exists():
 d=json.loads(prior.read_text());files=d['protected']
 for p in [Path('Docs/ProjectIndex/97_MonoLRAProposals.md'),Path('Docs/ProjectIndex/98_MonoLRAProposalPublication.md')]:files[str(p.resolve())]=sha(p)
 snapshot={'main':run('git','rev-parse','main'),'tags':refs(),'protected':files,'prior_packages':d['prior_packages']}
else:snapshot=json.loads(before.read_text())
errors=[p for p,digest in snapshot['protected'].items() if not Path(p).is_file() or sha(p)!=digest]
assert not errors,errors
assert run('git','rev-parse','main')==snapshot['main']
assert refs()==snapshot['tags']
if not before.exists():before.write_text(json.dumps(snapshot,indent=2)+'\n')
result={'passed':True,'protected_files':len(snapshot['protected']),'tags':len(snapshot['tags']),'main':snapshot['main'],'prior_packages':snapshot['prior_packages']}
(root/sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
