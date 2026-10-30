"""Copy selected textual M17 receipts; never package audio/binaries/connection data."""
from pathlib import Path
import hashlib,json,shutil
p=Path(__file__).resolve().parent
out=p.parents[1]/'Tests/Reference/results/m17-qualification-2026-10-30'
out.mkdir(parents=True,exist_ok=True)
files=[]
for pattern in ['preflight.json','preservation.json','remote-preservation.txt','reuse-audit.json','cross-build-candidate-*.json','candidate-*/identity.json','sequoia-candidate-*/*.json','sequoia-candidate-*/*.log','sequoia-candidate-*/marker-files/marker-tests.json','ui-candidate-*/observations.json','ui-candidate-*/independent-report-audit.json','ui-candidate-*/session-audit.json','ui-candidate-*/report-*/report.json','ui-candidate-*/failed-job-report/report.json','ui-candidate-*/verified-job-export/report.json','final-*.json','dev-*.json','dev-*.log','runner-development*.log','configure*.log','development-build-*.log','marker-assertion-build.log','complete_candidate02.py','collect_evidence.py']:
 for f in sorted(p.glob(pattern)):
  if f.is_file() and f not in files:files.append(f)
index={}
for f in files:
 relative=f.relative_to(p);t=out/relative;t.parent.mkdir(parents=True,exist_ok=True);digest=hashlib.sha256(f.read_bytes()).hexdigest()
 if t.exists():assert hashlib.sha256(t.read_bytes()).hexdigest()==digest, 'Never overwrite an earlier receipt '+str(t)
 else:shutil.copy2(f,t)
 index[str(relative)]={'sha256':digest,'bytes':f.stat().st_size}
(out/'evidence-index.json').write_text(json.dumps(index,indent=2)+'\n')
print(len(index),'text receipts',sum(x['bytes'] for x in index.values()),'bytes')
