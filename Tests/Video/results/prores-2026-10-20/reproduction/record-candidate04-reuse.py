from pathlib import Path
import hashlib,json
R=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation');B=R/'build/prores-20261020';dirs=[B/'candidate-03',B/'candidate-04'];ids=[json.loads((p/'identity.json').read_text()) for p in dirs]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for p,i in zip(dirs,ids):
 for name,want in i['artifacts_sha256'].items():assert sha(p/name)==want,name
source_diff=sorted(k for k in ids[0]['source_files_sha256'].keys()|ids[1]['source_files_sha256'].keys() if ids[0]['source_files_sha256'].get(k)!=ids[1]['source_files_sha256'].get(k))
assert source_diff==['Tests/VideoPlaybackQualification.cpp','Tests/run_prores_playback_qualification.py','Tests/verify_prores_playback_evidence.py']
unchanged={k:v for k,v in ids[0]['artifacts_sha256'].items() if k.startswith('bin/') or k.endswith('/MacOS/LevelScope Standalone')}
assert len(unchanged)==8 and all(ids[1]['artifacts_sha256'][k]==v for k,v in unchanged.items())
evidence=['sequoia-native-candidate03','sequoia-regressions-candidate03']
for scope in evidence:
 summary=json.loads((B/'evidence'/scope/'summary.json').read_text());assert summary.get('passed',summary.get('automated_passed'))
ui=json.loads((B/'evidence/product-ui/sequoia-candidate03.json').read_text());assert ui['passed'] and ui['binary_sha256']==unchanged['LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone']
files={}
for scope in evidence:
 for p in sorted((B/'evidence'/scope).rglob('*')):
  if p.is_file():files[str(p.relative_to(B/'evidence'))]=sha(p)
files['product-ui/sequoia-candidate03.json']=sha(B/'evidence/product-ui/sequoia-candidate03.json')
result=dict(passed=True,scope='Exact candidate03 native/M1-M5/UI binary evidence reused by candidate04; no playback-harness evidence reused',candidate_identities=[sha(p/'identity.json') for p in dirs],production_source_unchanged=True,source_differences=source_diff,unchanged_executables=unchanged,evidence_sha256=files,required_fresh_playback='All matrix and six endurance runs use final candidate04 harness')
p=B/'evidence/candidate04-native-regression-ui-reuse.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
print('PASS verified reuse of native, M1-M5 and Sequoia UI evidence with unchanged production and 8 identical binaries')
