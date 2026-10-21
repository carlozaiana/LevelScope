from pathlib import Path
import sys,subprocess,json,hashlib,datetime,platform
which=sys.argv[1]
assert which in ('sequoia','high-sierra')
if which=='sequoia':
 root=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020')
 out=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-integration-20261021/sequoia')
 actual=Path('/Users/fuzz/Developer/LevelScope/Tests/51_SyncTest.mov')
else:
 root=Path('/Users/fuzz/LevelScope-ProRes-qualification-2026-10-20')
 out=Path('/Users/fuzz/LevelScope-ProRes-integration-2026-10-21/high-sierra')
 actual=root/'51_SyncTest.mov'
p=root/'candidate-04';out.mkdir(parents=True,exist_ok=False)
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as stream:
  for b in iter(lambda:stream.read(4194304),b''):h.update(b)
 return h.hexdigest()
i=json.loads((p/'identity.json').read_text())
bad=[n for n,h in i['artifacts_sha256'].items() if sha(p/n)!=h]
report=dict(start_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),platform=which,machine=platform.machine(),candidate=str(p),identity_sha256=sha(p/'identity.json'),verified_artifacts=len(i['artifacts_sha256']),artifact_mismatches=bad,actual_sha256=sha(actual),fixture_manifest_sha256=sha(root/'fixtures/manifest.json'),audio_devices_opened=False,runs=[])
assert not bad and report['identity_sha256']=='45219de9006b845230c5e1296859bff2c483beb49c4b827b8c33c552f8381285'
for name,args in [
 ('native',[str(p/'Tests/run_prores_native_qualification.py'),'--binary',str(p/'bin/LevelScopeVideoTests'),'--fixtures',str(root/'fixtures'),'--actual',str(actual),'--output',str(out/'native')]),
 ('regressions',[str(p/'Tests/run_m5_qualification.py'),'--platform',which,'--output',str(out/'regressions')])]:
 print('START',which,name,flush=True)
 with (out/(name+'.log')).open('w') as f:
  rc=subprocess.run([sys.executable,'-B',*args],cwd=p,stdout=f,stderr=subprocess.STDOUT).returncode
 report['runs'].append(dict(name=name,argv=args,exit_code=rc))
 (out/'integration.json').write_text(json.dumps(report,indent=2)+'\n')
 print('FINISH',which,name,rc,flush=True)
 assert rc==0,name+' failed'
report['artifact_mismatches_after']=[n for n,h in i['artifacts_sha256'].items() if sha(p/n)!=h]
assert not report['artifact_mismatches_after']
report['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();report['passed']=True
(out/'integration.json').write_text(json.dumps(report,indent=2)+'\n')
print(which,'integration regressions PASS',flush=True)
