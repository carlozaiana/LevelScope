from pathlib import Path
import subprocess,sys,json,time,ast
R=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation');B=R/'build/prores-20261020';P=B/'candidate-04';E=B/'evidence'
# Candidate03's native and regression binaries are identical to candidate04.
# Let that already-running sequence finish before timed UI playback.
end=time.monotonic()+1800
while not (E/'sequoia-regressions-candidate03/summary.json').exists():
 assert time.monotonic()<end,'Previous regression sequence did not finish'
 time.sleep(2)
assert json.loads((E/'sequoia-native-candidate03/summary.json').read_text())['passed']
assert json.loads((E/'sequoia-regressions-candidate03/summary.json').read_text())['automated_passed']
# Load only the independent verifier function; do not execute its CLI.
module=ast.parse((P/'Tests/verify_prores_playback_evidence.py').read_text()); ns={}
exec('import math',ns)
exec(compile(ast.Module(body=[x for x in module.body if isinstance(x,ast.FunctionDef) and x.name=='verify_memory'],type_ignores=[]),'memory-verifier','exec'),ns)
exe=P/'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
report=E/'sequoia-memory-v2-300s-candidate04.json';assert not report.exists()
cmd=[str(exe),'--allow-audio','--muted','--device','Mac Studio Speakers','--seconds','300','--audio',str(B/'fixtures/selected-clock-48000.wav'),'--video',str(B/'fixtures/proxy-1800.mov'),'--report',str(report)]
print('START final candidate04 stabilization',flush=True)
with (B/'stabilization-candidate04.log').open('w') as log:
 rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,timeout=750).returncode
r=json.loads(report.read_text());errors=ns['verify_memory'](r)
print('Stabilization exit',rc,'independent memory errors',errors,flush=True)
assert rc==0 and r['passed'] and not errors, 'Final stabilization failed'
for phase in ['matrix','endurance']:
 print('START',phase,'candidate04',flush=True)
 with (B/(phase+'-candidate04.log')).open('w') as log:
  rc=subprocess.run([sys.executable,'-B',str(P/'Tests/run_prores_playback_qualification.py'),'--binary',str(exe),'--fixtures',str(B/'fixtures'),'--actual','/Users/fuzz/Developer/LevelScope/Tests/51_SyncTest.mov','--output',str(E/('sequoia-'+phase+'-candidate04')),'--device','Mac Studio Speakers','--sample-rate','48000','--phase',phase],stdout=log,stderr=subprocess.STDOUT).returncode
 print('FINISH',phase,'exit',rc,flush=True)
 assert rc==0,phase+' failed'
 for p in (E/('sequoia-'+phase+'-candidate04')).glob('*.json'):
  rr=json.loads(p.read_text())
  if 'memory_method' in rr:assert not ns['verify_memory'](rr),(p,ns['verify_memory'](rr))
print('Final Sequoia playback sequence completed',flush=True)
