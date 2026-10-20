from pathlib import Path
import subprocess,sys,json,ast
R=Path('/Users/fuzz/LevelScope-ProRes-qualification-2026-10-20');P=R/'candidate-04';E=R/'evidence'
assert json.loads((R/'candidate04-verification.json').read_text())['passed']
assert json.loads((E/'product-ui/high-sierra-candidate04.json').read_text())['passed']
for stage,cmd in [
 ('native',[sys.executable,'-B',str(P/'Tests/run_prores_native_qualification.py'),'--binary',str(P/'bin/LevelScopeVideoTests'),'--fixtures',str(R/'fixtures'),'--actual',str(R/'51_SyncTest.mov'),'--output',str(E/'high-sierra-native-candidate04'),'--loading']),
 ('regressions',[sys.executable,'-B',str(P/'Tests/run_m5_qualification.py'),'--platform','high-sierra','--output',str(E/'high-sierra-regressions-candidate04')])]:
 print('START',stage,flush=True)
 with (R/(stage+'-candidate04.log')).open('w') as log:rc=subprocess.run(cmd,cwd=P,stdout=log,stderr=subprocess.STDOUT).returncode
 print('FINISH',stage,rc,flush=True)
 assert rc==0,stage+' failed'
for phase in ['matrix','endurance']:
 print('START',phase,flush=True)
 with (R/(phase+'-candidate04.log')).open('w') as log:
  rc=subprocess.run([sys.executable,'-B',str(P/'Tests/run_prores_playback_qualification.py'),'--binary',str(P/'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'),'--fixtures',str(R/'fixtures'),'--actual',str(R/'51_SyncTest.mov'),'--output',str(E/('high-sierra-'+phase+'-candidate04')),'--device','Built-in Output','--sample-rate','44100','--phase',phase],stdout=log,stderr=subprocess.STDOUT).returncode
 print('FINISH',phase,rc,flush=True)
 assert rc==0,phase+' failed'
print('Golem final qualification sequence completed',flush=True)
