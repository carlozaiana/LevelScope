#!/usr/bin/env python3
"""M6 candidate-01 orchestration only; explicitly coordinated muted devices."""
import argparse, hashlib, json, pathlib, subprocess, sys, time
p=argparse.ArgumentParser()
for key in ('root','fixtures','actual','vfr','device','sample-rate'):p.add_argument('--'+key,required=True)
a=p.parse_args();root=pathlib.Path(a.root);candidate=root/'candidate-01';fixtures=pathlib.Path(a.fixtures)
binary=candidate/'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
runs=[]
def run(name,cmd):
 start=time.monotonic();print('START',name,flush=True)
 with (root/(name+'.log')).open('x') as log:rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
 runs.append(dict(name=name,argv=list(map(str,cmd)),exit_code=rc,seconds=time.monotonic()-start))
 (root/'muted-video-driver-runs.json').write_text(json.dumps(runs,indent=2)+'\n')
 print('FINISH',name,rc,flush=True)
 if rc:raise SystemExit(rc)
for domain in ('source','current'):
 cmd=[sys.executable,'-B',str(candidate/'Tests/run_prores_playback_qualification.py'),'--binary',str(binary),'--fixtures',a.fixtures,'--actual',a.actual,'--output',str(root/('prores-matrix-'+domain+'-01')),'--device',a.device,'--sample-rate',a.sample_rate,'--phase','matrix','--audition-state',domain]
 if domain=='current':cmd+=['--processing-gain','-3.305074118944']
 run('prores-matrix-'+domain+'-01',cmd)
for domain in ('source','current'):
 for kind,video in [('cfr',candidate/'video-fixtures/half-hour.mov'),('vfr',pathlib.Path(a.vfr))]:
  name='h264-'+kind+'-'+domain+'-01';out=root/(name+'.json')
  cmd=[str(binary),'--allow-audio','--muted','--device',a.device,'--audio',str(fixtures/('selected-clock-'+a.sample_rate+'.wav')),'--video',str(video),'--seconds','18','--report',str(out),'--audition-state',domain]
  if domain=='current':cmd+=['--processing-gain','-3.305074118944']
  run(name,cmd)
  if not json.loads(out.read_text())['passed']:raise SystemExit('H264 report failed')
cmd=[sys.executable,'-B',str(candidate/'Tests/run_prores_playback_qualification.py'),'--binary',str(binary),'--fixtures',a.fixtures,'--actual',a.actual,'--output',str(root/'prores-endurance-current-01'),'--device',a.device,'--sample-rate',a.sample_rate,'--phase','endurance','--audition-state','current','--processing-gain','-3.305074118944']
run('prores-endurance-current-01',cmd)
(root/'muted-video-driver-summary.json').write_text(json.dumps(dict(passed=True,script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),runs=runs),indent=2)+'\n')
