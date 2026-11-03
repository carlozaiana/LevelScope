"""Four explicitly coordinated, muted physical video cases. Never chooses a default."""
import argparse
import json
from pathlib import Path
import subprocess
import time
from run_short import sha
p=argparse.ArgumentParser();p.add_argument('candidate',type=Path);p.add_argument('session',type=Path);p.add_argument('video',type=Path);p.add_argument('output',type=Path)
p.add_argument('--device',required=True,choices=['Mac Studio Speakers','Built-in Output']);p.add_argument('--coordinated',action='store_true',required=True);p.add_argument('--seconds',type=int,default=10,choices=range(5,61));a=p.parse_args()
root=a.candidate.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);identity=json.loads((root/'identity.json').read_text())
for f,h in identity['artifacts_sha256'].items():assert sha(root/f)==h
binary=root/'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification';rows=[]
for name,file in [('h264','h264-cfr.mov'),('prores','422-25-1.mov')]:
    for domain in ('source','current'):
        key=name+'-'+domain;report=out/(key+'.json')
        cmd=[str(binary),'--allow-audio','--muted','--device',a.device,'--programme',str(a.session.resolve()),'--audition-state',domain,'--video',str((a.video/file).resolve()),'--seconds',str(a.seconds),'--report',str(report)]
        if domain=='current':cmd.append('--session-processing')
        t=time.monotonic()
        with (out/(key+'.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
        d=json.loads(report.read_text());rows.append(dict(case=key,command=cmd,exit_code=r.returncode,seconds=time.monotonic()-t,report_sha256=sha(report),passed=d['passed']))
        (out/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
        if r.returncode or not d['passed']:raise RuntimeError('Physical case failed; preserve report and review '+key)
(out/'summary.json').write_text(json.dumps(dict(passed=True,device=a.device,all_routes_muted=True,build=identity['measurement_build'],identity_sha256=sha(root/'identity.json'),cases=rows),indent=2)+'\n')
