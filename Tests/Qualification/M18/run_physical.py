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
        d=json.loads(report.read_text())
        # Preserve the old zero-picture-skip verdict separately. M18 allows
        # occasional picture skips, subject to review of the recorded intervals.
        # Never mask audio, clock, attachment, memory or processing failures.
        core=(r.returncode in (0,1) and not d['error'] and d['memory_passed']
              and d['attachment_within_budget'] and d['bound_failures']==0
              and d['blank_paints']==0 and d['underrun_blocks']==0
              and d['source_seconds']>=a.seconds and d['paints']>0
              and d['pre_monitor_channel_count_valid'] and d['pre_monitor_frames_observed']>0
              and d['all_routes_muted'] and d['hal_output_timestamp']
              and d['device']==a.device and d['sample_rate']==48000
              and d['measurement_build']==identity['measurement_build']
              and d['qualification_eligible'] and d['selected_audio_channels']==2)
        review=d['actual_paint_skipped_frames']>0 or d['paint_deadline_failures']>0
        rows.append(dict(case=key,command=cmd,exit_code=r.returncode,seconds=time.monotonic()-t,
                         report_sha256=sha(report),strict_legacy_passed=d['passed'],
                         non_picture_checks_passed=core,picture_review_required=review))
        (out/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
passed=all(r['non_picture_checks_passed'] and r['strict_legacy_passed'] for r in rows)
summary=dict(passed=passed,non_picture_checks_passed=all(r['non_picture_checks_passed'] for r in rows),
             picture_review_required=any(r['picture_review_required'] for r in rows),
             device=a.device,all_routes_muted=True,build=identity['measurement_build'],
             identity_sha256=sha(root/'identity.json'),cases=rows)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2));raise SystemExit(0 if passed else 1)
