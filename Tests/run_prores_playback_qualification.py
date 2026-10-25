#!/usr/bin/env python3
"""Explicitly muted, named-device ProRes display qualification; no default-device selection."""
from pathlib import Path
import argparse, subprocess, hashlib, json, platform, time
p=argparse.ArgumentParser()
for name in ['binary','fixtures','actual','output']:p.add_argument('--'+name,type=Path,required=True)
p.add_argument('--device',required=True);p.add_argument('--sample-rate',type=int,choices=[44100,48000],required=True)
p.add_argument('--phase',choices=['matrix','endurance'],required=True)
p.add_argument('--audition-state',choices=['source','current'],default='source');p.add_argument('--processing-gain',type=float)
p.add_argument('--audio',type=Path,help='Explicit selected PCM source; defaults to the preserved stereo clock fixture')
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
identity=dict(binary_sha256=sha(a.binary),architecture=platform.machine(),os=subprocess.check_output(['/usr/bin/sw_vers'],text=True),
              fixture_manifest_sha256=sha(a.fixtures/'manifest.json'),actual_sha256=sha(a.actual),device=a.device,
              sample_rate=a.sample_rate,all_routes_muted=True,audition_state=a.audition_state,processing_gain_db=a.processing_gain)
audio=a.audio or a.fixtures/('selected-clock-'+str(a.sample_rate)+'.wav')
identity['selected_audio_sha256']=sha(audio)
(a.output/'identity.json').write_text(json.dumps(identity,indent=2)+'\n')
manifest=json.loads((a.fixtures/'manifest.json').read_text())
cases=[(a.fixtures/name,18) for name,m in manifest.items() if m.get('required')]
cases+=[(a.actual,100)]
if a.phase=='endurance':cases=[(a.fixtures/(profile+'-1800.mov'),1800) for profile in ['proxy','lt','422']]
results=[]
for movie,duration in cases:
 name=movie.stem;report=a.output/(name+'.json')
 cmd=[str(a.binary.resolve()),'--allow-audio','--muted','--device',a.device,'--seconds',str(duration),
      '--audio',str(audio.resolve()),
      '--video',str(movie.resolve()),'--report',str(report.resolve())]
 cmd += ['--audition-state',a.audition_state]
 if a.processing_gain is not None:cmd += ['--processing-gain',str(a.processing_gain)]
 start=time.monotonic()
 print('START',name,'requested_seconds',duration,flush=True)
 with (a.output/(name+'.log')).open('w') as log:
  proc=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
  while proc.poll() is None:
   if time.monotonic()-start>duration+300+movie.stat().st_size/(64*1024*1024):
    proc.terminate();proc.wait();raise RuntimeError('Playback watchdog timeout: '+name)
   time.sleep(1)
 rc=proc.wait();data=json.loads(report.read_text()) if report.exists() else {}
 result=dict(name=name,exit_code=rc,passed=rc==0 and data.get('passed',False) and data.get('qualification_eligible') is True
                 and data.get('memory_method')=='prores-memory-v2-user-approved-2026-10-20',
             elapsed_seconds=time.monotonic()-start,argv=cmd)
 results.append(result);(a.output/'runs.json').write_text(json.dumps(results,indent=2)+'\n')
 print('FINISH',name,'PASS' if result['passed'] else 'FAIL','seconds',round(result['elapsed_seconds'],2),
       'skips',data.get('actual_paint_skipped_frames'),'peak_rss',data.get('peak_sampled_playback_rss_bytes'),
       'bound_failures',data.get('bound_failures'),flush=True)
ok=all(r['passed'] for r in results)
(a.output/'summary.json').write_text(json.dumps(dict(passed=ok,identity=identity,test_count=len(results)),indent=2)+'\n')
raise SystemExit(0 if ok else 1)
