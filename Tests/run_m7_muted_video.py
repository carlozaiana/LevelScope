#!/usr/bin/env python3
"""M7 named-device timing checks, only after explicit user approval. All routes muted."""
import argparse,json,pathlib,subprocess,sys,time,hashlib
p=argparse.ArgumentParser(description=__doc__)
for name in ('candidate','fixtures','actual','vfr','output'):p.add_argument('--'+name,type=pathlib.Path,required=True)
p.add_argument('--device',required=True);p.add_argument('--sample-rate',type=int,choices=(44100,48000),required=True)
p.add_argument('--allow-audio',action='store_true');a=p.parse_args()
if not a.allow_audio:raise SystemExit('Explicit coordinated --allow-audio is required')
a.output.mkdir(parents=True,exist_ok=False)
binary=a.candidate/'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
audio=a.candidate/('m7-fixtures/clock-71-'+str(a.sample_rate)+'.flac')
runs=[]
def run(name,cmd):
    start=time.monotonic();print('START',name,flush=True)
    with (a.output/(name+'.log')).open('x') as log:rc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT).returncode
    runs.append(dict(name=name,argv=list(map(str,cmd)),exit_code=rc,seconds=time.monotonic()-start))
    (a.output/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
    print('FINISH',name,rc,flush=True)
    if rc:raise SystemExit(rc)
for domain in ('source','current'):
    cmd=[sys.executable,'-B',str(a.candidate/'Tests/run_prores_playback_qualification.py'),'--binary',str(binary),'--fixtures',str(a.fixtures),'--actual',str(a.actual),'--output',str(a.output/('prores-matrix-'+domain)),'--device',a.device,'--sample-rate',str(a.sample_rate),'--phase','matrix','--audition-state',domain,'--audio',str(audio)]
    if domain=='current':cmd+=['--processing-gain','-3.305074118944']
    run('prores-matrix-'+domain,cmd)
for domain in ('source','current'):
    for kind,movie in [('cfr',a.candidate/'video-fixtures/half-hour.mov'),('vfr',a.vfr)]:
        name='h264-'+kind+'-'+domain;report=a.output/(name+'.json')
        cmd=[str(binary),'--allow-audio','--muted','--device',a.device,'--audio',str(audio),'--video',str(movie),'--seconds','18','--report',str(report),'--audition-state',domain]
        if domain=='current':cmd+=['--processing-gain','-3.305074118944']
        run(name,cmd)
        data=json.loads(report.read_text())
        if not data['passed'] or data['selected_audio_channels']!=8 or data['pre_monitor_nonzero_channel_mask']!=255:raise SystemExit('H264 eight-channel evidence failed')
cmd=[sys.executable,'-B',str(a.candidate/'Tests/run_prores_playback_qualification.py'),'--binary',str(binary),'--fixtures',str(a.fixtures),'--actual',str(a.actual),'--output',str(a.output/'prores-endurance-current'),'--device',a.device,'--sample-rate',str(a.sample_rate),'--phase','endurance','--audition-state','current','--processing-gain','-3.305074118944','--audio',str(audio)]
run('prores-endurance-current',cmd)
(a.output/'summary.json').write_text(json.dumps(dict(passed=True,physical_71_listening='Not assessed; all routes muted on two-output device',script_sha256=hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),runs=runs),indent=2)+'\n')
