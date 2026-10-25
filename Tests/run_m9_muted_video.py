#!/usr/bin/env python3
"""Approved M9 bounded physical block: ten short cases and one 30-minute run; no retries."""
import argparse
import datetime
import hashlib
import json
import platform
import struct
import subprocess
import time
from pathlib import Path
from collect_m1_environment import sha256
from m7_video_acceptance import evaluate, POLICY
from VideoPaint.run_comparison import audio_snapshot, check_device


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--platform',choices=('sequoia','high-sierra'),required=True)
    parser.add_argument('--phase',choices=('short','endurance'),required=True)
    parser.add_argument('--allow-audio',action='store_true')
    args=parser.parse_args()
    if not args.allow_audio:parser.error('Coordinated --allow-audio is required')
    arch,prefix,device,rate=('arm64','15.','Mac Studio Speakers',48000) if args.platform=='sequoia' else ('x86_64','10.13.','Built-in Output',44100)
    os_version=subprocess.check_output(['/usr/bin/sw_vers','-productVersion'],universal_newlines=True).strip()
    if platform.machine()!=arch or not os_version.startswith(prefix):raise RuntimeError('Actual platform guard failed')
    candidate=args.candidate.resolve();output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    identity=json.loads((candidate/'identity.json').read_text())
    for name,expected in identity['artifacts_sha256'].items():
        if sha256(candidate/name)!=expected:raise RuntimeError('Frozen candidate changed: '+name)
    programme=output/'programme.levelscope'
    command=[str(candidate/'bin/LevelScopeProgrammeTests'),'--only','--prepare-session',str(candidate/('m9-fixtures/clock-'+str(rate)+'/programme.levelscope')),'--destination',str(programme),'--output',str(output/'preparation.json')]
    with (output/'preparation.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    with programme.open('rb') as stream:
        header=stream.read(80);core=json.loads(stream.read(struct.unpack('<Q',header[8:16])[0]))
    group_identity=core['programme']['assembly']['confirmation']
    before=audio_snapshot();check_device(before,dict(device=device,rate=rate))
    binary=candidate/'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
    videos=candidate/'m9-fixtures/physical-video'
    if args.phase=='short':
        cases=[(name+'-'+domain, videos/movie,18,domain,fps) for domain in ('source','current') for name,movie,fps in (
            ('h264-cfr','h264-cfr.mov',30000/1001),('h264-vfr','h264-vfr.mov',None),
            ('prores-proxy','proxy-25-1.mov',25),('prores-lt','lt-25-1.mov',25),('prores-422','422-25-1.mov',25))]
    else:cases=[('prores-current-30min',videos/'proxy-1800.mov',1800,'current',25)]
    record=dict(policy=POLICY,source_digest=identity['source_digest'],identity_sha256=sha256(candidate/'identity.json'),
        binary_sha256=sha256(binary),architecture=arch,os=subprocess.check_output(['/usr/bin/sw_vers'],universal_newlines=True),
        cpu=subprocess.check_output(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string'],universal_newlines=True).strip(),
        device=device,sample_rate=rate,all_routes_muted=True,phase=args.phase,script_sha256=sha256(Path(__file__)),programme_identity=group_identity,
        started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),audio_configuration_before_sha256=hashlib.sha256(json.dumps(before,sort_keys=True).encode()).hexdigest())
    (output/'identity.json').write_text(json.dumps(record,indent=2)+'\n');runs=[]
    for name,movie,seconds,domain,fps in cases:
        report=output/(name+'.json');gain=0. if domain=='source' else -3.305074118944
        argv=[str(binary),'--allow-audio','--muted','--device',device,'--programme',str(programme),'--video',str(movie),'--seconds',str(seconds),'--report',str(report),'--audition-state',domain]
        if domain=='current':argv+=['--processing-gain',str(gain)]
        print('START',name,seconds,'seconds; all routes muted',flush=True);start=time.monotonic()
        with (output/(name+'.log')).open('x') as log:run=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=seconds+420)
        data=json.loads(report.read_text());expected=dict(architecture=arch,device=device,sample_rate=rate,audition_state='Source' if domain=='source' else 'Current State',processing_gain_db=gain,measurement_build=identity['measurement_build'],video_sha256=sha256(movie),audio_sha256=group_identity,requested_seconds=seconds,programme_assembly='mono-channel-assembly-v2',alignment_confirmation=group_identity)
        result=evaluate(data,expected,run.returncode,fps,channels=6,layout='mono-channel-assembly-v2')
        runs.append(dict(name=name,argv=argv,seconds=time.monotonic()-start,report_sha256=sha256(report),expected=expected,fps=fps,**result));(output/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
        print('FINISH',name,'ACCEPTED' if result['accepted'] else 'FAILED',result['failures'],'skips',data.get('actual_paint_skipped_frames'),flush=True)
        if not result['accepted']:break
    unchanged=before==audio_snapshot();passed=len(runs)==len(cases) and all(r['accepted'] for r in runs) and unchanged
    (output/'summary.json').write_text(json.dumps(dict(accepted=passed,policy=POLICY,expected_cases=len(cases),completed_cases=len(runs),audio_configuration_unchanged=unchanged,physical_multichannel_listening='Not assessed; every route muted on two-output device'),indent=2)+'\n')
    return 0 if passed else 1


if __name__=='__main__':raise SystemExit(main())
