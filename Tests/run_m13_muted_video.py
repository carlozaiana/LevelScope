#!/usr/bin/env python3
"""Four explicitly coordinated, all-muted physical cases. No endurance/retry loop."""
import argparse,datetime,hashlib,json,pathlib,platform,subprocess,time
from collect_m1_environment import sha256
from m7_video_acceptance import evaluate,POLICY
from VideoPaint.run_comparison import audio_snapshot,check_device

def measured_audio_identity(model):
    source=model['source']
    key='measured_programme_definition_digest' if 'composition' in model else 'measured_content_sha256'
    value=source.get(key)
    if not source.get('measurement_available') or not isinstance(value,str) or len(value)!=64 or any(c not in '0123456789abcdef' for c in value):
        raise ValueError('Missing verified measured audio identity')
    if 'composition' in model and value!=model['composition']['definition_digest']:
        raise ValueError('Composition measurement definition mismatch')
    return value

def validate_physical(root,identity,arch):
    summary=json.loads((root/'summary.json').read_text());record=json.loads((root/'identity.json').read_text());rows=json.loads((root/'runs.json').read_text())
    wanted={f'{codec}-{domain}' for codec in ('h264-vfr','prores-422') for domain in ('source','current')}
    binary='LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
    if (summary.get('accepted') is not True or summary.get('source_digest')!=identity['source_digest'] or summary.get('audio_configuration_unchanged') is not True
        or record.get('binary_sha256')!=identity['artifacts_sha256'][binary] or record.get('architecture')!=arch or len(rows)!=4 or {r['name'] for r in rows}!=wanted):raise ValueError('Missing or foreign physical subset')
    for row in rows:
        file=root/(row['name']+'.json')
        if sha256(file)!=row['report_sha256']:raise ValueError('Physical report changed')
        data=json.loads(file.read_text());result=evaluate(data,row['expected'],row['strict_exit_code'],row['fps'],channels=2,layout='stereo')
        if not result['accepted'] or not row['accepted'] or row['expected']['measurement_build']!=identity['measurement_build'] or row['expected']['architecture']!=arch:raise ValueError('Physical audio/timing/identity requirement failed')
        if row['name'].endswith('-current') and not data.get('section_plan_key'):raise ValueError('Current physical case omitted section plan')
    return summary

def main():
    p=argparse.ArgumentParser();p.add_argument('--candidate',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--platform',choices=('sequoia','high-sierra'),required=True);p.add_argument('--case-folder',type=pathlib.Path,required=True);p.add_argument('--allow-audio',action='store_true');a=p.parse_args()
    if not a.allow_audio:p.error('Coordinated --allow-audio is required')
    arch,prefix,device,rate=('arm64','15.','Mac Studio Speakers',48000) if a.platform=='sequoia' else ('x86_64','10.13.','Built-in Output',44100)
    osv=subprocess.check_output(['/usr/bin/sw_vers','-productVersion'],text=True).strip()
    if platform.machine()!=arch or not osv.startswith(prefix):raise ValueError('Actual platform required')
    candidate=a.candidate.resolve();identity=json.loads((candidate/'identity.json').read_text())
    for name,digest in identity['artifacts_sha256'].items():
        if sha256(candidate/name)!=digest:raise ValueError('Candidate changed: '+name)
    folder=a.case_folder.resolve();session=folder/'applied.levelscope';model=json.loads((folder/'measurement-only.json').read_text())
    if model['build']!=identity['measurement_build'] or model['sample_rate']!=rate or model['ordered_input_roles']!=['L','R'] or model['current']['evaluation']['verdict']!='PASS':raise ValueError('Matching measured stereo session required')
    plan=model['current_processing_plan'];gain=plan['global_gain_db'];audio_identity=measured_audio_identity(model);key=plan['audio_key']
    if not key or not session.is_file():raise ValueError('Applied section session required')
    before=audio_snapshot();check_device(before,dict(device=device,rate=rate));out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    binary=candidate/'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
    record=dict(policy=POLICY,source_digest=identity['source_digest'],identity_sha256=sha256(candidate/'identity.json'),binary_sha256=sha256(binary),architecture=arch,os=osv,
        cpu=subprocess.check_output(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),device=device,sample_rate=rate,all_routes_muted=True,
        session_sha256=sha256(session),measurement_report_sha256=sha256(folder/'measurement-only.json'),started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),audio_configuration_before_sha256=hashlib.sha256(json.dumps(before,sort_keys=True).encode()).hexdigest())
    (out/'identity.json').write_text(json.dumps(record,indent=2)+'\n');runs=[]
    for domain in ('source','current'):
        for codec,movie,fps in (('h264-vfr','h264-vfr.mov',None),('prores-422','422-25-1.mov',25)):
            name=codec+'-'+domain;video=candidate/'video-fixtures'/movie;report=out/(name+'.json')
            argv=[str(binary),'--allow-audio','--muted','--device',device,'--programme',str(session),'--video',str(video),'--seconds','8','--report',str(report),'--audition-state',domain]
            if domain=='current':argv+=['--session-processing']
            print('START',name,'8 seconds; every route muted',flush=True);start=time.monotonic()
            with (out/(name+'.log')).open('x') as log:run=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=428)
            data=json.loads(report.read_text());expected=dict(architecture=arch,device=device,sample_rate=rate,audition_state='Source' if domain=='source' else 'Current State',processing_gain_db=0. if domain=='source' else gain,section_plan_key='' if domain=='source' else key,measurement_build=identity['measurement_build'],video_sha256=sha256(video),audio_sha256=audio_identity,requested_seconds=8,programme_assembly='aligned-two-component-sum-v1' if 'composition' in model else 'interleaved')
            if 'composition' in model:expected['composition_identity']=model['composition']['definition_digest']
            result=evaluate(data,expected,run.returncode,fps,channels=2,layout='stereo');runs.append(dict(name=name,argv=argv,seconds=time.monotonic()-start,report_sha256=sha256(report),expected=expected,fps=fps,**result));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
            print('FINISH',name,'ACCEPTED' if result['accepted'] else 'FAILED',result['failures'],'skips',data.get('actual_paint_skipped_frames'),flush=True)
            if not result['accepted']:break
        if runs and not runs[-1]['accepted']:break
    unchanged=before==audio_snapshot();passed=len(runs)==4 and all(r['accepted'] for r in runs) and unchanged
    (out/'summary.json').write_text(json.dumps(dict(accepted=passed,source_digest=identity['source_digest'],policy=POLICY,expected_cases=4,completed_cases=len(runs),audio_configuration_unchanged=unchanged,physical_multichannel_listening='Not assessed; every route muted on two-output device'),indent=2)+'\n')
    if passed:validate_physical(out,identity,arch)
    return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
