#!/usr/bin/env python3
"""Explicit M11 impact-selected checks; no historical milestone runners."""
import argparse,json,math,pathlib,platform,shutil,subprocess,sys,time
from collections import Counter
from collect_m1_environment import sha256
from verify_m11_pcm import verify as verify_pcm

UNITS={'measurement':('LevelScopeMeasurementTests',None),'scaffold':('LevelScopeScaffoldTests',None),
 'compliance':('LevelScopeComplianceTests','--report'),'proposals':('LevelScopeProposalTests','--report'),
 'diagnostics':('LevelScopeDiagnosticsTests','--report'),'sections':('LevelScopeSectionsTests','--report'),
 'audition':('LevelScopeAuditionTests','--report'),'multichannel':('LevelScopeMultichannelTests','--report'),
 'sessions':('LevelScopeSessionTests','--output'),'programmes':('LevelScopeProgrammeTests','--output'),
 'export':('LevelScopeExportTests','--output'),'video-timing':('LevelScopeVideoTests',None)}
REQUIRED_SHORT=set(UNITS)|{'m11-core','m11-video','ebu','itu','ebu-runner','m11-runner','legacy-rejection','independent-pcm'}
REQUIRED={'short':REQUIRED_SHORT,'long':{'m11-long','independent-rf64'},'exchange':{'m11-exchange'}}
UI_CASES={'three-presets','confirmation-gate','component-roles-and-files','source-current-m3',
 'diagnostics-sections','session-reopen-relink','cancel-invalid-preserves-session','stopped-only-edits',
 'verified-export-report','device-closed-muted','older-session-save-as'}
MANDATORY_CHECKS={'invalid mapping purpose alignment or identity cannot replace valid session',
 'labels preserve audio revision and valid measurements','cancelled verification leaves applied gain unverified without rollback',
 'already-passes preserves measured Current and exact gain','missing one component member recovers settings but cannot retain a current pass',
 'old composition completion cannot become current','sum overload silences entire callback and preserves presented audio position',
 'unsupported or mismatched composition core rejected safely','damaged Current cache preserves composition settings and valid Source cache',
 'nonfinite component rejected before it can replace confirmed composition','active composition copy cancellation returns no partial session',
 'mixed WAV AIFF FLAC PCM16/24/32 composed samples equal independent decoded sum',
 'muting composed peak owner leaves measured rules and peak-constrained proposal unchanged'}

def validate_report(r,identity,arch,phase):
    if r.get('passed') is not True or r.get('build')!=identity['measurement_build'] or r.get('architecture')!=arch or r.get('phase')!=phase or r.get('physical_devices_opened') is not False:raise ValueError('Failed/mismatched M11 '+phase)
    checks=r.get('checks');
    if not isinstance(checks,list) or not checks or not all(x.get('passed') is True for x in checks):raise ValueError('Missing/failed assertions')
    if phase=='core':
        expected={(rate,p,group,gain) for rate in (44100,48000,96000,192000) for p in range(3) for group in ('interleaved','grouped') for gain in (0,3.125,-6.02059991328)}
        actual={(x['rate'],int(x['kind'].split('-')[1]),x['kind'].split('-')[2],x['gain_db']) for x in r.get('cases',[])}
        if actual!=expected or len(r['cases'])!=72:raise ValueError('Missing/duplicate composition matrix case')
        if not MANDATORY_CHECKS.issubset({x['name'] for x in checks}):raise ValueError('Missing required negative tests')
    if phase=='video':
        selections=r.get('video_selections',[])
        if r.get('video_cases')!=6 or len(selections)!=60:raise ValueError('Incomplete H264/ProRes cases')
        expected=Counter((file,domain,offset,frame) for file in ('cfr-exact.mov','vfr-exact.mov','422-with-tone.mov') for domain in ('Source','Current State') for offset in ('-1001/30000','7/3') for frame in ('256','72256','144256','256','192000'))
        actual=Counter((row['file'],row['domain'],row['offset'],row['audio_frame']) for row in selections)
        if actual!=expected:raise ValueError('Missing/duplicate native selection case')
        for row in selections:
            if not 0<=row['selection_ms']<10000 or row['player_skipped_frames']<0 or row['physical_paint_test'] is not False:raise ValueError('Invalid native selection evidence')
    if phase=='exchange':
        rows=r.get('exchange_cases',[])
        if len(rows)!=6 or {(x['preset'],x['representation']) for x in rows}!={(p,k) for p in range(3) for k in ('grouped','interleaved')}:raise ValueError('Incomplete cross-platform exchange inventory')
        if any(x['origin_platform']==x['receiving_platform'] or x['foreign_cache_reused'] is not False or x['remeasured_current_pass'] is not True for x in rows):raise ValueError('Invalid cross-platform exchange evidence')
    if phase=='long':
        v=r['long_composition']
        for k,cap in (('total_ms',10800000),('cancel_ms',1000),('ui_p95_ms',100),('rss_growth_bytes',256*1024*1024),('session_bytes',128*1024*1024),('report_bytes',16*1024*1024)):
            if type(v.get(k)) not in (int,float) or not math.isfinite(v[k]) or not 0<=v[k]<=cap:raise ValueError('Long resource budget '+k)
        for k in ('save_ms','parse_ms'):
            if not 0<=v['stages_ms'][k]<=5000:raise ValueError('Long persistence budget '+k)
        if (v['frames'],v['sample_rate'],v['channels'],v['components'],v['annotations'],v['capture_cases'],v['artifact_bytes'])!=('518400000',48000,6,2,1000,6,9331200104) or v['ui_callbacks']<10:raise ValueError('Incomplete real long workflow')

def validate_runs(runs,phase):
    names=[r['name'] for r in runs]
    if set(names)!=REQUIRED[phase] or len(names)!=len(REQUIRED[phase]) or any(r.get('exit_code')!=0 for r in runs):raise ValueError('Missing, duplicate, unexpected or failed required run')

def validate_ui(r,identity,arch):
    if r.get('passed') is not True or r.get('source_digest')!=identity['source_digest'] or r.get('architecture')!=arch or r.get('device_closed') is not True or r.get('routes_muted') is not True:raise ValueError('UI identity/state incomplete')
    rows=r.get('cases',[])
    if {x['id'] for x in rows}!=UI_CASES or len(rows)!=len(UI_CASES) or not all(x.get('passed') is True and x.get('evidence') for x in rows):raise ValueError('UI checklist incomplete')

def main():
    p=argparse.ArgumentParser();p.add_argument('--platform',choices=('sequoia','high-sierra'),required=True);p.add_argument('--phase',choices=tuple(REQUIRED),required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--decoder-library',type=pathlib.Path,required=True);p.add_argument('--exchange',type=pathlib.Path);a=p.parse_args()
    root=pathlib.Path.cwd();identity=json.loads((root/'identity.json').read_text());arch='arm64' if a.platform=='sequoia' else 'x86_64';version=subprocess.check_output(['/usr/bin/sw_vers','-productVersion'],universal_newlines=True).strip()
    if platform.machine()!=arch or not version.startswith('15.' if arch=='arm64' else '10.13.'):raise SystemExit('Actual platform required')
    changed=[name for name,digest in identity['artifacts_sha256'].items() if not (root/name).is_file() or sha256(root/name)!=digest]
    if changed:raise SystemExit('Frozen artifact mismatch: '+str(changed))
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);shutil.copy2(root/'identity.json',out/'identity.json');runs=[]
    env=dict(architecture=arch,os=version,os_build=subprocess.check_output(['/usr/bin/sw_vers','-buildVersion'],universal_newlines=True).strip(),cpu=subprocess.check_output(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string'],universal_newlines=True).strip(),source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],identity_sha256=sha256(root/'identity.json'),decoder_library_sha256=sha256(a.decoder_library),python=sys.version,physical_devices_opened=False,artifact_integrity_passed=True)
    (out/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
    def run(name,args,timeout=3600,expected=0):
        start=time.monotonic()
        with (out/(name+'.log')).open('w') as log:code=subprocess.run([str(x) for x in args],stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
        runs.append(dict(name=name,argv=[str(x) for x in args],exit_code=0 if code==expected else code or -1,actual_exit_code=code,expected_exit_code=expected,seconds=time.monotonic()-start,log_sha256=sha256(out/(name+'.log'))));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(name,'PASS' if code==expected else 'FAIL',flush=True)
        if code!=expected:raise ValueError('Failed required run: '+name)
    def exe(name):return ['/usr/bin/arch','-'+arch,root/'bin'/name]
    def internal(name,call):
        start=time.monotonic();result=call();(out/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n');runs.append(dict(name=name,exit_code=0,seconds=time.monotonic()-start,evidence_sha256=sha256(out/(name+'.json'))));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
    passed=False;error=''
    try:
        if a.phase=='short':
            run('m11-core',exe('LevelScopeCompositionTests')+['--root',out/'core-files','--output',out/'core.json','--samples',out/'samples']);validate_report(json.loads((out/'core.json').read_text()),identity,arch,'core')
            run('m11-video',exe('LevelScopeCompositionTests')+['--phase','video','--root',out/'video-files','--output',out/'video.json','--video',root/'video-fixtures','--prores',root/'m7-video-fixtures']);validate_report(json.loads((out/'video.json').read_text()),identity,arch,'video')
            for name,(target,flag) in UNITS.items():
                args=exe(target)+([flag,out/(name+'.json')] if flag else [])
                if name=='programmes':args+=['--m8-session',root/'m9-fixtures/m8-sessions']
                if name=='sessions':args+=['--short']
                if name=='export':args+=['--artifacts',out/'export-artifacts']
                run(name,args)
                if flag:
                    r=json.loads((out/(name+'.json')).read_text())
                    if r.get('passed') is not True or r.get('build')!=identity['measurement_build'] or r.get('architecture')!=arch:raise ValueError('Regression identity/failure: '+name)
            for name,script in (('ebu','run_ebu_measurements.py'),('itu','run_reference_measurements.py')):
                run(name,[sys.executable,'-B',root/'Tests'/script,'--binary',root/'bin/LevelScopeMeasurementTests','--vectors',root/'vectors'/name,'--report',out/(name+'.json'),'--compare',root/'baseline'/(name+'.json'),'--arch',arch])
            for name,script in (('ebu-runner','EbuReferenceRunnerTests.py'),('m11-runner','M11QualificationRunnerTests.py')):run(name,[sys.executable,'-B',root/'Tests'/script])
            session=out/'core-files/48000-0-files/composition.levelscope'
            run('legacy-rejection',['/usr/bin/arch','-'+arch,root/'legacy/M10ProgrammeTests','--only','--prepare-session',session,'--destination',out/'must-not-exist.levelscope','--output',out/'legacy-rejection.json'],expected=1)
            old=json.loads((out/'legacy-rejection.json').read_text())
            if old.get('passed') is not False or old.get('failure')!='prepare grouped qualification session' or (out/'must-not-exist.levelscope').exists():raise ValueError('M10 reader rejection not established')
            def pcm():
                rows=[dict(case=f.name,**verify_pcm(f,a.decoder_library)) for f in sorted((out/'core-files').glob('*-files'))]+[dict(case=f.name,**verify_pcm(f,a.decoder_library)) for f in sorted((out/'core-files').glob('*-groups'))]
                if len(rows)!=24:raise ValueError('Missing independent exported composition cases')
                return dict(passed=True,cases=rows)
            internal('independent-pcm',pcm)
        elif a.phase=='long':
            run('m11-long',exe('LevelScopeCompositionTests')+['--phase','long','--root',out/'long-files','--output',out/'long.json'],timeout=11100);validate_report(json.loads((out/'long.json').read_text()),identity,arch,'long');internal('independent-rf64',lambda:verify_pcm(out/'long-files',a.decoder_library,True))
        else:
            if a.exchange is None:raise ValueError('Foreign exchange directory required')
            run('m11-exchange',exe('LevelScopeCompositionTests')+['--phase','exchange','--root',out/'exchange-files','--exchange',a.exchange.resolve(),'--output',out/'exchange.json']);validate_report(json.loads((out/'exchange.json').read_text()),identity,arch,'exchange')
        validate_runs(runs,a.phase);passed=True
    except (AssertionError,ValueError,KeyError,OSError,TypeError,subprocess.TimeoutExpired) as e:error=str(e)
    (out/'summary.json').write_text(json.dumps(dict(automated_passed=passed,phase=a.phase,error=error,source_digest=identity['source_digest'],qualification_complete=False,required_runs=sorted(REQUIRED[a.phase]),remaining='Both actual platforms, short/long/exchange, cross-platform comparison and explicit UI receipts are required'),indent=2)+'\n');return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
