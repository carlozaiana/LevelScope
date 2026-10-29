#!/usr/bin/env python3
"""Explicit Sequoia audio-only M14 phases. Deferred coverage cannot become a pass."""
import argparse,json,math,pathlib,platform,subprocess,sys,time
from collect_m1_environment import sha256
from verify_m14_pcm import verify as verify_pcm,raw as verify_raw
UNITS={'measurement':('LevelScopeMeasurementTests',None),'scaffold':('LevelScopeScaffoldTests',None),
 'compliance':('LevelScopeComplianceTests','--report'),'proposals':('LevelScopeProposalTests','--report'),
 'diagnostics':('LevelScopeDiagnosticsTests','--report'),'sections':('LevelScopeSectionsTests','--report'),
 'audition':('LevelScopeAuditionTests','--report'),'multichannel':('LevelScopeMultichannelTests','--report'),
 'sessions':('LevelScopeSessionTests','--output'),'programmes':('LevelScopeProgrammeTests','--output'),
 'export':('LevelScopeExportTests','--output')}
REQUIRED={'short':set(UNITS)|{'m14-core','m12-regression','m13-regression','composition-smoke','composition-safety','ebu','itu','ebu-runner','m14-runner','legacy-rejection','independent-pcm'},
          'long':{'m14-long','independent-long-pcm'},'cinema':{'m14-cinema'},'legacy':{'legacy-baseline','legacy-rejection'}}
UI_CASES={'mono-import-roles','source-measure-profile','manual-draft-apply','independent-current','proposal-apply-verify','already-pass','session-reader-warning-reopen','export-report','rf64-compatible-input','invalid-input-preserves','stopped-muted-device-closed'}
DEFERRED=['actual-high-sierra-runtime-ui','cross-mac-exchange','video-execution','three-hour-mono-rf64-end-to-end']
MANDATORY={'existing bounded method finds measured mono candidate within pass cap','mono Apply preserves Source and makes Current unverified',
 'independent Current mono measurement and sample oracle pass','already passing mono makes no change or private scan',
 'export remains mono with exact confirmed role and explicit mapping','odd PCM24 padding is not audio and artifact has zero extra processing',
 'input compatibility never relaxes exported artifact validation','every short compatible RF64 decoded sample matches independent PCM bytes',
 'bounded seek beyond four GiB decodes expected last PCM sample','all mono WAV AIFF FLAC manual samples equal independent oracle',
 'one-channel loudness differs from explicit duplicated stereo by calculated energy sum','mono role change invalidates Source and Current'}
def validate_runs(rows,phase):
    names=[r['name'] for r in rows]
    if len(names)!=len(REQUIRED[phase]) or set(names)!=REQUIRED[phase] or any(r.get('exit_code')!=0 for r in rows):raise ValueError('Missing duplicated or failed required M14 checks')
def validate_legacy(applied,baseline):
    prefix='M8 core migration remains interleaved: '
    # This old helper requires schema 1 after reading. A valid schema-4 baseline
    # fails that helper assertion with the successful core-read message. The applied mono
    # plan must fail the parser itself, with the specific capability error.
    for r in (applied,baseline):
        if r.get('passed') is not False or r.get('architecture')!='arm64' or r.get('build')!='5fd55bb1a81c80b4b556715936b6873bccae29a5d04c00b0890a9b43296a20d2':raise ValueError('Wrong legacy reader identity')
    if applied.get('failure')!=prefix+'Applied section plan does not match this programme' or baseline.get('failure')!=prefix+'Session core read; media verification required':raise ValueError('Explicit old-reader rejection/control not established')
def validate_report(r,identity,name):
    if r.get('passed') is not True or r.get('build')!=identity['measurement_build'] or r.get('architecture')!='arm64':raise ValueError('Wrong/failed Sequoia candidate report')
    checks=r.get('checks',[])
    if not checks or any(x.get('passed') is not True for x in checks):raise ValueError('Missing or failed assertions')
    if r.get('physical_devices_opened') is not False or r.get('video_loaded') is not False:raise ValueError('Deferred media/device execution')
    if name=='m14-core':
        rows=r.get('cases',[])
        if len(rows)!=12 or {(x['rate'],x['role']) for x in rows}!={(rate,role) for rate in (44100,48000,96000,192000) for role in (1,2,3)} or not MANDATORY.issubset({x['name'] for x in checks}):raise ValueError('Incomplete mono matrix/safeguards')
        if any(x['frames']!=x['rate']*10+1 or x['mask']!={1:1,2:2,3:4}[x['role']] or not 0<x['private_passes']<=3 for x in rows):raise ValueError('Invalid mono matrix identity')
    if name=='m14-long':
        v=r.get('long_workflow',{})
        if (v.get('frames'),v.get('rate'),v.get('channels'),v.get('intervals'),v.get('capture_cases'))!=('1036800000',192000,1,1000,8):raise ValueError('Incomplete 90-minute workload')
        if v.get('three_hour_qualified') is not False or v.get('rf64_end_to_end_qualified') is not False or not 1<=v.get('private_passes',0)<=3:raise ValueError('Unsupported duration/pass claim')
        for k,cap in (('total_ms',120*60000),('cancel_ms',1000),('ui_p95_ms',100),('additional_peak_rss_bytes',256*1024*1024),('session_bytes',128*1024*1024),('report_bytes',16*1024*1024)):
            x=v.get(k)
            if type(x) not in (int,float) or not math.isfinite(x) or not 0<=x<=cap:raise ValueError('Resource budget '+k)
        if any(not 0<=v['stages_ms'][k]<=5000 for k in ('save_ms','parse_ms')):raise ValueError('Save/parse budget')
    if name=='m14-cinema':
        c=r.get('cinema',{})
        if (c.get('frames'),c.get('rate'),c.get('channels'))!=('269875200',48000,6) or c.get('original_sha256')!='1fb7f93a33d600bf0a62d63e642e5703f03aec8e6920f8a6e510cf65f61fdeaf':raise ValueError('Wrong cinema identity')
        if len(c.get('bounded_reads',[]))!=6 or c.get('full_proposal_export_repeated') is not False:raise ValueError('Incomplete bounded cinema comparison')
def validate_ui(r,identity):
    if r.get('passed') is not True or r.get('source_digest')!=identity['source_digest'] or r.get('architecture')!='arm64' or r.get('device_closed') is not True or r.get('routes_muted') is not True or r.get('video_loaded') is not False:raise ValueError('Wrong UI candidate/state')
    if r.get('app_sha256')!=identity['artifacts_sha256']['LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone']:raise ValueError('UI binary mismatch')
    rows=r.get('cases',[])
    if len(rows)!=len(UI_CASES) or {x['id'] for x in rows}!=UI_CASES or any(x.get('passed') is not True or not x.get('evidence') for x in rows):raise ValueError('Missing silent UI evidence')
def main():
    p=argparse.ArgumentParser();p.add_argument('--phase',choices=tuple(REQUIRED),required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--decoder-library',type=pathlib.Path,required=True);p.add_argument('--cinema',type=pathlib.Path);p.add_argument('--copy',type=pathlib.Path);p.add_argument('--mono-fixtures',type=pathlib.Path);a=p.parse_args()
    root=pathlib.Path.cwd();identity=json.loads((root/'identity.json').read_text());osv=subprocess.check_output(['/usr/bin/sw_vers','-productVersion'],text=True).strip()
    if platform.machine()!='arm64' or not osv.startswith('15.'):raise ValueError('Actual arm64 Sequoia phase required; High Sierra deferred')
    for name,h in identity['artifacts_sha256'].items():
        if sha256(root/name)!=h:raise ValueError('Frozen artifact changed: '+name)
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);runs=[]
    env=dict(os=osv,os_build=subprocess.check_output(['/usr/bin/sw_vers','-buildVersion'],text=True).strip(),architecture='arm64',cpu=subprocess.check_output(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],identity_sha256=sha256(root/'identity.json'),decoder_sha256=sha256(a.decoder_library),physical_devices_opened=False,video_loaded=False)
    (out/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
    def run(name,args,timeout=3600,expected=0):
        start=time.monotonic()
        with (out/(name+'.log')).open('x') as log:code=subprocess.run([str(x) for x in args],stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
        runs.append(dict(name=name,argv=[str(x) for x in args],exit_code=0 if code==expected else code or -1,actual_exit_code=code,expected_exit_code=expected,seconds=time.monotonic()-start,log_sha256=sha256(out/(name+'.log'))));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(name,'PASS' if code==expected else 'FAIL',flush=True)
        if code!=expected:raise ValueError('Required run failed: '+name)
    def exe(name):return ['/usr/bin/arch','-arm64',root/'bin'/name]
    def internal(name,call):
        start=time.monotonic();value=call();(out/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n');runs.append(dict(name=name,exit_code=0,seconds=time.monotonic()-start,evidence_sha256=sha256(out/(name+'.json'))));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(name,'PASS',flush=True)
    def core(name,target,extra):
        run(name,exe(target)+extra+['--output',out/(name+'.json')],timeout=125*60 if a.phase=='long' else 3600)
        r=json.loads((out/(name+'.json')).read_text())
        if name.startswith('m14-'):validate_report(r,identity,name)
        elif r.get('passed') is not True or r.get('build')!=identity['measurement_build'] or r.get('architecture')!='arm64':raise ValueError('Failed regression identity: '+name)
    passed=False;error=''
    try:
        if a.phase=='short':
            core('m14-core','LevelScopeMonoTests',['--root',out/'mono-files'])
            core('m12-regression','LevelScopeSectionProcessingTests',['--smoke','--root',out/'m12-files'])
            core('m13-regression','LevelScopeSectionProposalTests',['--smoke','--root',out/'m13-files'])
            core('composition-smoke','LevelScopeCompositionTests',['--smoke','--root',out/'composition-files'])
            core('composition-safety','LevelScopeCompositionTests',['--phase','safety','--root',out/'composition-safety-files'])
            for name,(target,flag) in UNITS.items():
                args=exe(target)+([flag,out/(name+'.json')] if flag else [])
                if name=='sessions':args+=['--short']
                if name=='programmes':args+=['--m8-session',root/'audio-only-m8-sessions']
                if name=='export':args+=['--artifacts',out/'export-artifacts']
                run(name,args)
                if flag:
                    r=json.loads((out/(name+'.json')).read_text())
                    if r.get('passed') is not True or r.get('build')!=identity['measurement_build'] or r.get('architecture')!='arm64':raise ValueError('Failed regression identity: '+name)
            for name,script in (('ebu','run_ebu_measurements.py'),('itu','run_reference_measurements.py')):
                run(name,[sys.executable,'-B',root/'Tests'/script,'--binary',root/'bin/LevelScopeMeasurementTests','--vectors',root/'vectors'/name,'--report',out/(name+'.json'),'--compare',root/'baseline'/(name+'.json'),'--arch','arm64'])
            for name,script in (('ebu-runner','EbuReferenceRunnerTests.py'),('m14-runner','M14QualificationRunnerTests.py')):run(name,[sys.executable,'-B',root/'Tests'/script])
            run('legacy-rejection',['/usr/bin/arch','-arm64',root/'legacy/M13ProgrammeTests','--only','--m8-session',out/'mono-files/48000-3/applied.levelscope','--output',out/'legacy-rejection.json'],expected=1)
            old=json.loads((out/'legacy-rejection.json').read_text())
            if old.get('passed') is not False or old.get('failure')!='M8 core migration remains interleaved: Applied section plan does not match this programme':raise ValueError('Old reader did not safely reject mono plan')
            def independent():
                folders=sorted((out/'mono-files').glob('[0-9]*-*'))
                if len(folders)!=12:raise ValueError('Missing mono export matrix')
                return dict(passed=True,raw=verify_raw(out/'mono-files/raw-pcm'),cases=[dict(case=f.name,**verify_pcm(f,a.decoder_library)) for f in folders])
            internal('independent-pcm',independent)
        elif a.phase=='long':
            core('m14-long','LevelScopeMonoTests',['--long','--root',out/'long-files']);internal('independent-long-pcm',lambda:verify_pcm(out/'long-files',a.decoder_library,True))
        elif a.phase=='cinema':
            if not a.cinema or not a.copy:raise ValueError('Explicit preserved original and M13 copy required')
            core('m14-cinema','LevelScopeMonoTests',['--cinema',a.cinema.resolve(),'--copy',a.copy.resolve(),'--root',out/'cinema-files'])
        else:
            if not a.mono_fixtures:raise ValueError('Explicit frozen mono fixtures required')
            for name,file in (('legacy-baseline','baseline.levelscope'),('legacy-rejection','applied.levelscope')):
                run(name,['/usr/bin/arch','-arm64',root/'legacy/M13ProgrammeTests','--only','--m8-session',a.mono_fixtures.resolve()/file,'--output',out/(name+'.json')],expected=1)
            validate_legacy(json.loads((out/'legacy-rejection.json').read_text()),json.loads((out/'legacy-baseline.json').read_text()))
        validate_runs(runs,a.phase);passed=True
    except (AssertionError,ValueError,KeyError,OSError,TypeError,subprocess.TimeoutExpired) as e:error=str(e)
    hashes={f.name:sha256(f) for f in out.glob('*.json')}
    (out/'summary.json').write_text(json.dumps(dict(passed=passed,phase=a.phase,platform='sequoia',error=error,source_digest=identity['source_digest'],required_runs=sorted(REQUIRED[a.phase]),report_hashes=hashes,qualification_complete=False,deferred=DEFERRED),indent=2)+'\n');return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
