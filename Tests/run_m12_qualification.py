#!/usr/bin/env python3
"""Explicit Sequoia audio-only phase. Never access Golem or launch video tests."""
import argparse, json, math, pathlib, platform, shutil, subprocess, sys, time
from collect_m1_environment import sha256
from verify_m12_pcm import verify as verify_pcm

UNITS={'measurement':('LevelScopeMeasurementTests',None),'scaffold':('LevelScopeScaffoldTests',None),
 'compliance':('LevelScopeComplianceTests','--report'),'proposals':('LevelScopeProposalTests','--report'),
 'diagnostics':('LevelScopeDiagnosticsTests','--report'),'sections':('LevelScopeSectionsTests','--report'),
 'audition':('LevelScopeAuditionTests','--report'),'multichannel':('LevelScopeMultichannelTests','--report'),
 'sessions':('LevelScopeSessionTests','--output'),'programmes':('LevelScopeProgrammeTests','--output'),
 'export':('LevelScopeExportTests','--output')}
REQUIRED={'short':set(UNITS)|{'m12-core','m11-safety','ebu','itu','ebu-runner','m12-runner','legacy-rejection','independent-pcm'},
          'long':{'m12-long','independent-rf64'}}
UI_CASES={'manual-sections','draft-applied-curves','explicit-apply-unverified','remeasure-actual',
 'annotation-independence','stale-draft','already-pass-noop','session-roundtrip','draft-save-discard',
 'export-report','device-closed-muted'}
MANDATORY={'all-zero section draft preserves legacy identity and measurement',
 'geometry edit stales draft without mutating audio','cancelled verification retains applied gains and unverified state',
 'identical applied section plan is an exact no-op','damaged Current cache discarded with exact applied plan and Source cache recovered',
 'missing group member retains processing settings and blocks partial audio/current pass',
 'finite over-full-scale samples measured without clipping; actual failure displayed and gain retained',
 'one representable gain step changes measurement identity before display rounding',
 'actual Current worker cancellation after section samples publishes no result',
 'editable gain text round-trips full double precision including scientific notation',
 'already-pass global trim retains exact plan and measurement'}

def reuse_compatible(old,new,target):
    production=lambda i:{k:v for k,v in i['source_files_sha256'].items() if k=='CMakeLists.txt' or k.startswith('Source/')}
    key='bin/'+target
    return bool(production(old)) and production(old)==production(new) and old['measurement_build']==new['measurement_build'] and key in old['artifacts_sha256'] and old['artifacts_sha256'][key]==new['artifacts_sha256'].get(key)

def validate_report(r,identity,phase):
    if r.get('passed') is not True or r.get('build')!=identity['measurement_build'] or r.get('architecture')!='arm64' or r.get('phase')!=phase or r.get('physical_devices_opened') is not False or r.get('video_loaded') is not False:raise ValueError('Failed/wrong M12 candidate or phase')
    checks=r.get('checks',[])
    if not checks or any(x.get('passed') is not True for x in checks):raise ValueError('Missing/failed assertions')
    if phase=='sequoia-audio':
        expected={(rate,ch,'M12 matrix input-'+str(mode)) for rate in (44100,48000,96000,192000) for ch in (2,6) for mode in range(4)}
        actual={(x['rate'],x['channels'],x['kind']) for x in r.get('cases',[])}
        if expected!=actual or len(r['cases'])!=32 or not MANDATORY.issubset({x['name'] for x in checks}):raise ValueError('Incomplete matrix or required negative tests')
    elif phase=='sequoia-long':
        v=r.get('long_workflow',{})
        if (v.get('frames'),v.get('rate'),v.get('channels'),v.get('intervals'),v.get('capture_cases'))!=('518400000',48000,6,1000,8):raise ValueError('Incomplete real long workflow')
        for key,cap in (('total_ms',10800000),('cancel_ms',1000),('ui_p95_ms',100),('additional_peak_rss_bytes',256*1024*1024),('session_bytes',128*1024*1024),('report_bytes',16*1024*1024)):
            x=v.get(key)
            if type(x) not in (int,float) or not math.isfinite(x) or not 0<=x<=cap:raise ValueError('Long resource budget: '+key)
        if any(not 0<=v['stages_ms'][key]<=5000 for key in ('save_ms','parse_ms')):raise ValueError('Persistence resource budget')
    else:raise ValueError('Unknown phase')

def validate_runs(runs,phase):
    names=[x['name'] for x in runs]
    if set(names)!=REQUIRED[phase] or len(names)!=len(REQUIRED[phase]) or any(x.get('exit_code')!=0 for x in runs):raise ValueError('Missing/duplicate/failed required Sequoia run')

def validate_ui(r,identity):
    if r.get('passed') is not True or r.get('source_digest')!=identity['source_digest'] or r.get('architecture')!='arm64' or r.get('device_closed') is not True or r.get('routes_muted') is not True or r.get('video_loaded') is not False:raise ValueError('Wrong/missing silent UI identity/state')
    if r.get('app_sha256')!=identity['artifacts_sha256']['LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone']:raise ValueError('Wrong/missing silent UI executable identity')
    rows=r.get('cases',[])
    if {x['id'] for x in rows}!=UI_CASES or len(rows)!=len(UI_CASES) or any(x.get('passed') is not True or not x.get('evidence') for x in rows):raise ValueError('Incomplete UI checklist')

def main():
    p=argparse.ArgumentParser();p.add_argument('--phase',choices=tuple(REQUIRED),required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--decoder-library',type=pathlib.Path,required=True);p.add_argument('--long-session',type=pathlib.Path);p.add_argument('--reuse-short',type=pathlib.Path);p.add_argument('--reuse-package',type=pathlib.Path);a=p.parse_args()
    root=pathlib.Path.cwd();identity=json.loads((root/'identity.json').read_text());os_version=subprocess.check_output(['/usr/bin/sw_vers','-productVersion'],text=True).strip()
    if platform.machine()!='arm64' or not os_version.startswith('15.'):raise SystemExit('Actual arm64 Sequoia required')
    for name,digest in identity['artifacts_sha256'].items():
        if not (root/name).is_file() or sha256(root/name)!=digest:raise ValueError('Frozen artifact mismatch: '+name)
    for name,asset in identity['referenced_assets'].items():
        if sha256(pathlib.Path(asset['path']))!=asset['sha256']:raise ValueError('Historical input asset mismatch: '+name)
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);runs=[]
    (out/'environment.json').write_text(json.dumps(dict(os=os_version,os_build=subprocess.check_output(['/usr/bin/sw_vers','-buildVersion'],text=True).strip(),architecture=platform.machine(),cpu=subprocess.check_output(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],identity_sha256=sha256(root/'identity.json'),decoder_sha256=sha256(a.decoder_library),physical_devices_opened=False,video_loaded=False,high_sierra_tested=False),indent=2)+'\n')
    reused={}
    if a.reuse_short or a.reuse_package:
        if a.phase!='short' or not a.reuse_short or not a.reuse_package:raise ValueError('Short reuse requires an explicit result root and original package')
        old=json.loads((a.reuse_package/'identity.json').read_text());environment=json.loads((a.reuse_short/'environment.json').read_text());now=json.loads((out/'environment.json').read_text())
        if environment.get('identity_sha256')!=sha256(a.reuse_package/'identity.json') or any(environment.get(k)!=now[k] for k in ('os','os_build','architecture','cpu','decoder_sha256','physical_devices_opened','video_loaded')):raise ValueError('Reuse environment mismatch')
        binaries={k:v[0] for k,v in UNITS.items()};binaries.update({'m12-core':'LevelScopeSectionProcessingTests','m11-safety':'LevelScopeCompositionTests'})
        for row in json.loads((a.reuse_short/'runs.json').read_text()):
            name=row['name']
            if name in binaries and row.get('exit_code')==0 and reuse_compatible(old,identity,binaries[name]):
                if sha256(a.reuse_package/'bin'/binaries[name])!=old['artifacts_sha256']['bin/'+binaries[name]] or sha256(a.reuse_short/(name+'.log'))!=row['log_sha256']:raise ValueError('Changed reusable binary/log')
                reused[name]=dict(row,reused=True,reused_from=str(a.reuse_short.resolve()),original_source_digest=row.get('original_source_digest',old['source_digest']),original_identity_sha256=row.get('original_identity_sha256',environment['identity_sha256']),reuse_chain=row.get('reuse_chain',[])+[dict(result_root=str(a.reuse_short.resolve()),source_digest=old['source_digest'],identity_sha256=environment['identity_sha256'])],reuse_reason='Exact production source/build and executable SHA-256; identical actual Sequoia environment. Original results retained.')
    def run(name,args,timeout=3600,expected=0):
        if name in reused:
            shutil.copy2(a.reuse_short/(name+'.log'),out/(name+'.log'))
            report={'m12-core':'core','m11-safety':'composition-safety'}.get(name,name)+'.json'
            if (a.reuse_short/report).exists():shutil.copy2(a.reuse_short/report,out/report)
            if name=='m12-core':(out/'core-files').symlink_to((a.reuse_short/'core-files').resolve(),target_is_directory=True)
            runs.append(reused[name]);(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(name,'REUSED exact executable',flush=True);return
        start=time.monotonic()
        with (out/(name+'.log')).open('w') as log:code=subprocess.run([str(x) for x in args],stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
        runs.append(dict(name=name,argv=[str(x) for x in args],exit_code=0 if code==expected else code or -1,actual_exit_code=code,expected_exit_code=expected,seconds=time.monotonic()-start,log_sha256=sha256(out/(name+'.log'))));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(name,'PASS' if code==expected else 'FAIL',flush=True)
        if code!=expected:raise ValueError('Failed required run: '+name)
    def exe(name):return ['/usr/bin/arch','-arm64',root/'bin'/name]
    def asset(name):return pathlib.Path(identity['referenced_assets'][name]['path'])
    def internal(name,call):
        start=time.monotonic();result=call();(out/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n');runs.append(dict(name=name,exit_code=0,seconds=time.monotonic()-start,evidence_sha256=sha256(out/(name+'.json'))));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
    passed=False;error=''
    try:
        if a.phase=='short':
            run('m12-core',exe('LevelScopeSectionProcessingTests')+['--root',out/'core-files','--output',out/'core.json']);validate_report(json.loads((out/'core.json').read_text()),identity,'sequoia-audio')
            run('m11-safety',exe('LevelScopeCompositionTests')+['--phase','safety','--root',out/'composition-safety','--output',out/'composition-safety.json'])
            safety=json.loads((out/'composition-safety.json').read_text())
            if safety.get('passed') is not True or safety.get('build')!=identity['measurement_build'] or safety.get('architecture')!='arm64' or safety.get('phase')!='safety' or not safety.get('checks') or safety.get('video_cases') or safety.get('physical_devices_opened') is not False:raise ValueError('Missing/wrong-candidate audio-only composition safety evidence')
            for name,(target,flag) in UNITS.items():
                args=exe(target)+([flag,out/(name+'.json')] if flag else [])
                if name=='programmes':args+=['--m8-session',root/'audio-only-m8-sessions']
                if name=='sessions':args+=['--short']
                if name=='export':args+=['--artifacts',out/'export-artifacts']
                run(name,args)
                if flag:
                    r=json.loads((out/(name+'.json')).read_text())
                    if r.get('passed') is not True or r.get('build')!=identity['measurement_build'] or r.get('architecture')!='arm64':raise ValueError('Regression identity/failure: '+name)
            for name,script,example in (('ebu','run_ebu_measurements.py','readme.txt'),('itu','run_reference_measurements.py',None)):
                paths=[v['path'] for k,v in identity['referenced_assets'].items() if k.startswith('vectors/'+name+'/')]
                run(name,[sys.executable,'-B',root/'Tests'/script,'--binary',root/'bin/LevelScopeMeasurementTests','--vectors',pathlib.Path(paths[0]).parent,'--report',out/(name+'.json'),'--compare',asset('baseline/'+name+'.json'),'--arch','arm64'])
            for name,script in (('ebu-runner','EbuReferenceRunnerTests.py'),('m12-runner','M12QualificationRunnerTests.py')):run(name,[sys.executable,'-B',root/'Tests'/script])
            run('legacy-rejection',['/usr/bin/arch','-arm64',root/'legacy/M11ProgrammeTests','--only','--prepare-session',out/'core-files/48000-2-0/section.levelscope','--destination',out/'must-not-exist.levelscope','--output',out/'legacy-rejection.json'],expected=1)
            old=json.loads((out/'legacy-rejection.json').read_text())
            if old.get('passed') is not False or old.get('failure')!='prepare grouped qualification session' or (out/'must-not-exist.levelscope').exists():raise ValueError('Old reader rejection not established')
            def pcm():
                rows=[dict(case=f.name,**verify_pcm(f,a.decoder_library)) for f in sorted((out/'core-files').glob('48000-*-*'))]
                if len(rows)!=8:raise ValueError('Missing independent export cases')
                return dict(passed=True,cases=rows)
            internal('independent-pcm',pcm)
        else:
            if a.long_session is None:raise ValueError('Explicit immutable M11 input session required')
            run('m12-long',exe('LevelScopeSectionProcessingTests')+['--long-session',a.long_session.resolve(),'--root',out/'long-files','--output',out/'long.json'],11100)
            validate_report(json.loads((out/'long.json').read_text()),identity,'sequoia-long');internal('independent-rf64',lambda:verify_pcm(out/'long-files',a.decoder_library,True))
        validate_runs(runs,a.phase);passed=True
    except (AssertionError,ValueError,KeyError,OSError,TypeError,subprocess.TimeoutExpired) as e:error=str(e)
    (out/'summary.json').write_text(json.dumps(dict(sequoia_automated_passed=passed,phase=a.phase,error=error,source_digest=identity['source_digest'],qualification_complete=False,required_runs=sorted(REQUIRED[a.phase]),deferred=['Actual High Sierra after M13','All video loading/decoding/playback checks after M13'],remaining='Required Sequoia short, long and silent UI receipts; no publication or M13 authorization'),indent=2)+'\n');return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
