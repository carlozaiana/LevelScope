#!/usr/bin/env python3
"""Explicit impact-selected M13 phases. No recursive historical runners or device opening."""
import argparse,json,math,pathlib,platform,subprocess,sys,time
from collections import Counter
from collect_m1_environment import sha256
from verify_m13_pcm import verify as verify_pcm
from run_m12_qualification import MANDATORY as M12_REQUIRED_ASSERTIONS
UNITS={'measurement':('LevelScopeMeasurementTests',None),'scaffold':('LevelScopeScaffoldTests',None),
 'compliance':('LevelScopeComplianceTests','--report'),'proposals':('LevelScopeProposalTests','--report'),
 'diagnostics':('LevelScopeDiagnosticsTests','--report'),'sections':('LevelScopeSectionsTests','--report'),
 'audition':('LevelScopeAuditionTests','--report'),'multichannel':('LevelScopeMultichannelTests','--report'),
 'sessions':('LevelScopeSessionTests','--output'),'programmes':('LevelScopeProgrammeTests','--output'),'export':('LevelScopeExportTests','--output')}
REQUIRED={'short':set(UNITS)|{'m13-core','m12-core','m11-safety','ebu','itu','ebu-runner','m13-runner','m7-video-policy','legacy-rejection','independent-pcm'},
 'long':{'m13-long','independent-rf64'},'video':{'m13-video','video-timing'},'exchange':{'m13-exchange'},'cinema':{'cinema'}}
UI_CASES={'generate-noop','global-first','ready-private','explicit-apply-unverified','independent-verify','policy-failure','stale-rejection','cancel-preserves','proposal-curve-highlight','save-discard-reopen','report-export','stopped-device-closed',
 'm12-manual-draft-apply','m12-annotation-independence','m12-full-precision-entry','m12-stale-draft','m12-session-cache-recovery'}
MANDATORY={'exact twenty percent including transition support accepted','one frame beyond twenty percent rejects without epsilon',
 'representable additional attenuation never exceeds exact three dB cap','private candidate domain cannot become Current through measurement acceptance',
 'actual verification failure displayed without rollback retry or new correction','active private measurement cancellation acknowledged within one second with no draft publication',
 'already passing returns genuine no-op with no private passes','feasible global-only decision precedes section work with zero private passes',
 'capped envelope measured once and fails honestly with uncovered deficit and no second candidate',
 'already-acceptable exact local candidate skips identical final measurement without centering gain',
 'applied proposal pre-monitor samples equal independent section/global oracle at seek and EOF','actual proposal PCM24 artifact independently remeasured and passes'}
def validate_runs(rows,phase):
    names=[x['name'] for x in rows]
    if len(names)!=len(REQUIRED[phase]) or set(names)!=REQUIRED[phase] or any(x.get('exit_code')!=0 for x in rows):raise ValueError('Missing, duplicated or failed required phase checks')
def validate_report(r,identity,arch,phase):
    if r.get('passed') is not True or r.get('build')!=identity['measurement_build'] or r.get('architecture')!=arch:raise ValueError('Wrong candidate/architecture or failed report')
    checks=r.get('checks',[])
    if not checks or any(c.get('passed') is not True for c in checks):raise ValueError('Missing/failed assertions')
    if r.get('physical_devices_opened') is not False:raise ValueError('Audio-only test opened a physical device')
    if phase.startswith('m13-') or phase=='cinema':
        if r.get('video_loaded') is not (phase=='m13-video'):raise ValueError('Unexpected video execution')
    if phase=='m12-core':
        expected={(rate,ch,'M12 matrix input-'+str(mode)) for rate in (44100,48000,96000,192000) for ch in (2,6) for mode in range(4)}
        rows=r.get('cases',[])
        if len(rows)!=32 or {(x['rate'],x['channels'],x['kind']) for x in rows}!=expected or not M12_REQUIRED_ASSERTIONS.issubset({c['name'] for c in checks}):raise ValueError('Incomplete combined M12 regression')
    if phase=='m13-core':
        expected={(rate,ch,mode) for rate in (44100,48000,96000,192000) for ch in (2,6) for mode in range(4 if ch==2 else 6)}
        rows=r.get('cases',[])
        if len(rows)!=40 or {(x['rate'],x['channels'],x['mode']) for x in rows}!=expected or not MANDATORY.issubset({x['name'] for x in checks}):raise ValueError('Incomplete M13 rate/representation matrix or safeguards')
    if phase=='m13-long':
        v=r.get('long_workflow',{})
        if (v.get('frames'),v.get('rate'),v.get('channels'),v.get('intervals'),v.get('capture_cases'),v.get('private_passes'))!=('518400000',48000,6,1000,8,3):raise ValueError('Incomplete three-hour proposal workflow')
        for key,cap in (('total_ms',(240 if arch=='arm64' else 300)*60000),('cancel_ms',1000),('ui_p95_ms',100),('additional_peak_rss_bytes',256*1024*1024),('session_bytes',128*1024*1024),('report_bytes',16*1024*1024)):
            x=v.get(key)
            if type(x) not in (int,float) or not math.isfinite(x) or not 0<=x<=cap:raise ValueError('Resource budget '+key)
        if any(not 0<=v['stages_ms'][k]<=5000 for k in ('save_ms','parse_ms')):raise ValueError('Save/parse budget')
    if phase=='m13-video':
        rows=r.get('video_cases',[])
        expected=Counter((d,f,o,frame) for d in ('Source','Current State') for f in ('h264-cfr.mov','h264-vfr.mov','422-25-1.mov','lt-30000-1001.mov') for o in ('-1001/30000','7/3') for frame in ('256','240128','288256','256','480000'))
        if len(rows)!=80 or Counter((x['domain'],x['file'],x['offset'],x['audio_frame']) for x in rows)!=expected:raise ValueError('Missing eight video cases/selections')
        if any(not 0<=x['selection_ms']<10000 or x['skipped_frames']<0 or x['late_polls']<0 or x['physical_paint_test'] is not False for x in rows):raise ValueError('Invalid video timing evidence')
    if phase=='m13-exchange':
        rows=r.get('exchange_cases',[])
        if len(rows)!=6 or {x['case'] for x in rows}!={f'48000-{ch}-{mode}' for ch in (2,6) for mode in (0,1,3)}:raise ValueError('Missing actual-platform exchange cases')
        if any(x['origin_platform']==x['receiver'] or x['cache_reused'] is not False for x in rows):raise ValueError('Invalid cross-platform cache evidence')
def validate_ui(r,identity,arch):
    if r.get('passed') is not True or r.get('source_digest')!=identity['source_digest'] or r.get('architecture')!=arch or r.get('device_closed') is not True or r.get('routes_muted') is not True:raise ValueError('Wrong UI candidate/state')
    if r.get('app_sha256')!=identity['artifacts_sha256']['LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone']:raise ValueError('UI binary mismatch')
    rows=r.get('cases',[])
    if len(rows)!=len(UI_CASES) or {x['id'] for x in rows}!=UI_CASES or any(x.get('passed') is not True or not x.get('evidence') for x in rows):raise ValueError('Missing silent UI evidence')
def main():
    p=argparse.ArgumentParser();p.add_argument('--platform',choices=('sequoia','high-sierra'),required=True);p.add_argument('--phase',choices=tuple(REQUIRED),required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--decoder-library',type=pathlib.Path,required=True);p.add_argument('--long-session',type=pathlib.Path);p.add_argument('--applied-session',type=pathlib.Path);p.add_argument('--exchange',type=pathlib.Path);p.add_argument('--cinema',type=pathlib.Path);a=p.parse_args()
    root=pathlib.Path.cwd();identity=json.loads((root/'identity.json').read_text());osv=subprocess.check_output(['/usr/bin/sw_vers','-productVersion'],text=True).strip()
    arch,prefix=('arm64','15.') if a.platform=='sequoia' else ('x86_64','10.13.')
    if platform.machine()!=arch or not osv.startswith(prefix):raise ValueError('Actual requested operating system/architecture required')
    if a.phase=='cinema' and arch!='arm64':raise ValueError('Cinema evaluation is Sequoia only')
    for name,h in identity['artifacts_sha256'].items():
        if sha256(root/name)!=h:raise ValueError('Frozen artifact changed: '+name)
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);runs=[]
    env=dict(os=osv,os_build=subprocess.check_output(['/usr/bin/sw_vers','-buildVersion'],text=True).strip(),architecture=arch,cpu=subprocess.check_output(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],identity_sha256=sha256(root/'identity.json'),decoder_sha256=sha256(a.decoder_library),physical_devices_opened=False,video_loaded=a.phase=='video')
    (out/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
    def run(name,args,timeout=3600,expected=0):
        start=time.monotonic()
        with (out/(name+'.log')).open('x') as log:code=subprocess.run([str(x) for x in args],stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
        runs.append(dict(name=name,argv=[str(x) for x in args],exit_code=0 if code==expected else code or -1,actual_exit_code=code,expected_exit_code=expected,seconds=time.monotonic()-start,log_sha256=sha256(out/(name+'.log'))));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(name,'PASS' if code==expected else 'FAIL',flush=True)
        if code!=expected:raise ValueError('Required run failed: '+name)
    def exe(name):return ['/usr/bin/arch','-'+arch,root/'bin'/name]
    def internal(name,call):
        start=time.monotonic();v=call();(out/(name+'.json')).write_text(json.dumps(v,indent=2)+'\n');runs.append(dict(name=name,exit_code=0,seconds=time.monotonic()-start,evidence_sha256=sha256(out/(name+'.json'))));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
    def core(name,target,extra):
        run(name,exe(target)+extra+['--output',out/(name+'.json')],(305*60 if a.phase=='long' else 7200));validate_report(json.loads((out/(name+'.json')).read_text()),identity,arch,name)
    passed=False;error=''
    try:
        if a.phase=='short':
            core('m13-core','LevelScopeSectionProposalTests',['--root',out/'m13-files'])
            core('m12-core','LevelScopeSectionProcessingTests',['--root',out/'m12-files'])
            core('m11-safety','LevelScopeCompositionTests',['--phase','safety','--root',out/'m11-safety-files'])
            for name,(target,flag) in UNITS.items():
                args=exe(target)+([flag,out/(name+'.json')] if flag else [])
                if name=='programmes':args+=['--m8-session',root/'audio-only-m8-sessions']
                if name=='sessions':args+=['--short']
                if name=='export':args+=['--artifacts',out/'export-artifacts']
                run(name,args)
                if flag:
                    r=json.loads((out/(name+'.json')).read_text())
                    if r.get('passed') is not True or r.get('build')!=identity['measurement_build'] or r.get('architecture')!=arch:raise ValueError('Regression identity/failure: '+name)
            for name,script in (('ebu','run_ebu_measurements.py'),('itu','run_reference_measurements.py')):
                run(name,[sys.executable,'-B',root/'Tests'/script,'--binary',root/'bin/LevelScopeMeasurementTests','--vectors',root/'vectors'/name,'--report',out/(name+'.json'),'--compare',root/'baseline'/(name+'.json'),'--arch',arch])
            for name,script in (('ebu-runner','EbuReferenceRunnerTests.py'),('m13-runner','M13QualificationRunnerTests.py'),('m7-video-policy','M7VideoAcceptanceTests.py')):run(name,[sys.executable,'-B',root/'Tests'/script])
            run('legacy-rejection',['/usr/bin/arch','-'+arch,root/'legacy/M11ProgrammeTests','--only','--prepare-session',out/'m13-files/48000-2-0/applied.levelscope','--destination',out/'must-not-exist.levelscope','--output',out/'legacy-rejection.json'],expected=1)
            old=json.loads((out/'legacy-rejection.json').read_text())
            if old.get('passed') is not False or old.get('failure')!='prepare grouped qualification session' or (out/'must-not-exist.levelscope').exists():raise ValueError('Older reader rejection not established')
            def pcm():
                folders=sorted((out/'m13-files').glob('48000-*'))
                if len(folders)!=10:raise ValueError('Missing 48 kHz export matrix')
                return dict(passed=True,cases=[dict(case=f.name,**verify_pcm(f,a.decoder_library)) for f in folders])
            internal('independent-pcm',pcm)
        elif a.phase=='long':
            if not a.long_session:raise ValueError('Explicit preserved long session required')
            core('m13-long','LevelScopeSectionProposalTests',['--root',out/'long-files','--long-session',a.long_session.resolve()]);internal('independent-rf64',lambda:verify_pcm(out/'long-files',a.decoder_library,True))
        elif a.phase=='video':
            if not a.applied_session:raise ValueError('Same-candidate applied session required')
            core('m13-video','LevelScopeSectionProposalTests',['--root',out/'video-files','--video',root/'video-fixtures','--session',a.applied_session.resolve()]);run('video-timing',exe('LevelScopeVideoTests'))
        elif a.phase=='exchange':
            if not a.exchange:raise ValueError('Other actual platform fixtures required')
            core('m13-exchange','LevelScopeSectionProposalTests',['--root',out/'exchange-files','--exchange',a.exchange.resolve()])
        else:
            if not a.cinema:raise ValueError('Explicit preserved cinema evaluation copy required')
            core('cinema','LevelScopeSectionProposalTests',['--root',out/'cinema-files','--cinema',a.cinema.resolve()])
        validate_runs(runs,a.phase);passed=True
    except (AssertionError,ValueError,KeyError,OSError,TypeError,subprocess.TimeoutExpired) as e:error=str(e)
    hashes={p.name:sha256(p) for p in out.glob('*.json')}
    (out/'summary.json').write_text(json.dumps(dict(passed=passed,phase=a.phase,platform=a.platform,error=error,source_digest=identity['source_digest'],required_runs=sorted(REQUIRED[a.phase]),report_hashes=hashes,qualification_complete=False),indent=2)+'\n');return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
