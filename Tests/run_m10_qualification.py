#!/usr/bin/env python3
"""Frozen M10 actual-platform checks. Physical/UI gates are explicitly separate."""
import argparse,json,pathlib,platform,subprocess,sys,time
from collect_m1_environment import sha256
from verify_m10_decoder import verify_bundle
from verify_m10_pcm import verify as verify_oracle

def validate_core(report,identity,arch):
    if not report.get('passed') or report.get('build')!=identity['measurement_build'] or report.get('architecture')!=arch or report.get('physical_devices_opened') is not False:
        raise ValueError('M10 failed or mismatched core identity')
    if not report.get('checks') or not all(row.get('passed') is True for row in report['checks']):raise ValueError('M10 failed/missing assertions')
    if len(report.get('cases',[]))!=16 or report.get('video_cases')!=10:raise ValueError('Incomplete export/reference matrix')
    expected={(r,c,k) for r in (44100,48000,96000,192000) for c in (2,6) for k in ('interleaved','mono-group')}
    if {(x['sample_rate'],len(x['per_channel_peaks']),x['input_kind']) for x in report['cases']}!=expected:raise ValueError('Incomplete rate/layout/input matrix')
    for artifact in report['cases']:
        if artifact.get('domain')!='Exported artifact' or artifact.get('verification_gain_db')!=0 or artifact.get('dither_applications')!=1 or artifact['evaluation']['verdict']!='PASS':raise ValueError('Invalid actual artifact evidence')

def validate_long(report,identity,arch):
    if not report.get('passed') or report.get('build')!=identity['measurement_build'] or report.get('architecture')!=arch or report.get('physical_devices_opened') is not False:raise ValueError('M10 long identity failure')
    r=report['long_export']
    for key,cap in [('cancel_ms',1000),('ui_p95_ms',100),('inspection_p95_ms',100),('rss_growth_bytes',256*1024*1024),('report_bytes',16*1024*1024)]:
        if type(r.get(key)) not in (int,float) or not 0<=r[key]<=cap:raise ValueError('M10 long budget: '+key)
    if r['frames']!='518400000' or r['artifact_bytes']!=9331200104 or r['annotations']!=1000 or r['ui_callbacks']<10 or r['rf64_audition_cases']!=6:raise ValueError('Missing real three-hour RF64 evidence')
    if sum(r[k] for k in ('prerequisite_ms','render_verify_ms','publish_ms','rf64_reimport_audition_ms'))>10800000:raise ValueError('M10 long time budget')
    if not report.get('checks') or not all(x['passed'] for x in report['checks']):raise ValueError('Long assertions failed')

def main():
    p=argparse.ArgumentParser();p.add_argument('--platform',required=True,choices=('sequoia','high-sierra'));p.add_argument('--output',required=True,type=pathlib.Path);p.add_argument('--phase',choices=('core','long','all'),default='all');p.add_argument('--decoder-library',required=True,type=pathlib.Path);a=p.parse_args()
    root=pathlib.Path.cwd();identity=json.loads((root/'identity.json').read_text());arch='arm64' if a.platform=='sequoia' else 'x86_64';prefix='15.' if arch=='arm64' else '10.13.'
    version=subprocess.check_output(['/usr/bin/sw_vers','-productVersion'],universal_newlines=True).strip()
    if platform.machine()!=arch or not version.startswith(prefix):raise SystemExit('Actual platform guard failed')
    changed=[n for n,d in identity['artifacts_sha256'].items() if not (root/n).is_file() or sha256(root/n)!=d]
    if changed:raise SystemExit('Candidate integrity failed: '+str(changed))
    out=a.output.resolve();out.mkdir(parents=True,exist_ok=False);runs=[];passed=True
    evidence=dict(source_digest=identity['source_digest'],identity_sha256=sha256(root/'identity.json'),build=identity['measurement_build'],architecture=arch,os=version,os_build=subprocess.check_output(['/usr/bin/sw_vers','-buildVersion'],universal_newlines=True).strip(),cpu=subprocess.check_output(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string'],universal_newlines=True).strip(),decoder_library_sha256=sha256(a.decoder_library),physical_devices_opened=False)
    (out/'environment.json').write_text(json.dumps(evidence,indent=2)+'\n')
    def run(name,argv,timeout=14400):
        start=time.monotonic()
        with (out/(name+'.log')).open('w') as log:code=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=timeout).returncode
        runs.append(dict(name=name,argv=argv,exit_code=code,seconds=time.monotonic()-start,log_sha256=sha256(out/(name+'.log'))));(out/'runs.json').write_text(json.dumps(runs,indent=2)+'\n');print(name,'PASS' if code==0 else 'FAIL',flush=True);return code==0
    exe=['/usr/bin/arch','-'+arch,str(root/'bin/LevelScopeExportTests')]
    try:
        if a.phase in ('core','all'):
            passed &= run('regressions',[sys.executable,'-B','Tests/run_m8_qualification.py','--platform',a.platform,'--output',str(out/'regressions')])
            programme=['/usr/bin/arch','-'+arch,str(root/'bin/LevelScopeProgrammeTests')]
            passed &= run('programmes',programme+['--output',str(out/'programmes.json'),'--samples',str(out/'programme-samples'),'--video',str(root/'video-fixtures'),'--prores',str(root/'m7-video-fixtures'),'--exchange-out',str(out/'programme-exchange'),'--m8-session',str(root/'m9-fixtures/m8-sessions')])
            passed &= run('programme-extended',programme+['--only','--extended-oracle','--output',str(out/'programme-extended.json')])
            passed &= run('m9-runner',[sys.executable,'-B','Tests/M9QualificationRunnerTests.py'])
            passed &= run('export',exe+['--output',str(out/'export.json'),'--artifacts',str(out/'artifacts'),'--video',str(root/'video-fixtures'),'--prores',str(root/'m7-video-fixtures')])
            report=json.loads((out/'export.json').read_text());validate_core(report,identity,arch)
            oracle=verify_oracle(out/'artifacts');(out/'independent-oracle.json').write_text(json.dumps(oracle,indent=2)+'\n')
            decoders=[dict(folder=f.name,**verify_bundle(f,a.decoder_library)) for f in sorted((out/'artifacts').glob('verified-*'))]
            if len(decoders)!=16:raise ValueError('Missing independent decoder cases')
            (out/'independent-decoders.json').write_text(json.dumps(decoders,indent=2)+'\n')
            passed &= run('m10-runner',[sys.executable,'-B','Tests/M10QualificationRunnerTests.py'])
        if a.phase in ('long','all'):
            passed &= run('long-export',exe+['--output',str(out/'long-export.json'),'--artifacts',str(out/'long-artifacts'),'--long-session',str(root/'m9-fixtures/long-48000/programme.levelscope'),'--long-destination',str(out/'long-export')])
            validate_long(json.loads((out/'long-export.json').read_text()),identity,arch)
            decoder=verify_bundle(out/'long-export',a.decoder_library,True);(out/'long-independent-decoder.json').write_text(json.dumps(decoder,indent=2)+'\n')
        error=''
    except (AssertionError,ValueError,KeyError,OSError,TypeError) as e:error=str(e);passed=False
    (out/'summary.json').write_text(json.dumps(dict(automated_passed=bool(passed),error=error,phase=a.phase,qualification_complete=False,interactive_ui='Separate evidence required',physical_video='Separate coordinated evidence required',cross_platform='Separate comparison and exchange required',source_digest=identity['source_digest']),indent=2)+'\n');return 0 if passed else 1
if __name__=='__main__':raise SystemExit(main())
