"""Explicit sequential M14–M18 qualification gate on either actual Mac."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

REQUIRED={'measurement','compliance','global-proposal','diagnostics','sections','audition',
          'multichannel','sessions','session-playback','export','section-processing',
          'section-proposal','preferences','proposal-workflow','mono','markers',
          'requirements','waveform-core','waveform-service','timeline','peak-parity',
          'video-core','ebu','itu','native-stereo','native-surround','audit-stereo','audit-surround'}

def main():
    p=argparse.ArgumentParser();p.add_argument('candidate',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--ebu',type=Path,required=True);p.add_argument('--itu',type=Path,required=True)
    a=p.parse_args();root=a.candidate.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    identity=json.loads((root/'identity.json').read_text());arch=platform.machine()
    assert arch in ('arm64','x86_64') and platform.system()=='Darwin'
    for path,digest in identity['artifacts_sha256'].items():assert sha(root/path)==digest,path
    rows=[];start=time.monotonic();failure=''
    def run(name,command,report=None,text_pass=None):
        command=list(map(str,command));t=time.monotonic();print('START',name,flush=True)
        with (out/(name+'.log')).open('w') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        row=dict(name=name,command=command,exit_code=r.returncode,seconds=time.monotonic()-t)
        if report:
            data=json.loads(Path(report).read_text())
            row.update(report_sha256=sha(report),passed=data.get('passed'),build=data.get('build',data.get('measurement_build')))
            # EBU has a documented official-check result rather than a top-level passed.
            if name=='ebu':row['passed']=data.get('passed',data.get('minimum_requirements_passed'))
            if name=='itu':row['passed']=data.get('passed',all(x['passed'] for x in data['results']))
        elif text_pass:row['passed']=text_pass in (out/(name+'.log')).read_text()
        rows.append(row);(out/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
        assert r.returncode==0 and row.get('passed') is True,'Failed or missing result: '+name
        if row.get('build'):assert row['build']==identity['measurement_build'],'Mismatched build: '+name
        print('PASS',name,round(row['seconds'],2),flush=True)
    def binary(name):return ['/usr/bin/arch','-'+arch,root/'bin'/name]
    try:
        for name,target in [('measurement','Measurement'),('compliance','Compliance'),('global-proposal','Proposal'),
                            ('diagnostics','Diagnostics'),('sections','Sections'),('audition','Audition'),('multichannel','Multichannel')]:
            run(name,binary('LevelScope'+target+'Tests')+['--report',out/(name+'.json')],out/(name+'.json'))
        run('sessions',binary('LevelScopeSessionTests')+['--short','--output',out/'sessions.json'],out/'sessions.json')
        run('session-playback',binary('LevelScopeSessionPlaybackTests')+['--report',out/'session-playback.json'],out/'session-playback.json')
        run('export',binary('LevelScopeExportTests')+['--artifacts',out/'export-files','--output',out/'export.json'],out/'export.json')
        for name,target,flags in [('section-processing','SectionProcessing',['--smoke']),('section-proposal','SectionProposal',['--smoke']),
                                  ('preferences','SectionProposal',['--preferences']),('proposal-workflow','SectionProposal',['--m18-workflow']),
                                  ('mono','Mono',['--smoke']),('requirements','JobRequirements',[]),
                                  ('waveform-service','Waveform',[]),('timeline','TimelineProgress',[])]:
            run(name,binary('LevelScope'+target+'Tests')+flags+['--root',out/(name+'-files'),'--output',out/(name+'.json')],out/(name+'.json'))
        run('markers',binary('LevelScopeMarkerTests')+[out/'marker-files',root/'source/Tests/Markers/nuendo-15.0.30-points.csv'],out/'marker-files/marker-tests.json')
        run('waveform-core',binary('LevelScopeWaveformCoreTests'),text_pass='PASS ')
        run('peak-parity',binary('LevelScopePeakChannelParity'),text_pass='PASS')
        run('video-core',binary('LevelScopeVideoTests'),text_pass='PASS')
        for name,script,vectors in [('ebu','run_ebu_measurements.py',a.ebu),('itu','run_reference_measurements.py',a.itu)]:
            run(name,[sys.executable,root/'source/Tests'/script,'--binary',root/'bin/LevelScopeMeasurementTests',
                      '--vectors',vectors,'--report',out/(name+'.json'),'--arch',arch,'--workers','1'],out/(name+'.json'))
        for channels,name in [(2,'stereo'),(6,'surround')]:
            run('native-'+name,binary('short-acceptance')+[root/'fixtures'/('programme-'+str(channels)+'.wav'),out/('native-'+name),str(channels)],out/('native-'+name)/'result.json')
            run('audit-'+name,[sys.executable,root/'source/Tests/Profiling/ShortAcceptance/check_export.py',out/('native-'+name),
                              root/'fixtures'/('programme-'+str(channels)+'.wav')],out/('audit-'+name+'.log'))
        assert len(rows)==len(REQUIRED) and {r['name'] for r in rows}==REQUIRED,'Missing or unexpected required group'
    except Exception as e:failure=str(e)
    report=dict(passed=not failure,failure=failure,build=identity['measurement_build'],identity_sha256=sha(root/'identity.json'),
                architecture=arch,os=subprocess.check_output(['sw_vers','-productVersion'],text=True).strip(),
                seconds=time.monotonic()-start,commands=rows,physical_devices_opened=False,full_qualification_complete=False)
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='commands'},indent=2))
    return 1 if failure else 0
if __name__=='__main__':raise SystemExit(main())
