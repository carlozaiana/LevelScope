"""Sequential peak decisions, exchange preparation and selected picture regressions."""
import argparse
import json
from pathlib import Path
import platform
import subprocess
import sys
import time
from run_short import sha

p=argparse.ArgumentParser();p.add_argument('candidate',type=Path);p.add_argument('short',type=Path);p.add_argument('output',type=Path);p.add_argument('--video',type=Path,required=True);a=p.parse_args()
root,short,out=a.candidate.resolve(),a.short.resolve(),a.output.resolve();out.mkdir(parents=True,exist_ok=False)
identity=json.loads((root/'identity.json').read_text());assert json.loads((short/'summary.json').read_text())['passed']
for f,h in identity['artifacts_sha256'].items():assert sha(root/f)==h
arch=platform.machine();rows=[];failure='';start=time.monotonic()
def run(name,cmd,report):
    cmd=list(map(str,cmd));t=time.monotonic();print('START',name,flush=True)
    with (out/(name+'.log')).open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
    d=json.loads(report.read_text());row=dict(name=name,command=cmd,exit_code=r.returncode,seconds=time.monotonic()-t,report_sha256=sha(report));rows.append(row);(out/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
    assert r.returncode==0 and d.get('passed',True) is True,name
    print('PASS',name,round(row['seconds'],2),flush=True);return d

def binary(name):return ['/usr/bin/arch','-'+arch,root/'bin'/name]
def signature(d):
    decisions={name:{k:v for k,v in row.items() if k not in ('main_preflight_ms','main_accept_ms')} for name,row in d['outcomes'].items()}
    return d['identities'],decisions,d['file_sha256'],d['source_lufs'],d['current_lufs'],d['base_build']
try:
    peaks=[]
    for ch in (2,6):
        expected=None
        for mode in ('serial','bounded','app-mode'):
            name='peaks-'+str(ch)+'-'+mode
            d=run(name,binary('peaks/'+mode)+[root/'peak-fixtures'/('programme-'+str(ch)+'.wav'),out/(name+'.json'),'--decisions'],out/(name+'.json'))
            assert d['base_build']==identity['measurement_build'] and len(d['outcomes'])==7
            sig=signature(d)
            if expected is None:expected=sig
            assert sig==expected,'Exact serial proposal/measurement mismatch: '+name
            peaks.append(dict(layout=ch,mode=mode,operations=d['operations']))
    (out/'peak-policy.json').write_text(json.dumps(dict(passed=True,architecture=arch,exact_match=True,cases=peaks),indent=2)+'\n')
    run('exchange-prepare',binary('session-exchange')+['prepare',short,out/'exchange-source'],out/'exchange-source/result.json')
    run('video-matrix',binary('LevelScopeSectionProposalTests')+['--video',a.video,'--session',short/'section-proposal-files/48000-6-0/applied.levelscope',
        '--root',out/'video-matrix-files','--output',out/'video-matrix.json'],out/'video-matrix.json')
    run('video-native',binary('video-pressure')+[short/'native-surround/job.levelscope',a.video,out/'video-native'],out/'video-native/result.json')
    assert len(rows)==9 and {r['name'] for r in rows}=={'peaks-'+str(c)+'-'+m for c in (2,6) for m in ('serial','bounded','app-mode')}|{'exchange-prepare','video-matrix','video-native'}
except Exception as e:failure=str(e)
summary=dict(passed=not failure,failure=failure,architecture=arch,build=identity['measurement_build'],identity_sha256=sha(root/'identity.json'),seconds=time.monotonic()-start,commands=rows,physical_devices_opened=False)
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='commands'},indent=2));raise SystemExit(0 if not failure else 1)
