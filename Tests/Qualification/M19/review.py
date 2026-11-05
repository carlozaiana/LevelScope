"""Review explicit numerical evidence; missing native/listening evidence never qualifies.

Candidate-01 numerical checks may be reused by the UI-label-only candidate-02.
The caller supplies separately verified package integrity and that source diff.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct

p = argparse.ArgumentParser()
p.add_argument('root', type=Path)
p.add_argument('--native-sequoia', required=True)
p.add_argument('--native-high-sierra', required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
required = ['core', 'stereo', 'surround', 'rates', 'failures', 'inputs', 'gating', 'real', 'measurement', 'section', 'requirements']
result = dict(qualified=False, checks=[], cross_platform=[])
def check(ok, name):
    result['checks'].append(dict(name=name, passed=bool(ok)))
    if not ok: raise AssertionError(name)
def read(path): return json.loads(path.read_text())
def core(path):
    with path.open('rb') as f:
        h=f.read(80); data=f.read(struct.unpack('<Q',h[8:16])[0])
    check(h[:8]==b'LSSESS01' and hashlib.sha256(data).hexdigest()==h[16:80].decode(), 'Session core checksum: '+str(path))
    return json.loads(data)
def number(hexadecimal): return struct.unpack('>d',bytes.fromhex(hexadecimal))[0]
try:
    for host, architecture, os in [('sequoia','arm64','15.7.9'),('high-sierra','x86_64','10.13.6')]:
        run=read(a.root/(host+'-01/run.json'))
        check(run['machine']==architecture and os in run['os'], 'Actual platform '+host)
        check(run['required_automated_complete'] and sorted(v['group'] for v in run['cases'])==sorted(required), 'Every required direct group present: '+host)
        for row in run['cases']:
            evidence=read(a.root/(host+'-01')/(row['group']+'.json'))
            check(row['exit_code']==0 and evidence['passed'] and row['build']==evidence['build'], 'Required passed receipt '+host+'/'+row['group'])
        native=read(a.root/(a.native_sequoia if host=='sequoia' else a.native_high_sierra)/'result.json')
        check(native['passed'] and not native['physical_devices_opened'] and not native['video_loaded'], 'Native silent UI '+host)
        check(native['passages_p95_ms']<=100 and native['cancel_ms']<1000, 'Native responsiveness/cancellation '+host)
    for group in ['stereo','surround','rates']:
        x=read(a.root/'sequoia-01'/(group+'.json'))['cases'];y=read(a.root/'high-sierra-01'/(group+'.json'))['cases']
        check(len(x)==len(y), 'Matching workflow count '+group)
        for index,(left,right) in enumerate(zip(x,y)):
            check(left['passes']==right['passes'] and left['channels']==right['channels'], 'Same assessed decision/pass budget '+group+str(index))
            delta={k:abs(left[k]-right[k]) for k in ['before_lra','after_lra','gain_db']}
            check(delta['before_lra']<=.1 and delta['after_lra']<=.1 and delta['gain_db']<=.01, 'Cross-platform numerical bounds '+group+str(index))
            result['cross_platform'].append(dict(case=group+str(index), differences=delta))
    for case in ['stereo','surround','rates/44100','rates/96000','rates/192000']:
        x=core(a.root/'sequoia-01'/case/'verified.levelscope');y=core(a.root/'high-sierra-01'/case/'verified.levelscope')
        check(x['programme']['files'][0]['media']['sha256']==y['programme']['files'][0]['media']['sha256'], 'Pinned equivalent generated audio '+case)
        check(x['programme']['layout']==y['programme']['layout'] and x['jobRequirements']==y['jobRequirements'] and x['proposalPolicy']==y['proposalPolicy'], 'Same roles, requirements and policy '+case)
        check(x['current']['sectionGains']==y['current']['sectionGains'] and x['current']['gain']==y['current']['gain'], 'Exact selected envelope and global gain '+case)
    x=read(a.root/'sequoia-01/real/real-programme.json');y=read(a.root/'high-sierra-01/real/real-programme.json')
    check(x['excerpt_sha256']==y['excerpt_sha256'] and x['outcome']==y['outcome'] and x['private_passes']==y['private_passes'], 'Same bounded cinema content/outcome/pass count')
    delta={k:abs(x[k]-y[k]) for k in ['source_integrated_lufs','source_lra_lu','source_true_peak_dbtp','p10','p95']}
    check(all(v<= (.1 if 'lra' in k else .01) for k,v in delta.items()), 'Bounded cinema metric comparison')
    result['cross_platform'].append(dict(case='cinema excerpt',differences=delta))
    result['numerical_and_native_complete']=True
    result['remaining']='Review frozen-package integrity, UI-only reuse proof, supplemental resources, and separately coordinated listening observations before qualification.'
except Exception as e:
    result['numerical_and_native_complete']=False;result['failure']=str(e)
a.output.write_text(json.dumps(result,indent=2)+'\n')
raise SystemExit(0 if result['numerical_and_native_complete'] else 1)
