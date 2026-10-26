#!/usr/bin/env python3
"""Bounded M7/M9 core-capture comparisons; does not claim their long qualification."""
import hashlib,json,math,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'Tests'))
from compare_m7_reports import compare_capture
r=Path(__file__).parent
roots=[r/'sequoia-candidate-03',r/'high-sierra-candidate-03']
env=[json.loads((p/'environment.json').read_text()) for p in roots]
assert {e['architecture'] for e in env}=={'arm64','x86_64'}
assert env[0]['source_digest']==env[1]['source_digest']
assert env[0]['os'].startswith('15.') and env[1]['os'].startswith('10.13.')
result=dict(passed=True,source_digest=env[0]['source_digest'],scope='M7/M9 core captured samples only; long measurements are separate evidence',suites=[])
for suite,relative,key,folder,count in [('M7','regressions/regressions/multichannel.json','sample_cases','regressions/regressions/m7-samples',218),('M9','programmes.json','cases','programme-samples',56)]:
    reports=[json.loads((p/relative).read_text()) for p in roots]
    assert all(x['passed'] for x in reports) and reports[0]['build']==reports[1]['build']
    assert len(reports[0][key])==len(reports[1][key])==count
    ulp=0;delta=0.
    for a,b in zip(reports[0][key],reports[1][key]):
        fields=['rate','channels','gain_db','frames','kind','file','lra_stable']
        if suite=='M9':fields+=['sample_peak_frames','true_peak_eighth_frames','ordered_roles','integrated_included','lra_included']
        for k in fields:assert a[k]==b[k],(suite,k)
        data=[]
        for root,row in zip(roots,(a,b)):
            assert Path(row['file']).name==row['file']
            raw=(root/folder/row['file']).read_bytes()
            assert len(raw)==row['frames']*row['channels']*4 and hashlib.sha256(raw).hexdigest()==row['sha256']
            data.append(raw)
        ulp=max(ulp,compare_capture(*data,a['gain_db']==0))
        for k,m in a['metrics'].items():
            n=b['metrics'][k];assert m['status']==n['status']
            if m['value'] is None:assert n['value'] is None
            else:
                assert math.isfinite(m['value']) and math.isfinite(n['value'])
                d=abs(m['value']-n['value']);assert d<=(.1 if k=='lra' else .01)
                delta=max(delta,d)
    result['suites'].append(dict(suite=suite,cases=count,maximum_float32_ulp=ulp,maximum_metric_difference=delta))
result['comparator_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
(r/'m7-m9-core-captures-candidate-03.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
