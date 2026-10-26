#!/usr/bin/env python3
"""Compare M10 actual-platform measures and fixed-seed sample oracles."""
import argparse,json,math,pathlib,struct

def compare_peaks(left,right,identical_audio):
    assert len(left)==len(right)
    maximum=0.
    for x,y in zip(left,right):
        assert x['role']==y['role']
        for name in ('sample_peak','true_peak'):
            m,n=x[name],y[name]
            assert m['status']==n['status'] and m['unit']==n['unit']
            if m['value'] is None:assert n['value'] is None
            else:
                assert math.isfinite(m['value']) and math.isfinite(n['value'])
                difference=abs(m['value']-n['value']);assert difference<=.01
                maximum=max(maximum,difference)
        if identical_audio:
            for name in ('sample_peak_frame','true_peak_eighth_frame'):assert x[name]==y[name]
    return maximum

def ordered_bits(x):
    n=struct.unpack('<I',struct.pack('<f',x))[0]
    return (0x80000000-n) if n&0x80000000 else (n+0x80000000)

def compare(a,b):
    ra=json.loads((a/'export.json').read_text());rb=json.loads((b/'export.json').read_text())
    assert ra['passed'] and rb['passed'] and ra['build']==rb['build'] and {ra['architecture'],rb['architecture']}=={'arm64','x86_64'}
    assert len(ra['cases'])==len(rb['cases'])==16
    max_diff=0.;different_hashes=0
    for x,y in zip(ra['cases'],rb['cases']):
        for k in ('domain','frames','sample_rate','input_kind','container','bits','channel_mask','verification_gain_db','dither_applications','output_to_input_channels'):assert x[k]==y[k],k
        assert x['evaluation']['verdict']==y['evaluation']['verdict']=='PASS'
        different_hashes+=x['sha256']!=y['sha256']
        for name,m in x['metrics'].items():
            n=y['metrics'][name];assert m['status']==n['status']
            if m['value'] is None:assert n['value'] is None
            else:
                diff=abs(m['value']-n['value']);max_diff=max(diff,max_diff);assert diff<=(.1 if name=='lra' else .01),(name,diff)
        # Quantized artifacts may have different seeds when normalized gain differs.
        # Locations are exact for identical content. Peak dB values retain the
        # documented cross-platform metric budget, including log/interpolation rounding.
        max_diff=max(max_diff,compare_peaks(x['per_channel_peaks'],y['per_channel_peaks'],x['sha256']==y['sha256']))
    rows=json.loads((a/'artifacts/raw-oracles.json').read_text());other=json.loads((b/'artifacts/raw-oracles.json').read_text());assert rows==other
    max_ulp=0;max_code=0
    for row in rows:
        x=(a/'artifacts'/row['capture']).read_bytes();y=(b/'artifacts'/row['capture']).read_bytes();assert len(x)==len(y)
        for (u,),(v,) in zip(struct.iter_unpack('<f',x),struct.iter_unpack('<f',y)):
            d=abs(ordered_bits(u)-ordered_bits(v));max_ulp=max(max_ulp,d);assert d<=(0 if row['gain_db']==0 else 1)
        x=(a/'artifacts'/row['file']).read_bytes()[68:];y=(b/'artifacts'/row['file']).read_bytes()[68:];assert len(x)==len(y)
        for i in range(0,len(x),3):
            d=abs(int.from_bytes(x[i:i+3],'little',signed=True)-int.from_bytes(y[i:i+3],'little',signed=True));max_code=max(max_code,d);assert d<=1
    return dict(passed=True,cases=16,maximum_metric_difference=max_diff,distinct_artifact_hashes=different_hashes,maximum_prequant_float_ulp=max_ulp,maximum_same_seed_pcm24_code_difference=max_code,allowances='Metric comparisons only; never extra profile acceptance allowance')
def compare_long(a,b):
    ra=json.loads((a/'long-export.json').read_text());rb=json.loads((b/'long-export.json').read_text())
    assert ra['passed'] and rb['passed'] and ra['build']==rb['build'] and {ra['architecture'],rb['architecture']}=={'arm64','x86_64'}
    a=json.loads((a/'long-export/report.json').read_text());b=json.loads((b/'long-export/report.json').read_text())
    maxima={}
    for domain in ('source','current','artifact'):
        x,y=a[domain],b[domain];assert x['frames']==y['frames']=='518400000' and x['sample_rate']==y['sample_rate']==48000
        assert x['evaluation']['verdict']==y['evaluation']['verdict']
        if domain!='source':assert x['evaluation']['verdict']=='PASS'
        delta=0.
        for name,m in x['metrics'].items():
            n=y['metrics'][name];assert m['status']==n['status']
            if m['value'] is None:assert n['value'] is None
            else:
                d=abs(m['value']-n['value']);delta=max(d,delta);assert d<=(.1 if name=='lra' else .01),(domain,name,d)
        assert [p['role'] for p in x['per_channel_peaks']]==[p['role'] for p in y['per_channel_peaks']]
        maxima[domain]=delta
    return dict(passed=True,frames=518400000,channels=6,sample_rate=48000,maximum_metric_difference_by_domain=maxima,artifact_sha256=[a['artifact']['sha256'],b['artifact']['sha256']])
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sequoia',type=pathlib.Path);p.add_argument('high_sierra',type=pathlib.Path);p.add_argument('--long',action='store_true');p.add_argument('--output',required=True,type=pathlib.Path);a=p.parse_args();r=compare_long(a.sequoia,a.high_sierra) if a.long else compare(a.sequoia,a.high_sierra);a.output.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
