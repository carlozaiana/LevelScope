"""Compare fresh actual-platform M18 receipts using established comparison allowances."""
import argparse
import json
import math
from pathlib import Path
import struct
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from compare_m11_reports import metric_delta,bits
from run_ebu_measurements import compare_results
from run_short import REQUIRED,sha

def read(path):return json.loads(path.read_text())
def compare(left,right,include_long):
    roots=[left,right];summary=[read(p/'short/summary.json') for p in roots]
    assert [s['architecture'] for s in summary]==['arm64','x86_64']
    assert summary[0]['os'].startswith('15.') and summary[1]['os'].startswith('10.13.')
    assert all(s['passed'] and {r['name'] for r in s['commands']}==REQUIRED for s in summary)
    assert summary[0]['identity_sha256']==summary[1]['identity_sha256'] and summary[0]['build']==summary[1]['build']
    ebu=[read(p/'short/ebu.json') for p in roots];assert ebu[0]['manifest_sha256']==ebu[1]['manifest_sha256']
    assert not compare_results(ebu[0]['results'],ebu[1])
    itu=[read(p/'short/itu.json') for p in roots];assert itu[0]['manifest_sha256']==itu[1]['manifest_sha256']
    assert len(itu[0]['results'])==len(itu[1]['results'])
    maximum=0.;ulp=0;cases=[]
    for a,b in zip(itu[0]['results'],itu[1]['results']):
        assert (a['file'],a['sha256'],a['passed'])==(b['file'],b['sha256'],True)
        for key in ('I','M','S','TP','LRA'):
            x,y=a['measured'][key],b['measured'][key];assert (x is None)==(y is None)
            if x is not None:assert math.isfinite(x) and math.isfinite(y) and abs(x-y)<=(.1 if key=='LRA' else .01);maximum=max(maximum,abs(x-y))
    names=['native-stereo','native-surround']
    pairs=[(name,[read(p/'short'/name/'export/report.json') for p in roots]) for name in names]
    if include_long:
        for p in roots:
            r=read(p/'long/result.json');audit=read(p/'long-audit.json');assert r['passed'] and audit['passed'] and audit['three_hour'] and audit['rf64'] and audit['samples_beyond_4_gib']>0
        pairs.append(('three-hour-mono',[read(p/'long/export/report.json') for p in roots]))
    for name,reports in pairs:
        a,b=reports
        assert a['ordered_input_roles']==b['ordered_input_roles'] and a['source']['measured_content_sha256']==b['source']['measured_content_sha256']
        for domain in ('source','current','artifact'):
            x,y=a[domain],b[domain];assert x['evaluation']['verdict']==y['evaluation']['verdict'] and x['lra_stable']==y['lra_stable']
            maximum=max(maximum,metric_delta(x['metrics'],y['metrics']))
            assert len(x['per_channel_peaks'])==len(y['per_channel_peaks'])
            for xp,yp in zip(x['per_channel_peaks'],y['per_channel_peaks']):
                assert all(xp[k]==yp[k] for k in (('channel','role') if domain=='artifact' else ('channel','role','sample_peak_frame','true_peak_eighth_frame')))
                maximum=max(maximum,metric_delta({k:xp[k] for k in ('sample_peak','true_peak')},{k:yp[k] for k in ('sample_peak','true_peak')}))
        x,y=a['current_processing_plan'],b['current_processing_plan'];assert x['input_identity']==y['input_identity'] and x['transitions']==y['transitions']
        assert len(x['applied_intervals'])==len(y['applied_intervals'])
        assert abs(float.fromhex(x['global_gain_hex'])-float.fromhex(y['global_gain_hex']))<=.01
        for xp,yp in zip(x['applied_intervals'],y['applied_intervals']):
            assert (xp['begin_frame'],xp['end_frame'])==(yp['begin_frame'],yp['end_frame'])
            assert abs(float.fromhex(xp['gain_hex'])-float.fromhex(yp['gain_hex']))<=.01
        cases.append(dict(name=name,source_current_peak_positions_equal=True,verdicts_equal=True,section_geometry_equal=True))
    # Existing M13 short captures retain samplewise processed audition parity.
    for channels in (2,6):
        files=[p/'short/section-proposal-files'/('48000-'+str(channels)+'-0')/'audition-capture.f32' for p in roots]
        data=[f.read_bytes() for f in files];assert len(data[0])==len(data[1])==3*4096*channels*4
        for (x,),(y,) in zip(struct.iter_unpack('<f',data[0]),struct.iter_unpack('<f',data[1])):
            assert math.isfinite(x) and math.isfinite(y);distance=abs(bits(x)-bits(y));assert distance<=1;ulp=max(ulp,distance)
    return dict(passed=True,build=summary[0]['build'],identity_sha256=summary[0]['identity_sha256'],maximum_metric_difference=maximum,
                maximum_capture_float32_ulp=ulp,cases=cases,three_hour_included=include_long,
                allowances='Cross-Mac comparisons only: 0.01 LU/dB/gain, LRA 0.1 LU, processed float32 <=1 ULP; no additional compliance allowance')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sequoia',type=Path);p.add_argument('high_sierra',type=Path);p.add_argument('--long',action='store_true');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    try:report=compare(a.sequoia,a.high_sierra,a.long)
    except Exception as e:report=dict(passed=False,error=repr(e))
    a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2));raise SystemExit(0 if report['passed'] else 1)
