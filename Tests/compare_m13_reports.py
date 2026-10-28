#!/usr/bin/env python3
"""M13 actual-platform comparison; tolerances are comparison-only, never QC allowances."""
import argparse,json,math,pathlib,struct
from collect_m1_environment import sha256
from compare_m11_reports import metric_delta,bits
from run_m13_qualification import validate_report
from m13_evidence_identity import identities,select

def load(path):return json.loads(path.read_text())
def compare(candidate,left,right,reference=None):
    identity=load(candidate/'identity.json');cores=[]
    records,reuse=identities(candidate,reference)
    for root,arch,prefix in ((left,'arm64','15.'),(right,'x86_64','10.13.')):
        env=load(root/'short/environment.json')
        evidence_identity,evidence_hash=select(records,env.get('source_digest'))
        if env.get('architecture')!=arch or not env.get('os','').startswith(prefix) or env.get('identity_sha256')!=evidence_hash:raise ValueError('Wrong actual-platform candidate')
        report=load(root/'short/m13-core.json');validate_report(report,identity,arch,'m13-core');cores.append({(r['rate'],r['channels'],r['mode']):r for r in report['cases']})
    cases=[];maximum=0;largest=0.
    for key,a in sorted(cores[0].items()):
        b=cores[1][key];rate,channels,mode=key
        if a['private_passes']!=b['private_passes']:raise ValueError('Different decision path')
        for metric in ('gain','integrated','true_peak'):
            delta=abs(a[metric]-b[metric])
            if not math.isfinite(delta) or delta>.01:raise ValueError('Metric/gain comparison exceeded .01')
            largest=max(largest,delta)
        row=dict(rate=rate,channels=channels,mode=mode,maximum_ulp=None)
        if rate==48000:
            folder=f'{rate}-{channels}-{mode}';captures=[];reports=[]
            for root,core in ((left,a),(right,b)):
                f=root/'short/m13-files'/folder/'audition-capture.f32'
                if sha256(f)!=core['capture_sha256'] or f.stat().st_size!=3*4096*channels*4:raise ValueError('Missing/tampered pre-monitor capture')
                captures.append(f.read_bytes());reports.append(load(root/'short/m13-files'/folder/'export/report.json'))
            peak=0
            for (x,),(y,) in zip(struct.iter_unpack('<f',captures[0]),struct.iter_unpack('<f',captures[1])):
                if not math.isfinite(x) or not math.isfinite(y):raise ValueError('Nonfinite capture')
                peak=max(peak,abs(bits(x)-bits(y)))
            if peak>1:raise ValueError('Processed audition capture exceeds one float32 ULP')
            maximum=max(maximum,peak);row['maximum_ulp']=peak
            for domain in ('source','current','artifact'):
                if reports[0][domain]['evaluation']['verdict']!=reports[1][domain]['evaluation']['verdict']:raise ValueError('Measured verdict differs')
                largest=max(largest,metric_delta(reports[0][domain]['metrics'],reports[1][domain]['metrics']))
                if reports[0][domain]['lra_stable']!=reports[1][domain]['lra_stable']:raise ValueError('LRA stability differs')
                px,py=reports[0][domain]['per_channel_peaks'],reports[1][domain]['per_channel_peaks']
                if len(px)!=channels or len(py)!=channels:raise ValueError('Missing all-channel peaks')
                for x,y in zip(px,py):
                    if x['role']!=y['role']:raise ValueError('Ordered roles differ')
                    largest=max(largest,metric_delta({k:x[k] for k in ('sample_peak','true_peak')},{k:y[k] for k in ('sample_peak','true_peak')}))
            pa,pb=reports[0]['current_processing_plan'],reports[1]['current_processing_plan']
            if pa['input_identity']!=pb['input_identity']:raise ValueError('Different programme inputs')
            if len(pa['applied_intervals'])!=len(pb['applied_intervals']):raise ValueError('Different envelope selection')
            for x,y in zip(pa['applied_intervals'],pb['applied_intervals']):
                if (x['begin_frame'],x['end_frame'])!=(y['begin_frame'],y['end_frame']) or abs(float.fromhex(x['gain_hex'])-float.fromhex(y['gain_hex']))>.01:raise ValueError('Envelope differs across actual platforms')
        cases.append(row)
    long_reports=[]
    for root,arch in ((left,'arm64'),(right,'x86_64')):
        validate_report(load(root/'long/m13-long.json'),identity,arch,'m13-long');long_reports.append(load(root/'long/long-files/export/report.json'))
    if long_reports[0]['composition']['definition_digest']!=long_reports[1]['composition']['definition_digest']:raise ValueError('Different sustained workload')
    for domain in ('source','current','artifact'):largest=max(largest,metric_delta(long_reports[0][domain]['metrics'],long_reports[1][domain]['metrics']))
    return dict(passed=True,source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],cases=cases,maximum_float32_ulp=maximum,maximum_metric_difference=largest,
        allowances='Cross-platform comparison only: I/M/S/TP/SP and gain .01; LRA .1; processed float32 <=1 ULP. No additional profile allowance.',physical_multichannel_listening='Separately pending',tool_only_equivalence=reuse)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--candidate',type=pathlib.Path,required=True);p.add_argument('--sequoia',type=pathlib.Path,required=True);p.add_argument('--high-sierra',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--equivalent-candidate',type=pathlib.Path);a=p.parse_args()
    try:result=compare(a.candidate,a.sequoia,a.high_sierra,a.equivalent_candidate)
    except (ValueError,KeyError,TypeError,OSError) as e:result=dict(passed=False,error=str(e))
    a.output.write_text(json.dumps(result,indent=2)+'\n');raise SystemExit(0 if result['passed'] else 1)
