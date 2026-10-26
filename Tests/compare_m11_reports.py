#!/usr/bin/env python3
"""Actual-platform M11 identity, capture, metric and required-phase comparison."""
import argparse,json,math,pathlib,struct
from collect_m1_environment import sha256
from run_m11_qualification import validate_report,validate_runs,validate_ui

def bits(x):
    n=struct.unpack('<I',struct.pack('<f',x))[0];return 0x80000000-n if n&0x80000000 else n+0x80000000
def metric_delta(a,b):
    largest=0.
    if set(a)!=set(b):raise ValueError('Metric inventory mismatch')
    for key,x in a.items():
        y=b[key]
        if x['status']!=y['status'] or (x['value'] is None)!=(y['value'] is None):raise ValueError('Availability mismatch')
        if x['value'] is not None:
            d=abs(x['value']-y['value'])
            if not math.isfinite(d) or d>(.1 if key=='lra' else .01):raise ValueError('Cross-platform metric tolerance exceeded')
            largest=max(largest,d)
    return largest
def compare(left,right,complete=False):
    roots=(left,right);ids=[json.loads((p/'short/identity.json').read_text()) for p in roots]
    if ids[0]!=ids[1]:raise ValueError('Different frozen candidates')
    cores=[]
    for root,arch in zip(roots,('arm64','x86_64')):
        for phase in ('short','long','exchange') if complete else ('short',):
            p=root/phase;env=json.loads((p/'environment.json').read_text());summary=json.loads((p/'summary.json').read_text());runs=json.loads((p/'runs.json').read_text())
            if env['architecture']!=arch or not env['os'].startswith('15.' if arch=='arm64' else '10.13.') or env['artifact_integrity_passed'] is not True or env['source_digest']!=ids[0]['source_digest'] or summary['automated_passed'] is not True:raise ValueError('Actual-platform or required phase failed')
            validate_runs(runs,phase)
            if phase=='short':
                core=json.loads((p/'core.json').read_text());validate_report(core,ids[0],arch,'core');cores.append(core);validate_report(json.loads((p/'video.json').read_text()),ids[0],arch,'video')
            else:validate_report(json.loads((p/(phase+'.json')).read_text()),ids[0],arch,phase)
        if complete:validate_ui(json.loads((root/'ui.json').read_text()),ids[0],arch)
    rows=[];maximum=0;metric=0.
    for x,y in zip(cores[0]['cases'],cores[1]['cases']):
        for key in ('kind','rate','channels','gain_db','frames','file','lra_stable'):
            if x[key]!=y[key]:raise ValueError('Capture case mismatch: '+key)
        samples=[]
        for p,row in zip(roots,(x,y)):
            file=p/'short/samples'/row['file']
            if pathlib.Path(row['file']).name!=row['file'] or sha256(file)!=row['sha256'] or file.stat().st_size!=row['frames']*row['channels']*4:raise ValueError('Capture integrity mismatch')
            samples.append(file.read_bytes())
        peak=0
        for (a,),(b,) in zip(struct.iter_unpack('<f',samples[0]),struct.iter_unpack('<f',samples[1])):
            if not math.isfinite(a) or not math.isfinite(b):raise ValueError('Non-finite capture')
            d=abs(bits(a)-bits(b));peak=max(peak,d)
            if d>(0 if x['gain_db']==0 else 1):raise ValueError('Pre-monitor capture exceeds float32 parity allowance')
        delta=metric_delta(x['metrics'],y['metrics']);metric=max(metric,delta);maximum=max(maximum,peak);rows.append(dict(case=x['file'],maximum_ulp=peak,maximum_metric_difference=delta))
    if complete:
        longReports=[json.loads((p/'long/long-files/export/report.json').read_text()) for p in roots]
        if longReports[0]['composition']['definition_digest']!=longReports[1]['composition']['definition_digest']:raise ValueError('Different long composition')
        for domain in ('source','current','artifact'):metric=max(metric,metric_delta(longReports[0][domain]['metrics'],longReports[1][domain]['metrics']))
    return dict(passed=True,qualification_complete=complete,source_digest=ids[0]['source_digest'],maximum_float32_ulp=maximum,maximum_metric_difference=metric,cases=rows,physical_multichannel_listening='separately pending',allowances='Cross-platform comparison only; no extra profile acceptance allowance')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('sequoia',type=pathlib.Path);p.add_argument('high_sierra',type=pathlib.Path);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--complete',action='store_true');a=p.parse_args()
    try:r=compare(a.sequoia,a.high_sierra,a.complete)
    except (AssertionError,ValueError,KeyError,TypeError,OSError) as e:r=dict(passed=False,qualification_complete=False,error=str(e))
    a.output.write_text(json.dumps(r,indent=2)+'\n');raise SystemExit(0 if r['passed'] else 1)
