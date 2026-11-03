"""Compare the two actual service policies on identical PCM, with fixed runtime identities."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import statistics
import subprocess

parser=argparse.ArgumentParser()
parser.add_argument('stage',choices=['decisions','timings','enabled'])
for arg in ('binaries','fixtures','destination'): parser.add_argument(arg,type=Path)
parser.add_argument('--serial-reference',type=Path,help='optional fresh spot checks of the preserved serial executable')
args=parser.parse_args()
root=args.destination.resolve()
dest=root/args.stage
dest.mkdir(parents=True,exist_ok=False)
expected={}

def signature(data):
    decisions={name:{k:v for k,v in row.items() if k not in ('main_preflight_ms','main_accept_ms')}
               for name,row in data['outcomes'].items()}
    return (data['identities'],decisions,data['file_sha256'],data['source_lufs'],data['current_lufs'],data['base_build'])

def run(ch,mode,label,extended=False):
    f=dest/f'{label}-{ch}-{mode}.json'
    cmd=[str(args.binaries.resolve()/mode),str(args.fixtures.resolve()/f'programme-{ch}.wav'),str(f)]
    if extended: cmd+=['--decisions']
    with f.with_suffix('.log').open('w') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    d=json.loads(f.read_text())
    s=signature(d)
    assert s==expected.setdefault(ch,s), f'Serial parity mismatch in {f}'
    print(f.stem,*(f'{k}={r["wall_ms"]/1000:.3f}s' for k,r in d['operations'].items()),flush=True)
    return d

def p95(items):
    return sorted(items)[math.ceil(len(items)*.95)-1]

if args.stage in ('decisions','enabled'):
    for ch in (2,6):
        a=run(ch,'serial','compare',True)
        b=run(ch,'app-mode' if args.stage=='enabled' else 'bounded','compare',True)
        assert len(a['outcomes'])==7
        if args.stage=='enabled':
            for k in ('source','fresh','reused','current'):
                assert b['operations'][k]['wall_ms'] <= .90*a['operations'][k]['wall_ms'],f'app path did not accelerate {k}'
            source=b['operations']['source']
            assert source['process_cpu_ms']/source['wall_ms']>1.15, 'normal service did not show concurrent execution'
    (dest/'passed.json').write_text(json.dumps({'exact_match':True,'channels':[2,6],'cases':list(a['outcomes'])},indent=2)+'\n')
else:
    assert json.loads((root/'decisions/passed.json').read_text())['exact_match']
    data={n:{mode:[] for mode in ('serial','bounded')} for n in (2,6)}
    for n in (2,6):
        for mode in ('serial','bounded'): run(n,mode,'warmup')
    for i in range(1,4):
        for n in ((2,6) if i%2 else (6,2)):
            for mode in (('serial','bounded') if (n==2)==(i%2==1) else ('bounded','serial')):
                data[n][mode].append(run(n,mode,f'trial-{i}'))
    summary={}
    for n in (2,6):
        rows={}; ui={}; limits=True
        ref=None
        if args.serial_reference:
            ref=json.loads((args.serial_reference/f'prior-{n}.json').read_text())
            baseline=data[n]['serial'][0]
            for field in ('file_sha256','source_lufs','current_lufs'):
                assert ref[field]==baseline[field],field
            for plan in ('fresh_plan','reused_plan'): assert ref['identities'][plan]==baseline['identities'][plan]
        for k in data[n]['serial'][0]['operations']:
            entry={}
            for mode,reps in data[n].items():
                vals=[d['operations'][k]['wall_ms'] for d in reps]
                if k=='global':
                    vals=[v+d['outcomes'][k]['main_preflight_ms']+d['outcomes'][k]['main_accept_ms'] for v,d in zip(vals,reps)]
                core=[d['operations'][k]['process_cpu_ms']/d['operations'][k]['wall_ms'] for d in reps]
                entry[mode]={'ms':statistics.median(vals),'range_ms':[min(vals),max(vals)],'process_cores':statistics.median(core)}
            entry['speedup']=entry['serial']['ms']/entry['bounded']['ms']
            rows[k]=entry
            if k in ('source','fresh','reused','current'):limits &= entry['bounded']['ms'] <= .90*entry['serial']['ms']
            if ref and k in ('source','fresh','reused','current'):
                entry['prior_binary_ms']=ref['operations'][k]['wall_ms']
                limits &= entry['serial']['ms']<=1.10*entry['prior_binary_ms'] and entry['bounded']['ms']<=.90*entry['prior_binary_ms']
        for mode,reps in data[n].items():
            jobs=[d['operations'][k] for d in reps for k in ('source','current','fresh','reused')]
            gaps=[x for j in jobs for x in j['heartbeat_gaps_ms']]
            ui[mode]={'p95_ms':p95(gaps),'max_ms':max(max(gaps),max(j['heartbeat_tail_ms'] for j in jobs))}
        limits &= rows['global']['bounded']['ms']<10 and ui['bounded']['p95_ms']<=ui['serial']['p95_ms']+20 and ui['bounded']['max_ms']<=200
        summary[n]={'operations':rows,'ui':ui,'gates_pass':limits,'identities_match':True}
    summary['binaries']={mode:hashlib.sha256((args.binaries/mode).read_bytes()).hexdigest() for mode in ('serial','bounded')}
    summary['build']=data[2]['serial'][0]['base_build']
    (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    assert all(summary[n]['gates_pass'] for n in (2,6)), 'Review performance gate before enabling'
