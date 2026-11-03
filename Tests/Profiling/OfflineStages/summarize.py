"""Summarize warmed, complete trials only. All durations are milliseconds."""
import argparse,hashlib,json,math,statistics as st
from pathlib import Path
def percentile(values,p):
    return sorted(values)[math.ceil(len(values)*p)-1]
def stats(values):
    return dict(median=st.median(values),min=min(values),max=max(values))
def summarize(root,n):
    groups={mode:[json.loads(p.read_text()) for p in sorted((root/'runs').glob(f'trial-*-{n}-{mode}.json'))] for mode in ('control','profile')}
    base=groups['control'][0]
    for rows in groups.values():
        assert len(rows)==3
        assert all(d['identities']==base['identities'] and d['base_build']==base['base_build'] for d in rows)
    operations={}
    for name in base['operations']:
        ctrl, probes=([d['operations'][name] for d in groups[k]] for k in ('control','profile'))
        total=stats([o['wall_ms'] for o in ctrl]); walls=stats([o['wall_ms'] for o in probes])
        stages={}
        for phase, fields in probes[0]['timings_ms'].items():
            stages[phase]={key:st.median(o['timings_ms'][phase][key] for o in probes) for key in fields}
        item=dict(control_wall_ms=total,profile_wall_ms=walls,instrumented_delta_pct=100*(walls['median']/total['median']-1),
                  stages_ms=stages,control_worker_core_equivalents=st.median(o['worker_cpu_ms']/o['wall_ms'] for o in ctrl))
        if name in base['outcomes']:
            item['passes']=base['outcomes'][name]['passes']
            item['with_preflight_accept_ms']=stats([d['operations'][name]['wall_ms']+d['outcomes'][name]['main_preflight_ms']+
                                                      d['outcomes'][name]['main_accept_ms'] for d in groups['control']])
        operations[name]=item
    long_jobs=[d['operations'][j] for d in groups['control'] for j in ('source','fresh','reused','current')]
    gaps=[g for o in long_jobs for g in o['heartbeat_gaps_ms']]
    paints=[v for o in long_jobs for v in o['timeline_paint_ms']]
    ui=dict(heartbeat_count=len(gaps),gap_p95_ms=percentile(gaps,.95),
            max_gap_or_tail_ms=max(max(gaps),max(o['heartbeat_tail_ms'] for o in long_jobs)),paint_p95_ms=percentile(paints,.95))
    return dict(operations=operations,ui=ui,base_build=base['base_build'],input_sha256=base['file_sha256'],roles=base['roles'],
                sample_stacks=json.loads((root/'samples'/f'hotspots-{n}.json').read_text()))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('evidence',type=Path);p.add_argument('bin',type=Path);args=p.parse_args()
    report={str(n):summarize(args.evidence,n) for n in (2,6)}
    report['binary_sha256']={name:hashlib.sha256((args.bin/name).read_bytes()).hexdigest() for name in ('control','profile','sample-control')}
    (args.evidence/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
