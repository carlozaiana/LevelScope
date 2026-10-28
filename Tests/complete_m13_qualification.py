#!/usr/bin/env python3
"""Final M13 gate: both actual platforms plus distinct physical/UI/cinema receipts."""
import argparse,json,pathlib
from collect_m1_environment import sha256
from run_m13_qualification import validate_runs,validate_report,validate_ui
from run_m13_muted_video import validate_physical
from compare_m13_reports import compare
from m13_evidence_identity import identities,select

def validate_receipt_files(folder,rows):
    for row in rows:
        if 'log_sha256' in row:
            if sha256(folder/(row['name']+'.log'))!=row['log_sha256']:raise ValueError('Changed run log')
        elif 'evidence_sha256' in row:
            if sha256(folder/(row['name']+'.json'))!=row['evidence_sha256']:raise ValueError('Changed internal check evidence')
        else:raise ValueError('Missing receipt integrity')

def complete(candidate,results,reference=None):
    identity=json.loads((candidate/'identity.json').read_text());evidence={}
    records,reuse=identities(candidate,reference)
    for platform,arch,prefix in (('sequoia','arm64','15.'),('high-sierra','x86_64','10.13.')):
        root=results/platform
        for phase in ('short','long','video','exchange'):
            folder=root/phase;s=json.loads((folder/'summary.json').read_text());env=json.loads((folder/'environment.json').read_text())
            evidence_identity,evidence_hash=select(records,s.get('source_digest'))
            if s.get('passed') is not True or s.get('phase')!=phase:raise ValueError('Missing/failed '+platform+' '+phase)
            if env.get('architecture')!=arch or not env.get('os','').startswith(prefix) or env.get('identity_sha256')!=evidence_hash or env.get('source_digest')!=s.get('source_digest'):raise ValueError('Wrong actual-platform identity')
            rows=json.loads((folder/'runs.json').read_text());validate_runs(rows,phase);validate_receipt_files(folder,rows)
            if s.get('report_hashes')!={p.name:sha256(p) for p in folder.glob('*.json') if p.name!='summary.json'}:raise ValueError('Phase evidence changed after completion')
            key={'short':'m13-core','long':'m13-long','video':'m13-video','exchange':'m13-exchange'}[phase]
            validate_report(json.loads((folder/(key+'.json')).read_text()),identity,arch,key)
            evidence[platform+'/'+phase]=sha256(folder/'summary.json')
        ui=json.loads((root/'ui-checklist.json').read_text());ui_identity,_=select(records,ui.get('source_digest'));validate_ui(ui,ui_identity,arch);evidence[platform+'/ui']=sha256(root/'ui-checklist.json')
        validate_physical(root/'physical-video',identity,arch)
        evidence[platform+'/physical-video']=sha256(root/'physical-video/summary.json')
        if reuse:
            repair=json.loads((root/'tool-repair.json').read_text())
            if repair.get('passed') is not True or repair.get('source_digest')!=identity['source_digest'] or repair.get('architecture')!=arch or repair.get('runner_tests_sha256')!=identity['artifacts_sha256']['Tests/M13QualificationRunnerTests.py'] or repair.get('log_sha256')!=sha256(root/'tool-repair.log'):raise ValueError('Missing fresh qualification-tool regression')
            evidence[platform+'/tool-repair']=sha256(root/'tool-repair.json')
    comparison=json.loads((results/'cross-platform-comparison.json').read_text())
    if comparison!=compare(candidate,results/'sequoia',results/'high-sierra',reference):raise ValueError('Missing or changed cross-platform comparison')
    cinema=results/'sequoia/cinema';s=json.loads((cinema/'summary.json').read_text())
    select(records,s.get('source_digest'))
    if s.get('passed') is not True:raise ValueError('Missing bounded cinema evaluation')
    validate_runs(json.loads((cinema/'runs.json').read_text()),'cinema')
    validate_receipt_files(cinema,json.loads((cinema/'runs.json').read_text()))
    validate_report(json.loads((cinema/'cinema.json').read_text()),identity,'arm64','cinema')
    if s.get('report_hashes')!={p.name:sha256(p) for p in cinema.glob('*.json') if p.name!='summary.json'}:raise ValueError('Cinema evidence changed')
    return dict(qualified=True,source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],identity_sha256=sha256(candidate/'identity.json'),evidence=evidence,tool_only_equivalence=reuse,
        publication_authorized=False,main_integration_authorized=False,main_remains='cfb1033aa7fb65a2655db276bee696161196bdeb',physical_full_71_listening='Separately pending; outside M13',picture_skips='Accepted development limitation under retained audio/drift/lifecycle conditions')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--candidate',type=pathlib.Path,required=True);p.add_argument('--results',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--equivalent-candidate',type=pathlib.Path);a=p.parse_args();a.output.write_text(json.dumps(complete(a.candidate,a.results,a.equivalent_candidate),indent=2)+'\n')
