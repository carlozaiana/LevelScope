#!/usr/bin/env python3
"""Close only M12's approved Sequoia phase; never declare full qualification."""
import argparse, json, pathlib
from collect_m1_environment import sha256
from run_m12_qualification import validate_runs, validate_report, validate_ui

def complete(package,short,long,ui):
    identity=json.loads((package/'identity.json').read_text())
    for name,digest in identity['artifacts_sha256'].items():
        if sha256(package/name)!=digest:raise ValueError('Frozen artifact changed: '+name)
    evidence={}
    for phase,folder,report in (('short',short,'core.json'),('long',long,'long.json')):
        summary=json.loads((folder/'summary.json').read_text())
        if summary.get('sequoia_automated_passed') is not True or summary.get('source_digest')!=identity['source_digest'] or summary.get('phase')!=phase:raise ValueError('Missing/failed/wrong-candidate '+phase)
        runs=json.loads((folder/'runs.json').read_text());validate_runs(runs,phase)
        for run in runs:
            suffix,key=('.log','log_sha256') if 'log_sha256' in run else ('.json','evidence_sha256')
            if sha256(folder/(run['name']+suffix))!=run[key]:raise ValueError('Missing/changed required receipt: '+run['name'])
        validate_report(json.loads((folder/report).read_text()),identity,'sequoia-audio' if phase=='short' else 'sequoia-long')
        environment=json.loads((folder/'environment.json').read_text())
        if environment.get('identity_sha256')!=sha256(package/'identity.json') or environment.get('architecture')!='arm64' or not environment.get('os','').startswith('15.') or environment.get('video_loaded') is not False or environment.get('physical_devices_opened') is not False:raise ValueError('Invalid Sequoia environment')
        for p in folder.glob('*.json'):evidence[str(p.resolve())]=sha256(p)
    validate_ui(json.loads(ui.read_text()),identity);evidence[str(ui.resolve())]=sha256(ui)
    return dict(sequoia_phase_complete=True,qualification_complete=False,
      status='Sequoia tested; High Sierra and video qualification pending',
      source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],
      identity_sha256=sha256(package/'identity.json'),evidence_sha256=evidence,
      deferred=['Actual Intel High Sierra runtime and cross-platform qualification after M13','All selected video loading, decoding and playback qualification after M13'],
      publication_authorized=False,main_integration_authorized=False,m13_implementation_authorized=False)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('package','short','long','ui','output'):p.add_argument('--'+name,required=True,type=pathlib.Path)
    a=p.parse_args();result=complete(a.package,a.short,a.long,a.ui);a.output.write_text(json.dumps(result,indent=2)+'\n')
