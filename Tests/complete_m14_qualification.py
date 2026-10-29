#!/usr/bin/env python3
"""Sequoia-phase gate; never claims full M14 platform/video/duration qualification."""
import argparse,json,pathlib
from collect_m1_environment import sha256
from run_m14_qualification import validate_runs,validate_report,validate_ui,validate_legacy,DEFERRED
def complete(candidate,results,tested_candidate=None):
    identity=json.loads((candidate/'identity.json').read_text());evidence={}
    for name,h in identity['artifacts_sha256'].items():
        if sha256(candidate/name)!=h:raise ValueError('Changed frozen artifact '+name)
    tested_candidate=tested_candidate or candidate
    tested=json.loads((tested_candidate/'identity.json').read_text())
    allowed={'Tests/run_m14_qualification.py','Tests/complete_m14_qualification.py','Tests/M14QualificationRunnerTests.py'}
    changed={k for k in set(identity['source_files_sha256'])|set(tested['source_files_sha256']) if identity['source_files_sha256'].get(k)!=tested['source_files_sha256'].get(k)}
    if not changed<=allowed or identity['measurement_build']!=tested['measurement_build']:raise ValueError('Reuse requires identical production and C++ test implementation')
    for name,h in tested['artifacts_sha256'].items():
        if sha256(tested_candidate/name)!=h:raise ValueError('Changed original tested artifact')
        if name.startswith(('bin/','legacy/','LevelScope Standalone.app/')) and identity['artifacts_sha256'].get(name)!=h:raise ValueError('Reuse binary mismatch')
    for phase in ('short','long','cinema','legacy'):
        folder=results/phase;s=json.loads((folder/'summary.json').read_text());env=json.loads((folder/'environment.json').read_text());rows=json.loads((folder/'runs.json').read_text())
        expected=identity if phase=='legacy' else tested;package=candidate if phase=='legacy' else tested_candidate
        if s.get('passed') is not True or s.get('source_digest')!=expected['source_digest'] or s.get('phase')!=phase or s.get('deferred')!=DEFERRED:raise ValueError('Missing/failed phase or deferral')
        if env.get('architecture')!='arm64' or not env.get('os','').startswith('15.') or env.get('identity_sha256')!=sha256(package/'identity.json') or env.get('video_loaded') is not False or env.get('physical_devices_opened') is not False:raise ValueError('Wrong Sequoia environment/candidate')
        validate_runs(rows,phase)
        for row in rows:
            key,ext=('log_sha256','.log') if 'log_sha256' in row else ('evidence_sha256','.json')
            if row.get(key)!=sha256(folder/(row['name']+ext)):raise ValueError('Changed/missing run evidence')
        if s.get('report_hashes')!={p.name:sha256(p) for p in folder.glob('*.json') if p.name!='summary.json'}:raise ValueError('Phase evidence changed')
        if phase=='legacy':validate_legacy(json.loads((folder/'legacy-rejection.json').read_text()),json.loads((folder/'legacy-baseline.json').read_text()))
        else:
            key={'short':'m14-core','long':'m14-long','cinema':'m14-cinema'}[phase];validate_report(json.loads((folder/(key+'.json')).read_text()),expected,key)
        evidence[phase]=sha256(folder/'summary.json')
    ui=results/'ui-checklist.json';validate_ui(json.loads(ui.read_text()),tested);evidence['ui']=sha256(ui)
    return dict(sequoia_phase_complete=True,qualified=False,status='Sequoia tested; High Sierra, video and full-duration mono/RF64 qualification pending',source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],identity_sha256=sha256(candidate/'identity.json'),evidence=evidence,reused_from_source=tested['source_digest'],tool_only_differences=sorted(changed),production_and_test_binaries_identical=True,deferred=DEFERRED,combined_window='After M18; no later than M19; later milestones require separate approval',main_remains='a9be070ab9eb5cb4f591ff9d991d8999c3a016d0',publication_authorized=False,main_integration_authorized=False)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--candidate',type=pathlib.Path,required=True);p.add_argument('--tested-candidate',type=pathlib.Path);p.add_argument('--results',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args();a.output.write_text(json.dumps(complete(a.candidate,a.results,a.tested_candidate),indent=2)+'\n')
