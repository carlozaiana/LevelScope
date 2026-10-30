#!/usr/bin/env python3
"""Explicit short audio-only M16 checks; no recursive runners/video/device/Golem."""
import argparse,hashlib,json,pathlib,platform,struct,subprocess,time
LEGACY_BUILD='8a7a2866d48f7ba6d4db1c7f532eaffbfce03fa78fbca02211d3f6b9fb79ef47'
REQUIRED={'preferences','global-proposals','sessions-short','sections','m12-smoke','m13-smoke','markers','export-short','m15-reader-rejects-schema6','runner-self-tests'}
GROUPS={'policy-validation','support-and-gain-limits','capacity-and-complete-inspection','ui-controls-and-responsiveness','identity-session-report','large-measured-proposal','mono-policy-path'}
def sha(f):
    h=hashlib.sha256()
    with open(f,'rb') as s:
        for b in iter(lambda:s.read(1048576),b''):h.update(b)
    return h.hexdigest()
def validate(rows,build):
    if len(rows)!=len(REQUIRED) or {r['name'] for r in rows}!=REQUIRED:raise ValueError('Missing/duplicate/unexpected required M16 check')
    for r in rows:
        if r['exit_code']!=r.get('expected_exit',0):raise ValueError('Failed command '+r['name'])
        if r['name']=='m15-reader-rejects-schema6':
            if (r.get('passed') is not False or r.get('reason')!='foreign session has exact applied plan and provenance'
                or r.get('build')!=LEGACY_BUILD or r.get('legacy_input_schema')!=6 or len(r.get('legacy_input_sha256',''))!=64):raise ValueError('Legacy reader rejection not established')
        elif r['name']!='runner-self-tests' and (r.get('passed') is not True or r.get('build')!=build):raise ValueError('Failed/mismatched receipt '+r['name'])
        if r['name']=='preferences':
            groups=r.get('groups',[])
            if len(groups)!=len(GROUPS) or set(groups)!=GROUPS:raise ValueError('Missing/duplicate M16 acceptance group')
def main():
    p=argparse.ArgumentParser();p.add_argument('--candidate',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args();root=a.candidate.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    identity=json.loads((root/'identity.json').read_text());assert platform.system()=='Darwin' and platform.machine()=='arm64'
    for f,h in identity['artifacts_sha256'].items():assert sha(root/f)==h,f
    rows=[]
    def run(name,args,report=None,expected=0):
        started=time.monotonic();command=[str(x) for x in args]
        with (out/(name+'.log')).open('w') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        row=dict(name=name,command=command,exit_code=r.returncode,expected_exit=expected,seconds=time.monotonic()-started)
        if report and pathlib.Path(report).is_file():
            d=json.loads(pathlib.Path(report).read_text());row.update(passed=d.get('passed'),build=d.get('build',d.get('source_build')),reason=d.get('message',d.get('failure',d.get('error',''))),report_sha256=sha(report))
            if name=='preferences':row['groups']=d.get('policy_groups',[])
        if name=='m15-reader-rejects-schema6':
            fixture=out/'m13-files/48000-2-0/applied.levelscope';raw=fixture.read_bytes();assert raw[:8]==b'LSSESS01'
            size=struct.unpack_from('<Q',raw,8)[0];core=raw[80:80+size];assert hashlib.sha256(core).hexdigest()==raw[16:80].decode()
            document=json.loads(core);assert document['schema']==6 and document['current']['sectionGains'] and document['caches'] and len(raw)>80+size
            row.update(legacy_input_schema=6,legacy_input_sha256=sha(fixture))
        rows.append(row);(out/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
        if r.returncode!=expected:raise RuntimeError('Required command failed: '+name)
    def binary(name):return ['/usr/bin/arch','-arm64',root/'bin'/name]
    start=time.monotonic();failure=''
    try:
        run('preferences',binary('LevelScopeSectionProposalTests')+['--preferences','--root',out/'preferences-files','--output',out/'preferences.json'],out/'preferences.json')
        run('global-proposals',binary('LevelScopeProposalTests')+['--report',out/'global-proposals.json'],out/'global-proposals.json')
        run('sessions-short',binary('LevelScopeSessionTests')+['--short','--output',out/'sessions.json'],out/'sessions.json')
        run('sections',binary('LevelScopeSectionsTests')+['--report',out/'sections.json'],out/'sections.json')
        run('m12-smoke',binary('LevelScopeSectionProcessingTests')+['--smoke','--root',out/'m12-files','--output',out/'m12.json'],out/'m12.json')
        run('m13-smoke',binary('LevelScopeSectionProposalTests')+['--smoke','--root',out/'m13-files','--output',out/'m13.json'],out/'m13.json')
        run('markers',binary('LevelScopeMarkerTests')+[out/'marker-files',root/'source/Tests/Markers/nuendo-15.0.30-points.csv'],out/'marker-files/marker-tests.json')
        run('export-short',binary('LevelScopeExportTests')+['--artifacts',out/'export-files','--output',out/'export.json'],out/'export.json')
        run('m15-reader-rejects-schema6',['/usr/bin/arch','-arm64',root/'legacy/M15SectionProposalTests','--exchange',out/'m13-files','--root',out/'legacy-read-files','--output',out/'legacy-rejection.json'],out/'legacy-rejection.json',1)
        run('runner-self-tests',['python3','-B',root/'Tests/M16QualificationRunnerTests.py'])
        validate(rows,identity['measurement_build'])
    except Exception as e:failure=str(e)
    summary=dict(short_automated_passed=not failure,failure=failure,source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],
      os=subprocess.check_output(['sw_vers'],text=True).strip(),cpu=subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),architecture=platform.machine(),seconds=time.monotonic()-start,commands=rows,
      silent_UI='separate required candidate receipt',physical_devices_opened=False,video_loaded=False,high_sierra='deferred after M18, no later than M19',three_hour='deferred; short structural test is not full-duration audio',qualification_complete=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='commands'},indent=2));return 1 if failure else 0
if __name__=='__main__':raise SystemExit(main())
