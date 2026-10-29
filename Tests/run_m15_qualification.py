#!/usr/bin/env python3
"""Explicit short M15 Sequoia selection. No recursive milestone runners/video/Golem."""
import argparse,hashlib,json,pathlib,platform,shutil,subprocess,time
REQUIRED={'markers','sessions-short','sections','m12-smoke','m13-smoke','export-short','m14-reader-rejects-schema5','runner-self-tests'}
REQUIRED_MARKER_GROUPS={"exact-mapping","bounded-readers","nuendo-actual-structure","pro-tools-actual-structure","resource-limits","ui-resource","legacy schemas / cache recovery","transaction/session/report 1","transaction/session/report 2","transaction/session/report 6"}
def validate_marker_groups(groups):
    if len(groups)!=len(REQUIRED_MARKER_GROUPS) or set(groups)!=REQUIRED_MARKER_GROUPS:raise ValueError("Missing, duplicate or unexpected M15 marker acceptance group")
REUSABLE={"sessions-short":"LevelScopeSessionTests","sections":"LevelScopeSectionsTests","m12-smoke":"LevelScopeSectionProcessingTests","m13-smoke":"LevelScopeSectionProposalTests","export-short":"LevelScopeExportTests"}
def validate_reuse(current,previous,target):
    code=lambda d:{k:v for k,v in d["source_files_sha256"].items() if k.startswith("Source/") or k=="CMakeLists.txt"}
    key="bin/"+target
    if code(current)!=code(previous) or current["measurement_build"]!=previous["measurement_build"] or current["artifacts_sha256"][key]!=previous["artifacts_sha256"][key]:raise ValueError("Reuse requires identical implementation, build and executable")
def sha(f):
    h=hashlib.sha256()
    with open(f,'rb') as stream:
        for b in iter(lambda:stream.read(1048576),b''):h.update(b)
    return h.hexdigest()
def validate(rows,build):
    names=[r['name'] for r in rows]
    if len(names)!=len(REQUIRED) or set(names)!=REQUIRED:raise ValueError('Missing, duplicate or unexpected required check')
    for r in rows:
        if r['exit_code']!=r.get('expected_exit',0):raise ValueError('Failed required command: '+r['name'])
        if r['name']=='m14-reader-rejects-schema5':
            if r.get('passed') is not False or not r.get('reason'):raise ValueError('Legacy rejection not established')
        elif r['name']!='runner-self-tests' and (r.get('passed') is not True or r.get('build')!=build):raise ValueError('Failed/mismatched result '+r['name'])
def main():
    p=argparse.ArgumentParser();p.add_argument('--candidate',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--pro-tools-original',type=pathlib.Path);p.add_argument('--reuse-output',type=pathlib.Path);p.add_argument('--reuse-candidate',type=pathlib.Path);a=p.parse_args();root=a.candidate.resolve();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    identity=json.loads((root/'identity.json').read_text());assert platform.system()=='Darwin' and platform.machine()=='arm64'
    for name,h in identity['artifacts_sha256'].items():assert sha(root/name)==h,name
    reused_rows={};previous=None
    if a.reuse_output or a.reuse_candidate:
        assert a.reuse_output and a.reuse_candidate
        previous=json.loads((a.reuse_candidate/"identity.json").read_text());old=json.loads((a.reuse_output/"summary.json").read_text());assert old["short_automated_passed"] and old["source_digest"]==previous["source_digest"]
        reused_rows={r["name"]:r for r in old["commands"]}
    rows=[]
    def run(name,args,report=None,expected=0):
        if previous and name in REUSABLE:
            validate_reuse(identity,previous,REUSABLE[name]);row=dict(reused_rows[name]);old_report=a.reuse_output/pathlib.Path(report).name
            assert sha(old_report)==row['report_sha256'];shutil.copy2(old_report,report);shutil.copy2(a.reuse_output/(name+'.log'),out/(name+'.log'))
            row['reused_from']={'source_digest':previous['source_digest'],'output':str(a.reuse_output.resolve()),'reason':'Identical production source/build and this exact executable; only marker-test coverage/reporting and qualification tooling changed'}
            rows.append(row);(out/'commands.json').write_text(json.dumps(rows,indent=2)+'\n');return
        started=time.monotonic();command=[str(x) for x in args]
        with (out/(name+'.log')).open('w') as log:r=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT)
        row=dict(name=name,command=command,exit_code=r.returncode,expected_exit=expected,seconds=time.monotonic()-started)
        if report and pathlib.Path(report).is_file():
            d=json.loads(pathlib.Path(report).read_text());row.update(passed=d.get('passed'),build=d.get('build',d.get('source_build')),reason=d.get('message',d.get('failure',d.get('error',''))),report_sha256=sha(report))
        if name=="markers" and report and pathlib.Path(report).is_file():validate_marker_groups(d.get("completed_groups",[]))
        rows.append(row);(out/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
        if r.returncode!=expected:raise RuntimeError('Required command failed: '+name)
    def binary(name):return ['/usr/bin/arch','-arm64',root/'bin'/name]
    started=time.monotonic();failure=''
    try:
        run('markers',binary('LevelScopeMarkerTests')+[out/'marker-files',root/'source/Tests/Markers/nuendo-15.0.30-points.csv']+([a.pro_tools_original.resolve()] if a.pro_tools_original else []),out/'marker-files/marker-tests.json')
        run('sessions-short',binary('LevelScopeSessionTests')+['--short','--output',out/'sessions.json'],out/'sessions.json')
        run('sections',binary('LevelScopeSectionsTests')+['--report',out/'sections.json'],out/'sections.json')
        run('m12-smoke',binary('LevelScopeSectionProcessingTests')+['--smoke','--root',out/'m12-files','--output',out/'m12.json'],out/'m12.json')
        run('m13-smoke',binary('LevelScopeSectionProposalTests')+['--smoke','--root',out/'m13-files','--output',out/'m13.json'],out/'m13.json')
        run('export-short',binary('LevelScopeExportTests')+['--artifacts',out/'export-files','--output',out/'export.json'],out/'export.json')
        run('m14-reader-rejects-schema5',['/usr/bin/arch','-arm64',root/'legacy/M14ProgrammeTests','--only','--m8-session',out/'marker-files/transaction-2/markers.levelscope','--output',out/'legacy-rejection.json'],out/'legacy-rejection.json',1)
        run('runner-self-tests',['python3','-B',root/'Tests/M15QualificationRunnerTests.py'])
        validate(rows,identity['measurement_build'])
    except Exception as e:failure=str(e)
    result=dict(short_automated_passed=not failure,failure=failure,source_digest=identity['source_digest'],measurement_build=identity['measurement_build'],
      runtime_os=subprocess.check_output(['sw_vers'],text=True).strip(),cpu=subprocess.check_output(['sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),architecture=platform.machine(),seconds=time.monotonic()-started,commands=rows,
      real_daw_interoperability='Nuendo 15.0.30 actual CSV and programme-zero origin tested at user-confirmed 25 fps; exact settings in fixture record. Pro Tools 10.3.10 actual table tested with explicit 24 fps LOCATION and matching-rate sample references; programme-origin confirmation recorded separately',silent_UI='separate required receipt',high_sierra='deferred after M18, no later than M19',video='deferred; no loading/decoding/playback',three_hour='deferred; synthetic section-population/index tests are not full-duration audio',physical_devices_opened=False,qualification_complete=False)
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='commands'},indent=2));return 1 if failure else 0
if __name__=='__main__':raise SystemExit(main())
