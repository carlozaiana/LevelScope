import importlib.util,json,pathlib,subprocess,time,struct,hashlib
p=pathlib.Path(__file__).resolve().parent;c=p/'candidate-02';old=p/'candidate-01';out=p/'sequoia-candidate-02';out.mkdir(exist_ok=False)
s=importlib.util.spec_from_file_location('runner',c/'Tests/run_m17_qualification.py');m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
a=json.loads((old/'identity.json').read_text());b=json.loads((c/'identity.json').read_text())
for root,i in [(old,a),(c,b)]:
 for f,h in i['artifacts_sha256'].items():assert m.sha(root/f)==h,f
sd=[f for f in a['source_files_sha256'] if a['source_files_sha256'][f]!=b['source_files_sha256'].get(f)]
ad=[f for f in a['artifacts_sha256'] if a['artifacts_sha256'][f]!=b['artifacts_sha256'].get(f)]
assert sd==['Tests/MarkerTests.cpp'] and set(ad)=={'source/Tests/MarkerTests.cpp','bin/LevelScopeMarkerTests'}
assert a['measurement_build']==b['measurement_build']
prior=p/'sequoia-candidate-01';summary=json.loads((prior/'summary.json').read_text());rows=[]
for r in summary['commands']:
 if r['name']=='markers':continue
 assert r['exit_code']==0 and r['passed'] is True and r['build']==b['measurement_build']
 r.update(evidence_origin_candidate='candidate-01',evidence_origin_source_digest=a['source_digest'],reuse_reason='Only marker test assertion and marker test binary changed; application and this test executable are byte-identical.')
 rows.append(r)
assert len(rows)==8
reuse=dict(source_difference=sd,artifact_difference=ad,original_source=a['source_digest'],accepted_source=b['source_digest'],original_summary_sha256=m.sha(prior/'summary.json'),reason='Correct obsolete M16-reader marker test assertion to M17 reader. No production change; historical failed run retained.')
(out/'equivalence.json').write_text(json.dumps(reuse,indent=2)+'\n')
def run(name,args,report=None,expected=0):
 start=time.monotonic();cmd=[str(x) for x in args]
 with (out/(name+'.log')).open('w') as log:r=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
 row=dict(name=name,command=cmd,exit_code=r.returncode,expected_exit=expected,seconds=time.monotonic()-start,evidence_origin_candidate='candidate-02')
 if report and report.is_file():
  d=json.loads(report.read_text());row.update(passed=d.get('passed'),build=d.get('build',d.get('source_build')),reason=d.get('message',d.get('failure',d.get('error',''))),report_sha256=m.sha(report))
 if name=='m16-reader-rejects-schema7':
  f=prior/'m13-files/48000-2-0/applied.levelscope';raw=f.read_bytes();assert raw[:8]==b'LSSESS01';n=struct.unpack_from('<Q',raw,8)[0];core=raw[80:80+n];assert hashlib.sha256(core).hexdigest()==raw[16:80].decode();d=json.loads(core);assert d['schema']==7 and d['current']['sectionGains'] and d['caches'] and len(raw)>80+n
  row.update(legacy_input_schema=7,legacy_input_sha256=m.sha(f))
 rows.append(row);(out/'commands.json').write_text(json.dumps(rows,indent=2)+'\n')
 print(name,r.returncode,flush=True)
 if r.returncode!=expected:raise RuntimeError('Required command failed: '+name)
def binary(n):return ['/usr/bin/arch','-arm64',c/'bin'/n]
failure='';start=time.monotonic()
try:
 run('markers',binary('LevelScopeMarkerTests')+[out/'marker-files',c/'source/Tests/Markers/nuendo-15.0.30-points.csv'],out/'marker-files/marker-tests.json')
 run('export-short',binary('LevelScopeExportTests')+['--artifacts',out/'export-files','--output',out/'export.json'],out/'export.json')
 run('m16-reader-rejects-schema7',['/usr/bin/arch','-arm64',c/'legacy/M16SectionProposalTests','--exchange',prior/'m13-files','--root',out/'legacy-read-files','--output',out/'legacy-rejection.json'],out/'legacy-rejection.json',1)
 run('runner-self-tests',['python3','-B',c/'Tests/M17QualificationRunnerTests.py'])
 m.validate(rows,b['measurement_build'])
except Exception as e:failure=str(e)
result=dict(short_automated_passed=not failure,failure=failure,source_digest=b['source_digest'],measurement_build=b['measurement_build'],os=summary['os'],cpu=summary['cpu'],architecture=summary['architecture'],seconds_fresh=time.monotonic()-start,commands=rows,equivalence=reuse,silent_UI='separate required receipt',physical_devices_opened=False,video_loaded=False,high_sierra='deferred after M18, no later than M19',three_hour='deferred',qualification_complete=False)
(out/'summary.json').write_text(json.dumps(result,indent=2)+'\n');print(failure or 'All twelve required short commands validated; eight reused with exact binary equivalence.')
raise SystemExit(bool(failure))
