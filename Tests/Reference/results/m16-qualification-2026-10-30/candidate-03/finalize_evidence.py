import pathlib,json,hashlib,subprocess,shutil,re,datetime
root=pathlib.Path.cwd(); work=root/'build/m16-qualification-2026-10-30'; candidate=work/'candidate-03'; result=work/'sequoia-candidate-03'; evidence=root/'Tests/Reference/results/m16-qualification-2026-10-30/candidate-03'
evidence.mkdir(exist_ok=False)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(p.read_text())
def write(name,value):
 p=evidence/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(value,indent=2)+'\n')
def cmd(*a):return subprocess.check_output(a,text=True)
def cp(src,name):
 dst=evidence/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
def verify_package(p):
 d=read(p/'identity.json')
 for f,h in d['artifacts_sha256'].items():assert sha(p/f)==h,(p,f)
 return d
now=datetime.datetime.now(datetime.timezone.utc).isoformat(); ident=verify_package(candidate);first=verify_package(work/'candidate-02')
for f,h in ident['source_files_sha256'].items():assert sha(root/f)==h,f
assert first['measurement_build']!=ident['measurement_build']
source_diff=[f for f,h in ident['source_files_sha256'].items() if first['source_files_sha256'].get(f)!=h]
assert set(source_diff)=={'Source/Core/Proposal/Methods/SectionTrimMethod.cpp','Source/Standalone/StandaloneMainComponent.cpp','Tests/ProposalPreferencesTests.h','Tests/run_m16_qualification.py','Tests/M16QualificationRunnerTests.py'},source_diff
appkey='LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone'
write('candidate-03-change-impact.json',dict(from_candidate='candidate-03',to_candidate='candidate-03',from_source=first['source_digest'],to_source=ident['source_digest'],changed_source=source_diff,production_change='Display only: affected frames/seconds/percentage in accepted summaries; remove irrelevant empty global-draft notice when section proposal exists. No selection/sample arithmetic change.',test_change='Assert duration display; strengthen old-reader proof with exact expected build/failure and valid schema6 input checksum.',qualification='Full explicit short manifest rerun on candidate-03; silent UI repeated on candidate-03.'))
cp(candidate/'identity.json','candidate-03-identity.json');cp(work/'candidate-01/identity.json','candidate-01-identity.json')
summary=read(result/'summary.json');assert summary['short_automated_passed'] and not summary['failure'];assert summary['source_digest']==ident['source_digest']
import importlib.util
spec=importlib.util.spec_from_file_location('runner',candidate/'Tests/run_m16_qualification.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);mod.validate(summary['commands'],ident['measurement_build'])
for f in result.iterdir():
 if f.is_file() and f.suffix in ('.json','.log'):cp(f,pathlib.Path('sequoia-candidate-03')/f.name)
cp(result/'marker-files/marker-tests.json','sequoia-candidate-03/marker-tests.json')
cp(result/'preferences-files/large-measured-proposal/large-proposal-result.txt','sequoia-candidate-03/large-proposal-result.txt')
for f in ['summary.json','m13.json','m13-smoke.log']:
 cp(work/'sequoia-candidate-01'/f,pathlib.Path('historical-candidate-01')/f)
for f in sorted(work.glob('dev-policy-*.json')):cp(f,pathlib.Path('development')/f.name)
cp(work/'dev-large-01.json','development/dev-large-01.json')
prior_records=[]
for folder,name in [('m13-qualification-2026-10-28','candidate-02'),('m14-qualification-2026-10-29','candidate-02'),('m15-qualification-2026-10-29','candidate-05')]:
 p=root/'build'/folder/name;d=verify_package(p);prior_records.append(dict(path=str(p),source=d['source_digest'],identity_sha256=sha(p/'identity.json'),artifacts_verified=len(d['artifacts_sha256'])))
last=read(root/'build/m15-qualification-2026-10-29/candidate-05/identity.json')
unchanged={f:h for f,h in ident['source_files_sha256'].items() if f.startswith('Source/') and last['source_files_sha256'].get(f)==h}
protected=[f for f in ident['source_files_sha256'] if f.startswith(('Source/Core/Measurement/','Source/Core/Processing/')) or any(f.startswith('Source/Standalone/'+s) for s in ('Pcm24Wave','ProgrammeSource','ProgrammeComposition','StandaloneAudition','StandaloneMeasurement','StandaloneMedia','StandaloneVideo','VideoTiming','AudioPresentationClock','MeasurementTimeline')) or f.startswith(('Source/Core/Proposal/Model/TargetProfile','Source/Core/Proposal/VerificationEngine'))]
assert all(f in unchanged for f in protected)
write('reuse-audit.json',dict(previous_candidates=prior_records,unchanged_dependency_hashes={f:unchanged[f] for f in protected},changed_production=ident['changed_production_files'],cmake_change='Only adds StandaloneProposalPreferences.cpp to existing DiagnosticsUI; compiler flags, deployment target, dependencies and audio/video linkage unchanged.',compiler_and_dependencies='Apple clang / JUCE 7.0.12 revision recorded in candidate identity, same established universal Release configuration.',reused_evidence=[dict(source_manifest='50eba19a2d66096607ec878c8eaa05e610580c7f5e8854d7b375d846a89c5f27',report='Docs/ProjectIndex/49_M13QualificationAndHandover.md',scope='Historical pinned EBU/ITU and timing/video lower-level baseline. Not M16 runtime or new-feature qualification.'),dict(source_manifest='9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013',report='Docs/ProjectIndex/52_M14QualificationAndHandover.md',scope='Unchanged mono/RF64 arithmetic/decoder and prior bounded cinema / sustained 90-minute workflow; original binary-equivalent candidate-01 provenance retained.'),dict(source_manifest=last['source_digest'],report='Docs/ProjectIndex/55_M15QualificationAndHandover.md',scope='Actual Nuendo/Pro Tools export setup provenance and unchanged mapping; focused marker/session checks freshly rerun.')],fresh_shared_callers=[r['name'] for r in summary['commands']],deferred='Actual High Sierra, cross-Mac, video and three-hour testing after M18, no later than M19. No earlier-triggering arithmetic/clock/buffering/dependency change found.'))
arches=[]
for name in sorted(ident['artifacts_sha256']):
 if name.startswith('bin/') or name==appkey:
  p=candidate/name;arch=cmd('/usr/bin/lipo','-archs',str(p)).split();assert set(arch)=={'arm64','x86_64'}
  loads=cmd('/usr/bin/otool','-arch','x86_64','-l',str(p));match=re.search(r'LC_VERSION_MIN_MACOSX\s+cmdsize \d+\s+version ([\d.]+)',loads)
  assert match and match.group(1)=='10.13',(name,match.group(1) if match else 'missing')
  arches.append(dict(file=name,sha256=sha(p),architectures=arch,intel_minimum_macos=match.group(1)))
write('cross-build.json',dict(actual_runtime=summary['os'],cpu=summary['cpu'],executed_architecture='arm64',compiler=ident['compiler'],build_configuration=ident['build_configuration'],binaries=arches,actual_intel_high_sierra_executed=False,plugins_built_or_installed=False))
before=read(work/'before.json');head=cmd('git','rev-parse','HEAD');main=cmd('git','rev-parse','main');tags=cmd('git','for-each-ref','--format=%(refname) %(objectname) %(*objectname)','refs/tags');trees=cmd('git','worktree','list','--porcelain')
assert head==before['head'] and main==before['main'] and tags==before['tags'] and trees==before['worktrees']
for line in before['status'].splitlines():
 assert line.startswith('?? ') and (root/line[3:]).exists()
write('preservation.json',dict(checked_at=now,head=head.strip(),main=main.strip(),tag_count=len(tags.splitlines()),tag_objects_and_targets=tags,worktrees=trees,all_preexisting_tags_unchanged=True,prior_frozen_artifacts_verified=prior_records,preexisting_unrelated_files_preserved=before['status'].splitlines(),network_publication_performed=False,main_integration_performed=False))
write('ui-candidate-03.json',dict(date='2026-10-30',method='Native CUA accessibility and screenshot observation; no physical device opened, no video loaded.',actually_executed_app=str(work/'candidate-03/LevelScope Standalone.app'),app_sha256=ident['artifacts_sha256'][appkey],fresh_candidate03_checks=True,session=str(work/'sequoia-candidate-03/preferences-files/large-measured-proposal/large-proposal-input.levelscope'),checks=[
 dict(name='Safe reopen',passed=True,observed='EBU R128 profile/applicability restored; Source measurements current; Current gain0dB; stopped, device closed, both routes muted.'),
 dict(name='Large proposal generation',passed=True,observed='READY FOR EXPLICIT APPLY; 134 corrections / 1 connected region / 66720 of480000frames; displayed1.390000seconds and13.900000%; 3 private passes. No processing applied. Irrelevant no-global-draft notice is absent after completion.'),
 dict(name='Complete paginated inspection',passed=True,observed='Navigated pages1–4; final Section134 and complete final private measured PASS visible by scrolling. Current explicitly unchanged/unverified after future Apply; narrow selected-rule wording.'),
 dict(name='Draft cancellation',passed=True,observed='Entered100% then Cancel edits; active20% unchanged.'),
 dict(name='Invalid reserve',passed=True,observed='Additional0.2dB / reserve0.5dB rejected with explicit error; no clamp or audio change.'),
 dict(name='Explicit count above64',passed=True,observed='Confirmed65regions and100.0%; active digest91393bac03c8d8c16031433a31402305d1c67f5609a47b6ebdce30e9875626e4.'),
 dict(name='Stale rejection and preserved measurement',passed=True,observed='STALE PROPOSAL visible; Apply unavailable; Source measurements current and failing original IL still shown; Current gain0dB. Native automated identity test separately verifies both snapshot pointers/evaluations unchanged.'),
 dict(name='SaveAs and reopen',passed=True,observed='ui-candidate-03-preferences.levelscope saved, New session then reopened.65/100.0/3.0/12.0/0.5exact; no transient proposal restored; measurements current, deviceclosed/routesmuted. Save response completed in0.464s including UI overhead; reopen response0.578s for1000annotations, not maximum-session benchmark.'),
 dict(name='Explicit defaults',passed=True,observed='Draft defaults then Use restores unrestrictedcount/20%/3dB/12dB/0.5dB, digest9020f36cdf7584536fc2708e4f57196cbf3ca8ef1d72e954ae8ad2e88bb48e27. Audio unchanged.'),
],physical_devices_opened=False,video_loaded=False,passed=True))
resources=read(result/'preferences.json')['policy_resources'];testbytes=sum(p.stat().st_size for p in result.rglob('*') if p.is_file());assert testbytes<2*1024**3
ui_saved=work/'sequoia-candidate-03/preferences-files/large-measured-proposal/ui-candidate-03-preferences.levelscope';assert ui_saved.stat().st_size<128*1024**2
resources.update(test_files_bytes=testbytes,silent_ui_session_bytes=ui_saved.stat().st_size,silent_ui_session_sha256=sha(ui_saved),maximum_cache_session_not_rebenchmarked=True,scope='Short synthetic audio and10000-interval structural stress. Synthetic10800-second diagnostic coordinates are not a three-hour audio workflow.')
write('resources.json',resources)
write('status.json',dict(milestone='M16',candidate='candidate-03',source_manifest=ident['source_digest'],measurement_build=ident['measurement_build'],app_sha256=ident['artifacts_sha256'][appkey],base_commit=ident['base_commit'],implementation_commit='Not yet committed; exact implementation identified by frozen source manifest',status='Sequoia tested; High Sierra, video and three-hour checks pending.',sequoia_phase_complete=True,actual_sequoia_automated_passed=True,silent_ui_passed=True,universal_cross_build_passed=True,full_cross_platform_qualification_complete=False,ready_for_development_checkpoint_review=True,published=False,main=main.strip(),deferred=['Actual Intel High Sierra and cross-Mac exchange','Video loading/decoding/playback and selected synchronization tests','Three-hour/full-duration qualification','Combined qualification after M18, no later than M19; integration requires subsequent approval'],known_limitations=['Existing one-envelope section-led attenuation method; no coordinated neighbours or new strategies','Conservative connected support includes complete transition endpoints','Full private proposal explanations are transient; save measurement report to retain them; sessions preserve active/historical policy and exact applied processing','No custom LRA rules, progressive analysis, interactive editing or DSP introduced','Physical full7.1 listening separately pending'],originals_installations_and_audio_configuration_untouched=True))
cp(pathlib.Path(__file__),'finalize_evidence.py')
cp(evidence/'status.json','../status.json')
print(json.dumps(dict(source=ident['source_digest'],checks=len(summary['commands']),ui=True,test_bytes=testbytes,tags=len(tags.splitlines()),evidence=str(evidence)),indent=2))
