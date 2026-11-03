"""Focused measured M/S, timeline, transport and Target development selection; no devices or video."""
import hashlib, json, pathlib, re, shutil, subprocess, sys, time
# This script lives under Tests/Reference/results/<slice>/.
repo = pathlib.Path(__file__).resolve().parents[4]
build = repo / 'build/m18-dev-01'
out = pathlib.Path(sys.argv[1]).resolve()
out.mkdir(parents=True, exist_ok=False)
fixture = pathlib.Path('/tmp') / ('levelscope-m18-progress-navigation-' + out.name)
assert not fixture.exists(), fixture
fixture.mkdir()
def binary(target):
    return build / (target + '_artefacts/Release') / target
selection = [
 ('timeline',binary('LevelScopeTimelineProgressTests'),['--root',str(fixture/'timeline'),'--output',str(out/'timeline.json')]),
 ('measurement',binary('LevelScopeMeasurementTests'),['--report',str(out/'measurement.json')]),
 ('diagnostics',binary('LevelScopeDiagnosticsTests'),['--report',str(out/'diagnostics.json')]),
 ('waveforms',binary('LevelScopeWaveformTests'),['--root',str(fixture/'waveforms'),'--output',str(out/'waveforms.json')]),
 ('targets',binary('LevelScopeJobRequirementsTests'),['--root',str(fixture/'targets'),'--output',str(out/'targets.json')]),
]
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()
flags=(build/'CMakeFiles/LevelScopeStandaloneSupport.dir/flags.make').read_text()
identity=re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([0-9a-f]{64})',flags).group(1)
source={str(p.relative_to(repo)):sha(p) for p in sorted((repo/'Source').rglob('*')) if p.suffix in ('.h','.cpp','.mm')}
inputs={str(p.relative_to(repo)):sha(p) for p in sorted((repo/'Tests').glob('*')) if p.suffix in ('.h','.cpp')}
inputs['CMakeLists.txt']=sha(repo/'CMakeLists.txt')
for extra in [pathlib.Path(__file__).resolve(),repo/'Tests/run_ebu_measurements.py',repo/'Tests/Reference/ebu-v05-manifest.json']:
    inputs[str(extra.relative_to(repo))]=sha(extra)
(out/'source-manifest.json').write_text(json.dumps(source,indent=2)+'\n')
(out/'validation-manifest.json').write_text(json.dumps({'measurement_build':identity,'inputs':inputs,'juce_commit':'4f43011b96eb0636104cb3e433894cda98243626','build_configuration':[s for s in (build/'CMakeCache.txt').read_text().splitlines() if s.startswith(('CMAKE_BUILD_TYPE:','CMAKE_OSX_ARCHITECTURES:','CMAKE_OSX_DEPLOYMENT_TARGET:','CMAKE_CXX_FLAGS:','CMAKE_CXX_FLAGS_RELEASE:','FETCHCONTENT_SOURCE_DIR_JUCE:'))]},indent=2)+'\n')
executables=[item[1] for item in selection]+[build/'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone']
binaries=[]
for exe in executables:
    archs=subprocess.check_output(['/usr/bin/lipo','-archs',str(exe)],text=True).strip()
    loads=subprocess.check_output(['/usr/bin/otool','-arch','x86_64','-l',str(exe)],text=True)
    minimum=re.search(r'cmd LC_VERSION_MIN_MACOSX\s+cmdsize \d+\s+version (\S+)',loads)
    if not minimum: minimum=re.search(r'cmd LC_BUILD_VERSION\s+cmdsize \d+\s+platform \S+\s+minos (\S+)',loads)
    assert set(archs.split())=={'arm64','x86_64'} and minimum and minimum.group(1)=='10.13', (exe,archs,minimum)
    binaries.append({'path':str(exe.relative_to(repo)),'sha256':sha(exe),'architectures':archs,'intel_minimum':minimum.group(1)})
(out/'binaries.json').write_text(json.dumps(binaries,indent=2)+'\n')
rows=[]
for name,exe,args in selection:
    print('Running '+name,flush=True); started=time.monotonic()
    with (out/(name+'.log')).open('w') as log: result=subprocess.run([str(exe)]+args,cwd=repo,stdout=log,stderr=subprocess.STDOUT)
    passed=result.returncode==0
    if passed:
        report=json.loads((out/(name+'.json')).read_text()) if (out/(name+'.json')).exists() else {}
        passed=passed and report.get('passed',False) and report.get('build')==identity
    row={'name':name,'command':[str(exe)]+args,'exit_code':result.returncode,'passed':passed,'seconds':time.monotonic()-started}
    rows.append(row);(out/'summary.json').write_text(json.dumps({'measurement_build':identity,'all_passed':False,'runs':rows},indent=2)+'\n')
    print(name+': '+('PASS' if passed else 'FAIL'),flush=True)
# Run only the four pinned EBU cases whose acceptance rules exercise measured
# window curves, using the existing independent reference runner/checks.
sys.path.insert(0,str(repo/'Tests'))
from run_ebu_measurements import run_case
vectors=pathlib.Path('/Users/fuzz/Developer/LevelScope/Tests/Ebu_Loudness_Testset_V05/ebu-loudness-test-setv05')
manifest=json.loads((repo/'Tests/Reference/ebu-v05-manifest.json').read_text())
curves=[c for c in manifest['files'] if any(rule['kind']=='curve' for rule in c['checks'])]
assert len(curves)==4
print('Running pinned EBU M/S curve cases',flush=True);started=time.monotonic()
references=[run_case(case,binary('LevelScopeMeasurementTests'),vectors,'arm64') for case in curves]
ref_pass=all(item['passed'] and item['measured']['build_identity']==identity for item in references)
(out/'reference-curves.json').write_text(json.dumps({'passed':ref_pass,'results':references},indent=2)+'\n')
rows.append({'name':'pinned EBU M/S curves','passed':ref_pass,'seconds':time.monotonic()-started,'files':[c['file'] for c in curves]})
for name in ('partial.png','navigation.png','compact.png'):
    shutil.copy2(fixture/'timeline'/name,out/name)
assert all(sha(repo/p)==digest for p,digest in {**source,**inputs}.items()),'inputs changed during run'
(out/'summary.json').write_text(json.dumps({'measurement_build':identity,'all_passed':len(rows)==6 and all(r['passed'] for r in rows),'runs':rows,'app_path':str(executables[-1].parents[2]),'physical_devices_opened':False},indent=2)+'\n')
print('Reference curves: '+('PASS' if ref_pass else 'FAIL'),flush=True)
sys.exit(0 if all(r['passed'] for r in rows) else 1)
