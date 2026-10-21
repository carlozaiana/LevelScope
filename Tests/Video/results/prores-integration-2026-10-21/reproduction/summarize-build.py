from pathlib import Path
import json,hashlib,subprocess,xml.etree.ElementTree as ET
R=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation');B=R/'build/prores-integration-20261021';C=B/'build'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
xml=ET.parse(B/'ctest.xml').getroot();tests=list(xml.iter('testcase'))
assert len(tests)==7 and not list(xml.iter('failure')) and not list(xml.iter('error'))
targets=['LevelScopeStandalone','LevelScopeVideoTests','LevelScopeVideoQualification','LevelScopeMeasurementTests','LevelScopeScaffoldTests','LevelScopeComplianceTests','LevelScopeProposalTests','LevelScopeDiagnosticsTests','LevelScopeSectionsTests']
products={'LevelScopeStandalone':'LevelScope Standalone','LevelScopeVideoQualification':'LevelScope Video Qualification'}
binaries=[]
for t in targets:
 p=C/(t+'_artefacts')/'Release'
 p=p/(products[t]+'.app')/'Contents/MacOS'/products[t] if t in products else p/t
 assert p.is_file(),p
 arch=subprocess.check_output(['lipo','-archs',str(p)],text=True).strip()
 assert set(arch.split())=={'arm64','x86_64'}
 binaries.append(dict(target=t,path=str(p),sha256=sha(p),architectures=arch))
plugins=list(C.rglob('*.vst3'))
plugin_binaries=[str(f) for p in plugins for f in (p/'Contents/MacOS').glob('*') if f.is_file()]
assert not plugin_binaries
report=dict(passed=True,source_commit='6662091f1dfaff948a4a3eb3f000f47e37cfdc21',integration_revision=subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),source_changed=False,configuration='Release; universal x86_64/arm64; CMAKE_OSX_DEPLOYMENT_TARGET=10.13; all seven standalone test options ON',dependency='/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12',dependency_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd='/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12',text=True).strip(),targets=targets,executable_identities=binaries,ctest_count=len(tests),ctest_tests=[dict(t.attrib) for t in tests],ctest_architecture='arm64 on actual Sequoia',physical_audio_devices_opened=False,vst3_executables=plugin_binaries,vst3_configuration_placeholders=[str(p) for p in plugins],vst3_note='CMake generated an Info.plist and empty Contents/MacOS directory only; no plugin target, executable or installation',build_log=str(B/'build.log'),build_log_sha256=sha(B/'build.log'),configure_log_sha256=sha(B/'configure.log'),workflow_sha256=sha(R/'.github/workflows/build.yml'),scope_note='Fresh CI build and CTest regression only. Does not replace or inherit the frozen candidate04 two-platform playback/endurance qualification.')
(B/'build-summary.json').write_text(json.dumps(report,indent=2)+'\n');print('Fresh build: 9 universal targets, 7 CTest passes, no VST3 product.')
