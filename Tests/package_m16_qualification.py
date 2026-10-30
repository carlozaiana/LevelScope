#!/usr/bin/env python3
"""Freeze the explicit M16 Sequoia candidate; never rebuild/test/install implicitly."""
import argparse,hashlib,json,pathlib,re,shutil,subprocess
MAIN='a9be070ab9eb5cb4f591ff9d991d8999c3a016d0'
BASE='63eafac9acc210d762dfd5e3dcc70defbe5bbdf3'
PRIOR='38676a45d47ddba6eb54f7ab299f52f622f6be26554256781b8960749c609564'
TARGETS=('LevelScopeSectionProposalTests','LevelScopeProposalTests','LevelScopeSectionsTests','LevelScopeSessionTests','LevelScopeSectionProcessingTests','LevelScopeMarkerTests','LevelScopeExportTests')
def sha(f):
    h=hashlib.sha256()
    with open(f,'rb') as s:
        for b in iter(lambda:s.read(1048576),b''):h.update(b)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--build',type=pathlib.Path,required=True);p.add_argument('--m15',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
    for ref,want in [('main',MAIN),('HEAD',BASE)]:assert subprocess.check_output(['git','rev-parse',ref],text=True).strip()==want
    assert subprocess.check_output(['git','branch','--show-current'],text=True).strip()=='codex/logic-module'
    prior=json.loads((a.m15/'identity.json').read_text());assert prior['source_digest']==PRIOR
    for f,h in prior['artifacts_sha256'].items():assert sha(a.m15/f)==h,('M15 preservation',f)
    a.output.mkdir(parents=True,exist_ok=False)
    def copy(src,name):
        dst=a.output/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
    for target in TARGETS:
        src=a.build/(target+'_artefacts/Release')/target
        assert set(subprocess.check_output(['/usr/bin/lipo','-archs',str(src)],text=True).split())=={'arm64','x86_64'}
        copy(src,pathlib.Path('bin')/target)
    app=a.build/'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app'
    assert set(subprocess.check_output(['/usr/bin/lipo','-archs',str(app/'Contents/MacOS/LevelScope Standalone')],text=True).split())=={'arm64','x86_64'}
    shutil.copytree(app,a.output/app.name,symlinks=True,ignore=shutil.ignore_patterns('.DS_Store'))
    copy(a.m15/'bin/LevelScopeSectionProposalTests','legacy/M15SectionProposalTests')
    paths=[pathlib.Path('CMakeLists.txt')]+[f for f in pathlib.Path('Source').rglob('*') if f.is_file() and f.name!='.DS_Store']
    paths += [f for f in pathlib.Path('Tests').glob('*') if f.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths += [f for f in pathlib.Path('Tests/Markers').rglob('*') if f.suffix in ('.csv','.txt')]
    source={str(f):sha(f) for f in sorted(paths)}
    for f in paths:copy(f,pathlib.Path('source')/f)
    for f in ('run_m16_qualification.py','M16QualificationRunnerTests.py'):copy(pathlib.Path('Tests')/f,pathlib.Path('Tests')/f)
    flags=(a.build/'CMakeFiles/LevelScopeStandaloneSupport.dir/flags.make').read_text();build=re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([0-9a-f]{64})',flags).group(1)
    changed=[f for f,h in source.items() if f.startswith('Source/') and prior['source_files_sha256'].get(f)!=h]
    identity=dict(milestone='M16 configurable proposal preferences and limits',base_commit=BASE,main_preserved=MAIN,branch='codex/logic-module',
      source_state='Unpublished exact source manifest; base commit does not identify this implementation',qualification_status='Sequoia candidate checks pending; actual High Sierra, video and three-hour checks deferred after M18, no later than M19',
      source_files_sha256=source,source_digest=hashlib.sha256(json.dumps(source,sort_keys=True,separators=(',',':')).encode()).hexdigest(),measurement_build=build,
      build_configuration=[s for s in (a.build/'CMakeCache.txt').read_text().splitlines() if s.startswith(('CMAKE_BUILD_TYPE:','CMAKE_OSX_ARCHITECTURES:','CMAKE_OSX_DEPLOYMENT_TARGET:','CMAKE_CXX_FLAGS:','CMAKE_CXX_FLAGS_RELEASE:','FETCHCONTENT_SOURCE_DIR_JUCE:'))],
      compiler=subprocess.check_output(['/usr/bin/clang','--version'],text=True).strip(),juce_commit='4f43011b96eb0636104cb3e433894cda98243626',
      previous_candidate=PRIOR,changed_production_files=changed,unchanged_m15_files={f:h for f,h in source.items() if prior['source_files_sha256'].get(f)==h},
      artifacts_sha256={str(f.relative_to(a.output)):sha(f) for f in sorted(a.output.rglob('*')) if f.is_file()})
    (a.output/'identity.json').write_text(json.dumps(identity,indent=2)+'\n');print(identity['source_digest'])
if __name__=='__main__':main()
