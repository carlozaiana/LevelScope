#!/usr/bin/env python3
"""Freeze an M15 development candidate; no media/video/device/network operations."""
import argparse, hashlib, json, pathlib, re, shutil, subprocess
MAIN='a9be070ab9eb5cb4f591ff9d991d8999c3a016d0'
BASE='be4bf23c8bb9f661ae1ede30f3b6e79a067aea37'
TARGETS=('LevelScopeMarkerTests','LevelScopeSectionsTests','LevelScopeSessionTests',
         'LevelScopeSectionProcessingTests','LevelScopeSectionProposalTests','LevelScopeExportTests')
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    p=argparse.ArgumentParser();p.add_argument('--build',type=pathlib.Path,required=True);p.add_argument('--m14',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
    assert subprocess.check_output(['git','rev-parse','main'],text=True).strip()==MAIN
    assert subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip()==BASE
    assert subprocess.check_output(['git','branch','--show-current'],text=True).strip()=='codex/logic-module'
    prior=json.loads((a.m14/'identity.json').read_text());assert prior['source_digest']=='9dd22a598281e83431039f210202a8a24fe49517294cc6bb8d83032c4975e013'
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
    old=a.m14/'bin/LevelScopeProgrammeTests';assert sha(old)==prior['artifacts_sha256']['bin/LevelScopeProgrammeTests'];copy(old,'legacy/M14ProgrammeTests')
    paths=[pathlib.Path('CMakeLists.txt')]+[f for f in pathlib.Path('Source').rglob('*') if f.is_file() and f.name!='.DS_Store']
    paths += [f for f in pathlib.Path('Tests').glob('*') if f.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths += [f for f in pathlib.Path('Tests/Markers').rglob('*') if f.suffix in ('.csv','.txt')]
    source={str(f):sha(f) for f in sorted(paths)}
    for f in paths:copy(f,pathlib.Path('source')/f)
    for f in ('run_m15_qualification.py','M15QualificationRunnerTests.py'):copy(pathlib.Path('Tests')/f,pathlib.Path('Tests')/f)
    flags=(a.build/'CMakeFiles/LevelScopeStandaloneSupport.dir/flags.make').read_text();build=re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([0-9a-f]{64})',flags).group(1)
    unchanged={name:h for name,h in prior['source_files_sha256'].items() if name in source and source[name]==h}
    identity=dict(milestone='M15 point-marker import',base_commit=BASE,main_preserved=MAIN,branch='codex/logic-module',
      source_state='Unpublished source manifest; base commit is not this implementation',
      qualification_status='Pending Sequoia candidate checks; Nuendo and Pro Tools actual exports available; current-candidate interoperability checks pending; High Sierra, video and three-hour checks deferred',
      source_files_sha256=source,source_digest=hashlib.sha256(json.dumps(source,sort_keys=True,separators=(',',':')).encode()).hexdigest(),measurement_build=build,
      build_configuration=[s for s in (a.build/'CMakeCache.txt').read_text().splitlines() if s.startswith(('CMAKE_BUILD_TYPE:','CMAKE_OSX_ARCHITECTURES:','CMAKE_OSX_DEPLOYMENT_TARGET:','CMAKE_CXX_FLAGS:','CMAKE_CXX_FLAGS_RELEASE:','FETCHCONTENT_SOURCE_DIR_JUCE:'))],
      compiler=subprocess.check_output(['/usr/bin/clang','--version'],text=True).strip(),juce_commit='4f43011b96eb0636104cb3e433894cda98243626',
      reused_evidence_candidate=prior['source_digest'],unchanged_m14_files=unchanged,
      artifacts_sha256={str(f.relative_to(a.output)):sha(f) for f in sorted(a.output.rglob('*')) if f.is_file()})
    (a.output/'identity.json').write_text(json.dumps(identity,indent=2)+'\n');print(identity['source_digest'])
if __name__=='__main__':main()
