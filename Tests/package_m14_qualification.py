#!/usr/bin/env python3
"""Freeze M14 Sequoia candidate and universal binaries; no devices or video assets."""
import argparse,hashlib,json,pathlib,re,shutil,subprocess
from collect_m1_environment import sha256
BASE='a9be070ab9eb5cb4f591ff9d991d8999c3a016d0'
TARGETS=('LevelScopeMonoTests','LevelScopeSectionProcessingTests','LevelScopeSectionProposalTests',
 'LevelScopeMeasurementTests','LevelScopeScaffoldTests','LevelScopeComplianceTests','LevelScopeProposalTests',
 'LevelScopeDiagnosticsTests','LevelScopeSectionsTests','LevelScopeAuditionTests','LevelScopeMultichannelTests',
 'LevelScopeSessionTests','LevelScopeProgrammeTests','LevelScopeExportTests','LevelScopeCompositionTests')
def main():
    p=argparse.ArgumentParser();p.add_argument('--build',type=pathlib.Path,required=True);p.add_argument('--m13',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
    assert subprocess.check_output(['git','rev-parse','main'],text=True).strip()==BASE
    assert subprocess.check_output(['git','branch','--show-current'],text=True).strip()=='codex/logic-module'
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
    prior=json.loads((a.m13/'identity.json').read_text());assert prior['source_digest']=='50eba19a2d66096607ec878c8eaa05e610580c7f5e8854d7b375d846a89c5f27'
    reused={}
    for name,h in prior['artifacts_sha256'].items():
        if name.startswith(('vectors/','baseline/','audio-only-m8-sessions/')):
            f=a.m13/name;assert sha256(f)==h;copy(f,name);reused[name]={'original_sha256':h,'original_candidate':prior['source_digest']}
    old=a.m13/'bin/LevelScopeProgrammeTests';assert sha256(old)==prior['artifacts_sha256']['bin/LevelScopeProgrammeTests'];copy(old,'legacy/M13ProgrammeTests')
    for f in pathlib.Path('Tests').glob('*.py'):copy(f,pathlib.Path('Tests')/f.name)
    for f in pathlib.Path('Tests/Reference').glob('*-manifest.json'):copy(f,f)
    paths=[pathlib.Path('CMakeLists.txt')]+[f for f in pathlib.Path('Source').rglob('*') if f.is_file() and f.name!='.DS_Store']
    paths += [f for f in pathlib.Path('Tests').glob('*') if f.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths += [f for f in pathlib.Path('Tests/VideoPaint').glob('*') if f.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths += list(pathlib.Path('Tests/Reference').glob('*-manifest.json'))
    source={str(f):sha256(f) for f in sorted(paths)}
    for f in paths:copy(f,pathlib.Path('source')/f)
    flags=(a.build/'CMakeFiles/LevelScopeStandaloneSupport.dir/flags.make').read_text();build=re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([0-9a-f]{64})',flags).group(1)
    identity=dict(milestone='M14 mono processing and bounded RF64 input',base_commit=BASE,branch='codex/logic-module',
      source_state='Unpublished source manifest; base commit is not the candidate implementation',
      qualification_status='Pending Sequoia phase; High Sierra, video and three-hour mono/RF64 qualification deferred',
      source_files_sha256=source,source_digest=hashlib.sha256(json.dumps(source,sort_keys=True,separators=(',',':')).encode()).hexdigest(),measurement_build=build,
      build_configuration=[s for s in (a.build/'CMakeCache.txt').read_text().splitlines() if s.startswith(('CMAKE_BUILD_TYPE:','CMAKE_OSX_ARCHITECTURES:','CMAKE_OSX_DEPLOYMENT_TARGET:','CMAKE_CXX_FLAGS:','CMAKE_CXX_FLAGS_RELEASE:','FETCHCONTENT_SOURCE_DIR_JUCE:'))],
      compiler=subprocess.check_output(['/usr/bin/clang','--version'],text=True).strip(),juce_commit='4f43011b96eb0636104cb3e433894cda98243626',historical_assets=reused,
      artifacts_sha256={str(f.relative_to(a.output)):sha256(f) for f in sorted(a.output.rglob('*')) if f.is_file()})
    (a.output/'identity.json').write_text(json.dumps(identity,indent=2)+'\n');print(identity['source_digest'])
if __name__=='__main__':main()
