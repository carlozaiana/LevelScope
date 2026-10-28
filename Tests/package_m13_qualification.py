#!/usr/bin/env python3
"""Freeze named M12/M13 standalone/test targets. Never overwrite a package."""
import argparse,hashlib,json,pathlib,re,shutil,subprocess
from collect_m1_environment import sha256
from package_m12_qualification import audio_only_session
TARGETS=('LevelScopeSectionProposalTests','LevelScopeSectionProcessingTests','LevelScopeMeasurementTests','LevelScopeScaffoldTests','LevelScopeComplianceTests','LevelScopeProposalTests','LevelScopeDiagnosticsTests','LevelScopeSectionsTests','LevelScopeAuditionTests','LevelScopeMultichannelTests','LevelScopeSessionTests','LevelScopeProgrammeTests','LevelScopeExportTests','LevelScopeCompositionTests','LevelScopeVideoTests')
def main():
    p=argparse.ArgumentParser();p.add_argument('--build',type=pathlib.Path,required=True);p.add_argument('--m11',type=pathlib.Path,required=True);p.add_argument('--m12',type=pathlib.Path,required=True);p.add_argument('--physical-video',type=pathlib.Path,required=True);p.add_argument('--fractional-video',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
    assert subprocess.check_output(['git','rev-parse','main'],text=True).strip()=='cfb1033aa7fb65a2655db276bee696161196bdeb'
    a.output.mkdir(parents=True,exist_ok=False)
    def copy(src,relative):
        dest=a.output/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    for target in TARGETS:
        src=a.build/(target+'_artefacts/Release')/target
        assert set(subprocess.check_output(['/usr/bin/lipo','-archs',str(src)],text=True).split())=={'arm64','x86_64'}
        copy(src,pathlib.Path('bin')/target)
    for target,name in [('LevelScopeStandalone','LevelScope Standalone.app'),('LevelScopeVideoQualification','LevelScope Video Qualification.app')]:
        src=a.build/(target+'_artefacts/Release')/name;shutil.copytree(src,a.output/name,symlinks=True,ignore=shutil.ignore_patterns('.DS_Store'))
    baseline=json.loads((a.m11/'identity.json').read_text());assert baseline['source_digest']=='669b350652dbc3fb09d0ff724758208fec21bba40e0f6bc0188b89e6404c2c61'
    assets={}
    for name,digest in baseline['artifacts_sha256'].items():
        if name.startswith(('vectors/','baseline/')):
            path=a.m11/name;assert sha256(path)==digest;copy(path,name);assets[name]=dict(original=str(path.resolve()),sha256=digest)
    for name in ('stereo','5.1','7.1'):
        key='m9-fixtures/m8-sessions/'+name+'.levelscope';path=a.m11/key;assert sha256(path)==baseline['artifacts_sha256'][key]
        dest=a.output/'audio-only-m8-sessions'/path.name;dest.parent.mkdir(exist_ok=True);dest.write_bytes(audio_only_session(path.read_bytes()))
        assets[str(dest.relative_to(a.output))]=dict(original=str(path.resolve()),original_sha256=sha256(path),sha256=sha256(dest),change='Historical attachment removed only; not video migration evidence')
    for name in ('h264-cfr.mov','h264-vfr.mov','422-25-1.mov'):copy(a.physical_video/name,pathlib.Path('video-fixtures')/name)
    copy(a.fractional_video,'video-fixtures/lt-30000-1001.mov')
    for old,name in ((a.m11,'M11ProgrammeTests'),(a.m12,'M12ProgrammeTests')):
        old_id=json.loads((old/'identity.json').read_text());src=old/'bin/LevelScopeProgrammeTests';assert sha256(src)==old_id['artifacts_sha256']['bin/LevelScopeProgrammeTests'];copy(src,'legacy/'+name)
    for path in pathlib.Path('Tests').glob('*.py'):copy(path,pathlib.Path('Tests')/path.name)
    for path in pathlib.Path('Tests/Reference').glob('*-manifest.json'):copy(path,path)
    for path in pathlib.Path('Tests/VideoPaint').glob('*.py'):copy(path,path)
    paths=[pathlib.Path('CMakeLists.txt')]+[x for x in pathlib.Path('Source').rglob('*') if x.is_file() and x.name!='.DS_Store']
    paths+=[x for x in pathlib.Path('Tests').glob('*') if x.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths+=[x for x in pathlib.Path('Tests/VideoPaint').glob('*') if x.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths+=list(pathlib.Path('Tests/Reference').glob('*-manifest.json'))
    source={str(x):sha256(x) for x in sorted(paths)}
    for path in paths:copy(path,pathlib.Path('source')/path)
    flags=(a.build/'CMakeFiles/LevelScopeStandaloneSupport.dir/flags.make').read_text();build=re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([0-9a-f]{64})',flags).group(1)
    identity=dict(milestone='M13 bounded section-led proposals with M12 processing',base_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),branch=subprocess.check_output(['git','branch','--show-current'],text=True).strip(),
        source_state='Unpublished candidate; source manifest identifies implementation, not base commit',qualification_status='Pending actual-platform qualification',source_files_sha256=source,source_digest=hashlib.sha256(json.dumps(source,sort_keys=True,separators=(',',':')).encode()).hexdigest(),measurement_build=build,
        build_configuration=[x for x in (a.build/'CMakeCache.txt').read_text().splitlines() if x.startswith(('CMAKE_BUILD_TYPE:','CMAKE_OSX_ARCHITECTURES:','CMAKE_OSX_DEPLOYMENT_TARGET:','CMAKE_CXX_FLAGS:','CMAKE_CXX_FLAGS_RELEASE:','FETCHCONTENT_SOURCE_DIR_JUCE:'))],compiler=subprocess.check_output(['/usr/bin/clang','--version'],text=True).strip(),juce_commit='4f43011b96eb0636104cb3e433894cda98243626',historical_assets=assets,
        artifacts_sha256={str(x.relative_to(a.output)):sha256(x) for x in sorted(a.output.rglob('*')) if x.is_file()})
    (a.output/'identity.json').write_text(json.dumps(identity,indent=2)+'\n');print(identity['source_digest'])
if __name__=='__main__':main()
