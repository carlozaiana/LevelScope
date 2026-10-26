#!/usr/bin/env python3
"""Freeze only explicit M11 targets/assets; never invoke historical packagers."""
import argparse,hashlib,json,pathlib,re,shutil,subprocess
from collect_m1_environment import sha256

TARGETS=('LevelScopeCompositionTests','LevelScopeMeasurementTests','LevelScopeScaffoldTests',
 'LevelScopeComplianceTests','LevelScopeProposalTests','LevelScopeDiagnosticsTests','LevelScopeSectionsTests',
 'LevelScopeAuditionTests','LevelScopeMultichannelTests','LevelScopeSessionTests','LevelScopeProgrammeTests','LevelScopeExportTests','LevelScopeVideoTests')
SCRIPTS=('run_m11_qualification.py','M11QualificationRunnerTests.py','verify_m11_pcm.py','compare_m11_reports.py',
 'verify_m10_pcm.py','verify_m10_decoder.py','collect_m1_environment.py','run_ebu_measurements.py','run_reference_measurements.py','EbuReferenceRunnerTests.py')

def main():
    p=argparse.ArgumentParser();p.add_argument('--build',type=pathlib.Path,required=True);p.add_argument('--m10',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    def copy(source,relative):
        destination=a.output/relative;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination)
    for target in TARGETS:
        source=a.build/(target+'_artefacts/Release')/target
        if set(subprocess.check_output(['/usr/bin/lipo','-archs',str(source)],universal_newlines=True).split())!={'arm64','x86_64'}:raise ValueError('Universal binary required: '+target)
        copy(source,pathlib.Path('bin')/target)
    app=a.build/'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app';shutil.copytree(app,a.output/app.name,symlinks=True,ignore=shutil.ignore_patterns('.DS_Store'))
    old=json.loads((a.m10/'identity.json').read_text())
    if old['source_digest']!='4ff89967e4927600d0b6defafabb28c372f4068870a0f1a9adb1ee84fe1781d8':raise ValueError('Pinned M10 candidate required')
    assets=['video-fixtures/cfr-exact.mov','video-fixtures/vfr-exact.mov','m7-video-fixtures/422-with-tone.mov']
    assets+=['m9-fixtures/m8-sessions/'+name+'.levelscope' for name in ('stereo','5.1','7.1')]
    for suite,directory in (('itu-bs2217','itu'),('ebu-v05','ebu')):
        manifest=pathlib.Path('Tests/Reference')/(suite+'-manifest.json');data=json.loads(manifest.read_text());copy(manifest,manifest)
        assets+=['vectors/'+directory+'/'+row['file'] for row in data['files']]
        if directory=='ebu':assets.append('vectors/ebu/readme.txt')
        assets.append('baseline/'+directory+'.json')
    for name in assets:
        if sha256(a.m10/name)!=old['artifacts_sha256'][name]:raise ValueError('Pinned asset changed: '+name)
        copy(a.m10/name,pathlib.Path(name))
    legacy='bin/LevelScopeProgrammeTests'
    if sha256(a.m10/legacy)!=old['artifacts_sha256'][legacy]:raise ValueError('M10 reader changed')
    copy(a.m10/legacy,'legacy/M10ProgrammeTests')
    for name in SCRIPTS:copy(pathlib.Path('Tests')/name,pathlib.Path('Tests')/name)
    paths=[pathlib.Path('CMakeLists.txt')]+[x for x in pathlib.Path('Source').rglob('*') if x.is_file() and x.name!='.DS_Store']
    paths += [x for x in pathlib.Path('Tests').glob('*') if x.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths += [x for x in pathlib.Path('Tests/VideoPaint').glob('*') if x.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths += list(pathlib.Path('Tests/Reference').glob('*-manifest.json'))
    source={str(x):sha256(x) for x in sorted(paths)}
    for x in paths:copy(x,pathlib.Path('source')/x)
    flags=(a.build/'CMakeFiles/LevelScopeStandaloneSupport.dir/flags.make').read_text();build=re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([0-9a-f]{64})',flags).group(1)
    identity=dict(milestone='M11 aligned M&E and dialogue composition',base_commit=subprocess.check_output(['git','rev-parse','HEAD'],universal_newlines=True).strip(),
      branch=subprocess.check_output(['git','branch','--show-current'],universal_newlines=True).strip(),source_state='Unpublished working candidate; source manifest, not base commit, identifies implementation',
      source_files_sha256=source,source_digest=hashlib.sha256(json.dumps(source,sort_keys=True,separators=(',',':')).encode()).hexdigest(),measurement_build=build,
      build_configuration=[x for x in (a.build/'CMakeCache.txt').read_text().splitlines() if x.startswith(('CMAKE_BUILD_TYPE:','CMAKE_OSX_ARCHITECTURES:','CMAKE_OSX_DEPLOYMENT_TARGET:','CMAKE_CXX_FLAGS:','CMAKE_CXX_FLAGS_RELEASE:','FETCHCONTENT_SOURCE_DIR_JUCE:'))],
      compiler=subprocess.check_output(['/usr/bin/clang','--version'],universal_newlines=True).strip(),juce_version='7.0.12',juce_commit='4f43011b96eb0636104cb3e433894cda98243626',
      asset_provenance=dict(m10_source_digest=old['source_digest'],m10_identity_sha256=sha256(a.m10/'identity.json'),legacy_reader_sha256=sha256(a.m10/legacy)),
      artifacts_sha256={str(x.relative_to(a.output)):sha256(x) for x in sorted(a.output.rglob('*')) if x.is_file()})
    (a.output/'identity.json').write_text(json.dumps(identity,indent=2)+'\n');print(identity['source_digest'])
if __name__=='__main__':main()
