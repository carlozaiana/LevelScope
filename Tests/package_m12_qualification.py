#!/usr/bin/env python3
"""Freeze explicitly named M12 audio-only targets; preserve all older packages."""
import argparse, hashlib, json, pathlib, re, shutil, struct, subprocess
from collect_m1_environment import sha256

TARGETS=('LevelScopeSectionProcessingTests','LevelScopeMeasurementTests','LevelScopeScaffoldTests',
 'LevelScopeComplianceTests','LevelScopeProposalTests','LevelScopeDiagnosticsTests','LevelScopeSectionsTests',
 'LevelScopeAuditionTests','LevelScopeMultichannelTests','LevelScopeSessionTests','LevelScopeProgrammeTests',
 'LevelScopeExportTests','LevelScopeCompositionTests')
SCRIPTS=('package_m12_qualification.py','run_m12_qualification.py','complete_m12_sequoia.py','M12QualificationRunnerTests.py','verify_m12_pcm.py','verify_m11_pcm.py',
 'verify_m10_pcm.py','verify_m10_decoder.py','collect_m1_environment.py','run_ebu_measurements.py',
 'run_reference_measurements.py','EbuReferenceRunnerTests.py')

def audio_only_session(original):
    """Replace only attachment values, retaining the reader's canonical JSON."""
    length=struct.unpack_from('<Q',original,8)[0]
    encoded=original[80:80+length]
    pattern=rb'("video": \{"locator": )"(?:\\.|[^"\\])*"(, "sha256": )"[0-9a-f]*"(, "stream": )[0-9]+'
    encoded,count=re.subn(pattern,lambda m:m[1]+b'""'+m[2]+b'""'+m[3]+b'0',encoded)
    if count!=1:raise ValueError('Expected one canonical historical video descriptor')
    return original[:8]+struct.pack('<Q',len(encoded))+hashlib.sha256(encoded).hexdigest().encode()+encoded+original[80+length:]

def main():
    p=argparse.ArgumentParser();p.add_argument('--build',required=True,type=pathlib.Path);p.add_argument('--m11',required=True,type=pathlib.Path);p.add_argument('--output',required=True,type=pathlib.Path);a=p.parse_args()
    if subprocess.check_output(['git','rev-parse','main'],text=True).strip()!='cfb1033aa7fb65a2655db276bee696161196bdeb':raise ValueError('Main must remain at M11')
    a.output.mkdir(parents=True,exist_ok=False)
    def copy(src,relative):
        dest=a.output/relative;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dest)
    for target in TARGETS:
        src=a.build/(target+'_artefacts/Release')/target
        if set(subprocess.check_output(['/usr/bin/lipo','-archs',str(src)],text=True).split())!={'arm64','x86_64'}:raise ValueError('Universal build required: '+target)
        copy(src,pathlib.Path('bin')/target)
    app=a.build/'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app'
    shutil.copytree(app,a.output/app.name,symlinks=True,ignore=shutil.ignore_patterns('.DS_Store'))
    baseline=json.loads((a.m11/'identity.json').read_text())
    if baseline['source_digest']!='669b350652dbc3fb09d0ff724758208fec21bba40e0f6bc0188b89e6404c2c61':raise ValueError('Pinned M11 package required')
    assets={}
    for suite,directory in (('itu-bs2217','itu'),('ebu-v05','ebu')):
        manifest=pathlib.Path('Tests/Reference')/(suite+'-manifest.json');copy(manifest,manifest)
        names=['vectors/'+directory+'/'+row['file'] for row in json.loads(manifest.read_text())['files']]
        if directory=='ebu':names.append('vectors/ebu/readme.txt')
        names.append('baseline/'+directory+'.json')
        for name in names:
            path=a.m11/name;expected=baseline['artifacts_sha256'][name]
            if sha256(path)!=expected:raise ValueError('Historical asset changed: '+name)
            assets[name]=dict(path=str(path.resolve()),sha256=expected)
    migration_inputs=[]
    for name in ('stereo','5.1','7.1'):
        key='m9-fixtures/m8-sessions/'+name+'.levelscope';path=a.m11/key;expected=baseline['artifacts_sha256'][key]
        if sha256(path)!=expected:raise ValueError('Historical session changed')
        assets[key]=dict(path=str(path.resolve()),sha256=expected)
        # A historical stereo session has an attached video. Derive a clearly
        # identified audio-only test document; never open that attachment in M12.
        original=path.read_bytes();length=struct.unpack_from('<Q',original,8)[0]
        core=json.loads(original[80:80+length]);old_video=core['video'].copy()
        derived=audio_only_session(original)
        relative=pathlib.Path('audio-only-m8-sessions')/(name+'.levelscope');dest=a.output/relative;dest.parent.mkdir(exist_ok=True);dest.write_bytes(derived)
        migration_inputs.append(dict(original_sha256=expected,derived_path=str(relative),derived_sha256=sha256(dest),original_video=old_video,change='Remove video attachment only; preserve rational offset. Not video migration evidence.'))
    legacy=a.m11/'bin/LevelScopeProgrammeTests';copy(legacy,'legacy/M11ProgrammeTests')
    if sha256(legacy)!=baseline['artifacts_sha256']['bin/LevelScopeProgrammeTests']:raise ValueError('M11 reader changed')
    for name in SCRIPTS:copy(pathlib.Path('Tests')/name,pathlib.Path('Tests')/name)
    paths=[pathlib.Path('CMakeLists.txt')]+[x for x in pathlib.Path('Source').rglob('*') if x.is_file() and x.name!='.DS_Store']
    paths += [x for x in pathlib.Path('Tests').glob('*') if x.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths += [x for x in pathlib.Path('Tests/VideoPaint').glob('*') if x.suffix in ('.cpp','.h','.mm','.py','.sh')]
    paths += list(pathlib.Path('Tests/Reference').glob('*-manifest.json'))
    source={str(x):sha256(x) for x in sorted(paths)}
    for x in paths:copy(x,pathlib.Path('source')/x)
    flags=(a.build/'CMakeFiles/LevelScopeStandaloneSupport.dir/flags.make').read_text();build=re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([0-9a-f]{64})',flags).group(1)
    identity=dict(milestone='M12 section gains — Sequoia phase',base_commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
      branch=subprocess.check_output(['git','branch','--show-current'],text=True).strip(),source_state='Unpublished candidate; manifest identifies implementation, not base commit',
      qualification_status='Not yet Sequoia tested; High Sierra and video qualification pending',source_files_sha256=source,
      source_digest=hashlib.sha256(json.dumps(source,sort_keys=True,separators=(',',':')).encode()).hexdigest(),measurement_build=build,
      build_configuration=[x for x in (a.build/'CMakeCache.txt').read_text().splitlines() if x.startswith(('CMAKE_BUILD_TYPE:','CMAKE_OSX_ARCHITECTURES:','CMAKE_OSX_DEPLOYMENT_TARGET:','CMAKE_CXX_FLAGS:','CMAKE_CXX_FLAGS_RELEASE:','FETCHCONTENT_SOURCE_DIR_JUCE:'))],
      compiler=subprocess.check_output(['/usr/bin/clang','--version'],text=True).strip(),juce_commit='4f43011b96eb0636104cb3e433894cda98243626',
      referenced_assets=assets,audio_only_migration_inputs=migration_inputs,m11_identity_sha256=sha256(a.m11/'identity.json'),
      artifacts_sha256={str(x.relative_to(a.output)):sha256(x) for x in sorted(a.output.rglob('*')) if x.is_file()})
    (a.output/'identity.json').write_text(json.dumps(identity,indent=2)+'\n');print(identity['source_digest'])
if __name__=='__main__':main()
