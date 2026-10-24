#!/usr/bin/env python3
"""Create a new, exactly identified diagnostic package; never modify preserved candidates."""
import argparse
import hashlib
import json
import pathlib
import shutil
import subprocess
from run_comparison import sha256, write_new

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['build','output','control-package','fixture-directory','actual']:
        parser.add_argument('--'+name,type=pathlib.Path,required=True)
    args=parser.parse_args(); root=pathlib.Path(__file__).resolve().parents[2]
    build=args.build.resolve(); output=args.output.resolve(); control=args.control_package.resolve()
    original=json.loads((control/'identity.json').read_text())
    original_exe='LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
    expected='c8efc1a6d0b760206c8f99105e553542a4d32fdb9722dc3e19b2172f60bb36e9'
    if sha256(control/original_exe)!=expected:raise RuntimeError('Preserved M6 harness identity mismatch')
    if sha256(args.actual)!='0870d5ab93d1492f23524ed1fff9319af1a6a9cf71882543f2190929afeb8e24':
        raise RuntimeError('Actual reference movie identity mismatch')
    output.mkdir(parents=True,exist_ok=False)
    def copy(source,relative):
        target=output/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
    app='LevelScope Video Paint Diagnostic.app'
    shutil.copytree(build/'LevelScopeVideoPaintDiagnostic_artefacts/Release'/app,output/app,symlinks=True,ignore=shutil.ignore_patterns('.DS_Store'))
    shutil.copytree(control/'LevelScope Video Qualification.app',output/'control/LevelScope Video Qualification.app',symlinks=True,ignore=shutil.ignore_patterns('.DS_Store'))
    copy(control/'identity.json','control/original-m6-identity.json')
    for target in ['LevelScopeVideoPaintTraceTests','LevelScopeVideoTests','LevelScopeMeasurementTests','LevelScopeAuditionTests']:
        source=build/target if target=='LevelScopeVideoPaintTraceTests' else build/(target+'_artefacts/Release')/target
        copy(source,'bin/'+target)
    for rate in [44100,48000]:
        name='selected-clock-%d.wav'%rate;copy(args.fixture_directory/name,'fixtures/'+name)
    copy(args.actual,'fixtures/51_SyncTest.mov')
    paths=[root/'CMakeLists.txt']
    paths += [p for p in (root/'Source').rglob('*') if p.is_file() and p.name!='.DS_Store']
    paths += [p for p in (root/'Tests').glob('*') if p.suffix in ['.cpp','.h','.mm']]
    paths += [p for p in (root/'Tests/VideoPaint').glob('*') if p.suffix in ['.cpp','.h','.mm','.py']]
    source_hashes={p.relative_to(root).as_posix():sha256(p) for p in sorted(paths)}
    for p in paths:copy(p,pathlib.Path('Source')/p.relative_to(root))
    for p in (root/'Tests/VideoPaint').glob('*.py'):copy(p,'Tests/'+p.name)
    copy(build/'CMakeCache.txt','build-metadata/CMakeCache.txt')
    for target in ['LevelScopeVideoPaintDiagnostic','LevelScopeStandaloneSupport']:
        copy(build/'CMakeFiles'/(target+'.dir')/'flags.make','build-metadata/'+target+'-flags.make')
    executables={'control':'control/'+original_exe,'diagnostic':app+'/Contents/MacOS/LevelScope Video Paint Diagnostic'}
    identity={
        'scope':'Diagnostic comparison only; M7 remains frozen and no qualification is claimed',
        'base_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,universal_newlines=True).strip(),
        'source_state':'Uncommitted diagnostic candidate identified by exact source hashes',
        'source_files_sha256':source_hashes,
        'source_digest':hashlib.sha256(json.dumps(source_hashes,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
        'juce_commit':'4f43011b96eb0636104cb3e433894cda98243626','juce_version':'7.0.12',
        'compiler':subprocess.check_output(['/usr/bin/clang','--version'],universal_newlines=True).strip(),
        'frozen_control':{'source_digest':original['source_digest'],'identity_sha256':sha256(control/'identity.json'),'harness_sha256':expected},
        'executables':executables,'seconds_per_run':100,'fixed_order':['control','diagnostic']*3,
        'playback_requires_fresh_coordination':True,
        'artifacts_sha256':{p.relative_to(output).as_posix():sha256(p) for p in sorted(output.rglob('*')) if p.is_file()}}
    write_new(output/'identity.json',identity)
    print(json.dumps({'package':str(output),'source_digest':identity['source_digest'],
        'identity_sha256':sha256(output/'identity.json'),'diagnostic_sha256':sha256(output/executables['diagnostic']),
        'artifact_count':len(identity['artifacts_sha256'])},indent=2))
if __name__=='__main__':main()
