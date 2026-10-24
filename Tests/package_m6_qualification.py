#!/usr/bin/env python3
"""Package an exact M6 candidate without installing or requiring a Git commit.

Source hashes identify the candidate; base_commit is ancestry, not a claim that
uncommitted changes are in that commit. Output includes audio/binaries: ignore it.
"""
import argparse
import hashlib
import json
import pathlib
import shutil
import subprocess

from collect_m1_environment import sha256
from package_m1_qualification import git

VIDEO_FIXTURES = ('cfr-exact.mov', 'vfr-exact.mov', 'two-streams-audio.mov',
                  'timecode.mov', 'cfr-exact.mp4', 'cfr-exact.m4v', 'half-hour.mov')
TARGETS = ('LevelScopeMeasurementTests', 'LevelScopeScaffoldTests',
           'LevelScopeComplianceTests', 'LevelScopeVideoTests', 'LevelScopeProposalTests', 'LevelScopeDiagnosticsTests', 'LevelScopeSectionsTests', 'LevelScopeAuditionTests')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=pathlib.Path)
    parser.add_argument('--build', required=True, type=pathlib.Path)
    parser.add_argument('--m1-bundle', required=True, type=pathlib.Path)
    parser.add_argument('--video-fixtures', required=True, type=pathlib.Path)
    parser.add_argument('--ui-fixtures', required=True, type=pathlib.Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)

    def copy(source, relative):
        destination = args.output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    app = args.build / 'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app'
    shutil.copytree(app, args.output / 'LevelScope Standalone.app', symlinks=True,
                    ignore=shutil.ignore_patterns('.DS_Store'))
    video_app = args.build / 'LevelScopeVideoQualification_artefacts/Release/LevelScope Video Qualification.app'
    shutil.copytree(video_app, args.output / video_app.name, symlinks=True, ignore=shutil.ignore_patterns('.DS_Store'))
    for target in TARGETS:
        copy(args.build / (target + '_artefacts/Release') / target, pathlib.Path('bin') / target)
    for name in ('run_prores_native_qualification.py', 'run_prores_playback_qualification.py', 'run_reference_measurements.py', 'run_ebu_measurements.py',
                 'EbuReferenceRunnerTests.py', 'collect_m1_environment.py', 'run_m6_qualification.py',
                 'compare_m2_reports.py', 'compare_m3_reports.py', 'M3QualificationRunnerTests.py', 'compare_m4_reports.py', 'M4QualificationRunnerTests.py', 'compare_m5_reports.py', 'M5QualificationRunnerTests.py', 'compare_m6_reports.py', 'M6QualificationRunnerTests.py'):
        copy(pathlib.Path('Tests') / name, pathlib.Path('Tests') / name)
    for suite, directory in [('itu-bs2217', 'itu'), ('ebu-v05', 'ebu')]:
        manifest_path = pathlib.Path('Tests/Reference') / (suite + '-manifest.json')
        manifest = json.loads(manifest_path.read_text())
        copy(manifest_path, manifest_path)
        files = {case['file']: case['sha256'] for case in manifest['files']}
        if directory == 'ebu':
            files['readme.txt'] = manifest['readme_sha256']
        for name, expected in files.items():
            relative = pathlib.Path('vectors') / directory / name
            source = args.m1_bundle / relative
            if sha256(source) != expected:
                raise SystemExit('Pinned fixture mismatch: ' + str(source))
            copy(source, relative)
        copy(args.m1_bundle / 'baseline' / (directory + '.json'), 'baseline/' + directory + '.json')
    for name in VIDEO_FIXTURES:
        copy(args.video_fixtures / name, pathlib.Path('video-fixtures') / name)
    for name in ('stereo.wav', 'stereo.aiff', 'stereo.flac', 'long-stereo.wav', 'unsupported-12ch.wav', 'trim-required.wav', 'gate-crossing.wav', 'diagnostic-scenes.wav', 'diagnostic-silence.wav', 'diagnostic-long-silence.wav', 'diagnostic-short.wav', 'diagnostic-three-hour.wav'):
        copy(args.ui_fixtures / name, pathlib.Path('ui-fixtures') / name)

    # Capture all implementation inputs, including untracked new source/tests.
    paths = [pathlib.Path('CMakeLists.txt')]
    paths += [p for p in pathlib.Path('Source').rglob('*') if p.is_file() and p.name != '.DS_Store']
    paths += [p for p in pathlib.Path('Tests').glob('*') if p.suffix in ('.cpp', '.h', '.py', '.sh')]
    paths += list(pathlib.Path('Tests/Reference').glob('*-manifest.json'))
    source_hashes = {str(p): sha256(p) for p in sorted(paths)}
    for path in paths:
        copy(path, pathlib.Path('source') / path)
    source_digest = hashlib.sha256(json.dumps(source_hashes, sort_keys=True,
                                             separators=(',', ':')).encode()).hexdigest()
    build_metadata = {}
    for line in (args.build / 'CMakeCache.txt').read_text().splitlines():
        if line.startswith(('CMAKE_BUILD_TYPE:', 'CMAKE_OSX_ARCHITECTURES:', 'CMAKE_OSX_DEPLOYMENT_TARGET:',
                            'CMAKE_OSX_SYSROOT:', 'CMAKE_CXX_FLAGS:', 'CMAKE_CXX_FLAGS_RELEASE:', 'FETCHCONTENT_SOURCE_DIR_JUCE:')):
            key, value = line.split('=', 1)
            build_metadata[key] = value
    identity = dict(base_commit=git('rev-parse', 'HEAD'), branch=git('branch', '--show-current'),
                    source_state='Working candidate; exact files below, not base commit alone',
                    source_files_sha256=source_hashes, source_digest=source_digest,
                    build_configuration=build_metadata,
                    compiler=subprocess.check_output(['/usr/bin/clang', '--version'], universal_newlines=True).strip(),
                    juce_version='7.0.12', juce_commit='4f43011b96eb0636104cb3e433894cda98243626',
                    artifacts_sha256={str(p.relative_to(args.output)): sha256(p)
                                      for p in sorted(args.output.rglob('*')) if p.is_file()})
    (args.output / 'identity.json').write_text(json.dumps(identity, indent=2) + '\n')
    print('Source digest:', source_digest)
    print('identity.json SHA256:', sha256(args.output / 'identity.json'))


if __name__ == '__main__':
    main()
