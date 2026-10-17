#!/usr/bin/env python3
"""Make an isolated, non-installing M1 test folder from a committed source tree.

Output contains third-party reference audio and binaries: keep it outside Git.
Transfer this folder to the older Mac only after inspecting its available runtimes.
"""
import argparse
import hashlib
import json
import pathlib
import shutil
import subprocess

from collect_m1_environment import sha256


def git(*args):
    return subprocess.check_output(['git'] + list(args), universal_newlines=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=pathlib.Path)
    parser.add_argument('--build', type=pathlib.Path, default=pathlib.Path('build/scaffold'))
    parser.add_argument('--itu', type=pathlib.Path, default=pathlib.Path('build/measurement-reference/itu'))
    parser.add_argument('--ebu', type=pathlib.Path, default=pathlib.Path('Tests/Ebu_Loudness_Testset_V05/ebu-loudness-test-setv05'))
    parser.add_argument('--fixtures', required=True, type=pathlib.Path)
    parser.add_argument('--ebu-baseline', required=True, type=pathlib.Path)
    parser.add_argument('--itu-baseline', required=True, type=pathlib.Path)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit('Output already exists; choose a fresh isolated directory')
    status = git('status', '--porcelain', '--untracked-files=all', '--', 'Source', 'Tests', 'CMakeLists.txt')
    changes = [line for line in status.splitlines() if not line.endswith('/.DS_Store')]
    if changes:
        raise SystemExit('Commit source/tests/evidence before packaging:\n' + '\n'.join(changes))
    args.output.mkdir(parents=True)

    def copy(source, relative):
        destination = args.output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    app = args.build / 'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app'
    shutil.copytree(app, args.output / 'LevelScope Standalone.app', symlinks=True,
                    ignore=shutil.ignore_patterns('.DS_Store'))
    for target in ('LevelScopeMeasurementTests', 'LevelScopeScaffoldTests'):
        copy(args.build / (target + '_artefacts/Release') / target, pathlib.Path('bin') / target)
    for name in ('run_reference_measurements.py', 'run_ebu_measurements.py', 'EbuReferenceRunnerTests.py', 'collect_m1_environment.py', 'run_m1_bundle.sh'):
        copy(pathlib.Path('Tests') / name, pathlib.Path('Tests') / name)
    for suite, vectors in [('itu-bs2217', args.itu), ('ebu-v05', args.ebu)]:
        manifest_path = pathlib.Path('Tests/Reference') / (suite + '-manifest.json')
        manifest = json.loads(manifest_path.read_text())
        copy(manifest_path, manifest_path)
        destination = pathlib.Path('vectors') / ('itu' if suite == 'itu-bs2217' else 'ebu')
        for case in manifest['files']:
            source = vectors / case['file']
            if sha256(source) != case['sha256']:
                raise SystemExit('Fixture mismatch: ' + str(source))
            copy(source, destination / case['file'])
        if suite == 'ebu-v05':
            if sha256(vectors / 'readme.txt') != manifest['readme_sha256']:
                raise SystemExit('EBU README mismatch')
            copy(vectors / 'readme.txt', destination / 'readme.txt')
    for name in ('stereo.wav', 'stereo.aiff', 'stereo.flac', 'long-stereo.wav', 'unsupported-12ch.wav'):
        source = args.fixtures / name
        if source.stat().st_size <= 44:
            raise SystemExit('Missing/incomplete UI fixture: ' + name)
        copy(source, pathlib.Path('ui-fixtures') / name)
    copy(args.ebu_baseline, 'baseline/ebu.json')
    copy(args.itu_baseline, 'baseline/itu.json')
    source_paths = subprocess.check_output(['git', 'ls-files', '-z', 'Source', 'Tests', 'CMakeLists.txt']).decode().split('\0')
    source_hashes = {p: sha256(pathlib.Path(p)) for p in source_paths if p and not p.endswith('.DS_Store')}
    build_metadata = {}
    for line in (args.build / 'CMakeCache.txt').read_text().splitlines():
        if line.startswith(('CMAKE_BUILD_TYPE:', 'CMAKE_OSX_ARCHITECTURES:', 'CMAKE_OSX_DEPLOYMENT_TARGET:', 'CMAKE_OSX_SYSROOT:', 'CMAKE_CXX_FLAGS:', 'CMAKE_CXX_FLAGS_RELEASE:')):
            key, value = line.split('=', 1)
            build_metadata[key] = value
    compiler = subprocess.check_output(['/usr/bin/clang', '--version'], universal_newlines=True).strip()
    artifact_hashes = {str(p.relative_to(args.output)): sha256(p)
                       for p in sorted(args.output.rglob('*')) if p.is_file()}
    identity = dict(source_commit=git('rev-parse', 'HEAD'), source_tree=git('rev-parse', 'HEAD^{tree}'),
                    source_files_sha256=source_hashes, build_configuration=build_metadata,
                    compiler=compiler, juce_version='7.0.12', artifacts_sha256=artifact_hashes)
    (args.output / 'identity.json').write_text(json.dumps(identity, indent=2) + '\n')
    print('Prepared', args.output, 'from', identity['source_commit'])
    print('identity.json SHA256:', sha256(args.output / 'identity.json'))


if __name__ == '__main__':
    main()
