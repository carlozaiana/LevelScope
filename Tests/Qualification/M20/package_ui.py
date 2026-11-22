"""Freeze a new combined M20 UI candidate; never publish, install or overwrite."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def sha(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('compile', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--targets', nargs='+', required=True)
    parser.add_argument('--logic-handover', default='294f9b52abc1ae21754d5cf42c879b9de674244e')
    parser.add_argument('--published-baseline', default='22961f03de7919e1de631f8b13eb4b852108ac0b')
    parser.add_argument('--kind', default='M20 combined UI review candidate')
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    assert re.fullmatch(r'[a-f0-9]{40}', args.logic_handover), 'Exact Logic handover commit required'
    subprocess.run(['git', 'merge-base', '--is-ancestor', args.logic_handover, 'HEAD'], cwd=repo, check=True)
    build, destination = args.compile.resolve(), args.destination.resolve()
    assert not destination.exists(), 'Preserved packages are immutable'
    assert not subprocess.check_output(['git', 'diff', 'HEAD', '--', 'Source', 'Tests', 'CMakeLists.txt'], cwd=repo), 'Commit tested inputs before freezing'
    names = subprocess.check_output(['git', 'ls-files', '-z', 'Source', 'Tests', 'CMakeLists.txt'], cwd=repo).decode().split('\0')
    destination.mkdir(parents=True)
    for name in names:
        if not name or 'results' in Path(name).parts or Path(name).name == '.DS_Store':
            continue
        source = repo / name
        if source.is_file():
            target = destination / 'source' / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    shutil.copytree(build / 'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app', destination / 'LevelScope Standalone.app', symlinks=True)
    (destination / 'bin').mkdir()
    (destination / 'configuration').mkdir()
    shutil.copy2(build / 'CMakeCache.txt', destination / 'configuration/CMakeCache.txt')
    identity_flags = (build / 'CMakeFiles/LevelScopeStandalone.dir/flags.make').read_text()
    identity = re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([a-f0-9]{64})', identity_flags)[1]
    for target in ['LevelScopeStandalone', 'LevelScopeStandaloneSupport', *args.targets]:
        assert identity in (build / ('CMakeFiles/' + target + '.dir/flags.make')).read_text(), target
        for name in ('flags.make', 'link.txt'):
            source = build / ('CMakeFiles/' + target + '.dir/' + name)
            if source.exists():
                shutil.copy2(source, destination / 'configuration' / (target + '-' + name))
    for target in args.targets:
        shutil.copy2(build / (target + '_artefacts/Release/' + target), destination / 'bin' / target)
    binaries = [destination / 'LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone', *sorted((destination / 'bin').iterdir())]
    binary_checks = {}
    for binary in binaries:
        architectures = subprocess.check_output(['lipo', '-archs', str(binary)], text=True).split()
        commands = subprocess.check_output(['otool', '-arch', 'x86_64', '-l', str(binary)], text=True)
        minimum = re.search(r'cmd LC_VERSION_MIN_MACOSX\s+cmdsize \d+\s+version ([\d.]+)', commands)
        assert set(architectures) == {'arm64', 'x86_64'} and minimum and minimum[1] == '10.13'
        binary_checks[str(binary.relative_to(destination))] = {'architectures': architectures, 'intel_minimum': minimum[1]}
    sources = {str(p.relative_to(destination / 'source')): sha(p) for p in sorted((destination / 'source').rglob('*')) if p.is_file()}
    result = {'kind': args.kind, 'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
              'published_baseline': args.published_baseline, 'logic_handover': args.logic_handover,
              'measurement_build': identity, 'source_sha256': sources,
              'source_manifest_sha256': hashlib.sha256(json.dumps(sources, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
              'binary_checks': binary_checks, 'artifacts_sha256': {str(p.relative_to(destination)): sha(p) for p in sorted(destination.rglob('*')) if p.is_file()},
              'combined_qualification_complete': False, 'publication_authorized': False}
    (destination / 'identity.json').write_text(json.dumps(result, indent=2) + '\n')
    (destination / 'identity.sha256').write_text(sha(destination / 'identity.json') + '  identity.json\n')
    for p in destination.rglob('*'):
        if p.is_file():
            p.chmod(p.stat().st_mode & ~0o222)
    print(json.dumps({k: result[k] for k in ('source_commit', 'measurement_build', 'source_manifest_sha256')}, indent=2))


if __name__ == '__main__':
    main()
