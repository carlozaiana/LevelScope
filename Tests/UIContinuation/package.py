"""Preserve an exact focused UI candidate without overwriting tested builds."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('build', type=Path)
    parser.add_argument('destination', type=Path)
    parser.add_argument('--slice', required=True)
    parser.add_argument('--targets', nargs='+', required=True)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[2]
    build, dest = args.build.resolve(), args.destination.resolve()
    assert not dest.exists(), 'A preserved candidate must never be overwritten'
    assert not subprocess.check_output(['git', 'diff', 'HEAD', '--', 'Source', 'CMakeLists.txt', 'Tests/UIContinuation'], cwd=repo), 'Commit candidate inputs before freezing'
    dest.mkdir(parents=True)
    inputs = subprocess.check_output(['git', 'ls-files', '-z', 'Source', 'Tests', 'CMakeLists.txt'], cwd=repo).decode().split('\0')
    for name in inputs:
        if not name or 'results' in Path(name).parts or Path(name).name == '.DS_Store':
            continue
        source = repo / name
        if source.is_file():
            target = dest / 'source' / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
    shutil.copytree(build / 'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app', dest / 'LevelScope Standalone.app', symlinks=True)
    (dest / 'bin').mkdir()
    (dest / 'configuration').mkdir()
    shutil.copy2(build / 'CMakeCache.txt', dest / 'configuration/CMakeCache.txt')
    for target in ['LevelScopeStandalone', 'LevelScopeStandaloneSupport', *args.targets]:
        for name in ['flags.make', 'link.txt']:
            source = build / 'CMakeFiles' / (target + '.dir') / name
            if source.exists():
                shutil.copy2(source, dest / 'configuration' / (target + '-' + name))
    for target in args.targets:
        shutil.copy2(build / (target + '_artefacts/Release') / target, dest / 'bin' / target)
    binaries = [dest / 'LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone', *sorted((dest / 'bin').iterdir())]
    binary_checks = {}
    for binary in binaries:
        architectures = subprocess.check_output(['lipo', '-archs', str(binary)], text=True).split()
        commands = subprocess.check_output(['otool', '-arch', 'x86_64', '-l', str(binary)], text=True)
        minimum = re.search(r'cmd LC_VERSION_MIN_MACOSX\s+cmdsize \d+\s+version ([\d.]+)', commands)
        if not minimum:
            minimum = re.search(r'cmd LC_BUILD_VERSION\s+cmdsize \d+\s+platform \d+\s+minos ([\d.]+)', commands)
        assert set(architectures) == {'arm64', 'x86_64'} and minimum and minimum[1] == '10.13'
        binary_checks[str(binary.relative_to(dest))] = {'architectures': architectures, 'intel_minimum': minimum[1]}
    flags = (build / 'CMakeFiles/LevelScopeStandalone.dir/flags.make').read_text()
    measurement_build = re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([a-f0-9]{64})', flags)[1]
    for target in args.targets:
        assert measurement_build in (build / 'CMakeFiles' / (target + '.dir') / 'flags.make').read_text()
    identity = {
        'schema': 1, 'slice': args.slice,
        'source_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
        'baseline': '5fc8c6ee95e705077cf996cecbc6543dbc20d257',
        'measurement_build': measurement_build, 'binary_checks': binary_checks,
        'source_sha256': {str(p.relative_to(dest / 'source')): sha(p) for p in sorted((dest / 'source').rglob('*')) if p.is_file()},
        'artifacts_sha256': {str(p.relative_to(dest)): sha(p) for p in sorted(dest.rglob('*')) if p.is_file()},
        'qualification_scope': 'Focused UI/connection tests; reuse unchanged M19/Intel qualification. Physical opening requires separate coordination.',
    }
    (dest / 'identity.json').write_text(json.dumps(identity, indent=2) + '\n')
    (dest / 'identity.sha256').write_text(sha(dest / 'identity.json') + '  identity.json\n')
    for path in dest.rglob('*'):
        if path.is_file():
            path.chmod(path.stat().st_mode & ~0o222)
    print(json.dumps({k: v for k, v in identity.items() if k not in ('source_sha256', 'artifacts_sha256')}, indent=2))


if __name__ == '__main__':
    main()
