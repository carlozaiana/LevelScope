"""Freeze a local M20 Logic integration package; no publication or installation."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def sha(file):
    value = hashlib.sha256()
    with file.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('compile', type=Path)
    parser.add_argument('destination', type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[3]
    dest = args.destination.resolve()
    dest.mkdir(parents=True, exist_ok=False)
    names = subprocess.check_output(['git', 'ls-files', '-z', '--cached', '--others',
                                     '--exclude-standard', 'Source', 'Tests', 'CMakeLists.txt'], cwd=repo).decode().split('\0')
    sources = {}
    for name in sorted(set(names)):
        f = Path(name)
        if not name or not (repo / f).is_file() or 'results' in f.parts or f.name == '.DS_Store':
            continue
        if f.suffix not in ('.cpp', '.h', '.mm', '.py', '.md', '.txt'):
            continue
        to = dest / 'source' / f
        to.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(repo / f, to)
        sources[name] = sha(to)
    compile_dir = args.compile.resolve()
    shutil.copytree(compile_dir / 'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app',
                    dest / 'LevelScope Standalone.app', symlinks=True)
    (dest / 'bin').mkdir()
    for name in ('LevelScopeRegionPlanningTests', 'LevelScopeLraProposalTests',
                 'LevelScopeJobRequirementsTests', 'LevelScopeSectionProposalTests'):
        shutil.copy2(compile_dir / (name + '_artefacts/Release/' + name), dest / 'bin' / name)
    (dest / 'configuration').mkdir()
    shutil.copy2(compile_dir / 'CMakeCache.txt', dest / 'configuration/CMakeCache.txt')
    for target in ('LevelScopeStandalone', 'LevelScopeStandaloneSupport', 'LevelScopeRegionPlanningTests'):
        for name in ('flags.make', 'link.txt'):
            src = compile_dir / ('CMakeFiles/' + target + '.dir') / name
            shutil.copy2(src, dest / 'configuration' / (target + '-' + name))
    flags = (compile_dir / 'CMakeFiles/LevelScopeStandalone.dir/flags.make').read_text()
    build = re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([a-f0-9]{64})', flags)[1]
    identity = dict(kind='M20 Logic integration package, not combined qualification',
                    source_base_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
                    source_files_sha256=sources,
                    source_manifest_sha256=hashlib.sha256(json.dumps(sources, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
                    measurement_build=build, architectures=['arm64', 'x86_64'], intel_minimum='10.13',
                    compiler=subprocess.check_output(['/usr/bin/c++', '--version'], text=True),
                    os=subprocess.check_output(['/usr/bin/sw_vers'], text=True),
                    combined_qualification_complete=False,
                    artifacts_sha256={str(f.relative_to(dest)): sha(f) for f in sorted(dest.rglob('*')) if f.is_file() and 'source' not in f.relative_to(dest).parts})
    (dest / 'identity.json').write_text(json.dumps(identity, indent=2) + '\n')
    for file in dest.rglob('*'):
        if file.is_file():
            file.chmod(file.stat().st_mode & ~0o222)
    print(json.dumps({k: identity[k] for k in ('measurement_build', 'source_manifest_sha256')}))


if __name__ == '__main__':
    main()
