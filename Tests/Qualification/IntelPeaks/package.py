"""Freeze this focused candidate without touching any existing package."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


p = argparse.ArgumentParser()
p.add_argument('build', type=Path)
p.add_argument('destination', type=Path)
a = p.parse_args()
repo = Path(__file__).resolve().parents[3]
base, dest = a.build.resolve(), a.destination.resolve()
dest.mkdir(parents=True, exist_ok=False)
compile_dir = base / 'compile'
files = subprocess.check_output(['git', 'ls-files', '-z', '--cached', '--others', '--exclude-standard',
                                 'Source', 'Tests', 'CMakeLists.txt'], cwd=repo).decode().split('\0')
for name in sorted(set(files)):
    f = Path(name)
    if not name or not (repo / f).is_file() or 'results' in f.parts or f.name == '.DS_Store':
        continue
    if f.suffix not in ('.cpp', '.h', '.mm', '.py', '.md', '.txt'):
        continue
    to = dest / 'source' / f
    to.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(repo / f, to)
sources = {str(f.relative_to(dest / 'source')): sha(f)
           for f in sorted((dest / 'source').rglob('*')) if f.is_file()}
shutil.copytree(compile_dir / 'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app',
                dest / 'LevelScope Standalone.app', symlinks=True)
(dest / 'bin/peaks').mkdir(parents=True)
for name in ('LevelScopeMeasurementTests', 'LevelScopeSectionProposalTests'):
    shutil.copy2(compile_dir / (name + '_artefacts/Release/' + name), dest / 'bin' / name)
shutil.copy2(compile_dir / 'LevelScopePeakChannelParity', dest / 'bin/LevelScopePeakChannelParity')
shutil.copy2(base / 'drivers/native', dest / 'bin/native')
for name in ('serial', 'app-mode'):
    shutil.copy2(base / 'drivers/peaks' / name, dest / 'bin/peaks' / name)
shutil.copytree(base / 'fixtures', dest / 'fixtures')
(dest / 'configuration').mkdir()
shutil.copy2(compile_dir / 'CMakeCache.txt', dest / 'configuration/CMakeCache.txt')
for target in ('LevelScopeStandalone', 'LevelScopeStandaloneSupport', 'LevelScopeMeasurementCore'):
    for name in ('flags.make', 'link.txt'):
        f = compile_dir / ('CMakeFiles/' + target + '.dir') / name
        if f.exists():
            shutil.copy2(f, dest / 'configuration' / (target + '-' + name))
flags = (compile_dir / 'CMakeFiles/LevelScopeStandalone.dir/flags.make').read_text()
build = re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([a-f0-9]{64})', flags)[1]
artifacts = {str(f.relative_to(dest)): sha(f) for f in sorted(dest.rglob('*')) if f.is_file()}
identity = dict(schema=1, name='Intel bounded true peak ' + dest.name,
    source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=repo, text=True).strip(),
    source_manifest_sha256=hashlib.sha256(json.dumps(sources, sort_keys=True, separators=(',', ':')).encode()).hexdigest(),
    source_files_sha256=sources, measurement_build=build, artifacts_sha256=artifacts,
    architectures=['arm64', 'x86_64'], intel_minimum='10.13',
    serial_oracle='Private same-build service object with LEVELSCOPE_ENABLE_BOUNDED_TRUE_PEAK=0; not the product app',
    app_policy='Bounded per-channel stereo/5.1 on both macOS architectures; existing serial fallbacks',
    compiler=subprocess.check_output(['/usr/bin/c++', '--version'], text=True),
    juce_commit=subprocess.check_output(['git', '-C', '/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12', 'rev-parse', 'HEAD'], text=True).strip(),
    qualification_complete=False)
(dest / 'identity.json').write_text(json.dumps(identity, indent=2) + '\n')
for f in dest.rglob('*'):
    if f.is_file():
        f.chmod(f.stat().st_mode & ~0o222)
print(json.dumps({k: v for k, v in identity.items() if k not in ('source_files_sha256', 'artifacts_sha256')}, indent=2))
