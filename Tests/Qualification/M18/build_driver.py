"""Compile a universal qualification launcher against the frozen product objects."""
import argparse
from pathlib import Path
import shlex
import subprocess

p = argparse.ArgumentParser()
p.add_argument('base', type=Path)
p.add_argument('source', type=Path)
p.add_argument('output', type=Path)
p.add_argument('--target', default='LevelScopeStandalone')
p.add_argument('--replace', default='StandaloneApp.cpp.o')
p.add_argument('--private', action='store_true')
a = p.parse_args()
base, source, output = a.base.resolve(), a.source.resolve(), a.output.resolve()
output.parent.mkdir(parents=True, exist_ok=True)
target = base / ('CMakeFiles/' + a.target + '.dir')
flags = (target / 'flags.make').read_text().splitlines()
options = []
for key in ('CXX_DEFINES', 'CXX_INCLUDES', 'CXX_FLAGS'):
    options += shlex.split(next(l.split(' = ', 1)[1] for l in flags if l.startswith(key + ' = ')))
obj = output.with_suffix('.o')
root = Path(__file__).resolve().parents[3]
subprocess.run(['/usr/bin/c++', *options, '-x', 'objective-c++',
                *(['-fno-access-control'] if a.private else []), '-I' + str(root / 'Tests'),
                '-I' + str(root / 'Source'), '-c', str(source), '-o', str(obj)], check=True)
link = shlex.split((target / 'link.txt').read_text())
matches = [i for i, x in enumerate(link) if x.endswith('/' + a.replace)]
assert len(matches) == 1, matches
link[matches[0]] = str(obj)
link[link.index('-o') + 1] = str(output)
subprocess.run(link, cwd=base, check=True)
