import hashlib
import json
import shlex
import subprocess
import sys
from pathlib import Path

out = Path(__file__).resolve().parent
ui = Path('/Users/fuzz/.codex/worktrees/a411/LevelScope')
candidate = ui / 'build/m20-competing-boundary-ui-2026-11-21/candidate-01'
compile_dir = ui / 'build/m20-ui-integration/compile'
target = 'LevelScopeEventPlanningTests'
manifest = json.loads((candidate / 'identity.json').read_text())
def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for name, value in manifest['source_sha256'].items():
    assert digest(candidate/'source'/name) == value, name
link = shlex.split((candidate/'configuration'/f'{target}-link.txt').read_text())
link[link.index('-o')+1] = str(out/target)
subprocess.run(link, cwd=compile_dir, check=True)
assert digest(out/target)==manifest['artifacts_sha256'][f'bin/{target}'], 'Linked original differs from frozen executable'
if '--proof' in sys.argv:
    print('All frozen sources verified; original compiled link inputs reproduce exact frozen executable.')
    sys.exit()
settings = {}
for row in (candidate/'configuration'/f'{target}-flags.make').read_text().splitlines():
    if ' = ' in row:
        name, value = row.split(' = ',1)
        settings[name] = shlex.split(value)
args = settings['CXX_DEFINES'] + settings['CXX_INCLUDES'] + settings['CXX_FLAGS']
args = [str(arg).replace(f'-I{ui}/Source',f'-I{candidate}/source/Source') for arg in args]
source = out/'qualified.cpp'
obj = out/'qualified.o'
subprocess.run(['/usr/bin/c++',*args,f'-I{candidate}/source/Tests','-c',source,'-o',obj],cwd=out,check=True)
link[link.index('-o')+1] = str(out/'m20-offline')
idx = next(i for i,a in enumerate(link) if a.endswith('/EventPlanningTests.cpp.o'))
link[idx] = str(obj)
subprocess.run(link,cwd=compile_dir,check=True)
(out/'binary-identity.json').write_text(json.dumps(dict(source_commit=manifest['source_commit'],measurement_build=manifest['measurement_build'],protocol=digest(out/'protocol.txt'),harness_cpp=digest(source),build_driver=digest(__file__),supplemental_universal=digest(out/'m20-offline'),original_executable_relinked=digest(out/target),original_frozen_executable=digest(candidate/'bin'/target)),indent=2)+'\n')
