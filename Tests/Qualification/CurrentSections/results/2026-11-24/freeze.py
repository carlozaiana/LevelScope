"""Freeze a clean committed implementation and the exact Ninja-built native binaries used in qualification."""
import hashlib, json, os, re, shutil, subprocess, sys
from pathlib import Path
root=Path(__file__).resolve().parent
repo=root.parent.parent
build=root/'compile'
dest=root/sys.argv[1]
targets=sys.argv[2:]
def run(args,cwd=repo):return subprocess.check_output(args,cwd=cwd,text=True).strip()
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
assert not dest.exists()
assert not run(['git','status','--short','--untracked-files=no'])
tip=run(['git','rev-parse','HEAD'])
sources=run(['git','ls-files','Source','Tests','CMakeLists.txt']).splitlines()
dest.mkdir()
for name in sources:
 p=Path(name)
 if 'results' in p.parts or p.name=='.DS_Store':continue
 f=dest/'source'/p;f.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(repo/p,f)
shutil.copytree(build/'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app',dest/'LevelScope Standalone.app',symlinks=True)
(dest/'bin').mkdir();(dest/'configuration').mkdir()
shutil.copy2(build/'CMakeCache.txt',dest/'configuration/CMakeCache.txt')
commands = []
identity=None
for target in ['LevelScopeStandaloneSupport','LevelScopeStandalone',*targets]:
 lines=run(['ninja','-t','commands',target],build).splitlines()
 own=[line for line in lines if ('CMakeFiles/'+target+'.dir/') in line]
 ids={m for line in own for m in re.findall(r'LEVELSCOPE_MEASUREMENT_BUILD[^a-f0-9]+([0-9a-f]{64})',line)}
 assert len(ids)==1,(target,ids)
 which=next(iter(ids))
 if identity is None:identity=which
 assert identity==which
 (dest/'configuration'/(target+'-commands.txt')).write_text('\n'.join(own)+'\n')
 if target in targets:shutil.copy2(build/(target+'_artefacts/Release')/target,dest/'bin'/target)
binaries={}
for path in [dest/'LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone',*sorted((dest/'bin').iterdir())]:
 arch=run(['lipo','-archs',path]).split()
 assert set(arch)=={'arm64','x86_64'}
 load=run(['otool','-arch','x86_64','-l',path])
 assert re.search(r'LC_VERSION_MIN_MACOSX\s+cmdsize \d+\s+version 10\.13\b|platform \d+\s+minos 10\.13\b',load)
 binaries[str(path.relative_to(dest))]={'sha256':sha(path),'architectures':arch,'intel_minimum':'10.13'}
app=Path('LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone')
assert sha(dest/app)==sha(build/'LevelScopeStandalone_artefacts/Release'/app)
# Preserve the exact local development link product, including its normal
# linker signature. This is a local test package, with no distribution signing.

files={str(p.relative_to(dest)):sha(p) for p in sorted(dest.rglob('*')) if p.is_file()}
manifest={k.removeprefix('source/'):v for k,v in files.items() if k.startswith('source/')}
identity_doc={'kind':'Initialize absent Current on Propose and create sections from applied regions','source_commit':tip,'source_tree':run(['git','rev-parse','HEAD:Source']),'qualified_main':'c29b19940253d615ba4c3d3909eafa9b32bcc2f2','service_commit':run(['git','rev-list','-1','HEAD','--','Source/Standalone/VerifiedSource.cpp']),'measurement_build':identity,'source_manifest_sha256':hashlib.sha256(json.dumps(manifest,sort_keys=True,separators=(',',':')).encode()).hexdigest(),'binaries':binaries,'artifacts_sha256':files}
(dest/'identity.json').write_text(json.dumps(identity_doc,indent=2)+'\n')
receipt=sha(dest/'identity.json')
(dest/'identity.sha256').write_text(receipt+'  identity.json\n')
for path in dest.rglob('*'):
 if path.is_file():path.chmod(path.stat().st_mode & ~0o222)
print(tip,identity,receipt,len(files))
