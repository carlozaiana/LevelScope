"""Lossless transport recipe for repeated native-encoded fixtures; no transcoding."""
from pathlib import Path
import hashlib,json,struct,shutil,tarfile
F=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020/fixtures')
P=Path('/private/tmp/levelscope-prores-20261020/fixture-transport')
P.mkdir(exist_ok=False)
manifest=json.loads((F/'manifest.json').read_text())
def parts(path):
 with path.open('rb') as f:
  p=0;mdat=None
  while p<path.stat().st_size:
   f.seek(p);v=f.read(16);n,k=struct.unpack('>I4s',v[:8]);header=8
   if n==1:n=struct.unpack('>Q',v[8:16])[0];header=16
   if n==0:n=path.stat().st_size-p
   assert n>=header
   if k==b'mdat':
    assert mdat is None
    mdat=(p+header,n-header)
   p+=n
  if not mdat:return None
  start,size=mdat;f.seek(0);prefix=f.read(start);f.seek(start+size);suffix=f.read()
  return prefix,start,size,suffix
cycles={}
for profile in ['proxy','lt','422']:
 p=F/(profile+'-seed.mov');pre,start,n,post=parts(p)
 with p.open('rb') as f:f.seek(start);cycles[profile]=f.read(n)
 (P/(profile+'.cycle')).write_bytes(cycles[profile])
recipes={}
for name,m in manifest.items():
 p=F/name
 if p.suffix=='.mov' and m.get('profile') in cycles:
  pre,start,n,post=parts(p);cycle=cycles[m['profile']]
  copies,remainder=divmod(n,len(cycle))
  h=hashlib.sha256(pre)
  for _ in range(copies):h.update(cycle)
  h.update(cycle[:remainder]);h.update(post)
  if h.hexdigest()==m['sha256']:
   (P/(name+'.prefix')).write_bytes(pre);(P/(name+'.suffix')).write_bytes(post)
   recipes[name]=dict(cycle=m['profile']+'.cycle',copies=copies,remainder=remainder)
   print('EXACT_RECIPE',name,flush=True);continue
 shutil.copy2(p,P/name);recipes[name]=dict(full=name)
 print('FULL_BYTES',name,flush=True)
(P/'recipes.json').write_text(json.dumps(recipes,indent=2)+'\n')
shutil.copy2(F/'manifest.json',P/'manifest.json')
with tarfile.open('/private/tmp/levelscope-prores-20261020/fixture-transport.tar.gz','w:gz',compresslevel=1) as archive:
 archive.add(P,arcname='fixture-transport')
print('Compact transport ready',flush=True)
