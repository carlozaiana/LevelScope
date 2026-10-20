from pathlib import Path
import subprocess, hashlib, json, struct
R=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation')
F=R/'build/prores-20261020/fixtures'
G='/private/tmp/levelscope-prores-20261020/ProResFixtureBuilder'
manifest=json.loads((F/'manifest.json').read_text())
def record(p,**metadata):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 manifest[p.name]=dict(bytes=p.stat().st_size,sha256=h.hexdigest(),**metadata)
 (F/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
def transform(data,variant):
 def atoms(data,start=0):
  out=bytearray(data[:start]); p=start
  while p<len(data):
   if len(data)-p<8:
    out+=data[p:];break
   n,kind=struct.unpack_from('>I4s',data,p); assert n>=8 and p+n<=len(data)
   body=data[p+8:p+n]
   if kind in (b'moov',b'trak',b'mdia',b'minf',b'stbl'): body=atoms(body)
   elif kind==b'stsd': body=atoms(body,8)
   elif kind in (b'apco',b'apcs',b'apcn'):body=atoms(body,78)
   elif kind==b'colr':
    assert body[:4]==b'nclc'
    if variant=='hdr':body=b'nclc'+struct.pack('>HHH',1,16,1)
    elif variant=='nclx-full-range':body=b'nclx'+struct.pack('>HHH',1,1,1)+b'\x80'
    elif variant=='nclx-709':body=b'nclx'+struct.pack('>HHH',1,1,1)+b'\0'
    elif variant=='duplicate-colour':
     out+=data[p:p+n]
   elif kind==b'pasp' and variant=='nonsquare':body=struct.pack('>II',2,1)
   elif kind==b'clap' and variant=='cropped':
    body=struct.pack('>I',1918)+body[4:]
   out+=struct.pack('>I4s',len(body)+8,kind)+body;p+=n
  return bytes(out)
 # The seed's moov is after media data, so growing it leaves all media offsets intact.
 assert data.rfind(b'moov')>data.find(b'mdat')
 return atoms(data)
for profile in ('proxy','lt','422'):
 seed=(F/(profile+'-seed.mov')).read_bytes()
 for variant,reason in [('hdr','Rec.709'),('nclx-full-range','Rec.709'),
                        ('duplicate-colour','Rec.709'),('nonsquare','square pixels'),('cropped','full, centred aperture')]:
  p=F/(profile+'-'+variant+'.mov');assert not p.exists()
  p.write_bytes(transform(seed,variant));record(p,profile=profile,expect='reject',reason=reason,variant=variant)
 p=F/(profile+'-nclx-709.mov');assert not p.exists();p.write_bytes(transform(seed,'nclx-709'))
 record(p,profile=profile,frames=256,scale=30000,step=1001,generated=True,positive_policy=True)
for name,code in [('excluded-hq','apch'),('excluded-4444','ap4h')]:
 p=F/(name+'.mov');assert not p.exists()
 subprocess.run([G,'generate',str(p),code,'1920','1080','4','25','1','normal'],check=True)
 record(p,expect='reject',reason='Supported codecs',variant=name)
for name,m in manifest.items():
 if m.get('expect')=='reject' and 'reason' not in m:
  m['reason']='progressive' if m.get('variant')=='interlaced' else 'Rec.709'
(F/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print('Extra colour/profile boundary fixtures complete.')
