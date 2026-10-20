from pathlib import Path
import hashlib,json,struct,shutil
F=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020/fixtures')
manifest=json.loads((F/'manifest.json').read_text());before=hashlib.sha256((F/'manifest.json').read_bytes()).hexdigest()
shutil.copy2(F/'manifest.json',F/'manifest-before-crop-correction.json')
def atoms(data,start=0):
 out=bytearray(data[:start]);p=start
 while p<len(data):
  if len(data)-p<8:out+=data[p:];break
  n,k=struct.unpack_from('>I4s',data,p);assert n>=8 and p+n<=len(data)
  b=data[p+8:p+n]
  if k in (b'moov',b'trak',b'mdia',b'minf',b'stbl'):b=atoms(b)
  elif k==b'stsd':b=atoms(b,8)
  elif k in (b'apco',b'apcs',b'apcn'):
   # Insert an explicit 1918x1080 clean aperture before the optional trailing pad.
   assert b'clap' not in b
   b=b[:78]+struct.pack('>I4s8I',40,b'clap',1918,1,1080,1,0,1,0,1)+b[78:]
  out+=struct.pack('>I4s',len(b)+8,k)+b;p+=n
 return bytes(out)
for profile in ('proxy','lt','422'):
 old=manifest[profile+'-cropped.mov']
 assert old['sha256']==manifest[profile+'-seed.mov']['sha256']
 old.pop('expect');old['invalid_negative_fixture']=True
 old['note']='Superseded test data: no clap atom existed to mutate, so this is identical to the valid seed.'
 p=F/(profile+'-cropped-v2.mov');assert not p.exists()
 data=atoms((F/(profile+'-seed.mov')).read_bytes());p.write_bytes(data)
 manifest[p.name]=dict(bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),profile=profile,
                       expect='reject',reason='full, centred aperture',variant='explicit-1918-clean-aperture')
(F/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
record=dict(previous_manifest_sha256=before,corrected_manifest_sha256=hashlib.sha256((F/'manifest.json').read_bytes()).hexdigest(),
            implementation_changed=False,binaries_changed=False,affected_checks='Three cropped-aperture negative fixture checks; all previous fixture bytes preserved.')
(F/'crop-fixture-correction.json').write_text(json.dumps(record,indent=2)+'\n')
print(record)
