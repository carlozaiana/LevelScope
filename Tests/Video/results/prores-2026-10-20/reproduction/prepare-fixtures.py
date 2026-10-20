from pathlib import Path
import subprocess,hashlib,json,math,shutil,struct,wave
R=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation')
B=R/'build/prores-20261020'; F=B/'fixtures'; F.mkdir(exist_ok=False)
TOOL=Path('/private/tmp/levelscope-prores-20261020/ProResFixtureBuilder')
MUX=Path('/Users/fuzz/Developer/LevelScope/build/video-implementation/probe/LevelScopeVideoFixtureBuilder')
def run(args):
 print(' '.join(map(str,args)),flush=True)
 subprocess.run(list(map(str,args)),check=True)
manifest={}
def record(p,**metadata):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 manifest[p.name]=dict(bytes=p.stat().st_size,sha256=h.hexdigest(),**metadata)
 (F/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
audio=F/'embedded-distinctive-1973hz.wav'
with wave.open(str(audio),'wb') as w:
 w.setparams((2,2,48000,48000*12,'NONE','not compressed'))
 block=bytearray()
 for n in range(48000*12):
  v=round(24000*math.sin(2*math.pi*1973*n/48000));block+=struct.pack('<hh',v,-v)
 w.writeframes(block)
record(audio,role='Distinctive embedded 1973 Hz tone; never physically played')
for profile,code in [('proxy','apco'),('lt','apcs'),('422','apcn')]:
 seed=F/(profile+'-seed.mov')
 if profile=='proxy':shutil.copy2('/private/tmp/levelscope-prores-20261020/proxy-seed.mov',seed)
 else:run([TOOL,'generate',seed,code,1920,1080,256,30000,1001,'normal'])
 record(seed,profile=profile,frames=256,scale=30000,step=1001,generated=True)
 for scale,step in [(24000,1001),(24,1),(25,1),(30000,1001),(30,1)]:
  count=math.ceil(20*scale/step);p=F/f'{profile}-{scale}-{step}.mov'
  run([TOOL,'repeat',seed,p,count,scale,step]);record(p,profile=profile,frames=count,scale=scale,step=step,generated=True,required=True)
 p=F/(profile+'-with-tone.mov')
 run([MUX,'mux',seed,audio,p]);record(p,profile=profile,frames=256,embedded_tone_hz=1973,generated=True)
 p=F/(profile+'-vfr.mov')
 run([TOOL,'generate',p,code,1920,1080,120,30000,1001,'vfr']);record(p,profile=profile,frames=120,scale=30000,step=1001,generated=True,vfr=True)
 # Explicit policy negatives preserve a valid compressed video payload.
 data=seed.read_bytes()
 for variant in ['missing-colour','conflicting-colour','interlaced']:
  altered=bytearray(data)
  if variant.startswith('missing'):
   pos=altered.rfind(b'colrnclc');assert pos>=0;altered[pos:pos+4]=b'free'
  elif variant.startswith('conflicting'):
   pos=altered.rfind(b'colrnclc');assert pos>=0;altered[pos+13]=6
  else:
   pos=altered.rfind(b'fiel');assert pos>=0;altered[pos+4]=2
  p=F/(profile+'-'+variant+'.mov');p.write_bytes(altered);record(p,profile=profile,expect='reject',variant=variant)
for rate in [44100,48000]:
 p=F/f'selected-clock-{rate}.wav'
 run(['/usr/bin/python3','-B',R/'Tests/generate_video_sync_audio.py',p,'--sample-rate',rate]);record(p,role='Selected source clock fixture, routes physically muted',sample_rate=rate)
print('Short fixture preparation complete',flush=True)
