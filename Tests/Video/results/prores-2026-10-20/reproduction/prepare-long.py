from pathlib import Path
import subprocess,json,hashlib,time
F=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020/fixtures')
tool='/private/tmp/levelscope-prores-20261020/ProResFixtureBuilder'
for profile in ['proxy','lt','422']:
 p=F/(profile+'-1800.mov');assert not p.exists()
 start=time.monotonic(); print('Preparing',p.name,flush=True)
 subprocess.run([tool,'repeat',str(F/(profile+'-seed.mov')),str(p),'54000','30000','1001'],check=True)
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 manifest=json.loads((F/'manifest.json').read_text())
 manifest[p.name]=dict(profile=profile,frames=54000,scale=30000,step=1001,generated=True,
     sha256=h.hexdigest(),bytes=p.stat().st_size,generation_seconds=time.monotonic()-start,
     content='256 detailed moving frames repeated with exact unique PTS and 8-bit modulo-256 burn-in')
 (F/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print('Ready',p.name,p.stat().st_size,flush=True)
