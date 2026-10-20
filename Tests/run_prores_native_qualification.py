"""Inaudible native ProRes qualification; every invocation uses a new output directory."""
from pathlib import Path
import argparse,subprocess,json,time,hashlib,platform
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--fixtures',type=Path,required=True);p.add_argument('--actual',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--loading',action='store_true');args=p.parse_args()
args.output.mkdir(parents=True,exist_ok=False)
manifest=json.loads((args.fixtures/'manifest.json').read_text());results=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
identity=dict(binary_sha256=sha(args.binary),architecture=platform.machine(),os=subprocess.check_output(['/usr/bin/sw_vers'],text=True),fixtures_manifest_sha256=sha(args.fixtures/'manifest.json'),actual_sha256=sha(args.actual))
(args.output/'identity.json').write_text(json.dumps(identity,indent=2)+'\n')
def run(name,cmd):
 start=time.monotonic();peak=0;samples=[]
 with (args.output/(name+'.log')).open('w') as log:
  process=subprocess.Popen([str(args.binary.resolve()),*map(str,cmd)],stdout=log,stderr=subprocess.STDOUT)
  while process.poll() is None:
   try:
    r=subprocess.check_output(['/bin/ps','-o','rss=','-p',str(process.pid)],text=True).strip()
    if r:
     rss=int(r)*1024;peak=max(peak,rss);samples.append([time.monotonic()-start,rss])
   except subprocess.CalledProcessError:pass
   if time.monotonic()-start>1800:process.terminate();process.wait();raise RuntimeError('Native test timeout')
   time.sleep(.5)
 rc=process.wait()
 result=dict(name=name,argv=list(map(str,cmd)),exit_code=rc,seconds=time.monotonic()-start,peak_rss_bytes=peak,rss_samples=samples)
 results.append(result);(args.output/'runs.json').write_text(json.dumps(results,indent=2)+'\n')
 print(name,'PASS' if rc==0 else 'FAIL',round(result['seconds'],2),'s',round(peak/1024**2,1),'MiB',flush=True)
 if rc:print((args.output/(name+'.log')).read_text(),flush=True)
 return rc==0
ok=True
for name,m in manifest.items():
 if m.get('required') or m.get('vfr') or m.get('positive_policy'):
  ok=run(name,[ '--prores',args.fixtures/name,m['frames'],m['scale'],m['step'],'vfr' if m.get('vfr') else 'normal']) and ok
 elif m.get('expect')=='reject':
  ok=run(name,['--reject',args.fixtures/name,m.get('reason','')]) and ok
 elif m.get('embedded_tone_hz'):
  ok=run(name,['--isolation',args.fixtures/name]) and ok
ok=run('actual-file',['--actual',args.actual]) and ok
ok=run('actual-file-isolation',['--isolation',args.actual]) and ok
if args.loading:ok=run('large-file-loading',['--loading',args.fixtures/'422-1800.mov']) and ok
(args.output/'summary.json').write_text(json.dumps(dict(passed=ok,identity=identity,test_count=len(results)),indent=2)+'\n')
raise SystemExit(0 if ok else 1)
