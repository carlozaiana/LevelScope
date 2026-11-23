"""Run the finite device-free suite against an immutable candidate (Python 3.8+)."""
import argparse,hashlib,json,os,platform,subprocess,sys,time
from pathlib import Path
def sha(path):
 h=hashlib.sha256()
 with open(str(path),'rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def verify(package,expected):
 identity=package/'identity.json'
 assert sha(identity)==expected, 'candidate identity mismatch'
 m=json.loads(identity.read_text()); files=m['artifacts_sha256']
 for name,digest in files.items():
  assert sha(package/name)==digest,name
 return m,len(files)
p=argparse.ArgumentParser()
p.add_argument('package',type=Path);p.add_argument('expected');p.add_argument('out',type=Path)
p.add_argument('--benchmark',type=Path);p.add_argument('--timing',nargs=2,type=Path)
p.add_argument('--only',nargs='+',choices=['service','playback-ui','feedback','clock','propose','timing-stereo','timing-surround','audition','devices','timeline','ordinary-ui'])
args=p.parse_args(); package=args.package.resolve(); dest=args.out.resolve()
m,count=verify(package,args.expected); dest.mkdir(parents=True,exist_ok=False)
def runner(label,binary,cmd,receipt=None):
 t=time.monotonic()
 with open(str(dest/(label+'.log')),'w') as log:
  proc=subprocess.run([str(package/'bin'/binary)]+[str(s) for s in cmd],stdout=log,stderr=subprocess.STDOUT,timeout=300)
 row={'group':label,'binary':binary,'exit':proc.returncode,'seconds':time.monotonic()-t}
 if receipt and receipt.is_file():
  d=json.loads(receipt.read_text());row.update(assertions=len(d.get('checks',[])),result=d.get('passed') is True,receipt=receipt.name)
 elif not receipt:
  import re
  match=re.search(r'PASS:\s*(\d+) video timing/transport checks',(dest/(label+'.log')).read_text())
  row.update(assertions=int(match.group(1)) if match else 0,result=bool(match))
 else:row['result']=False
 print(label,proc.returncode,row.get('assertions'),row['result'],flush=True)
 return row
only=args.only
def enabled(label):return only is None or label in only
rows=[]
if enabled('service'):rows.append(runner('service','LevelScopeSharedSourceTests',[dest/'service.json'],dest/'service.json'))
if enabled('feedback'):rows.append(runner('feedback','LevelScopePreparationFeedbackTests',[dest/'feedback'],dest/'feedback/result.json'))
if enabled('audition'):rows.append(runner('audition','LevelScopeAuditionTests',['--report',dest/'audition.json'],dest/'audition.json'))
if enabled('devices'):rows.append(runner('devices','LevelScopeDeviceSettingsTests',[dest/'devices'],dest/'devices/result.json'))
if enabled('timeline'):rows.append(runner('timeline','LevelScopeTimelineProgressTests',['--root',dest/'timeline','--output',dest/'timeline.json'],dest/'timeline.json'))
if enabled('clock'):rows.append(runner('clock','LevelScopeVideoTests',[]))
if enabled('ordinary-ui'):rows.append(runner('ordinary-ui','LevelScopeRegionUITests',[dest/'ordinary-ui'],dest/'ordinary-ui/result.json'))
if enabled('playback-ui'):rows.append(runner('playback-ui','LevelScopeRegionUITests',[dest/'playback-ui','--playback-ui'],dest/'playback-ui/result.json'))
if args.benchmark and enabled('propose'):rows.append(runner('propose','LevelScopeSharedSourceTests',['--benchmark',args.benchmark,dest/'propose.json'],dest/'propose.json'))
if args.timing:
 for name,file in zip(('timing-stereo','timing-surround'),args.timing):
  if enabled(name):rows.append(runner(name,'LevelScopePlaybackTiming',[file,dest/(name+'.json')],dest/(name+'.json')))
_,post=verify(package,args.expected)
d={'os':platform.mac_ver()[0],'build':subprocess.check_output(['sw_vers','-buildVersion'],universal_newlines=True).strip(),'architecture':platform.machine(),'candidate':args.expected,'implementation':m['source_commit'],'measurement_build':m['measurement_build'],'verified_package_files_before':count,'verified_package_files_after':post,'groups':rows,'passed':bool(rows) and all(x['exit']==0 and x['result'] for x in rows),'physical_outputs_opened':False}
(dest/'receipt.json').write_text(json.dumps(d,indent=2)+'\n')
sys.exit(0 if d['passed'] else 1)
