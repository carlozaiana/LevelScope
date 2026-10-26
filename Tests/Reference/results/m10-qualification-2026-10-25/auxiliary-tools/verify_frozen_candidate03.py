#!/usr/bin/env python3
import argparse,pathlib,json,hashlib,datetime
p=argparse.ArgumentParser();p.add_argument('--root',type=pathlib.Path,required=True);p.add_argument('--workspace',type=pathlib.Path);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args();cache={}
def sha(p):
 st=p.stat();k=(st.st_dev,st.st_ino,st.st_size,st.st_mtime_ns)
 if k not in cache:
  h=hashlib.sha256()
  with p.open('rb') as f:
   for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
  cache[k]=h.hexdigest()
 return cache[k]
records=[]
for name in ('candidate-01','candidate-02','candidate-03'):
 folder=a.root/name;identity=json.loads((folder/'identity.json').read_text());files=identity['artifacts_sha256'];selected=files if name=='candidate-03' else {n:h for n,h in files.items() if n.startswith(('source/','Tests/','bin/')) or '.app/' in n};errors=[]
 for n,h in selected.items():
  try:ok=sha(folder/n)==h
  except OSError:ok=False
  if not ok:errors.append(n)
 records.append({'candidate':name,'source_digest':identity['source_digest'],'identity_sha256':sha(folder/'identity.json'),'scope':'entire recorded package' if name=='candidate-03' else 'source, scripts, executables and app bundles; historical media not rehashed here','files_checked':len(selected),'errors':errors,'passed':not errors})
 if name=='candidate-03' and a.workspace:
  errors=[n for n,h in identity['source_files_sha256'].items() if sha(a.workspace/n)!=h];records.append({'scope':'current implementation against frozen manifest','passed':not errors,'errors':errors,'files_checked':len(identity['source_files_sha256'])})
result={'passed':all(x['passed'] for x in records),'checked_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'records':records}
with a.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] else 1)
