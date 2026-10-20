from pathlib import Path
import json,ast,copy,hashlib,statistics
R=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation');B=R/'build/prores-20261020';P=B/'candidate-04'
p=B/'evidence/sequoia-memory-v2-300s-candidate04.json';base=json.loads(p.read_text());assert base['passed']
src=(P/'Tests/verify_prores_playback_evidence.py').read_text();tree=ast.parse(src);ns={};exec('import math',ns);exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='verify_memory'],type_ignores=[]),'verifier','exec'),ns);check=ns['verify_memory']
def rollup(d):
 first=[x for x in d['memory_samples'] if x['seconds']<62];last=[x for x in d['memory_samples'] if x['seconds']>d['requested_seconds']-60]
 d['first_minute_memory_sample_count']=len(first);d['last_minute_memory_sample_count']=len(last)
 for k,peak in [('rss_bytes','peak_sampled_playback_rss_bytes'),('non_reusable_resident_bytes','peak_non_reusable_resident_bytes'),('physical_footprint_bytes','peak_physical_footprint_bytes')]:
  d[peak]=max(x[k] for x in d['memory_samples']);d['first_minute_mean_'+k]=statistics.mean(x[k] for x in first) if first else 0;d['last_minute_mean_'+k]=statistics.mean(x[k] for x in last) if last else 0
rows=[]
def trial(name,change,want_pass):
 d=copy.deepcopy(base);change(d);errors=check(d);passed=not errors
 assert passed==want_pass,(name,errors)
 rows.append(dict(name=name,expected_pass=want_pass,observed_pass=passed,errors=errors))
trial('unaltered final 300-second report',lambda d:None,True)
trial('diagnostic report cannot qualify',lambda d:d.update(qualification_eligible=False),False)
trial('allocator relief cannot qualify',lambda d:d.update(allocator_relief_requested=True),False)
trial('VM query failure fails closed',lambda d:d.update(memory_query_failures=1),False)
trial('missing raw samples fails closed',lambda d:d.update(memory_samples=[]),False)
trial('falsified mean fails independent recomputation',lambda d:d.update(last_minute_mean_rss_bytes=1),False)
def delta(d,key,reusable=False):
 for x in d['memory_samples']:
  if x['seconds']>d['requested_seconds']-60:
   x[key]+=64*1024**2
   if key=='rss_bytes':x['reusable_bytes' if reusable else 'non_reusable_resident_bytes']+=64*1024**2
 rollup(d)
trial('non-reusable growth above limit rejected',lambda d:delta(d,'rss_bytes'),False)
trial('physical footprint growth above limit rejected',lambda d:delta(d,'physical_footprint_bytes'),False)
trial('reusable-only growth allowed below unchanged RSS ceiling',lambda d:delta(d,'rss_bytes',True),True)
def ceiling(d):
 d['memory_samples'][-1]['rss_bytes']=769*1024**2;d['memory_samples'][-1]['reusable_bytes']=d['memory_samples'][-1]['rss_bytes']-d['memory_samples'][-1]['non_reusable_resident_bytes'];rollup(d)
trial('total RSS ceiling remains enforced for reusable memory',ceiling,False)
def short(d):
 d['memory_samples']=[x for x in d['memory_samples'] if x['seconds']>=8];rollup(d)
trial('insufficient first-window samples rejected',short,False)
out=B/'evidence/memory-verifier-negative-checks-candidate04.json';assert not out.exists();out.write_text(json.dumps(dict(passed=True,source_report_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),verifier_sha256=hashlib.sha256(src.encode()).hexdigest(),tests=rows),indent=2)+'\n');print('PASS',len(rows),'independent acceptance checks')
