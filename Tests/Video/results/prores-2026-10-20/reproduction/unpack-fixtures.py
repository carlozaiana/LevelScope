"""Expand exact bytes and verify every file against the original fixture manifest."""
from pathlib import Path
import json,hashlib,shutil,time
R=Path('/Users/fuzz/LevelScope-ProRes-qualification-2026-10-20')
P=R/'fixture-transport';F=R/'fixtures';F.mkdir(exist_ok=True)
manifest=json.loads((P/'manifest.json').read_text());recipes=json.loads((P/'recipes.json').read_text())
cycles={p.name:p.read_bytes() for p in P.glob('*.cycle')}
results={}
for name,m in manifest.items():
 target=F/name;start=time.monotonic()
 if target.exists():
  archived=target.with_name(target.name+'.interrupted-scp-20261020')
  assert not archived.exists();target.rename(archived)
 h=hashlib.sha256();recipe=recipes[name]
 with target.open('xb') as out:
  def write(b):out.write(b);h.update(b)
  if 'full' in recipe:
   with (P/recipe['full']).open('rb') as f:
    for b in iter(lambda:f.read(4*1024*1024),b''):write(b)
  else:
   write((P/(name+'.prefix')).read_bytes());cycle=cycles[recipe['cycle']]
   for _ in range(recipe['copies']):write(cycle)
   write(cycle[:recipe['remainder']]);write((P/(name+'.suffix')).read_bytes())
 actual=h.hexdigest();ok=actual==m['sha256'] and target.stat().st_size==m['bytes']
 results[name]=dict(passed=ok,sha256=actual,bytes=target.stat().st_size,seconds=time.monotonic()-start)
 (R/'fixture-expansion.json').write_text(json.dumps(results,indent=2)+'\n')
 print(name,'PASS' if ok else 'FAIL',round(time.monotonic()-start,2),flush=True)
 if not ok:raise SystemExit('Fixture integrity failure')
shutil.copy2(P/'manifest.json',F/'manifest.json')
print('All fixture bytes match original identities.')
