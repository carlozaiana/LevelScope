import json,sys
from pathlib import Path
base,new,out=map(Path,sys.argv[1:])
a,b=[json.loads(p.read_text()) for p in (base,new)]
assert a['passed'] and b['passed']
keys=['media_sha256','source_numeric_sha256','applied_plan','applied_choice','processed_current_numeric_sha256','rendered_numeric_sha256','rendered_pcm_sha256']
matches={k:a[k]==b[k] for k in keys}
runs=[]
for x,y in zip(a['runs'],b['runs']):
 for k in ['outcome','retained_attempt','attempts']:assert x[k]==y[k],k
 assert x['copy_passes']==x['validation_passes']==x['validation_hash_passes']==1
 assert y['copy_passes']==y['validation_passes']==y['validation_hash_passes']==0
 runs.append(dict(baseline_ms=x['elapsed_ms'],candidate_ms=y['elapsed_ms'],baseline_preparation_ms=x['preparation_ms'],candidate_preparation_ms=y['preparation_ms'],baseline_copy_passes=1,candidate_copy_passes=0,private_candidates=len(y['attempts']),retained_attempt=y['retained_attempt']))
d={'passed':all(matches.values()) and len(runs)==2,'qualified_measurement_build':a['build'],'candidate_measurement_build':b['build'],'input_content_sha256':b['media_sha256'],'exact_payloads_and_selection':matches,'exported_pcm_sha256':b['rendered_pcm_sha256'],'runs':runs}
assert d['passed'];assert not out.exists();out.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d))
