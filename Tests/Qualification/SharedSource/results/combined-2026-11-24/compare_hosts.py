"""Compare the two actual runs and both qualified references. Preserve binary64 differences."""
import base64
from collections import defaultdict
import hashlib
import io
import json
import struct
import sys
from pathlib import Path

root, local, remote, destination = map(Path, sys.argv[1:])
read = lambda file: json.loads(file.read_text())
a, b = [read(path / 'propose.json') for path in (local, remote)]
old_a = read(root / 'comparison/qualified/propose-03.json')
old_b = read(remote.parent / 'baseline-01/propose.json')

def decode(proposal):
    blob = base64.b64decode(proposal['rendered_numeric_base64'], validate=True)
    assert hashlib.sha256(blob).hexdigest() == proposal['rendered_numeric_sha256']
    stream = io.BytesIO(blob)
    def get(fmt):
        size = struct.calcsize('<' + fmt)
        return struct.unpack('<' + fmt, stream.read(size))[0]
    def opt():
        return get('d') if get('?') else None
    fields = {}
    fields['rate'] = get('i')
    fields['frames'] = get('q')
    for metric in ('integrated','lra','momentaryMax','shortTermMax','samplePeak','truePeak'):
        fields[metric + '.status'] = get('i')
        fields[metric] = opt()
    for name in ('integratedRelativeGate','lraRelativeGate','p10','p95'):
        fields[name] = opt()
    fields['integratedIncluded'] = get('q')
    fields['lraIncluded'] = get('q')
    fields['lraStable'] = get('?')
    for kind in ('momentary','shortTerm','lraWindows'):
        fields[kind+'.count'] = get('q')
        for i in range(fields[kind+'.count']):
            fields[f'{kind}.endFrame[{i}]'] = get('q')
            fields[f'{kind}.energy[{i}]'] = get('d')
    for kind, fmt in [('channelSamplePeaks','d'),('channelTruePeaks','d'),('samplePeakFrames','q'),('truePeakEighthFrames','q')]:
        fields[kind+'.count'] = get('q')
        for i in range(fields[kind+'.count']):
            fields[f'{kind}[{i}]'] = get(fmt)
    assert not stream.read()
    return fields

def exact(v):
    return v.hex() if isinstance(v,float) else v

left, right = decode(a), decode(b)
assert left.keys() == right.keys()
differences = defaultdict(list)
for name in left:
    x, y = left[name], right[name]
    if exact(x) != exact(y):
        family = name.split('[')[0]
        differences[family].append((name, x, y))
summary = {}
for family, entries in differences.items():
    name, x, y = max(entries, key=lambda item: abs(item[1] - item[2]))
    summary[family] = dict(count=len(entries), max_abs_delta=abs(x-y), field=name,
                           sequoia=x, high_sierra=y, sequoia_exact=exact(x), high_sierra_exact=exact(y))
keys = ('media_sha256','source_numeric_sha256','applied_plan','applied_choice','processed_current_numeric_sha256','rendered_numeric_sha256','rendered_pcm_sha256')
cross = {key: a[key] == b[key] for key in keys}
same_host = {name: {key: current[key] == old[key] for key in keys}
             for name,current,old in [('sequoia',a,old_a),('high_sierra',b,old_b)]}
attempt_keys = ('attempts','outcome','retained_attempt')
all_runs = [proposal['runs'] for proposal in (a,b,old_a,old_b)]
assert all(len(runs)==2 for runs in all_runs)
attempts = all(all(run[key] == a['runs'][i][key] for key in attempt_keys)
               for runs in all_runs for i,run in enumerate(runs))
assert all(all(row.values()) for row in same_host.values()) and attempts
assert all(match for key,match in cross.items() if key!='rendered_numeric_sha256')
assert a['build'] == b['build'] and old_a['build'] == old_b['build']
assert all(x['passed'] for x in (a,b,old_a,old_b))
result = dict(same_host_exact_against_qualified=same_host,cross_host_exact=cross,
              outcomes_order_and_all_attempt_receipts_equal=attempts,
              same_host_reference_also_differs_on_render_measurement=(old_a['rendered_numeric_sha256'] != old_b['rendered_numeric_sha256']),
              fields_compared=len(left),fields_differing=sum(len(v) for v in differences.values()),
              rendered_measurement_difference=summary,exported_file_sha256=a['rendered_pcm_sha256'])
assert not destination.exists()
destination.write_text(json.dumps(result,indent=2)+'\n')
