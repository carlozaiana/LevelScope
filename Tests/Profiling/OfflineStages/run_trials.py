"""Warm-up, then sequential matched profile/control trials, with alternating layout and binary order."""
import argparse,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('bin',type=Path);p.add_argument('fixtures',type=Path);p.add_argument('destination',type=Path)
a=p.parse_args()
dest=a.destination.resolve();dest.mkdir(parents=True,exist_ok=False)
expected={}
def run(channel,binary,label):
    path=dest/f'{label}-{channel}-{binary}.json'
    cmd=[str(a.bin.resolve()/binary),str(a.fixtures.resolve()/f'programme-{channel}.wav'),str(path)]
    with path.with_suffix('.log').open('w') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    d=json.loads(path.read_text())
    assert d['identities']==expected.setdefault(channel,d['identities']), f'changed results: {path}'
    assert all(bool(r['timings_ms'])==(binary=='profile') for r in d['operations'].values() if r is not d['operations']['global'])
    print(label,channel,binary,*(f'{k}={r["wall_ms"]/1000:.3f}s' for k,r in d['operations'].items()),flush=True)
for n in (2,6):run(n,'profile','warmup')
for i in range(1,4):
    for n in ((2,6) if i%2 else (6,2)):
        for mode in (('control','profile') if (n==2)==(i%2==1) else ('profile','control')):
            run(n,mode,f'trial-{i}')
