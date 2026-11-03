"""Focused M18 proposal workflow and independent historical attachment checks."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import time

repo = Path(__file__).resolve().parents[4]
build = repo / 'build/m18-dev-01'
parser = argparse.ArgumentParser()
parser.add_argument('destination', type=Path)
parser.add_argument('--pilot', action='store_true')
args = parser.parse_args()
out = args.destination.resolve()
out.mkdir(parents=True, exist_ok=False)
tmp = Path(tempfile.mkdtemp(prefix='levelscope-m18-proposals-', dir='/tmp'))
def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()
def write(name, value):
    (out/name).write_text(json.dumps(value, indent=2)+'\n')
def executable(name):
    return build/(name+'_artefacts/Release')/name
def manifest():
    paths = [p for folder in ['Source','Tests'] for p in (repo/folder).rglob('*') if p.suffix in ('.cpp','.h','.mm')]
    paths += [repo/'CMakeLists.txt', Path(__file__).resolve()]
    return {str(p.relative_to(repo)):sha(p) for p in sorted(paths)}
inputs = manifest();write('inputs.json',inputs)
identity = re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([0-9a-f]{64})',(build/'CMakeFiles/LevelScopeStandaloneSupport.dir/flags.make').read_text())[1]
runs = [('proposal-workflow','LevelScopeSectionProposalTests',['--m18-workflow']),
        ('section-smoke','LevelScopeSectionProposalTests',['--smoke']),
        ('preferences','LevelScopeSectionProposalTests',['--preferences']),
        ('global','LevelScopeProposalTests',[]),('session-short','LevelScopeSessionTests',['--short']),
        ('timeline','LevelScopeTimelineProgressTests',[]),('requirements','LevelScopeJobRequirementsTests',[])]
if args.pilot:runs=runs[:1]
records=[]
for name,target,options in runs:
    command=[str(executable(target)),*options]
    if target not in ['LevelScopeProposalTests','LevelScopeSessionTests']:command+=['--root',str(tmp/name)]
    command+=['--report' if target=='LevelScopeProposalTests' else '--output',str(out/(name+'.json'))]
    started=time.monotonic()
    with (out/(name+'.log')).open('w') as log:
        process=subprocess.run(command,cwd=repo,stdout=log,stderr=subprocess.STDOUT)
    path=out/(name+'.json');result=json.loads(path.read_text()) if path.is_file() else {}
    record=dict(group=name,passed=process.returncode==0 and result.get('passed',False) and result.get('build')==identity,
                checks=(result.get('checks',0) if isinstance(result.get('checks'),int) else len(result.get('checks',[]))),seconds=time.monotonic()-started,command=command,exit_code=process.returncode)
    records.append(record);print(name,record['passed'],record['checks'],result.get('failure',''),flush=True)
    if not record['passed']:break
for capture in tmp.rglob('*.png'):shutil.copy2(capture,out/capture.name)
app=build/'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app';binaries=[]
if not args.pilot:
    for exe in [app/'Contents/MacOS/LevelScope Standalone']+sorted({executable(t) for _,t,_ in runs}):
        archs=subprocess.check_output(['lipo','-archs',str(exe)],text=True).split()
        load=subprocess.check_output(['otool','-arch','x86_64','-l',str(exe)],text=True)
        assert set(archs)=={'arm64','x86_64'} and re.search(r'(version|minos) 10\.13\s',load)
        assert identity.encode() in exe.read_bytes()
        binaries.append(dict(path=str(exe),sha256=sha(exe),architectures=archs,intel_minimum='10.13'))
    assert inputs==manifest(),'input changed during final run'
write('summary.json',dict(all_passed=len(records)==len(runs) and all(r['passed'] for r in records),pilot=args.pilot,
    measurement_build=identity,records=records,app_path=str(app),binaries=binaries,fixtures=str(tmp),
    physical_device_opened=False,video_loaded=False))
raise SystemExit(0 if len(records)==len(runs) and all(r['passed'] for r in records) else 1)
