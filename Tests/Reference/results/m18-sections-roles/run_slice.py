"""Run the M18 section, role review and preset interaction checks with private fixtures."""
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
tmp = Path(tempfile.mkdtemp(prefix='levelscope-m18-section-import-', dir='/tmp'))
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
    paths = [p for p in (repo/'Source').rglob('*') if p.suffix in ('.cpp','.h','.mm')]
    paths += [repo/'CMakeLists.txt', repo/'Tests/TimelineProgressTests.cpp', repo/'Tests/TimelineGainLaneChecks.h',
              repo/'Tests/AuditionTestSupport.h', repo/'Tests/JobRequirementsTests.cpp',repo/'Tests/TimelineSectionImportChecks.h', Path(__file__).resolve()]
    return {str(p.relative_to(repo)):sha(p) for p in sorted(paths)}
inputs = manifest()
write('inputs.json',inputs)
identity = re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([0-9a-f]{64})',(build/'CMakeFiles/LevelScopeStandaloneSupport.dir/flags.make').read_text())[1]
runs = [('timeline','LevelScopeTimelineProgressTests',[])]
if not args.pilot: runs.append(('target-presets','LevelScopeJobRequirementsTests',[]))
records = []
for name,target,options in runs:
    command = [str(executable(target)), *options, '--root', str(tmp/name), '--output', str(out/(name+'.json'))]
    started = time.monotonic()
    with (out/(name+'.log')).open('w') as log:
        process = subprocess.run(command, cwd=repo, stdout=log, stderr=subprocess.STDOUT)
    p = out/(name+'.json')
    result = json.loads(p.read_text()) if p.is_file() else {}
    if name=='timeline' and not args.pilot:
        labels={c['name'] for c in result.get('checks',[]) if c['passed']}
        assert {'AIFF CHAN description order survives decoding','Edit order enables the selected file-channel editor',
                'Sections Space toggles playback and stop',
                'Save and Cancel both destroy preset dialog and release modality'}.issubset(labels)
    record = dict(group=name, passed=process.returncode==0 and result.get('passed',False) and (args.pilot or result.get('build')==identity),
                  checks=len(result.get('checks',[])), seconds=time.monotonic()-started, command=command, exit_code=process.returncode)
    records.append(record); print(record['group'],record['passed'],record['checks'],result.get('failure',''),flush=True)
for capture in (tmp/'timeline').glob('*.png'):
    shutil.copy2(capture,out/capture.name)
app = build/'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app'
binaries=[]
if not args.pilot:
    for exe in [app/'Contents/MacOS/LevelScope Standalone']+[executable(t) for _,t,_ in runs]:
        archs=subprocess.check_output(['lipo','-archs',str(exe)],text=True).split()
        load=subprocess.check_output(['otool','-arch','x86_64','-l',str(exe)],text=True)
        assert set(archs)=={'arm64','x86_64'} and re.search(r'(version|minos) 10\.13\s',load)
        assert identity.encode() in exe.read_bytes()
        binaries.append(dict(path=str(exe),sha256=sha(exe),architectures=archs,intel_minimum='10.13'))
    assert inputs==manifest(), 'input changed during run'
write('summary.json',dict(all_passed=all(r['passed'] for r in records),pilot=args.pilot,measurement_build=identity,
                          records=records,app_path=str(app),binaries=binaries,physical_device_opened=False,video_loaded=False))
raise SystemExit(0 if all(r['passed'] for r in records) else 1)
