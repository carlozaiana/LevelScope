"""Freeze an owned universal package; preserve all development and historical builds."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
from run_short import REQUIRED,sha

p=argparse.ArgumentParser();p.add_argument('build',type=Path);p.add_argument('destination',type=Path);a=p.parse_args()
repo=Path(__file__).resolve().parents[3];base=a.build.resolve();dest=a.destination.resolve();assert not dest.exists();dest.mkdir(parents=True)
compile=base/'compile';drivers=base/'drivers'
source_files=[Path(f) for f in subprocess.check_output(['git','ls-files','-z','--cached','--others','--exclude-standard','Source','Tests','CMakeLists.txt','Docs/ProjectIndex/63_M18InteractiveTimelineAndProgressiveAnalysis.md','Docs/ProjectIndex/64_M18UITimelineHandover.md'],cwd=repo).decode().split('\0') if f]
source_files=sorted(set(f for f in source_files if f.name!='.DS_Store' and '__pycache__' not in f.parts and 'results' not in f.parts))
for f in source_files:
    if not (repo/f).is_file():continue
    to=dest/'source'/f;to.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(repo/f,to)
source_hashes={str(p.relative_to(dest/'source')):sha(p) for p in sorted((dest/'source').rglob('*')) if p.is_file()}
source_digest=hashlib.sha256(json.dumps(source_hashes,sort_keys=True,separators=(',',':')).encode()).hexdigest()
shutil.copytree(compile/'LevelScopeStandalone_artefacts/Release/LevelScope Standalone.app',dest/'LevelScope Standalone.app',symlinks=True)
(dest/'bin').mkdir()
for target in ('Measurement','Compliance','Proposal','Diagnostics','Sections','Audition','Multichannel','Session','SessionPlayback','Export','SectionProcessing','SectionProposal','Mono','Marker','JobRequirements','Waveform','TimelineProgress','Video'):
    name='LevelScope'+target+'Tests';shutil.copy2(compile/(name+'_artefacts/Release/'+name),dest/'bin'/name)
for name in ('LevelScopeWaveformCoreTests','LevelScopePeakChannelParity'):shutil.copy2(compile/name,dest/'bin'/name)
for name in ('short-acceptance','long-workflow','session-exchange','video-pressure'):shutil.copy2(drivers/name,dest/'bin'/name)
(dest/'bin/peaks').mkdir()
for name in ('serial','bounded','app-mode'):shutil.copy2(drivers/'peaks'/name,dest/'bin/peaks'/name)
shutil.copytree(compile/'LevelScopeVideoQualification_artefacts/Release/LevelScope Video Qualification.app',dest/'LevelScope Video Qualification.app',symlinks=True)
(dest/'configuration').mkdir();shutil.copy2(compile/'CMakeCache.txt',dest/'configuration/CMakeCache.txt')
for target in ('LevelScopeStandalone','LevelScopeStandaloneSupport','LevelScopeMeasurementCore','LevelScopeSectionProposalTests','LevelScopeMonoTests'):
    for name in ('flags.make','link.txt'):
        path=compile/('CMakeFiles/'+target+'.dir')/name
        if path.exists():shutil.copy2(path,dest/'configuration'/(target+'-'+name))
compiler=next((compile/'CMakeFiles').glob('*/CMakeCXXCompiler.cmake'));shutil.copy2(compiler,dest/'configuration/compiler.cmake')
subprocess.run([str(Path('/usr/bin/python3')),str(dest/'source/Tests/Profiling/ShortAcceptance/fixtures.py'),str(dest/'fixtures')],check=True)
subprocess.run([str(Path('/usr/bin/python3')),str(dest/'source/Tests/Profiling/OfflineStages/fixtures.py'),str(dest/'peak-fixtures')],check=True)
flags=(compile/'CMakeFiles/LevelScopeStandalone.dir/flags.make').read_text();build=re.search(r'LEVELSCOPE_MEASUREMENT_BUILD=\\"([a-f0-9]{64})',flags)[1]
artifacts={str(f.relative_to(dest)):sha(f) for f in sorted(dest.rglob('*')) if f.is_file()}
identity=dict(schema=1,name='M18 combined qualification '+dest.name,source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip(),
              source_digest=source_digest,source_files_sha256=source_hashes,measurement_build=build,
              bounded_true_peak='ON; app service uses bounded stereo/5.1 on Apple silicon, serial on Intel',
              architectures=['arm64','x86_64'],intel_minimum='10.13',
              juce_commit=subprocess.check_output(['git','-C','/Users/fuzz/Developer/LevelScope-dependencies/JUCE-7.0.12','rev-parse','HEAD'],text=True).strip(),
              compiler=subprocess.check_output(['/usr/bin/c++','--version'],text=True),
              required_short_groups=sorted(REQUIRED),required_other_groups=['peak-policy-exact','cross-mac-measurements','exchange-both-directions','video-native-pressure','video-injected-matrix','mono-three-hour-both-macs','rf64-independent-oracle','physical-muted-coordinated','manual-native-review'],
              artifacts_sha256=artifacts,qualification_complete=False)
(dest/'identity.json').write_text(json.dumps(identity,indent=2)+'\n')
(dest/'identity.sha256').write_text(sha(dest/'identity.json')+'  identity.json\n')
for f in dest.rglob('*'):
    if f.is_file():f.chmod(f.stat().st_mode & ~0o222)
print(json.dumps({k:v for k,v in identity.items() if k not in ('artifacts_sha256','source_files_sha256')},indent=2))
