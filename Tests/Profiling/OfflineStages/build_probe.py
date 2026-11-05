"""Build two native harnesses against the same Release objects without rebuilding or editing the app."""
import argparse
import hashlib
from pathlib import Path
import re
import shlex
import subprocess

here=Path(__file__).resolve().parent
root=here.parents[2]
p=argparse.ArgumentParser()
p.add_argument('base',type=Path);p.add_argument('output',type=Path)
p.add_argument('--universal',action='store_true',help='Retain both macOS architectures for actual-platform qualification')
modes=p.add_mutually_exclusive_group()
modes.add_argument('--sampling-only',action='store_true')
modes.add_argument('--compare-peaks',action='store_true')
args=p.parse_args()
base=args.base.resolve(); out=args.output.resolve();out.mkdir(parents=True,exist_ok=True)
support='LevelScopeStandaloneSupport'; target='LevelScopeSectionProposalTests'
source=root/'Source/Standalone'
# The overlays must come from the same production sources as the linked objects.
cache=(base/'CMakeCache.txt').read_text()
compiler=next((base/'CMakeFiles').glob('*/CMakeCXXCompiler.cmake')).read_text()
def compiler_value(key):
    return re.search(r'set\(CMAKE_CXX_COMPILER_'+key+r' "([^"]*)"\)',compiler)[1]
def cache_value(key):
    return re.search('^'+key+r':[^=]+=([^\n]*)$',cache,re.M)[1]
identity_text=';'.join((compiler_value('ID'),compiler_value('VERSION'),cache_value('CMAKE_BUILD_TYPE'),
                        cache_value('CMAKE_CXX_FLAGS'),cache_value('CMAKE_CXX_FLAGS_RELEASE'),'no-fast-math','fp-contract-off'))
if 'LEVELSCOPE_ENABLE_BOUNDED_TRUE_PEAK:BOOL=' in cache:
    identity_text+=';bounded-peaks='+cache_value('LEVELSCOPE_ENABLE_BOUNDED_TRUE_PEAK')
for file in sorted((f for f in (root/'Source').rglob('*') if f.suffix in ('.cpp','.h','.mm')),key=str):
    identity_text+=hashlib.sha256(file.read_bytes()).hexdigest()+';'
expected=hashlib.sha256(identity_text.encode()).hexdigest()
assert expected in (base/'CMakeFiles'/f'{support}.dir'/'flags.make').read_text(), 'Rebuild the base Release target for these sources first'
def compile_one(src, tgt, extra=()):
    flags=(base/'CMakeFiles'/f'{tgt}.dir'/'flags.make').read_text().splitlines()
    opts=[]
    for key in ('CXX_DEFINES','CXX_INCLUDES','CXX_FLAGS' if args.universal else 'CXX_FLAGSarm64'):
        opts+=shlex.split(next(line.split(' = ',1)[1] for line in flags if line.startswith(key+' = ')))
    obj=out/(src.stem+'.o')
    subprocess.run(['/usr/bin/c++',*opts,*extra,'-I'+str(here),'-I'+str(source),'-c',str(src),'-o',str(obj)],check=True)
    return str(obj)
def insert(text, marker, addition, before=False):
    assert text.count(marker)==1, marker
    return text.replace(marker,addition+'\n'+marker if before else marker+'\n'+addition)
overlay=[]
for filename in ('StandaloneMedia.cpp','VerifiedSource.cpp','StandaloneMeasurementService.cpp','StandaloneSectionProposal.cpp'):
    if args.sampling_only or args.compare_peaks: break
    text=(source/filename).read_text()
    if filename=='StandaloneMedia.cpp':
        text=insert(text,'constexpr int blockSize = 4096;', '    lsprofile::Lap scan;', True)
        text=insert(text,'    if (flacChecksum.isNotEmpty())','    scan.mark(lsprofile::decodeCheck);',True)
        text=insert(text,'    CancellableInput stream (*input, control);','    lsprofile::Lap hash;')
        text=insert(text,'    media.sha256 = juce::SHA256 (stream).toHexString();','    hash.mark(lsprofile::validationHash);')
    if filename=='VerifiedSource.cpp':
        text=insert(text,'    try\n    {','        lsprofile::Span verificationScope{lsprofile::verified};')
        text=insert(text,'        auto input = context.sourceFile.createInputStream();','        lsprofile::Lap prep;',True)
        text=insert(text,'        output.reset();','        prep.mark(lsprofile::copy);')
        text=insert(text,'        const auto hash = juce::SHA256 (hashInput).toHexString();','        prep.mark(lsprofile::copyHash);')
    if filename=='StandaloneMeasurementService.cpp':
        text, marker, tail = text.partition('template <typename Fn>')
        assert marker
        text=insert(text,'        if (! reader->read (audio.getArrayOfWritePointers(), audio.getNumChannels(), frame, count))','        lsprofile::Lap block;',True)
        text=insert(text,'        if (! processing.process (audio.getArrayOfWritePointers(), audio.getNumChannels(), count, frame))','        block.mark(lsprofile::read);',True)
        text=insert(text,'        if (observer) observer->observe (audio.getArrayOfReadPointers(), audio.getNumChannels(), count, frame);','        block.mark(lsprofile::gain);',True)
        text=insert(text,'        engine.push (audio.getArrayOfReadPointers(), count);','        block.mark(lsprofile::wave);',True)
        text=insert(text,'        engine.push (audio.getArrayOfReadPointers(), count);','        block.mark(lsprofile::meter);')
        text=insert(text,'        if(loudness)loudness->publish(frame);','        block.mark(lsprofile::publish);')
        text=insert(text,'    result.detail = engine.finish();','    lsprofile::Lap tail;',True)
        text=insert(text,'    result.detail = engine.finish();','    tail.mark(lsprofile::finish);')
        text+=marker+tail
    if filename=='StandaloneSectionProposal.cpp':
        text=insert(text,'            return worker.analyzeCandidate(context,**copy.value,pass,observer);',
                    '            lsprofile::Pass tagged(observer ? 1 : context.processing.staticGainDb == r.currentContext.processing.staticGainDb ? 2 : 3);',True)
    path=out/filename
    path.write_text('#include "Trace.h"\n'+text)
    overlay.append(compile_one(path,support))
if args.sampling_only:
    # Let macOS load symbols before the short import and Source worker have ended.
    # Only this third launcher waits. Its actual trial links the unmodified services.
    sampling_source=out/'SampleDriver.cpp'
    sampling_source.write_text('#include <chrono>\n#include <thread>\n'
        'namespace { struct Delay { Delay(){std::this_thread::sleep_for(std::chrono::seconds(7));} } delay; }\n'
        '#include "'+str(here/'Profile.cpp')+'"\n')
    driver=compile_one(sampling_source,target)
else:
    driver=compile_one(here/'Profile.cpp',target)
link=shlex.split((base/'CMakeFiles'/f'{target}.dir'/'link.txt').read_text())
if not args.universal:
    i=link.index('x86_64');del link[i-1:i+1]
variants=[('sample-control',[])] if args.sampling_only else [('profile',overlay),('control',[])]
if args.compare_peaks:
    # Only the service's execution choice differs. Keep the same DSP, build
    # binding, driver and optimization flags so exact identities are comparable.
    variants=[('app-mode',[])]
    for enable,name in enumerate(('serial','bounded')):
        src=out/(name+'-measurement.cpp')
        # Both supported macOS architectures now use the production policy.
        # Only this private oracle overrides the build switch; app-mode links
        # the unchanged production service object from the ordinary app build.
        src.write_text('#include "'+str(source/'StandaloneMeasurementService.cpp')+'"\n')
        obj=compile_one(src,support,['-ULEVELSCOPE_ENABLE_BOUNDED_TRUE_PEAK',f'-DLEVELSCOPE_ENABLE_BOUNDED_TRUE_PEAK={enable}'])
        variants.append((name,[obj]))
for name,objects in variants:
    cmd=link[:]; idx=next(i for i,a in enumerate(cmd) if a.endswith('/Tests/SectionProposalTests.cpp.o'))
    cmd[idx:idx+1]=[driver,*objects]
    cmd.remove(next(a for a in cmd if a.endswith('/VideoReaderAudit_mac.mm.o')))
    cmd[cmd.index('-o')+1]=str(out/name)
    subprocess.run(cmd,cwd=base,check=True)
print('Built native profiling harness:',out)
