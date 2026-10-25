#!/usr/bin/env python3
"""Fixed remaining M7 device programme; no retries, no diagnostic tracing."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time
from m7_video_acceptance import evaluate, POLICY
from VideoPaint.run_comparison import audio_snapshot, check_device


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(4*1024*1024), b''): h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('candidate', 'fixtures', 'actual', 'vfr', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--platform', choices=('sequoia', 'high-sierra'), required=True)
    p.add_argument('--phase', choices=('short', 'endurance'), required=True)
    p.add_argument('--allow-audio', action='store_true')
    a = p.parse_args()
    if not a.allow_audio: p.error('Coordinated --allow-audio required; routes remain muted')
    arch, prefix, device, rate = ('arm64','15.','Mac Studio Speakers',48000) if a.platform == 'sequoia' else ('x86_64','10.13.','Built-in Output',44100)
    os_version = subprocess.check_output(['/usr/bin/sw_vers','-productVersion'], text=True).strip()
    if platform.machine() != arch or not os_version.startswith(prefix):
        raise RuntimeError('Actual platform guard failed')
    audio_before = audio_snapshot()
    check_device(audio_before, dict(device=device, rate=rate))
    identity = json.loads((a.candidate/'identity.json').read_text())
    for name, digest in identity['artifacts_sha256'].items():
        if sha(a.candidate/name) != digest: raise RuntimeError('Candidate changed: '+name)
    binary_key = 'LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification'
    binary = (a.candidate/binary_key).resolve()
    audio = (a.candidate/('m7-fixtures/clock-71-'+str(rate)+'.flac')).resolve()
    manifest = json.loads((a.fixtures/'manifest.json').read_text())
    cases = []
    if a.phase == 'short':
        for name, m in manifest.items():
            if m.get('required'):
                cases.append(('prores-current-'+Path(name).stem,a.fixtures/name,18,'current',m['scale']/m['step'],m['sha256']))
        if len(cases) != 15: raise RuntimeError('ProRes matrix must have 15 generated cases')
        cases.append(('prores-current-reference',a.actual,100,'current',25,sha(a.actual)))
        for domain in ('source','current'):
            cases += [('h264-cfr-'+domain,a.candidate/'video-fixtures/half-hour.mov',18,domain,30000/1001,identity['artifacts_sha256']['video-fixtures/half-hour.mov']),
                      ('h264-vfr-'+domain,a.vfr,18,domain,None,sha(a.vfr))]
    else:
        for profile in ('proxy','lt','422'):
            name = profile+'-1800.mov'; m = manifest[name]
            cases.append(('prores-current-'+profile+'-1800',a.fixtures/name,1800,'current',m['scale']/m['step'],m['sha256']))
    a.output.mkdir(parents=True, exist_ok=False)
    record = dict(policy=POLICY, source_digest=identity['source_digest'], identity_sha256=sha(a.candidate/'identity.json'),
        binary_sha256=sha(binary), architecture=arch, os=subprocess.check_output(['/usr/bin/sw_vers'],text=True),
        cpu=subprocess.check_output(['/usr/sbin/sysctl','-n','machdep.cpu.brand_string'],text=True).strip(),
        fixture_manifest_sha256=sha(a.fixtures/'manifest.json'), device=device, sample_rate=rate,
        all_routes_muted=True, phase=a.phase, script_sha256=sha(__file__),
        audio_configuration_before_sha256=hashlib.sha256(json.dumps(audio_before,sort_keys=True).encode()).hexdigest(),
        started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (a.output/'identity.json').write_text(json.dumps(record,indent=2)+'\n')
    runs = []
    for name, movie, duration, domain, fps, movie_hash in cases:
        if sha(movie) != movie_hash: raise RuntimeError('Fixture changed: '+str(movie))
        report = (a.output/(name+'.json')).resolve()
        gain = 0.0 if domain == 'source' else -3.305074118944
        argv = [str(binary),'--allow-audio','--muted','--device',device,'--audio',str(audio),
                '--video',str(movie.resolve()),'--seconds',str(duration),'--report',str(report),'--audition-state',domain]
        if domain == 'current': argv += ['--processing-gain',str(gain)]
        print('START',name,duration,'seconds; all routes muted',flush=True)
        start=time.monotonic()
        with (a.output/(name+'.log')).open('x') as log:
            completed=subprocess.run(argv,stdout=log,stderr=subprocess.STDOUT,timeout=duration+420)
        data=json.loads(report.read_text())
        expected=dict(architecture=arch,device=device,sample_rate=rate,
            audition_state='Source' if domain=='source' else 'Current State',processing_gain_db=gain,
            measurement_build=identity['measurement_build'], video_sha256=movie_hash,
            audio_sha256=identity['artifacts_sha256'][str(audio.relative_to(a.candidate.resolve()))],requested_seconds=duration)
        result=evaluate(data,expected,completed.returncode,fps)
        runs.append(dict(name=name,argv=argv,seconds=time.monotonic()-start,report_sha256=sha(report),
                         expected=expected,fps=fps,**result))
        (a.output/'runs.json').write_text(json.dumps(runs,indent=2)+'\n')
        print('FINISH',name,'ACCEPTED' if result['accepted'] else 'FAILED',result['failures'],
              'skips',data.get('actual_paint_skipped_frames'),flush=True)
        if not result['accepted']: break
    audio_unchanged = audio_before == audio_snapshot()
    passed=len(runs)==len(cases) and all(r['accepted'] for r in runs) and audio_unchanged
    (a.output/'summary.json').write_text(json.dumps(dict(accepted=passed,policy=POLICY,
        expected_cases=len(cases),completed_cases=len(runs),physical_71_listening='Pending; not assessed',
        audio_configuration_unchanged=audio_unchanged,
        total_picture_skips=sum(r.get('accepted_picture_skips') or 0 for r in runs)),indent=2)+'\n')
    return 0 if passed else 1


if __name__=='__main__': raise SystemExit(main())
