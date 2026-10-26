#!/usr/bin/env python3
"""Re-evaluate preserved physical reports, retaining strict and accepted verdicts."""
import hashlib
import importlib.util
import json
import pathlib
import sys

r = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(r / 'candidate-03/Tests'))
from m7_video_acceptance import evaluate, POLICY

candidate = json.loads((r / 'candidate-03/identity.json').read_text())
rows, errors = [], []
for platform in ('sequoia', 'high-sierra'):
    for phase in ('short', 'endurance'):
        folder = r / (platform + '-physical-' + phase + '-candidate-03')
        required = [folder / name for name in ('identity.json', 'runs.json', 'summary.json')]
        if not all(p.is_file() for p in required):
            errors.append('Incomplete physical block: ' + folder.name)
            continue
        identity, runs, summary = [json.loads(p.read_text()) for p in required]
        if identity['source_digest'] != candidate['source_digest']:
            errors.append('Source identity: ' + folder.name)
        if identity['binary_sha256'] != candidate['artifacts_sha256']['LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification']:
            errors.append('Executable identity: ' + folder.name)
        if not summary['accepted'] or not summary['audio_configuration_unchanged']:
            errors.append('Block not accepted: ' + folder.name)
        if len(runs) != (4 if phase == 'short' else 1):
            errors.append('Incomplete case count: ' + folder.name)
        for run in runs:
            report = folder / (run['name'] + '.json')
            raw = json.loads(report.read_text())
            digest = hashlib.sha256(report.read_bytes()).hexdigest()
            if digest != run['report_sha256']:
                errors.append('Changed raw report: ' + str(report))
            verdict = evaluate(raw, run['expected'], run['strict_exit_code'], run['fps'], channels=2, layout='stereo')
            if not verdict['accepted']:
                errors.append('Unaccepted case: ' + str(report))
            rows.append(dict(platform=platform, phase=phase, name=run['name'],
                report=str(report.relative_to(r)), report_sha256=digest,
                strict_harness_passed=raw['passed'], accepted=verdict['accepted'],
                skipped_pictures=raw['actual_paint_skipped_frames'],
                skip_events=raw['skip_events'],
                source_seconds=raw['source_seconds'], underrun_blocks=raw['underrun_blocks'],
                maximum_paint_interval_ms=raw['maximum_paint_interval_ms'],
                maximum_frame_pts_skew_seconds=raw['maximum_frame_pts_skew_seconds'],
                first_to_last_mean_skew_difference_seconds=abs(raw['last_minute_mean_skew_seconds']-raw['first_minute_mean_skew_seconds']),
                embedded_audio_policy=raw['embedded_audio_policy'],
                memory_passed=raw['memory_passed'], memory_growth_assessed=raw['memory_growth_assessed'],
                first_minute_mean_rss_bytes=raw['first_minute_mean_rss_bytes'],
                last_minute_mean_rss_bytes=raw['last_minute_mean_rss_bytes'],
                peak_physical_footprint_bytes=raw['peak_physical_footprint_bytes'],
                peak_non_reusable_resident_bytes=raw['peak_non_reusable_resident_bytes']))
result = dict(passed=not errors and len(rows)==10, policy=POLICY, cases=rows, errors=errors,
    source_digest=candidate['source_digest'],
    scope='Approved muted two-output devices, internal audio-clock/UI evidence. Not optical synchronization or full multichannel listening qualification.')
(r / 'physical-video-summary-candidate-03.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(dict(passed=result['passed'], completed_cases=len(rows), errors=errors), indent=2))
