"""Separate M7 development acceptance; never edits the strict harness report.

The 2026-10-25 decision accepts occasional picture skips, not freezes or drift.
The 250 ms observation guard is a conservative stop/review guard for sustained
UI stalls, not a new synchronization tolerance or a promise of optical sync.
"""
import ast
import math
from pathlib import Path

POLICY = 'm7-development-picture-skips-2026-10-25'
MAX_PAINT_INTERVAL_MS = 250.0


def memory_verifier():
    path = Path(__file__).with_name('verify_prores_playback_evidence.py')
    function = next(n for n in ast.parse(path.read_text()).body
                    if isinstance(n, ast.FunctionDef) and n.name == 'verify_memory')
    namespace = {'math': math}
    exec(compile(ast.Module(body=[function], type_ignores=[]), str(path), 'exec'), namespace)
    return namespace['verify_memory']


def evaluate(report, expected, exit_code, fps=None):
    """Fail closed on missing/invalid evidence; preserve original strict verdict."""
    failures = [key for key, value in expected.items() if report.get(key) != value]
    required = dict(error='', all_routes_muted=True, selected_audio_channels=8,
        layout_capability='pcm-7.1-bs2051-system-i-v1',
        pre_monitor_nonzero_channel_mask=255, pre_monitor_channel_count_valid=True,
        blank_paints=0, bound_failures=0, paint_deadline_failures=0, underrun_blocks=0,
        attachment_within_budget=True, memory_passed=True, hal_output_timestamp=True,
        skipped_frame_limit=0)
    failures += [key for key, value in required.items() if report.get(key) != value]
    finite_fields = ('source_seconds', 'requested_seconds', 'maximum_paint_interval_ms',
        'maximum_frame_pts_skew_seconds', 'maximum_cpu_paint_ms', 'mean_skew_seconds',
        'first_minute_mean_skew_seconds', 'last_minute_mean_skew_seconds')
    finite = all(type(report.get(k)) in (int, float) and math.isfinite(report[k]) for k in finite_fields)
    if not finite:
        failures.append('missing/non-finite timing')
    else:
        if report['requested_seconds'] <= 0 or report['source_seconds'] < report['requested_seconds']:
            failures.append('incomplete duration')
        if not 0 < report['maximum_paint_interval_ms'] <= MAX_PAINT_INTERVAL_MS:
            failures.append('sustained paint stall / missing interval')
        if report['requested_seconds'] >= 120 and report.get('memory_growth_assessed') is not True:
            failures.append('missing memory growth assessment')
    for key in ('actual_paint_skipped_frames', 'unique_painted_frames', 'paints', 'pre_monitor_frames_observed'):
        if type(report.get(key)) is not int or report[key] < (0 if key == 'actual_paint_skipped_frames' else 1):
            failures.append('invalid counter ' + key)
    if not isinstance(report.get('audition_identity'), str) or not report['audition_identity']:
        failures.append('missing audition identity')
    skipped = report.get('actual_paint_skipped_frames', 0)
    events = report.get('skip_events')
    event_count = 0
    if not isinstance(events, list):
        failures.append('missing skip event timing')
    else:
        previous_time = -1.0
        for event in events:
            try:
                count = event['frame'] - event['previous_frame'] - 1
                stamp = event['audio_seconds']
                interval = event['paint_interval_ms']
                if (type(count) is not int or count <= 0 or not math.isfinite(stamp)
                    or stamp < previous_time or not math.isfinite(interval)
                    or not 0 < interval <= MAX_PAINT_INTERVAL_MS):
                    raise ValueError()
                previous_time = stamp
                event_count += count
            except (KeyError, TypeError, ValueError):
                failures.append('invalid skip event')
        if event_count != skipped:
            failures.append('incomplete skip event accounting')
    if type(skipped) is int and skipped >= 0:
        # A raw failure is acceptable ONLY when the preserved zero-skip gate is
        # the sole difference. The harness exit and strict flag must agree.
        if report.get('passed') is not (skipped == 0) or exit_code != (0 if skipped == 0 else 1):
            failures.append('strict verdict/exit unexplained by picture skips')
        if finite and fps:
            if report.get('unique_painted_frames', 0) + skipped < math.floor((report['requested_seconds'] - 2) * fps) - 1:
                failures.append('painted plus skipped coverage')
            drift = abs(report['last_minute_mean_skew_seconds'] - report['first_minute_mean_skew_seconds'])
            if drift > 1.0 / fps + 1e-9:
                failures.append('accumulating synchronization drift')
    try:
        failures += memory_verifier()(report)
    except (KeyError, TypeError, ValueError, ZeroDivisionError, OverflowError):
        failures.append('invalid memory evidence')
    return dict(policy=POLICY, accepted=not failures, failures=failures,
        strict_harness_passed=report.get('passed'), strict_exit_code=exit_code,
        accepted_picture_skips=skipped if not failures else None,
        maximum_paint_interval_ms=report.get('maximum_paint_interval_ms'),
        maximum_frame_pts_skew_seconds=report.get('maximum_frame_pts_skew_seconds'),
        scope='Internal audio-clock / UI selection; not optical sync or physical 7.1 listening')
