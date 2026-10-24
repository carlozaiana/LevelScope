#!/usr/bin/env python3
"""Read completed immutable comparison outputs; never launch a test or alter input."""
import argparse
import bisect
import hashlib
import gzip
import json
import pathlib
import statistics


def read_bytes(path):
    return path.read_bytes() if path.exists() else gzip.decompress(path.with_name(path.name + '.gz').read_bytes())


def sha(path):
    return hashlib.sha256(read_bytes(path)).hexdigest()


def timing(events):
    gaps = [(b[0] - a[0]) / 1e6 for a, b in zip(events, events[1:])]
    return {'events': len(events), 'maximum_interval_ms': max(gaps, default=0),
            'mean_interval_ms': statistics.mean(gaps) if gaps else 0}


def trace_analysis(path, report):
    trace = json.loads(read_bytes(path))
    assert trace['valid'] and not trace['ui_dropped'] and not trace['display_dropped']
    ui, ticks = trace['ui_events'], trace['display_events']
    spans, active = [], []
    for event in ui:
        if event[1] in (1, 5):
            active.append({'begin': event, 'selections': []})
        elif event[1] == 3 and active:
            active[-1]['selections'].append(event)
        elif event[1] in (2, 6):
            begin_kind = 1 if event[1] == 2 else 5
            assert active and active[-1]['begin'][1] == begin_kind, 'Unexpected nested callback structure'
            span = active.pop(); span['end'] = event
            span['duration_ms'] = (event[0] - span['begin'][0]) / 1e6
            spans.append(span)
    # finish() stops the trace from inside the final timer, intentionally before
    # serialization. It is not a completed full-timer duration sample.
    assert all(s['begin'][1] == 1 for s in active)
    def qualified(span):
        return span['selections'] and span['selections'][-1][3] / span['selections'][-1][4] >= 2
    paints = [s for s in spans if s['begin'][1] == 5 and qualified(s)]
    timers = [s for s in spans if s['begin'][1] == 1 and qualified(s)]
    assert paints
    lo, hi = paints[0]['begin'][0], paints[-1]['begin'][0]
    in_window = lambda events: [e for e in events if lo <= e[0] <= hi]
    qualified_ui = in_window(ui)
    by_kind = {kind: [e for e in qualified_ui if e[1] == kind] for kind in [1, 3, 4, 5, 6, 7, 9]}
    native = in_window(ticks)
    qualified_selections = [e for e in ui if e[1] == 3 and lo <= e[0] <= paints[-1]['end'][0]]
    gaps = [(a, b, (b['begin'][0] - a['begin'][0]) / 1e6) for a, b in zip(paints, paints[1:])]
    largest = []
    for a, b, gap in sorted(gaps, key=lambda item: item[2], reverse=True)[:5]:
        start, stop = a['begin'][0], b['begin'][0]
        between = [e for e in ui if start < e[0] < stop]
        selected = [e for e in between if e[1] == 3]
        prev = a['selections'][-1]; current = b['selections'][-1]
        largest.append({
            'gap_ms': gap, 'previous_frame': prev[8], 'frame': current[8],
            'audio_seconds': current[3] / current[4],
            'audio_elapsed_ms': 1000 * (current[3] / current[4] - prev[3] / prev[4]),
            'timer_entries_between': sum(e[1] == 1 for e in between),
            'repaint_requests_between': sum(e[1] == 4 for e in between),
            'main_vblanks_between': sum(e[1] == 7 for e in between),
            'independent_ticks_between': sum(start < e[0] < stop for e in ticks),
            'selected_frames_between': sorted(set(e[8] for e in selected)),
            'missing_picture_selections_between': sum(e[8] < 0 for e in selected),
        })
    observed_skips = sum(max(0, b['selections'][-1][8] - a['selections'][-1][8] - 1)
                         for a, b in zip(paints, paints[1:]))
    result = {
        'scope': 'Derived callback observations within the >=2 second playback paint span; not physical presentation',
        'trace_sha256': sha(path), 'ui_events_total': len(ui), 'independent_ticks_total': len(ticks),
        'trace_loss': {'ui': trace['ui_dropped'], 'display': trace['display_dropped']},
        'nominal_hz': trace['nominal_hz'], 'display_id': trace['display_id'],
        'qualified_paints_derived': len(paints), 'qualified_paints_reported': report['paints'],
        'paint_count_matches_report': len(paints) == report['paints'],
        'paint_skips_derived': observed_skips,
        'paint_skips_match_report': observed_skips == report['actual_paint_skipped_frames'],
        'maximum_full_paint_ms': max(s['duration_ms'] for s in paints),
        'maximum_full_timer_ms': max((s['duration_ms'] for s in timers), default=0),
        'paint': timing([s['begin'] for s in paints]),
        'timer': timing(by_kind[1]), 'main_vblank': timing(by_kind[7]), 'independent_display': timing(native),
        'transport_generations': sorted(set(e[5] for e in qualified_selections)),
        'view_epochs': sorted(set(e[6] for e in qualified_selections)),
        'window_flags_at_paint': sorted(set(s['begin'][2] for s in paints)),
        'display_ids_at_paint': sorted(set(s['begin'][12] for s in paints)),
        'missing_picture_selections': sum(e[8] < 0 for e in qualified_selections),
        'selection_mismatches': sum(e[8] >= 0 and e[7] != e[8] for e in qualified_selections),
        'maximum_queue_count': max(e[11] for e in qualified_selections),
        'largest_paint_gaps': largest,
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('comparison', type=pathlib.Path)
    parser.add_argument('output', type=pathlib.Path)
    args = parser.parse_args()
    root = args.comparison.resolve()
    summary = json.loads((root / 'summary.json').read_text())
    runs = []
    keys = ['passed', 'device', 'sample_rate', 'all_routes_muted', 'source_seconds', 'paints',
            'actual_paint_skipped_frames', 'blank_paints', 'bound_failures', 'paint_deadline_failures',
            'underrun_blocks', 'maximum_paint_interval_ms', 'maximum_cpu_paint_ms', 'hal_output_timestamp',
            'memory_passed', 'memory_growth_assessed', 'peak_sampled_playback_rss_bytes',
            'audio_sha256', 'video_sha256']
    for outcome in summary['runs']:
        folder = root / ('%02d-%s' % (outcome['number'], outcome['role']))
        report_path = folder / 'report.json'; report = json.loads(report_path.read_text())
        assert sha(report_path) == outcome['report_sha256']
        run = {'number': outcome['number'], 'role': outcome['role'], 'report_sha256': sha(report_path),
               **{k: report[k] for k in keys}}
        if outcome['role'] == 'diagnostic':
            assert sha(folder / 'trace.json') == outcome['trace']['trace_sha256']
            run['trace'] = trace_analysis(folder / 'trace.json', report)
        runs.append(run)
    data = {'scope': 'Read-only derived comparison analysis; no additional test runs',
            'comparison_path': str(root), 'summary_sha256': sha(root / 'summary.json'),
            'candidate_identity_sha256': summary['identity_sha256'],
            'audio_configuration_unchanged': summary['audio_configuration_unchanged'], 'runs': runs}
    with args.output.open('x') as f: json.dump(data, f, indent=2); f.write('\n')
    print(json.dumps([{'run': r['number'], 'role': r['role'], 'passed': r['passed'],
                      'skips': r['actual_paint_skipped_frames'],
                      'maximum_qualified_paint_gap_ms': r['maximum_paint_interval_ms'],
                      'trace': {k: v for k, v in r.get('trace', {}).items() if k not in ['largest_paint_gaps']}}
                     for r in runs], indent=2))


if __name__ == '__main__':
    main()
