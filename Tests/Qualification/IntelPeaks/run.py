"""Explicit, sequential scheduling-change qualification. Never invokes historical runners."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess
import time


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def signature(d):
    outcomes = {name: {k: v for k, v in row.items()
                       if k not in ('main_preflight_ms', 'main_accept_ms')}
                for name, row in d['outcomes'].items()}
    return [d['identities'], outcomes, d['file_sha256'], d['source_lufs'],
            d['current_lufs'], d['base_build']]


def main():
    p = argparse.ArgumentParser()
    p.add_argument('candidate', type=Path)
    p.add_argument('output', type=Path)
    a = p.parse_args()
    root, out = a.candidate.resolve(), a.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    identity = json.loads((root / 'identity.json').read_text())
    arch = platform.machine()
    rows, failure, comparisons = [], '', {}

    def run(name, args, report=None, text_pass=None):
        started = time.monotonic()
        command = ['/usr/bin/arch', '-' + arch, *map(str, args)]
        print('START', name, flush=True)
        with (out / (name + '.log')).open('w') as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        row = dict(name=name, command=command, exit_code=result.returncode,
                   seconds=time.monotonic() - started)
        rows.append(row)
        (out / 'commands.json').write_text(json.dumps(rows, indent=2) + '\n')
        assert result.returncode == 0, name + ': nonzero exit'
        data = None
        if report:
            data = json.loads(report.read_text())
            row['report_sha256'] = sha(report)
            assert data.get('passed', True) is True, name + ': failed report'
            assert data.get('build', data.get('base_build')) == identity['measurement_build']
        if text_pass:
            assert text_pass in (out / (name + '.log')).read_text()
        row['passed'] = True
        (out / 'commands.json').write_text(json.dumps(rows, indent=2) + '\n')
        print('PASS', name, round(row['seconds'], 2), flush=True)
        return data

    try:
        assert arch in ('arm64', 'x86_64') and platform.system() == 'Darwin'
        if arch == 'x86_64':
            assert platform.mac_ver()[0] == '10.13.6', 'Intel must run on actual High Sierra'
        for path, digest in identity['artifacts_sha256'].items():
            assert sha(root / path) == digest, path
        run('core-parity', [root / 'bin/LevelScopePeakChannelParity'], text_pass='PASS: all result bits')
        run('measurement', [root / 'bin/LevelScopeMeasurementTests', '--report', out / 'measurement.json'], out / 'measurement.json')
        run('section-smoke', [root / 'bin/LevelScopeSectionProposalTests', '--smoke', '--root', out / 'section-files',
                              '--output', out / 'section-smoke.json'], out / 'section-smoke.json')
        for channels in (2, 6):
            fixture = root / 'fixtures' / ('programme-%d.wav' % channels)
            expected = None
            timing = {mode: [] for mode in ('serial', 'app-mode')}
            # One decision/untimed-warmup comparison followed by three alternating
            # ordinary repetitions. No repeated runs to hide a failing result.
            for trial in range(4):
                for mode in (('serial', 'app-mode') if trial % 2 == 0 else ('app-mode', 'serial')):
                    name = '%d-%s-%d' % (channels, mode, trial)
                    report = out / (name + '.json')
                    data = run(name, [root / 'bin/peaks' / mode, fixture, report,
                                      *(['--decisions'] if trial == 0 else [])], report)
                    if trial == 0:
                        assert len(data['outcomes']) == 7
                        if expected is None:
                            expected = signature(data)
                        assert signature(data) == expected, 'exact decision parity'
                    else:
                        timing[mode].append(data)
                        # The five regular operation identities are compared
                        # between repeats with identical runtime/job sequence.
                        assert signature(data) == signature(timing['serial'][0] if timing['serial'] else data)
            assert all(signature(d) == signature(timing['serial'][0])
                       for ds in timing.values() for d in ds), 'exact repeated-workflow parity'
            operations = {}
            for op in ('source', 'fresh', 'reused', 'current'):
                values = {mode: [d['operations'][op]['wall_ms'] for d in ds]
                          for mode, ds in timing.items()}
                medians = {mode: statistics.median(v) for mode, v in values.items()}
                operations[op] = dict(milliseconds=values, median_ms=medians,
                    speedup=medians['serial'] / medians['app-mode'])
                assert medians['app-mode'] <= .90 * medians['serial'], 'less than 10% acceleration: ' + op
            ui = {}
            for mode, ds in timing.items():
                gaps = sorted(v for d in ds for op in ('source', 'fresh', 'reused', 'current')
                              for v in d['operations'][op]['heartbeat_gaps_ms'])
                ui[mode] = dict(p95_ms=gaps[math.ceil(.95 * len(gaps)) - 1], maximum_ms=max(gaps))
            assert ui['app-mode']['p95_ms'] <= ui['serial']['p95_ms'] + 20
            assert ui['app-mode']['maximum_ms'] <= 200
            cores = [d['operations']['source']['process_cpu_ms'] / d['operations']['source']['wall_ms']
                     for d in timing['app-mode']]
            assert statistics.median(cores) > 1.15, 'no observed service concurrency'
            comparisons[str(channels)] = dict(exact=True, operations=operations, heartbeat=ui,
                                              source_process_cores=cores)
            run('native-%d' % channels, [root / 'bin/native', fixture, out / ('native-%d' % channels)],
                out / ('native-%d' % channels) / 'result.json')
        assert len(rows) == 21 and all(r.get('passed') for r in rows), 'required selection incomplete'
    except Exception as error:
        failure = repr(error)
    summary = dict(passed=not failure, failure=failure, comparisons=comparisons,
                   architecture=arch, os=platform.mac_ver()[0],
                   os_build=subprocess.check_output(['sw_vers', '-buildVersion'], text=True).strip(),
                   cpu=subprocess.check_output(['sysctl', '-n', 'machdep.cpu.brand_string'], text=True).strip(),
                   measurement_build=identity['measurement_build'], identity_sha256=sha(root / 'identity.json'),
                   commands=rows, physical_devices_opened=False, video_loaded=False,
                   ordinary_app_stack_check='Separate required receipt on Golem')
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('PASS' if not failure else failure, flush=True)
    return 1 if failure else 0


if __name__ == '__main__':
    raise SystemExit(main())
