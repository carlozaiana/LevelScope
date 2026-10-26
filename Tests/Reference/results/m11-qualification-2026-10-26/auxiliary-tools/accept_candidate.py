#!/usr/bin/env python3
"""Post-freeze evidence gate; never builds, renders, opens devices or publishes Git."""
import argparse
import hashlib
import json
import pathlib
import subprocess
import sys


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--workspace', type=pathlib.Path, required=True)
    parser.add_argument('--work', type=pathlib.Path, required=True)
    parser.add_argument('--output', type=pathlib.Path, required=True)
    args = parser.parse_args()
    candidate = args.work / 'candidate-02'
    identity = json.loads((candidate / 'identity.json').read_text())
    for name, expected in identity['source_files_sha256'].items():
        assert digest(args.workspace / name) == expected, 'Workspace source changed: ' + name
    for name, expected in identity['artifacts_sha256'].items():
        assert digest(candidate / name) == expected, 'Frozen package changed: ' + name
    results = args.work / 'candidate-02-results'
    comparison = results / 'complete-cross-platform.json'
    subprocess.run([sys.executable, '-B', str(candidate / 'Tests/compare_m11_reports.py'),
                    str(results / 'sequoia'), str(results / 'high-sierra'),
                    '--output', str(comparison), '--complete'], check=True)
    cross = json.loads(comparison.read_text())
    assert cross['passed'] is True and cross['qualification_complete'] is True
    platforms = {}
    for platform in ('sequoia', 'high-sierra'):
        root = results / platform
        phases = {}
        for phase in ('short', 'long', 'exchange'):
            summary = root / phase / 'summary.json'
            assert json.loads(summary.read_text())['automated_passed'] is True
            phases[phase] = dict(summary_sha256=digest(summary), runs_sha256=digest(root / phase / 'runs.json'))
        independent = root / 'long/independent-rf64.json'
        decoded = json.loads(independent.read_text())
        assert decoded['passed'] is True and decoded['independent_sum_and_single_final_dither'] is True
        assert decoded['composition_oracle_codes'] == 98304
        long_result = json.loads((root / 'long/long.json').read_text())['long_composition']
        assert decoded['artifact_sha256'] == long_result['artifact_sha256']
        platforms[platform] = dict(phases=phases, environment=json.loads((root / 'long/environment.json').read_text()),
                                   ui_sha256=digest(root / 'ui.json'), long_workflow=long_result,
                                   independent_rf64_sha256=digest(independent))
    receipt = dict(candidate='candidate-02', qualified=True, ready_for_checkpoint_review=True,
                   source_digest=identity['source_digest'], base_commit=identity['base_commit'],
                   measurement_build=identity['measurement_build'],
                   app_sha256=identity['artifacts_sha256']['LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone'],
                   identity_sha256=digest(candidate / 'identity.json'), platforms=platforms,
                   complete_comparison_sha256=digest(comparison),
                   evidence_policy='Explicit M11 impact plan, ProjectIndex 42; scoped M10 reuse recorded separately',
                   required_checks_outstanding=[], physical_devices_opened=False,
                   physical_multichannel_listening='Separately pending; not part of approved offline/injected M11 scope',
                   accepted_limitations=['Occasional picture skips; no acceptance of audio interruption, drift or sustained freezes',
                                         'Unity component gains and sample-zero alignment only; no dialogue analysis or component DSP'],
                   checkpoint_published=False, main_merged=False)
    assert not args.output.exists(), 'Never replace an existing acceptance receipt'
    args.output.write_text(json.dumps(receipt, indent=2) + '\n')
    print('All required M11 gates passed; ready for user checkpoint review. No Git publication performed.')


if __name__ == '__main__':
    main()
