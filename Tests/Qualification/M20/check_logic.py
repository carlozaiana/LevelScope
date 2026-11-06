"""Explicit Logic handover acceptance; never launches historical runners or devices."""
import argparse
import hashlib
import json
from pathlib import Path

REQUIRED = (
    'core', 'stereo', 'surround', 'heterogeneous', 'failures', 'audio', 'population',
    'gating', 'cinema', 'legacy-core', 'legacy-stereo', 'legacy-surround',
    'legacy-inputs', 'legacy-failures', 'legacy-gating', 'requirements', 'section-smoke',
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('receipts', type=Path)
    parser.add_argument('identity', type=Path)
    parser.add_argument('--repository', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    identity = json.loads(args.identity.read_text())
    errors, results = [], {}
    for name in REQUIRED:
        path = args.receipts / (name + '.json')
        if not path.is_file():
            errors.append(name + ': missing required receipt')
            continue
        value = json.loads(path.read_text())
        if value.get('passed') is not True:
            errors.append(name + ': did not pass')
        if value.get('build') != identity['measurement_build']:
            errors.append(name + ': belongs to a different product source identity')
        if any(value.get(key) for key in ('devices_opened', 'physical_devices_opened', 'video_loaded')):
            errors.append(name + ': exceeded silent Logic-stage scope')
        results[name] = {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                         'checks': len(value.get('checks', [])), 'cases': value.get('cases', [])}
    for name, expected in identity['source_files_sha256'].items():
        file = args.repository / name
        if not file.is_file() or hashlib.sha256(file.read_bytes()).hexdigest() != expected:
            errors.append('source changed: ' + name)
    summary = dict(logic_handover_checks_passed=not errors, m20_combined_qualified=False,
                   errors=errors, required=list(REQUIRED), results=results,
                   measurement_build=identity['measurement_build'],
                   source_manifest_sha256=identity['source_manifest_sha256'],
                   pending=['UI implementation/integration and native responsiveness',
                            'actual Intel High Sierra and cross-platform checks',
                            'independent real holdouts/control variants and coordinated listening'],
                   physical_devices_opened=False, video_loaded=False)
    args.output.write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps({'logic_handover_checks_passed': not errors, 'errors': errors}))
    return 0 if not errors else 1


if __name__ == '__main__':
    raise SystemExit(main())
