#!/usr/bin/env python3
"""Freeze M9 over the preserved M8/M7 regression assets; no installation."""
import argparse
import json
import pathlib
import shutil
import sys
import package_m8_qualification
from collect_m1_environment import sha256


def main():
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--m9-fixtures', required=True, type=pathlib.Path)
    extras, remaining = parser.parse_known_args()
    sys.argv = [sys.argv[0]] + remaining
    package_m8_qualification.main()
    output = pathlib.Path(remaining[remaining.index('--output') + 1])
    build = pathlib.Path(remaining[remaining.index('--build') + 1])
    target = 'LevelScopeProgrammeTests'
    shutil.copy2(build / (target + '_artefacts/Release') / target, output / 'bin' / target)
    shutil.copytree(extras.m9_fixtures, output / 'm9-fixtures', ignore=shutil.ignore_patterns('.DS_Store'))
    for name in ('run_m9_qualification.py', 'compare_m9_reports.py', 'compare_m9_extended_oracles.py', 'M9QualificationRunnerTests.py', 'run_m9_muted_video.py'):
        shutil.copy2(pathlib.Path('Tests') / name, output / 'Tests' / name)
    path = output / 'identity.json'
    identity = json.loads(path.read_text())
    identity['milestone'] = 'M9 separate-file stereo/5.1 programme import'
    identity['artifacts_sha256'] = {str(p.relative_to(output)): sha256(p)
        for p in sorted(output.rglob('*')) if p.is_file() and p != path}
    path.write_text(json.dumps(identity, indent=2) + '\n')
    print('Final M9 identity SHA256:', sha256(path))


if __name__ == '__main__':
    main()
