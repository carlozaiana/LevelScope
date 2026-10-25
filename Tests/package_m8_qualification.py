#!/usr/bin/env python3
"""Freeze M8 using the pinned M7 regression media; never modify the input bundle."""
import json
import pathlib
import shutil
import sys

import package_m7_qualification
from collect_m1_environment import sha256


def main():
    # Keep the established package's explicit media, source and build provenance.
    package_m7_qualification.main()
    output = pathlib.Path(sys.argv[sys.argv.index('--output') + 1])
    build = pathlib.Path(sys.argv[sys.argv.index('--build') + 1])
    for target in ('LevelScopeSessionTests', 'LevelScopeSessionPlaybackTests'):
        shutil.copy2(build / (target + '_artefacts/Release') / target, output / 'bin' / target)
    for name in ('run_m8_qualification.py', 'compare_m8_reports.py', 'M8QualificationRunnerTests.py'):
        shutil.copy2(pathlib.Path('Tests') / name, output / 'Tests' / name)
    path = output / 'identity.json'
    identity = json.loads(path.read_text())
    identity['milestone'] = 'M8 reproducible sessions'
    identity['artifacts_sha256'] = {str(p.relative_to(output)): sha256(p)
        for p in sorted(output.rglob('*')) if p.is_file() and p != path}
    path.write_text(json.dumps(identity, indent=2) + '\n')
    print('Final M8 identity SHA256:', sha256(path))


if __name__ == '__main__':
    main()
