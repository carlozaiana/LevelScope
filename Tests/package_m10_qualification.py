#!/usr/bin/env python3
"""Freeze an M10 candidate over pinned M9 regression assets, without installation."""
import json, pathlib, shutil, sys
import package_m9_qualification
from collect_m1_environment import sha256

def main():
    package_m9_qualification.main()
    output=pathlib.Path(sys.argv[sys.argv.index('--output')+1]);build=pathlib.Path(sys.argv[sys.argv.index('--build')+1])
    target='LevelScopeExportTests';shutil.copy2(build/(target+'_artefacts/Release')/target,output/'bin'/target)
    for name in ('run_m10_qualification.py','verify_m10_pcm.py','verify_m10_decoder.py','compare_m10_reports.py','M10QualificationRunnerTests.py','run_m10_muted_video.py'):
        shutil.copy2(pathlib.Path('Tests')/name,output/'Tests'/name)
    path=output/'identity.json';identity=json.loads(path.read_text());identity['milestone']='M10 verified PCM export and reporting'
    identity['artifacts_sha256']={str(p.relative_to(output)):sha256(p) for p in sorted(output.rglob('*')) if p.is_file() and p!=path}
    path.write_text(json.dumps(identity,indent=2)+'\n');print('Final M10 identity SHA256:',sha256(path))

if __name__=='__main__':main()
