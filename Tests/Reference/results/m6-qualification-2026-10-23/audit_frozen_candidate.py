#!/usr/bin/env python3
"""Read-only source/package/Git preservation audit after qualification."""
import argparse
import hashlib
import json
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--repository', type=Path, required=True)
p.add_argument('--candidate', type=Path, required=True)
p.add_argument('--baseline', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
identity = json.loads((a.candidate / 'identity.json').read_text())
baseline = json.loads(a.baseline.read_text())

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def mismatches(root, manifest):
    return [name for name, expected in manifest.items()
            if not (root / name).is_file() or sha(root / name) != expected]

def git(*args):
    return subprocess.check_output(['git', '-C', str(a.repository), *args], text=True).strip()

tags = git('for-each-ref', '--format=%(refname:short) %(objectname) %(*objectname)', 'refs/tags')
worktrees = git('worktree', 'list', '--porcelain')
result = dict(source_digest=identity['source_digest'],
              identity_sha256=sha(a.candidate / 'identity.json'),
              source_files=len(identity['source_files_sha256']),
              source_mismatches=mismatches(a.repository, identity['source_files_sha256']),
              artifacts=len(identity['artifacts_sha256']),
              artifact_mismatches=mismatches(a.candidate, identity['artifacts_sha256']),
              head=git('rev-parse', 'HEAD'), main=git('rev-parse', 'refs/heads/main'),
              head_unchanged=git('rev-parse', 'HEAD') == baseline['head'],
              tags_unchanged=tags == baseline['tags'],
              worktrees_unchanged=worktrees == baseline['worktrees'],
              tags=tags, worktrees=worktrees,
              publication='No M6 commit, tag, push or merge performed in this implementation task')
result['passed'] = (not result['source_mismatches'] and not result['artifact_mismatches']
                    and result['head_unchanged'] and result['tags_unchanged'] and result['worktrees_unchanged'])
a.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k not in ('tags', 'worktrees')}, indent=2))
raise SystemExit(0 if result['passed'] else 1)
