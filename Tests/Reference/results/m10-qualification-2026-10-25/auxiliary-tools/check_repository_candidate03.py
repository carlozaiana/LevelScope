#!/usr/bin/env python3
import datetime
import hashlib
import json
import pathlib
import subprocess

r = pathlib.Path(__file__).resolve().parent
w = r.parents[1]
before = json.loads((r / 'preservation-before.json').read_text())
identity = json.loads((r / 'candidate-03/identity.json').read_text())
def git(*args):
    return subprocess.check_output(['git', '-C', str(w), *args], text=True)
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
source_errors = [name for name, digest in identity['source_files_sha256'].items() if sha(w / name) != digest]
tags = git('for-each-ref', '--format=%(refname) %(objectname) %(*objectname)', 'refs/tags').splitlines()
changed_tags = [line for line in before['tags'] if line not in tags]
head = git('rev-parse', 'HEAD').strip()
unrelated_errors = [name for name, digest in before['unrelated'].items() if sha(w / name) != digest]
result = dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    current_source_digest=identity['source_digest'], source_manifest_matches=not source_errors,
    source_files_checked=len(identity['source_files_sha256']), source_errors=source_errors,
    existing_tags_unchanged=not changed_tags, existing_tag_count=len(before['tags']),
    changed_existing_tags=changed_tags, extra_tags=[line for line in tags if line not in before['tags']],
    worktrees_unchanged=git('worktree', 'list', '--porcelain') == before['worktrees'],
    head_unchanged=head == before['head'], head=head,
    branch=git('branch', '--show-current').strip(),
    unrelated_unchanged=not unrelated_errors, unrelated_errors=unrelated_errors)
(r / 'repository-preservation-final-candidate-03.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
raise SystemExit(0 if all(result[k] for k in ('source_manifest_matches','existing_tags_unchanged','worktrees_unchanged','head_unchanged','unrelated_unchanged')) else 1)
