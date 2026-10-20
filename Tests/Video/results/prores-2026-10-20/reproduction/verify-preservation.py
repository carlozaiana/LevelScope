from pathlib import Path
import hashlib,json,subprocess,shlex,argparse

parser=argparse.ArgumentParser()
parser.add_argument("--checkpoint-tag",required=True)
args=parser.parse_args()
root=Path("/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation")
original=Path("/Users/fuzz/Developer/LevelScope")
build=root/"build/prores-20261020"
before=json.loads((build/"evidence/preservation-before.json").read_text())
def git(*args,cwd=root):
    return subprocess.check_output(["git",*args],cwd=cwd,text=True).strip()
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(4*1024*1024),b""):h.update(block)
    return h.hexdigest()

tag_records=[]
remote={line.split()[1]:line.split()[0] for line in git("ls-remote","origin","refs/heads/main","refs/heads/codex/video-implementation","refs/tags/*").splitlines()}
for line in before["tags"].splitlines():
    ref,obj,commit=line.split()
    actual_obj=git("rev-parse",ref)
    actual_commit=git("rev-parse",ref+"^{commit}")
    matches=actual_obj==obj==remote.get(ref) and actual_commit==commit==remote.get(ref+"^{}")
    tag_records.append(dict(ref=ref,tag_object=obj,commit=commit,unchanged_local_and_remote=matches))
assert all(t["unchanged_local_and_remote"] for t in tag_records)
checkpoint_ref="refs/tags/"+args.checkpoint_tag
checkpoint_obj=git("rev-parse",checkpoint_ref)
checkpoint_commit=git("rev-parse",checkpoint_ref+"^{commit}")
assert remote.get(checkpoint_ref)==checkpoint_obj and remote.get(checkpoint_ref+"^{}")==checkpoint_commit
assert checkpoint_commit==git("rev-parse","HEAD")==remote["refs/heads/codex/video-implementation"]

inventory_changes=[]
for name,digest in before["original"].items():
    p=original/name
    current=sha(p) if p.is_file() else None
    if current!=digest:
        inventory_changes.append(dict(path=name,before=digest,after=current))
actual_movie=sha(original/"Tests/51_SyncTest.mov")
assert actual_movie==before["actual_fixture_sha256"]
candidate=json.loads((build/"candidate-04/identity.json").read_text())
source_mismatches=[name for name,digest in candidate["source_files_sha256"].items() if sha(root/name)!=digest]
assert not source_mismatches
artifact_mismatches=[name for name,digest in candidate["artifacts_sha256"].items() if sha(build/"candidate-04"/name)!=digest]
assert not artifact_mismatches
remote_code = """from pathlib import Path
import json,hashlib
p=Path('/Users/fuzz/LevelScope-ProRes-qualification-2026-10-20/candidate-04')
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as stream:
  for b in iter(lambda:stream.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
i=json.loads((p/'identity.json').read_text())
bad=[n for n,h in i['artifacts_sha256'].items() if sha(p/n)!=h]
print(json.dumps(dict(verified_artifact_count=len(i['artifacts_sha256']),mismatches=bad,identity_sha256=sha(p/'identity.json'))))
"""
remote_candidate=json.loads(subprocess.check_output(['ssh','-S','/tmp/levelscope-prores-golem','-o','BatchMode=yes','fuzz@192.168.2.2','/usr/local/bin/python3 -c '+shlex.quote(remote_code)],text=True))
assert not remote_candidate['mismatches']
assert remote_candidate['identity_sha256']==sha(build/'candidate-04/identity.json')
report=dict(
    checkpoint_tags=tag_records,
    published_checkpoint=dict(ref=checkpoint_ref,tag_object=checkpoint_obj,commit=checkpoint_commit,verified_local_and_remote=True),
    original_head_before=before["original_head"],original_head_now=git("rev-parse","HEAD",cwd=original),
    original_branch_now=git("branch","--show-current",cwd=original),
    original_inventory_changes=inventory_changes,
    original_tracked_diff=git("diff","--stat",cwd=original),
    original_index_diff=git("diff","--cached","--stat",cwd=original),
    actual_movie_sha256=actual_movie,
    tested_source_matches_candidate=True,
    local_verified_artifact_count=len(candidate["artifacts_sha256"]),
    golem_final_candidate_verification=remote_candidate,
    video_branch=git("branch","--show-current"),video_head=git("rev-parse","HEAD"),
    remote_main=remote["refs/heads/main"],remote_video=remote["refs/heads/codex/video-implementation"],
    worktrees=git("worktree","list","--porcelain"),
    preservation_note="No metadata differences are reverted. No installed products, Pro40 settings or historical qualification bundles were modified by this task.",
)
out=build/"evidence/preservation-after.json"
assert not out.exists()
out.write_text(json.dumps(report,indent=2)+"\n")
print("Checkpoint tags preserved; exact source and actual-movie hashes match.")
print("Original inventory differences:",[r["path"] for r in inventory_changes])
