import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path('/Users/fuzz/.codex/worktrees/a411/LevelScope')
WORK = ROOT / 'build/boundaries-workflow-publication-2026-11-23'
INTEGRATION = ROOT / 'build/m18-publication-2026-11-03/integration'
BRANCH = 'codex/propose-verify-workflow-2026-11-22'
TAG = 'm20-boundary-editing-propose-verify-2026-11-23'
CHECKPOINT = 'd5f6e86801c675794bb0b37281afc44bd30aa961'
TESTED = '79b55eb177a01e56fdbe86f843e03140df356bcd'
MAIN_BEFORE = '1056ad4ace536abfd68904f083c7fc56a1a8999f'
QUAL = ROOT / 'Tests/Qualification/ProposalWorkflow/results/2026-11-23'

def git(*args, cwd=ROOT):
    return subprocess.check_output(['git', *args], cwd=cwd, text=True).strip()

def read(path):
    return json.loads(path.read_text())

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')

def refs(path):
    return {line.split('\t')[1]: line.split('\t')[0]
            for line in path.read_text().splitlines()
            if '\t' in line and not line.startswith('ref:')}

def preflight():
    previous = read(ROOT / 'Tests/Qualification/ProposalWorkflow/results/2026-11-22/preservation-before.json')
    assert read(QUAL / 'acceptance.json')['focused_actual_platform_qualification_passed']
    protected = dict(previous['protected'])
    for name, digest in protected.items():
        assert sha(Path(name)) == digest, name
    definitions = previous['packages'] + [
        {'path': str(ROOT / 'build/applied-region-boundaries-2026-11-22/candidate-01'),
         'identity_sha256': '64d5f870f5382f332f58564b945f0ebf61012218570dc1bb970bb06d65171942'},
        {'path': str(ROOT / 'build/propose-verify-workflow-2026-11-22/candidate-01'),
         'identity_sha256': '4eb921611c8cec84e98b79ef7cb94aa9fd09da6b27bfbbeed958e35706e68ce3'}]
    packages = []
    for package in definitions:
        path = Path(package['path'])
        assert sha(path / 'identity.json') == package['identity_sha256']
        manifest = read(path / 'identity.json')
        for name, digest in manifest['artifacts_sha256'].items():
            assert sha(path / name) == digest, str(path / name)
            protected[str(path / name)] = digest
        for name in ['identity.json', 'identity.sha256']:
            if (path / name).exists():
                protected[str(path / name)] = sha(path / name)
        packages.append(dict(package, indexed_artifacts=len(manifest['artifacts_sha256'])))
    for folder in [ROOT / 'Tests/Qualification/ProposalWorkflow/results', ROOT / 'Tests/Qualification/AppliedBoundaries/results']:
        for path in folder.rglob('*'):
            if path.is_file():
                protected[str(path)] = sha(path)
    remote = refs(WORK / 'remote-before.txt')
    tags = [line.split(' ') for line in git('for-each-ref', '--format=%(refname) %(objectname) %(*objectname)', 'refs/tags').splitlines()]
    old_tags = [line.split(' ') for line in previous['tags'].splitlines()]
    assert all(row in tags for row in old_tags)
    for ref, obj, peeled in tags:
        assert remote[ref] == obj and remote.get(ref + '^{}') == peeled
    assert remote['refs/heads/main'] == git('rev-parse', 'main') == git('rev-parse', 'origin/main') == MAIN_BEFORE
    assert git('branch', '--show-current') == BRANCH and git('rev-parse', 'HEAD') == CHECKPOINT
    assert 'refs/heads/' + BRANCH not in remote and 'refs/tags/' + TAG not in remote
    assert not git('diff', 'HEAD', '--')
    subprocess.run(['git', 'merge-base', '--is-ancestor', MAIN_BEFORE, CHECKPOINT], cwd=ROOT, check=True)
    assert git('rev-parse', TESTED + ':Source') == git('rev-parse', CHECKPOINT + ':Source') == '0a96e9ab0b47847d301c7b273577b945a4c1921a'
    assert not git('diff', TESTED, CHECKPOINT, '--', 'Source', 'CMakeLists.txt')
    manifest = read(Path(definitions[-1]['path']) / 'identity.json')
    for name, digest in manifest['source_sha256'].items():
        assert hashlib.sha256(subprocess.check_output(['git', 'show', CHECKPOINT + ':' + name], cwd=ROOT)).hexdigest() == digest, name
    worktrees = []
    for block in git('worktree', 'list', '--porcelain').split('\n\n'):
        path = Path(block.splitlines()[0][9:])
        untracked = {}
        for name in subprocess.check_output(['git', 'ls-files', '--others', '--exclude-standard', '-z'], cwd=path).decode().split('\0'):
            if name and (path / name).is_file():
                untracked[name] = sha(path / name)
        worktrees.append({'path': str(path), 'head': git('rev-parse', 'HEAD', cwd=path),
                          'branch': git('branch', '--show-current', cwd=path),
                          'tracked_diff_sha256': hashlib.sha256(git('diff', 'HEAD', '--', cwd=path).encode()).hexdigest(),
                          'untracked': untracked})
    assert not git('diff', 'HEAD', '--', cwd=INTEGRATION)
    result = {'authorization': {'source_chat': '01a13c8f-90a9-7150-813f-cb057626ef7e',
              'text': 'Please tell the UI chat to publish the combined boundary-editing and workflow checkpoint to GitHub and merge it into main, preserving all existing tags, tested builds, evidence and the agreed qualification deferrals.'},
              'repository': git('remote', 'get-url', 'origin'), 'branch': BRANCH, 'checkpoint': CHECKPOINT,
              'tested': TESTED, 'source_tree': git('rev-parse', CHECKPOINT + ':Source'), 'tag': TAG,
              'main_before': MAIN_BEFORE, 'remote_before': remote, 'local_tags_before': tags,
              'packages': packages, 'protected': protected, 'worktrees': worktrees,
              'candidate_source_entries_verified': len(manifest['source_sha256']),
              'qualification_receipt_sha256': sha(QUAL / 'acceptance.json'), 'passed': True}
    save(WORK / 'preflight.json', result)
    print('Preflight PASS:', len(packages), 'packages;', len(protected), 'protected files;', len(tags), 'existing tags;', len(manifest['source_sha256']), 'source entries; fast-forward available.')

def verify(label, expected):
    before = read(WORK / 'preflight.json')
    remote = refs(WORK / ('remote-' + label + '.txt'))
    checks = {}
    def check(name, ok):
        checks[name] = bool(ok)
        assert ok, name
    for name, digest in before['protected'].items():
        check('protected:' + name, sha(Path(name)) == digest)
    tags = {row.split(' ')[0]: row.split(' ')[1:] for row in git('for-each-ref', '--format=%(refname) %(objectname) %(*objectname)', 'refs/tags').splitlines()}
    for ref, obj, peeled in before['local_tags_before']:
        check('old_local_tag:' + ref, tags.get(ref) == [obj, peeled])
        check('old_remote_tag:' + ref, remote.get(ref) == obj and remote.get(ref + '^{}') == peeled)
    new_ref = 'refs/tags/' + TAG
    tag_object = git('rev-parse', new_ref)
    check('annotated_tag', git('cat-file', '-t', tag_object) == 'tag')
    check('tag_identity', git('rev-parse', new_ref + '^{}') == CHECKPOINT and remote.get(new_ref) == tag_object and remote.get(new_ref + '^{}') == CHECKPOINT)
    check('tag_inventory', set(tags) == {x[0] for x in before['local_tags_before']} | {new_ref})
    for ref, value in before['remote_before'].items():
        if ref not in ['HEAD', 'refs/heads/main']:
            check('unrelated_remote_ref:' + ref, remote.get(ref) == value)
    check('remote_inventory', set(remote) == set(before['remote_before']) | {'refs/heads/' + BRANCH, new_ref, new_ref + '^{}'})
    for ref in ['HEAD', 'refs/heads/main', 'refs/heads/' + BRANCH]:
        check('published:' + ref, remote.get(ref) == expected)
    check('local_heads', git('rev-parse', 'main') == git('rev-parse', BRANCH) == git('rev-parse', 'origin/main') == expected)
    check('tested_source', git('rev-parse', expected + ':Source') == before['source_tree'] and not git('diff', TESTED, expected, '--', 'Source', 'CMakeLists.txt'))
    subprocess.run(['git', 'merge-base', '--is-ancestor', CHECKPOINT, expected], cwd=ROOT, check=True)
    for name in ['Docs/ProjectIndex/91_AppliedRegionBoundaryEditing.md', 'Docs/ProjectIndex/92_ProposeApplyVerifyWorkflow.md', 'Tests/Qualification/ProposalWorkflow/results/2026-11-23/acceptance.json']:
        check('qualified_evidence:' + name, subprocess.check_output(['git', 'show', expected + ':' + name], cwd=ROOT) == subprocess.check_output(['git', 'show', CHECKPOINT + ':' + name], cwd=ROOT))
    finder_changes = []
    for state in before['worktrees']:
        path = Path(state['path'])
        check('worktree_head:' + str(path), git('rev-parse', 'HEAD', cwd=path) == (expected if path in [ROOT, INTEGRATION] else state['head']))
        check('worktree_branch:' + str(path), git('branch', '--show-current', cwd=path) == state['branch'])
        check('worktree_diff:' + str(path), hashlib.sha256(git('diff', 'HEAD', '--', cwd=path).encode()).hexdigest() == state['tracked_diff_sha256'])
        for name, digest in state['untracked'].items():
            actual = sha(path / name) if (path / name).is_file() else None
            if actual != digest and Path(name).name == '.DS_Store':
                finder_changes.append(str(path / name))
            else:
                check('unrelated_work:' + str(path / name), actual == digest)
    save(WORK / (label + '-verification.json'), {'passed': all(checks.values()), 'stage': label,
         'main': expected, 'checkpoint': CHECKPOINT, 'tag': TAG, 'tag_object': tag_object,
         'branch': BRANCH, 'source_tree': before['source_tree'], 'protected_files': len(before['protected']),
         'old_tags_preserved': len(before['local_tags_before']), 'packages_preserved': len(before['packages']),
         'candidate_source_entries_verified': before['candidate_source_entries_verified'],
         'finder_metadata_changes': finder_changes, 'checks': checks, 'runtime_tests_repeated': False,
         'rebuild_performed': False, 'physical_devices_opened': False, 'audible_tests_run': False})
    print(label, 'PASS:', len(checks), 'checks; main', expected, 'tag object', tag_object)

if __name__ == '__main__':
    if sys.argv[1] == 'preflight':
        preflight()
    else:
        verify(sys.argv[1], sys.argv[2])
