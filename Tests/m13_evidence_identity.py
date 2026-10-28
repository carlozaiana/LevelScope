"""Explicit equivalence for M13 qualification-tool corrections; never product fixes."""
import hashlib,json,pathlib
from collect_m1_environment import sha256

TOOL_ONLY_PATHS={
    'Tests/run_m13_muted_video.py','Tests/M13QualificationRunnerTests.py',
    'Tests/complete_m13_qualification.py','Tests/compare_m13_reports.py',
    'Tests/m13_evidence_identity.py'}

def equivalent(previous,current):
    for key in ('measurement_build','build_configuration','compiler','juce_commit'):
        if previous.get(key)!=current.get(key):raise ValueError('Different production build: '+key)
    left,right=previous['source_files_sha256'],current['source_files_sha256']
    changed={p for p in set(left)|set(right) if left.get(p)!=right.get(p)}
    if not changed.issubset(TOOL_ONLY_PATHS):raise ValueError('Non-tool source changes require affected qualification')
    allowed=TOOL_ONLY_PATHS|{'source/'+p for p in TOOL_ONLY_PATHS}
    a,b=previous['artifacts_sha256'],current['artifacts_sha256']
    differences={p for p in set(a)|set(b) if a.get(p)!=b.get(p)}
    if not differences.issubset(allowed):raise ValueError('Binary or non-tool artifact changed')
    return dict(previous_source_digest=previous['source_digest'],current_source_digest=current['source_digest'],
        changed_tool_sources=sorted(changed),production_and_test_binaries_unchanged=True,
        reason='Qualification-tool repair only. Original phase receipts remain unchanged; corrected runner is tested freshly.')

def identities(candidate,reference=None):
    candidates=[pathlib.Path(candidate)]
    if reference is not None:candidates.append(pathlib.Path(reference))
    records=[]
    for folder in candidates:
        file=folder/'identity.json';record=json.loads(file.read_text())
        for name,digest in record['artifacts_sha256'].items():
            if sha256(folder/name)!=digest:raise ValueError('Frozen artifact changed: '+str(folder/name))
        expected=hashlib.sha256(json.dumps(record['source_files_sha256'],sort_keys=True,separators=(',',':')).encode()).hexdigest()
        if record['source_digest']!=expected:raise ValueError('Invalid source manifest digest')
        records.append((record,sha256(file)))
    reuse=equivalent(records[1][0],records[0][0]) if len(records)==2 else None
    return {record['source_digest']:(record,digest) for record,digest in records},reuse

def select(records,source_digest):
    if source_digest not in records:raise ValueError('Unknown or non-equivalent evidence candidate')
    return records[source_digest]
