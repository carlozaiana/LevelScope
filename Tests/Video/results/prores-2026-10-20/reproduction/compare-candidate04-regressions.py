from pathlib import Path
import hashlib,json,sys
R=Path('/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation');B=R/'build/prores-20261020'
sys.path.insert(0,str(R/'Tests'));from compare_m5_reports import compare_m5
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
left=B/'evidence/sequoia-regressions-candidate03';right=B/'evidence/high-sierra-regressions-candidate04';packages=[B/'candidate-03',B/'candidate-04'];ids=[json.loads((p/'identity.json').read_text()) for p in packages]
reuse=json.loads((B/'evidence/candidate04-native-regression-ui-reuse.json').read_text());assert reuse['passed']
for ev,p,id in zip([left,right],packages,ids):
 env=json.loads((ev/'environment.json').read_text());assert env['identity_sha256']==sha(p/'identity.json') and env['source_digest']==id['source_digest']
 for name,want in id['artifacts_sha256'].items():assert sha(p/name)==want,name
srcdiff=sorted(k for k in ids[0]['source_files_sha256'].keys()|ids[1]['source_files_sha256'].keys() if ids[0]['source_files_sha256'].get(k)!=ids[1]['source_files_sha256'].get(k))
assert srcdiff==reuse['source_differences']
diff=sorted(k for k in ids[0]['artifacts_sha256'].keys()|ids[1]['artifacts_sha256'].keys() if ids[0]['artifacts_sha256'].get(k)!=ids[1]['artifacts_sha256'].get(k))
assert diff==['LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification','Tests/run_prores_playback_qualification.py','Tests/verify_prores_playback_evidence.py','source/Tests/VideoPlaybackQualification.cpp','source/Tests/run_prores_playback_qualification.py','source/Tests/verify_prores_playback_evidence.py']
raw=compare_m5(left,right);expected=['environment: identity_sha256','environment: source_digest'];assert sorted(raw['errors'])==sorted(expected),raw
p=B/'evidence/m5-cross-platform-candidate04-ordinary.json';assert not p.exists();p.write_text(json.dumps(raw,indent=2)+'\n')
result=dict(passed=True,numeric_and_regression_comparison_passed=True,ordinary_comparator=raw,accepted_identity_mismatches=expected,source_differences=srcdiff,artifact_differences=diff,reuse_proof_sha256=sha(B/'evidence/candidate04-native-regression-ui-reuse.json'),evidence_directories=[str(left),str(right)],candidate_identity_sha256=[sha(p/'identity.json') for p in packages],reuse_reason='All production source and 8 Standalone/native/M1-M5 binaries are identical; only the unused playback harness and its runner/verifier differ. No historical report is rewritten.')
p=B/'evidence/m5-cross-platform-candidate04-with-reuse.json';assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n');print('PASS exact-binary cross-platform M1–M5/H264/ITU/EBU comparison')
