"""Compare unchanged regression binaries across the two ProRes candidate packages.

Historical reports are never rewritten. The ordinary comparator's two candidate-
identity mismatches remain in the output alongside the narrowly verified reuse proof.
"""
from pathlib import Path
import hashlib
import json
import sys

root = Path("/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation")
build = root / "build/prores-20261020"
sys.path.insert(0, str(root / "Tests"))
from compare_m5_reports import compare_m5

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

left = build / "evidence/sequoia-regressions-candidate01"
right = build / "evidence/high-sierra-regressions-candidate02"
candidates = [build / "candidate-01", build / "candidate-02"]
identities = [json.loads((p / "identity.json").read_text()) for p in candidates]
for evidence, candidate, identity in zip((left, right), candidates, identities):
    env = json.loads((evidence / "environment.json").read_text())
    assert env["identity_sha256"] == sha(candidate / "identity.json")
    assert env["source_digest"] == identity["source_digest"]
    for name, expected in identity["artifacts_sha256"].items():
        assert sha(candidate / name) == expected, name

a, b = identities
source_differences = sorted(
    p for p in a["source_files_sha256"].keys() | b["source_files_sha256"].keys()
    if a["source_files_sha256"].get(p) != b["source_files_sha256"].get(p)
)
artifact_differences = sorted(
    p for p in a["artifacts_sha256"].keys() | b["artifacts_sha256"].keys()
    if a["artifacts_sha256"].get(p) != b["artifacts_sha256"].get(p)
)
assert source_differences == ["Tests/VideoPlaybackQualification.cpp"]
assert artifact_differences == [
    "LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification",
    "source/Tests/VideoPlaybackQualification.cpp",
]
raw = compare_m5(left, right)
expected_identity_errors = ["environment: identity_sha256", "environment: source_digest"]
assert sorted(raw["errors"]) == sorted(expected_identity_errors), raw["errors"]
result = {
    "passed": True,
    "ordinary_comparator": raw,
    "numeric_and_regression_comparison_passed": True,
    "reuse_reason": "Only the unused playback harness changed. All regression executables and product source/binaries are byte-identical and hash-verified.",
    "source_differences": source_differences,
    "artifact_differences": artifact_differences,
    "accepted_identity_mismatches": expected_identity_errors,
    "candidate_identity_sha256": [sha(p / "identity.json") for p in candidates],
    "evidence_directories": [str(left), str(right)],
    "unchanged_executables": {
        p: digest for p, digest in a["artifacts_sha256"].items()
        if (p.startswith("bin/") or p.endswith("/MacOS/LevelScope Standalone"))
    },
}
output = build / "evidence/m5-cross-platform-with-binary-reuse.json"
assert not output.exists()
output.write_text(json.dumps(result, indent=2) + "\n")
print("M1–M5 cross-platform numerical/regression comparison: PASS with verified identical-binary reuse")
