"""Additional conservative audit: include loading in the OS RSS high-water ceiling."""
from pathlib import Path
import hashlib
import json
import math

B = Path("/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020")
rows = []
for platform in ["sequoia", "high-sierra"]:
    for phase, count in [("matrix", 16), ("endurance", 3)]:
        root = B / "evidence" / f"{platform}-{phase}-candidate04"
        files = [p for p in root.glob("*.json") if p.stem not in ["identity", "runs", "summary"]]
        assert len(files) == count, (platform, phase, len(files))
        for p in sorted(files):
            r = json.loads(p.read_text())
            value = r.get("process_peak_rss_including_load_bytes")
            passed = (isinstance(value, (int, float)) and math.isfinite(value)
                      and 0 < value <= 768 * 1024**2)
            rows.append(dict(platform=platform, phase=phase, name=p.stem,
                             process_highwater_bytes=value, passed=passed,
                             report_sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
out = B / "evidence/process-highwater-candidate04.json"
assert not out.exists()
result = dict(passed=all(r["passed"] for r in rows), limit_bytes=768 * 1024**2,
              method="getrusage(RUSAGE_SELF).ru_maxrss, macOS bytes; includes loading and playback",
              scope="Additional conservative check alongside unchanged coherent-snapshot memory gates",
              runs=rows)
out.write_text(json.dumps(result, indent=2) + "\n")
print("Process RSS high-water including loading:", "PASS" if result["passed"] else "FAIL")
raise SystemExit(0 if result["passed"] else 1)
