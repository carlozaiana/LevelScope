"""Read-only before/after audio device configuration check after qualification."""
from pathlib import Path
import json
import subprocess
import hashlib

T = Path("/private/tmp/levelscope-prores-20261020")
B = Path("/Users/fuzz/Developer/LevelScope-worktrees/VideoImplementation/build/prores-20261020")
before_local = (T / "sequoia-audio-before.txt").read_text()
before_remote = (T / "golem-preflight.txt").read_text().split("/Users/fuzz/", 1)[0]
after_local = subprocess.check_output(
    ["/usr/sbin/system_profiler", "SPAudioDataType"], text=True)
after_remote = subprocess.check_output(
    ["ssh", "-S", "/tmp/levelscope-prores-golem", "-o", "BatchMode=yes",
     "fuzz@192.168.2.2", "/usr/sbin/system_profiler SPAudioDataType"], text=True)
def normalized(s):
    return "\n".join(line.strip() for line in s.splitlines() if line.strip())
rows = []
for platform, before, after in [
    ("sequoia", before_local, after_local),
    ("high-sierra", before_remote, after_remote),
]:
    rows.append(dict(platform=platform, unchanged=normalized(before) == normalized(after),
                     before_audio_devices=before.strip(), after_audio_devices=after.strip()))
out = B / "evidence/audio-configuration-after.json"
assert not out.exists()
report = dict(
    passed=all(r["unchanged"] for r in rows),
    scope="Read-only system_profiler device/default/input/output/source/sample-rate snapshot comparison",
    limitation="This comparison covers exposed device settings; no claim to inventory every CoreAudio property.",
    observations=rows,
)
out.write_text(json.dumps(report, indent=2) + "\n")
print("Audio configuration snapshot comparison:", "PASS" if report["passed"] else "DIFFERENCE")
raise SystemExit(0 if report["passed"] else 1)
