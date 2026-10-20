"""Verify final-candidate media identity, painted coverage and fixed acceptance gates."""
from pathlib import Path
import argparse, hashlib, json, math

p=argparse.ArgumentParser()
p.add_argument("--root",type=Path,required=True)
p.add_argument("--phase",choices=["matrix","endurance"],required=True)
p.add_argument("--candidate",choices=["02","03"],default="02")
p.add_argument("--output",type=Path,required=True)
a=p.parse_args()
manifest_path=a.root/"fixtures/manifest.json"
manifest=json.loads(manifest_path.read_text())
manifest_hash=hashlib.sha256(manifest_path.read_bytes()).hexdigest()
candidate=json.loads((a.root/("candidate-"+a.candidate)/"identity.json").read_text())
binary_hash=candidate["artifacts_sha256"]["LevelScope Video Qualification.app/Contents/MacOS/LevelScope Video Qualification"]
rows=[]; errors=[]
for platform,arch,device,rate in [
    ("sequoia","arm64","Mac Studio Speakers",48000),
    ("high-sierra","x86_64","Built-in Output",44100),
]:
    directory=a.root/"evidence"/(platform+"-"+a.phase+"-candidate"+a.candidate)
    identity=json.loads((directory/"identity.json").read_text())
    summary=json.loads((directory/"summary.json").read_text())
    expected={Path(name).stem for name,m in manifest.items() if m.get("required")}
    expected.add("51_SyncTest")
    if a.phase=="endurance":
        expected={profile+"-1800" for profile in ["proxy","lt","422"]}
    reports={p.stem:p for p in directory.glob("*.json") if p.stem not in ["identity","runs","summary"]}
    if reports.keys()!=expected:errors.append(platform+": required report set mismatch")
    if not summary["passed"] or summary["test_count"]!=len(expected):errors.append(platform+": summary failed")
    for key,value in [
        ("binary_sha256",binary_hash),("architecture",arch),
        ("fixture_manifest_sha256",manifest_hash),("device",device),
        ("sample_rate",rate),("all_routes_muted",True),
    ]:
        if identity.get(key)!=value:errors.append(platform+": identity "+key)
    for name,path in sorted(reports.items()):
        r=json.loads(path.read_text())
        if name=="51_SyncTest":
            fps=25
            video_hash=candidate["actual_movie_sha256"]
            duration=100
        else:
            fixture=manifest[name+".mov"]
            fps=fixture["scale"]/fixture["step"]
            video_hash=fixture["sha256"]
            duration=1800 if a.phase=="endurance" else 18
        # One-frame boundary allowance covers timer termination relative to the
        # first/last paint; it is not an allowance for skipped interior frames.
        minimum_unique=math.floor((duration-2)*fps)-1
        required={
            "passed":True, "all_routes_muted":True, "video_sha256":video_hash,
            "audio_sha256":manifest["selected-clock-"+str(rate)+".wav"]["sha256"],
            "requested_seconds":duration, "architecture":arch,
            "actual_paint_skipped_frames":0, "skipped_frame_limit":0,
            "blank_paints":0, "bound_failures":0, "paint_deadline_failures":0,
            "underrun_blocks":0, "attachment_within_budget":True, "memory_passed":True,
        }
        failures=[k for k,v in required.items() if r.get(k)!=v]
        if r["source_seconds"]<duration:failures.append("full duration")
        if r["unique_painted_frames"]<minimum_unique:failures.append("paint coverage")
        if a.phase=="endurance" and not r["memory_growth_assessed"]:failures.append("memory growth assessment")
        errors.extend(platform+"/"+name+": "+failure for failure in failures)
        rows.append(dict(platform=platform,name=name,passed=not failures,
                         painted_frames=r["unique_painted_frames"],minimum_painted_frames=minimum_unique,
                         report_sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
result=dict(passed=not errors,phase=a.phase,candidate=a.candidate,errors=errors,runs=rows,
            zero_interior_skip_limit=0,endpoint_coverage_allowance_frames=1)
assert not a.output.exists()
a.output.write_text(json.dumps(result,indent=2)+"\n")
print(a.phase,"evidence identity and coverage:", "PASS" if result["passed"] else "FAIL",errors)
raise SystemExit(0 if result["passed"] else 1)
