"""Verify final-candidate media identity, painted coverage and fixed acceptance gates."""
from pathlib import Path
import argparse, hashlib, json, math

def verify_memory(report):
    """Recompute the approved gate from raw coherent snapshots, independently of C++."""
    failures = []
    if report.get("memory_method") != "prores-memory-v2-user-approved-2026-10-20":
        failures.append("memory method")
    if report.get("qualification_eligible") is not True or report.get("allocator_relief_requested") is not False:
        failures.append("diagnostic/allocator exclusion")
    if report.get("memory_query_failures") != 0:
        failures.append("VM query failure")
    samples = report.get("memory_samples", [])
    if not samples:
        return failures + ["missing memory samples"]
    fields = ["rss_bytes", "reusable_bytes", "non_reusable_resident_bytes", "physical_footprint_bytes"]
    seconds = []
    for sample in samples:
        if not all(isinstance(sample.get(key), (float, int)) and math.isfinite(sample[key]) and sample[key] >= 0 for key in fields):
            return failures + ["invalid memory counter"]
        second = sample.get("seconds")
        if not isinstance(second, int) or second < 2 or second > math.ceil(report["source_seconds"]):
            failures.append("memory sample timestamp")
        seconds.append(second)
        if sample["rss_bytes"] <= 0 or sample["physical_footprint_bytes"] <= 0 or sample["reusable_bytes"] > sample["rss_bytes"]:
            failures.append("inconsistent VM snapshot")
        if sample["non_reusable_resident_bytes"] != sample["rss_bytes"] - sample["reusable_bytes"]:
            failures.append("non-reusable calculation")
    if seconds != sorted(set(seconds)):
        failures.append("duplicate/unordered memory snapshots")
    if max(s["rss_bytes"] for s in samples) > 768 * 1024**2:
        failures.append("raw RSS peak ceiling")
    duration = report["requested_seconds"]
    first = [s for s in samples if s["seconds"] < 62]
    last = [s for s in samples if s["seconds"] > duration - 60]
    if duration >= 120 and (len(first) < 55 or len(last) < 55):
        failures.append("insufficient growth-window samples")
    if report.get("first_minute_memory_sample_count") != len(first) or report.get("last_minute_memory_sample_count") != len(last):
        failures.append("sample count rollup")
    for key, peak_key, first_key, last_key in [
        ("rss_bytes", "peak_sampled_playback_rss_bytes", "first_minute_mean_rss_bytes", "last_minute_mean_rss_bytes"),
        ("non_reusable_resident_bytes", "peak_non_reusable_resident_bytes", "first_minute_mean_non_reusable_resident_bytes", "last_minute_mean_non_reusable_resident_bytes"),
        ("physical_footprint_bytes", "peak_physical_footprint_bytes", "first_minute_mean_physical_footprint_bytes", "last_minute_mean_physical_footprint_bytes"),
    ]:
        f = sum(s[key] for s in first) / len(first) if first else 0
        l = sum(s[key] for s in last) / len(last) if last else 0
        for label, value in [(peak_key, max(s[key] for s in samples)), (first_key, f), (last_key, l)]:
            if not math.isclose(report.get(label, -1), value, rel_tol=1e-12, abs_tol=0.01):
                failures.append("memory rollup " + label)
        if duration >= 120 and key != "rss_bytes" and l - f > 32 * 1024**2:
            failures.append("growth gate " + key)
    return failures

p=argparse.ArgumentParser()
p.add_argument("--root",type=Path,required=True)
p.add_argument("--phase",choices=["matrix","endurance"],required=True)
p.add_argument("--candidate",choices=["04"],default="04")
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
        failures.extend(verify_memory(r))
        if r.get("hal_output_timestamp") is not True: failures.append("hardware presentation clock")
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
