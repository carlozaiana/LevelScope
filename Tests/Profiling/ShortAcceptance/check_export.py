"""Independently check the published bundle and representative PCM24 samples."""
import argparse
import hashlib
import json
import pathlib
import struct

MASK64 = 2**64 - 1
STEP = 0x9E3779B97F4A7C15


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def chunks(path):
    with path.open("rb") as f:
        tag, size, kind = struct.unpack("<4sI4s", f.read(12))
        assert (tag, kind, size + 8) == (b"RIFF", b"WAVE", path.stat().st_size)
        found = {}
        while f.tell() < size + 8:
            tag, length = struct.unpack("<4sI", f.read(8))
            assert tag not in found
            found[tag] = (f.tell(), length)
            f.seek(length + length % 2, 1)
        assert f.tell() == size + 8
        return found


def mix(value):
    value &= MASK64
    value = ((value ^ (value >> 30)) * 0xBF58476D1CE4E5B9) & MASK64
    value = ((value ^ (value >> 27)) * 0x94D049BB133111EB) & MASK64
    return value ^ (value >> 31)


def noise(seed, frame, ch):
    key = seed ^ mix((ch + 1) * STEP)
    counter = key + frame * 2 + STEP
    return ((mix(counter) >> 11) - (mix(counter + 1) >> 11)) / 2**53


def f32(x):
    return struct.unpack("<f", struct.pack("<f", x))[0]


def local_db(plan, frame):
    rows = plan["applied_intervals"]
    for i, ramp in enumerate(plan["transitions"]):
        begin, end = map(int, (ramp["begin_frame"], ramp["end_frame"]))
        if begin <= frame < end:
            a, b = rows[i]["gain_db"], rows[i+1]["gain_db"]
            return a + (b-a) * (frame-begin) / (end-begin-1)
    return next((r["gain_db"] for r in rows if int(r["begin_frame"]) <= frame < int(r["end_frame"])), 0)


def run(directory, source):
    result = json.loads((directory / "result.json").read_text())
    bundle = directory / "export"
    report = json.loads((bundle / "report.json").read_text())
    manifest = json.loads((bundle / "manifest.json").read_text())
    assert result["passed"] and all(c["passed"] for c in result["checks"])
    names = {"programme.wav", "report.json", "report.html"}
    assert set(manifest["files_sha256"]) == names
    for name in names:
        assert digest(bundle / name) == manifest["files_sha256"][name]
    assert report["publication"] == "published-verified-pcm" and not report["session_unsaved"]
    art = report["artifact"]
    src_hash = digest(source)
    assert src_hash == result["input_sha256"] == report["source"]["measured_content_sha256"] == report["current"]["measured_content_sha256"]
    assert art["sha256"] == result["artifact_sha256"] == manifest["files_sha256"]["programme.wav"]
    assert art["evaluation"]["verdict"] == report["current"]["evaluation"]["verdict"] == "PASS"
    channels = len(report["ordered_input_roles"])
    assert report["ordered_input_roles"] == (["L", "R"] if channels == 2 else ["L", "R", "C", "Lfe", "Ls", "Rs"])
    assert art["output_to_input_channels"] == list(range(1, channels + 1))
    path = bundle / "programme.wav"
    wav = chunks(path)
    src = chunks(source)
    frames = 2880000
    assert wav[b"data"][1] == src[b"data"][1] == frames * channels * 3
    with path.open("rb") as out, source.open("rb") as original:
        out.seek(wav[b"fmt "][0])
        fmt = out.read(wav[b"fmt "][1])
        assert struct.unpack("<HHIIHHHHI", fmt[:24]) == (65534, channels, 48000, 48000*channels*3, channels*3, 24, 22, 24, 3 if channels == 2 else 63)
        assert fmt[24:] == bytes.fromhex("0100000000001000800000aa00389b71")
        plan = report["current_processing_plan"]
        selected = {0, 1, frames-1, 31*48000+24000}
        selected.update(range(0, frames, 48001))
        for t in plan["transitions"]:
            a, b = int(t["begin_frame"]), int(t["end_frame"])
            for middle in (a, (a+b)//2, b):
                selected.update((middle-1, middle, middle+1))
        seed = int(art["dither_seed_hex"], 16)
        for frame in sorted(selected):
            original.seek(src[b"data"][0] + frame * channels * 3)
            out.seek(wav[b"data"][0] + frame * channels * 3)
            for ch in range(channels):
                sample = int.from_bytes(original.read(3), "little", signed=True) / 2**23
                local = f32(sample * 10**(local_db(plan, frame)/20))
                global_sample = f32(local * 10**(plan["global_gain_db"]/20))
                expected = round(global_sample * 2**23 + noise(seed, frame, ch))
                actual = int.from_bytes(out.read(3), "little", signed=True)
                assert actual == expected, (frame, ch, actual, expected)
    # The diagnostic JSON writer uses 15 fractional decimal places. The native
    # driver checks exact plan identity at Apply. The artifact carries exact hex.
    assert float.fromhex(plan["global_gain_hex"]) == plan["global_gain_db"]
    assert abs(result["proposed_global_db"] - plan["global_gain_db"]) <= 5.1e-16
    return {"passed": True, "bundle_hashes_checked": sorted(names), "frames_per_channel": frames,
            "ordered_roles": report["ordered_input_roles"], "pcm_positions_per_channel_checked": len(selected),
            "exact_independent_pcm24_samples": len(selected)*channels, "session_sha256": digest(directory / "job.levelscope")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=pathlib.Path)
    parser.add_argument("source", type=pathlib.Path)
    args = parser.parse_args()
    print(json.dumps(run(args.directory, args.source), indent=2))
