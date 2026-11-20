"""Read the archived baseline observations. Run with Python 3, from any directory."""

import csv
import hashlib
import io
import json
import math
from itertools import groupby
from pathlib import Path
from statistics import median, median_high
from zipfile import ZipFile


def read_rows(archive, name):
    with archive.open(name) as raw:
        return list(csv.DictReader(io.TextIOWrapper(raw, encoding="utf-8")))


def coarse_edges(windows, diagnostic, allowed_begin, allowed_end):
    """The two existing centre-based nomination tests, prior to ±1.5 s expansion."""
    settings = diagnostic["settings"]
    rate = diagnostic["sample_rate"]
    step = rate // 10
    halo = 3 * rate
    persistence = settings["persistenceMs"] * rate // 1000
    threshold = settings["needTenths"] / 10
    # The diagnostic used the 18 LU requirement recorded in the frozen report.
    limit = diagnostic["p10"] + 18
    samples = []
    for row in windows:
        centre = (int(row["begin"]) + int(row["end"])) // 2
        if row["lufs"] and allowed_begin - halo <= centre < allowed_end + halo:
            samples.append((centre, max(0, float(row["lufs"]) - limit)))

    def values(begin, end):
        return [q for frame, q in samples if begin <= frame < end]

    edges = []
    for f in range(allowed_begin, allowed_end, step):
        left = median_high(values(f - halo, f))
        right = median_high(values(f, f + halo))
        middle = (left + right) / 2
        before = values(f - persistence, f)
        after = values(f, f + persistence)
        if len(before) < persistence // step or len(after) < persistence // step:
            continue
        if abs(left - right) < threshold:
            continue
        if left < right and max(before) <= middle <= min(after):
            edges.append((f / rate, "rising median", right - left))
        if left > right and min(before) >= middle >= max(after):
            edges.append((f / rate, "falling median", left - right))

    for is_low, group in groupby(samples, key=lambda item: item[1] <= threshold / 2):
        run = list(group)
        begin, end = run[0][0], run[-1][0] + step
        if not is_low or end - begin < persistence:
            continue
        previous = values(begin - halo, begin)
        following = values(end, end + halo)
        left = median_high(previous) if previous else -1
        right = median_high(following) if following else -1
        strength = max(left, right)
        if strength > threshold / 2:
            for f, label in ((begin, "start of low run"), (end, "end of low run")):
                if allowed_begin < f < allowed_end:
                    edges.append((f / rate, label, strength))
    return sorted(edges)


def first_quiet_run(trace, decision, after, before, rate, hop):
    """Conditional 300 ms return to the background at one recorded decision."""
    bed = trace[decision]
    channels = (("fast_l", "bg_l"), ("fast_r", "bg_r"),
                ("low_l", "bg_low_l"), ("low_r", "bg_low_r"))
    limits = [bed[bg] * 10 ** 0.3 if bed[bg] > 1e-8 else 1e-8 for _, bg in channels]
    run_begin = None
    for end in range(after + hop, before + 1, hop):
        if all(trace[end][name] <= limit for (name, _), limit in zip(channels, limits)):
            if run_begin is None:
                run_begin = end - hop
            if end - run_begin >= 300 * rate // 1000:
                return run_begin, end
        else:
            run_begin = None
    return None


def main():
    directory = Path(__file__).resolve().parent
    manifest = json.loads((directory / "provenance.json").read_text())
    with ZipFile(directory / "trace-evidence.zip") as archive:
        for name, expected in manifest["files_sha256"].items():
            if hashlib.sha256(archive.read(name)).hexdigest() != expected:
                raise ValueError(f"Evidence digest mismatch for {name}")
        diagnostic = json.loads(archive.read("diagnostic.json"))
        windows = read_rows(archive, "windows.csv")
        trace = {int(row["frame"]): {k: float(v) for k, v in row.items()}
                 for row in read_rows(archive, "trace.csv")}

    rate = diagnostic["sample_rate"]
    hop = diagnostic["settings"]["eventFeatures"]["hopMs"] * rate // 1000
    frame = lambda seconds: round(seconds * rate)
    momentary = {int(r["begin"]): float(r["energy"]) for r in windows if r["kind"] == "momentary"}
    high = [r for r in windows if r["kind"] == "lra"]
    short_term = {int(r["begin"]): float(r["energy"]) for r in high}

    def tiled_energy(start, end):
        return sum(momentary[f] * 0.4 for f in range(frame(start), frame(end), frame(0.4)))

    shares = (tiled_energy(52.5, 52.9) / (3 * short_term[frame(49.9)]),
              tiled_energy(57.2, 58.0) / (3 * short_term[frame(57.2)]),
              tiled_energy(52.4, 58.0) / tiled_energy(50, 60))
    print("Later passage, baseline energy shares (first window, last window, ten-second sample)")
    print(*(f"{100 * share:.4f}%" for share in shares))
    print("Coarse later-root nominations before support offsets (seconds, origin, strength)")
    for location, kind, strength in coarse_edges(high, diagnostic, frame(48.9), frame(61.2)):
        print(f"{location:.2f}  {kind:17}  {strength:.6f}")

    print("Selected recorded attacks and their 100 ms median fast / low references (dB)")
    for t in (11.80, 12.07, 12.34, 12.86, 13.39, 29.18, 29.71, 30.26, 31.29,
              51.91, 52.00, 52.52, 57.34, 58.17):
        start = frame(t)
        cells = [trace[f] for f in range(start, start + frame(0.1) + 1, hop)]
        if not cells[0]["attack"]:
            raise ValueError(f"No recorded attack at {t}")
        refs = [median((c[band + "_l"] + c[band + "_r"]) / 2 for c in cells)
                for band in ("fast", "low")]
        print(f"{t:5.2f}  {10 * math.log10(refs[0]):6.2f}  {10 * math.log10(refs[1]):6.2f}")

    print("First 300 ms local-bed return in each ending inspection interval (source frames)")
    for decision, after, before in ((12.34, 30.35, 31.26), (52.52, 57.78, 59.70)):
        print(decision, first_quiet_run(trace, frame(decision), frame(after), frame(before), rate, hop))


if __name__ == "__main__":
    main()
