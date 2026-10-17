#!/usr/bin/env python3
"""Generate paired stereo clicks for the 30000/1001 burn-in movie.

Flash every 30 video frames. Click positions round the exact rational time
to the nearest audio sample; errors never accumulate.
The file is test media, stays outside Git, and is never played by this script.
"""
import argparse
import hashlib
import json
import math
import pathlib
import struct
import wave

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('output', type=pathlib.Path)
parser.add_argument('--video-frames', type=int, default=54000)
parser.add_argument('--sample-rate', type=int, choices=[44100, 48000], default=48000)
parser.add_argument('--peak-amplitude', type=float, choices=[0.1, 0.4], default=0.1)
args = parser.parse_args()
if args.output.exists() or args.video_frames <= 0:
    raise SystemExit('Choose a new output path and a positive frame count')
rate = args.sample_rate
frames = (args.video_frames * 1001 * rate + 29999) // 30000
pulse_frames = rate // 200
pulse = bytearray(pulse_frames * 4)
for index in range(pulse_frames):
    value = round(32767 * args.peak_amplitude * math.sin(2 * math.pi * 1000 * index / rate)
                  * math.sin(math.pi * index / pulse_frames) ** 2)
    struct.pack_into('<hh', pulse, index * 4, value, value)
with wave.open(str(args.output), 'wb') as writer:
    writer.setparams((2, 2, rate, frames, 'NONE', 'not compressed'))
    position = 0
    event = 1
    while position < frames:
        end = min(frames, (event * 1001 * rate + 500) // 1000)
        cycle = bytearray((end - position) * 4)
        cycle[:min(len(cycle), len(pulse))] = pulse[:min(len(cycle), len(pulse))]
        writer.writeframesraw(cycle)
        position = end
        event += 1
digest = hashlib.sha256()
with args.output.open('rb') as handle:
    for block in iter(lambda: handle.read(1024 * 1024), b''):
        digest.update(block)
print(json.dumps(dict(file=str(args.output), sha256=digest.hexdigest(),
    sample_rate=rate, channels=2, frames=frames, click_period_seconds=[1001, 1000],
    click_rounding="nearest source sample",
    flash_period_video_frames=30, video_rate=[30000, 1001],
    peak_amplitude=args.peak_amplitude, pulse_frames=pulse_frames), indent=2))
