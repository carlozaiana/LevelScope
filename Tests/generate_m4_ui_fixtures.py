#!/usr/bin/env python3
"""Generate M4-only silent UI inputs in a new folder; never open an audio device.

Copy earlier M3 UI fixtures with --prior-fixtures when packaging a full candidate.
These PCM16 stimuli are descriptive UI fixtures, not normative reference vectors.
The separate three-hour measurement-sized benchmark lives in DiagnosticsTests.cpp.
"""
import argparse
import math
from pathlib import Path
import shutil
import struct
import wave


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--prior-fixtures', type=Path)
    p.add_argument('--three-hour', action='store_true', help='Also generate a 3-hour stereo PCM16 UI endurance input (about 2.1 GB)')
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    if args.prior_fixtures:
        for name in ('stereo.wav', 'stereo.aiff', 'stereo.flac', 'long-stereo.wav',
                     'unsupported-12ch.wav', 'trim-required.wav', 'gate-crossing.wav'):
            shutil.copy2(args.prior_fixtures / name, args.output / name)
    rate = 48000
    cycle = [math.sin(2 * math.pi * n / 48) for n in range(48)]
    for name, seconds in (('scenes', 120), ('silence', 4), ('long-silence', 120), ('short', 2)):
        with wave.open(str(args.output / ('diagnostic-' + name + '.wav')), 'wb') as out:
            out.setparams((2, 2, rate, 0, 'NONE', 'not compressed'))
            for second in range(seconds):
                db = (-80 if second < 10 else -40 if second < 30 else
                      -23 if second < 80 else -16 if second < 100 else -35)
                amplitude = 0 if 'silence' in name else 32767 * 10 ** ((-23 if name == 'short' else db) / 20)
                block = bytearray()
                for n in range(rate):
                    left = round(amplitude * cycle[n % 48])
                    right = 30000 if name == 'scenes' and second == 90 and n == 2400 else left
                    block.extend(struct.pack('<hh', left, right))
                out.writeframesraw(block)
    if args.three_hour:
        # Repeat the finite 120-second scene stimulus without holding hours of
        # PCM in memory. This supplies an actual long-file UI measurement in
        # addition to the independent generated-snapshot performance benchmark.
        with wave.open(str(args.output / 'diagnostic-scenes.wav'), 'rb') as source:
            block = source.readframes(source.getnframes())
        with wave.open(str(args.output / 'diagnostic-three-hour.wav'), 'wb') as out:
            out.setparams((2, 2, rate, 0, 'NONE', 'not compressed'))
            for _ in range(90):
                out.writeframesraw(block)


if __name__ == '__main__':
    main()
