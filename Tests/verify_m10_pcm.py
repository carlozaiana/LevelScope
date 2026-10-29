#!/usr/bin/env python3
"""Independent PCM24/TPDF oracle; standard-library only, Python 3.8 compatible."""
import argparse, hashlib, json, math, pathlib, struct
MASK = (1 << 64) - 1
STEP = 0x9e3779b97f4a7c15

def mix(v):
    v = ((v ^ (v >> 30)) * 0xbf58476d1ce4e5b9) & MASK
    v = ((v ^ (v >> 27)) * 0x94d049bb133111eb) & MASK
    return v ^ (v >> 31)

def noise(seed, frame, channel):
    key = seed ^ mix(((channel + 1) * STEP) & MASK)
    return ((mix((key + frame * 2 + STEP) & MASK) >> 11) -
            (mix((key + frame * 2 + 1 + STEP) & MASK) >> 11)) / 2**53

def inspect(path):
    size = path.stat().st_size
    with path.open('rb') as f:
        head = f.read(12)
        magic, length, wave = struct.unpack('<4sI4s', head)
        assert magic in (b'RIFF', b'RF64') and wave == b'WAVE'
        assert length == (0xffffffff if magic == b'RF64' else size - 8)
        ds = None; fmt = None; data = None
        while f.tell() < size:
            chunk, n = struct.unpack('<4sI', f.read(8)); pos = f.tell()
            if chunk == b'ds64':
                assert magic == b'RF64' and ds is None and pos == 20 and n == 28
                ds = struct.unpack('<QQQI', f.read(28)); assert ds[0] == size - 8 and ds[3] == 0
            elif chunk == b'fmt ':
                assert fmt is None and n == 40
                raw = f.read(40); fmt = struct.unpack('<HHIIHHHHI', raw[:24])
                assert raw[24:] == bytes.fromhex('0100000000001000800000aa00389b71')
            elif chunk == b'data':
                assert data is None and fmt is not None
                if magic == b'RF64':
                    assert ds and n == 0xffffffff; n = ds[1]
                data = (pos, n)
            assert pos + n + (n & 1) <= size
            f.seek(pos + n + (n & 1))
        assert data and fmt
        tag, channels, rate, byte_rate, align, bits, cb, valid, mask = fmt
        assert tag == 0xfffe and channels in (1, 2, 6) and bits == valid == 24 and cb == 22
        assert (mask in (1,2,4) if channels == 1 else mask == (3 if channels == 2 else 63))
        assert align == channels * 3 and byte_rate == rate * align
        assert size == data[0] + data[1] + (data[1] & 1)
        assert data[1] % align == 0
        frames = data[1] // align
        assert not ds or ds[2] == frames
        return dict(container=magic.decode(), channels=channels, rate=rate, frames=frames,
                    offset=data[0], bytes=data[1], mask=mask)

def verify(root):
    rows = json.loads((root / 'raw-oracles.json').read_text()); count = 0; zero_changed = False
    for row in rows:
        path = root / row['file']; info = inspect(path)
        assert all(info[k] == row[k] for k in ('channels','rate','frames'))
        frames, channels = row['frames'], row['channels']
        raw = (root / row['capture']).read_bytes(); samples = struct.unpack('<' + 'f' * (frames * channels), raw)
        with path.open('rb') as f:
            f.seek(info['offset']); pcm = f.read(info['bytes'])
        for frame in range(frames):
            for ch, src in enumerate(row['map']):
                sample = samples[src * frames + frame]
                expected = round(sample * 2**23 + noise(int(row['seed_hex'],16),frame,ch))
                start = (frame * channels + ch) * 3
                actual = int.from_bytes(pcm[start:start+3], 'little', signed=True)
                assert actual == expected, (row['file'], frame, ch, actual, expected)
                assert abs(actual / 2**23 - sample) <= 1.5 / 2**23
                if row['gain_db'] == 0 and actual / 2**23 != sample: zero_changed = True
                count += 1
    assert len(rows) == 48 and zero_changed
    return dict(passed=True, cases=len(rows), exact_codes=count, zero_gain_is_not_bit_identical=zero_changed,
                oracle='independent Python integer container parser; binary64 TPDF and ties-even round')

if __name__ == '__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=pathlib.Path);p.add_argument('--output',type=pathlib.Path);a=p.parse_args()
    result=verify(a.root);text=json.dumps(result,indent=2)+'\n'
    if a.output:a.output.write_text(text)
    print(text)
