"""Read-only extraction of a bounded PCM24 excerpt; no dither or resampling.

This fixture utility is not a product importer. The original's known content
hash and previously qualified RF64 metadata/layout must match before extraction.
"""
import argparse
import hashlib
import json
from pathlib import Path
import struct

p = argparse.ArgumentParser()
p.add_argument('source', type=Path)
p.add_argument('destination', type=Path)
a = p.parse_args()
assert not a.destination.exists()
before = a.source.stat()
h = hashlib.sha256()
with a.source.open('rb') as f:
    for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
        h.update(block)
original = h.hexdigest()
assert original == '1fb7f93a33d600bf0a62d63e642e5703f03aec8e6920f8a6e510cf65f61fdeaf'
with a.source.open('rb') as f:
    assert f.read(12) == b'RF64\xff\xff\xff\xffWAVE'
    fmt = None
    while True:
        name, size = struct.unpack('<4sI', f.read(8))
        if name == b'ds64':
            assert 28 <= size <= 4096
            ds = f.read(size)
            _, data_bytes, samples, entries = struct.unpack('<QQQI', ds[:28])
            assert samples == entries == 0
        elif name == b'fmt ':
            assert size == 40
            fmt = f.read(size)
            assert struct.unpack('<HHIIHH', fmt[:16]) == (65534, 6, 48000, 864000, 18, 24)
            assert struct.unpack('<HHI', fmt[16:24]) == (22, 24, 63)
        elif name == b'data':
            assert size == 0xffffffff and fmt and data_bytes == 4857753600
            offset = f.tell()
            break
        else:
            f.seek(size + (size & 1), 1)
    begin, end = 27 * 60 * 48000, 35 * 60 * 48000
    size = (end - begin) * 18
    assert end <= data_bytes // 18 and size < 2**32 - 60
    a.destination.parent.mkdir(parents=True, exist_ok=True)
    with a.destination.open('xb') as out:
        out.write(b'RIFF' + struct.pack('<I', 4 + 8 + len(fmt) + 8 + size) + b'WAVE')
        out.write(b'fmt ' + struct.pack('<I', len(fmt)) + fmt + b'data' + struct.pack('<I', size))
        payload = hashlib.sha256()
        f.seek(offset + begin * 18)
        remaining = size
        while remaining:
            data = f.read(min(4 * 1024 * 1024, remaining))
            assert data
            out.write(data)
            payload.update(data)
            remaining -= len(data)
    # Byte comparison of the entire bounded payload is stronger than a few probes.
    with a.destination.open('rb') as check:
        check.seek(68)
        actual = hashlib.sha256()
        for block in iter(lambda: check.read(4 * 1024 * 1024), b''):
            actual.update(block)
        assert actual.digest() == payload.digest()
after = a.source.stat()
assert (before.st_size, before.st_mtime_ns, before.st_ino) == (after.st_size, after.st_mtime_ns, after.st_ino)
output = hashlib.sha256()
with a.destination.open('rb') as f:
    for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
        output.update(block)
record = dict(original_path=str(a.source), original_sha256=original,
    original_size=before.st_size, original_mtime_ns=before.st_mtime_ns,
    sample_rate=48000, channels=['L', 'R', 'C', 'LFE', 'Ls', 'Rs'], bits=24,
    begin_source_frame=begin, end_source_frame=end, seconds=480,
    excerpt_sha256=output.hexdigest(), pcm_payload_sha256=payload.hexdigest(),
    payload_byte_identical=True, dither=False, resampling=False,
    selection='27:00–35:00 audio-source time, encompassing the prior 28:15–31:31 marker cluster; measured quiet/loud supports will be recorded from this excerpt.',
    interpretation='Bounded excerpt, not full-film measurement or Nuendo comparison.')
a.destination.with_suffix('.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
