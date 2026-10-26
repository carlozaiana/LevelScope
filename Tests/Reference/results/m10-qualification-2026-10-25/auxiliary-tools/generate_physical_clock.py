#!/usr/bin/env python3
"""M10 auxiliary media only: aligned stereo continuous tones, no endpoint impulse."""
import ctypes as C
import hashlib
import json
import math
from pathlib import Path

root = Path(__file__).parent
libpath = root / 'independent-decoder/_soundfile_data/libsndfile_arm64.dylib'
lib = C.CDLL(str(libpath))
class Info(C.Structure):
    _fields_ = [('frames', C.c_int64), ('samplerate', C.c_int), ('channels', C.c_int), ('format', C.c_int), ('sections', C.c_int), ('seekable', C.c_int)]
lib.sf_open.argtypes = [C.c_char_p, C.c_int, C.POINTER(Info)]
lib.sf_open.restype = C.c_void_p
lib.sf_writef_short.argtypes = [C.c_void_p, C.POINTER(C.c_short), C.c_int64]
lib.sf_writef_short.restype = C.c_int64
lib.sf_close.argtypes = [C.c_void_p]
lib.sf_close.restype = C.c_int
out = root / 'physical-clock-fixtures'
out.mkdir(exist_ok=False)
rows = []
for rate in (48000, 44100):
    folder = out / str(rate)
    folder.mkdir()
    for channel, amplitude in enumerate((1200, 1400)):
        table = [round(amplitude * math.sin(i * 2 * math.pi / 128)) for i in range(128)]
        block = (C.c_short * 131072)(*(table * 1024))
        path = folder / ('mono-%d.flac' % channel)
        info = Info(0, rate, 1, 0x170002, 0, 0)
        handle = lib.sf_open(str(path).encode(), 0x20, C.byref(info))
        assert handle
        frames = rate * 1810
        try:
            for first in range(0, frames, len(block)):
                n = min(len(block), frames-first)
                assert lib.sf_writef_short(handle, block, n) == n
        finally:
            assert lib.sf_close(handle) == 0
        rows.append(dict(file=str(path.relative_to(out)), rate=rate, frames=frames,
                         role=('Left','Right')[channel], pcm16_period=table,
                         sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
manifest = dict(purpose='Muted timing qualification only; independent generated aligned stereo',
                alignment='Sample zero aligned; exact frame counts; no automatic alignment',
                seconds=1810, generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                decoder_library_sha256=hashlib.sha256(libpath.read_bytes()).hexdigest(), files=rows)
(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps({k:v for k,v in manifest.items() if k!='files'},indent=2))
