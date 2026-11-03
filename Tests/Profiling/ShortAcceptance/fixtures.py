"""Generate two short programmes with distinct global and section proposal outcomes."""
from pathlib import Path
import runpy
import sys

folder = Path(sys.argv[1])
runpy.run_path(str(Path(__file__).resolve().parents[1] / "OfflineStages" / "fixtures.py"))
# That writer produces ordinary, 44-byte RIFF/PCM24 headers. Keep the original
# 5.1 crest in LFE. Lower the stereo crest enough for global trim alone to pass.
sample = (31 * 48000 + 24000) * 2
with (folder / "programme-2.wav").open("r+b") as audio:
    audio.seek(44 + 3 * sample)
    assert int.from_bytes(audio.read(3), "little", signed=True) == round(.8 * 2**23)
    audio.seek(-3, 1)
    audio.write(round(.3 * 2**23).to_bytes(3, "little", signed=True))
