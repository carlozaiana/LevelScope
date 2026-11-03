"""Generate matching 60-second interleaved PCM24 fixtures outside the repository."""
from pathlib import Path
import math,sys,wave,hashlib
root=Path(sys.argv[1]);root.mkdir(parents=True,exist_ok=False)
rate=48000
for channels in (2,6):
    # Periodic finite PCM, with one high peak in the centre of a section.
    cycle=bytearray()
    for n in range(960):
        for c in range(channels):
            lfe=c==3
            amp=.012 if lfe else .06 if channels==2 else .034
            f=50 if lfe else 1000
            x=round(32768*amp*math.sin(2*math.pi*f*n/rate))
            cycle.extend((x*256).to_bytes(3,'little',signed=True))
    second=bytes(cycle)*50
    path=root/f'programme-{channels}.wav'
    with wave.open(str(path),'wb') as wav:
        wav.setparams((channels,3,rate,0,'NONE','not compressed'))
        for sec in range(60):
            payload=second
            if sec==31:
                at=(rate//2*channels+(3 if channels==6 else 0))*3
                payload=second[:at]+round(.8*2**23).to_bytes(3,'little',signed=True)+second[at+3:]
            wav.writeframesraw(payload)
    print(path,hashlib.sha256(path.read_bytes()).hexdigest())
