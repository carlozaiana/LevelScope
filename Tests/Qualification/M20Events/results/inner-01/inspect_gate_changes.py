"""Read the verified listening sessions and compare measured gate membership by end frame."""
import hashlib
import io
import json
from pathlib import Path
import struct
import sys


def populations(file):
    data=file.read_bytes()
    assert data[:8]==b'LSSESS01'
    length=struct.unpack_from('<Q',data,8)[0]
    core=data[80:80+length]
    assert hashlib.sha256(core).hexdigest()==data[16:80].decode()
    result=[]
    for entry in json.loads(core)['caches']:
        begin=80+length+int(entry['offset'])
        raw=data[begin:begin+int(entry['size'])]
        assert hashlib.sha256(raw).hexdigest()==entry['sha256']
        stream=io.BytesIO(raw)
        def u(): return struct.unpack('<Q',stream.read(8))[0]
        def number(): return struct.unpack('<d',stream.read(8))[0]
        header=json.loads(stream.read(u()))
        rate,frames=u(),u()
        assert (rate,frames)==(48000,13228800)
        for _ in range(6):
            u()
            if u(): number()
        for _ in range(4):
            if u(): number()
        counts=[u(),u()]
        u()
        membership={}
        for kind in ('integrated','shortTerm','lra'):
            windows=[(u(),number()) for _ in range(u())]
            if kind=='shortTerm': continue
            lra=kind=='lra'
            def above(value,gate): return value>=gate if lra else value>gate
            absolute=10**((-70+.691)/10)
            ungated=[energy for _,energy in windows if above(energy,absolute)]
            total=0.0
            for e in ungated: total+=e
            relative=total/len(ungated)/(100 if lra else 10)
            included={end for end,e in windows if above(e,absolute) and above(e,relative)}
            assert len(included)==counts[int(lra)]
            membership[kind]=included
        result.append((header['identity'],membership))
    return result


if __name__=='__main__':
    output={}
    for file in map(Path,sys.argv[1:]):
        (sid,source),(cid,current)=populations(file)
        output[file.name]=dict(source=sid,current=cid,groups={k:dict(
            before=len(source[k]),after=len(current[k]),entered_window_ends=sorted(current[k]-source[k]),
            left_window_ends=sorted(source[k]-current[k])) for k in source})
    print(json.dumps(output,indent=2))
