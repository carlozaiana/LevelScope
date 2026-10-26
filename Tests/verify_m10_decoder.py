#!/usr/bin/env python3
"""Independent libsndfile interoperability and report integrity (Python 3.8+)."""
import argparse, ctypes as C, hashlib, json, math, pathlib, struct
from html.parser import HTMLParser
from verify_m10_pcm import inspect, noise

class Info(C.Structure):
    _fields_=[('frames',C.c_int64),('samplerate',C.c_int),('channels',C.c_int),('format',C.c_int),('sections',C.c_int),('seekable',C.c_int)]
class Html(HTMLParser):
    def __init__(self):super().__init__(convert_charrefs=True);self.inside=False;self.text=''
    def handle_starttag(self,tag,attrs):
        assert tag not in ('script','iframe','img','link','audio','video')
        assert not any(k.startswith('on') or k in ('src','href') for k,v in attrs)
        if tag=='pre':self.inside=True
    def handle_endtag(self,tag):
        if tag=='pre':self.inside=False
    def handle_data(self,data):
        if self.inside:self.text+=data

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        while True:
            b=f.read(1024*1024)
            if not b:break
            h.update(b)
    return h.hexdigest()

def decode_spans(frames):
    assert frames>0
    if frames<1000000:return [(0,frames)]
    starts=sorted({0,min(123456789,frames-4096),frames//2,frames-4096})
    return [(first,min(4096,frames-first)) for first in starts]

def verify_bundle(folder, library, long_fixture=False):
    lib=C.CDLL(str(library));lib.sf_version_string.restype=C.c_char_p
    lib.sf_open.argtypes=[C.c_char_p,C.c_int,C.POINTER(Info)];lib.sf_open.restype=C.c_void_p
    lib.sf_readf_int.argtypes=[C.c_void_p,C.POINTER(C.c_int),C.c_int64];lib.sf_readf_int.restype=C.c_int64
    lib.sf_seek.argtypes=[C.c_void_p,C.c_int64,C.c_int];lib.sf_seek.restype=C.c_int64
    lib.sf_close.argtypes=[C.c_void_p];lib.sf_close.restype=C.c_int
    model=json.loads((folder/'report.json').read_text());manifest=json.loads((folder/'manifest.json').read_text())
    html=Html();html.feed((folder/'report.html').read_text());assert json.loads(html.text)==model
    assert manifest['kind']=='verified-pcm-export' and set(manifest['files_sha256'])=={'programme.wav','report.json','report.html'}
    for name,digest in manifest['files_sha256'].items():
        assert name in ('programme.wav','report.json','report.html') and sha(folder/name)==digest
    artifact=model['artifact'];path=folder/'programme.wav';parsed=inspect(path)
    assert int(artifact['file_bytes'])==path.stat().st_size and int(artifact['frames'])==parsed['frames'] and artifact['sample_rate']==parsed['rate'] and int(artifact['channel_mask'],16)==parsed['mask']
    assert artifact['sha256']==manifest['files_sha256']['programme.wav'] and artifact['verification_gain_db']==0
    assert artifact['dither_applications']==1 and artifact['measurement_available']
    assert model['current']['evaluation']['verdict']=='PASS' and artifact['evaluation']['verdict']=='PASS'
    info=Info();handle=lib.sf_open(str(path).encode(),0x10,C.byref(info));assert handle
    frames_tested=0;codes_tested=0
    try:
        assert (info.frames,info.channels,info.samplerate)==(parsed['frames'],parsed['channels'],parsed['rate'])
        assert info.format&0xffff==3 # SF_FORMAT_PCM_24
        spans=decode_spans(info.frames)
        with path.open('rb') as f:
            for first,count in spans:
                assert lib.sf_seek(handle,first,0)==first
                f.seek(parsed['offset']+first*info.channels*3)
                for start in range(first,first+count,4096):
                    n=min(4096,first+count-start);buffer=(C.c_int*(n*info.channels))();assert lib.sf_readf_int(handle,buffer,n)==n
                    raw=f.read(n*info.channels*3);assert len(raw)==n*info.channels*3
                    for idx,value in enumerate(buffer):
                        actual=int.from_bytes(raw[idx*3:idx*3+3],'little',signed=True)
                        assert value==actual*256
                        if long_fixture:
                            frame=start+idx//info.channels;ch=idx%info.channels
                            pcm16=(10000+ch*1000) if frame==info.frames-1 else round(500*(ch+1)*math.sin((frame%128)*2*math.pi/128))
                            processed=struct.unpack('<f',struct.pack('<f',pcm16/32768*10**(model['processing_gain_db']/20)))[0]
                            expected=round(processed*2**23+noise(int(artifact['dither_seed_hex'],16),frame,ch))
                            assert actual==expected,(frame,ch,actual,expected)
                        codes_tested+=1
                    frames_tested+=n
    finally:assert lib.sf_close(handle)==0
    return dict(passed=True,decoder=lib.sf_version_string().decode(),library_sha256=sha(library),artifact_sha256=artifact['sha256'],
                format=parsed,frames_decoded=frames_tested,codes_exact=codes_tested,long_known_source_oracle=long_fixture,
                report_manifest_html_agree=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--library',required=True,type=pathlib.Path);p.add_argument('--bundle',required=True,type=pathlib.Path);p.add_argument('--long-fixture',action='store_true');p.add_argument('--output',required=True,type=pathlib.Path);a=p.parse_args()
    result=verify_bundle(a.bundle,a.library,a.long_fixture);a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
