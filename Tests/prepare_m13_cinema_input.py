#!/usr/bin/env python3
"""Qualification-only lossless rewrap, not a product RF64 importer/repair claim."""
import argparse,hashlib,json,pathlib,struct,time

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser();p.add_argument('--original',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);a=p.parse_args()
    original=a.original.resolve();before=original.stat();a.output.mkdir(parents=True,exist_ok=False)
    with original.open('rb') as source:
        assert source.read(12)==b'RF64\xff\xff\xff\xffWAVE'
        ds=None;fmt=None
        while True:
            header=source.read(8);assert len(header)==8
            kind,size=struct.unpack('<4sI',header);start=source.tell()
            if kind==b'ds64':
                assert ds is None and 28<=size<=4096
                value=source.read(size);ds=struct.unpack_from('<QQQI',value);assert ds[3]==0
            elif kind==b'fmt ':
                assert fmt is None and size==40;fmt=source.read(size)
            elif kind==b'data':
                assert ds and fmt and size==0xffffffff;offset=start;data_bytes=ds[1];break
            assert start+size+(size&1)<=before.st_size
            source.seek(start+size+(size&1))
        assert fmt==bytes.fromhex('feff060080bb0000002f0d0012001800160018003f0000000100000000001000800000aa00389b71')
        assert data_bytes>0 and data_bytes%18==0 and offset+data_bytes<=before.st_size
        frames=data_bytes//18;assert frames<=48000*10800
        header=b'RF64'+struct.pack('<I',0xffffffff)+b'WAVEds64'+struct.pack('<IQQQI',28,data_bytes+96,data_bytes,frames,0)
        header+=b'fmt '+struct.pack('<I',40)+fmt+b'data'+struct.pack('<I',0xffffffff)
        assert len(header)==104
        output=a.output/'cinema-pcm24-payload-identical.wav';full=hashlib.sha256();payload=hashlib.sha256()
        source.seek(0);full.update(source.read(offset));remaining=data_bytes
        with output.open('xb') as dest:
            dest.write(header)
            while remaining:
                block=source.read(min(1024*1024,remaining));assert block
                full.update(block);payload.update(block);dest.write(block);remaining-=len(block)
            for block in iter(lambda:source.read(1024*1024),b''):full.update(block)
        after=original.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
    verified=hashlib.sha256()
    with output.open('rb') as f:
        assert f.read(104)==header
        for block in iter(lambda:f.read(1024*1024),b''):verified.update(block)
    assert verified.hexdigest()==payload.hexdigest()
    evidence=dict(original=str(original),original_sha256=full.hexdigest(),original_bytes=before.st_size,original_mtime_ns=before.st_mtime_ns,
        output=str(output.resolve()),output_sha256=digest(output),payload_sha256=payload.hexdigest(),payload_verified_identical=True,
        source_payload_offset=offset,source_payload_bytes=data_bytes,original_ds64=list(ds),trailing_nonpayload_bytes=before.st_size-offset-data_bytes,
        channels=['L','R','C','LFE','Ls','Rs'],sample_rate=48000,bits=24,frames=frames,range=[0,frames],
        transformation='Canonical qualified RF64 header around identical PCM payload; no resampling, gain, decode/reencode, dither or original modification.',
        limitation='Original extended ds64/container metadata is outside the product-qualified subset. This test-only rewrap does not expand importer support.',
        external_reference='User-reported Nuendo statistics, 28 October 2026; exact version/settings/range/content equivalence unconfirmed. Historical Cubase 10 screenshot remains supplementary.')
    (a.output/'provenance.json').write_text(json.dumps(evidence,indent=2)+'\n');print(output)
if __name__=='__main__':main()
