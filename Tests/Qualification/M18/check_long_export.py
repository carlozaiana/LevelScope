"""Independent streaming hashes, RF64 parser and PCM24 oracle for the combined job."""
import argparse
import importlib.util
import json
import math
from pathlib import Path
import struct

helper=Path(__file__).resolve().parents[2]/'Profiling/ShortAcceptance/check_export.py'
spec=importlib.util.spec_from_file_location('pcm_oracle',helper)
q=importlib.util.module_from_spec(spec);spec.loader.exec_module(q)

def chunks(path):
    with path.open('rb') as f:
        tag,size,kind=struct.unpack('<4sI4s',f.read(12));assert kind==b'WAVE' and tag in (b'RIFF',b'RF64')
        rf64=tag==b'RF64';assert size==(0xffffffff if rf64 else path.stat().st_size-8)
        end=path.stat().st_size;found={};ds=None
        while f.tell()<end:
            name,length=struct.unpack('<4sI',f.read(8));assert name not in found
            if name==b'ds64':
                assert rf64 and length==28
                ds=struct.unpack('<QQQI',f.read(28));assert ds[0]+8==end and ds[3]==0
                found[name]=(f.tell()-28,28);continue
            if name==b'data' and rf64:assert length==0xffffffff and ds;length=ds[1]
            assert f.tell()+length+(length&1)<=end
            found[name]=(f.tell(),length);f.seek(length+(length&1),1)
        assert f.tell()==end and b'fmt ' in found and b'data' in found
        return found,rf64,ds

def run(root):
    result=json.loads((root/'result.json').read_text());assert result['passed'] and all(c['passed'] for c in result['checks'])
    long=result['long_workflow'];frames=int(long['frames']);full=frames==2073600000
    bundle=root/'export';manifest=json.loads((bundle/'manifest.json').read_text());report=json.loads((bundle/'report.json').read_text());art=report['artifact']
    assert set(manifest['files_sha256'])=={'programme.wav','report.json','report.html'}
    for name,digest in manifest['files_sha256'].items():assert q.digest(bundle/name)==digest
    original=root/'mono-duration-192k.wav';path=bundle/'programme.wav'
    source_hash=q.digest(original);assert source_hash==report['source']['measured_content_sha256']==report['current']['measured_content_sha256']
    assert art['sha256']==long['artifact_sha256']==manifest['files_sha256']['programme.wav']
    assert report['ordered_input_roles']==['C'] and art['output_to_input_channels']==[1]
    assert art['evaluation']['verdict']==report['current']['evaluation']['verdict']=='PASS'
    wav,rf,ds=chunks(path);src,srf,sds=chunks(original);assert rf==srf==full
    assert wav[b'data'][1]==src[b'data'][1]==frames*3
    if full:assert ds[2]==sds[2]==frames and wav[b'data'][1]>2**32
    plan=report['current_processing_plan'];plan['global_gain_db']=float.fromhex(plan['global_gain_hex'])
    for row in plan['applied_intervals']:row['gain_db']=float.fromhex(row['gain_hex'])
    impulse=frames//2+192000;selected={0,1,frames-1,impulse-1,impulse,impulse+1}
    selected.update(range(0,frames,max(1,frames//1000)))
    crossing=(2**32-wav[b'data'][0])//3
    if full:selected.update(range(crossing-4,crossing+5))
    for t in plan['transitions']:
        a,b=int(t['begin_frame']),int(t['end_frame'])
        for middle in (a,(a+b)//2,b):selected.update((middle-1,middle,middle+1))
    beyond=0;seed=int(art['dither_seed_hex'],16)
    with path.open('rb') as out,original.open('rb') as source:
        out.seek(wav[b'fmt '][0]);fmt=out.read(wav[b'fmt '][1])
        assert struct.unpack('<HHIIHHHHI',fmt[:24])==(65534,1,192000,576000,3,24,22,24,4)
        assert fmt[24:]==bytes.fromhex('0100000000001000800000aa00389b71')
        for frame in sorted(selected):
            if not 0<=frame<frames:continue
            source.seek(src[b'data'][0]+frame*3);code=int.from_bytes(source.read(3),'little',signed=True)
            expected_source=6710886 if frame==impulse else round(.08*8388608*math.sin((frame%192)*2*math.pi/192))
            assert code==expected_source,(frame,code,expected_source)
            local=q.f32(code/2**23*10**(q.local_db(plan,frame)/20));value=q.f32(local*10**(plan['global_gain_db']/20))
            expected=round(value*2**23+q.noise(seed,frame,0));out.seek(wav[b'data'][0]+frame*3);actual=int.from_bytes(out.read(3),'little',signed=True)
            assert actual==expected,(frame,actual,expected)
            beyond+=wav[b'data'][0]+frame*3>2**32
    assert beyond>0 if full else beyond==0
    return dict(passed=True,frames=frames,rf64=rf,independent_samples=len(selected),samples_beyond_4_gib=beyond,source_sha256=source_hash,artifact_sha256=art['sha256'],three_hour=full)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    report=run(a.directory);a.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
