#!/usr/bin/env python3
"""Independent composition export oracle plus libsndfile/report integrity checks."""
import argparse,json,math,pathlib,struct
from verify_m10_pcm import inspect,noise
from verify_m10_decoder import verify_bundle,decode_spans

def f32(value):return struct.unpack('<f',struct.pack('<f',value))[0]
def operand(component,channel,frame,frames):
    if frame==frames-1:return 1000+channel*100 if component==0 else 500+channel*50
    period=128 if component==0 else 64
    x=(400 if component==0 else 180)*(channel+1)*math.sin((frame%period)*2*math.pi/period)
    return math.floor(x+.5) if x>=0 else math.ceil(x-.5)

def verify(folder,library,long=False):
    bundle=folder/'export';result=verify_bundle(bundle,library)
    model=json.loads((bundle/'report.json').read_text());composition=model['composition']
    assert len(composition['components'])==2 and composition['confirmation']==composition['definition_digest']
    artifact=model['artifact'];parsed=inspect(bundle/'programme.wav');channels=parsed['channels'];frames=parsed['frames']
    seed=int(artifact['dither_seed_hex'],16);gain=10**(model['processing_gain_db']/20);checked=0
    original=None if long else (folder/'independent-composed-oracle.wav').open('rb')
    try:
        with (bundle/'programme.wav').open('rb') as output:
            for first,count in decode_spans(frames):
                output.seek(parsed['offset']+first*channels*3)
                if original:original.seek(44+first*channels*4)
                for start in range(first,first+count,4096):
                    n=min(4096,first+count-start);raw=output.read(n*channels*3)
                    samples=None if long else struct.unpack('<'+'f'*(n*channels),original.read(n*channels*4))
                    for i in range(n*channels):
                        frame=start+i//channels;channel=i%channels
                        summed=f32(operand(0,channel,frame,frames)/32768+operand(1,channel,frame,frames)/32768) if long else samples[i]
                        expected=round(f32(summed*gain)*2**23+noise(seed,frame,channel))
                        actual=int.from_bytes(raw[i*3:i*3+3],'little',signed=True)
                        assert actual==expected,(frame,channel,actual,expected)
                        checked+=1
    finally:
        if original:original.close()
    result.update(composition_oracle_codes=checked,independent_sum_and_single_final_dither=True,long_fixture=long)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=pathlib.Path,required=True);p.add_argument('--library',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--long',action='store_true');a=p.parse_args()
    a.output.write_text(json.dumps(verify(a.folder,a.library,a.long),indent=2)+'\n')
