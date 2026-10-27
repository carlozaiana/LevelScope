#!/usr/bin/env python3
"""Independent section/global/PCM24 oracle and libsndfile/report verification."""
import argparse, json, math, pathlib, struct
from verify_m10_pcm import inspect, noise
from verify_m10_decoder import verify_bundle, decode_spans
from verify_m11_pcm import operand, f32

def section_db(rows, rate, frame):
    for n in range(1, len(rows)):
        a, b = rows[n-1], rows[n]
        half = min(rate//20, (a[1]-a[0])//2, (b[1]-b[0])//2)
        start, end = b[0]-half, b[0]+half
        if frame == start: return a[2]
        if frame == end-1: return b[2]
        if start < frame < end-1: return a[2]+(b[2]-a[2])*((frame-start)/(end-start-1))
    for begin,end,gain in rows:
        if begin <= frame < end: return gain
    raise ValueError('Missing section frame')

def verify(folder, library, long=False):
    bundle=folder/'export';result=verify_bundle(bundle,library)
    model=json.loads((bundle/'report.json').read_text());assert model['schema']=='levelscope-measurement-report-v3'
    artifact=model['artifact'];parsed=inspect(bundle/'programme.wav');channels=parsed['channels'];frames=parsed['frames']
    plan=model['current_processing_plan'];rows=[(int(x['begin_frame']),int(x['end_frame']),float.fromhex(x['gain_hex'])) for x in plan['applied_intervals']]
    assert rows and rows[0][0]==0 and rows[-1][1]==frames and plan['engine']=='section-gain-v1'
    assert all(a[1]==b[0] for a,b in zip(rows,rows[1:]))
    seed=int(artifact['dither_seed_hex'],16);global_gain=10**(float.fromhex(plan['global_gain_hex'])/20)
    raw=None if long else (folder/'independent-input.f32').read_bytes()
    if raw is not None: assert len(raw)==frames*channels*4
    spans=decode_spans(frames)
    if long:
        spans += [(rows[n][0]-2500,5000) for n in (1,500,999)]
    checked=0
    with (bundle/'programme.wav').open('rb') as output:
        for first,count in spans:
            output.seek(parsed['offset']+first*channels*3)
            for start in range(first,first+count,4096):
                n=min(4096,first+count-start);pcm=output.read(n*channels*3)
                for frame in range(start,start+n):
                    db=section_db(rows,parsed['rate'],frame);coefficient=10**(db/20)
                    for ch,source in enumerate(artifact['output_to_input_channels']):
                        source-=1
                        original=f32(operand(0,source,frame,frames)/32768+operand(1,source,frame,frames)/32768) if long else struct.unpack_from('<f',raw,(source*frames+frame)*4)[0]
                        section=original if db==0 else f32(original*coefficient)
                        processed=f32(section*global_gain);expected=round(processed*2**23+noise(seed,frame,ch))
                        offset=((frame-start)*channels+ch)*3;actual=int.from_bytes(pcm[offset:offset+3],'little',signed=True)
                        assert actual==expected,(frame,ch,actual,expected)
                        checked+=1
    result.update(section_global_single_dither_oracle=True,codes_checked=checked,long=long)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=pathlib.Path,required=True);p.add_argument('--library',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--long',action='store_true');a=p.parse_args()
    a.output.write_text(json.dumps(verify(a.folder,a.library,a.long),indent=2)+'\n')
