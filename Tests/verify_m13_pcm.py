#!/usr/bin/env python3
"""M13 independent PCM24 check; the long fixture changes only two LFE impulse samples."""
import argparse, json, math, pathlib, struct
from verify_m10_pcm import inspect, noise
from verify_m10_decoder import verify_bundle, decode_spans
from verify_m11_pcm import operand, f32
from verify_m12_pcm import section_db, verify as verify_short

def verify(folder, library, long=False):
    if not long:return verify_short(folder,library)
    bundle=folder/'export';result=verify_bundle(bundle,library)
    model=json.loads((bundle/'report.json').read_text());assert model['schema']=='levelscope-measurement-report-v3'
    artifact=model['artifact'];parsed=inspect(bundle/'programme.wav');channels=parsed['channels'];frames=parsed['frames']
    assert frames==518400000 and channels==6 and parsed['rate']==48000
    plan=model['current_processing_plan'];rows=[(int(x['begin_frame']),int(x['end_frame']),float.fromhex(x['gain_hex'])) for x in plan['applied_intervals']]
    assert rows and rows[0][0]==0 and rows[-1][1]==frames and plan['engine']=='section-gain-v1'
    assert all(a[1]==b[0] for a,b in zip(rows,rows[1:]))
    seed=int(artifact['dither_seed_hex'],16);global_gain=10**(float.fromhex(plan['global_gain_hex'])/20)
    impulse=frames//2+48000
    spans=decode_spans(frames)+[(impulse-2048,4096)]+[(r[0]-2500,5000) for r in rows[1:]]
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
                        original=1.5 if source==3 and frame==impulse else f32(operand(0,source,frame,frames)/32768+operand(1,source,frame,frames)/32768)
                        section=original if db==0 else f32(original*coefficient)
                        processed=f32(section*global_gain);expected=round(processed*2**23+noise(seed,frame,ch))
                        offset=((frame-start)*channels+ch)*3;actual=int.from_bytes(pcm[offset:offset+3],'little',signed=True)
                        assert actual==expected,(frame,ch,actual,expected)
                        checked+=1
    result.update(section_global_single_dither_oracle=True,codes_checked=checked,long=True,impulse_frame=impulse)
    return result
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=pathlib.Path,required=True);p.add_argument('--library',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--long',action='store_true');a=p.parse_args()
    a.output.write_text(json.dumps(verify(a.folder,a.library,a.long),indent=2)+'\n')
