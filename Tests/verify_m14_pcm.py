#!/usr/bin/env python3
"""Independent mono container/PCM/dither and report checks; no historical runner."""
import argparse,json,math,pathlib,struct
from verify_m10_pcm import inspect,noise
from verify_m10_decoder import verify_bundle,decode_spans
from verify_m12_pcm import section_db

def f32(x):return struct.unpack('<f',struct.pack('<f',x))[0]

def raw(root):
    rows=json.loads((root/'raw-oracles.json').read_text());expected={(r,m,n) for r in (44100,48000,96000,192000) for m in (1,2,4) for n in (1,2,101,102)}
    assert len(rows)==48 and {(r['rate'],r['mask'],r['frames']) for r in rows}==expected
    checked=0;changed=False
    for row in rows:
        info=inspect(root/row['file']);assert info['channels']==1 and all(info[k]==row[k] for k in ('rate','mask','frames'))
        samples=struct.unpack('<'+'f'*row['frames'],(root/row['capture']).read_bytes())
        with (root/row['file']).open('rb') as file:file.seek(info['offset']);pcm=file.read(info['bytes']);pad=file.read()
        assert pad==(b'\0' if info['bytes']&1 else b'')
        for frame,sample in enumerate(samples):
            actual=int.from_bytes(pcm[frame*3:frame*3+3],'little',signed=True)
            assert actual==round(sample*2**23+noise(int(row['seed_hex'],16),frame,0))
            changed|=actual/2**23!=sample;checked+=1
    assert changed
    return dict(passed=True,cases=len(rows),codes_checked=checked,zero_gain_is_not_bit_identical=True)

def verify(folder,library,long=False):
    bundle=folder/'export';result=verify_bundle(bundle,library,full_decode=not long)
    model=json.loads((bundle/'report.json').read_text());artifact=model['artifact'];info=inspect(bundle/'programme.wav')
    assert info['channels']==1 and artifact['output_to_input_channels']==[1]
    plan=model['current_processing_plan'];rows=[(int(x['begin_frame']),int(x['end_frame']),float.fromhex(x['gain_hex'])) for x in plan['applied_intervals']]
    assert rows and rows[0][0]==0 and rows[-1][1]==info['frames'] and plan['engine']=='section-gain-v1'
    assert all(a[1]==b[0] for a,b in zip(rows,rows[1:]))
    seed=int(artifact['dither_seed_hex'],16);gain=10**(float.fromhex(plan['global_gain_hex'])/20)
    samples=None if long else (folder/'independent-input.f32').read_bytes()
    if samples is not None:assert len(samples)==info['frames']*4
    if long:
        assert info['rate']==192000 and info['frames']==1036800000 and info['container']=='RIFF'
        impulse=info['frames']//2+192000
        spans=decode_spans(info['frames'])+[(impulse-2048,4096)]+[(r[0]-10000,20000) for r in rows[1:]]
    else:spans=[(0,info['frames'])]
    table=[round(.08*8388608*math.sin(n*2*math.pi/192)) for n in range(192)]
    checked=0
    with (bundle/'programme.wav').open('rb') as output:
        for first,count in spans:
            assert first>=0 and first+count<=info['frames']
            output.seek(info['offset']+first*3)
            for start in range(first,first+count,4096):
                n=min(4096,first+count-start);pcm=output.read(n*3);assert len(pcm)==n*3
                for frame in range(start,start+n):
                    original=f32((6710886 if frame==impulse else table[frame%192])/2**23) if long else struct.unpack_from('<f',samples,frame*4)[0]
                    db=section_db(rows,info['rate'],frame);section=original if db==0 else f32(original*10**(db/20))
                    processed=f32(section*gain);expected=round(processed*2**23+noise(seed,frame,0))
                    offset=(frame-start)*3;actual=int.from_bytes(pcm[offset:offset+3],'little',signed=True)
                    assert actual==expected,(frame,actual,expected)
                    checked+=1
    result.update(section_global_single_dither_oracle=True,codes_checked=checked,long=long,
                  coverage='bounded long locations and complete hash' if long else 'every short decoded/output PCM sample',
                  three_hour_qualified=False,rf64_end_to_end_qualified=False)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--folder',type=pathlib.Path,required=True);p.add_argument('--library',type=pathlib.Path,required=True);p.add_argument('--output',type=pathlib.Path,required=True);p.add_argument('--long',action='store_true');p.add_argument('--raw',action='store_true');a=p.parse_args()
    result=raw(a.folder) if a.raw else verify(a.folder,a.library,a.long)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
