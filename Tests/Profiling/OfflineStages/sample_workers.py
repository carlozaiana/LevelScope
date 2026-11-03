"""Separate sampling run, excluded from wall-time statistics. Store host-specific traces in scratch."""
from pathlib import Path
import argparse, json, re, subprocess

def extract(report):
    graph=report.split('Call graph:',1)[1].split('Total number in stack',1)[0]
    counts={}
    worker=None
    for line in graph.splitlines():
        if re.match(r'^\s*\d+ Thread_',line):
            match=re.search(r'^\s*(\d+) Thread_\d+.*?\bls-(source|current|fresh|reused|import)\b',line)
            worker=match[2] if match else None
            if worker:
                values=counts.setdefault(worker,{'worker_stacks':0,'meter_push':0,'add_peak':0,'collector_true_peak':0,'sha_full_block':0})
                values['worker_stacks']+=int(match[1])
            continue
        if not worker: continue
        frame=re.search(r'\b(\d+)\s+(.+)\s+\(in ',line)
        if frame:
            for key,pattern in [('meter_push','OfflineMeasurement::push('),('add_peak','OfflineMeasurement::addPeak('),('collector_true_peak','SectionPeakCollector::truePeak('),('sha_full_block','juce::SHA256Processor::processFullBlock(')]:
                if pattern in frame[2]:counts[worker][key]+=int(frame[1])
    for job in ('source','current','fresh','reused'):
        assert counts[job]['worker_stacks']>100 and 0<counts[job]['add_peak']<=counts[job]['worker_stacks'],(job,counts)
    return counts

def main():
    p=argparse.ArgumentParser()
    for name in ('bin','fixtures','scratch','destination'):p.add_argument(name,type=Path)
    a=p.parse_args()
    dest=a.destination.resolve();dest.mkdir(exist_ok=True)
    assert not any(dest.iterdir()), 'use a fresh results directory'
    for n in (2,6):
        path=a.scratch.resolve()/f'worker-{n}.sample.txt'
        child_report=a.scratch.resolve()/f'sampled-{n}.json'
        with (a.scratch/f'sampled-{n}.log').open('w') as log:
            child=subprocess.Popen([str(a.bin.resolve()/'sample-control'),str(a.fixtures.resolve()/f'programme-{n}.wav'),str(child_report)],stdout=log,stderr=subprocess.STDOUT)
            result=subprocess.run(['/usr/bin/sample',str(child.pid),'120','3','-mayDie','-file',str(path)],capture_output=True,text=True)
            assert child.wait()==0 and path.is_file(),result.stderr
        timed=json.loads((dest.parent/'runs'/f'trial-1-{n}-control.json').read_text())
        sampled=json.loads(child_report.read_text())
        assert sampled['identities']==timed['identities'] and sampled['base_build']==timed['base_build']
        counts=extract(path.read_text(errors='replace'))
        (dest/f'hotspots-{n}.json').write_text(json.dumps(counts,indent=2)+'\n')
        print(n,counts,flush=True)
if __name__=='__main__':main()
