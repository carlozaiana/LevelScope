#!/usr/bin/env python3
import json
import datetime
import pathlib
import tempfile
import unittest
from unittest.mock import patch
import run_comparison as runner

class RunnerTests(unittest.TestCase):
    def test_device_guard(self):
        p = runner.PROFILES['high-sierra']
        device = {'_name':p['device'],'coreaudio_device_srate':44100,'coreaudio_device_output':2}
        runner.check_device([device],p)
        for invalid in [[],[device,device],[dict(device,coreaudio_device_srate=48000)]]:
            with self.assertRaises(RuntimeError):runner.check_device(invalid,p)
    def test_profiler_metadata_is_not_audio_configuration(self):
        import plistlib
        device={'_name':'Built-in Output','coreaudio_device_srate':44100,'coreaudio_device_output':2}
        a=plistlib.dumps([{'_items':[device],'_SPResponseTime':datetime.datetime(2026,10,24)}])
        b=plistlib.dumps([{'_items':[device],'_SPResponseTime':datetime.datetime(2026,10,25)}])
        with patch.object(runner.subprocess,'check_output',side_effect=[a,b]):
            self.assertEqual(runner.audio_snapshot(),runner.audio_snapshot())
    def test_fixed_muted_commands(self):
        package=pathlib.Path('/package'); identity={'executables':{'control':'control','diagnostic':'diagnostic'}}
        for p in runner.PROFILES.values():
            for role in runner.ORDER:
                argv=runner.command(package,identity,p,role,pathlib.Path('/new-output'))
                self.assertIn('--muted',argv)
                self.assertEqual(argv[argv.index('--seconds')+1],'100')
                self.assertEqual('--trace' in argv,role=='diagnostic')
    def test_preview_never_launches(self):
        with patch.object(runner,'verify_package',return_value={'source_digest':'test'}),patch.object(runner,'sha256',return_value='test'),patch.object(runner,'audio_snapshot',return_value=[{'recorded': datetime.datetime(2026,10,24)}]),patch.object(runner,'check_device'),patch.object(runner.platform,'machine',return_value='arm64'),patch.object(runner.subprocess,'check_output',return_value='15.7.9'),patch.object(runner.subprocess,'run') as launch:
            self.assertEqual(runner.main(['--package','/package','--platform','sequoia']),0)
            with self.assertRaises(RuntimeError):runner.main(['--package','/package','--platform','sequoia','--run'])
            launch.assert_not_called()
    def test_integrity_and_exclusive_output(self):
        with tempfile.TemporaryDirectory() as directory:
            p=pathlib.Path(directory); (p/'binary').write_bytes(b'original')
            runner.write_new(p/'identity.json',{'artifacts_sha256':{'binary':runner.sha256(p/'binary')}})
            runner.verify_package(p)
            with self.assertRaises(FileExistsError):runner.write_new(p/'identity.json',{})
            (p/'binary').write_bytes(b'changed')
            with self.assertRaises(RuntimeError):runner.verify_package(p)
    def test_trace_gap_and_overflow(self):
        with tempfile.TemporaryDirectory() as directory:
            p=pathlib.Path(directory)/'trace.json'
            def event(ns,kind):return [ns,kind,3,0,1,0,0,-1,-1,-1,-1,0,1,0]
            data={'schema':'levelscope-paint-trace-v1','valid':True,'ui_dropped':0,'display_dropped':0,
                'allocated_buffer_bytes':100,'nominal_hz':60,'display_id':1,
                'ui_events':[event(0,5),event(1000000,6),event(50000000,5),event(52000000,6),event(53000000,10)],
                'display_events':[event(0,8),event(16666667,8)]}
            p.write_text(json.dumps(data)); summary=runner.trace_summary(p)
            self.assertEqual(summary['paint_max_interval_ms'],50)
            self.assertEqual(summary['full_paint_max_ms'],2)
            data['ui_dropped']=1;p.write_text(json.dumps(data))
            with self.assertRaises(RuntimeError):runner.trace_summary(p)

if __name__=='__main__':unittest.main()
