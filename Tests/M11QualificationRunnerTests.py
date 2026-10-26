#!/usr/bin/env python3
"""Fail-closed coverage tests for the explicit M11 qualification inventory."""
import copy,unittest
from run_m11_qualification import REQUIRED,UI_CASES,MANDATORY_CHECKS,validate_report,validate_runs,validate_ui

class RequiredEvidence(unittest.TestCase):
    def core(self):
        return dict(passed=True,build='build',architecture='arm64',phase='core',physical_devices_opened=False,
            checks=[dict(name=n,passed=True) for n in MANDATORY_CHECKS],cases=[dict(rate=r,kind='composition-%d-%s'%(p,k),gain_db=g) for r in (44100,48000,96000,192000) for p in range(3) for k in ('interleaved','grouped') for g in (0,3.125,-6.02059991328)])
    def test_core_matrix(self):
        r=self.core();validate_report(r,{'measurement_build':'build'},'arm64','core')
        for key,value in [('passed',False),('build','other'),('architecture','x86_64'),('physical_devices_opened',True),('cases',r['cases'][:-1]),('checks',[]),('phase','development-smoke')]:
            with self.subTest(key=key),self.assertRaises(ValueError):validate_report(dict(r,**{key:value}),{'measurement_build':'build'},'arm64','core')
        d=copy.deepcopy(r);d['cases'][-1]=d['cases'][0]
        with self.assertRaises(ValueError):validate_report(d,{'measurement_build':'build'},'arm64','core')
        d=copy.deepcopy(r);d['checks'].pop()
        with self.assertRaises(ValueError):validate_report(d,{'measurement_build':'build'},'arm64','core')
    def test_inventory(self):
        for phase,names in REQUIRED.items():
            rows=[dict(name=n,exit_code=0) for n in names];validate_runs(rows,phase)
            for removed in range(len(rows)):
                with self.assertRaises(ValueError):validate_runs(rows[:removed]+rows[removed+1:],phase)
            with self.assertRaises(ValueError):validate_runs(rows+rows[:1],phase)
            with self.assertRaises(ValueError):validate_runs(rows+[dict(name='recursive-historical-runner',exit_code=0)],phase)
            d=copy.deepcopy(rows);d[0]['exit_code']=1
            with self.assertRaises(ValueError):validate_runs(d,phase)
    def test_ui(self):
        r=dict(passed=True,architecture='arm64',source_digest='source',device_closed=True,routes_muted=True,cases=[dict(id=n,passed=True,evidence='observed receipt') for n in UI_CASES]);identity={'source_digest':'source'};validate_ui(r,identity,'arm64')
        for key,value in [('source_digest','wrong'),('device_closed',False),('cases',r['cases'][:-1])]:
            with self.assertRaises(ValueError):validate_ui(dict(r,**{key:value}),identity,'arm64')
        r['cases'][0]['evidence']=''
        with self.assertRaises(ValueError):validate_ui(r,identity,'arm64')
    def test_long_and_video(self):
        r=dict(passed=True,build='build',architecture='arm64',phase='long',physical_devices_opened=False,checks=[dict(name='long',passed=True)],long_composition=dict(frames='518400000',sample_rate=48000,channels=6,components=2,annotations=1000,capture_cases=6,artifact_bytes=9331200104,ui_callbacks=10,stages_ms=dict(save_ms=1,parse_ms=1),total_ms=1,cancel_ms=1,ui_p95_ms=1,rss_growth_bytes=1,session_bytes=1,report_bytes=1));identity={'measurement_build':'build'};validate_report(r,identity,'arm64','long')
        for key,value in [('frames','100'),('artifact_bytes',104),('capture_cases',5),('total_ms',10800001),('ui_p95_ms',101),('rss_growth_bytes',float('nan'))]:
            d=copy.deepcopy(r);d['long_composition'][key]=value
            with self.assertRaises(ValueError):validate_report(d,identity,'arm64','long')
        v=dict(r,phase='video',video_cases=6,video_selections=[dict(file=file,domain=domain,offset=offset,audio_frame=frame,selection_ms=1,player_skipped_frames=0,physical_paint_test=False) for file in ('cfr-exact.mov','vfr-exact.mov','422-with-tone.mov') for domain in ('Source','Current State') for offset in ('-1001/30000','7/3') for frame in ('256','72256','144256','256','192000')]);validate_report(v,identity,'arm64','video')
        with self.assertRaises(ValueError):validate_report(dict(v,video_cases=5),identity,'arm64','video')
        d=copy.deepcopy(v);d['video_selections'][-1]=d['video_selections'][0]
        with self.assertRaises(ValueError):validate_report(d,identity,'arm64','video')
        e=dict(r,phase='exchange',exchange_cases=[dict(preset=p,representation=k,origin_platform='x86_64',receiving_platform='arm64',foreign_cache_reused=False,remeasured_current_pass=True) for p in range(3) for k in ('grouped','interleaved')]);validate_report(e,identity,'arm64','exchange')
        for key,value in [('origin_platform','arm64'),('foreign_cache_reused',True),('remeasured_current_pass',False)]:
            d=copy.deepcopy(e);d['exchange_cases'][0][key]=value
            with self.assertRaises(ValueError):validate_report(d,identity,'arm64','exchange')
        with self.assertRaises(ValueError):validate_report(dict(e,exchange_cases=e['exchange_cases'][:-1]),identity,'arm64','exchange')
if __name__=='__main__':unittest.main()
