#!/usr/bin/env python3
"""Negative qualification-gate tests: missing evidence must never qualify M14."""
import copy, unittest
from run_m14_qualification import REQUIRED, MANDATORY, UI_CASES, validate_runs, validate_report, validate_ui, validate_legacy

class Gates(unittest.TestCase):
    identity={'source_digest':'source','measurement_build':'build','artifacts_sha256':{'LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone':'binary'}}
    def report(self):
        return dict(passed=True,build='build',architecture='arm64',physical_devices_opened=False,video_loaded=False,
                    checks=[dict(name=n,passed=True) for n in MANDATORY],
                    cases=[dict(rate=r,role=c,frames=r*10+1,mask={1:1,2:2,3:4}[c],private_passes=3) for r in (44100,48000,96000,192000) for c in (1,2,3)])
    def test_required_runs(self):
        for phase,names in REQUIRED.items():
            rows=[dict(name=n,exit_code=0) for n in names];validate_runs(rows,phase)
            for invalid in (rows[:-1],rows+[rows[0]],[dict(r,exit_code=1) for r in rows]):
                with self.assertRaises(ValueError):validate_runs(invalid,phase)
    def test_old_group_helper_is_not_rejection_evidence(self):
        baseline=dict(passed=False,architecture='arm64',build='5fd55bb1a81c80b4b556715936b6873bccae29a5d04c00b0890a9b43296a20d2',failure='M8 core migration remains interleaved: Session core read; media verification required')
        applied=dict(baseline,failure='M8 core migration remains interleaved: Applied section plan does not match this programme');validate_legacy(applied,baseline)
        for bad in (dict(applied,failure='prepare grouped qualification session'),baseline,dict(applied,build='wrong')):
            with self.assertRaises(ValueError):validate_legacy(bad,baseline)
    def test_core_identity_and_matrix(self):
        r=self.report();validate_report(r,self.identity,'m14-core')
        for key,value in [('passed',False),('build','old'),('architecture','x86_64'),('physical_devices_opened',True),('video_loaded',True),('checks',[]),('cases',r['cases'][:-1])]:
            with self.assertRaises(ValueError):validate_report(dict(r,**{key:value}),self.identity,'m14-core')
        for key,value in [('mask',4),('frames',0),('private_passes',4)]:
            bad=copy.deepcopy(r);bad['cases'][0][key]=value
            with self.assertRaises(ValueError):validate_report(bad,self.identity,'m14-core')
    def test_long_budgets_and_no_overclaim(self):
        r=self.report();v=dict(frames='1036800000',rate=192000,channels=1,intervals=1000,capture_cases=8,three_hour_qualified=False,rf64_end_to_end_qualified=False,private_passes=3,total_ms=1,cancel_ms=1,ui_p95_ms=1,additional_peak_rss_bytes=1,session_bytes=1,report_bytes=1,stages_ms={'save_ms':1,'parse_ms':1});r['long_workflow']=v
        validate_report(r,self.identity,'m14-long')
        for key,value in [('frames','2073600000'),('capture_cases',7),('three_hour_qualified',True),('rf64_end_to_end_qualified',True),('private_passes',4),('total_ms',7200001),('cancel_ms',1001),('ui_p95_ms',101),('additional_peak_rss_bytes',268435457),('session_bytes',134217729),('report_bytes',16777217),('total_ms',float('nan'))]:
            with self.assertRaises(ValueError):validate_report(dict(r,long_workflow=dict(v,**{key:value})),self.identity,'m14-long')
    def test_cinema_bounded_identity(self):
        r=self.report();c=dict(frames='269875200',rate=48000,channels=6,original_sha256='1fb7f93a33d600bf0a62d63e642e5703f03aec8e6920f8a6e510cf65f61fdeaf',bounded_reads=[{}]*6,full_proposal_export_repeated=False);r['cinema']=c
        validate_report(r,self.identity,'m14-cinema')
        for key,value in [('frames','1'),('original_sha256','wrong'),('bounded_reads',[{}]*5),('full_proposal_export_repeated',True)]:
            with self.assertRaises(ValueError):validate_report(dict(r,cinema=dict(c,**{key:value})),self.identity,'m14-cinema')
    def test_ui_required_observations(self):
        r=dict(passed=True,source_digest='source',architecture='arm64',device_closed=True,routes_muted=True,video_loaded=False,app_sha256='binary',cases=[dict(id=n,passed=True,evidence='observed') for n in UI_CASES]);validate_ui(r,self.identity)
        for key,value in [('cases',r['cases'][:-1]),('app_sha256','wrong'),('source_digest','old'),('device_closed',False),('routes_muted',False),('video_loaded',True)]:
            with self.assertRaises(ValueError):validate_ui(dict(r,**{key:value}),self.identity)
        r['cases'][0]['evidence']=''
        with self.assertRaises(ValueError):validate_ui(r,self.identity)

if __name__=='__main__':unittest.main()
