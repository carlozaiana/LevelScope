#!/usr/bin/env python3
import copy,json,pathlib,tempfile,unittest
from run_m13_qualification import REQUIRED,UI_CASES,MANDATORY,validate_runs,validate_report,validate_ui
from complete_m13_qualification import validate_receipt_files
from collect_m1_environment import sha256
from run_m13_muted_video import measured_audio_identity
from m13_evidence_identity import equivalent,select
class Qualification(unittest.TestCase):
    def test_physical_composition_identity(self):
        digest='a'*64
        single=dict(source=dict(measurement_available=True,measured_content_sha256=digest))
        composed=dict(source=dict(measurement_available=True,measured_programme_definition_digest=digest),composition=dict(definition_digest=digest))
        self.assertEqual(measured_audio_identity(single),digest);self.assertEqual(measured_audio_identity(composed),digest)
        for bad in (dict(source={}),dict(source=dict(measurement_available=False,measured_content_sha256=digest)),dict(source=composed['source'],composition=dict(definition_digest='b'*64)),dict(source=single['source'],composition=composed['composition'])):
            with self.assertRaises(ValueError):measured_audio_identity(bad)
    def test_equivalence_never_accepts_product_or_binary_changes(self):
        old=dict(source_digest='old',measurement_build='build',build_configuration=['same'],compiler='same',juce_commit='same',source_files_sha256={'Source/app.cpp':'same','Tests/run_m13_muted_video.py':'old'},artifacts_sha256={'bin/app':'same','Tests/run_m13_muted_video.py':'old'})
        new=copy.deepcopy(old);new['source_digest']='new';new['source_files_sha256']['Tests/run_m13_muted_video.py']='new';new['artifacts_sha256']['Tests/run_m13_muted_video.py']='new'
        self.assertTrue(equivalent(old,new)['production_and_test_binaries_unchanged'])
        for field,key in (('source_files_sha256','Source/app.cpp'),('source_files_sha256','Tests/SectionProposalTests.cpp'),('artifacts_sha256','bin/app')):
            bad=copy.deepcopy(new);bad[field][key]='changed'
            with self.assertRaises(ValueError):equivalent(old,bad)
        bad=copy.deepcopy(new);bad['measurement_build']='other'
        with self.assertRaises(ValueError):equivalent(old,bad)
        with self.assertRaises(ValueError):select({'new':(new,'hash')},'foreign')
    def test_required_receipts(self):
        for phase,names in REQUIRED.items():
            rows=[dict(name=n,exit_code=0) for n in names];validate_runs(rows,phase)
            for bad in (rows[:-1],rows+rows[:1],[dict(name=r['name'],exit_code=1) for r in rows]):
                with self.assertRaises(ValueError):validate_runs(bad,phase)
    def test_empty_or_foreign_report(self):
        identity=dict(measurement_build='candidate')
        for report in ({},dict(passed=True,build='other',architecture='arm64',checks=[dict(passed=True)]),dict(passed=True,build='candidate',architecture='x86_64',checks=[dict(passed=True)]),dict(passed=True,build='candidate',architecture='arm64',checks=[])):
            with self.assertRaises(ValueError):validate_report(report,identity,'arm64','m13-core')
    def test_incomplete_matrix(self):
        r=dict(passed=True,build='candidate',architecture='arm64',physical_devices_opened=False,video_loaded=False,
            checks=[dict(name=k,passed=True) for k in MANDATORY],cases=[dict(rate=rate,channels=ch,mode=mode) for rate in (44100,48000,96000,192000) for ch in (2,6) for mode in range(4 if ch==2 else 6)])
        identity=dict(measurement_build='candidate');validate_report(r,identity,'arm64','m13-core')
        for key,value in (('cases',r['cases'][:-1]),('checks',r['checks'][:-1]),('video_loaded',True),('physical_devices_opened',True)):
            with self.assertRaises(ValueError):validate_report(dict(r,**{key:value}),identity,'arm64','m13-core')
    def test_ui_identity_and_completeness(self):
        identity=dict(source_digest='source',artifacts_sha256={'LevelScope Standalone.app/Contents/MacOS/LevelScope Standalone':'app'})
        r=dict(passed=True,source_digest='source',app_sha256='app',architecture='arm64',device_closed=True,routes_muted=True,cases=[dict(id=k,passed=True,evidence='observed') for k in UI_CASES]);validate_ui(r,identity,'arm64')
        for field,value in (('cases',r['cases'][:-1]),('app_sha256','wrong'),('device_closed',False),('architecture','x86_64')):
            changed=copy.deepcopy(r);changed[field]=value
            with self.assertRaises(ValueError):validate_ui(changed,identity,'arm64')
    def test_missing_long_metrics(self):
        with self.assertRaises(ValueError):validate_report(dict(passed=True,build='candidate',architecture='arm64',checks=[dict(passed=True)]),dict(measurement_build='candidate'),'arm64','m13-long')
    def test_long_budget_failure(self):
        v=dict(frames='518400000',rate=48000,channels=6,intervals=1000,capture_cases=8,private_passes=3,total_ms=100,cancel_ms=10,ui_p95_ms=15,additional_peak_rss_bytes=1,session_bytes=1,report_bytes=1,stages_ms=dict(save_ms=1,parse_ms=1))
        r=dict(passed=True,build='candidate',architecture='arm64',physical_devices_opened=False,video_loaded=False,checks=[dict(passed=True)],long_workflow=v);identity=dict(measurement_build='candidate');validate_report(r,identity,'arm64','m13-long')
        for key,value in (('intervals',999),('private_passes',4),('ui_p95_ms',101),('cancel_ms',float('nan')),('total_ms',14400001)):
            bad=copy.deepcopy(r);bad['long_workflow'][key]=value
            with self.assertRaises(ValueError):validate_report(bad,identity,'arm64','m13-long')
    def test_receipt_tampering(self):
        with tempfile.TemporaryDirectory() as name:
            root=pathlib.Path(name);file=root/'test.log';file.write_text('passed')
            rows=[dict(name='test',log_sha256=sha256(file))];validate_receipt_files(root,rows)
            file.write_text('changed')
            with self.assertRaises(ValueError):validate_receipt_files(root,rows)
            with self.assertRaises(ValueError):validate_receipt_files(root,[dict(name='test')])
if __name__=='__main__':unittest.main()
